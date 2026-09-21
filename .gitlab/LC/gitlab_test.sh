#!/usr/bin/env bash

set -o errexit
set -o nounset
set -o pipefail

# ===========================================================================
# Phase timing instrumentation
# ===========================================================================
# Each phase is wrapped in phase_begin / phase_end. Durations accumulate and
# are printed as a summary table by the EXIT trap, so the table appears even
# when a phase fails under `set -o errexit`.
#
#   CI_TIMING_SECTIONS=false   disable GitLab collapsible-section markers
#   NINJA_SLOWEST_N=15         how many slow build steps to report (0 = off)
#   CTEST_SLOWEST_N=20         how many slow tests to report (0 = off)
# ===========================================================================

CI_TIMING_SECTIONS=${CI_TIMING_SECTIONS:-"true"}
NINJA_SLOWEST_N=${NINJA_SLOWEST_N:-15}
CTEST_SLOWEST_N=${CTEST_SLOWEST_N:-20}

PHASE_NAMES=()
PHASE_SECONDS=()
CURRENT_PHASE=""
CURRENT_PHASE_START=0
SCRIPT_START=$(date +%s)

# ANSI colors for the job log. GitLab renders standard SGR codes, and its own
# built-in sections ("Preparing environment" etc.) use bold green, so C_HEAD
# matches them. Set CI_COLOR=false for plain output.
if [[ "${CI_COLOR:-true}" == "true" ]]; then
    C_HEAD=$'\e[32;1m'   # bold green  - section headers
    C_WARN=$'\e[31;1m'   # bold red    - warnings
    C_DIM=$'\e[36;1m'    # bold cyan   - phase end / summary
    C_OFF=$'\e[0m'
else
    C_HEAD="" C_WARN="" C_DIM="" C_OFF=""
fi

fmt_hms() {
    local t=${1}
    printf '%02d:%02d:%02d' $((t / 3600)) $(((t % 3600) / 60)) $((t % 60))
}

# GitLab renders these as collapsible sections with their own duration badge.
section_marker() {
    # $1 = start|end, $2 = section id, $3 = header text (start only)
    [[ "${CI_TIMING_SECTIONS}" == "true" ]] || return 0
    local header=""
    [[ -n "${3:-}" ]] && header="${C_HEAD}${3}${C_OFF}"
    printf '\e[0Ksection_%s:%s:%s\r\e[0K%s\n' "${1}" "$(date +%s)" "${2}" "${header}"
}

phase_id() {
    echo "${1}" | tr '[:upper:] ()/-' '[:lower:]____'
}

phase_begin() {
    CURRENT_PHASE="${1}"
    CURRENT_PHASE_START=$(date +%s)
    section_marker start "$(phase_id "${CURRENT_PHASE}")" "${CURRENT_PHASE}"
    echo "======================================================================"
    echo ">>> PHASE START: ${CURRENT_PHASE}   ($(date '+%F %T'))"
    echo "======================================================================"
}

phase_end() {
    local now elapsed
    now=$(date +%s)
    elapsed=$((now - CURRENT_PHASE_START))
    PHASE_NAMES+=("${CURRENT_PHASE}")
    PHASE_SECONDS+=("${elapsed}")
    echo "${C_DIM}<<< PHASE END:   ${CURRENT_PHASE}   elapsed $(fmt_hms "${elapsed}") (${elapsed}s)${C_OFF}"
    section_marker end "$(phase_id "${CURRENT_PHASE}")"
    CURRENT_PHASE=""
}

# Record a phase that was deliberately not run, so the table stays readable.
phase_skip() {
    PHASE_NAMES+=("${1} [skipped]")
    PHASE_SECONDS+=(0)
    echo "--- PHASE SKIPPED: ${1}"
}

timing_report() {
    local rc=${1:-0}
    local now total width i s pct
    now=$(date +%s)
    total=$((now - SCRIPT_START))

    if [[ -n "${CURRENT_PHASE}" ]]; then
        PHASE_NAMES+=("${CURRENT_PHASE} [FAILED/INCOMPLETE]")
        PHASE_SECONDS+=($((now - CURRENT_PHASE_START)))
        # A phase that died never reached phase_end, so its section_start has
        # no matching section_end. Close it here, otherwise GitLab folds the
        # failure output -- and this summary -- into the still-open section.
        section_marker end "$(phase_id "${CURRENT_PHASE}")"
        CURRENT_PHASE=""
    fi

    width=20
    if ((${#PHASE_NAMES[@]} > 0)); then
        for i in "${!PHASE_NAMES[@]}"; do
            ((${#PHASE_NAMES[i]} > width)) && width=${#PHASE_NAMES[i]}
        done
    fi

    echo
    echo "======================================================================"
    echo "Phase timing summary (script exit status: ${rc})"
    echo "======================================================================"
    printf "%-${width}s  %10s  %9s  %7s\n" "PHASE" "ELAPSED" "SECONDS" "% TOTAL"

    if ((${#PHASE_NAMES[@]} > 0)); then
        for i in "${!PHASE_NAMES[@]}"; do
            s=${PHASE_SECONDS[i]}
            pct=0
            ((total > 0)) && pct=$((100 * s / total))
            printf "%-${width}s  %10s  %9d  %6d%%\n" \
                   "${PHASE_NAMES[i]}" "$(fmt_hms "${s}")" "${s}" "${pct}"
        done
    fi
    printf "%-${width}s  %10s  %9d  %6d%%\n" "TOTAL" "$(fmt_hms "${total}")" "${total}" 100
    echo
    echo "End: $(date)"
}

trap 'timing_report $?' EXIT

# Per-target timings straight out of Ninja's own log. Useful for spotting a
# single expensive link step hiding inside an otherwise healthy build phase.
report_slowest_targets() {
    local dir=${1}
    local n=${2:-${NINJA_SLOWEST_N}}
    local log="${dir}/.ninja_log"

    ((n > 0)) || return 0
    [[ -f "${log}" ]] || return 0

    echo
    echo "--- ${n} slowest build steps ---"
    # .ninja_log is TAB separated, one line per output:
    #   start_ms <TAB> end_ms <TAB> mtime <TAB> output path <TAB> command hash
    # Split on tabs (paths may contain spaces) and keep the longest duration
    # seen per output, in case a target was built more than once.
    awk -F'\t' '
        NR > 1 && NF >= 4 {
            d = ($2 - $1) / 1000.0
            if (d > best[$4]) best[$4] = d
        }
        END { for (o in best) printf "%.1f\t%s\n", best[o], o }
    ' "${log}" \
        | sort -k1,1 -rn \
        | head -n "${n}" \
        | awk -F'\t' '{ printf "  %9.1f s  %s\n", $1, $2 }' \
        || true
    echo "--- (full log: ${log}) ---"
}

# Per-test timings, the ctest counterpart of report_slowest_targets. ctest
# --output-junit writes one <testcase> per test carrying its duration and
# status, which is steadier to parse than the console output. If that file is
# missing, fall back to Testing/Temporary/CTestCostData.txt, which ctest writes
# unprompted for its own scheduling -- no status there, and the figure is an
# average over runs, which equals the single run in a fresh build directory.
report_slowest_tests() {
    local junit=${1}
    local n=${2:-${CTEST_SLOWEST_N}}
    local costfile=${3:-}

    ((n > 0)) || return 0

    echo
    echo "--- ${n} slowest tests ---"

    if [[ -f "${junit}" ]]; then
        # Match attributes by name, not position, and keep names with spaces
        # intact by separating fields with tabs.
        grep -o '<testcase [^>]*' "${junit}" \
            | awk '
                {
                    name = ""; t = -1; st = ""
                    if (match($0, /name="[^"]*"/))   { name = substr($0, RSTART+6, RLENGTH-7) }
                    if (match($0, /time="[^"]*"/))   { t    = substr($0, RSTART+6, RLENGTH-7) + 0 }
                    if (match($0, /status="[^"]*"/)) { st   = substr($0, RSTART+8, RLENGTH-9) }
                    if (name != "" && t >= 0) { printf "%.2f\t%s\t%s\n", t, st, name }
                }' \
            | sort -k1,1 -rn \
            | head -n "${n}" \
            | awk -F'\t' '{ mark = ($2 == "run" ? "" : "  [" $2 "]");
                            printf "  %9.2f s  %s%s\n", $1, $3, mark }' \
            || true
        echo "--- (full results: ${junit}) ---"
    elif [[ -n "${costfile}" && -f "${costfile}" ]]; then
        # "name numRuns averageSeconds" lines, terminated by a "---" line
        awk '$1 == "---" { exit } NF >= 3 { printf "%.2f\t%s\n", $3, $1 }' "${costfile}" \
            | sort -k1,1 -rn \
            | head -n "${n}" \
            | awk -F'\t' '{ printf "  %9.2f s  %s\n", $1, $2 }' \
            || true
        echo "--- (full log: ${costfile}) ---"
    else
        echo "  (no timing data found)"
    fi
}

# Resolve build parallelism and check the step actually has the CPUs it is
# about to use. `nproc` honours the affinity mask the launcher imposed;
# `nproc --all` ignores it. If the launcher confined the step, the build
# oversubscribes and runs far slower than the node allows, with no error
# message anywhere.
detect_parallelism() {
    cores_avail=$(nproc)
    cores_node=$(nproc --all)
    # CI_BUILD_JOBS / CI_LINK_JOBS are set per-machine in the GitLab CI files.
    build_jobs=${CI_BUILD_JOBS:-${cores_avail}}
    link_jobs=${CI_LINK_JOBS:-4}

    echo "CPUs usable / on node: ${cores_avail} / ${cores_node}"
    echo "Build jobs: ${build_jobs}   Link jobs: ${link_jobs}"

    if ((cores_avail < build_jobs)); then
        echo "${C_WARN}WARNING: ${build_jobs} build jobs but only ${cores_avail} of"
        echo "         ${cores_node} CPUs are usable by this step;"
        echo "         this step was confined by the launcher."
        echo "         Flux: run this script as the INITIAL PROGRAM of a"
        echo "         whole-node subinstance (flux alloc -N 1 --exclusive"
        echo "         SCRIPT), not as a job inside the allocation -- an"
        echo "         initial program consumes no resources and is not bound."
        echo "         Slurm: check --cpus-per-task/--exact and --mpibind.${C_OFF}"
    fi
}

# ===========================================================================

echo "Start: $(date)"

echo "========="
echo "GitLab CI"
echo "========="

modules=${MODULE_LIST:-""}
mpiexec_executable=${MPIEXEC_EXECUTABLE:-"srun"}
# If using flux, append "run" after the flux executable path
if [[ "${mpiexec_executable}" == "flux" ]]
then
    mpiexec_executable="$(which ${mpiexec_executable}) run"
    flux jobs
    flux resource list
else
    mpiexec_executable="$(which ${mpiexec_executable})"
fi

mpiexec_preflags=${MPIEXEC_PREFLAGS:-""}
host=$(hostname)
build_type=${BUILD_TYPE:-"Debug"}

RUN_CTEST=${RUN_CTEST:-"true"}
RUN_CTEST=${RUN_CTEST,,}

ERF_ENABLE_CUDA=${ERF_ENABLE_CUDA:-"OFF"}
ERF_ENABLE_HIP=${ERF_ENABLE_HIP:-"OFF"}

echo "HOST: ${host}"
src_dir="${PWD}"
echo "Source directory: ${src_dir}"
build_dir="$(realpath -- "${src_dir}/../build_${host}_${CI_PIPELINE_ID}_${CI_JOB_ID}_$(date +%F_%H_%M_%S)")"
echo "Build directory: ${build_dir}"

phase_begin "Setup modules"

if [[ -n ${modules} ]]
then
    module load ${modules}
fi
module list

phase_end

echo "==========================="
echo "Available build parallelism"
echo "==========================="

detect_parallelism

# Default fcompare executable
FCOMPARE_EXE="${build_dir}/Submodules/AMReX/Tools/Plotfile/amrex_fcompare"

if [[ $RUN_CTEST == "true" ]]; then
    # For GPU builds we use a CPU version of fcompare to compare output files as it
    # can be faster than the GPU version because data does not need to migrate to
    # device memory.
    if [[ "${ERF_ENABLE_CUDA}" == "ON" || "${ERF_ENABLE_HIP}" == "ON" ]]
    then
        phase_begin "Configure fcompare (CPU)"
        cmake \
             -G Ninja \
             -S "${src_dir}" \
             -B "${build_dir}_cpu" \
             -D CMAKE_INSTALL_PREFIX:PATH=./install \
             -D CMAKE_CXX_COMPILER:STRING=${CMAKE_CXX_COMPILER:-"mpicxx"} \
             -D CMAKE_C_COMPILER:STRING=${CMAKE_C_COMPILER:-"mpicc"} \
             -D CMAKE_Fortran_COMPILER:STRING=${CMAKE_Fortran_COMPILER:-"mpifort"} \
             -D CMAKE_BUILD_TYPE:STRING=Release \
             -D ERF_DIM:STRING=3 \
             -D ERF_ENABLE_MPI:BOOL=ON \
             -D ERF_ENABLE_CUDA:BOOL=OFF \
             -D ERF_ENABLE_TESTS:BOOL=OFF \
             -D ERF_ENABLE_UNIT_TESTS:BOOL=OFF \
             -D ERF_ENABLE_FCOMPARE:BOOL=ON \
             -D ERF_ENABLE_DOCUMENTATION:BOOL=OFF \
             -D CMAKE_EXPORT_COMPILE_COMMANDS:BOOL=ON \
             -D ERF_PRECISION:STRING="${ERF_PRECISION^^}" \
             -D ERF_PARTICLES_PRECISION:STRING="${ERF_PARTICLES_PRECISION^^}" \
             -D CMAKE_JOB_POOLS:STRING="link=${link_jobs}" \
             -D CMAKE_JOB_POOL_LINK:STRING=link \
             -D ERF_ENABLE_CRAY_AUTO_FIXES=OFF
        phase_end

        phase_begin "Build fcompare (CPU)"
        cmake --build "${build_dir}_cpu" --parallel "${build_jobs}" --target fcompare
        phase_end

        FCOMPARE_EXE="${build_dir}_cpu/Submodules/AMReX/Tools/Plotfile/amrex_fcompare"
    else
        phase_skip "Configure fcompare (CPU)"
        phase_skip "Build fcompare (CPU)"
    fi
else
    phase_skip "Configure fcompare (CPU)"
    phase_skip "Build fcompare (CPU)"
fi

# Default gold files directory
ERF_TEST_GOLD_FILES_DIRECTORY="${src_dir}/Tests/ERFGoldFiles"
ERF_TEST_ENABLE_EXTRA_SDM_TESTS="OFF"

if [[ $RUN_CTEST == "true" ]]; then
    phase_begin "Clone LC gold files repo"

    # Clone LC gold files repo -- note that we need to grant this repo job
    # token permissions to the gold file repo
    rm -rf erf-llnl-gold-files
    git clone \
        --branch ${CI_GOLD_FILES_GIT_REF:-"main"} --depth 1 \
        https://gitlab-ci-token:${CI_JOB_TOKEN}@lc.llnl.gov/gitlab/erf-model/erf-llnl-gold-files.git
    cd erf-llnl-gold-files
    git log -1
    if [[ -d ${CI_MACHINE} ]]
    then
        ERF_TEST_GOLD_FILES_DIRECTORY="$(pwd)/${CI_MACHINE}"
        ERF_TEST_ENABLE_EXTRA_SDM_TESTS="ON"
        if [[ "${ERF_ENABLE_CUDA}" == "ON" || "${ERF_ENABLE_HIP}" == "ON" ]]
        then
            if [[ -d "${CI_MACHINE}/gpu" ]]; then
                ERF_TEST_GOLD_FILES_DIRECTORY+="/gpu"
            fi
        else
            if [[ -d "${CI_MACHINE}/cpu" ]]; then
                ERF_TEST_GOLD_FILES_DIRECTORY+="/cpu"
            fi
        fi
    fi
    cd -

    phase_end
else
    phase_skip "Clone LC gold files repo"
fi

phase_begin "Configure ERF"

cmake \
     -G Ninja \
     -S "${src_dir}" \
     -B "${build_dir}" \
     -D CMAKE_INSTALL_PREFIX:PATH=./install \
     -D CMAKE_CXX_COMPILER:STRING=${CMAKE_CXX_COMPILER:-"mpicxx"} \
     -D CMAKE_C_COMPILER:STRING=${CMAKE_C_COMPILER:-"mpicc"} \
     -D CMAKE_Fortran_COMPILER:STRING=${CMAKE_Fortran_COMPILER:-"mpifort"} \
     -D MPIEXEC_EXECUTABLE="${mpiexec_executable}" \
     -D MPIEXEC_PREFLAGS:STRING="${mpiexec_preflags}" \
     -D CMAKE_BUILD_TYPE:STRING="${build_type}" \
     -D ERF_DIM:STRING=3 \
     -D ERF_PRECISION:STRING="${ERF_PRECISION^^}" \
     -D ERF_ENABLE_PARTICLES:BOOL=ON \
     -D ERF_PARTICLES_PRECISION:STRING="${ERF_PARTICLES_PRECISION^^}" \
     -D ERF_ENABLE_MPI:BOOL=ON \
     -D ERF_ENABLE_CUDA:BOOL="${ERF_ENABLE_CUDA}" \
     -D CMAKE_CUDA_ARCHITECTURES:STRING="${CUDA_ARCH:-""}" \
     -D ERF_ENABLE_HIP:BOOL="${ERF_ENABLE_HIP:-"OFF"}" \
     -D AMReX_AMD_ARCH:STRING="${AMD_ARCH:-""}" \
     -D ERF_ENABLE_FFT:BOOL="${ERF_ENABLE_FFT:-"OFF"}" \
     -D ERF_ENABLE_FCOMPARE:BOOL=ON \
     -D FCOMPARE_EXE="${FCOMPARE_EXE}" \
     -D ERF_ENABLE_DOCUMENTATION:BOOL=OFF \
     -D ERF_ENABLE_TESTS:BOOL=${RUN_CTEST} \
     -D ERF_ENABLE_UNIT_TESTS:BOOL=${RUN_CTEST} \
     -D ERF_TEST_NRANKS:STRING=${ERF_TEST_NRANKS:-"4"} \
     -D ERF_TEST_GOLD_FILES_DIRECTORY="${ERF_TEST_GOLD_FILES_DIRECTORY}" \
     -D ERF_TEST_ENABLE_EXTRA_SDM_TESTS="${ERF_TEST_ENABLE_EXTRA_SDM_TESTS}" \
     -D ERF_TEST_FCOMPARE_RTOL="${ERF_TEST_FCOMPARE_RTOL:-"5.0e-9"}" \
     -D ERF_TEST_FCOMPARE_ATOL="${ERF_TEST_FCOMPARE_ATOL:-"2.0e-10"}" \
     -D CMAKE_EXPORT_COMPILE_COMMANDS:BOOL=ON \
     -D CMAKE_JOB_POOLS:STRING="link=${link_jobs}" \
     -D CMAKE_JOB_POOL_LINK:STRING=link \
     -D ERF_ENABLE_CRAY_AUTO_FIXES=OFF

phase_end

phase_begin "Build ERF"

cmake --build "${build_dir}" --parallel "${build_jobs}"
report_slowest_targets "${build_dir}"

phase_end

if [[ "$RUN_CTEST" == "true" ]]; then
    phase_begin "Test ERF"

    echo "fcompare executable: ${FCOMPARE_EXE}"
    echo "Gold files directory: ${ERF_TEST_GOLD_FILES_DIRECTORY}"
    echo "Extra SDM tests enabled: ${ERF_TEST_ENABLE_EXTRA_SDM_TESTS}"

    # errexit would abort here on a test failure, before the timings are
    # printed, so take the status by hand and re-raise it below.
    ctest_junit="${build_dir}/ctest_results.xml"
    ctest_rc=0
    ctest --test-dir "${build_dir}" --extra-verbose --output-on-failure \
          --no-tests=error --output-junit "${ctest_junit}" || ctest_rc=$?

    report_slowest_tests "${ctest_junit}" "${CTEST_SLOWEST_N}" \
                         "${build_dir}/Testing/Temporary/CTestCostData.txt"

    # Leave the phase open on failure so the EXIT trap still labels it
    # FAILED/INCOMPLETE, exactly as it did before this reporting was added.
    ((ctest_rc == 0)) || exit "${ctest_rc}"

    phase_end
else
    phase_skip "Test ERF"
fi

# The EXIT trap prints the timing summary and the closing "End:" line.
