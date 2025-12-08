#!/usr/bin/env bash

set -o errexit
set -o nounset
set -o pipefail

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

ERF_ENABLE_CUDA=${ERF_ENABLE_CUDA:-"OFF"}
ERF_ENABLE_HIP=${ERF_ENABLE_HIP:-"OFF"}

basehost=${host//[[:digit:]]/}

echo "HOST: ${host}"
src_dir="${PWD}"
echo "Source directory: ${src_dir}"
build_dir="$(realpath -- "${src_dir}/../build_${host}_${CI_PIPELINE_ID}_$(date +%F_%H_%M_%S)")"
echo "Build directory: ${build_dir}"

echo "============="
echo "Setup modules"
echo "============="

if [[ -n ${modules} ]]
then
    module load ${modules}
fi
module list

# Default fcompare executable
FCOMPARE_EXE="${build_dir}/Submodules/AMReX/Tools/Plotfile/amrex_fcompare"

# For GPU builds we use a CPU version of fcompare to compare output files as it
# can be faster than the GPU version because data does not need to migrate to
# device memory.
if [[ "${ERF_ENABLE_CUDA}" == "ON" || "${ERF_ENABLE_HIP}" == "ON" ]]
then
    echo "======================="
    echo "Build CPU amrex_fcompre"
    echo "======================="
    time cmake \
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
         -D ERF_ENABLE_FCOMPARE:BOOL=ON \
         -D ERF_ENABLE_DOCUMENTATION:BOOL=OFF \
         -D CMAKE_EXPORT_COMPILE_COMMANDS:BOOL=ON \
         -D ERF_ENABLE_CRAY_AUTO_FIXES=OFF
    time cmake --build "${build_dir}_cpu" --target fcompare
    FCOMPARE_EXE="${build_dir}_cpu/Submodules/AMReX/Tools/Plotfile/amrex_fcompare"
fi
echo "fcompare executable: ${FCOMPARE_EXE}"

echo "========================"
echo "Clone LC gold files repo"
echo "========================"

# Default gold files directory
ERF_TEST_GOLD_FILES_DIRECTORY="${src_dir}/Tests/ERFGoldFiles"

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

echo "Gold files directory: ${ERF_TEST_GOLD_FILES_DIRECTORY}"

echo "============="
echo "Configure ERF"
echo "============="

time cmake \
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
     -D ERF_ENABLE_MPI:BOOL=ON \
     -D ERF_ENABLE_CUDA:BOOL="${ERF_ENABLE_CUDA}" \
     -D AMReX_CUDA_ARCH:STRING="${CUDA_ARCH:-""}" \
     -D ERF_ENABLE_HIP:BOOL="${ERF_ENABLE_HIP:-"OFF"}" \
     -D AMReX_AMD_ARCH:STRING="${AMD_ARCH:-""}" \
     -D ERF_ENABLE_FCOMPARE:BOOL=ON \
     -D FCOMPARE_EXE="${FCOMPARE_EXE}" \
     -D ERF_ENABLE_DOCUMENTATION:BOOL=OFF \
     -D ERF_ENABLE_TESTS:BOOL=ON \
     -D ERF_TEST_NRANKS:STRING=${ERF_TEST_NRANKS:-"4"} \
     -D ERF_TEST_GOLD_FILES_DIRECTORY="${ERF_TEST_GOLD_FILES_DIRECTORY}" \
     -D ERF_TEST_FCOMPARE_RTOL="${ERF_TEST_FCOMPARE_RTOL:-"5.0e-9"}" \
     -D ERF_TEST_FCOMPARE_ATOL="${ERF_TEST_FCOMPARE_ATOL:-"2.0e-10"}" \
     -D CMAKE_EXPORT_COMPILE_COMMANDS:BOOL=ON \
     -D ERF_ENABLE_CRAY_AUTO_FIXES=OFF

echo "========="
echo "Build ERF"
echo "========="

time cmake --build "${build_dir}"

echo "========"
echo "Test ERF"
echo "========"

time ctest --test-dir "${build_dir}" --extra-verbose --output-on-failure

echo "End: $(date)"
