#!/usr/bin/env bash

set -o errexit
set -o nounset
set -o pipefail

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

echo ${host}

build_dir=build_${host}_${CI_PIPELINE_ID}_$(date +%F_%H_%M_%S)

if [[ -n ${modules} ]]
then
    module load ${modules}
fi

# Temporary workaround for CUDA builds:
#  AMReX fcompare seems to not work as expected if compiled with CUDA.
#  This builds a CPU version first and uses that fcompare executable during the
#  testing for the CUDA build
if [[ "${ERF_ENABLE_CUDA}" == "ON" ]]
then
    echo "====================================================="
    echo "Building CPU version first to get fcompare executable"
    echo "====================================================="
    mkdir "${build_dir}_cpu"
    cd "${build_dir}_cpu"
    pwd

    time cmake -DCMAKE_INSTALL_PREFIX:PATH=./install \
         -DCMAKE_CXX_COMPILER:STRING=${CMAKE_CXX_COMPILER:-"mpicxx"} \
         -DCMAKE_C_COMPILER:STRING=${CMAKE_C_COMPILER:-"mpicc"} \
         -DCMAKE_Fortran_COMPILER:STRING=${CMAKE_Fortran_COMPILER:-"mpifort"} \
         -DCMAKE_BUILD_TYPE:STRING=Release \
         -DERF_DIM:STRING=3 \
         -DERF_ENABLE_MPI:BOOL=ON \
         -DERF_ENABLE_CUDA:BOOL=OFF \
         -DERF_ENABLE_TESTS:BOOL=OFF \
         -DERF_ENABLE_FCOMPARE:BOOL=ON \
         -DERF_ENABLE_DOCUMENTATION:BOOL=OFF \
         -DCMAKE_EXPORT_COMPILE_COMMANDS:BOOL=ON ..
    time make -j fcompare

    FCOMPARE_EXE="$(pwd)/Submodules/AMReX/Tools/Plotfile/amrex_fcompare"

    cd ../

    echo "====================================================="
    echo "Using fcompare executable at: ${FCOMPARE_EXE}"
    echo "====================================================="
fi

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
    echo "====================================================="
    echo "Using gold files at: ${ERF_TEST_GOLD_FILES_DIRECTORY}"
    echo "====================================================="
fi
cd -

mkdir ${build_dir}
cd ${build_dir}
pwd

time cmake \
     -DCMAKE_INSTALL_PREFIX:PATH=./install \
     -DCMAKE_CXX_COMPILER:STRING=${CMAKE_CXX_COMPILER:-"mpicxx"} \
     -DCMAKE_C_COMPILER:STRING=${CMAKE_C_COMPILER:-"mpicc"} \
     -DCMAKE_Fortran_COMPILER:STRING=${CMAKE_Fortran_COMPILER:-"mpifort"} \
     -DMPIEXEC_EXECUTABLE="${mpiexec_executable}" \
     -DMPIEXEC_PREFLAGS:STRING="${mpiexec_preflags}" \
     -DCMAKE_BUILD_TYPE:STRING="${build_type}" \
     -DERF_DIM:STRING=3 \
     -DERF_ENABLE_MPI:BOOL=ON \
     -DERF_ENABLE_CUDA:BOOL="${ERF_ENABLE_CUDA}" \
     -DAMReX_CUDA_ARCH:STRING="${CUDA_ARCH:-""}" \
     -DERF_ENABLE_HIP:BOOL="${ERF_ENABLE_HIP:-"OFF"}" \
     -DAMReX_AMD_ARCH:STRING="${AMD_ARCH:-""}" \
     -DERF_ENABLE_TESTS:BOOL=ON \
     -DERF_TEST_NRANKS:STRING=${ERF_TEST_NRANKS:-"4"} \
     -DERF_ENABLE_FCOMPARE:BOOL=ON \
     -DERF_ENABLE_DOCUMENTATION:BOOL=OFF \
     -DFCOMPARE_EXE="${FCOMPARE_EXE:-"$(pwd)/Submodules/AMReX/Tools/Plotfile/amrex_fcompare"}" \
     -DERF_TEST_GOLD_FILES_DIRECTORY="${ERF_TEST_GOLD_FILES_DIRECTORY:-"$(pwd)/../Tests/ERFGoldFiles"}" \
     -DERF_TEST_FCOMPARE_RTOL="${ERF_TEST_FCOMPARE_RTOL:-"5.0e-9"}" \
     -DERF_TEST_FCOMPARE_ATOL="${ERF_TEST_FCOMPARE_ATOL:-"2.0e-10"}" \
     -DCMAKE_EXPORT_COMPILE_COMMANDS:BOOL=ON \
     ..
time make -j ${OMP_NUM_THREADS:-16}
time ctest -VV --output-on-failure
