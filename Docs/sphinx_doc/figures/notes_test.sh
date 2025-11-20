#!/bin/bash

git clone --recursive git@github.com:erf-model/ERF
cd ERF

source Build/machines/perlmutter_erf.profile
BACKEND=CUDA  # or HIP or SYCL

git fetch origin pull/2725/head:pr-2725
git checkout pr-2725

# Replace make -j commands in cmake and build scripts
#find Build -name "cmake*.sh" -type f -exec sed -i 's/make -j[0-9]\+/make/g; s/make -j /make /g' {} \;
#find Build -name "build_*.sh" -type f -exec sed -i 's/make -j[0-9]\+/make/g; s/make -j /make /g' {} \;

#ERF_HOME puts things in build and install automatically, for testing we want to be more specific
ERF_HOME=$(pwd) ERF_BUILD_DIR=$(pwd)/build_gpu ERF_INSTALL_DIR=install_gpu ./Build/cmake_with_kokkos_many_${BACKEND,,}.sh
ERF_HOME=$(pwd) ERF_BUILD_DIR=$(pwd)/build_cpu ERF_INSTALL_DIR=install_cpu ./Build/cmake_with_kokkos_many.sh

# Create GPU version of cmake_with_shoc.sh
sed "/ERF_ENABLE_MPI/a\      -DERF_ENABLE_${BACKEND}:BOOL=ON \\\\" Build/cmake_with_shoc.sh > Build/cmake_with_shoc_${BACKEND,,}.sh
chmod +x Build/cmake_with_shoc_${BACKEND,,}.sh

# Create GPU version of build_erf_with_shoc.sh that calls the GPU cmake script
sed "s/cmake_with_shoc.sh/cmake_with_shoc_${BACKEND,,}.sh/" Build/build_erf_with_shoc.sh > Build/build_erf_with_shoc_${BACKEND,,}.sh
chmod +x Build/build_erf_with_shoc_${BACKEND,,}.sh

# Run the GPU build script
./Build/build_erf_with_shoc_${BACKEND,,}.sh

mv build build_shoc_gpu

# Run the CPU build script
./Build/build_erf_with_shoc.sh

mkdir build_tmp
cd build_tmp

# Replace make -j in locally generated cmake scripts
sed -i 's/make -j[0-9]\+/make/g; s/make -j /make /g' ../Build/cmake*.sh

../Build/cmake.sh
make install
make uninstall
make distclean
../Build/cmake_${BACKEND,,}.sh
make distclean
cd ../
