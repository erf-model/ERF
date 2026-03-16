#!/bin/bash

# Set path to ERF
ERF_DIR="${ERF_DIR:-/home/ERF}"

# Set the path to external in ERF
EXT_DIR=${ERF_DIR}/external

# NOTE: These git commands requires git version > 2.34
git clone --filter=blob:none https://github.com/E3SM-Project/E3SM.git ${EXT_DIR}/E3SM --sparse
cd ${EXT_DIR}/E3SM
git checkout 2eb52d9b3d452e593f3825f9408507cae3519bf6
git sparse-checkout set components/eamxx/src/physics
git sparse-checkout add components/eamxx/src/share
export EAMXX_HOME=${EXT_DIR}/E3SM/components/eamxx/src
