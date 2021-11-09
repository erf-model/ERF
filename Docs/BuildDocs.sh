#!/bin/bash

# exit immediately on error
set -e

# To avoid issues when calling the script from different directories
# sets the directory to the location of the script
cd $(dirname $0)

# This short script builds both the doxygen and sphinx documentation

# Define pretty colors
YEL='\033[0;33m'
GRN='\033[1;32m'
NC='\033[0m'


# Build doxygen documents -- configuration parameters contained in Doxyfile.in
# Errors can be viewed in the file: doxy.log
doxygen Doxyfile.in



# Build sphinx documents
cd sphinx_doc/
make clean          # fixes occasional unexpected behavior
make html


echo ""
echo -e "-- ${GRN}Building Docs Complete${NC} --"
echo -e "See ${YEL}doxygen_output/html/${NC} for doxygen html files"
echo -e "See ${YEL}sphinx_doc/_build/html/${NC} for sphinx html files"
