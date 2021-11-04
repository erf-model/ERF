#!/bin/bash


# This short script builds both the doxygen and sphinx documentation

# Define pretty colors
YEL='\033[0;33m'
GRN='\033[1;32m'
NC='\033[0m'


# Build doxygen documents -- configuration parameters contained in Doxyfile.in
doxygen Doxyfile.in

# After build, doxygen html files can be found in `doxygen_output/html/`


# Build sphinx documents
cd sphinx_doc/
make clean          # fixes occasional unexpected behavior
make html 

# After build, sphinx html files can be found in `sphinx_doc/_build/html/`

echo "" 
echo -e "-- ${GRN}Building Docs Complete${NC} --"
echo -e "See ${YEL}doxygen_output/html/${NC} for doxygen html files"
echo -e "See ${YEL}sphinx_doc/_build/html/${NC} for sphinx html files"
