#!/bin/bash


# This short script builds both the doxygen and sphinx documentation


# Build doxygen documents -- configuration parameters contained in Doxyfile.in
doxygen Doxyfile.in

# After build, doxygen html files can be found in `doxygen_output/html/`


# Build sphinx documents
cd sphinx_doc/
make clean          # fixes occasional unexpected behavior
make html 

# After build, sphinx html files can be found in `sphinx_doc/_build/html/`

