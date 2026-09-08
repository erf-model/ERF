# Pre-processing for Data Assimilation

This directory contains the pre-processing code that generates the background initial conditions for ensemble-based data assimilation.

## Overview

The code takes a fine-mesh "ground truth" plotfile and a coarse-mesh plotfile, then interpolates the fine-mesh data onto the coarse mesh and writes out a custom-written binary file with the coarse data, that can be read-in by the data assimilation code. 

**Concept:** Satellite observations capture coarse representations of flow features ie. the features are present, but at a lower resolution. To mimic this, this tool generates a coarsened version of the ground truth and exports it to a custom binary data file used as input for ensemble simulations.

## Build and Run

1. **Compile:**
   ```bash
   make -j8
   ```
2. `vars.txt` should contain the list of variables to be included in the coarsened data.

3. **Execute:**
   ```bash
   ./out -- <fine-mesh-plt-file> <coarse-mesh-plt-file>
   ```
Note: The coarse-mesh plotfile **must contain only 1 box**. Ensure it is run on a single MPI rank with a sufficiently large `max_grid_size`.
