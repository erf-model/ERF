 .. role:: cpp(code)
    :language: c++

 .. _Performance:

GPU weak scaling
================

The plot shows weak scaling of the ABL application with the Smagorinsky LES model using A100 GPUs on the Perlmutter system at NERSC.
The domain size is **amr.n_cell = 128 128 512** for a **single GPU**; this is progressively scaled up to **2048 1024 512** for **128 GPUs**.
This test uses all 4 GPUs per node with GPU-aware MPI communication and runs 100 time steps.

The times shown includes both initialization and time evolution. All I/O and diagnostic calculations are turned off.

.. image:: figures/scaling.png
  :width: 400
