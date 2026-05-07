.. _sec:kestrel_build_run:

Kestrel (NREL)
==============

The `Kestrel <https://nrel.github.io/HPC/Documentation/Systems/Kestrel/>`__ cluster is an HPE Cray system with Intel Xeon Sapphire Rapids CPU nodes (104 cores) and NVIDIA H100 GPUs (4 per node).

.. note::
   Kestrel has **separate login nodes for GPU work**. Access GPU login nodes via ``kestrel-gpu.hpc.nrel.gov``. GPU jobs should only be submitted from GPU login nodes.

**Building**

.. tab-set::

   .. tab-item:: GNU Make (CPU)

      .. code-block:: bash

         # Reset to default environment
         module restore

         # Build with Cray compilers
         cd ${ERF_HOME}/Exec
         make realclean
         make -j COMP=cray

   .. tab-item:: GNU Make (GPU)

      .. code-block:: bash

         # Load GPU modules
         module purge
         module load PrgEnv-gnu/8.5.0
         module load cuda/12.3
         module load craype-x86-milan

         # Build
         cd ${ERF_HOME}/Exec
         make realclean
         make -j COMP=gnu USE_CUDA=TRUE

   .. tab-item:: CMake (GPU)

      .. code-block:: bash

         # Load GPU modules
         module purge
         module load PrgEnv-gnu/8.5.0
         module load cuda/12.3
         module load craype-x86-milan

         # Configure and build
         mkdir build && cd build
         cmake -DCMAKE_BUILD_TYPE=Release \
               -DERF_ENABLE_MPI=ON \
               -DERF_ENABLE_CUDA=ON \
               ..
         make -j

.. warning::
   **System updates on Kestrel periodically change required modules.** Verify current module names with ``module avail`` before building.

**Memory Allocation**

Kestrel allows partial node allocations. For memory-intensive operations (e.g., CUDA compilation):

.. code-block:: bash

   # Option 1: Request exclusive node access
   #SBATCH --exclusive

   # Option 2: Request specific memory
   #SBATCH --mem=240G
   # or
   #SBATCH --mem-per-cpu=2G

Without these flags, CUDA compilation may fail due to insufficient memory.

**Performance and Cost Considerations**

GPU node hours on Kestrel are charged at **10× the rate** of CPU node hours. Understanding performance trade-offs is essential for efficient use of allocations.

**Typical Performance Characteristics:**

* GPU nodes (4× H100): **10-20× faster** than CPU nodes (96-104 cores)
* Best efficiency with **>1M cells per GPU**
* Smaller problems may not fully utilize GPU capability

**When to Use GPU vs CPU:**

.. list-table::
   :header-rows: 1
   :widths: 30 35 35

   * - Problem Size
     - GPU Nodes (10× cost)
     - CPU Nodes (1× cost)
   * - < 500K cells/GPU
     - May not justify 10× cost
     - **Recommended**
   * - 500K - 1M cells/GPU
     - Marginal benefit
     - Consider for development
   * - > 1M cells/GPU
     - **Recommended** - cost effective
     - Slower time to solution
   * - > 5M cells/GPU
     - **Excellent utilization**
     - May exceed wall-time limits

**Recommendations:**

1. **Profile your specific case** - Performance varies with physics packages and I/O frequency
2. **Development on CPU** - Use CPU nodes for code development and small test cases
3. **Production on GPU** - Use GPU nodes for production runs with well-sized domains
4. **Monitor utilization** - Check GPU usage with ``nvidia-smi`` to verify saturation

.. note::
   The 10-20× performance gain typically justifies the 10× cost increase for production runs, providing faster time-to-solution and 1-2× better overall cost efficiency measured in allocation units per simulation.
