 .. role:: cpp(code)
    :language: c++

 .. _Performance:

Scaling results
================

To assess the performance of the ERF code on CPUs and GPUs, strong and weak scaling studies were completed. The ABL simulation was chosen for the scaling studies since it is representative of a typical atmospheric flow simulation. The scaling was performed on the CPU and GPU nodes on Perlmutter -- a supercomputer hosted at the National Energy Research Scientific Computing Center (NERSC), which is located at Lawrence Berkeley National Laboratory (Berkeley Lab). A single CPU node consists of 128 AMD EPYC 7763 (Milan) cores while a GPU node has 64 AMD EPYC 7763 (Milan) cores and 4 NVIDIA A100 (Ampere) GPUs. In all simulations, the domain has lengths :math:`(L_x, L_y, L_z) = (2048, 2048, 1024)` [m].

For the strong scaling study, a fixed mesh size of :math:`(N_x, N_y, N_z) = (512, 512, 256)` is employed, and a 128-core simulation is taken as the reference. Subsequent simulations were completed by doubling the number of cores until 4096. For CPU-only simulations, the observed strong scaling timings are illustrated in the figure below, while the parallelization efficiency, defined as:

.. math::
    E = \left( \frac{T}{T_\text{128}} \right) \left( \frac{128}{N} \right) \times 100 \%,

is shown in the figure below. In the above, N is the number of cores and T is the time taken per time step. At 2048 cores, and ~32^3 cells per rank, the parallel efficiency is 69%. On GPU nodes, the number of CPU ranks is held constant at 4 --- i.e., each rank offloads its work to a single GPU. Therefore, the speed-up presented here between GPU and CPU is per node, rather than per rank. On node, speed-ups of 5--15x are achieved up to 16 GPU nodes, as shown in the figure below.

.. figure:: figures/StrongScaling_CPU.png
   :alt: Strong scaling on CPUs
   :name: strong_CPU
   :figwidth: 45%
   :align: left

   **Strong scaling on CPUs. The number of mesh cells per rank is shown in blue.**

.. figure:: figures/ParEff.png
   :alt: Parallelization efficiency for strong scaling
   :name: strong_pareff
   :figwidth: 45%
   :align: left

   **Parallelization efficiency for strong scaling.**

.. figure:: figures/WeakScaling_CPU.png
   :alt: Weak scaling on CPUs
   :name: weak_CPU
   :figwidth: 45%
   :align: left

   **Weak scaling on CPUs. The total number of mesh cells is shown in blue.**

.. figure:: figures/CPUvsGPU.png
   :alt: Comparison of timings on CPU and GPU showing the speed-up factor
   :name: CPUvsGPU
   :figwidth: 45%
   :align: left

   **Comparison of timings on CPU and GPU showing the speed-up factor. We compare a CPU node with 128 ranks to a GPU node with 4 ranks, so there are 32x more points per GPU than per CPU core. Points on the same vertical line represent the same number of nodes.**
