.. _Glossary:

Glossary
========

This page collects the acronyms and abbreviations that appear throughout the ERF
documentation.  Each entry gives the expansion, a one-line description, and --
where the documentation discusses the topic in more detail -- a link to the
relevant section.

.. glossary::
   :sorted:

   ABL
      Atmospheric Boundary Layer -- the lowest portion of the atmosphere, which is
      directly influenced by the surface.  See :ref:`sec:ABLDriverInputs`.

   AGL
      Above Ground Level -- a height measured from the local terrain surface
      rather than from sea level.  See :ref:`sec:Plotfile2DReference`.

   AMR
      Adaptive Mesh Refinement -- the use of locally finer grids to resolve
      selected regions of the domain.  See :ref:`MeshRefinement`.

   AMReX
      The block-structured adaptive mesh refinement (AMR) framework, developed for
      exascale architectures, on which ERF is built.  See :ref:`subsec:AMReX`.

   AMR-Wind
      An incompressible wind-energy flow solver, also built on AMReX, that ERF can
      drive through boundary-plane data exchange.  See :ref:`CouplingToAMRWind`.

   API
      Application Programming Interface -- here, the source-level interface of the
      ERF classes and functions documented by Doxygen.  See :ref:`doxygen_link`.

   CCN
      Cloud Condensation Nuclei -- aerosol particles on which water vapor condenses
      to form cloud droplets.  See :ref:`Microphysics`.

   CFL
      Courant-Friedrichs-Lewy -- the dimensionless number that limits the stable
      explicit time step.  See :ref:`inputs-time-step`.

   CI
      Continuous Integration -- the automated builds and tests that run on each
      pull request.  See :ref:`Testing`.

   CPM
      Cell Perturbation Method -- a technique that adds temperature perturbations
      near an inflow boundary to accelerate the development of turbulence.  See
      :ref:`sec:InflowTurbulenceGeneration`.

   CPU
      Central Processing Unit -- the conventional processor on which ERF runs when
      no accelerator backend is enabled.  See :ref:`sec:build:overview`.

   CSV
      Comma-Separated Values -- the plain-text tabular format used for several ERF
      sampling and time-series outputs.  See :ref:`inputs-data-sampling-outputs`.

   CUDA
      Compute Unified Device Architecture -- NVIDIA's programming model, used as
      one of the AMReX GPU backends.  See :ref:`building:configuration`.

   DNS
      Direct Numerical Simulation -- a mode in which all dynamically relevant
      scales of motion are resolved rather than modeled.  See :ref:`DNSvsLES`.

   E3SM
      Energy Exascale Earth System Model -- the DOE Earth system model from which
      ERF borrows several physics packages.  See :ref:`sec:build:library`.

   EAMxx
      E3SM Atmosphere Model in C++ -- the C++ rewrite of the E3SM atmosphere model
      that supplies the SHOC implementation used by ERF.  See
      :ref:`sec:build:library`.

   EB
      Embedded Boundary -- the cut-cell representation of solid geometry such as
      buildings.  See :ref:`inputs-embedded-boundary-eb-tuning`.

   EDMF
      Eddy-Diffusivity Mass-Flux -- a boundary-layer closure that combines local
      diffusion with a nonlocal mass-flux plume contribution.  See :ref:`MYNNEDMF`.

   EKAT
      E3SM Kokkos Application Toolkit -- the utility library, bundled as a
      submodule, required by the EAMxx-derived physics.  See
      :ref:`sec:build:library`.

   EOS
      Equation of State -- the thermodynamic relation between pressure, density and
      temperature.  See :ref:`GoverningEquations`.

   ERA5
      ECMWF Reanalysis version 5 -- the global reanalysis dataset commonly used to
      drive hindcast simulations.  See :ref:`sec:HindCast`.

   ERF
      Energy Research and Forecasting -- the atmospheric model documented here.
      See :ref:`GettingStarted`.

   EWP
      Explicit Wake Parametrization -- a wind-farm parameterization that represents
      turbine wakes through an analytically prescribed velocity deficit.  See
      :ref:`explicit-wake-parametrization-ewp-model`.

   FFT
      Fast Fourier Transform -- used by the spectral Poisson solvers available for
      the anelastic formulation.  See :ref:`subsec:LinearSolvers`.

   GAD
      Generalized Actuator Disk -- the most detailed of the ERF wind-turbine models,
      which computes blade-element forces around the rotor disk.  See
      :ref:`generalized_actuator_disk_model`.

   GCM
      General Circulation Model -- a global climate or weather model; the term
      appears in the name of the RRTMGP radiation package.  See :ref:`Radiation`.

   GPU
      Graphics Processing Unit -- the accelerator hardware targeted through the
      CUDA, HIP and SYCL backends.  See :ref:`sec:build:overview`.

   HDF5
      Hierarchical Data Format version 5 -- an optional output format for ERF
      plotfiles.  See :ref:`sec:Plotfiles`.

   HIP
      Heterogeneous-Compute Interface for Portability -- AMD's programming model,
      used as one of the AMReX GPU backends.  See :ref:`building:configuration`.

   HPC
      High-Performance Computing -- the large parallel machines on which production
      ERF simulations are run.  See :ref:`sec:build:hpc`.

   IF
      Immersed Forcing -- the representation of terrain or buildings in which large
      body forces drive the velocity to zero inside solid cells, selected with
      ``erf.terrain_type`` or ``erf.buildings_type`` = ``ImmersedForcing``.  See
      :ref:`sec:ImmersedForcingInputs`.

   INAS
      Ice Nucleation Active Surface site -- the density of sites on an insoluble
      aerosol particle that can trigger freezing.  See :ref:`sec:SuperDroplets`.

   INP
      Ice Nucleating Particle -- an aerosol particle capable of initiating ice
      formation.  See :ref:`Microphysics`.

   LAD
      Leaf Area Density -- the one-sided leaf area per unit volume used by the
      forest canopy drag model.  See :ref:`sec:Forest`.

   LAI
      Leaf Area Index -- the one-sided leaf area per unit ground area used by the
      land surface model.  See :ref:`SLM`.

   LES
      Large-Eddy Simulation -- a mode in which the large turbulent eddies are
      resolved and the smaller ones are modeled.  See :ref:`DNSvsLES`.

   LSM
      Land Surface Model -- the component that evolves soil and surface state and
      supplies energy and moisture fluxes at a land lower boundary.  See
      :ref:`inputs-land-surface-model`, :ref:`SLM` and :ref:`CouplingToNoahMP`.

   LW
      Longwave -- the thermal infrared part of the radiative spectrum.  See
      :ref:`Radiation`.

   MOST
      Monin-Obukhov Similarity Theory -- the surface-layer theory used to relate
      surface fluxes to the flow in the lowest grid cells.  See
      :ref:`sec:surface_layer`.

   MPI
      Message Passing Interface -- the library used for distributed-memory
      parallelism.  See :ref:`sec:build:overview`.

   MRF
      Medium Range Forecast -- a nonlocal PBL scheme that originated in the NCEP MRF
      model.  See :ref:`MRFPBL`.

   MYJ
      Mellor-Yamada-Janjic -- a local, TKE-based PBL scheme.  See :ref:`MYJ`.

   MYNN
      Mellor-Yamada-Nakanishi-Niino -- a family of TKE-based PBL closures, of which
      Level 2.5 is implemented in ERF.  See :ref:`MYNN25`.

   NERSC
      National Energy Research Scientific Computing Center -- the DOE facility that
      hosts Perlmutter.  See :ref:`sec:hpc:guides`.

   NetCDF
      Network Common Data Form -- the self-describing file format used for WRF
      inputs, boundary files and several ERF outputs.  See :ref:`sec:Initialization`.

   Noah-MP
      Noah Multi-Parameterization -- the community land surface model that ERF can
      couple to.  See :ref:`CouplingToNoahMP`.

   NREL
      National Renewable Energy Laboratory -- the source of the reference wind
      turbine specifications used by the wind farm models.  See
      :ref:`sec:WindFarmModels`.

   OpenMP
      Open Multi-Processing -- the directive-based interface used for shared-memory
      threading.  See :ref:`building:configuration`.

   PBL
      Planetary Boundary Layer -- the turbulent layer of the atmosphere adjacent to
      the surface, and by extension the schemes that parameterize it.  See
      :ref:`PBLschemes`.

   PBS
      Portable Batch System -- one of the job schedulers used on supported HPC
      systems.  See :ref:`sec:build:hpc`.

   RANS
      Reynolds-Averaged Navier-Stokes -- a modeling mode in which all turbulence is
      parameterized rather than resolved.  See :ref:`RANS`.

   RH
      Relative Humidity -- the ratio of the vapor pressure to its saturation value.
      See :ref:`sec:Plotfile3DReference`.

   RHS
      Right-Hand Side -- the collected source and flux-divergence terms advanced by
      the time integrator.  See :ref:`TimeAdvance`.

   RICO
      Rain In Cumulus Over the Ocean -- a precipitating shallow-cumulus benchmark
      case.  See :ref:`Microphysics`.

   RK2
   RK3
      Second- and third-order Runge-Kutta -- the explicit multi-stage time
      integrators used by ERF.  See :ref:`TimeAdvance`.

   RRTMGP
      Rapid Radiative Transfer Model for GCMs, Parallel -- the radiative transfer
      package used by ERF.  See :ref:`Radiation`.

   SAD
      Simplified Actuator Disk -- a wind-turbine model that applies thrust over the
      rotor disk without resolving the blades.  See
      :ref:`actuator_disk_model_simplified`.

   SAI
      Stem Area Index -- the one-sided stem area per unit ground area used by the
      land surface model.  See :ref:`SLM`.

   SAM
      System for Atmospheric Modeling -- the cloud-resolving model whose
      single-moment bulk microphysics scheme is available in ERF.  See
      :ref:`Microphysics`.

   SatAdj
      Saturation Adjustment -- the simplest ERF moisture model, which condenses or
      evaporates water to remove supersaturation.  See :ref:`Microphysics`.

   SBM
      Spectral-Bin Microphysics -- a scheme that resolves the hydrometeor size
      distribution over discrete bins.  See :ref:`sec:SpectralBinMicrophysics`.

   SDM
      Super-Droplet Method -- a particle-based, probabilistic microphysics model.
      See :ref:`sec:SuperDroplets`.

   SEB
      Surface Energy Balance -- the budget of radiative, sensible, latent and ground
      heat fluxes that sets the surface temperature.  See :ref:`sec:IBSEB`.

   SGS
      Sub-Grid Scale -- the motions smaller than the mesh spacing, which must be
      modeled rather than resolved.  See :ref:`DNSvsLES`.

   SHOC
      Simplified Higher-Order Closure -- a unified turbulence and cloud macrophysics
      parameterization taken from EAMxx.  See :ref:`SHOC`.

   SLM
      Simplified Land Model -- ERF's lightweight land surface model.  See :ref:`SLM`.

   SST
      Sea Surface Temperature -- the prescribed or coupled temperature of the ocean
      surface.  See :ref:`inputs-ocean-surface-model`.

   STF
      Smoothed Terrain Following -- a terrain-fitted mesh in which small-scale
      terrain features are progressively damped with height.  See :ref:`sec:Meshing`.

   STL
      Stereolithography -- the triangulated-surface file format used to describe
      building geometry for embedded boundaries.  See :ref:`sec:EBBuildingsSTL`.

   SW
      Shortwave -- the solar part of the radiative spectrum.  See :ref:`Radiation`.

   SYCL
      The Khronos C++ programming model for heterogeneous processors, used as the
      AMReX backend for Intel GPUs.  See :ref:`building:configuration`.

   TKE
      Turbulent Kinetic Energy -- the kinetic energy of the turbulent fluctuations,
      carried as a prognostic variable by several closures.  See :ref:`DNSvsLES`.

   TSK
      Skin Temperature -- the surface skin temperature field, named after the
      corresponding WRF variable.  See :ref:`sec:Plotfile2DReference`.

   UPP
      Unified Post Processor -- the NOAA post-processing package whose sea-level
      pressure reduction ERF reproduces.  See :ref:`sec:Plotfile2DReference`.

   UTC
      Coordinated Universal Time -- the time standard used for simulation start
      times and provenance metadata.  See :ref:`sec:Provenance`.

   UUID
      Universally Unique Identifier -- the identifier stamped into plotfiles and
      checkpoints to tie outputs to the run that produced them.  See
      :ref:`sec:Provenance`.

   WDM6
      WRF Double-Moment 6-class -- a bulk microphysics scheme that predicts both
      mass and number concentration for several hydrometeor species.  See
      :ref:`Microphysics`.

   WENO
      Weighted Essentially Non-Oscillatory -- a family of high-order advection
      schemes designed to limit spurious oscillations.  See
      :ref:`inputs-advection-schemes`.

   WoA
      Witch of Agnesi -- the analytic hill profile used as a canonical
      flow-over-terrain test case.  See :ref:`sec:Meshing`.

   WPS
      WRF Preprocessing System -- the toolchain that generates the ``wrfinput`` and
      ``wrfbdy`` files ERF can be initialized from.  See :ref:`sec:Initialization`.

   WRF
      Weather Research and Forecasting model -- the widely used mesoscale model that
      ERF is frequently compared to and can be initialized from.  See
      :ref:`ERFvsWRF`.

   WSM6
      WRF Single-Moment 6-class -- a bulk microphysics scheme that predicts the mass
      of six hydrometeor species.  See :ref:`Microphysics`.

   WW3
      WAVEWATCH III -- the ocean surface wave model that ERF can couple to.  See
      :ref:`CouplingToWW3`.

   YSU
      Yonsei University -- a nonlocal PBL scheme with countergradient mixing.  See
      :ref:`YSUPBL`.
