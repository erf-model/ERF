---
title: 'ERF: Energy Research and Forecasting'

tags:
  - C++
  - atmospheric modeling
  - mesoscale
  - microscale
  - wind energy

authors:
  - name: Ann Almgren
    orcid: 0000-0003-2103-312X
    corresponding: true
    affiliation: 1
  - name: Aaron Lattanzi
    orcid: 0000-0002-3372-406X
    affiliation: 1
  - name: Riyaz Haque
    orcid: 0000-0001-8930-5721
    affiliation: 4
  - name: Pankaj Jha
    orcid: 0000-0003-1476-6747
    affiliation: 4
  - name: Branko Kosovic
    orcid: 0000-0002-1746-0746
    affiliation: 6
  - name: Jeffrey Mirocha
    orcid: 0000-0001-5914-6168
    affiliation: 4
  - name: Bruce Perry
    orcid: 0000-0002-9150-8103
    affiliation: 2
  - name: Eliot Quon
    orcid: 0000-0002-8445-5840
    affiliation: 2
  - name: Michael Sanders
    orcid: 0000-0001-6644-2571
    affiliation: 5
  - name: David Wiersema
    orcid: 0000-0001-8452-4095
    affiliation: 4
  - name: Donald Willcox
    orcid: 0000-0003-2300-5165
    affiliation: 1
  - name: Xingqiu Yuan
    orcid: 0000-0002-6146-1260
    affiliation: 3
  - name: Weiqun Zhang
    orcid: 0000-0001-8092-1974
    affiliation: 1

affiliations:
 - name: Lawrence Berkeley National Laboratory
   index: 1
 - name: National Renewable Energy Laboratory
   index: 2
 - name: Argonne National Laboratory
   index: 3
 - name: Lawrence Livermore National Laboratory
   index: 4
 - name: San Diego State University
   index: 5
 - name: National Center for Atmospheric Research
   index: 6

date: 17 February 2023

bibliography: paper.bib
---

# Summary

The Energy Research and Forecasting (ERF) code is a new model that simulates the mesoscale and microscale
dynamics of the atmosphere using the latest high-performance computing architectures.  It employs
hierarchical parallelism using an MPI+X model, where X may be OpenMP on multicore CPU-only systems,
or CUDA, HIP, or SYCL on GPU-accelerated systems.
ERF is built on AMReX [@AMReX:JOSS; @AMReX:IJHPCA],
a block-structured adaptive mesh refinement (AMR) software framework that
provides the underlying performance-portable software infrastructure for block-structured mesh operations. 
The "energy" aspect of ERF indicates that the software has been developed with renewable energy applications in mind.
In addition to being a numerical weather prediction model, ERF is designed to provide a flexible
computational framework for the exploration and investigation of different physics parameterizations
and numerical strategies, and to characterize the flow field that impacts the
ability of wind turbines to extract wind energy.  The ERF development is part of a broader effort
led by the US Department of Energy's Wind Energy Technologies Office.

# ERF Features

### Hydrodynamics Models

ERF solves the fully compressible Navier-Stokes equations for
dry or moist air and includes a planetary boundary layer (PBL)
parameterization as well as subfilter flux parameterizations for
large-eddy simulations (LES). The PBL parameterization is based on
the work of Mellor and Yamada [@PBL:Mellor] and Nakanishi and Niino [@PBL:Nakanishi],
the so-called MYNN model for mesoscale simulations. LES parameterizations
are Smagorinsky-type [@LES:Smagorinsky; @LES:Lilly] and Deardorff [@LES:Deardorff].

### Microphysics Options

Microphysics options in ERF include a warm, non-precipitating model
that evolves cloud water and cloud vapor and a single-moment model [@SAMXX:marat] that evolves precipitating and
nonprecipitating tracers, such as water vapor, rain, ice, snow, and graupel. 
These prognostic variables can track particle evolution through all the important mechanisms of ice and water growth,
including vapor deposition, aggregation, autoconversion, and condensation.

### Time and Space Discretization and Terrain

The time discretization in ERF utilizes a third-order Runge-Kutta scheme with
substepping of perturbational quantities at the acoustic time scale [@FAST:Klemp].
(A non-substepping method is also available as a run-time option.)
The spatial discretization in ERF uses the classic Arakawa C-grid with 
scalar quantities at cell centers and normal velocities at cell faces.
For simulations over complex topography, a terrain-following, height-based
vertical coordinate is employed.  The model includes capability for application
of some common map projections (e.g., Lambert Conformal, Mercator).
The advection terms may be calculated using second- through sixth-order accurate
spatial discretizations, including both centered difference and upwind 
schemes.  Third- and fifth-order weighted essentially non-oscillatory (WENO) advection schemes
are also available for the cell-centered scalars.
ERF supports both static and dynamic (adaptive) mesh refinement,
with subcycling in time at finer levels of refinement.

### Physical Forcings and Boundary Conditions

Physical forcings include Coriolis and geostrophic forcing as well as
Rayleigh damping in the upper regions of the domain.  Lateral boundary
conditions can be specified as periodic, inflow/outflow, or time-varying
values read in from external files in netcdf format generated by the WRF
Preprocessing System (WPS) [@WRF:Skamarock]. The surface boundary condition may be
specified either as a simple wall or by using Monin-Obukhov similarity theory (MOST)
[@MOST:monin; @MOST:van] to model the surface layer. The initial data can
be read from WPS-generated files, reconstructed from 1-d input sounding
data, or specified by the user.

# Statement of need

Most widely used atmospheric modeling codes today do not have the 
ability to use GPU acceleration, which limits their ability to 
efficiently utilize current and next-generation high performance computing 
architectures.  ERF provides an atmospheric modeling capability that runs on the latest high-performance
computing architectures, from laptops to supercomputers, 
whether CPU-only or GPU-accelerated.  In addition, ERF is based on AMReX,
a modern, well-supported AMR library,
which provides a performance portable interface that shields ERF
from most of the detailed changes needed to adapt to new systems.
The active and large developer community contributing to AMReX helps ensure
that ERF will continue to run efficiently as architectures and operating systems
evolve.

To support renewable energy research and development, ERF provides an essential
resource characterization and forensic capability for terrestrial and offshore
applications. For wind energy, ERF includes a standard suite of physical process
parameterizations that supports simulation across weather (meso) and
turbulence-resolving (micro) scales, allowing for efficient downscaling of
flow field information that specifies realistic inflow, surface, and background
conditions for wind farm simulation.  Realistic conditions can include extreme
wind-shear events (e.g., low-level jets), thunderstorms, or tropical cyclones
(e.g., hurricanes). This modeling capability also captures the impacts of clouds
and precipitation, and is similarly applicable to solar farms and hybrid energy
systems.

# Acknowledgements

Funding for this work was provided by the U.S. Department of Energy
Office of Energy Efficiency and Renewable Energy Wind Energy Technologies Office.
We acknowledge the help of the AMReX team
in developing and supporting new AMReX features needed by ERF.
The work at LBNL was supported by the U.S. Department of Energy
under contract No. DE-AC02-05CH11231. 
The work at LLNL was supported by the U.S. Department of Energy
under contract No. DE-AC52-07NA27344.
The contribution of Branko Kosovic was supported by the National Center for Atmospheric Research,
which is a major facility sponsored by the National Science Foundation under Cooperative Agreement No. 1852977.
This work was authored in part by the
National Renewable Energy Laboratory, operated by Alliance for Sustainable Energy, LLC,
for the U.S. Department of Energy (DOE) under Contract No. DE-AC36-08GO28308.
The views expressed in the article do not necessarily represent the views of the DOE or the U.S. Government. The U.S. Government retains and the publisher, by accepting the article for publication, acknowledges that the U.S. Government retains a nonexclusive, paid-up, irrevocable, worldwide license to publish or reproduce the published form of this work, or allow others to do so, for U.S. Government purposes.
A portion of this code development was performed using computational resources sponsored by the
U.S. Department of Energy's Office of Energy Efficiency and Renewable Energy and located
at the National Renewable Energy Laboratory. This development also used resources of the 
National Energy Research Scientific Computing Center (NERSC), 
a U.S. Department of Energy Office of Science User Facility located at 
Lawrence Berkeley National Laboratory.

# References
