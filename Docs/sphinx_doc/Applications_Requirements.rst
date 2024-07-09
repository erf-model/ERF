.. highlight:: rst

#############################
Applications and Requirements
#############################

Last update: 2021-01-14

Overview
========

The Energy Research and Forecasting (ERF) model addresses a key component of the wind energy simulation supply chain, the downscaling of the energy contained within mesoscale atmospheric flows into the microscale wind plant environment, where turbines generate electricity. The efficient, robust and accurate downscaling of mesoscale flow information, which defines not only the energy available for conversion to electricity, but also many characteristics influencing the extraction of that energy, is crucial for the design and optimization of wind plant layout and operation in realistic, complicated operating conditions, including the more challenging complex terrain and offshore settings defining the majority of future deployment.

Much of the functionality ERF will provide is obtained today from weather prediction models such as the Weather Research and Forecasting (WRF) model. However, neither WRF nor similar models were designed to prioritize accurate prediction of low-level winds, or to be easily configurable for wind energy applications, focusing instead on larger-scale weather parameters such as storm dynamics, precipitation and temperature. While among such models WRF is most amenable to, and has been the most widely adopted into wind energy work flows, these and other limitations of its software architecture and design, including poor load balancing and inefficient data transfer in downscaling configurations, no ability to decompose in the vertical direction which limits parallelizability, and complications arising from its pressure-based vertical coordinate, render WRF inadequate as a long-term solution for wind energy applications.

Most importantly, WRF’s software architecture cannot effectively utilize Graphics Processing Units (GPUs), a requirement both for the increased computational expense of emerging applications, as well as for use of the next generation of high-performance computing (HPC) assets slated for rollout first at the DOE national laboratories, and eventually defining the majority of available scientific computing hardware. Moreover, no current code development efforts pursuing GPU-compatible architectures are addressing the downscaling functionality required of wind energy applications, with the other efforts pursuing either numerical weather prediction (NWP) on global-scale grids, or focusing on microscale simulation without including the critical mesoscale component defining the scales of energy input.

Beyond bridging the scale gap between mesoscale and microscale energy flow, ERF will provide functionality not contained within other mesoscale or microscale simulation codes, specifically the offshore and complex terrain environments that define key challenges for future wind energy expansion. To provide the required functionality, ERF will include advanced physical process models for mesoscale-to-microscale coupling, while seamlessly interfacing with external solvers, such as wave and sea-state models, and the ExaWind microscale wind plant simulation code. Coupling ERF with ExaWind is especially critical, as ExaWind relies on a larger-scale model such as ERF for energy inflow and boundary condition data, in order to provide accurate wind plant optimization guidance in complicated settings and conditions.

In short, ERF will provide a modern, flexible, and efficient GPU-capable software framework to supply critical atmospheric and environmental drivers of energy availability to the microscale wind plant environment, thereby enabling the significant advances in siting, design and operation required to support continued industry expansion into increasingly challenging environments and high penetration scenarios.


ERF User Base
=============

The ERF model is being designed for use by moderately skilled to advanced practitioners of computational fluid dynamics (CFD) codes such as OpenFOAM and WRF, working within the wind energy industry, national laboratories and university research groups, in applications involving numerical simulations of atmospheric flows, and the interactions of those flows with wind turbines, wind plants, and multiple interacting plants in regions of dense development. Beyond its applications, an additional goal of ERF is to serve as a bridge for users and developers of existing CFD and NWP codes to transition away from older legacy codes into a modern software architecture and programming paradigm that efficiently utilizes the next generation of GPU-accelerated HPC hardware, while providing opportunities for expanded applicability to modern wind energy research challenges that require the implementation and integration of new computational capabilities. ERF also targets user-developers who seek to contribute new code back to the ERF code base to improve and expand it, as has occurred within other open-source codes such as WRF and OpenFOAM.


ERF Applications
================

ERF is being designed for wind energy applications requiring simulation of the atmospheric flows that provide the energy available for conversion to power within wind plants. ERF will focus on capturing at high fidelity the atmospheric flow quantities of primary interest for wind turbine and plant design, operation, and optimization, at multiple time and space scales, and in variable operating conditions. The primary applications and phenomena that ERF will address, and how ERF will improve upon the present state of the art, are described below.

1. Resource Characterization
----------------------------
A key application targeted by ERF is wind resource characterization, the assessment of the potential of a given site to produce power over time. Today’s commonly applied resource characterization approaches used in both mesoscale and microscale settings are necessarily of restricted fidelity due to the computational expense of higher-fidelity techniques exceeding industry resources in typical workflows. The computational burden of high fidelity will be significantly mitigated using ERF’s more efficient code architecture, algorithms, and ability to use GPU-accelerated HPC hardware, enabling improved resource characterization at all scales.

At the largest mesoscale applications envisioned for ERF (domains of hundreds to thousands of km on a side; grid spacings of one to ten km), increased throughput will enable several advantages. One is the ability to perform larger ensembles for more robust quantification of uncertainty. A second is the ability to use finer grid spacings over the same geographical footprint in areas requiring increased resolution of landscape or surface features. A third will be the ability to use higher-fidelity physical process models, including the three-dimensional planetary boundary layer (PBL) scheme and machine-learning-based surface layer flux models, each developed within the portfolio of projects supported by the Wind Energy Technologies Office (WETO).
ERF will also support the implementation of more advanced offshore wave and sea-state models, including new capabilities developed under the WETO-supported Wind Forecast Improvement Project 3, which is focused on improving computational approaches for offshore environments. All of these new capabilities will enable much more accurate simulation and characterization of wind-energy-relevant flow phenomena within the atmospheric boundary layer (ABL) than is currently achievable within any existing code base.

At the smallest microscale applications envisioned for ERF (involving domains of a few to tens of km on a side, with grid spacings of a few to tens of m), higher computational throughput will enable improved characterization of the turbulence environment, as well as the impacts of smaller-scale terrain or surface heterogeneities on flow quantities of interest. More accurate treatment of complex terrain effects via immersed boundary methods, as well as the ability to integrate resolved wave and sea-state forcing for offshore applications will significantly enhance microscale resource assessment in these settings.

In addition to providing improved resource characterization in both mesoscale and microscale applications, ERF will enable much more efficient mesoscale-to-microscale coupling via efficient dynamic downscaling, interfacing the microscale turbulence field with the mesoscale forcing that drives it. Such multiscale coupling is especially critical in settings involving complex meteorology and landscape characteristics that supply forcing at scales larger than can be encompassed within even a very large single-domain microscale setup.

2. Forensics
------------
Another application ERF will support is the ability to simulate unique meteorological events of importance, such as those leading to damage or some other outcome for which improved understanding is desired. ERF will enable enhanced forensics abilities via higher-resolution, and higher-fidelity treatment of relevant physical processes impacting the flow, coupled with the flexibility to either ingest larger-scale forcing datasets from forecast models or analysis products, or to set up idealized process-level simulations with controlled forcing.

3. Wind Plant Inflow
--------------------
A primary use case for ERF is the downscaling of mesoscale atmospheric flows to microscale grid spacing, for which all of the relevant scales of motion, including turbulence, are sufficiently resolved to specify turbine-airflow interactions. While this information can be used to estimate resulting power, fatigue loading, and other data, ERF is also being designed to couple directly with the ExaWind microscale wind plant simulation code, within which turbine performance, loading, and controls models of various fidelity can operate directly within the ERF-generated inflow. These coupled ERF-ExaWind simulations will provide unprecedented levels of full-spectrum fidelity, information required to understand and optimize wind plant performance in general, complicated operating conditions and environments.

An additional aspect of wind plant inflow is the impact of entire wind plants on both their own inflow, via blockage effects, as well as on downstream plants via wind plant wakes, gravity waves, or other atmospheric disturbances that large wind plants generate. These issues are of particular importance in areas of dense development, and require a larger simulation footprint than is practical within even a very large single-resolution microscale domain. In addition, difficulties simulating gravity waves using the incompressible solvers that form the basis of many microscale wind plant simulation codes are ameliorated using a fully compressible solver such as ERF will employ. The multiple-resolution capability of ERF, coupled with the incorporation of wind plant wake models applicable at both mesoscale and microscale grid spacings, will provide a flexible framework to better understand wind plant interactions, regional wind power generation, and the regional integration of wind-generated power.

4. Offshore Development
-----------------------
A defining challenge of future wind energy development is the offshore environment, which presents unique operating conditions that require the creation of suitable simulation tools to understand and design for. Among the unique offshore conditions impacting wind energy are swell, wave and sea-state variability that impact the low-level atmospheric flow, hence the available energy. Sea surface temperature and roughness variability can also influence submesoscale motions that are important in offshore environments. The large thermal inertia of water can also support persistent static stability regimes with strong and long-lived impacts on flow and wake propagation, for example, via synoptic-scale advection of air masses with different thermal properties over the water, and due to variability in sea-surface temperature due to the existence of currents or bathymetric influences.

Improved parameterizations to represent these unique features of offshore environments at various scales will be implemented into ERF, along with abilities to explicitly specify wave characteristics in large-eddy simulation (LES) domains containing sufficient mesh resolution to capture the wave shape and impacts of moving wave surfaces on the flow. IBMs may provide a pathway to efficiently implement resolved wave impacts into LES domains as well.

5. Impacts of Complex Terrain
-----------------------------
ERF will be designed to improve the representation of complex terrain and its impacts on the flow, including gravity flows, gravity waves, mountain-valley circulations, and coastal jets, in mesoscale simulations, relative to other mesoscale models. Improved complex terrain capabilities will be incorporated via the use of higher-order numerical stencils for the evaluation of horizontal derivatives over moderately steep terrain, and immersed boundary methods (IBMs) for very steep terrain. These approaches will reduce numerical errors while extending ERF to much steeper slopes than the standard WRF model and similar codes can simulate. IBMs will also permit use of higher resolution, less smoothed terrain than with WRF, which will improve simulation fidelity in complex terrain, and thus improve the local wind accuracy around turbines and plants. IBMs can also stabilize numerical solutions over steep terrain, even if not strictly required, allowing for larger model time steps and therefore accelerating execution. ERF’s use of a vertical coordinate with a fixed height will lead to much more efficient use of IBMs in ERF than in WRF, where the changing heights require new interpolations and projections at every time step.

6. Impacts of Low-Level Jets
----------------------------
An important meteorological feature defining the energy resources in many geographic locations, including the US central great plains and offshore regions, is the low-level jet (LLJ), a narrow ribbon of fast moving air that occurs within the lowest several hundred meters above the surface. While LLJs provide a rich energy resource, LLJs characteristically contain strong sheer, veer and intermittent atmospheric wave and turbulence activity, all of which can increase fatigue loading. Moreover, details of their height and strength, as well as the timing of their onset and dissolution, which impact the integration of power produced, present numerous challenges to development within such regions. More efficient downscaling and higher-fidelity mesoscale and microscale turbulence models will provide enhanced understanding of LLJ impacts on wind power applications.

7. Impact of Clouds and Precipitation
-------------------------------------
Clouds are important modulators of the atmospheric flow, impacting turbulence intensity via shading of the surface, which also influences boundary layer growth and the vertical transfer of momentum, leading to changes in mean wind speed, shear and veer across the turbine rotor swept area. In the offshore environment, radiative cooling from liquid water in the stratocumulus field that often surmounts the marine ABL can drive increased turbulence within the flow below the cloud deck. Clouds also produce various species of precipitation which impact turbine and plant performance, including glaze and rime ice that reduce aerodynamic performance, raindrops that can accelerate leading edge erosion, and graupel and hail which can be particularly damaging. The greater computational expense of cloud microphysics models that can accurately depict these processes, as well as running these schemes at the fine mesh scales required, have hindered understanding and predictive capabilities for the impacts of clouds and precipitation on wind power generation thus far. ERF’s more efficient solution framework will facilitate addressing these gaps. Higher-fidelity representation of cloud impacts on surface shading will also improve ERF’s applicability to solar and hybrid plant operation.

8. Impacts of Thunderstorms and Tropical Disturbances
-----------------------------------------------------
Thunderstorms and tropical disturbances are multiscale phenomena with potentially high impacts on wind plant operation and reliability due to strong winds, high turbulence levels and thus gustiness, large raindrops and hail, and, in offshore environments, associated surface waves that can cause significant damage to turbine platforms. Current atmospheric simulation codes such as WRF provide some capability to investigate storm impacts on wind plants, but simulations are too expensive for routine application. ERF will provide a superior framework for investigating the impacts of strong storms on both onshore and especially offshore environments, where integration with high-fidelity wave and sea-state models will significantly extend predictive capabilities for ABL flow and potentially relevant sea-state characteristics.

9. Regional Integration and Hybrid Plants
-----------------------------------------
The integration of wind-generated electricity at high penetration scenarios is critical to the robust and efficient operation of the grid, and the ability to safely reduce baseload generation from fossil fuel-based plants. Integration relies on accurate predictions of both power generation and demand, each of which is impacted by atmospheric parameters computed by ERF. With a regional mesoscale footprint, ERF will be able to predict the regional distributions of atmospheric quantities impacting both production and demand. ERF can also improve solar forecasting, for many of the reasons it enhances wind predictions (greater throughput, more efficient downscaling, and higher-fidelity physics parameterizations), thereby facilitating integration of wind with another rapidly developing weather-dependent generation source. Moreover, the integration of solar and wind with data analytics from the grid will facilitate exploration of integration with energy storage, transmission, and grid services. As such, ERF will facilitate the development and operation of hybrid plants which balance multiple sources of production, storage and grid services, helping to transform meteorologically-dependent, intermittent sources of energy into robust, dispatchable power generation and delivery.

10. Development of Machine Learning and Artificial Intelligence Methods
-----------------------------------------------------------------------
Machine learning (ML) and artificial intelligence (AI) approaches represent significant opportunities to improve understanding of physical phenomena, enhance the fidelity of numerical simulations, improve the computational efficiency of such simulation codes, and to develop faster running reduced-order models for applications requiring increased throughput in industry workflows. ERF will address development of ML and AI methods by providing the ability to generate simulation-based training datasets via its integration of high-fidelity physical process parameterizations with efficient multiscale simulation. Coupling of ERF with the ExaWind code will provide unique datasets enabling the examination of relationships between mesoscale flow features and turbine-level response that can bypass the most expensive portion of full-physics numerical simulations, the downscaling of flow into and simulation within the microscale domain.

ERF's Role within the Spectrum of Geophysical and Wind Energy Applications
--------------------------------------------------------------------------
These above listed activities and phenomena define the key applications envisioned for the ERF model. However, there are two other flow simulation regimes of relevance to wind-energy that ERF is not intended to address. The first of these is weather forecasting. While ERF could, in principle, be extended to capture larger-scales of meteorological forcing, efforts are already underway elsewhere to create next-generation numerical weather NWP systems, operating at global scales, to capture the largest scales impacting weather system evolution, while also being designed to utilize GPUs for enhanced speed and efficiency. ERF will leverage these concurrent developments by interfacing with a data preprocessor to ingest forecast and analysis fields produced by these new larger-scale models. ERF will focus on the efficient downscaling of those solutions to footprints of relevance to wind energy applications, capturing the associated finer-scale mesoscale and turbulence features, as well as wind plant interactions, along the way.

The second application ERF will not address is very fine-scale microscale simulation and interaction with resolved turbine components. As is the case with NWP, other modeling tools, including the ExaWind code, are being developed to support those activities, with ERF functioning primarily as the provider of inflow and boundary information, downscaled from mesoscale or NWP scales, to those fine-scale application domains, through robust model coupling interfaces. ERF will also be designed to upscale information from finer-scale offline-coupled simulations back into its domain(s) for improved fidelity.


ERF Features and Requirements
=============================

Below is a list of the features that the ERF code must provide, and requirements of the code to support those features, in order to satisfy the above described user base and applications.

1. Excellent Performance on Both CPU- and GPU-Based HPC Platforms
-----------------------------------------------------------------
ERF must be able to run efficiently on both CPU- and GPU-based HPC platforms. This flexibility is required to support enhanced utilization of existing HPC architectures for which significant industry investments have been made, to serve as a vehicle to transition those users to next-generation platforms and programming paradigms, and to support current applications using GPU-accelerated hardware being rolled out today at leadership computing facilities (LCFs). To meet these use cases, ERF must compile and run on a variety of platforms, but also must be configurable for optimal performance on LCFs, including coupling with ExaWind tools on those platforms. Key metrics to assess adequate performance include superior scaling up to tens of thousands of cores on LCF systems, with several levels of mesh refinement, and solution accuracy that meets or exceeds that of legacy codes such as WRF and OpenFoam in similar applications. In addition to LCF machines at DOE labs, integration with new disk storage approaches (e.g., burst buffer at NERSC) should also be explored. Other platforms that would be desirable for ERF to utilize include emerging GPU-based small sized clusters, commodity desktops and laptops with GPU cards, and cloud resources, which are increasingly coming to replace industry-owned HPC resources at many wind energy companies.

2. C++ Base Code with Ability to Incorporate FORTRAN
----------------------------------------------------
For optimal utilization of GPU-based hardware, the ERF source code must be written in C++. However, ERF should also be able to incorporate legacy Fortran source code from other models. While Fortran code incorporated into the C++ code base will not result in optimal performance, it will allow for the rapid expansion of ERF’s capabilities, while providing a pathway to facilitate adoption of ERF by users and developers familiar with Fortran programming and legacy codes. Future development of ERF, including potential community development, can target the rewriting of desired Fortran modules into C++ for enhanced performance.

3. Configurability for Mesoscale, Microscale and Multiscale Simulations
-----------------------------------------------------------------------
ERF must provide flexible configurability that supports mesoscale or microscale simulations, each using a single mesh level, and for multiscale simulations, with as many levels of mesh refinement as required to span a given scale range, from horizontal grid spacings as large as several kilometers to as fine as several meters. It would be preferable to be able to use integer mesh refinement ratios from two to approximately ten or so in order to support the ability to avoid certain grid spacings depending on the application. The code must also support vertical refinement, however due to the nature of the vertical coordinate and required applications, vertical refinement should enable arbitrary height levels to be specified on all meshes.

Atmospheric downscaling applications will focus primarily on increased resolution over particular geographical areas, rather than tracking flow features that move in time, and therefore static refinement over rectangular volumes is sufficient for the near term (although eventually adaptive refinement will be useful to track features such as storms or turbine wakes). One-way coupling at the refinement interfaces, for which the fine mesh only receives information from, but does not transmit information back to the parent mesh, is acceptable at the beginning. Two-way nesting, for which the fine mesh provides information back to the parent mesh, is eventually needed to improve simulation also on those coarser domains.

The downscaling information exchange at mesh interfaces can follow procedures used in other codes, such as WRF, in which the parent mesh executes one time step, after which variables at the beginning and end of each parent time step are interpolated linearly in time to the refined-mesh time step to support integration of the finer mesh solution. For two-way coupling, the fine-mesh solution can be averaged to the coarse mesh after advancing to the coarse-mesh time step, where it can either replace the solution on that mesh, or provide a target for relaxation of the coarser-mesh solution.

4. Initial and Boundary Condition Preprocessing for Real-Data Simulations
-------------------------------------------------------------------------
To facilitate adoption of ERF for its primary intended workflows, those involving mesoscale energy simulation and its downscaling, ERF must support straightforward methods to prepare initial conditions and forcing for its lateral and surface boundaries. This process requires interpolation in time and space from the locations of state variables within the native datasets, to the locations of ERF’s mesh (or meshes in downscaling setups). The input data must also be transferred onto the ERF model’s map projection. Adopting some of the functionality of the WRF Preprocessing System (WPS) would provide an excellent pathway, as WPS has the ability to read in data from multiple larger-scale forecast and analysis products, and to project those data onto a number of standard map projections. As the WRF model supports global simulation, as well as configurability over arbitrary locations, including equatorial and polar latitudes, it contains several map projections. With ERF’s focus on mid-latitude locations, a Lambert Conformal projection would be ideal at the initial development phase. ERF can also initially focus on global analysis and forecast models such as NCEP’s Global Forecast System (GFS), as well as for North American applications, NOAA's High Resolution Rapid Refresh (HRRR) model.

5. Initial and Boundary Condition Support for Idealized Simulations
-------------------------------------------------------------------
To serve a variety of process-level applications and interaction studies, ERF must provide support for idealized setups, easily configurable into arbitrary (rectangular) domains with simplified, user-defined specification of input (vertical distributions of state variables), forcing (horizontal pressure gradient, geostrophic wind, and wind profile assimilation methods), and support for idealized boundary conditions, including open and periodic lateral boundaries, and appropriate top and surface boundary conditions such as free slip, no flux and no normal flow. Radiative fluxes at the domain top should also be easy to specify for use with low model tops, which is commonly done to reduce cost in high-resolution simulations under applicable atmospheric conditions.

6. Compressible Nonhydrostatic Dynamic Core
-------------------------------------------
To simulate atmospheric forcing arising from mesoscale processes, a fully compressible, non-hydrostatic equation set appropriate for dynamics involving vertical density variability and well-resolved vertical motions must be used. While the WRF model uses a pressure-based vertical coordinate system for which the height above the surface of model grid points change in time, ERF should adopt a system for which the heights remain constant over time, simplifying applicability to wind energy use cases, including the coupling with other codes such as ExaWind. The equation set and vertical coordinate used by the COSMO (Consortium for Small Scale Modeling) model would serve as an excellent template for ERF.

7. Second-Order Finite Difference Spatial Discretization
--------------------------------------------------------
The ERF equation set will require a discretization strategy and numerical solution procedures amenable to optimization for different mesh spacings and applications. For ease of implementation and familiarity with users of other code bases, as well as ability to incorporate modules from other codes such as WRF, a finite-difference spatial discretization strategy with second-order accuracy should be employed for ERF. Options for higher-order spatial differences can be included as well, however those methods may not scale as well on GPU-based hardware.

8. Fully Explicit and Mixed Implicit-Explicit Time Discretizations
------------------------------------------------------------------
For time advancement, a fully explicit method is required for general applications, however the flexibility to implement options for implicit treatment in the vertical, or all three directions, should be included. Different meshes must be able to use different time steps (selectable by the user based upon the spatial refinement ratio). Substepping in time for acoustic mode propagation should be included for the advancement of the pressure or density field in mesoscale domains where a time step sufficient to resolve those modes would be untenable.

9. Application-Relevant Physical Process Parameterizations
----------------------------------------------------------
To support mesoscale, microscale and downscaling simulations, ERF must contain physical process parameterizations appropriate to all relevant scales and processes, including:

*    Monin-Obukhov Similarity Theory (MOST) surface stress boundary condition
*    Wave damping treatments for the upper domain
*    Subgrid turbulence closure for large-eddy simulation
*    Subgrid turbulence closure for mesoscale simulation
*    Surface energy budget parameterization for surface energy and momentum fluxes
*    Vegetation canopy model for improved flow over tall vegetation
*    Shortwave and longwave radiative transfer to capture solar/diurnal forcing and cloud-induced radiation impacts
*    Cloud microphysics to capture cloud-radiative forcing and precipitation
*    Offshore wave, sea-state and ocean current models appropriate for various spatial scales
*    Wind plant wake models for microscale and mesoscale applications

Following the WRF model, the physical process parameterizations should be callable on user-adjustable time steps (for faster simulation execution), and if a Runge-Kutta time advancement scheme is chosen, computed on the first predictor step of the sequence.

Hindcasting and forensics applications will require methods to incorporate forcing data from either analysis datasets (“analysis” nudging) or from observations (“observation” nudging), using Newtonian relaxation, or spectral methods.

Many of ERF’s required parameterizations and other capabilities can be taken from existing models such as WRF and implemented into the ERF source code. Another pathway is to interface ERF with the community physics package (CPP), a repository of common physics modules developed for integration with various community models such as WRF, MPAS and FV3. As described previously, providing an ability to compile modules written in Fortran into the ERF executable will facilitate the incorporation of existing codes from other models, and from the CPP, as well as facilitating the transition of users of older legacy codes and modules to ERF.

10. Robust and Efficient Interfacing of ERF with Other Code Bases
-----------------------------------------------------------------
ERF must robustly and efficiently interface with other codes, such as ExaWind, modules within the CPP, wave, sea-state, and wind turbine aerodynamics and load models to support ERF’s applications. The interfacing must allow for the efficient exchange of required variable values or forcing terms at the domain boundaries or overlap regions of the coupled codes.

11. Efficient and Flexible Selection of Output and Internal Diagnostics
-----------------------------------------------------------------------
ERF should incorporate efficient parallel input/output (I/O) strategies to support runs over large processor counts, and frequent output, including options for asynchronous I/O to overlap I/O and computation. For variables on staggered grids, projection to cell-centered locations is desirable. For variables decomposed into perturbation and base states, reconstruction of the entire physical field before outputting is desirable to reduce file sizes.

To further reduce output file sizes and time spent writing output fields, ERF must provide abilities to select output for specific variables over domain subvolumes and slices of arbitrary size and orientation, and to compute diagnostic quantities such as spectra and Reynolds stresses, either over time or arbitrary spatial directions and volumes, during code execution.

ERF should adopt NetCDF climate and forecast (CF)-compliant metadata or a sufficiently close proximity thereof for output files to enable the use of common plotting and analysis software like xarray in Python, Ferret, NCL, etc.

12. Flexibility and Extensibility
---------------------------------
For the ERF code to remain viable on evolving computational hardware while maintaining abilities to address evolving wind energy applications, the code architecture must provide high levels of flexibility and extensibility. These attributes must encompass data and memory flow and management, as required to utilize future HPC hardware, ability to incorporate more efficient and accurate numerical solution methods and physical process models, and ability to interface with other codes, including not only those related to geophysical process and wind energy design, but other related applications such as grid, transmission, storage and others.

13. Detailed Documentation and Test Cases
-----------------------------------------
To attract new users and encourage development by the user community, ERF must contain detailed documentation and test cases describing configurability and pathways for extensibility. The documentation should take the form of a users’ guide for a high-level understanding of how to use the code for various applications, a detailed technical volume describing the equations, discretization strategy, solution methods, and data management, and extensive comments within the source code that describe the function of specific segments. The WRF model provides an excellent example of a documentation and code structure that facilitates adoption by new users, while encouraging community contributed extensions to its capabilities.

14. Balance between Performance, Robustness and Usability
---------------------------------------------------------
ERF development must strike an appropriate balance between achievement of optimal performance on priority HPC platforms
and other considerations, including robustness of the code to faithfully compile and run on different platforms,
tractability for use by a diverse set of practitioners in a diverse set of applications, and maintainability
to ensure a manageable workload for continued utility.

ERF Code Design
===============

Development of the ERF source code requires a detailed strategy for the implementation of specific components and capabilities, a vision for overall code topology to provide the required flexibility and extensibility, protocols for code development consistency and documentation, and robust testing to establish code performance and ensure continuity of performance and functionality under continuing development.

1. Development of Appropriate Code Development Protocols
--------------------------------------------------------
To ensure a robust record of code development history, ERF will be developed and managed within a GIT commit structure, with protocols for description and verification of implementations following the ExaWind project.

Documentation will be required on three levels:

*    Unit-level information within the source code to aid future users and developers
*    Higher-level design and implementation information in accompanying documentation on a readthedocs.org portal
*    A users’ guide detailing code functionality and describing accompanying test cases to demonstrate that functionality to new users; to learn by going through examples.
*    Shared analysis scripts to ensure consistency of results across the team

Development of the code will follow modern software project paradigms, including

*    Configurability of software packages

  * Different combinations of components constituting different applications

    - E.g., single- versus multi-resolution, real data versus idealized, microscale only versus mesoscale or multiscale (requiring physics packages), stand-alone versus coupled with other codes, …

  * Assessment of compile-time versus run-time configurability to guide optimization

* Encapsulation

  * Related functionality and data that can be grouped as a single class are organized into encapsulated code units that can have multiple alternative implementations.
  *    Code libraries to provide services such as discretization, data management, and orchestration of data movement for parallelization as well as I/O.
  *    Physics solvers are largely ignorant of the details managed by the framework (dynamic core and I/O), and operate in “plug-and-play” mode within the framework to enable streamlining for specific applications, or to perform ensembles for which configuration, dynamics and physics options can be modified at compile or run time

*    Separation of concerns so that, e.g., infrastructure and physics solvers development do not intrude into each other’s space

*    Hierarchy of parallelism within and across code units

2. Development of Robust Testing Strategies
-------------------------------------------
Robust testing at multiple levels is required for verification, validation, and characterization of code performance and simulation accuracy.

*    At the smallest unit level, test cases for unit-level commits will be provided by the originators of those code units and included in a master suite of unit tests which must exhibit acceptable performance during future code development.
*    At intermediate aggregation levels for which functional units group together to provide a moderately complex capability, those sub-model aggregations must exhibit acceptable performance against suitable test cases, under future code development.
*    At the highest aggregation levels for which functional groups combine together to provide a complex capability, those whole-model aggregations must exhibit acceptable performance against suitable aggregate-level test cases under future development.

Acceptable performance should consist of bit-for-bit agreement for the addition of unrelated code components, or components that operate on data management and movements, but are not expected to alter values.

Higher-level tests should consist of performance against standards, such as analytical or manufactured solutions, as well as convergence under spatial or temporal refinement, depending on the sub-model component. The highest-level or “whole-model” tests should consist of evaluation against data from field campaigns or previous appropriate whole-model results, such as obtained from WRF or OpenFoam.

These standards should be articulated at a sufficient level of detail to guide community contributors of requirements for code additions or improvements.


Implementation Details for Specific Code Components
===================================================

This section provides detailed information on core elements of the ERF code that enable its required functionality. Details of both the formulation of methods and their implementation and usability within the code will be added to the existing higher-level descriptions as the project progresses.

1. Governing Equation Set
-------------------------
Details of ERF’s governing equation set are currently still being formulated, however it will closely follow the formulation used by the COSMO model, with the following attributes:

*    Fully compressible non-hydrostatic formulation
*    Terrain-following computational mesh
*    Fixed-height vertical coordinate
*    Prognostic variables (three dimensional winds, a conserved temperature variable (potential or moist potential temperature, and pressure) cast in perturbation form, relative to hydrostatic base state.

This formulation does not conserve mass, but errors are negligible over intended timescales of simulations (hours to days).

2. Spatial Discretization
-------------------------
Spatial discretization of the ERF governing equation set will follow a second-order finite difference strategy, on an Arakawa C-grid. Horizontal grid spacing will be equal in each direction and constant within each mesh refinement level with respect to adjustments necessary for mapping factors. The vertical coordinate will support user specification of heights on each refinement level. Higher-order difference formulations, including upwinding schemes, will also be provided if needed and resources allow. A Lambert conformal map projection will be used.

Immersed boundary methods will also be implemented for a variety of applications, including reduced numerical errors over complex terrain, extension to very steep terrain up to and including vertical cliff walls, incorporation of embedded turbine and structure geometries, and potentially as a method to incorporate resolved wave impacts on atmospheric flows. Various methods developed within the WRF model can be straightforwardly extracted from that code and implemented into ERF, within which the ability to maintain constant heights of the model vertical coordinate will greatly enhance code performance when in use.

3. Temporal Discretization
--------------------------
The base temporal discretization strategy for ERF will be an explicit third-order Runge-Kutta predictor-corrector method, similar to that used in the WRF code and others. The scheme will permit users to specify time steps for each mesh refinement level. For improved performance, smaller time steps for acoustic modes within each Runge-Kutta time step will be enabled for applications with grid spacings large enough to improve performance. Implicit methods will also be investigated for further performance enhancements. Implicit methods for the vertical direction only will be explored for cases in which the vertical grid spacing is small relative to the horizontal, requiring a prohibitively small model time step, as well as to stabilize the code with a larger time step. Advantages of also using implicit methods in the horizontal direction relative to the increased memory requirements of such formulations will also be explored.

4. Numerical Solution Methods
-----------------------------
ERF should be formulated to support integration of a variety of numerical solution methods for different applications, and to incorporate future techniques that can provide superior performance or other advantages. The team will prioritize use of solvers available within the PETSc libraries that have already been incorporated into the mesh refinement framework that is also being used for ERF.

5. Mesh Refinement
------------------
ERF will utilize the AMReX adaptive mesh refinement framework for its computational mesh and refinement requirements. AMReX provides a flexible capability that can support all of ERF’s required mesh needs utilizing advanced data structures and memory management for robust and efficient data transfer and load balancing. Moreover, AMReX contains built-in abstractions to efficiently interface with a variety of CPU- and GPU-based HPC hardware. The continuing support of AMReX by the Exascale Computing Program makes AMReX an ideal choice for ERF.

While ERF will most likely never utilize all of the adaptive mesh refinement capabilities available within AMReX, requiring principally only static regions of refinement, there is no penalty for not utilizing those additional capabilities, and those capabilities might prove valuable to future applications.

6. Boundary Conditions
----------------------
List of top, bottom and lateral boundary condition implementations provided here.

The highest priority needed for initial code testing and to simplify additional development includes capability to support ABL flow, including free-slip boundaries at the top and bottom, with no flow normal to the domain top and surface, and no fluxes of energy or constituents, periodic lateral boundaries in each direction, and wave damping at the model top.

Once the basic code is functional, the ability to read time dependent boundary conditions from input files using modified WPS code will be implemented to enable real-world simulations.

7. Physical Process Parameterizations
-------------------------------------
List of physical process models described here. Priority development includes:

*    MOST surface stress boundary condition for surface momentum fluxes
*    Subgrid turbulence flux model for large-eddy simulation
*    Subgrid turbulence flux model (PBL scheme) for mesoscale simulation
*    Hooks to the CPP for radiation, cloud, surface, and boundary layer parameterizations

8. Model Coupling
-----------------
Explore pathways to interface both ERF and ExaWind within the same AMReX platform.

9. Input/Output
---------------
Short-term goal: Utilize native AMReX output which can be read by ParaView and VisIt. Add capability to write to NetCDF in a mostly-CF-compliant format, enabling integration with Xarray in Python, and potentially other plotting utilities.

Long-term goal: Implement capability to configure I/O at run time based, for example, on a YAML file that is separate from the configuration file used to run the code. The separate I/O file would support:

*    Adding and removing specific variables from the output
* Changing output frequency for different variables
*    Use of multiple output files with different sets of variables and output frequencies
