ERF 
----
*An atmospheric modeling code*

`ERF` is built on `PeleC` (https://pelec.readthedocs.io/en/latest/), a compressible hydrodynamics code for reacting flows built on the AMReX solution framework (https://amrex-codes.github.io/amrex/). 

Getting Started 
~~~~~~~~~~~~~~~

* To compile and run the `ERF` suite of codes, one needs a C++ compiler that supports the C++14 standard.  A hierarchical strategy for parallelism is supported, based MPI + OpenMP, or MPI + CUDA.  The codes work with all major MPI and OpenMP implementations.  The codes should build and run with no modifications to the `make` system if using a Linux system with the GNU compilers, version 4.9.4 and above.

To build `ERF` and run a sample problem:

1. One can have ERF use the default submodule for AMReX in its own repo by simply performing: ::

    git clone --recursive git@github.com:ERF-model.git
    cd ERF/ExecCpp/RegTests/Sod
    make
    ./ERF3d.xxx,yyy.ex inputs_ex

Alternatively, one can set environment variables to use AMReX repo from external locations: ::

1. Set the environment variable, AMREX_HOME, and clone a copy of `AMReX` there: ::

    export AMREX_HOME=<location for AMReX>    
    git clone git@github.com:AMReX-Codes/amrex.git ${AMREX_HOME}

2. Set the environment variable, ERF_HOME, and clone a copy of `ERF` there. You should be placed in the `development` branch: ::

    export ERF_HOME=<location for ERF code>
    git clone git@github.com:ERF-model.git ${ERF_HOME}

3. Move to an example build folder, build an executable, run a test case: ::

    cd ${ERF_HOME}/ExecCpp/RegTests/Sod
    make
    ./ERF3d.xxx,yyy.ex inputs_ex

* Notes

   A. In the exec line above, xxx.yyy is a tag identifying your compiler and various build options, and will vary across pltaform.  (Note that GNU compilers must be at least 4.8.4, and MPI should be at least version 3).
   B. In addition to informative output to the terminal, periodic plotfiles are written in the run folder.  These may be viewed with CCSE's Amrvis (<https://ccse.lbl.gov/Downloads/downloadAmrvis.html>) or Vis-It (<http://vis.lbl.gov/NERSC/Software/visit/>):

      1. In Vis-It, direct the File->Open dialogue to select the file named "Header" that is inside each plotfile folder..
      2. With Amrvis, "amrvis3d plt00030", for example.


Origin of ERF 
~~~~~~~~~~~~~

`ERF` was created as a renamed, stripped down version of `PeleC`
(<https://github.com/AMReX-combustion/PeleC>),
incorporates a modified RK3 compressible hydro integrator adapted from 
the FHDeX code base (<https://github.com/AMReX-FHD/FHDeX>), 
and is built on the `AMReX` library (<https://github.com/AMReX-codes/AMReX>).

Development model
~~~~~~~~~~~~~~~~~

To add a new feature to ERF, the procedure is:

1. Create a branch for the new feature (locally): ::

    git checkout -b AmazingNewFeature

2. Develop the feature, merging changes often from the development branch into your AmazingNewFeature branch: ::
   
    git commit -m "Developed AmazingNewFeature"
    git checkout development
    git pull                     [fix any identified conflicts between local and remote branches of "development"]
    git checkout AmazingNewFeature
    git merge development        [fix any identified conflicts between "development" and "AmazingNewFeature"]

3. Push feature branch to ERF repository: ::

    git push -u origin AmazingNewFeature [Note: -u option required only for the first push of new branch]

4. Submit a merge request through git@github.com:ERF-model.git, and make sure you are requesting a merge against the development branch

5. Check the CI status on Github and make sure the tests passed for merge request

.. note::

   Github CI uses the CMake build system and CTest to test the core source files of ERF. If you are adding source files, you will need to add them to the list of source files in the ``CMake`` directory for the tests to pass. Make sure to add them to the GNU make makefiles as well.


Test Status - ***UPDATE THIS***
~~~~~~~~~~~

Nightly test results for ERF against multiple compilers and machines can be seen on its CDash page `here <https://my.cdash.org/index.php?project=ERF>`_. Static analysis results for ERF can be seen in the notes of the newest GCC compiler on CDash. ERF is also tested using the Clang address sanitizer to detect memory leaks.

Test results for the GNU Make implementation of ERF can be seen `here <https://amrex-combustion.github.io/ERFRegressionTestResults>`_.


Documentation - ***UPDATE THIS***
~~~~~~~~~~~~~

The full documentation for ERF exists in the Docs directory; at present this is maintained inline using Doxygen
and Sphinx  `Sphinx <http://www.sphinx-doc.org>`_. With 
Sphinx, documentation is written in *Restructured Text*. reST is a markup language
similar to Markdown, but with somewhat greater capabilities (and idiosyncrasies). There
are several `primers <http://thomas-cokelaer.info/tutorials/sphinx/rest_syntax.html>`_
available to get started. One gotcha is that indentation matters.
To build the documentation, run Doxygen in the Docs directory then build the sphinx ::

    doxygen Doxyfile
    cd sphinx_doc
    make html


Acknowledgment - ***UPDATE THIS***
~~~~~~~~~~~~~~

This research was supported by the Exascale Computing Project (ECP), Project
Number: 17-SC-20-SC, a collaborative effort of two DOE organizations -- the
Office of Science and the National Nuclear Security Administration --
responsible for the planning and preparation of a capable exascale ecosystem --
including software, applications, hardware, advanced system engineering, and
early testbed platforms -- to support the nation's exascale computing
imperative.
