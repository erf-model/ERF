 .. role:: cpp(code)
    :language: c++

.. _AgenticWorkFlow:

ERF Agentic Workflow
====================

NOTE: This section is a work in progress and will be continually updated as appropriate.

`AMReX-Agent <https://github.com/AMReX-Codes/amrex-agent>`_ is an agentic AI workflow for amrex-based
simulation codes. Since ERF is based upon the AMReX library, one may utilize the AMReX-Agent software to
automate ERF simulations and postprocess results. Here, we provide a brief overview of how to utilize
the AMReX-Agent code but refer the interested readers to the `ERF agent demo <https://github.com/AMReX-Codes/amrex-agent/tree/development/demo/erf>`_ for more details.


Basic Set up
------------

1. Clone the ``AMReX-Agent`` software

.. code-block:: bash

   git clone --recursive https://github.com/AMReX-Codes/amrex-agent.git
   cd amrex-agent
   git checkout 6e77ca3

2. Set up your API key

.. code-block:: bash

   export CBORG_API_KEY=<your_key_name>

3. Clone ``ERF`` and export the path (recommended to add export to .bashrc)

.. code-block:: bash

    git clone https://github.com/erf-model/ERF.git --recursive
    export ERF_REPO_PATH=<path_to_ERF>

4. Set up the environment

.. code-block:: bash

   conda env create -f environment.yaml
   conda activate amrex-agent-dev

5. Build ERF schemas and indices

.. code-block:: bash

   bash demo/setup_demo_database.sh --code erf --force-rebuild

.. NOTE: Startup preflight does not block on ERF pin mismatch when local schema
   and FAISS artifacts are valid for the current ERF commit. The pinned commit in
   ``.dependencies.json`` remains the recommendation for reproducible shared runs.

6. Prompt the ``AMReX-Agent``. Example here requests a local simulation that runs a 2D squall line with 4 ranks and plots the cloud water

.. code-block:: bash

   python amrex_agent.py --run_ntasks 4 --indexing-strategy simple --inputs-file-strategy llm_compare \
   --json --prompt "Run a 2D squall line simulation with Kessler microphysics, open x boundaries, \
   and HO outflow aloft. Run the simulation for 10000 steps and visualize the cloud water."

