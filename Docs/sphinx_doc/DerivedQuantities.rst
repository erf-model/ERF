 .. role:: cpp(code)
    :language: c++

.. _sec:derived:

Derived Variables
=================

ERF has the ability to create new temporary variables derived from the state variables.

Access to the derived variable is through one of two amrex:AmrLevel functions
(which are inherited by ERF)

::

        /**
        * \brief Returns a MultiFab containing the derived data for this level.
        * The user is responsible for deleting this pointer when done
        * with it.  If ngrow>0 the MultiFab is built on the appropriately
        * grown BoxArray.
        */
        virtual std::unique_ptr<MultiFab> derive (const std::string& name,
                              Real               time,
                              int                ngrow);
        /**
        * \brief This version of derive() fills the dcomp'th component of mf
        * with the derived quantity.
        */
        virtual void derive (const std::string& name,
                             Real               time,
                             MultiFab&          mf,
                             int                dcomp);

Derived quantities as well as state variables can be output in the plotfiles.
