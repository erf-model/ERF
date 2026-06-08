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

MUCAPE
------

ERF provides ``mucape`` as an optional plotfile variable. The current ERF
implementation defines MUCAPE column-by-column using the following choices:

- Search the lowest 300 hPa of each column for the most unstable parcel.
- Define the most unstable parcel as the parcel in that search layer with the
  largest resulting CAPE.
- Dry-lift the parcel to the LCL, then continue with a pseudoadiabatic moist ascent.
- Compute buoyancy from virtual temperature.
- Integrate only the positive buoyancy contribution, so the result has units of J/kg.
- Store the single column value at every vertical level in that column, following
  the same column-diagnostic pattern used by quantities such as ``precipitable``
  and ``helicity``.

This definition is ERF-specific and may differ from other packages that choose a
different parcel search layer, parcel-selection rule, ascent model, or buoyancy correction.
