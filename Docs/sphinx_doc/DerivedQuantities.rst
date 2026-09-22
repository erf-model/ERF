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

Vorticity stretching
--------------------

ERF provides ``vort_stretching`` as an optional plotfile variable. It is the
stretching term in the evolution equation for the vertical component of
vorticity,

.. math::

   S = \zeta \, \frac{\partial w}{\partial z}, \qquad
   \zeta = \frac{\partial v}{\partial x} - \frac{\partial u}{\partial y},

and has units of :math:`\mathrm{s}^{-2}`. Writing the vertical vorticity
equation for an inviscid flow as

.. math::

   \frac{D\zeta}{Dt} = \underbrace{\zeta \frac{\partial w}{\partial z}}_{\rm stretching}
     + \underbrace{\omega_x \frac{\partial w}{\partial x}
                 + \omega_y \frac{\partial w}{\partial y}}_{\rm tilting}
     + \ \ldots

makes clear that ``vort_stretching`` is the first term only. The tilting term,
which converts horizontal vorticity into vertical vorticity, is not included in
this diagnostic and is not output separately, nor are the baroclinic, Coriolis,
or diffusive contributions hidden in the ellipsis.

Because the stretching term is proportional to the vertical vorticity already
present, it amplifies existing rotation rather than creating it: it is positive
where a rotating column is being stretched vertically and negative where it is
being compressed. Values are largest in regions of strong rotation coincident
with strong vertical convergence, which is why the field is useful in severe
convective storm and tornadogenesis studies.

The quantity is evaluated cell-by-cell from the cell-centered velocity using
centered differences. As for ``vorticity_x`` and ``vorticity_y``, the vertical
derivative is formed from the physical heights of the cell centers, so it is
correct on a vertically stretched mesh; on a terrain-fitted mesh it neglects the
horizontal metric terms of the mapping, exactly as the vorticity components
themselves do. Unlike ``helicity``, ``precipitable``, and ``mucape``, it is a
purely local quantity and places no restriction on how the grid is decomposed in
the vertical.
