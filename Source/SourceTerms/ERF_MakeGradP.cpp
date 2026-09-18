#include <AMReX_MultiFab.H>
#include <AMReX_ArrayLim.H>
#include "AMReX_BCRec.H"

#include "ERF.H"
#include "ERF_SrcHeaders.H"
#include "ERF_DataStruct.H"
#include "ERF_Utils.H"
#include "ERF_EB.H"
#include <ERF_EBSlopes.H>

using namespace amrex;

/**
 * Vertical derivative of a cell-centered quantity, taken with respect to physical
 * height.  Centered in the interior of the domain, one-sided at the bottom and top.
 *
 * @param[in] i x-index
 * @param[in] j y-index
 * @param[in] k z-index
 * @param[in] klo lowest k index of the domain
 * @param[in] khi highest k index of the domain
 * @param[in] q_arr cell-centered quantity to differentiate
 * @param[in] z_cc_arr physical height at cell centers
 * @return dq/dz at (i,j,k)
 */
AMREX_GPU_DEVICE
AMREX_FORCE_INLINE
Real
dqdz_cc (int i, int j, int k,
         int klo, int khi,
         const Array4<const Real>& q_arr,
         const Array4<const Real>& z_cc_arr)
{
    int km = (k == klo) ? k : k-1;
    int kp = (k == khi) ? k : k+1;
    if (kp == km) { return zero; }
    return (q_arr(i,j,kp) - q_arr(i,j,km)) / (z_cc_arr(i,j,kp) - z_cc_arr(i,j,km));
}

/**
 * Function for computing the pressure gradient
 *
 * @param[in]  level     level of resolution
 * @param[in]  geom      geometry container at this level
 * @param[in]  S_data    current solution
 * @param[in]  base_state base state (r0, p0, pi0, th0, qv0)
 * @param[in]  qt        total water mixing ratio (zero if there is no moisture)
 * @param[in]  z_phys_nd z on nodes
 * @param[in]  z_phys_cc z on cell centers
 * @param[in]  d_bcrec_ptr Boundary Condition Record
 * @param[in]  ebfact    EB factory container at this level
 * @param[out] gradp     pressure gradient
 */

void make_gradp_pert (int level,
                      const SolverChoice& solverChoice,
                      const Geometry& geom,
                      Vector<MultiFab>& S_data,
                      const MultiFab& base_state,
                      const MultiFab& qt,
                      const MultiFab& z_phys_nd,
                      const MultiFab& z_phys_cc,
                      Vector<std::unique_ptr<MultiFab>>& mapfac,
                      const eb_& ebfact,
                      Vector<MultiFab>& gradp)
{
    const bool l_use_moisture  = (solverChoice.moisture_type != MoistureType::None);
    const bool l_eb_terrain    = (solverChoice.terrain_type == TerrainType::EB);

    const bool l_use_pert_pres = (solverChoice.use_pert_pres_gradient);
    //
    // Note that we only recompute gradp if compressible;
    //      if anelastic then we have computed gradp in the projection
    //      and we can reuse it, no need to recompute it
    //
    if (solverChoice.anelastic[level] == 0)
    {
        if (solverChoice.gradp_type == 1) {
            AMREX_ASSERT_WITH_MESSAGE(solverChoice.terrain_type != TerrainType::EB,
                "gradp_type==1 not implemented for EB");
        }

        // gradp_type 2 and 3 carry p' to the height of the lateral face along the
        //    hydrostatic relation dp'/dz = -g rho', so they need rho'.
        //
        // NOTE: this must be an always-assert rather than a debug-only one, because with
        //       EB we grow the boxes by more ghost cells than qt carries.
        const bool l_need_rhopert = (solverChoice.gradp_type >= 2);
        if (l_need_rhopert) {
            AMREX_ALWAYS_ASSERT_WITH_MESSAGE(solverChoice.terrain_type != TerrainType::EB,
                "gradp_type 2 and 3 are not implemented for EB");
        }

        const int ngrow = (l_eb_terrain) ? 3 : 1;
        MultiFab p(S_data[Vars::cons].boxArray(), S_data[Vars::cons].DistributionMap(), 1, ngrow);

        MultiFab rhopert;
        if (l_need_rhopert) {
            rhopert.define(S_data[Vars::cons].boxArray(), S_data[Vars::cons].DistributionMap(), 1, ngrow);
        }

        // *****************************************************************************
        // Compute pressure
        // *****************************************************************************
        for ( MFIter mfi(S_data[Vars::cons]); mfi.isValid(); ++mfi)
        {
            Box gbx = mfi.tilebox();
            gbx.grow(IntVect(ngrow,ngrow,ngrow));

            if (gbx.smallEnd(2) < 0) gbx.setSmall(2,0);
            const Array4<const Real>& cell_data = S_data[Vars::cons].array(mfi);
            const Array4<      Real>& pp_arr = p.array(mfi);
            ParallelFor(gbx, [=] AMREX_GPU_DEVICE (int i, int j, int k) noexcept
            {
                Real qv_for_p = (l_use_moisture) ? cell_data(i,j,k,RhoQ1_comp)/cell_data(i,j,k,Rho_comp) : zero;
                pp_arr(i,j,k) = getPgivenRTh(cell_data(i,j,k,RhoTheta_comp),qv_for_p);
            });
        }

        // If we want to use the full pressure in the lateral gradients, call compute_gradp_xy here
        if (solverChoice.gradp_type == 0 && !l_use_pert_pres) {
            compute_gradp_xy(p,geom,z_phys_cc,mapfac,ebfact,gradp,solverChoice);
        }

        // *****************************************************************************
        // Compute perturbational pressure -- and, if needed, perturbational density
        //
        // NOTE: rho' here is the same quantity that buoyancy_rhopert forms, i.e. the
        //       total (moist) density minus the total base-state density.  Keeping the
        //       two definitions identical is what makes the hydrostatic reconstruction
        //       below consistent with the buoyancy term.
        // *****************************************************************************
        for ( MFIter mfi(S_data[Vars::cons]); mfi.isValid(); ++mfi)
        {
            Box gbx = mfi.tilebox();
            gbx.grow(IntVect(ngrow,ngrow,ngrow));

            if (gbx.smallEnd(2) < 0) gbx.setSmall(2,0);
            const Array4<const Real>& base_arr = base_state.const_array(mfi);
            const Array4<      Real>& pp_arr = p.array(mfi);
            ParallelFor(gbx, [=] AMREX_GPU_DEVICE (int i, int j, int k) noexcept
            {
                pp_arr(i,j,k) -= base_arr(i,j,k,BaseState::p0_comp);
            });

            if (l_need_rhopert) {
                const Array4<const Real>& cell_data = S_data[Vars::cons].const_array(mfi);
                const Array4<const Real>& qt_arr    = qt.const_array(mfi);
                const Array4<      Real>& rp_arr    = rhopert.array(mfi);
                ParallelFor(gbx, [=] AMREX_GPU_DEVICE (int i, int j, int k) noexcept
                {
                    Real qt_loc = (l_use_moisture) ? qt_arr(i,j,k) : zero;
                    rp_arr(i,j,k) = cell_data(i,j,k,Rho_comp) * (one + qt_loc)
                                  - base_arr(i,j,k,BaseState::r0_comp)
                                  * (one + base_arr(i,j,k,BaseState::qv0_comp));
                });
            }
        }

        // If we want to use the perturbational pressure in the lateral gradients, call compute_gradp_xy here
        if (solverChoice.gradp_type == 0 && l_use_pert_pres) {
            compute_gradp_xy(p,geom,z_phys_cc,mapfac,ebfact,gradp,solverChoice);
        }

        // gradp_type 2 and 3 replace the lateral gradients only; gpz is computed below
        //    exactly as it is for gradp_type 0
        if (l_need_rhopert) {
            compute_gradp_hse(p,rhopert,geom,z_phys_cc,mapfac,gradp,solverChoice);
        }

        if (solverChoice.gradp_type == 1) {
            compute_gradp_interpz(p,geom,z_phys_nd,z_phys_cc,mapfac,gradp,solverChoice);
        } else {
            compute_gradp_z(p,geom,z_phys_nd,ebfact,gradp,solverChoice);
        }

    } // not anelastic
}

/**
 * @brief Compute the full pressure gradient.
 * @param[in] p Pressure field.
 * @param[in] geom Geometry container.
 * @param[in] z_phys_nd Physical height on nodes.
 * @param[in] z_phys_cc Physical height on cell centers.
 * @param[in] mapfac Map factors.
 * @param[in] ebfact EB factory.
 * @param[out] gradp Pressure gradient components.
 * @param[in] solverChoice Solver options.
 */
void
compute_gradp (const MultiFab& p,
               const Geometry& geom,
               const MultiFab& z_phys_nd,
               const MultiFab& z_phys_cc,
               Vector<std::unique_ptr<MultiFab>>& mapfac,
               const eb_& ebfact,
               Vector<MultiFab>& gradp,
               const SolverChoice& solverChoice)
{
    compute_gradp_xy(p,geom,z_phys_cc,mapfac,ebfact,gradp,solverChoice);
    compute_gradp_z(p,geom,z_phys_nd,ebfact,gradp,solverChoice);
}

/**
 * @brief Compute the horizontal components of the pressure gradient.
 * @param[in] p Pressure field.
 * @param[in] geom Geometry container.
 * @param[in] z_phys_cc Physical height on cell centers.
 * @param[in] mapfac Map factors.
 * @param[in] ebfact EB factory.
 * @param[out] gradp Pressure gradient components.
 * @param[in] solverChoice Solver options.
 */
void
compute_gradp_xy (const MultiFab& p,
                  const Geometry& geom,
                  const MultiFab& z_phys_cc,
                  Vector<std::unique_ptr<MultiFab>>& mapfac,
                  const eb_& ebfact,
                  Vector<MultiFab>& gradp,
                  const SolverChoice& solverChoice)
{
    const bool l_use_terrain_fitted_coords = (solverChoice.mesh_type != MeshType::ConstantDz);

    const Box domain = geom.Domain();
    const int domain_klo = domain.smallEnd(2);
    const int domain_khi = domain.bigEnd(2);

    const GpuArray<Real, AMREX_SPACEDIM> dxInv = geom.InvCellSizeArray();

    // *****************************************************************************
    // Take gradient of relevant quantity (p0, pres, or pert_pres = pres - p0)
    // *****************************************************************************
    for ( MFIter mfi(p); mfi.isValid(); ++mfi)
    {
        Box tbx = mfi.nodaltilebox(0);
        Box tby = mfi.nodaltilebox(1);

        // Terrain metrics
        const Array4<const Real>& z_cc_arr = z_phys_cc.const_array(mfi);

        const Array4<const Real>& p_arr = p.const_array(mfi);

        const Array4<      Real>& gpx_arr = gradp[GpVars::gpx].array(mfi);
        const Array4<      Real>& gpy_arr = gradp[GpVars::gpy].array(mfi);

        const Array4<const Real>& mf_ux_arr = mapfac[MapFacType::u_x]->const_array(mfi);
        const Array4<const Real>& mf_vy_arr = mapfac[MapFacType::v_y]->const_array(mfi);

        if (solverChoice.terrain_type != TerrainType::EB) {

            ParallelFor(tbx, tby,
            [=] AMREX_GPU_DEVICE(int i, int j, int k) noexcept
            {
                //Note : mx/my == 1, so no map factor needed here
                Real gpx = dxInv[0] * (p_arr(i,j,k) - p_arr(i-1,j,k));

                if (l_use_terrain_fitted_coords) {
                    Real met_h_xi = (z_cc_arr(i,j,k) - z_cc_arr(i-1,j,k)) * dxInv[0];

                    Real dz_phys_hi, dz_phys_lo;
                    Real gpz_lo, gpz_hi;
                    if (k==domain_klo) {
                        dz_phys_hi = z_cc_arr(i  ,j,k+1) -   z_cc_arr(i  ,j,k  );
                        dz_phys_lo = z_cc_arr(i-1,j,k+1) -   z_cc_arr(i-1,j,k  );
                        gpz_hi  = (p_arr(i  ,j,k+1) - p_arr(i  ,j,k  )) / dz_phys_hi;
                        gpz_lo  = (p_arr(i-1,j,k+1) - p_arr(i-1,j,k  )) / dz_phys_lo;
                    } else if (k==domain_khi) {
                        dz_phys_hi = z_cc_arr(i  ,j,k  ) -   z_cc_arr(i  ,j,k-1);
                        dz_phys_lo = z_cc_arr(i-1,j,k  ) -   z_cc_arr(i-1,j,k-1);
                        gpz_hi  = (p_arr(i  ,j,k  ) - p_arr(i  ,j,k-1)) / dz_phys_hi;
                        gpz_lo  = (p_arr(i-1,j,k  ) - p_arr(i-1,j,k-1)) / dz_phys_lo;
                    } else {
                        dz_phys_hi = z_cc_arr(i  ,j,k+1) -   z_cc_arr(i  ,j,k-1);
                        dz_phys_lo = z_cc_arr(i-1,j,k+1) -   z_cc_arr(i-1,j,k-1);
                        gpz_hi  = (p_arr(i  ,j,k+1) - p_arr(i  ,j,k-1)) / dz_phys_hi;
                        gpz_lo  = (p_arr(i-1,j,k+1) - p_arr(i-1,j,k-1)) / dz_phys_lo;
                    }
                    Real gpx_metric = met_h_xi * myhalf * (gpz_hi + gpz_lo);
                    gpx -= gpx_metric;
                }
                gpx_arr(i,j,k) = gpx;

                // NOTE that the gradp array now carries the map factor!
                gpx_arr(i,j,k) *= mf_ux_arr(i,j,0);
            },
            [=] AMREX_GPU_DEVICE(int i, int j, int k) noexcept
            {
                //Note : mx/my == 1, so no map factor needed here
                Real gpy = dxInv[1] * (p_arr(i,j,k) - p_arr(i,j-1,k));

                if (l_use_terrain_fitted_coords) {
                    Real met_h_eta = (z_cc_arr(i,j,k) - z_cc_arr(i,j-1,k)) * dxInv[1];

                    Real dz_phys_hi, dz_phys_lo;
                    Real gpz_lo, gpz_hi;
                    if (k==domain_klo) {
                        dz_phys_hi = z_cc_arr(i,j  ,k+1) -   z_cc_arr(i,j  ,k  );
                        dz_phys_lo = z_cc_arr(i,j-1,k+1) -   z_cc_arr(i,j-1,k  );
                        gpz_hi  = (p_arr(i,j  ,k+1) - p_arr(i,j  ,k  )) / dz_phys_hi;
                        gpz_lo  = (p_arr(i,j-1,k+1) - p_arr(i,j-1,k  )) / dz_phys_lo;
                    } else if (k==domain_khi) {
                        dz_phys_hi = z_cc_arr(i,j  ,k  ) -   z_cc_arr(i,j  ,k-1);
                        dz_phys_lo = z_cc_arr(i,j-1,k  ) -   z_cc_arr(i,j-1,k-1);
                        gpz_hi  = (p_arr(i,j  ,k  ) - p_arr(i,j  ,k-1)) / dz_phys_hi;
                        gpz_lo  = (p_arr(i,j-1,k  ) - p_arr(i,j-1,k-1)) / dz_phys_lo;
                    } else {
                        dz_phys_hi = z_cc_arr(i,j  ,k+1) -   z_cc_arr(i,j  ,k-1);
                        dz_phys_lo = z_cc_arr(i,j-1,k+1) -   z_cc_arr(i,j-1,k-1);
                        gpz_hi  = (p_arr(i,j  ,k+1) - p_arr(i,j  ,k-1)) / dz_phys_hi;
                        gpz_lo  = (p_arr(i,j-1,k+1) - p_arr(i,j-1,k-1)) / dz_phys_lo;
                    }
                    Real gpy_metric = met_h_eta * myhalf * (gpz_hi + gpz_lo);
                    gpy -= gpy_metric;
                }
                gpy_arr(i,j,k) = gpy;

                // NOTE that the gradp array now carries the map factor!
                gpy_arr(i,j,k) *= mf_vy_arr(i,j,0);
            });

        } else {

            // Pressure gradients are fitted at the centroids of cut cells, if EB and Compressible.
            // Least-Squares Fitting: Compute slope using 3x3x3 stencil

            const bool l_fitting = false;

            const Real* dx_arr = geom.CellSize();
            const Real dx = dx_arr[0];
            const Real dy = dx_arr[1];
            const Real dz = dx_arr[2];

            // EB factory
            Array4<const EBCellFlag> cellflg = (ebfact.get_const_factory())->getMultiEBCellFlagFab()[mfi].const_array();

            // EB u-factory
            auto const* u_factory = ebfact.get_u_const_factory();
            Array4<const EBCellFlag> u_cellflg = u_factory->getMultiEBCellFlagFab()[mfi].const_array();
            Array4<const Real      > u_volfrac = u_factory->getVolFrac().const_array(mfi);
            bool u_is_cut = (u_factory->getMultiEBCellFlagFab()[mfi].getType() == FabType::singlevalued);
            Array4<const Real      > u_volcent = u_is_cut ? u_factory->getCentroid().const_array(mfi) : Array4<const Real>{};


            // EB v-factory
            auto const* v_factory = ebfact.get_v_const_factory();
            Array4<const EBCellFlag> v_cellflg = v_factory->getMultiEBCellFlagFab()[mfi].const_array();
            Array4<const Real      > v_volfrac = v_factory->getVolFrac().const_array(mfi);
            bool v_is_cut = (v_factory->getMultiEBCellFlagFab()[mfi].getType() == FabType::singlevalued);
            Array4<const Real      > v_volcent = v_is_cut ? v_factory->getCentroid().const_array(mfi) : Array4<const Real>{};

            if (l_fitting) {

                ParallelFor(tbx, tby,
                [=] AMREX_GPU_DEVICE(int i, int j, int k)
                {
                    if (u_volfrac(i,j,k) > zero) {

                        if (u_cellflg(i,j,k).isSingleValued()) {

                            GpuArray<Real,AMREX_SPACEDIM> slopes;
                            slopes = erf_calc_slopes_eb_staggered(Vars::xvel, Vars::cons, dx, dy, dz, i, j, k, p_arr, u_volcent, u_cellflg);

                            gpx_arr(i,j,k) = slopes[0];

                        } else {
                            gpx_arr(i,j,k) = dxInv[0] * (p_arr(i,j,k) - p_arr(i-1,j,k));
                        }

                    } else {
                        gpx_arr(i,j,k) = zero;
                    }
                },
                [=] AMREX_GPU_DEVICE(int i, int j, int k)
                {
                    if (v_volfrac(i,j,k) > zero) {

                        if (v_cellflg(i,j,k).isSingleValued()) {

                            GpuArray<Real,AMREX_SPACEDIM> slopes;
                            slopes = erf_calc_slopes_eb_staggered(Vars::yvel, Vars::cons, dx, dy, dz, i, j, k, p_arr, v_volcent, v_cellflg);

                            gpy_arr(i,j,k) = slopes[1];

                        } else {
                            gpy_arr(i,j,k) = dxInv[1] * (p_arr(i,j,k) - p_arr(i,j-1,k));
                        }
                    } else {
                        gpy_arr(i,j,k) = zero;
                    }
                });

            } else {

                // Simple calculation: assuming pressures at cell centers

                ParallelFor(tbx, tby,
                [=] AMREX_GPU_DEVICE(int i, int j, int k) noexcept
                {
                    if (u_volfrac(i,j,k) > zero) {
                        if (cellflg(i,j,k).isCovered()) {
                            gpx_arr(i,j,k) = dxInv[0] * (p_arr(i-3,j,k) - three*p_arr(i-2,j,k) + two*p_arr(i-1,j,k));
                        } else if (cellflg(i-1,j,k).isCovered()) {
                            gpx_arr(i,j,k) = dxInv[0] * (three*p_arr(i+1,j,k) - p_arr(i+2,j,k) - two*p_arr(i,j,k));
                        } else {
                            gpx_arr(i,j,k) = dxInv[0] * (p_arr(i,j,k) - p_arr(i-1,j,k));
                        }
                    } else {
                        gpx_arr(i,j,k) = zero;
                    }
                },
                [=] AMREX_GPU_DEVICE(int i, int j, int k) noexcept
                {
                    if (v_volfrac(i,j,k) > zero) {
                        if (cellflg(i,j,k).isCovered()) {
                            gpy_arr(i,j,k) = dxInv[1] * (p_arr(i,j-3,k) - three*p_arr(i,j-2,k) + two*p_arr(i,j-1,k));
                        } else if (cellflg(i,j-1,k).isCovered()) {
                            gpy_arr(i,j,k) = dxInv[1] * (three*p_arr(i,j+1,k) - p_arr(i,j+2,k) - two*p_arr(i,j,k));
                        } else {
                            gpy_arr(i,j,k) = dxInv[1] * (p_arr(i,j,k) - p_arr(i,j-1,k));
                        }
                    } else {
                        gpy_arr(i,j,k) = zero;
                    }
                });

            } // l_fitting

        } // TerrainType::EB

    } // mfi
}

/**
 * @brief Compute the vertical component of the pressure gradient.
 * @param[in] p Pressure field.
 * @param[in] geom Geometry container.
 * @param[in] z_phys_nd Physical height on nodes.
 * @param[in] ebfact EB factory.
 * @param[out] gradp Pressure gradient components.
 * @param[in] solverChoice Solver options.
 */
void
compute_gradp_z (const MultiFab& p,
                 const Geometry& geom,
                 const MultiFab& z_phys_nd,
                 const eb_& ebfact,
                 Vector<MultiFab>& gradp,
                 const SolverChoice& solverChoice)
{
    const bool l_use_terrain_fitted_coords = (solverChoice.mesh_type != MeshType::ConstantDz);

    const Box domain = geom.Domain();
    const int domain_klo = domain.smallEnd(2);
    const int domain_khi = domain.bigEnd(2);

    const GpuArray<Real, AMREX_SPACEDIM> dxInv = geom.InvCellSizeArray();

    // *****************************************************************************
    // Take gradient of relevant quantity (p0, pres, or pert_pres = pres - p0)
    // *****************************************************************************
    for ( MFIter mfi(p); mfi.isValid(); ++mfi)
    {
        Box tbz = mfi.nodaltilebox(2);

        // We don't compute gpz on the bottom or top domain boundary
        if (tbz.smallEnd(2) == domain_klo) {
            tbz.growLo(2,-1);
        }
        if (tbz.bigEnd(2) == domain_khi+1) {
            tbz.growHi(2,-1);
        }

        // Terrain metrics
        const Array4<const Real>& z_nd_arr = z_phys_nd.const_array(mfi);

        const Array4<const Real>& p_arr = p.const_array(mfi);

        const Array4<      Real>& gpz_arr = gradp[GpVars::gpz].array(mfi);

        if (solverChoice.terrain_type != TerrainType::EB) {

            ParallelFor(tbz, [=] AMREX_GPU_DEVICE(int i, int j, int k) noexcept
            {
                Real met_h_zeta = (l_use_terrain_fitted_coords) ? Compute_h_zeta_AtKface(i, j, k, dxInv, z_nd_arr) : 1;
                gpz_arr(i,j,k) = dxInv[2] * ( p_arr(i,j,k)-p_arr(i,j,k-1) )  / met_h_zeta;
            });

        } else {

            // Pressure gradients are fitted at the centroids of cut cells, if EB and Compressible.
            // Least-Squares Fitting: Compute slope using 3x3x3 stencil

            const bool l_fitting = false;

            const Real* dx_arr = geom.CellSize();
            const Real dx = dx_arr[0];
            const Real dy = dx_arr[1];
            const Real dz = dx_arr[2];

            // EB factory
            Array4<const EBCellFlag> cellflg = (ebfact.get_const_factory())->getMultiEBCellFlagFab()[mfi].const_array();

            // EB w-factory
            auto const* w_factory = ebfact.get_w_const_factory();
            Array4<const EBCellFlag> w_cellflg = w_factory->getMultiEBCellFlagFab()[mfi].const_array();
            Array4<const Real      > w_volfrac = w_factory->getVolFrac().const_array(mfi);
            bool w_is_cut = (w_factory->getMultiEBCellFlagFab()[mfi].getType() == FabType::singlevalued);
            Array4<const Real      > w_volcent = w_is_cut ? w_factory->getCentroid().const_array(mfi) : Array4<Real>{};

            if (l_fitting) {

                ParallelFor(tbz, [=] AMREX_GPU_DEVICE(int i, int j, int k)
                {
                    if (w_volfrac(i,j,k) > zero) {

                        if (w_cellflg(i,j,k).isSingleValued()) {

                            GpuArray<Real,AMREX_SPACEDIM> slopes;
                            slopes = erf_calc_slopes_eb_staggered(Vars::zvel, Vars::cons, dx, dy, dz, i, j, k, p_arr, w_volcent, w_cellflg);

                            gpz_arr(i,j,k) = slopes[2];

                        } else {
                            gpz_arr(i,j,k) = dxInv[2] * (p_arr(i,j,k) - p_arr(i,j,k-1));
                        }
                    } else {
                        gpz_arr(i,j,k) = zero;
                    }
                });

            } else {

                // Simple calculation: assuming pressures at cell centers

                ParallelFor(tbz, [=] AMREX_GPU_DEVICE(int i, int j, int k) noexcept
                {
                    if (w_volfrac(i,j,k) > zero) {
                        if (cellflg(i,j,k).isCovered()) {
                            gpz_arr(i,j,k) = dxInv[2] * ( p_arr(i,j,k-3) - three*p_arr(i,j,k-2) + two*p_arr(i,j,k-1) );
                        } else if (cellflg(i,j,k-1).isCovered()) {
                            gpz_arr(i,j,k) = dxInv[2] * ( three*p_arr(i,j,k+1) - p_arr(i,j,k+2) - two*p_arr(i,j,k) );
                        } else {
                            gpz_arr(i,j,k) = dxInv[2] * ( p_arr(i,j,k)-p_arr(i,j,k-1) );
                        }
                    } else {
                        gpz_arr(i,j,k) = zero;
                    }
                });

            } // l_fitting

        } // TerrainType::EB

    } // mfi
}

/**
 * @brief Compute the pressure gradient using vertical interpolation.
 * @param[in] p Pressure field.
 * @param[in] geom Geometry container.
 * @param[in] z_phys_nd Physical height on nodes.
 * @param[in] z_phys_cc Physical height on cell centers.
 * @param[in] mapfac Map factors.
 * @param[out] gradp Pressure gradient components.
 * @param[in] solverChoice Solver options.
 */
void
compute_gradp_interpz (const MultiFab& p,
                       const Geometry& geom,
                       const MultiFab& z_phys_nd,
                       const MultiFab& z_phys_cc,
                       Vector<std::unique_ptr<MultiFab>>& mapfac,
                       Vector<MultiFab>& gradp,
                       const SolverChoice& solverChoice)
{
    const bool l_use_terrain_fitted_coords = (solverChoice.mesh_type != MeshType::ConstantDz);

    const Box domain = geom.Domain();
    const int domain_klo = domain.smallEnd(2);
    const int domain_khi = domain.bigEnd(2);

    const GpuArray<Real, AMREX_SPACEDIM> dxInv = geom.InvCellSizeArray();

    // *****************************************************************************
    // Take gradient of relevant quantity (p0, pres, or pert_pres = pres - p0)
    // *****************************************************************************
    for ( MFIter mfi(p); mfi.isValid(); ++mfi)
    {
        Box tbx = mfi.nodaltilebox(0);
        Box tby = mfi.nodaltilebox(1);
        Box tbz = mfi.nodaltilebox(2);

        // We don't compute gpz on the bottom or top domain boundary
        if (tbz.smallEnd(2) == domain_klo) {
            tbz.growLo(2,-1);
        }
        if (tbz.bigEnd(2) == domain_khi+1) {
            tbz.growHi(2,-1);
        }

        // Terrain metrics
        const Array4<const Real>& z_nd_arr = z_phys_nd.const_array(mfi);
        const Array4<const Real>& z_cc_arr = z_phys_cc.const_array(mfi);

        const Array4<const Real>& p_arr = p.const_array(mfi);

        const Array4<      Real>& gpx_arr = gradp[GpVars::gpx].array(mfi);
        const Array4<      Real>& gpy_arr = gradp[GpVars::gpy].array(mfi);
        const Array4<      Real>& gpz_arr = gradp[GpVars::gpz].array(mfi);

        const Array4<const Real>& mf_ux_arr = mapfac[MapFacType::u_x]->const_array(mfi);
        const Array4<const Real>& mf_vy_arr = mapfac[MapFacType::v_y]->const_array(mfi);

        ParallelFor(tbx, tby, tbz,
        [=] AMREX_GPU_DEVICE(int i, int j, int k) noexcept
        {
            if (l_use_terrain_fitted_coords) {
                Real p_lo = p_arr(i-1,j,k);
                Real p_hi = p_arr(i,j,k);
                Real dz_int = myhalf * (z_cc_arr(i,j,k) - z_cc_arr(i-1,j,k));
                if (dz_int > 0) {
                    // Klemp 2011, Eqn. 16: s = 1/2
                    if (k==domain_klo) {
                        p_hi = quad_interp_1d(z_cc_arr(i,j,k) - dz_int,
                                              z_cc_arr(i,j,k  ), p_arr(i,j,k  ),
                                              z_cc_arr(i,j,k+1), p_arr(i,j,k+1),
                                              z_cc_arr(i,j,k+2), p_arr(i,j,k+2));
                    } else {
                        p_hi -= dz_int * ( (   p_arr(i  ,j,k  ) -    p_arr(i  ,j,k-1))
                                         / (z_cc_arr(i  ,j,k  ) - z_cc_arr(i  ,j,k-1)) );
                    }
                    if (k==domain_khi) {
                        p_lo = quad_interp_1d(z_cc_arr(i-1,j,k) + dz_int,
                                              z_cc_arr(i-1,j,k-2), p_arr(i-1,j,k-2),
                                              z_cc_arr(i-1,j,k-1), p_arr(i-1,j,k-1),
                                              z_cc_arr(i-1,j,k  ), p_arr(i-1,j,k  ));
                    } else {
                        p_lo += dz_int * ( (   p_arr(i-1,j,k+1) -    p_arr(i-1,j,k  ))
                                         / (z_cc_arr(i-1,j,k+1) - z_cc_arr(i-1,j,k  )) );
                    }
                } else if (dz_int < 0) {
                    // Klemp 2011, Eqn. 16: s = -1/2
                    if (k==domain_khi) {
                        p_hi = quad_interp_1d(z_cc_arr(i,j,k) - dz_int,
                                              z_cc_arr(i,j,k-2), p_arr(i,j,k-2),
                                              z_cc_arr(i,j,k-1), p_arr(i,j,k-1),
                                              z_cc_arr(i,j,k  ), p_arr(i,j,k  ));
                    } else {
                        p_hi -= dz_int * ( (   p_arr(i  ,j,k+1) -    p_arr(i  ,j,k  ))
                                         / (z_cc_arr(i  ,j,k+1) - z_cc_arr(i  ,j,k  )) );
                    }
                    if (k==domain_klo) {
                        p_lo = quad_interp_1d(z_cc_arr(i-1,j,k) + dz_int,
                                              z_cc_arr(i-1,j,k  ), p_arr(i-1,j,k  ),
                                              z_cc_arr(i-1,j,k+1), p_arr(i-1,j,k+1),
                                              z_cc_arr(i-1,j,k+2), p_arr(i-1,j,k+2));
                    } else {
                        p_lo += dz_int * ( (   p_arr(i-1,j,k  ) -    p_arr(i-1,j,k-1))
                                         / (z_cc_arr(i-1,j,k  ) - z_cc_arr(i-1,j,k-1)) );
                    }
                }
                gpx_arr(i,j,k) = dxInv[0] * (p_hi - p_lo);
            } else {
                gpx_arr(i,j,k) = dxInv[0] * (p_arr(i,j,k) - p_arr(i-1,j,k));
            }

            // NOTE that the gradp array now carries the map factor!
            gpx_arr(i,j,k) *= mf_ux_arr(i,j,0);
        },
        [=] AMREX_GPU_DEVICE(int i, int j, int k) noexcept
        {
            if (l_use_terrain_fitted_coords) {
                Real p_lo = p_arr(i,j-1,k);
                Real p_hi = p_arr(i,j,k);
                Real dz_int = myhalf * (z_cc_arr(i,j,k) - z_cc_arr(i,j-1,k));
                if (dz_int > 0) {
                    // Klemp 2011, Eqn. 16: s = 1/2
                    if (k==domain_klo) {
                        p_hi = quad_interp_1d(z_cc_arr(i,j,k) - dz_int,
                                              z_cc_arr(i,j,k  ), p_arr(i,j,k  ),
                                              z_cc_arr(i,j,k+1), p_arr(i,j,k+1),
                                              z_cc_arr(i,j,k+2), p_arr(i,j,k+2));
                    } else {
                        p_hi -= dz_int * ( (   p_arr(i,j  ,k  ) -    p_arr(i,j  ,k-1))
                                         / (z_cc_arr(i,j  ,k  ) - z_cc_arr(i,j  ,k-1)) );
                    }
                    if (k==domain_khi) {
                        p_lo = quad_interp_1d(z_cc_arr(i,j-1,k) + dz_int,
                                              z_cc_arr(i,j-1,k-2), p_arr(i,j-1,k-2),
                                              z_cc_arr(i,j-1,k-1), p_arr(i,j-1,k-1),
                                              z_cc_arr(i,j-1,k  ), p_arr(i,j-1,k  ));
                    } else {
                        p_lo += dz_int * ( (   p_arr(i,j-1,k+1) -    p_arr(i,j-1,k  ))
                                         / (z_cc_arr(i,j-1,k+1) - z_cc_arr(i,j-1,k  )) );
                    }
                } else if (dz_int < 0) {
                    // Klemp 2011, Eqn. 16: s = -1/2
                    if (k==domain_khi) {
                        p_hi = quad_interp_1d(z_cc_arr(i,j,k) - dz_int,
                                              z_cc_arr(i,j,k-2), p_arr(i,j,k-2),
                                              z_cc_arr(i,j,k-1), p_arr(i,j,k-1),
                                              z_cc_arr(i,j,k  ), p_arr(i,j,k  ));
                    } else {
                        p_hi -= dz_int * ( (   p_arr(i,j  ,k+1) -    p_arr(i,j  ,k  ))
                                         / (z_cc_arr(i,j  ,k+1) - z_cc_arr(i,j  ,k  )) );
                    }
                    if (k==domain_klo) {
                        p_lo = quad_interp_1d(z_cc_arr(i,j-1,k) + dz_int,
                                              z_cc_arr(i,j-1,k  ), p_arr(i,j-1,k  ),
                                              z_cc_arr(i,j-1,k+1), p_arr(i,j-1,k+1),
                                              z_cc_arr(i,j-1,k+2), p_arr(i,j-1,k+2));
                    } else {
                        p_lo += dz_int * ( (   p_arr(i,j-1,k  ) -    p_arr(i,j-1,k-1))
                                         / (z_cc_arr(i,j-1,k  ) - z_cc_arr(i,j-1,k-1)) );
                    }
                }
                gpy_arr(i,j,k) = dxInv[1] * (p_hi - p_lo);
            } else {
                gpy_arr(i,j,k) = dxInv[1] * (p_arr(i,j,k) - p_arr(i,j-1,k));
            }

            // NOTE that the gradp array now carries the map factor!
            gpy_arr(i,j,k) *= mf_vy_arr(i,j,0);
        },
        [=] AMREX_GPU_DEVICE(int i, int j, int k) noexcept
        {
            // Note: identical to gradp_type == 0
            Real met_h_zeta = (l_use_terrain_fitted_coords) ? Compute_h_zeta_AtKface(i, j, k, dxInv, z_nd_arr) : 1;
            gpz_arr(i,j,k) = dxInv[2] * ( p_arr(i,j,k)-p_arr(i,j,k-1) )  / met_h_zeta;
        });
    } // mfi
}

/**
 * @brief Compute the horizontal components of the pressure gradient, carrying the
 *        perturbational pressure to the height of the face along the hydrostatic
 *        relation dp'/dz = -g rho' rather than along a difference of p'.
 *
 * At an x-face (i,j,k) the two cell centers that straddle the face sit at different
 * physical heights, so each must be carried a signed distance
 *
 *     dz_int = 1/2 ( z_cc(i,j,k) - z_cc(i-1,j,k) )
 *
 * to reach the height of the face.  gradp_type 0 and 1 do this with a difference of
 * p' (centered, and one-sided respectively); here we instead integrate the hydrostatic
 * relation, so that
 *
 *     gpx = [p'(i,j,k) - p'(i-1,j,k)]/dx + g * met_h_xi * 1/2 * [rhot(i,j,k) + rhot(i-1,j,k)]
 *
 * where met_h_xi = 2*dz_int/dx and rhot is rho' evaluated at the midpoint of each of
 * the two extrapolation segments.  For gradp_type == 2 rho' is held constant over the
 * segment; for gradp_type == 3 it is reconstructed linearly, which makes the midpoint
 * rule exact and hence makes gpx vanish identically whenever rho' is linear in z and
 * p' is its exact hydrostatic integral, independently of the terrain slope.
 *
 * Note that the stencil needs no special casing at the bottom or top of the domain and
 * no branch on the sign of dz_int: rho' is defined in every cell, and the formula is
 * symmetric under dz_int -> -dz_int.  Only the piecewise-linear reconstruction used by
 * gradp_type == 3 reaches to k-1 and k+1, and there only for rho'.
 *
 * @param[in] p Perturbational pressure field.
 * @param[in] rhopert Perturbational (moist) density field.
 * @param[in] geom Geometry container.
 * @param[in] z_phys_cc Physical height on cell centers.
 * @param[in] mapfac Map factors.
 * @param[out] gradp Pressure gradient components.
 * @param[in] solverChoice Solver options.
 */
void
compute_gradp_hse (const MultiFab& p,
                   const MultiFab& rhopert,
                   const Geometry& geom,
                   const MultiFab& z_phys_cc,
                   Vector<std::unique_ptr<MultiFab>>& mapfac,
                   Vector<MultiFab>& gradp,
                   const SolverChoice& solverChoice)
{
    const bool l_use_terrain_fitted_coords = (solverChoice.mesh_type != MeshType::ConstantDz);

    // gradp_type == 3 reconstructs rho' linearly over the extrapolation segment;
    //                 gradp_type == 2 holds it constant
    const bool l_linear_rhopert = (solverChoice.gradp_type == 3);

    // Note this is the magnitude of gravity, i.e. it is positive
    const Real l_grav = solverChoice.gravity;

    const Box domain = geom.Domain();
    const int domain_klo = domain.smallEnd(2);
    const int domain_khi = domain.bigEnd(2);

    const GpuArray<Real, AMREX_SPACEDIM> dxInv = geom.InvCellSizeArray();

    for ( MFIter mfi(p); mfi.isValid(); ++mfi)
    {
        Box tbx = mfi.nodaltilebox(0);
        Box tby = mfi.nodaltilebox(1);

        // Terrain metrics
        const Array4<const Real>& z_cc_arr = z_phys_cc.const_array(mfi);

        const Array4<const Real>& p_arr  = p.const_array(mfi);
        const Array4<const Real>& rp_arr = rhopert.const_array(mfi);

        const Array4<      Real>& gpx_arr = gradp[GpVars::gpx].array(mfi);
        const Array4<      Real>& gpy_arr = gradp[GpVars::gpy].array(mfi);

        const Array4<const Real>& mf_ux_arr = mapfac[MapFacType::u_x]->const_array(mfi);
        const Array4<const Real>& mf_vy_arr = mapfac[MapFacType::v_y]->const_array(mfi);

        ParallelFor(tbx, tby,
        [=] AMREX_GPU_DEVICE(int i, int j, int k) noexcept
        {
            //Note : mx/my == 1, so no map factor needed here
            Real gpx = dxInv[0] * (p_arr(i,j,k) - p_arr(i-1,j,k));

            if (l_use_terrain_fitted_coords) {
                Real met_h_xi = (z_cc_arr(i,j,k) - z_cc_arr(i-1,j,k)) * dxInv[0];
                Real dz_int   = myhalf * (z_cc_arr(i,j,k) - z_cc_arr(i-1,j,k));

                // rho' at the midpoint of each of the two extrapolation segments
                Real rp_hi = rp_arr(i  ,j,k);
                Real rp_lo = rp_arr(i-1,j,k);
                if (l_linear_rhopert) {
                    rp_hi -= myhalf * dz_int * dqdz_cc(i  ,j,k,domain_klo,domain_khi,rp_arr,z_cc_arr);
                    rp_lo += myhalf * dz_int * dqdz_cc(i-1,j,k,domain_klo,domain_khi,rp_arr,z_cc_arr);
                }

                gpx += l_grav * met_h_xi * myhalf * (rp_hi + rp_lo);
            }
            gpx_arr(i,j,k) = gpx;

            // NOTE that the gradp array now carries the map factor!
            gpx_arr(i,j,k) *= mf_ux_arr(i,j,0);
        },
        [=] AMREX_GPU_DEVICE(int i, int j, int k) noexcept
        {
            //Note : mx/my == 1, so no map factor needed here
            Real gpy = dxInv[1] * (p_arr(i,j,k) - p_arr(i,j-1,k));

            if (l_use_terrain_fitted_coords) {
                Real met_h_eta = (z_cc_arr(i,j,k) - z_cc_arr(i,j-1,k)) * dxInv[1];
                Real dz_int    = myhalf * (z_cc_arr(i,j,k) - z_cc_arr(i,j-1,k));

                // rho' at the midpoint of each of the two extrapolation segments
                Real rp_hi = rp_arr(i,j  ,k);
                Real rp_lo = rp_arr(i,j-1,k);
                if (l_linear_rhopert) {
                    rp_hi -= myhalf * dz_int * dqdz_cc(i,j  ,k,domain_klo,domain_khi,rp_arr,z_cc_arr);
                    rp_lo += myhalf * dz_int * dqdz_cc(i,j-1,k,domain_klo,domain_khi,rp_arr,z_cc_arr);
                }

                gpy += l_grav * met_h_eta * myhalf * (rp_hi + rp_lo);
            }
            gpy_arr(i,j,k) = gpy;

            // NOTE that the gradp array now carries the map factor!
            gpy_arr(i,j,k) *= mf_vy_arr(i,j,0);
        });
    } // mfi
}
