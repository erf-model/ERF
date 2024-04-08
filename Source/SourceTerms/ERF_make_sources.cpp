#include <AMReX_MultiFab.H>
#include <AMReX_ArrayLim.H>
#include <AMReX_BCRec.H>
#include <AMReX_TableData.H>
#include <AMReX_GpuContainers.H>

#include <NumericalDiffusion.H>
#include <Src_headers.H>
#include <TI_slow_headers.H>

using namespace amrex;

/**
 * Function for computing the slow RHS for the evolution equations for the density, potential temperature and momentum.
 *
 * @param[in]  level level of resolution
 * @param[in]  nrk   which RK stage
 * @param[in]  dt    slow time step
 * @param[in]  S_data current solution
 * @param[in]  S_prim primitive variables (i.e. conserved variables divided by density)
 * @param[in] source source terms for conserved variables
 * @param[in]  geom   Container for geometric information
 * @param[in]  solverChoice  Container for solver parameters
 * @param[in] mapfac_u map factor at x-faces
 * @param[in] mapfac_v map factor at y-faces
 * @param[in] dptr_rhotheta_src  custom temperature source term
 * @param[in] dptr_rhoqt_src  custom moisture source term
 * @param[in] dptr_wbar_sub  subsidence source term
 * @param[in] d_rayleigh_ptrs_at_lev  Vector of {strength of Rayleigh damping, reference value of theta} used to define Rayleigh damping
 */

void make_sources (int level,
                   int /*nrk*/, Real dt,
                   Vector<MultiFab>& S_data,
                   const  MultiFab & S_prim,
                          MultiFab & source,
#ifdef ERF_USE_RRTMGP
                   const MultiFab& qheating_rates,
#endif
                   const Geometry geom,
                   const SolverChoice& solverChoice,
                   std::unique_ptr<MultiFab>& mapfac_u,
                   std::unique_ptr<MultiFab>& mapfac_v,
                   const Real* dptr_rhotheta_src,
                   const Real* dptr_rhoqt_src,
                   const Real* dptr_wbar_sub,
                   const Vector<Real*> d_rayleigh_ptrs_at_lev)
{
    BL_PROFILE_REGION("erf_make_sources()");

    // *****************************************************************************
    // Initialize source to zero since we re-compute it ever RK stage
    // *****************************************************************************
    source.setVal(0.0);

    const bool l_use_ndiff      = solverChoice.use_NumDiff;

    TurbChoice tc = solverChoice.turbChoice[level];
    const bool l_use_deardorff  = (tc.les_type == LESType::Deardorff);
    const bool l_use_QKE        = tc.use_QKE && tc.advect_QKE;

    const Box& domain = geom.Domain();

    const GpuArray<Real, AMREX_SPACEDIM> dxInv = geom.InvCellSizeArray();

    Real*      tau = d_rayleigh_ptrs_at_lev[Rayleigh::tau];
    Real* thetabar = d_rayleigh_ptrs_at_lev[Rayleigh::thetabar];

    // *****************************************************************************
    // Planar averages for subsidence terms
    // *****************************************************************************
    Table1D<Real>      dptr_t_plane, dptr_q_plane;
    TableData<Real, 1>  t_plane_tab,  q_plane_tab;
    if (dptr_wbar_sub)
    {
        PlaneAverage t_ave(&S_prim, geom, solverChoice.ave_plane, true);
        t_ave.compute_averages(ZDir(), t_ave.field());

        int ncell = t_ave.ncell_line();
        Gpu::HostVector<    Real> t_plane_h(ncell);
        Gpu::DeviceVector<  Real> t_plane_d(ncell);

        t_ave.line_average(PrimTheta_comp, t_plane_h);

        Gpu::copyAsync(Gpu::hostToDevice, t_plane_h.begin(), t_plane_h.end(), t_plane_d.begin());

        Real* dptr_t = t_plane_d.data();

        IntVect ng_c = S_prim.nGrowVect();
        Box tdomain = domain; tdomain.grow(2,ng_c[2]);
        t_plane_tab.resize({tdomain.smallEnd(2)}, {tdomain.bigEnd(2)});

        int offset = ng_c[2];
        dptr_t_plane = t_plane_tab.table();
        ParallelFor(ncell, [=] AMREX_GPU_DEVICE (int k) noexcept
        {
            dptr_t_plane(k-offset) = dptr_t[k];
        });

        if (solverChoice.moisture_type != MoistureType::None)
        {
            PlaneAverage q_ave(&S_prim, geom, solverChoice.ave_plane, true);
            q_ave.compute_averages(ZDir(), q_ave.field());

            Gpu::HostVector<  Real> q_plane_h(ncell);
            Gpu::DeviceVector<Real> q_plane_d(ncell);

            q_ave.line_average(PrimQ1_comp, q_plane_h);
            Gpu::copyAsync(Gpu::hostToDevice, q_plane_h.begin(), q_plane_h.end(), q_plane_d.begin());

            Real* dptr_q = q_plane_d.data();

            q_plane_tab.resize({tdomain.smallEnd(2)}, {tdomain.bigEnd(2)});

            dptr_q_plane = q_plane_tab.table();
            ParallelFor(ncell, [=] AMREX_GPU_DEVICE (int k) noexcept
            {
                dptr_q_plane(k-offset) = dptr_q[k];
            });
        }
    }

    // *****************************************************************************
    // Define source term for cell-centered conserved variables, from
    //    1. user-defined source terms for (rho theta) and (rho q_t)
    //    1. radiation           for (rho theta)
    //    2. Rayleigh damping    for (rho theta)
    //    3. custom forcing      for (rho theta) and (rho Q1)
    //    4. custom subsidence   for (rho theta) and (rho Q1)
    //    5. numerical diffusion for (rho theta)
    // *****************************************************************************

    // ***********************************************************************************************
    // Add remaining source terms
    // ***********************************************************************************************
#ifdef _OPENMP
#pragma omp parallel if (Gpu::notInLaunchRegion())
#endif
    {
    for ( MFIter mfi(S_data[IntVars::cons],TileNoZ()); mfi.isValid(); ++mfi)
    {
        Box bx  = mfi.tilebox();

        const Array4<const Real> & cell_data  = S_data[IntVars::cons].array(mfi);
        const Array4<const Real> & cell_prim  = S_prim.array(mfi);
        const Array4<Real>       & cell_src   = source.array(mfi);

#ifdef ERF_USE_RRTMGP
        // *************************************************************************************
        // Add radiation source terms to (rho theta)
        // *************************************************************************************
        {
            auto const& qheating_arr = qheating_rates.const_array(mfi);
            ParallelFor(bx, [=] AMREX_GPU_DEVICE (int i, int j, int k) noexcept
            {
                cell_src(i,j,k,RhoTheta_comp) += qheating_arr(i,j,k,0) + qheating_arr(i,j,k,1);
            });
        }

#endif

        // *************************************************************************************
        // Add Rayleigh damping for (rho theta)
        // *************************************************************************************
        if (solverChoice.use_rayleigh_damping && solverChoice.rayleigh_damp_T) {
            int n  = RhoTheta_comp;
            int nr = Rho_comp;
            int np = PrimTheta_comp;
            ParallelFor(bx, [=] AMREX_GPU_DEVICE (int i, int j, int k) noexcept
            {
                Real theta = cell_prim(i,j,k,np);
                cell_src(i, j, k, n) -= tau[k] * (theta - thetabar[k]) * cell_data(i,j,k,nr);
            });
        }

        // *************************************************************************************
        // Add custom forcing for (rho theta)
        // *************************************************************************************
        if (solverChoice.custom_rhotheta_forcing) {
            const int n = RhoTheta_comp;
            if (solverChoice.custom_forcing_prim_vars) {
                const int nr = Rho_comp;
                ParallelFor(bx, [=] AMREX_GPU_DEVICE (int i, int j, int k) noexcept
                {
                    cell_src(i, j, k, n) += cell_data(i,j,k,nr) * dptr_rhotheta_src[k];
                });
            } else {
                ParallelFor(bx, [=] AMREX_GPU_DEVICE (int i, int j, int k) noexcept
                {
                    cell_src(i, j, k, n) += dptr_rhotheta_src[k];
                });
            }
        }

        // *************************************************************************************
        // Add custom forcing for RhoQ1
        // *************************************************************************************
        if (solverChoice.custom_moisture_forcing) {
            const int n = RhoQ1_comp;
            if (solverChoice.custom_forcing_prim_vars) {
                const int nr = Rho_comp;
                ParallelFor(bx, [=] AMREX_GPU_DEVICE (int i, int j, int k) noexcept
                {
                    cell_src(i, j, k, n) += cell_data(i,j,k,nr) * dptr_rhoqt_src[k];
                });
            } else {
                ParallelFor(bx, [=] AMREX_GPU_DEVICE (int i, int j, int k) noexcept
                {
                    cell_src(i, j, k, n) += dptr_rhoqt_src[k];
                });
            }
        }

        // *************************************************************************************
        // Add custom subsidence for (rho theta)
        // *************************************************************************************
        if (solverChoice.custom_w_subsidence) {
            const int n = RhoTheta_comp;
            ParallelFor(bx, [=] AMREX_GPU_DEVICE (int i, int j, int k) noexcept
            {
                cell_src(i, j, k, n) += dptr_wbar_sub[k] *
                    0.5 * (dptr_t_plane(k+1) - dptr_t_plane(k-1)) * dxInv[2];

            });
        }

        // *************************************************************************************
        // Add custom subsidence for RhoQ1
        // *************************************************************************************
        if (solverChoice.custom_w_subsidence && (solverChoice.moisture_type != MoistureType::None)) {
            const int n = RhoQ1_comp;
            ParallelFor(bx, [=] AMREX_GPU_DEVICE (int i, int j, int k) noexcept
            {
                cell_src(i,j,k,n) += dptr_wbar_sub[k] *
                    0.5 * (dptr_q_plane(k+1) - dptr_q_plane(k-1)) * dxInv[2];
            });
        }

        // *************************************************************************************
        // Add numerical diffuion for rho and (rho theta)
        // *************************************************************************************
        if (l_use_ndiff) {
            int start_comp = 0;
            int   num_comp = 2;

            const Array4<const Real>& mf_u   = mapfac_u->const_array(mfi);
            const Array4<const Real>& mf_v   = mapfac_v->const_array(mfi);

            NumericalDiffusion(bx, start_comp, num_comp, dt, solverChoice.NumDiffCoeff,
                               cell_data, cell_src, mf_u, mf_v, false, false);

            if (l_use_deardorff) {
                int sc = RhoKE_comp;
                int nc = 1;
                NumericalDiffusion(bx, sc, nc, dt, solverChoice.NumDiffCoeff,
                                   cell_data, cell_src, mf_u, mf_v, false, false);
            }
            if (l_use_QKE) {
                int sc = RhoQKE_comp;
                int nc = 1;
                NumericalDiffusion(bx, sc, nc, dt, solverChoice.NumDiffCoeff,
                                   cell_data, cell_src, mf_u, mf_v, false, false);
            }
            {
                int sc = RhoScalar_comp;
                int nc = 1;
                NumericalDiffusion(bx, sc, nc, dt, solverChoice.NumDiffCoeff,
                                   cell_data, cell_src, mf_u, mf_v, false, false);
            }
        }

        // *************************************************************************************
        // Add sponging
        // *************************************************************************************
        ApplySpongeZoneBCsForCC(solverChoice.spongeChoice, geom, bx, cell_src, cell_data);

    } // mfi
    } // OMP
}
