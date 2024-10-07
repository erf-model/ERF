/**
 * \file ERF_init1d.cpp
 */
#include <ERF_EOS.H>
#include <ERF.H>
#include <ERF_TileNoZ.H>
#include <ERF_prob_common.H>
#include <ERF_ParFunctions.H>
#include <ERF_Utils.H>

#include <ERF_Interpolation_1D.H>

using namespace amrex;

/**
 * Initialize density and pressure base state in
 * hydrostatic equilibrium.
 */
void
ERF::initHSE (int lev)
{
    // This integrates up through column to update p_hse, pi_hse;
    // r_hse is not const b/c FillBoundary is called at the end for r_hse and p_hse

    MultiFab r_hse (base_state[lev], make_alias, 0, 1); // r_0  is first  component
    MultiFab p_hse (base_state[lev], make_alias, 1, 1); // p_0  is second component
    MultiFab pi_hse(base_state[lev], make_alias, 2, 1); // pi_0 is third  component

    bool all_boxes_touch_bottom = true;
    Box domain(geom[lev].Domain());

    int icomp = 0; int bccomp = BCVars::base_bc; int ncomp = 3;

    Real time = 0.;

    if (lev == 0) {
        BoxArray ba(base_state[lev].boxArray());
        for (int i = 0; i < ba.size(); i++) {
            if (ba[i].smallEnd(2) != domain.smallEnd(2)) {
                all_boxes_touch_bottom = false;
            }
        }
    }
    else
    {
        //
        // We need to do this interp from coarse level in order to set the values of
        // the base state inside the domain but outside of the fine region
        //
        base_state[lev-1].FillBoundary(geom[lev-1].periodicity());
        //
        // NOTE: this interpolater assumes that ALL ghost cells of the coarse MultiFab
        //       have been pre-filled - this includes ghost cells both inside and outside
        //       the domain
        //
        InterpFromCoarseLevel(base_state[lev], base_state[lev].nGrowVect(),
                              IntVect(0,0,0), // do not fill ghost cells outside the domain
                              base_state[lev-1], icomp, icomp, ncomp,
                              geom[lev-1], geom[lev],
                              refRatio(lev-1), &cell_cons_interp,
                              domain_bcs_type, BCVars::cons_bc);

         // We need to do this here because the interpolation above may leave corners unfilled
         //    when the corners need to be filled by, for example, reflection of the fine ghost
         //    cell outside the fine region but inide the domain.
         (*physbcs_base[lev])(base_state[lev],icomp,ncomp,base_state[lev].nGrowVect(),time,bccomp);
    }

    if (all_boxes_touch_bottom || lev > 0) {

        // Initial r_hse may or may not be in HSE -- defined in ERF_prob.cpp
        if (solverChoice.use_moist_background){
            prob->erf_init_dens_hse_moist(r_hse, z_phys_nd[lev], geom[lev]);
        } else {
            prob->erf_init_dens_hse(r_hse, z_phys_nd[lev], z_phys_cc[lev], geom[lev]);
        }

        erf_enforce_hse(lev, r_hse, p_hse, pi_hse, z_phys_cc[lev]);

    } else {

        BoxArray ba_new(domain);

        ChopGrids2D(ba_new, domain, ParallelDescriptor::NProcs());

        DistributionMapping dm_new(ba_new);

        MultiFab new_base_state(ba_new, dm_new, 3, 1);
        new_base_state.ParallelCopy(base_state[lev],0,0,3,1,1);

        MultiFab new_r_hse (new_base_state, make_alias, 0, 1); // r_0  is first  component
        MultiFab new_p_hse (new_base_state, make_alias, 1, 1); // p_0  is second component
        MultiFab new_pi_hse(new_base_state, make_alias, 2, 1); // pi_0 is third  component

        std::unique_ptr<MultiFab> new_z_phys_cc;
        std::unique_ptr<MultiFab> new_z_phys_nd;
        if (solverChoice.use_terrain) {
            new_z_phys_cc = std::make_unique<MultiFab>(ba_new,dm_new,1,1);
            new_z_phys_cc->ParallelCopy(*z_phys_cc[lev],0,0,1,1,1);

            BoxArray ba_new_nd(ba_new);
            ba_new_nd.surroundingNodes();
            new_z_phys_nd = std::make_unique<MultiFab>(ba_new_nd,dm_new,1,1);
            new_z_phys_nd->ParallelCopy(*z_phys_nd[lev],0,0,1,1,1);
        }

        // Initial r_hse may or may not be in HSE -- defined in ERF_prob.cpp
        if (solverChoice.use_moist_background){
            prob->erf_init_dens_hse_moist(new_r_hse, new_z_phys_nd, geom[lev]);
        } else {
            prob->erf_init_dens_hse(new_r_hse, new_z_phys_nd, new_z_phys_cc, geom[lev]);
        }

        erf_enforce_hse(lev, new_r_hse, new_p_hse, new_pi_hse, new_z_phys_cc);

        // Now copy back into the original arrays
        base_state[lev].ParallelCopy(new_base_state,0,0,3,1,1);
    }

    //
    // Impose physical bc's on the base state -- the values outside the fine region
    //   but inside the domain have already been filled in the call above to InterpFromCoarseLevel
    //
    (*physbcs_base[lev])(base_state[lev],icomp,ncomp,base_state[lev].nGrowVect(),time,bccomp);

    base_state[lev].FillBoundary(geom[lev].periodicity());
}

void
ERF::initHSE ()
{
    AMREX_ALWAYS_ASSERT(!init_sounding_ideal);
    for (int lev = 0; lev <= finest_level; lev++)
    {
        initHSE(lev);
    }
}

/**
 * Enforces hydrostatic equilibrium when using terrain.
 *
 * @param[in]  lev  Integer specifying the current level
 * @param[out] dens MultiFab storing base state density
 * @param[out] pres MultiFab storing base state pressure
 * @param[out] pi   MultiFab storing base state Exner function
 * @param[in]  z_cc Pointer to MultiFab storing cell centered z-coordinates
 */
void
ERF::erf_enforce_hse (int lev,
                      MultiFab& dens, MultiFab& pres, MultiFab& pi,
                      std::unique_ptr<MultiFab>& z_cc)
{
    Real l_gravity = solverChoice.gravity;
    bool l_use_terrain = solverChoice.use_terrain;

    const auto geomdata = geom[lev].data();
    const Real dz = geomdata.CellSize(2);

    const Box& domain = geom[lev].Domain();

    for ( MFIter mfi(dens, TileNoZ()); mfi.isValid(); ++mfi )
    {
        // Create a flat box with same horizontal extent but only one cell in vertical
        const Box& tbz = mfi.nodaltilebox(2);
        int klo = tbz.smallEnd(2);
        int khi = tbz.bigEnd(2);

        Box b2d = tbz; // Copy constructor
        b2d.grow(0,1); b2d.grow(1,1); // Grow by one in the lateral directions
        b2d.setRange(2,0);

        // We integrate to the first cell (and below) by using rho in this cell
        // If gravity == 0 this is constant pressure
        // If gravity != 0, hence this is a wall, this gives gp0 = dens[0] * gravity
        // (dens_hse*gravity would also be dens[0]*gravity because we use foextrap for rho at k = -1)
        // Note ng_pres_hse = 1

       // We start by assuming pressure on the ground is p_0 (in ERF_Constants.H)
       // Note that gravity is positive

        Array4<Real>  rho_arr = dens.array(mfi);
        Array4<Real> pres_arr = pres.array(mfi);
        Array4<Real>   pi_arr =   pi.array(mfi);
        Array4<Real> zcc_arr;
        if (l_use_terrain) {
           zcc_arr = z_cc->array(mfi);
        }

        const Real rdOcp = solverChoice.rdOcp;

        ParallelFor(b2d, [=] AMREX_GPU_DEVICE (int i, int j, int)
        {
            // Set value at surface from Newton iteration for rho
            if (klo == 0) {
                // Physical height of the terrain at cell center
                Real hz;
                if (l_use_terrain) {
                    hz = zcc_arr(i,j,klo);
                } else {
                    hz = 0.5*dz;
                }

                pres_arr(i,j,klo  ) = p_0 - hz * rho_arr(i,j,klo) * l_gravity;
                pi_arr(i,j,klo  ) = getExnergivenP(pres_arr(i,j,klo  ), rdOcp);

                // Set ghost cell with dz and rho at boundary
                pres_arr(i,j,klo-1) = p_0 + hz * rho_arr(i,j,klo) * l_gravity;
                pi_arr(i,j,klo-1) = getExnergivenP(pres_arr(i,j,klo-1), rdOcp);

            } else {
                // If level > 0 and klo > 0, we need to use the value of pres_arr(i,j,klo-1) which was
                //    filled from FillPatch-ing it.
                Real dz_loc;
                if (l_use_terrain) {
                    dz_loc = (zcc_arr(i,j,klo) - zcc_arr(i,j,klo-1));
                } else {
                    dz_loc = dz;
                }

                Real dens_interp = 0.5*(rho_arr(i,j,klo) + rho_arr(i,j,klo-1));
                pres_arr(i,j,klo) = pres_arr(i,j,klo-1) - dz_loc * dens_interp * l_gravity;

                pi_arr(i,j,klo  ) = getExnergivenP(pres_arr(i,j,klo  ), rdOcp);
                pi_arr(i,j,klo-1) = getExnergivenP(pres_arr(i,j,klo-1), rdOcp);
            }

            Real dens_interp;
            if (l_use_terrain) {
                for (int k = klo+1; k <= khi; k++) {
                    Real dz_loc = (zcc_arr(i,j,k) - zcc_arr(i,j,k-1));
                    dens_interp = 0.5*(rho_arr(i,j,k) + rho_arr(i,j,k-1));
                    pres_arr(i,j,k) = pres_arr(i,j,k-1) - dz_loc * dens_interp * l_gravity;
                    pi_arr(i,j,k) = getExnergivenP(pres_arr(i,j,k), rdOcp);
                }
            } else {
                for (int k = klo+1; k <= khi; k++) {
                    dens_interp = 0.5*(rho_arr(i,j,k) + rho_arr(i,j,k-1));
                    pres_arr(i,j,k) = pres_arr(i,j,k-1) - dz * dens_interp * l_gravity;
                    pi_arr(i,j,k) = getExnergivenP(pres_arr(i,j,k), rdOcp);
                }
            }
        });

        int domlo_x = domain.smallEnd(0); int domhi_x = domain.bigEnd(0);
        int domlo_y = domain.smallEnd(1); int domhi_y = domain.bigEnd(1);

        if (pres[mfi].box().smallEnd(0) < domlo_x)
        {
            Box bx = mfi.nodaltilebox(2);
            bx.setSmall(0,domlo_x-1);
            bx.setBig(0,domlo_x-1);
            ParallelFor(bx, [=] AMREX_GPU_DEVICE (int i, int j, int k) {
                pres_arr(i,j,k) = pres_arr(domlo_x,j,k);
                pi_arr(i,j,k) = getExnergivenP(pres_arr(i,j,k), rdOcp);
            });
        }

        if (pres[mfi].box().bigEnd(0) > domhi_x)
        {
            Box bx = mfi.nodaltilebox(2);
            bx.setSmall(0,domhi_x+1);
            bx.setBig(0,domhi_x+1);
            ParallelFor(bx, [=] AMREX_GPU_DEVICE (int i, int j, int k) {
                pres_arr(i,j,k) = pres_arr(domhi_x,j,k);
                pi_arr(i,j,k) = getExnergivenP(pres_arr(i,j,k), rdOcp);
            });
        }

        if (pres[mfi].box().smallEnd(1) < domlo_y)
        {
            Box bx = mfi.nodaltilebox(2);
            bx.setSmall(1,domlo_y-1);
            bx.setBig(1,domlo_y-1);
            ParallelFor(bx, [=] AMREX_GPU_DEVICE (int i, int j, int k) {
                pres_arr(i,j,k) = pres_arr(i,domlo_y,k);
                pi_arr(i,j,k) = getExnergivenP(pres_arr(i,j,k), rdOcp);
            });
        }

        if (pres[mfi].box().bigEnd(1) > domhi_y)
        {
            Box bx = mfi.nodaltilebox(2);
            bx.setSmall(1,domhi_y+1);
            bx.setBig(1,domhi_y+1);
            ParallelFor(bx, [=] AMREX_GPU_DEVICE (int i, int j, int k) {
                pres_arr(i,j,k) = pres_arr(i,domhi_y,k);
                pi_arr(i,j,k) = getExnergivenP(pres_arr(i,j,k), rdOcp);
            });
        }
    }

    dens.FillBoundary(geom[lev].periodicity());
    pres.FillBoundary(geom[lev].periodicity());
}
