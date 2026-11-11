/**
 * \file ERF_ImmersedProjection.cpp
 *
 * Implicit immersed-forcing drag in the anelastic projection (erf.if_implicit_projection,
 * erf-model/ERF#4919 for terrain, #4887 for buildings).
 *
 * The terrain immersed forcing drags the momentum of the solid faces toward zero at the rate
 * rate = beta_face * if_Cd_momentum / (dx dy dz)^(1/3). Applied explicitly, the drag is undone by
 * the projection, which treats the solid as fluid; the post-projection freeze of the solid
 * faces then leaves the surface faces (solid on one side only) unbalanced, so the solid row
 * under the surface becomes a mass sink and source. With the drag inside the projection,
 *
 *     m^{n+1} = (m* - grad phi) / (1 + dt rate) = sigma (m* - grad phi),
 *     div(sigma grad phi) = div(sigma m*),
 *
 * the projected momenta are divergence free everywhere, the solid included, and no freeze is
 * needed. The faces where the source term applies the wall law instead of the drag keep sigma = 1.
 * For buildings (erf.buildings_type = ImmersedForcing, erf-model/ERF#4887) the rate is the linear
 * drag of ImmersedForcingBuildings_{X,Y,Z}mom (partial-wall and interior faces, from the masks of
 * ERF_ImmersedBuildingMasks.H), and the faces between two solid cells take the full rate.
 */
#include "ERF.H"
#include "ERF_ImmersedWallCell.H"
#include "ERF_ImmersedBuildingMasks.H"

#ifdef ERF_USE_FFT
#include "ERF_ImmersedPoisson.H"
#include <AMReX_GMRES.H>
#endif

using namespace amrex;

/**
 * Face coefficients sigma = 1 / (1 + dt rate) of the implicit immersed drag in the projection,
 * with rate the drag rate of ImmersedForcingTerrain_{X,Y,Z}mom (or, for buildings,
 * ImmersedForcingBuildings_{X,Y,Z}mom) on that face, zero where those apply the wall law.
 * The faces on the boundary of the level's grids keep sigma = 1.
 *
 * @param lev   Level index
 * @param l_dt  Time step of the projection
 * @param sigma Face coefficients on the x-, y- and z-faces of the level's grids (defined here)
 */
void ERF::make_if_projection_sigma (int lev, Real l_dt, Array<MultiFab,AMREX_SPACEDIM>& sigma)
{
    BL_PROFILE("ERF::make_if_projection_sigma()");

    AMREX_ALWAYS_ASSERT(terrain_blanking[lev]);
    const MultiFab& tb = *terrain_blanking[lev];
    const MultiFab* tb_face[AMREX_SPACEDIM] = {terrain_blanking_xface[lev].get(),
                                               terrain_blanking_yface[lev].get(),
                                               terrain_blanking_zface[lev].get()};

    // 1 in the valid cells of this level (and their periodic images), 0 elsewhere: a face is
    // interior to the level's grids when both of its cells are 1. The other faces carry the
    // domain or c/f boundary data of the projection and keep sigma = 1.
    iMultiFab inside(grids[lev], dmap[lev], 1, 1);
    inside.setVal(0);
    inside.setVal(1, 0, 1, 0);
    inside.FillBoundary(geom[lev].periodicity());

    for (int dir = 0; dir < AMREX_SPACEDIM; ++dir) {
        sigma[dir].define(convert(grids[lev], IntVect::TheDimensionVector(dir)), dmap[lev], 1, 0);
    }

    // The drag rate of a solid face in the source terms, without their windspeed cap
    // (CdM * windspeed = min(1, |u|) * drag_coefficient there)
    const auto dx = geom[lev].CellSizeArray();
    const Real rate_solid    = solverChoice.if_Cd_momentum / std::pow(dx[0]*dx[1]*dx[2], one/three);
    const Real small_volfrac = solverChoice.if_small_volfrac;
    const bool use_most      = solverChoice.if_use_most;
    // the fraction-stress law (erf.if_wall_form = fraction_stress) has its own wall cell
    const bool fraction_stress = solverChoice.if_fraction_stress;
    // Buildings (erf.buildings_type = ImmersedForcing, erf-model/ERF#4887): the linear drag of
    // ImmersedForcingBuildings_{X,Y,Z}mom, from the masks of ERF_ImmersedBuildingMasks.H, read
    // with the blanking as those kernels read it (snapped with erf.if_snap_partial_cells)
    // (a level above the 50 m cutoff of the building forcing has no blanking, so every face
    // gets sigma = 1 there, as the kernels and the explicit path leave it alone)
    const bool buildings = (solverChoice.buildings_type == BuildingsType::ImmersedForcing);
    const bool snap      = solverChoice.if_snap_partial_cells;
    const int  dom_klo = geom[lev].Domain().smallEnd(2);
    const int  dom_khi = geom[lev].Domain().bigEnd(2);

#ifdef _OPENMP
#pragma omp parallel if (Gpu::notInLaunchRegion())
#endif
    for (MFIter mfi(tb, TilingIfNotGPU()); mfi.isValid(); ++mfi)
    {
        const Array4<const Real> t  = tb.const_array(mfi);
        const Array4<const int>  in = inside.const_array(mfi);

        for (int dir = 0; dir < AMREX_SPACEDIM; ++dir)
        {
            const int io = (dir == 0) ? 1 : 0;
            const int jo = (dir == 1) ? 1 : 0;
            const int ko = (dir == 2) ? 1 : 0;
            const Array4<const Real> tf = (tb_face[dir]) ? tb_face[dir]->const_array(mfi) : Array4<const Real>{};
            const Array4<Real>       sg = sigma[dir].array(mfi);

            if (buildings) {
                ParallelFor(mfi.nodaltilebox(dir), [=] AMREX_GPU_DEVICE (int i, int j, int k) noexcept
                {
                    if (in(i,j,k) == 0 || in(i-io,j-jo,k-ko) == 0) {
                        sg(i,j,k) = one;
                        return;
                    }
                    // the face as the building kernels read it, and the weight of the drag they drop
                    Real w;
                    if (dir == 0) {
                        const if_bld::HorizFace f = if_bld::horiz_face<0>(t, tf, i, j, k, dom_klo, snap, small_volfrac);
                        w = if_bld::implicit_weight(if_bld::horiz_masks(f, use_most, snap).linear(f.t), f.raw);
                    } else if (dir == 1) {
                        const if_bld::HorizFace f = if_bld::horiz_face<1>(t, tf, i, j, k, dom_klo, snap, small_volfrac);
                        w = if_bld::implicit_weight(if_bld::horiz_masks(f, use_most, snap).linear(f.t), f.raw);
                    } else {
                        const if_bld::VertFace f = if_bld::vert_face(t, tf, i, j, k, snap, small_volfrac);
                        w = if_bld::implicit_weight(if_bld::vert_masks(f, k >= dom_klo + 1, use_most, snap).linear(f.t), f.raw);
                    }
                    sg(i,j,k) = one / (one + l_dt * w * rate_solid);
                });
                continue;
            }

            ParallelFor(mfi.nodaltilebox(dir), [=] AMREX_GPU_DEVICE (int i, int j, int k) noexcept
            {
                if (in(i,j,k) == 0 || in(i-io,j-jo,k-ko) == 0) {
                    sg(i,j,k) = one;
                    return;
                }
                // The face blanking as the terrain kernels read it: face-centered if available,
                // otherwise the average of the two cells, thresholded at small_volfrac
                auto face_blank = [=] (int kk) -> Real {
                    const int kc = amrex::min(amrex::max(kk, dom_klo), dom_khi);
                    const Real b = (tf) ? tf(i,j,kc) : myhalf * (t(i,j,kc) + t(i-io,j-jo,kc-ko));
                    return (b < small_volfrac) ? zero : b;
                };
                const Real b = face_blank(k);
                // The horizontal kernels apply the wall law instead of the drag on a partial face
                // under a fluid face (a fully solid face there is held by the drag in the projection);
                // the fraction-stress law forces its wall cell (ib_wall::is_wall_cell), the fluid face
                // on a solid one included, and treats a face within small_volfrac of solid as solid
                const bool wall = (dir < 2) && use_most &&
                                  (fraction_stress ? ib_wall::is_wall_cell(face_blank(k-1), b, face_blank(k+1), small_volfrac)
                                                   : ((b > zero) && (b < one) && (face_blank(k+1) == zero)));
                const Real rate = wall ? zero : b * rate_solid;
                sg(i,j,k) = one / (one + l_dt * rate);
            });
        }
    }
}

/**
 * FFT-preconditioned GMRES solve of div(sigma grad phi) = rhs (erf.if_implicit_projection)
 *
 * @param lev       Level index for the solve
 * @param subdomain Box over which the solve is performed
 * @param rhs       Right-hand side, div(sigma m*)
 * @param phi       Solution field to fill
 * @param fluxes    Face fluxes to fill: -grad phi, without sigma
 * @param sigma     Face coefficients on the subdomain (make_if_projection_sigma)
 */
void ERF::solve_with_if_gmres (int lev, const Box& subdomain, MultiFab& rhs, MultiFab& phi,
                               Array<MultiFab,AMREX_SPACEDIM>& fluxes,
                               Array<MultiFab const*,AMREX_SPACEDIM> const& sigma)
{
#ifdef ERF_USE_FFT
    BL_PROFILE("ERF::solve_with_if_gmres()");

    auto const dom_lo = lbound(Geom(lev).Domain());
    auto const dom_hi = ubound(Geom(lev).Domain());
    auto const sub_lo = lbound(subdomain);
    auto const sub_hi = ubound(subdomain);
    auto dx = Geom(lev).CellSizeArray();

    // The same subdomain geometry as solve_with_gmres
    Array<int,AMREX_SPACEDIM> is_per; is_per[0] = 0; is_per[1] = 0; is_per[2] = 0;
    if (Geom(lev).isPeriodic(0) && sub_lo.x == dom_lo.x && sub_hi.x == dom_hi.x) { is_per[0] = 1;}
    if (Geom(lev).isPeriodic(1) && sub_lo.y == dom_lo.y && sub_hi.y == dom_hi.y) { is_per[1] = 1;}
    int coord_sys = 0;
    Geometry my_geom;
    if (subdomain == Geom(lev).Domain()) {
        my_geom.define(Geom(lev).Domain(), Geom(lev).ProbDomain(), coord_sys, is_per);
    } else {
        RealBox rb( sub_lo.x   *dx[0],  sub_lo.y   *dx[1],  sub_lo.z   *dx[2],
                   (sub_hi.x+1)*dx[0], (sub_hi.y+1)*dx[1], (sub_hi.z+1)*dx[2]);
        my_geom.define(subdomain, rb, coord_sys, is_per);
    }

    ImmersedPoisson ip(my_geom, Geom(lev), rhs.boxArray(), rhs.DistributionMap(), domain_bc_type,
                       sigma, solverChoice.use_real_bcs);
    ip.usePrecond(true);

    GMRES<MultiFab, ImmersedPoisson> gmsolver;
    gmsolver.define(ip);
    gmsolver.setVerbose(mg_verbose);
    gmsolver.setRestartLength(50);
    gmsolver.solve(phi, rhs, solverChoice.poisson_reltol, solverChoice.poisson_abstol);
    // GMRES returns quietly when it stops at its iteration cap; the projection then would not be
    // divergence free, so stop here (the FFT preconditioner keeps the count at about ten)
    if (gmsolver.getStatus() != 0) {
        amrex::Abort("erf.if_implicit_projection: the GMRES solve of div(sigma grad phi) = rhs did not converge "
                     "(residual " + std::to_string(gmsolver.getResidualNorm()) + " after " +
                     std::to_string(gmsolver.getNumIters()) + " iterations)");
    }
    if (mg_verbose > 0) {
        amrex::Print() << "erf.if_implicit_projection: GMRES converged in " << gmsolver.getNumIters()
                       << " iterations, residual " << gmsolver.getResidualNorm() << "\n";
    }

    ip.getFluxes(phi, fluxes);
#else
    amrex::ignore_unused(lev, rhs, phi, fluxes, sigma);
    amrex::Abort("erf.if_implicit_projection needs a build with FFT");
#endif

    ImposeBCsOnPhi(lev, phi, subdomain);
}
