/**
 * \file ERF_CoarseFineFaces.cpp
 *
 * A refined level that stops inside the domain in z (a level refined over part of the depth)
 * has coarse/fine faces at its top and bottom.  On a terrain-fitted mesh the anelastic
 * projection works with the contravariant flux through each z-face,
 *
 *     Omega = rho0 w - rho0 u dz/dx - rho0 v dz/dy,
 *
 * which needs the horizontal momenta one face below and above the level and, for the solve to
 * have a solution, fluxes through the level's boundary that add up to the coarse level's.
 * fill_cf_zghost_momenta supplies the first; set_cf_zface_omega_projection and
 * set_cf_zface_omega_fill the second, the same way in the projection and in the fine-level
 * fills the time step uses.
 */
#include "ERF.H"
#include "ERF_Utils.H"

#include <AMReX_FillPatchUtil.H>
#include <AMReX_Interpolater.H>

using namespace amrex;

namespace {

/**
 * The coarse/fine w-faces of every box of a level: its bottom (top) face where the cells just
 * below (above) it are inside the domain and not in another box of the level.
 *
 * @param fba  BoxArray of the level
 * @param fdom Domain of the level
 * @return     For each box index, the nodal-in-z boxes of its coarse/fine w-faces
 */
std::map<int, Vector<Box>>
cf_zface_patches (const BoxArray& fba, const Box& fdom)
{
    std::map<int, Vector<Box>> out;
    for (int n = 0; n < fba.size(); ++n) {
        const Box& b = fba[n];
        for (int side = 0; side < 2; ++side) {
            const int kc = (side == 0) ? b.smallEnd(2) - 1 : b.bigEnd(2) + 1;
            const int kf = (side == 0) ? b.smallEnd(2)     : b.bigEnd(2) + 1;
            if (kc < fdom.smallEnd(2) || kc > fdom.bigEnd(2)) { continue; }
            Box outside(b); outside.setRange(2, kc);
            BoxList uncovered = fba.complementIn(outside);
            for (Box u : uncovered) {
                u.surroundingNodes(2);
                u.setRange(2, kf);
                out[n].push_back(u);
            }
        }
    }
    return out;
}

} // namespace

/**
 * Fill the (rho0 u) and (rho0 v) faces in the z-ghost layer of a refined level from the coarse
 * level, where the level stops inside the domain in z.
 *
 * OmegaFromW and WFromOmega read one horizontal-momentum face below the lowest and above the
 * highest w-face of each box.  FillBoundary fills that face only where a neighbouring fine box
 * covers it; at a coarse/fine face in z it would keep whatever the array held (at start-up,
 * memory never written), and the terrain slope multiplies it into Omega.  Here the coarse
 * velocity, interpolated in space and time as the fine-level fill does, times the density on the
 * face, is written there.  Faces covered by another fine box are overwritten by the FillBoundary
 * that follows; faces outside the domain are left alone.
 *
 * @param lev   Level index (> 0)
 * @param time  Time of the momenta, for the coarse time interpolation
 * @param rho0u x-face momentum in the scaling of its valid faces
 * @param rho0v y-face momentum in the scaling of its valid faces
 * @param rho0  Cell-centred density of that scaling (component 0), with ghost cells
 */
void ERF::fill_cf_zghost_momenta (int lev, double time_d,
                                  MultiFab& rho0u, MultiFab& rho0v, MultiFab const& rho0)
{
    BL_PROFILE("ERF::fill_cf_zghost_momenta()");
    AMREX_ALWAYS_ASSERT(lev > 0);

    const Box& fdom = geom[lev].Domain();
    const int klo_dom = fdom.smallEnd(2);
    const int khi_dom = fdom.bigEnd(2);

    // Nothing to do when every box of the level spans the whole depth
    bool any_inside = false;
    for (int n = 0; n < grids[lev].size(); ++n) {
        const Box& b = grids[lev][n];
        if (b.smallEnd(2) > klo_dom || b.bigEnd(2) < khi_dom) { any_inside = true; break; }
    }
    if (!any_inside) { return; }

    AMREX_ALWAYS_ASSERT(rho0.nGrowVect()[0] >= 1 && rho0.nGrowVect()[1] >= 1 && rho0.nGrowVect()[2] >= 1);

    // Coarse time levels and times, as in FillIntermediatePatch
    const TimeLevelSel csel = select_time_level(time_d, t_old[lev-1], t_new[lev-1]);
    const double tref = t_new[lev-1];
    const Real time_rel = static_cast<Real>(time_d - tref);
    Vector<Real> ctime = {static_cast<Real>(t_old[lev-1] - tref), Real(0)};
    Vector<Real> ftime = {time_rel, time_rel};

    ReduceOps<ReduceOpMin, ReduceOpMax> reduce_op;
    ReduceData<Real, Real> reduce_data(reduce_op);
    using ReduceTuple = typename decltype(reduce_data)::Type;

    for (int dir = 0; dir < 2; ++dir)
    {
        MultiFab& mom = (dir == 0) ? rho0u : rho0v;
        const int vel_comp = (dir == 0) ? Vars::xvel : Vars::yvel;
        const int bc_comp  = (dir == 0) ? BCVars::xvel_bc : BCVars::yvel_bc;

        // Coarse velocity with its ghost faces filled, which FillPatchTwoLevels needs
        Vector<MultiFab*> csrc = select_time_data(csel, vars_old[lev-1][vel_comp], vars_new[lev-1][vel_comp]);
        Vector<MultiFab>  ctmp(2);
        Vector<MultiFab*> cmf(2);
        for (int n = 0; n < 2; ++n) {
            if (n == 1 && csrc[1] == csrc[0]) { cmf[1] = cmf[0]; break; }
            ctmp[n].define(csrc[n]->boxArray(), csrc[n]->DistributionMap(), 1, IntVect(1));
            MultiFab::Copy(ctmp[n], *csrc[n], 0, 0, 1, IntVect(1));
            ctmp[n].FillBoundary(geom[lev-1].periodicity());
            cmf[n] = &ctmp[n];
        }

        // Interpolate into the z-ghost layer of a scratch copy of this level's faces
        MultiFab ftmp(mom.boxArray(), mom.DistributionMap(), 1, IntVect(0,0,1));
        ftmp.setVal(0.0);
        Vector<MultiFab*> fmf = {&ftmp, &ftmp};
        FillPatchTwoLevels(ftmp, IntVect(0,0,1), IntVect(0), time_rel, cmf, ctime, fmf, ftime,
                           0, 0, 1, geom[lev-1], geom[lev], refRatio(lev-1),
                           &face_cons_linear_interp, domain_bcs_type, bc_comp);

        for (MFIter mfi(mom); mfi.isValid(); ++mfi)
        {
            const Box& vbx = mfi.validbox();
            const Array4<Real      >& m_arr = mom.array(mfi);
            const Array4<Real const>& f_arr = ftmp.const_array(mfi);
            const Array4<Real const>& r_arr = rho0.const_array(mfi);
            const int ioff = (dir == 0) ? 1 : 0;
            const int joff = (dir == 1) ? 1 : 0;
            for (int side = 0; side < 2; ++side) {
                const int kg = (side == 0) ? vbx.smallEnd(2) - 1 : vbx.bigEnd(2) + 1;
                if (kg < klo_dom || kg > khi_dom) { continue; }
                Box slab(vbx); slab.setRange(2, kg);
                ParallelFor(slab, [=] AMREX_GPU_DEVICE (int i, int j, int k) noexcept
                {
                    const Real r_face = myhalf * (r_arr(i,j,k) + r_arr(i-ioff,j-joff,k));
                    m_arr(i,j,k) = r_face * f_arr(i,j,k);
                });
                reduce_op.eval(slab, reduce_data,
                [=] AMREX_GPU_DEVICE (int i, int j, int k) noexcept -> ReduceTuple
                {
                    const Real r_face = myhalf * (r_arr(i,j,k) + r_arr(i-ioff,j-joff,k));
                    return {r_face, r_face};
                });
            }
        } // mfi
    } // dir

    // The density in the z-ghost cells must be set (a base state, or a filled density)
    ReduceTuple hv = reduce_data.value(reduce_op);
    Real rmin = amrex::get<0>(hv);
    Real rmax = amrex::get<1>(hv);
    ParallelDescriptor::ReduceRealMin(rmin);
    ParallelDescriptor::ReduceRealMax(rmax);
    if (!(rmin > zero) || !(rmax < Real(1.e10))) {
        amrex::Abort("fill_cf_zghost_momenta: the density on the coarse/fine z-ghost faces at level "
                     + std::to_string(lev) + " is not set (min " + std::to_string(rmin) + ", max "
                     + std::to_string(rmax) + ")");
    }
}

/**
 * The coarse level's Omega at time, interpolated onto every w-face of level lev so that the
 * fine faces over each coarse face average to it (face_cons_linear_interp; the ratio in z is
 * one).  Omega is formed in the density the projection uses: rho when the density is fixed, as
 * in anelastic runs, else the base-state rho0.
 *
 * @param lev          Level index (> 0)
 * @param time         Time, for the coarse time interpolation
 * @param omega_target Defined here on the level's w-faces, no ghost cells
 */
void ERF::cf_zface_omega_target (int lev, double time_d, MultiFab& omega_target)
{
    BL_PROFILE("ERF::cf_zface_omega_target()");
    const int levc = lev - 1;
    const TimeLevelSel csel = select_time_level(time_d, t_old[levc], t_new[levc]);
    Real fac_new = one, fac_old = zero;
    if (csel == TimeLevelSel::Old) {
        fac_new = zero; fac_old = one;
    } else if (csel == TimeLevelSel::Interp) {
        fac_new = static_cast<Real>((time_d - t_old[levc]) / (t_new[levc] - t_old[levc]));
        fac_old = one - fac_new;
    }
    auto at_time = [&] (int comp, MultiFab& dst) {
        const MultiFab& mnew = vars_new[levc][comp];
        const MultiFab& mold = vars_old[levc][comp];
        dst.define(mnew.boxArray(), mnew.DistributionMap(), 1, IntVect(1));
        if (csel == TimeLevelSel::New) {
            MultiFab::Copy(dst, mnew, 0, 0, 1, IntVect(1));
        } else if (csel == TimeLevelSel::Old) {
            MultiFab::Copy(dst, mold, 0, 0, 1, IntVect(1));
        } else {
            MultiFab::LinComb(dst, fac_old, mold, 0, fac_new, mnew, 0, 0, 1, IntVect(1));
        }
        dst.FillBoundary(geom[levc].periodicity());
    };
    MultiFab cu, cv, cw, crho;
    at_time(Vars::xvel, cu);
    at_time(Vars::yvel, cv);
    at_time(Vars::zvel, cw);
    if (solverChoice.fixed_density[levc]) {
        at_time(Vars::cons, crho);   // component 0 is rho
    } else {
        crho.define(base_state[levc].boxArray(), base_state[levc].DistributionMap(), 1, IntVect(1));
        MultiFab::Copy(crho, base_state[levc], BaseState::r0_comp, 0, 1, IntVect(1));
    }

    const Box& cdom = geom[levc].Domain();
    const auto cdlo = lbound(cdom);
    const auto cdhi = ubound(cdom);
    const bool per_x = geom[levc].isPeriodic(0);
    const bool per_y = geom[levc].isPeriodic(1);
    auto dxInv_c = geom[levc].InvCellSizeArray();

    MultiFab cOmega(convert(grids[levc], IntVect(0,0,1)), dmap[levc], 1, IntVect(1,1,0));
    cOmega.setVal(0.0);

    for (MFIter mfi(cOmega); mfi.isValid(); ++mfi)
    {
        Box cc = mfi.validbox(); cc.enclosedCells();
        const Box gx = amrex::surroundingNodes(amrex::grow(cc, IntVect(0,0,1)), 0);
        const Box gy = amrex::surroundingNodes(amrex::grow(cc, IntVect(0,0,1)), 1);
        FArrayBox rux(gx, 1, The_Async_Arena());
        FArrayBox rvy(gy, 1, The_Async_Arena());
        const Array4<Real> ru = rux.array();
        const Array4<Real> rv = rvy.array();
        const Array4<Real const> u_a = cu.const_array(mfi);
        const Array4<Real const> v_a = cv.const_array(mfi);
        const Array4<Real const> w_a = cw.const_array(mfi);
        const Array4<Real const> d_a = crho.const_array(mfi);
        ParallelFor(gx, gy,
        [=] AMREX_GPU_DEVICE (int i, int j, int k) noexcept {
            ru(i,j,k) = myhalf * (d_a(i,j,k) + d_a(i-1,j,k)) * u_a(i,j,k);
        },
        [=] AMREX_GPU_DEVICE (int i, int j, int k) noexcept {
            rv(i,j,k) = myhalf * (d_a(i,j,k) + d_a(i,j-1,k)) * v_a(i,j,k);
        });
        const Array4<Real const> ru_c = rux.const_array();
        const Array4<Real const> rv_c = rvy.const_array();
        const Array4<Real> om = cOmega.array(mfi);
        const Array4<Real const> z_nd = z_phys_nd[levc]->const_array(mfi);
        const Array4<Real const> mf_u = mapfac[levc][MapFacType::u_x]->const_array(mfi);
        const Array4<Real const> mf_v = mapfac[levc][MapFacType::v_y]->const_array(mfi);
        ParallelFor(mfi.validbox(), [=] AMREX_GPU_DEVICE (int i, int j, int k) noexcept
        {
            if (k <= cdlo.z || k > cdhi.z) {
                om(i,j,k) = zero;   // ground and top: never a coarse/fine face of a refined level
            } else {
                Real rw = myhalf * (d_a(i,j,k) + d_a(i,j,k-1)) * w_a(i,j,k);
                om(i,j,k) = OmegaFromW(i,j,k,rw,ru_c,rv_c,mf_u,mf_v,z_nd,dxInv_c);
            }
        });
    }
    cOmega.FillBoundary(geom[levc].periodicity());

    // Outside a non-periodic side of the domain, repeat the nearest face, so the interpolation's
    // slopes there stay finite and one-sided
    if (!per_x || !per_y) {
        for (MFIter mfi(cOmega); mfi.isValid(); ++mfi)
        {
            const Array4<Real> om = cOmega.array(mfi);
            const Box gbx = mfi.growntilebox(IntVect(1,1,0));
            ParallelFor(gbx, [=] AMREX_GPU_DEVICE (int i, int j, int k) noexcept
            {
                const int ii = per_x ? i : amrex::min(amrex::max(i, cdlo.x), cdhi.x);
                const int jj = per_y ? j : amrex::min(amrex::max(j, cdlo.y), cdhi.y);
                if (ii != i || jj != j) { om(i,j,k) = om(ii,jj,k); }
            });
        }
    }

    PhysBCFunctNoOp null_bc;
    omega_target.define(convert(grids[lev], IntVect(0,0,1)), dmap[lev], 1, 0);
    InterpFromCoarseLevel(omega_target, IntVect(0), Real(0), cOmega, 0, 0, 1,
                          geom[levc], geom[lev], null_bc, 0, null_bc, 0,
                          refRatio(levc), &face_cons_linear_interp, domain_bcs_type, BCVars::zvel_bc);
}

/**
 * In the projection, after Omega is formed: set Omega on the level's coarse/fine faces in z
 * from the coarse Omega.  The coarse projection leaves no net inflow into the coarse cells under
 * the level, so the level's fluxes then add up to zero over its boundary, and its solve (which
 * applies no correction through these faces) converges to a divergence-free field.  ERF used to
 * keep the coarse (rho w) that FillSet puts on these faces; the fine level then rebuilt Omega
 * from its own u, v and slopes, the fluxes no longer added up to the coarse ones, and the
 * difference stayed in the level as a uniform divergence (2.8e-6 1/s on the test hill).
 *
 * @param lev   Level index (> 0)
 * @param time  Time of the momenta
 * @param omega w-face array that holds Omega
 */
void ERF::set_cf_zface_omega_projection (int lev, double time_d, MultiFab& omega)
{
    BL_PROFILE("ERF::set_cf_zface_omega_projection()");
    if (lev == 0 || !solverChoice.anelastic[lev]) { return; }
    auto patches = cf_zface_patches(grids[lev], geom[lev].Domain());
    if (patches.empty()) { return; }

    MultiFab tgt;
    cf_zface_omega_target(lev, time_d, tgt);
    for (MFIter mfi(omega); mfi.isValid(); ++mfi) {
        auto it = patches.find(mfi.index());
        if (it == patches.end()) { continue; }
        const Array4<Real      > o_a = omega.array(mfi);
        const Array4<Real const> t_a = tgt.const_array(mfi);
        for (const Box& fb : it->second) {
            ParallelFor(fb, [=] AMREX_GPU_DEVICE (int i, int j, int k) noexcept
            { o_a(i,j,k) = t_a(i,j,k); });
        }
    }
}

/**
 * In a fine-level fill, on velocities: give the level's coarse/fine faces in z the w whose Omega
 * is the coarse Omega, w = WFromOmega(Omega_coarse, rho u, rho v) / rho.  Called after the fill
 * has set the coarse/fine faces from the coarse (rho w) (FillSet) and filled every ghost cell,
 * before the momenta are formed from the velocities, so that the time step transports theta and
 * the other scalars with the face fluxes the projection balanced.  Setting them in the
 * projection alone leaves the fill's (rho w) in the advection, and theta in the layer next to the
 * face drifts (2.7 K in 300 s on a 100 m hill).
 *
 * @param lev  Level index
 * @param time Time of the fill
 * @param cons Conserved state with filled ghost cells (rho in component 0)
 * @param xvel x-velocity with filled ghost faces
 * @param yvel y-velocity with filled ghost faces
 * @param zvel z-velocity; its coarse/fine faces in z are set here
 */
void ERF::set_cf_zface_omega_fill (int lev, double time_d, MultiFab& cons,
                                   MultiFab& xvel, MultiFab& yvel, MultiFab& zvel)
{
    BL_PROFILE("ERF::set_cf_zface_omega_fill()");
    if (lev == 0 || !solverChoice.anelastic[lev]) { return; }
    if (solverChoice.mesh_type != MeshType::VariableDz) { return; }
    auto patches = cf_zface_patches(grids[lev], geom[lev].Domain());
    if (patches.empty()) { return; }

    // Momenta in rho scaling, with one ghost layer, from the filled velocities
    MultiFab ru(xvel.boxArray(), xvel.DistributionMap(), 1, IntVect(1,1,1));
    MultiFab rv(yvel.boxArray(), yvel.DistributionMap(), 1, IntVect(1,1,1));
    MultiFab rw(zvel.boxArray(), zvel.DistributionMap(), 1, IntVect(1,1,0));
    VelocityToMomentum(xvel, IntVect(1,1,1), yvel, IntVect(1,1,1), zvel, IntVect(1,1,0), cons,
                       ru, rv, rw, geom[lev].Domain(), domain_bcs_type);

    MultiFab tgt;
    cf_zface_omega_target(lev, time_d, tgt);
    auto dxInv = geom[lev].InvCellSizeArray();

    for (MFIter mfi(zvel); mfi.isValid(); ++mfi)
    {
        auto it = patches.find(mfi.index());
        if (it == patches.end()) { continue; }
        const Array4<Real const> ru_a = ru.const_array(mfi);
        const Array4<Real const> rv_a = rv.const_array(mfi);
        const Array4<Real      > w_a  = zvel.array(mfi);
        const Array4<Real const> r_a  = cons.const_array(mfi);
        const Array4<Real const> t_a  = tgt.const_array(mfi);
        const Array4<Real const> z_nd = z_phys_nd[lev]->const_array(mfi);
        const Array4<Real const> mf_u = mapfac[lev][MapFacType::u_x]->const_array(mfi);
        const Array4<Real const> mf_v = mapfac[lev][MapFacType::v_y]->const_array(mfi);
        for (const Box& fb : it->second) {
            ParallelFor(fb, [=] AMREX_GPU_DEVICE (int i, int j, int k) noexcept
            {
                const Real rw_new = WFromOmega(i,j,k,t_a(i,j,k),ru_a,rv_a,mf_u,mf_v,z_nd,dxInv);
                w_a(i,j,k) = rw_new / (myhalf * (r_a(i,j,k,Rho_comp) + r_a(i,j,k-1,Rho_comp)));
            });
        }
    }
}
