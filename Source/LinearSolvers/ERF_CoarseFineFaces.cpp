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
 *
 * cf_zface_coarse_data builds both from the coarse level, on thin planes around the level's top
 * and bottom only: the horizontal momenta in the layer just outside the level, and the coarse
 * Omega on its coarse/fine faces in z.  project_momenta uses them through
 * set_cf_zghost_projection and set_cf_zface_omega; every fine-level fill uses them through
 * FillSetMomentaCF, right after the coarse/fine faces are set from the coarse level and before
 * any ghost fill, so the time step carries theta and the other scalars with the fluxes the
 * projection balanced.
 *
 * ERF used to keep the coarse (rho w) that FillSet puts on these faces: the fine level then
 * rebuilt Omega from its own u, v and slopes, the fluxes no longer added up to the coarse ones,
 * and the difference stayed in the level as a uniform divergence (2.8e-6 1/s on the test hill).
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
 */
std::map<int, Vector<Box>>
make_cf_zface_patches (const BoxArray& fba, const Box& fdom)
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

/**
 * Outside a non-periodic side of the domain, repeat the nearest value inside it.  The coarse
 * planes extend past the domain as valid cells, so the coarse patch the interpolation copies
 * from them sees finite, one-sided data there.
 */
void
clamp_outside_domain (MultiFab& mf, const Box& cdom, bool per_x, bool per_y)
{
    if (per_x && per_y) { return; }
    const IndexType ixt = mf.ixType();
    const int xlo = cdom.smallEnd(0);
    const int ylo = cdom.smallEnd(1);
    const int xhi = cdom.bigEnd(0) + (ixt.nodeCentered(0) ? 1 : 0);
    const int yhi = cdom.bigEnd(1) + (ixt.nodeCentered(1) ? 1 : 0);
    for (MFIter mfi(mf); mfi.isValid(); ++mfi) {
        const Array4<Real> a = mf.array(mfi);
        const Box bx = mfi.validbox();
        const int ilo = bx.smallEnd(0), ihi = bx.bigEnd(0);
        const int jlo = bx.smallEnd(1), jhi = bx.bigEnd(1);
        ParallelFor(bx, [=] AMREX_GPU_DEVICE (int i, int j, int k) noexcept
        {
            int ii = per_x ? i : amrex::min(amrex::max(i, xlo), xhi);
            int jj = per_y ? j : amrex::min(amrex::max(j, ylo), yhi);
            ii = amrex::min(amrex::max(ii, ilo), ihi);
            jj = amrex::min(amrex::max(jj, jlo), jhi);
            if (ii != i || jj != j) { a(i,j,k) = a(ii,jj,k); }
        });
    }
}

/**
 * Copy fine slabs (one per coarse/fine face of a box, with the boxes' distribution) into the
 * level's face array: slab s goes into box slab_box[s].
 */
void
copy_slabs_into (MultiFab& dst, const MultiFab& slabs, const Vector<int>& slab_box)
{
    for (MFIter mfi(slabs); mfi.isValid(); ++mfi) {
        const Array4<Real      > d_a = dst.array(slab_box[mfi.index()]);
        const Array4<Real const> s_a = slabs.const_array(mfi);
        ParallelFor(mfi.validbox(), [=] AMREX_GPU_DEVICE (int i, int j, int k) noexcept
        { d_a(i,j,k) = s_a(i,j,k); });
    }
}

} // namespace

/**
 * The layout of level lev's coarse/fine z-face data, kept until the level's grids or
 * distribution change: its coarse/fine w-faces per box, and one slab per box face whose outside
 * layer is inside the domain (slab s belongs to box slab_box[s], and the slabs have the boxes'
 * distribution), with the coarse planes under each slab.  Every fill reuses these BoxArrays and
 * this DistributionMapping, and so the copy patterns AMReX caches for them.
 *
 * @param lev Level index (> 0)
 */
const ERF::CFZFaceLayout&
ERF::cf_zface_layout (int lev)
{
    CFZFaceLayout& L = m_cf_zface_cache[lev];
    if (!L.fba.empty() && L.fba == grids[lev] && L.fdm == dmap[lev]) { return L; }

    L = CFZFaceLayout{};
    L.fba = grids[lev];
    L.fdm = dmap[lev];
    const BoxArray& fba = grids[lev];
    const Box& fdom = geom[lev].Domain();
    L.patches = make_cf_zface_patches(fba, fdom);

    const IntVect rr = refRatio(lev-1);
    BoxList bl_f (IndexType(IntVect(0,0,1))), bl_x (IndexType(IntVect(1,0,0))), bl_y (IndexType(IntVect(0,1,0)));
    BoxList bl_cf(IndexType(IntVect(0,0,1))), bl_cx(IndexType(IntVect(1,0,0))), bl_cy(IndexType(IntVect(0,1,0)));
    BoxList bl_r;
    Vector<int> pmap;
    for (int n = 0; n < fba.size(); ++n) {
        const Box& b = fba[n];
        for (int side = 0; side < 2; ++side) {
            const int kout = (side == 0) ? b.smallEnd(2) - 1 : b.bigEnd(2) + 1;
            const int kf   = (side == 0) ? b.smallEnd(2)     : b.bigEnd(2) + 1;
            if (kout < fdom.smallEnd(2) || kout > fdom.bigEnd(2)) { continue; }
            Box f = amrex::surroundingNodes(b, 2); f.setRange(2, kf);
            Box x = amrex::surroundingNodes(b, 0); x.setRange(2, kout);
            Box y = amrex::surroundingNodes(b, 1); y.setRange(2, kout);
            bl_f.push_back(f); bl_x.push_back(x); bl_y.push_back(y);
            // The coarse cells under the slab in the two layers around the face (ci, with room
            // for the interpolation), and around them room for the metric stencil of OmegaFromW
            Box c = amrex::coarsen(b, rr); c.setRange(2, kf-1, 2);
            const Box ci = amrex::grow(c, IntVect(1,1,0));
            bl_r.push_back(amrex::grow(c, IntVect(3,3,0)));
            Box cf = amrex::surroundingNodes(ci, 2); cf.setRange(2, kf);
            Box cx = amrex::surroundingNodes(ci, 0); cx.setRange(2, kout);
            Box cy = amrex::surroundingNodes(ci, 1); cy.setRange(2, kout);
            bl_cf.push_back(cf); bl_cx.push_back(cx); bl_cy.push_back(cy);
            pmap.push_back(dmap[lev][n]);
            L.slab_box.push_back(n);
        }
    }
    if (L.slab_box.empty()) { return L; }

    L.sdm   = DistributionMapping(std::move(pmap));
    L.ba_r  = BoxArray(std::move(bl_r));
    L.ba_ru = amrex::convert(L.ba_r, IntVect(1,0,0));
    L.ba_rv = amrex::convert(L.ba_r, IntVect(0,1,0));
    L.ba_rw = amrex::convert(L.ba_r, IntVect(0,0,1));
    L.ba_rz = amrex::convert(L.ba_r, IntVect(1,1,1));
    BoxList bl_2d = L.ba_r.boxList();
    for (auto& bb : bl_2d) { bb.setRange(2, 0); }
    const BoxArray ba_2d(std::move(bl_2d));
    L.ba_mu = amrex::convert(ba_2d, IntVect(1,0,0));
    L.ba_mv = amrex::convert(ba_2d, IntVect(0,1,0));
    L.ba_cf = BoxArray(std::move(bl_cf));
    L.ba_cx = BoxArray(std::move(bl_cx));
    L.ba_cy = BoxArray(std::move(bl_cy));
    L.ba_f  = BoxArray(std::move(bl_f));
    L.ba_x  = BoxArray(std::move(bl_x));
    L.ba_y  = BoxArray(std::move(bl_y));
    return L;
}

/**
 * Build, from the coarse level at time, the data a refined level needs at its coarse/fine faces
 * in z, on one slab per such face of each of its boxes (slab s belongs to box slab_box[s], and
 * the slabs have the boxes' distribution):
 *
 *   rhou_g, rhov_g: (rho u) and (rho v) on the x- and y-faces of the layer just outside the box
 *                   (below its bottom, above its top);
 *   omega_t:        the coarse Omega on the box's top or bottom w-face.
 *
 * All three are interpolated with face_cons_linear_interp, so the fine faces over each coarse
 * face average to it (the ratio in z is one).  The density is the conserved rho when
 * conserved_rho is true (the level's momenta are rho u, as with a fixed density), else the base
 * state's rho0, so the result is in the scaling of the level's own valid faces and needs no
 * density from the fine level.  Only coarse planes around these faces are copied and computed on.
 *
 * @param lev           Level index (> 0)
 * @param time          Time, for the coarse time interpolation
 * @param conserved_rho Scaling: conserved rho (true) or base-state rho0 (false)
 * @param omega_t       Fine slabs on the w-faces (output)
 * @param rhou_g        Fine slabs on the x-faces of the outside layer (output)
 * @param rhov_g        Fine slabs on the y-faces of the outside layer (output)
 * @param slab_box      Box of level lev each slab belongs to (output)
 * @return              False when the level has no coarse/fine faces in z (outputs not defined)
 */
bool ERF::cf_zface_coarse_data (int lev, double time_d, bool conserved_rho,
                                MultiFab& omega_t, MultiFab& rhou_g, MultiFab& rhov_g,
                                Vector<int>& slab_box)
{
    BL_PROFILE("ERF::cf_zface_coarse_data()");
    AMREX_ALWAYS_ASSERT(lev > 0);
    const int levc = lev - 1;
    const IntVect rr = refRatio(levc);
    const CFZFaceLayout& L = cf_zface_layout(lev);
    AMREX_ALWAYS_ASSERT(rr[2] == 1 || L.patches.empty());

    slab_box = L.slab_box;
    if (slab_box.empty()) { return false; }
    const DistributionMapping& sdm = L.sdm;

    // ---------------------------------------------------------------------------------
    // Coarse data on the planes, at time
    // ---------------------------------------------------------------------------------
    const TimeLevelSel csel = select_time_level(time_d, t_old[levc], t_new[levc]);
    Real fac_new = one, fac_old = zero;
    if (csel == TimeLevelSel::Old) {
        fac_new = zero; fac_old = one;
    } else if (csel == TimeLevelSel::Interp) {
        fac_new = static_cast<Real>((time_d - t_old[levc]) / (t_new[levc] - t_old[levc]));
        fac_old = one - fac_new;
    }
    const Periodicity cper = geom[levc].periodicity();
    // Component comp of the coarse data at time on the planes (valid coarse cells only; cells
    // past a non-periodic side stay zero and are clamped after the computation)
    auto plane_at_time = [&] (const MultiFab& mold, const MultiFab& mnew, int comp, MultiFab& dst) {
        dst.setVal(0.0);
        if (csel != TimeLevelSel::Old) {
            dst.ParallelCopy(mnew, comp, 0, 1, IntVect(0), IntVect(0), cper);
        }
        if (csel != TimeLevelSel::New) {
            MultiFab tmp(dst.boxArray(), dst.DistributionMap(), 1, 0);
            tmp.setVal(0.0);
            tmp.ParallelCopy(mold, comp, 0, 1, IntVect(0), IntVect(0), cper);
            if (csel == TimeLevelSel::Old) {
                MultiFab::Copy(dst, tmp, 0, 0, 1, 0);
            } else {
                MultiFab::LinComb(dst, fac_old, tmp, 0, fac_new, dst, 0, 0, 1, 0);
            }
        }
    };

    MultiFab pu(L.ba_ru, sdm, 1, 0);
    MultiFab pv(L.ba_rv, sdm, 1, 0);
    MultiFab pw(L.ba_rw, sdm, 1, 0);
    MultiFab pr(L.ba_r , sdm, 1, 0);
    plane_at_time(vars_old[levc][Vars::xvel], vars_new[levc][Vars::xvel], 0, pu);
    plane_at_time(vars_old[levc][Vars::yvel], vars_new[levc][Vars::yvel], 0, pv);
    plane_at_time(vars_old[levc][Vars::zvel], vars_new[levc][Vars::zvel], 0, pw);
    if (conserved_rho) {
        plane_at_time(vars_old[levc][Vars::cons], vars_new[levc][Vars::cons], Rho_comp, pr);
    } else {
        plane_at_time(base_state[levc], base_state[levc], BaseState::r0_comp, pr);
    }
    clamp_outside_domain(pr, geom[levc].Domain(), geom[levc].isPeriodic(0), geom[levc].isPeriodic(1));

    // The mesh and the map factors, their ghosts outside the domain included
    MultiFab pz(L.ba_rz, sdm, 1, 0);
    pz.setVal(0.0);
    pz.ParallelCopy(*z_phys_nd[levc], 0, 0, 1, z_phys_nd[levc]->nGrowVect(), IntVect(0), cper);
    MultiFab pmu(L.ba_mu, sdm, 1, 0);
    MultiFab pmv(L.ba_mv, sdm, 1, 0);
    pmu.setVal(1.0);
    pmv.setVal(1.0);
    pmu.ParallelCopy(*mapfac[levc][MapFacType::u_x], 0, 0, 1, mapfac[levc][MapFacType::u_x]->nGrowVect(), IntVect(0), cper);
    pmv.ParallelCopy(*mapfac[levc][MapFacType::v_y], 0, 0, 1, mapfac[levc][MapFacType::v_y]->nGrowVect(), IntVect(0), cper);

    // ---------------------------------------------------------------------------------
    // Coarse (rho u), (rho v) on the outside layer and Omega on the faces
    // ---------------------------------------------------------------------------------
    MultiFab cOm(L.ba_cf, sdm, 1, 0);
    MultiFab cRu(L.ba_cx, sdm, 1, 0);
    MultiFab cRv(L.ba_cy, sdm, 1, 0);
    const Box& cdom = geom[levc].Domain();
    const int ck_lo = cdom.smallEnd(2);
    const int ck_hi = cdom.bigEnd(2);
    auto dxInv_c = geom[levc].InvCellSizeArray();

    for (MFIter mfi(pr); mfi.isValid(); ++mfi)
    {
        const Box rbx = mfi.validbox();
        const Box gx = amrex::grow(amrex::surroundingNodes(rbx, 0), IntVect(-1,0,0));
        const Box gy = amrex::grow(amrex::surroundingNodes(rbx, 1), IntVect(0,-1,0));
        FArrayBox rux(gx, 1, The_Async_Arena());
        FArrayBox rvy(gy, 1, The_Async_Arena());
        const Array4<Real> ru = rux.array();
        const Array4<Real> rv = rvy.array();
        const Array4<Real const> u_a = pu.const_array(mfi);
        const Array4<Real const> v_a = pv.const_array(mfi);
        const Array4<Real const> w_a = pw.const_array(mfi);
        const Array4<Real const> d_a = pr.const_array(mfi);
        // Faces with a plane cell on each side (the planes' outermost faces are not needed)
        ParallelFor(gx, gy,
        [=] AMREX_GPU_DEVICE (int i, int j, int k) noexcept {
            ru(i,j,k) = myhalf * (d_a(i,j,k) + d_a(i-1,j,k)) * u_a(i,j,k);
        },
        [=] AMREX_GPU_DEVICE (int i, int j, int k) noexcept {
            rv(i,j,k) = myhalf * (d_a(i,j,k) + d_a(i,j-1,k)) * v_a(i,j,k);
        });

        const Array4<Real> om = cOm.array(mfi);
        const Array4<Real> cu = cRu.array(mfi);
        const Array4<Real> cv = cRv.array(mfi);
        const Array4<Real const> ru_c = rux.const_array();
        const Array4<Real const> rv_c = rvy.const_array();
        const Array4<Real const> z_nd = pz.const_array(mfi);
        const Array4<Real const> mf_u = pmu.const_array(mfi);
        const Array4<Real const> mf_v = pmv.const_array(mfi);
        ParallelFor(cOm[mfi].box(), [=] AMREX_GPU_DEVICE (int i, int j, int k) noexcept
        {
            if (k <= ck_lo || k > ck_hi) {
                om(i,j,k) = zero;   // ground and top: never a coarse/fine face of a refined level
            } else {
                Real rw = myhalf * (d_a(i,j,k) + d_a(i,j,k-1)) * w_a(i,j,k);
                om(i,j,k) = OmegaFromW(i,j,k,rw,ru_c,rv_c,mf_u,mf_v,z_nd,dxInv_c);
            }
        });
        ParallelFor(cRu[mfi].box(), cRv[mfi].box(),
        [=] AMREX_GPU_DEVICE (int i, int j, int k) noexcept { cu(i,j,k) = ru_c(i,j,k); },
        [=] AMREX_GPU_DEVICE (int i, int j, int k) noexcept { cv(i,j,k) = rv_c(i,j,k); });
    }
    const bool per_x = geom[levc].isPeriodic(0);
    const bool per_y = geom[levc].isPeriodic(1);
    clamp_outside_domain(cOm, cdom, per_x, per_y);
    clamp_outside_domain(cRu, cdom, per_x, per_y);
    clamp_outside_domain(cRv, cdom, per_x, per_y);

    // ---------------------------------------------------------------------------------
    // Interpolate to the fine slabs
    // ---------------------------------------------------------------------------------
    PhysBCFunctNoOp null_bc;
    omega_t.define(L.ba_f, sdm, 1, 0);
    rhou_g.define(L.ba_x, sdm, 1, 0);
    rhov_g.define(L.ba_y, sdm, 1, 0);
    InterpFromCoarseLevel(omega_t, IntVect(0), Real(0), cOm, 0, 0, 1, geom[levc], geom[lev],
                          null_bc, 0, null_bc, 0, rr, &face_cons_linear_interp, domain_bcs_type, BCVars::zvel_bc);
    InterpFromCoarseLevel(rhou_g, IntVect(0), Real(0), cRu, 0, 0, 1, geom[levc], geom[lev],
                          null_bc, 0, null_bc, 0, rr, &face_cons_linear_interp, domain_bcs_type, BCVars::xvel_bc);
    InterpFromCoarseLevel(rhov_g, IntVect(0), Real(0), cRv, 0, 0, 1, geom[levc], geom[lev],
                          null_bc, 0, null_bc, 0, rr, &face_cons_linear_interp, domain_bcs_type, BCVars::yvel_bc);
    return true;
}

/**
 * In the projection: fill the layer of momenta just below and above the level from the coarse
 * level (OmegaFromW and WFromOmega read it; FillBoundary fills it only where another fine box
 * covers it, so this comes first), and return the coarse Omega for the level's coarse/fine faces
 * in z, which set_cf_zface_omega sets once Omega is formed.  The scaling is the projected one:
 * rho when the density is fixed, else rho0.
 *
 * @param lev      Level index (> 0)
 * @param time     Time of the momenta
 * @param xmom     x-face momenta of the level (its outside layer is written)
 * @param ymom     y-face momenta of the level (its outside layer is written)
 * @param omega_t  Coarse Omega on the level's coarse/fine w-faces (output)
 * @param slab_box Box each slab of omega_t belongs to (output; empty if the level has none)
 */
void ERF::set_cf_zghost_projection (int lev, double time_d, MultiFab& xmom, MultiFab& ymom,
                                    MultiFab& omega_t, Vector<int>& slab_box)
{
    BL_PROFILE("ERF::set_cf_zghost_projection()");
    MultiFab rhou_g, rhov_g;
    if (!cf_zface_coarse_data(lev, time_d, solverChoice.fixed_density[lev],
                              omega_t, rhou_g, rhov_g, slab_box)) {
        return;
    }
    AMREX_ALWAYS_ASSERT(xmom.nGrowVect()[2] >= 1 && ymom.nGrowVect()[2] >= 1);
    copy_slabs_into(xmom, rhou_g, slab_box);
    copy_slabs_into(ymom, rhov_g, slab_box);
}

/**
 * In the projection, after Omega is formed: set Omega on the level's coarse/fine faces in z from
 * the coarse Omega (anelastic levels, whose fills do the same).  The coarse projection leaves no
 * net inflow into the coarse cells under the level, so the level's fluxes then add up to zero over
 * its boundary, and its solve, which applies no correction through these faces, converges to a
 * divergence-free field.
 *
 * @param lev      Level index (> 0)
 * @param omega    w-face array that holds Omega
 * @param omega_t  Coarse Omega from set_cf_zghost_projection
 * @param slab_box Box each slab of omega_t belongs to
 */
void ERF::set_cf_zface_omega (int lev, MultiFab& omega, const MultiFab& omega_t,
                              const Vector<int>& slab_box)
{
    BL_PROFILE("ERF::set_cf_zface_omega()");
    if (slab_box.empty() || !solverChoice.anelastic[lev]) { return; }
    const auto& patches = cf_zface_layout(lev).patches;
    for (MFIter mfi(omega_t); mfi.isValid(); ++mfi) {
        const int n = slab_box[mfi.index()];
        auto it = patches.find(n);
        if (it == patches.end()) { continue; }
        const Array4<Real      > o_a = omega.array(n);
        const Array4<Real const> t_a = omega_t.const_array(mfi);
        const Box& slab = mfi.validbox();
        for (const Box& fb : it->second) {
            const Box bx = fb & slab;
            if (!bx.ok()) { continue; }
            ParallelFor(bx, [=] AMREX_GPU_DEVICE (int i, int j, int k) noexcept
            { o_a(i,j,k) = t_a(i,j,k); });
        }
    }
}

/**
 * Set a refined level's momenta on its coarse/fine faces from the coarse level (FillSet), and on
 * an anelastic level over terrain give its faces in z the (rho w) whose Omega is the coarse
 * Omega,  rho w = WFromOmega(Omega_coarse, rho u, rho v),  with the layer of (rho u), (rho v)
 * outside the level taken from the coarse level too.  Every fine-level fill of the momenta goes
 * through here, before any ghost fill or FillBoundary, so every copy of these faces (the box's own
 * and its neighbours' ghost copies) gets the same value, and the time step transports theta and
 * the other scalars with the fluxes the projection balanced.  Setting them in the projection
 * alone left the fill's (rho w) in the advection, and theta in the layer next to the face drifted
 * (2.7 K in 300 s on a 100 m hill).
 *
 * @param lev  Level index (> 0)
 * @param time Time of the fill
 * @param xmom x-face momenta, (rho u)
 * @param ymom y-face momenta, (rho v)
 * @param zmom z-face momenta, (rho w)
 */
void ERF::FillSetMomentaCF (int lev, double time_d, MultiFab& xmom, MultiFab& ymom, MultiFab& zmom)
{
    BL_PROFILE("ERF::FillSetMomentaCF()");
    PhysBCFunctNoOp null_bc;
    FPr_u[lev-1].FillSet(xmom, time_d, null_bc, domain_bcs_type);
    FPr_v[lev-1].FillSet(ymom, time_d, null_bc, domain_bcs_type);
    FPr_w[lev-1].FillSet(zmom, time_d, null_bc, domain_bcs_type);

    if (!solverChoice.anelastic[lev] || solverChoice.mesh_type != MeshType::VariableDz) { return; }

    MultiFab omega_t, rhou_g, rhov_g;
    Vector<int> slab_box;
    if (!cf_zface_coarse_data(lev, time_d, true, omega_t, rhou_g, rhov_g, slab_box)) { return; }

    const auto& patches = cf_zface_layout(lev).patches;
    const BoxArray& fba = grids[lev];
    auto dxInv = geom[lev].InvCellSizeArray();
    for (MFIter mfi(omega_t); mfi.isValid(); ++mfi) {
        const int n = slab_box[mfi.index()];
        auto it = patches.find(n);
        if (it == patches.end()) { continue; }
        const Box& slab = mfi.validbox();

        // WFromOmega reads (rho u), (rho v) in the two layers around the face: the level's own
        // inside it, and outside it the coarse layer from rhou_g, rhov_g.  Put both in a scratch
        // slab, so the level's momenta are only read here and need no ghost layer in z.
        const int kf   = slab.smallEnd(2);
        const int kout = (kf == fba[n].smallEnd(2)) ? kf - 1 : kf;
        Box box_u = rhou_g[mfi].box(); box_u.setRange(2, kf-1, 2);
        Box box_v = rhov_g[mfi].box(); box_v.setRange(2, kf-1, 2);
        FArrayBox ru_fab(box_u, 1, The_Async_Arena());
        FArrayBox rv_fab(box_v, 1, The_Async_Arena());
        const Array4<Real      > ru_s = ru_fab.array();
        const Array4<Real      > rv_s = rv_fab.array();
        const Array4<Real const> xm_a = xmom.const_array(n);
        const Array4<Real const> ym_a = ymom.const_array(n);
        const Array4<Real const> xg_a = rhou_g.const_array(mfi);
        const Array4<Real const> yg_a = rhov_g.const_array(mfi);
        ParallelFor(box_u, box_v,
        [=] AMREX_GPU_DEVICE (int i, int j, int k) noexcept {
            ru_s(i,j,k) = (k == kout) ? xg_a(i,j,k) : xm_a(i,j,k);
        },
        [=] AMREX_GPU_DEVICE (int i, int j, int k) noexcept {
            rv_s(i,j,k) = (k == kout) ? yg_a(i,j,k) : ym_a(i,j,k);
        });

        const Array4<Real      > w_a  = zmom.array(n);
        const Array4<Real const> ru_a = ru_fab.const_array();
        const Array4<Real const> rv_a = rv_fab.const_array();
        const Array4<Real const> t_a  = omega_t.const_array(mfi);
        const Array4<Real const> z_nd = z_phys_nd[lev]->const_array(n);
        const Array4<Real const> mf_u = mapfac[lev][MapFacType::u_x]->const_array(n);
        const Array4<Real const> mf_v = mapfac[lev][MapFacType::v_y]->const_array(n);
        for (const Box& fb : it->second) {
            const Box bx = fb & slab;
            if (!bx.ok()) { continue; }
            ParallelFor(bx, [=] AMREX_GPU_DEVICE (int i, int j, int k) noexcept
            {
                w_a(i,j,k) = WFromOmega(i,j,k,t_a(i,j,k),ru_a,rv_a,mf_u,mf_v,z_nd,dxInv);
            });
        }
    }
}
