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
#include "ERF_TI_slow_headers.H"

#include <AMReX_FillPatchUtil.H>
#include <AMReX_Interpolater.H>

#include <limits>

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
 * The coarse Omega on the coarse/fine w-faces of level lev, and the coarse (rho u), (rho v) on the
 * layer just outside them, on the coarse planes of cf_zface_layout, at coarse time tc (linear in
 * time between the coarse old and new states).  Not yet clamped outside the domain.
 *
 * @param lev           Level index (> 0)
 * @param tc            Time
 * @param conserved_rho Scaling: conserved rho (true) or base-state rho0 (false)
 * @param pz            Coarse z_phys_nd on the planes
 * @param pmu, pmv      Coarse map factors on the planes
 * @param cOm           Coarse Omega (output)
 * @param cRu, cRv      Coarse (rho u), (rho v) on the outside layer (output)
 */
void ERF::cf_zface_coarse_at (int lev, double tc, bool conserved_rho,
                              const MultiFab& pz, const MultiFab& pmu, const MultiFab& pmv,
                              MultiFab& cOm, MultiFab& cRu, MultiFab& cRv)
{
    const int levc = lev - 1;
    const CFZFaceLayout& L = cf_zface_layout(lev);
    const DistributionMapping& sdm = L.sdm;
    const Periodicity cper = geom[levc].periodicity();
    const Box& cdom = geom[levc].Domain();
    const int ck_lo = cdom.smallEnd(2);
    const int ck_hi = cdom.bigEnd(2);
    auto dxInv_c = geom[levc].InvCellSizeArray();
    const bool per_x = geom[levc].isPeriodic(0);
    const bool per_y = geom[levc].isPeriodic(1);

    // ---------------------------------------------------------------------------------
    // Coarse data on the planes, at tc
    // ---------------------------------------------------------------------------------
    const TimeLevelSel csel = select_time_level(tc, t_old[levc], t_new[levc]);
    Real fac_new = one, fac_old = zero;
    if (csel == TimeLevelSel::Old) {
        fac_new = zero; fac_old = one;
    } else if (csel == TimeLevelSel::Interp) {
        fac_new = static_cast<Real>((tc - t_old[levc]) / (t_new[levc] - t_old[levc]));
        fac_old = one - fac_new;
    }
    // Component comp of the coarse data at tc on the planes (valid coarse cells only; cells
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
    clamp_outside_domain(pr, cdom, per_x, per_y);

    // ---------------------------------------------------------------------------------
    // Coarse (rho u), (rho v) on the outside layer and Omega on the faces
    // ---------------------------------------------------------------------------------
    cOm.define(L.ba_cf, sdm, 1, 0);
    cRu.define(L.ba_cx, sdm, 1, 0);
    cRv.define(L.ba_cy, sdm, 1, 0);

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

    const Periodicity cper = geom[levc].periodicity();
    const Box& cdom = geom[levc].Domain();
    const int ck_lo = cdom.smallEnd(2);
    const int ck_hi = cdom.bigEnd(2);
    const bool per_x = geom[levc].isPeriodic(0);
    const bool per_y = geom[levc].isPeriodic(1);

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

    MultiFab cOm, cRu, cRv;
    cf_zface_coarse_at(lev, time_d, conserved_rho, pz, pmu, pmv, cOm, cRu, cRv);

    // In the first stage of a fine step, shift Omega so that over the coarse step the fine faces
    // apply the Omega of the coarse step's stages (see cf_anelastic_stage_registration)
    if (cf_anelastic_stage1(lev)) {
        MultiFab cOm_o, cOm_n, tu, tv;
        cf_zface_coarse_at(lev, t_old[levc], conserved_rho, pz, pmu, pmv, cOm_o, tu, tv);
        cf_zface_coarse_at(lev, t_new[levc], conserved_rho, pz, pmu, pmv, cOm_n, tu, tv);
        MultiFab rOm(L.ba_cf, sdm, 1, 0);
        rOm.setVal(0.0);
        rOm.ParallelCopy(m_cf_step_zflux[levc], 0, 0, 1, IntVect(0), IntVect(0), cper);
        const Real kfac = (solverChoice.anelastic_type[lev] == AnelasticType::MidPoint) ? one : Real(2);
        for (MFIter mfi(cOm); mfi.isValid(); ++mfi) {
            const Array4<Real      > om = cOm.array(mfi);
            const Array4<Real const> oo = cOm_o.const_array(mfi);
            const Array4<Real const> on = cOm_n.const_array(mfi);
            const Array4<Real const> ro = rOm.const_array(mfi);
            ParallelFor(cOm[mfi].box(), [=] AMREX_GPU_DEVICE (int i, int j, int k) noexcept
            {
                if (k > ck_lo && k <= ck_hi) {
                    om(i,j,k) += kfac * (ro(i,j,k) - myhalf * (oo(i,j,k) + on(i,j,k)));
                }
            });
        }
    }

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
 * Whether the coarse/fine faces of this level take the horizontal mass flux from the coarse level:
 * an anelastic level on a terrain-fitted mesh.  The face areas vary along the coarse/fine faces
 * there, so interpolating (rho u) would leave the fine fluxes not adding up to the coarse flux.
 * The mass flux is the one the advection forms, ax (rho u) / m_uy and ay (rho v) / m_vx, so the
 * map factors, which can vary along the faces as well, are included.
 *
 * @param lev Level index (> 0)
 */
bool ERF::cf_interp_mass_flux (int lev) const
{
    return solverChoice.anelastic[lev] && solverChoice.mesh_type == MeshType::VariableDz;
}

/**
 * The mass flux a (rho u) / m of a horizontal momentum component, on the faces of its valid
 * region and the ghost faces that the momentum, the area and the map factor all hold.
 *
 * @param m    Face momenta, (rho u) or (rho v)
 * @param area Face areas, ax or ay
 * @param mf   Map factor on those faces, m_uy or m_vx
 */
static MultiFab
weighted_mass_flux (const MultiFab& m, const MultiFab& area, const MultiFab& mf)
{
    IntVect ng = amrex::min(m.nGrowVect(), area.nGrowVect());
    ng[0] = amrex::min(ng[0], mf.nGrowVect()[0]);
    ng[1] = amrex::min(ng[1], mf.nGrowVect()[1]);
    MultiFab w(m.boxArray(), m.DistributionMap(), 1, ng);
    for (MFIter mfi(w, TilingIfNotGPU()); mfi.isValid(); ++mfi) {
        const Box bx = mfi.growntilebox(ng);
        const auto w_arr = w.array(mfi);
        const auto m_arr = m.const_array(mfi);
        const auto a_arr = area.const_array(mfi);
        const auto f_arr = mf.const_array(mfi);
        ParallelFor(bx, [=] AMREX_GPU_DEVICE (int i, int j, int k) noexcept
        {
            w_arr(i,j,k) = m_arr(i,j,k) * a_arr(i,j,k) / f_arr(i,j,0);
        });
    }
    return w;
}

/**
 * Register the coarse horizontal momenta for the coarse/fine face fills of the level above,
 * as ax (rho u) / m_uy and ay (rho v) / m_vx where cf_interp_mass_flux holds.
 *
 * @param levc  Coarse level index
 * @param xold  x-face momenta at the first time
 * @param xnew  x-face momenta at the second time
 * @param yold  y-face momenta at the first time
 * @param ynew  y-face momenta at the second time
 * @param times The two times
 */
void ERF::RegisterCoarseMomentumCF (int levc, const MultiFab& xold, const MultiFab& xnew,
                                    const MultiFab& yold, const MultiFab& ynew,
                                    const Vector<double>& times)
{
    if (!cf_interp_mass_flux(levc+1)) {
        FPr_u[levc].RegisterCoarseData({&xold,&xnew}, times);
        FPr_v[levc].RegisterCoarseData({&yold,&ynew}, times);
        return;
    }
    const MultiFab& mf_uy = *mapfac[levc][MapFacType::u_y];
    const MultiFab& mf_vx = *mapfac[levc][MapFacType::v_x];
    const MultiFab xo = weighted_mass_flux(xold, *ax[levc], mf_uy);
    const MultiFab xn = weighted_mass_flux(xnew, *ax[levc], mf_uy);
    const MultiFab yo = weighted_mass_flux(yold, *ay[levc], mf_vx);
    const MultiFab yn = weighted_mass_flux(ynew, *ay[levc], mf_vx);
    FPr_u[levc].RegisterCoarseData({&xo,&xn}, times);
    FPr_v[levc].RegisterCoarseData({&yo,&yn}, times);
}

/**
 * Average the momenta of level crse_lev+1 down onto the faces of level crse_lev it covers, as
 * mass fluxes (the vertical as Omega, see below).  Horizontally: ax (rho u) / m_uy and ay (rho v) / m_vx are averaged, and divided by the
 * coarse area and map factor.  Over terrain the areas vary along a coarse face, so the plain
 * average of (rho u) is not the coarse momentum whose flux is the mean of the fine fluxes: on a
 * coarse/fine face it would replace the projected coarse momentum with another one, and the
 * uncovered coarse cell next to it would start the next step with a divergence.  The fine
 * coarse/fine faces end each step on the coarse mass flux (cf_anelastic_stage_registration), so
 * here the coarse/fine faces get back the flux they had.  Faces the fine level does not cover
 * are left as they are.
 *
 * @param crse_lev Coarse level index
 */
void ERF::average_down_mass_flux (int crse_lev)
{
    const int lev = crse_lev + 1;
    const IntVect rr = refRatio(crse_lev);
    MultiFab* mom_f[2]       = {&rU_new[lev], &rV_new[lev]};
    MultiFab* mom_c[2]       = {&rU_new[crse_lev], &rV_new[crse_lev]};
    const MultiFab* area_f[2] = {ax[lev].get(), ay[lev].get()};
    const MultiFab* area_c[2] = {ax[crse_lev].get(), ay[crse_lev].get()};
    const MultiFab* mf_f[2]   = {mapfac[lev][MapFacType::u_y].get(), mapfac[lev][MapFacType::v_x].get()};
    const MultiFab* mf_c[2]   = {mapfac[crse_lev][MapFacType::u_y].get(), mapfac[crse_lev][MapFacType::v_x].get()};

    // The coarse cells the fine level covers, with a layer in z
    const Periodicity fper = geom[lev].periodicity();
    const Periodicity cper = geom[crse_lev].periodicity();
    const iMultiFab covered = amrex::makeFineMask(grids[crse_lev], dmap[crse_lev], IntVect(0,0,1),
                                                  grids[lev], rr, cper, 0, 1);

    // The coarse Omega on the fine level's coarse/fine faces in z, before the horizontal momenta
    //    around them change: those faces keep their Omega (see the vertical below)
    const auto dxInv_c = geom[crse_lev].InvCellSizeArray();
    const int ck_lo = geom[crse_lev].Domain().smallEnd(2);
    const int ck_hi = geom[crse_lev].Domain().bigEnd(2);
    rU_new[crse_lev].FillBoundary(cper);
    rV_new[crse_lev].FillBoundary(cper);
    MultiFab om_cf(rW_new[crse_lev].boxArray(), rW_new[crse_lev].DistributionMap(), 1, 0);
    for (MFIter mfi(om_cf, TilingIfNotGPU()); mfi.isValid(); ++mfi) {
        const Box bx = mfi.tilebox();
        const auto o_arr = om_cf.array(mfi);
        const auto c_arr = covered.const_array(mfi);
        const auto u_arr = rU_new[crse_lev].const_array(mfi);
        const auto v_arr = rV_new[crse_lev].const_array(mfi);
        const auto w_arr = rW_new[crse_lev].const_array(mfi);
        const auto z_nd  = z_phys_nd[crse_lev]->const_array(mfi);
        const auto mf_u  = mapfac[crse_lev][MapFacType::u_x]->const_array(mfi);
        const auto mf_v  = mapfac[crse_lev][MapFacType::v_y]->const_array(mfi);
        ParallelFor(bx, [=] AMREX_GPU_DEVICE (int i, int j, int k) noexcept
        {
            const bool cf = (c_arr(i,j,k-1) != c_arr(i,j,k)) && k > ck_lo && k <= ck_hi;
            o_arr(i,j,k) = cf ? OmegaFromW(i,j,k,w_arr(i,j,k),u_arr,v_arr,mf_u,mf_v,z_nd,dxInv_c)
                              : zero;
        });
    }

    for (int d = 0; d < 2; ++d) {
        const MultiFab wf = weighted_mass_flux(*mom_f[d], *area_f[d], *mf_f[d]);

        // The coarse faces under the fine level, on the fine level's distribution
        const BoxArray cba = amrex::coarsen(mom_f[d]->boxArray(), rr);
        const DistributionMapping& dm = mom_f[d]->DistributionMap();
        MultiFab wc(cba, dm, 1, 0);
        amrex::average_down_faces(wf, wc, rr);

        // The coarse area and map factor on those faces (the map factor has no extent in z)
        MultiFab ac(cba, dm, 1, 0);
        ac.ParallelCopy(*area_c[d], 0, 0, 1);
        BoxList bl2d = cba.boxList();
        for (auto& b : bl2d) { b.setRange(2, 0); }
        MultiFab mfc(BoxArray(std::move(bl2d)), dm, 1, 0);
        mfc.ParallelCopy(*mf_c[d], 0, 0, 1);

        for (MFIter mfi(wc, TilingIfNotGPU()); mfi.isValid(); ++mfi) {
            const Box bx = mfi.tilebox();
            const auto w_arr = wc.array(mfi);
            const auto a_arr = ac.const_array(mfi);
            const auto f_arr = mfc.const_array(mfi);
            ParallelFor(bx, [=] AMREX_GPU_DEVICE (int i, int j, int k) noexcept
            {
                w_arr(i,j,k) = w_arr(i,j,k) * f_arr(i,j,0) / a_arr(i,j,k);
            });
        }
        mom_c[d]->ParallelCopy(wc, 0, 0, 1, IntVect(0), IntVect(0), geom[crse_lev].periodicity());
    }

    // ---------------------------------------------------------------------------------------
    // The vertical, as the flux az Omega / (m_x m_y).  Averaging (rho w) leaves the covered
    // coarse cells with a divergence over terrain, and their scalars, advected with it, reach
    // the uncovered cells through the advection stencil.  Each coarse face with covered cells
    // on both sides takes the mean of the fine fluxes over it, and its (rho w) is rebuilt from
    // that Omega and the coarse (rho u), (rho v) averaged down above; the fine level's
    // coarse/fine faces in z keep their coarse Omega, which the fine faces there reproduce on
    // average (cf_zface_coarse_data): their (rho w) is rebuilt from it, since the horizontal
    // momenta it depends on over terrain have changed.  All other coarse faces are left as they are.
    // ---------------------------------------------------------------------------------------
    rU_new[lev].FillBoundary(fper);
    rV_new[lev].FillBoundary(fper);
    rU_new[crse_lev].FillBoundary(cper);
    rV_new[crse_lev].FillBoundary(cper);

    // The fine fluxes, on the faces with fine cells on both sides (a coarse/fine face in z reads
    //    momenta outside the level, and its coarse face is not set from it below)
    iMultiFab in_fine(grids[lev], dmap[lev], 1, IntVect(0,0,1));
    in_fine.setVal(0);
    in_fine.setVal(1, 0);
    in_fine.FillBoundary(fper);
    const int fk_lo = geom[lev].Domain().smallEnd(2);
    const int fk_hi = geom[lev].Domain().bigEnd(2);
    const auto dxInv_f = geom[lev].InvCellSizeArray();
    MultiFab ff(rW_new[lev].boxArray(), rW_new[lev].DistributionMap(), 1, 0);
    for (MFIter mfi(ff, TilingIfNotGPU()); mfi.isValid(); ++mfi) {
        const Box bx = mfi.tilebox();
        const auto f_arr = ff.array(mfi);
        const auto m_arr = in_fine.const_array(mfi);
        const auto u_arr = rU_new[lev].const_array(mfi);
        const auto v_arr = rV_new[lev].const_array(mfi);
        const auto w_arr = rW_new[lev].const_array(mfi);
        const auto z_nd  = z_phys_nd[lev]->const_array(mfi);
        const auto mf_u  = mapfac[lev][MapFacType::u_x]->const_array(mfi);
        const auto mf_v  = mapfac[lev][MapFacType::v_y]->const_array(mfi);
        const auto mx    = mapfac[lev][MapFacType::m_x]->const_array(mfi);
        const auto my    = mapfac[lev][MapFacType::m_y]->const_array(mfi);
        const auto a_arr = az[lev]->const_array(mfi);
        ParallelFor(bx, [=] AMREX_GPU_DEVICE (int i, int j, int k) noexcept
        {
            if (k == fk_lo || k == fk_hi + 1 || m_arr(i,j,k-1) == 0 || m_arr(i,j,k) == 0) {
                f_arr(i,j,k) = zero;
            } else {
                const Real om = OmegaFromW(i,j,k,w_arr(i,j,k),u_arr,v_arr,mf_u,mf_v,z_nd,dxInv_f);
                f_arr(i,j,k) = a_arr(i,j,k) * om / (mx(i,j,0) * my(i,j,0));
            }
        });
    }

    // Their means over the coarse faces, on the coarse level's faces under the fine level
    const BoxArray cba_w = amrex::coarsen(ff.boxArray(), rr);
    MultiFab fc(cba_w, ff.DistributionMap(), 1, 0);
    amrex::average_down_faces(ff, fc, rr);
    MultiFab fcc(rW_new[crse_lev].boxArray(), rW_new[crse_lev].DistributionMap(), 1, 0);
    fcc.setVal(zero);
    fcc.ParallelCopy(fc, 0, 0, 1, IntVect(0), IntVect(0), cper);

    for (MFIter mfi(rW_new[crse_lev], TilingIfNotGPU()); mfi.isValid(); ++mfi) {
        const Box bx = mfi.tilebox();
        const auto w_arr = rW_new[crse_lev].array(mfi);
        const auto f_arr = fcc.const_array(mfi);
        const auto o_arr = om_cf.const_array(mfi);
        const auto c_arr = covered.const_array(mfi);
        const auto u_arr = rU_new[crse_lev].const_array(mfi);
        const auto v_arr = rV_new[crse_lev].const_array(mfi);
        const auto z_nd  = z_phys_nd[crse_lev]->const_array(mfi);
        const auto mf_u  = mapfac[crse_lev][MapFacType::u_x]->const_array(mfi);
        const auto mf_v  = mapfac[crse_lev][MapFacType::v_y]->const_array(mfi);
        const auto mx    = mapfac[crse_lev][MapFacType::m_x]->const_array(mfi);
        const auto my    = mapfac[crse_lev][MapFacType::m_y]->const_array(mfi);
        const auto a_arr = az[crse_lev]->const_array(mfi);
        ParallelFor(bx, [=] AMREX_GPU_DEVICE (int i, int j, int k) noexcept
        {
            if (c_arr(i,j,k-1) == 1 && c_arr(i,j,k) == 1) {
                const Real om = f_arr(i,j,k) * mx(i,j,0) * my(i,j,0) / a_arr(i,j,k);
                w_arr(i,j,k) = WFromOmega(i,j,k,om,u_arr,v_arr,mf_u,mf_v,z_nd,dxInv_c);
            } else if (c_arr(i,j,k-1) != c_arr(i,j,k) && k > ck_lo && k <= ck_hi) {
                w_arr(i,j,k) = WFromOmega(i,j,k,o_arr(i,j,k),u_arr,v_arr,mf_u,mf_v,z_nd,dxInv_c);
            }
        });
    }
}

//
// Anelastic coarse/fine faces in time.
//
// Over a coarse step the coarse level advects with the mass flux R its stages apply:
// (F_1 + F_2)/2 with the RK2 (Heun) integrator, F_2 with MidPoint, where F_s is the flux of
// the momenta at the start of stage s (avg_xmom and avg_ymom, set by erf_slow_rhs_pre).  The
// fine level's coarse/fine faces are filled by interpolating linearly in time between the
// coarse momenta at the start and end of the step, w^n and w^{n+1} (in mass-flux form), and
// over the fine steps of a coarse step its faces apply (w^n + w^{n+1})/2 -- the trapezoid and
// midpoint rules are exact for a linear function -- which is not R.  Scalars are refluxed and
// rho is not, so the difference changes a uniform scalar next to the fine level (#4873).
//
// So the coarse level records R on its faces (cf_anelastic_record), and in the first stage of
// each fine step the fine level fills its coarse/fine faces from the coarse momenta shifted by
// k D, D = R - (w^n + w^{n+1})/2, k = 2 for RK2 and 1 for MidPoint
// (cf_anelastic_stage_registration).  The faces that stage fills are those of the projection
// that ends it, whose momenta the second stage advects with; the first stage advects with the
// momenta the step started with, filled at its start.  A fine RK2 step then applies
// (I(t_a) + I(t_b) + 2D)/2 and a MidPoint step I(t_a + dt/2) + D, I the linear interpolation in
// time, and over the coarse step these add up to R.  The spatial interpolation keeps the mean of
// the fine faces under each coarse face, so this holds face by face.  The second stage fills
// from the unshifted momenta, so the step ends, and the next one starts, on the linear
// interpolation.
//
// erf.cf_flux_history = false turns this off, as it does the compressible flux history.
//

/**
 * Whether the coarse/fine faces of level lev apply, over a coarse step, the mass flux the coarse
 * level's stages applied: both levels anelastic, on a terrain-fitted mesh (cf_interp_mass_flux)
 *
 * @param lev Fine level index
 */
bool ERF::cf_anelastic_flux_match (int lev) const
{
    return cf_flux_history_enabled() && lev > 0 && lev <= finest_level &&
           solverChoice.anelastic[lev] && solverChoice.anelastic[lev-1] &&
           solverChoice.terrain_type != TerrainType::EB &&
           solverChoice.terrain_type != TerrainType::MovingFittedMesh;
}

/**
 * Add stage nrk's horizontal mass fluxes (avg_xmom, avg_ymom) to the flux R the coarse level
 * levc applies over its step.  Called once erf_slow_rhs_pre has set them.
 *
 * @param levc Coarse level index
 * @param nrk  The stage
 */
void ERF::cf_anelastic_record (int levc, int nrk)
{
    if (!cf_anelastic_flux_match(levc+1)) { return; }
    if (static_cast<int>(m_cf_step_flux.size()) <= levc) { m_cf_step_flux.resize(levc+1); }
    const Real b = StageFluxWeight(solverChoice, levc, nrk);
    const MultiFab* avg[2] = {&avg_xmom[levc], &avg_ymom[levc]};
    for (int d = 0; d < 2; ++d) {
        MultiFab& R = m_cf_step_flux[levc][d];
        if (R.boxArray() != avg[d]->boxArray() || R.DistributionMap() != avg[d]->DistributionMap()) {
            R.define(avg[d]->boxArray(), avg[d]->DistributionMap(), 1, 0);
        }
        if (nrk == 0) { R.setVal(zero); }
        MultiFab::Saxpy(R, b, *avg[d], 0, 0, 1, 0);
    }

    // The vertical, as Omega (avg_zmom is az Omega / (m_x m_y)); Omega is (rho w) without terrain
    if (static_cast<int>(m_cf_step_zflux.size()) <= levc) { m_cf_step_zflux.resize(levc+1); }
    MultiFab& Rz = m_cf_step_zflux[levc];
    if (Rz.boxArray() != avg_zmom[levc].boxArray() || Rz.DistributionMap() != avg_zmom[levc].DistributionMap()) {
        Rz.define(avg_zmom[levc].boxArray(), avg_zmom[levc].DistributionMap(), 1, 0);
    }
    if (nrk == 0) { Rz.setVal(zero); }
    for (MFIter mfi(Rz, TilingIfNotGPU()); mfi.isValid(); ++mfi) {
        const Box bx = mfi.tilebox();
        const auto r_arr = Rz.array(mfi);
        const auto f_arr = avg_zmom[levc].const_array(mfi);
        const auto a_arr = az[levc]->const_array(mfi);
        const auto mx = mapfac[levc][MapFacType::m_x]->const_array(mfi);
        const auto my = mapfac[levc][MapFacType::m_y]->const_array(mfi);
        ParallelFor(bx, [=] AMREX_GPU_DEVICE (int i, int j, int k) noexcept
        {
            r_arr(i,j,k) += b * f_arr(i,j,k) * mx(i,j,0) * my(i,j,0) / a_arr(i,j,k);
        });
    }
}

/**
 * Whether fine level lev is in the first stage of its step with its coarse/fine faces taking the
 * coarse step's flux (cf_anelastic_stage_registration): its fills then shift the coarse Omega on
 * its coarse/fine faces in z as they shift the horizontal momenta (cf_zface_coarse_data)
 *
 * @param lev Fine level index
 */
bool ERF::cf_anelastic_stage1 (int lev) const
{
    return lev < static_cast<int>(m_cf_stage1.size()) && m_cf_stage1[lev] != 0 &&
           cf_anelastic_flux_match(lev) && cf_interp_mass_flux(lev) &&
           lev-1 < static_cast<int>(m_cf_step_zflux.size()) && m_cf_step_zflux[lev-1].ok();
}

/**
 * Register the coarse horizontal momenta for the coarse/fine faces of fine level lev in stage
 * nrk of its step: shifted so that the fine faces apply the coarse step's flux in the first
 * stage, as they are in the second (see above)
 *
 * @param lev Fine level index
 * @param nrk The stage
 */
void ERF::cf_anelastic_stage_registration (int lev, int nrk)
{
    if (!cf_anelastic_flux_match(lev)) { return; }
    const int levc = lev - 1;
    if (static_cast<int>(m_cf_step_flux.size()) <= levc ||
        m_cf_step_flux[levc][0].boxArray() != avg_xmom[levc].boxArray() ||
        FPr_u[levc].crse_times().size() != 2) {
        Abort("cf_anelastic_stage_registration: the coarse level has not recorded its step");
    }

    if (static_cast<int>(m_cf_stage1.size()) <= lev) { m_cf_stage1.resize(lev+1, 0); }
    m_cf_stage1[lev] = (nrk == 0) ? 1 : 0;

    const Real k = (solverChoice.anelastic_type[lev] == AnelasticType::MidPoint) ? one : Real(2);
    // On a terrain-fitted mesh the coarse momenta are registered as mass fluxes
    //    (RegisterCoarseMomentumCF); on a flat mesh as momenta, and the shift is converted
    const bool mass_flux = cf_interp_mass_flux(lev);
    const MultiFab* mold[2] = {&rU_old[levc], &rV_old[levc]};
    const MultiFab* mnew[2] = {&rU_new[levc], &rV_new[levc]};
    const MultiFab* area[2] = {ax[levc].get(), ay[levc].get()};
    const MultiFab* mf[2]   = {mapfac[levc][MapFacType::u_y].get(), mapfac[levc][MapFacType::v_x].get()};
    ERFFillPatcher* fpr[2]  = {&FPr_u[levc], &FPr_v[levc]};

    for (int d = 0; d < 2; ++d) {
        if (!mass_flux && nrk != 0) {
            fpr[d]->RegisterCoarseData({mold[d], mnew[d]}, fpr[d]->crse_times());
            continue;
        }
        MultiFab wo = weighted_mass_flux(*mold[d], *area[d], *mf[d]);
        MultiFab wn = weighted_mass_flux(*mnew[d], *area[d], *mf[d]);
        if (nrk == 0) {
            // k D on the coarse level's faces, zero on the ghost faces outside the domain;
            //    as a momentum, k D m / a, where the momenta are registered as such
            const IntVect ng = wo.nGrowVect();
            MultiFab shift(wo.boxArray(), wo.DistributionMap(), 1, ng);
            shift.setVal(zero);
            const MultiFab& R = m_cf_step_flux[levc][d];
            for (MFIter mfi(shift, TilingIfNotGPU()); mfi.isValid(); ++mfi) {
                const Box bx = mfi.tilebox();
                const auto s_arr  = shift.array(mfi);
                const auto r_arr  = R.const_array(mfi);
                const auto wo_arr = wo.const_array(mfi);
                const auto wn_arr = wn.const_array(mfi);
                const auto a_arr  = area[d]->const_array(mfi);
                const auto f_arr  = mf[d]->const_array(mfi);
                ParallelFor(bx, [=] AMREX_GPU_DEVICE (int i, int j, int kk) noexcept
                {
                    const Real D = k * (r_arr(i,j,kk) - myhalf * (wo_arr(i,j,kk) + wn_arr(i,j,kk)));
                    s_arr(i,j,kk) = mass_flux ? D : D * f_arr(i,j,0) / a_arr(i,j,kk);
                });
            }
            shift.FillBoundary(geom[levc].periodicity());
            if (!mass_flux) {
                MultiFab::Copy(wo, *mold[d], 0, 0, 1, ng);
                MultiFab::Copy(wn, *mnew[d], 0, 0, 1, ng);
            }
            MultiFab::Add(wo, shift, 0, 0, 1, ng);
            MultiFab::Add(wn, shift, 0, 0, 1, ng);
        }
        // the times the coarse momenta were registered with (RegisterCoarseMomentumCF)
        fpr[d]->RegisterCoarseData({&wo,&wn}, fpr[d]->crse_times());
    }

    // The coarse/fine faces in z of a level that stops inside the domain in z.  On a flat mesh
    //    they are filled from the coarse (rho w), which is Omega there: shift it the same way.
    //    On a terrain-fitted mesh they take the coarse Omega instead (cf_zface_coarse_data).
    if (!mass_flux) {
        const MultiFab& wold = rW_old[levc];
        const MultiFab& wnew = rW_new[levc];
        if (nrk != 0) {
            FPr_w[levc].RegisterCoarseData({&wold, &wnew}, FPr_w[levc].crse_times());
        } else {
            const IntVect ng = wold.nGrowVect();
            MultiFab zo(wold.boxArray(), wold.DistributionMap(), 1, ng);
            MultiFab zn(wold.boxArray(), wold.DistributionMap(), 1, ng);
            MultiFab shift(wold.boxArray(), wold.DistributionMap(), 1, ng);
            shift.setVal(zero);
            const MultiFab& Rz = m_cf_step_zflux[levc];
            for (MFIter mfi(shift, TilingIfNotGPU()); mfi.isValid(); ++mfi) {
                const Box bx = mfi.tilebox();
                const auto s_arr  = shift.array(mfi);
                const auto r_arr  = Rz.const_array(mfi);
                const auto wo_arr = wold.const_array(mfi);
                const auto wn_arr = wnew.const_array(mfi);
                ParallelFor(bx, [=] AMREX_GPU_DEVICE (int i, int j, int kk) noexcept
                {
                    s_arr(i,j,kk) = k * (r_arr(i,j,kk) - myhalf * (wo_arr(i,j,kk) + wn_arr(i,j,kk)));
                });
            }
            shift.FillBoundary(geom[levc].periodicity());
            MultiFab::Copy(zo, wold, 0, 0, 1, ng);
            MultiFab::Copy(zn, wnew, 0, 0, 1, ng);
            MultiFab::Add(zo, shift, 0, 0, 1, ng);
            MultiFab::Add(zn, shift, 0, 0, 1, ng);
            FPr_w[levc].RegisterCoarseData({&zo, &zn}, FPr_w[levc].crse_times());
        }
    }
}

/**
 * Fill the coarse/fine faces of one horizontal momentum component.  Where cf_interp_mass_flux
 * holds, the fill gives the mass flux, which is multiplied by the fine map factor and divided by
 * the fine area on the faces it set; all other faces keep their values exactly.
 */
static void
fill_cf_momentum (ERFFillPatcher& patcher, MultiFab& mom, const MultiFab& area, const MultiFab& mapfac,
                  bool mass_flux, double time_d, const Vector<BCRec>& bcs)
{
    PhysBCFunctNoOp null_bc;
    if (!mass_flux) {
        patcher.FillSet(mom, time_d, null_bc, bcs);
        return;
    }
    constexpr Real not_set = std::numeric_limits<Real>::max();
    MultiFab flux(mom.boxArray(), mom.DistributionMap(), 1, 0);
    flux.setVal(not_set);
    patcher.FillSet(flux, time_d, null_bc, bcs);
    for (MFIter mfi(flux, TilingIfNotGPU()); mfi.isValid(); ++mfi) {
        const Box bx = mfi.tilebox();
        const auto f = flux.const_array(mfi);
        const auto a = area.const_array(mfi);
        const auto mf = mapfac.const_array(mfi);
        const auto m = mom.array(mfi);
        const Real flag = not_set;
        ParallelFor(bx, [=] AMREX_GPU_DEVICE (int i, int j, int k) noexcept
        {
            if (f(i,j,k) != flag) { m(i,j,k) = f(i,j,k) * mf(i,j,0) / a(i,j,k); }
        });
    }
}

/**
 * Fill the horizontal momenta on the coarse/fine faces of a level (see RegisterCoarseMomentumCF).
 *
 * @param lev    Level index (> 0)
 * @param time_d Time of the fill
 * @param xmom   x-face momenta, (rho u)
 * @param ymom   y-face momenta, (rho v)
 */
void ERF::FillSetHorizontalMomentaCF (int lev, double time_d, MultiFab& xmom, MultiFab& ymom)
{
    const bool mass_flux = cf_interp_mass_flux(lev);
    fill_cf_momentum(FPr_u[lev-1], xmom, *ax[lev], *mapfac[lev][MapFacType::u_y], mass_flux, time_d,
                     domain_bcs_type);
    fill_cf_momentum(FPr_v[lev-1], ymom, *ay[lev], *mapfac[lev][MapFacType::v_x], mass_flux, time_d,
                     domain_bcs_type);
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
    FillSetHorizontalMomentaCF(lev, time_d, xmom, ymom);
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
