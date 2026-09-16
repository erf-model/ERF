#include <AMReX_MultiFab.H>
#include <AMReX_GpuContainers.H>
#include <AMReX_Reduce.H>

#include <ERF_Diffusion.H>

#include <gtest/gtest.h>

#include <algorithm>
#include <cmath>
#include <limits>

using namespace amrex;

namespace {

// Column: theta = a + b z + c z^2 and qv = qa + qb z at the cell centres z = (k + 1/2) dz, rho = 1,
// no TKE, and a constant vertical eddy diffusivity rho K for theta and qv (no molecular diffusion,
// no strain, no dissipation, no surface layer).  The theta diffusion then writes the face flux
// -K (theta_k - theta_{k-1}) / dz = -K (b + 2 c z_face) exactly, and the buoyancy production of k in
// cell k is |g| / theta_ref times the average of the fluxes at its two faces, -K (b + 2 c z_k),
// whether the vertical diffusion is explicit or implicit.
constexpr Real dz      = Real(0.5);
constexpr Real K       = Real(0.8);
constexpr Real a       = Real(300.0);
constexpr Real b       = Real(0.004);
constexpr Real qa      = Real(0.01);
constexpr Real qb      = Real(-0.0003);
constexpr Real gabs    = Real(9.81);
constexpr Real theta0  = Real(300.0);

struct BuoyancyResult {
    Real max_hfx_error = 0.0;   // hfx_z against -K (b + 2 c z_face) over interior faces
    Real max_qfx_error = 0.0;   // qfx1_z against -K qb over interior faces
    Real max_src_error = 0.0;   // RhoKE rhs against |g|/theta0 (-K (b + 2 c z_k)) over interior cells
    Real flux_scale    = 0.0;
    Real src_scale     = 0.0;
};

Real max_abs_error (const MultiFab& mf, int comp, const Box& region, int kind, Real c)
{
    ReduceOps<ReduceOpMax> reduce_op;
    ReduceData<Real> reduce_data(reduce_op);
    using ReduceTuple = typename decltype(reduce_data)::Type;
    for (MFIter mfi(mf); mfi.isValid(); ++mfi) {
        const Box bx = region & mfi.fabbox();
        if (!bx.ok()) { continue; }
        const auto arr = mf.const_array(mfi);
        reduce_op.eval(bx, reduce_data, [=] AMREX_GPU_DEVICE (int i, int j, int k) -> ReduceTuple
        {
            Real expected = Real(0.0);
            if (kind == 0) {        // theta face flux at z_face = k dz
                expected = -K * (b + Real(2.0) * c * static_cast<Real>(k) * dz);
            } else if (kind == 1) { // qv face flux
                expected = -K * qb;
            } else {                // buoyancy production at z_k = (k + 1/2) dz
                expected = gabs / theta0 * (-K * (b + Real(2.0) * c * (static_cast<Real>(k) + Real(0.5)) * dz));
            }
            return std::abs(arr(i,j,k,comp) - expected);
        });
    }
    return get<0>(reduce_data.value(reduce_op));
}

void fill_column (MultiFab& cell_data, MultiFab& cell_prim, MultiFab& mu_turb, Real c)
{
    for (MFIter mfi(cell_data); mfi.isValid(); ++mfi) {
        auto data = cell_data.array(mfi);
        auto prim = cell_prim.array(mfi);
        auto mu   = mu_turb.array(mfi);
        ParallelFor(mfi.fabbox(), [=] AMREX_GPU_DEVICE (int i, int j, int k) noexcept
        {
            const Real z = (static_cast<Real>(k) + Real(0.5)) * dz;
            const Real theta = a + b * z + c * z * z;
            const Real qv = qa + qb * z;
            data(i,j,k,Rho_comp)      = Real(1.0);
            data(i,j,k,RhoTheta_comp) = theta;
            data(i,j,k,RhoKE_comp)    = Real(0.0);
            data(i,j,k,RhoQ1_comp)    = qv;
            prim(i,j,k,PrimTheta_comp) = theta;
            prim(i,j,k,PrimKE_comp)    = Real(0.0);
            prim(i,j,k,PrimQ1_comp)    = qv;
            mu(i,j,k,EddyDiff::Theta_v) = K;
            mu(i,j,k,EddyDiff::Q_v)     = K;
        });
    }
    Gpu::streamSynchronize();
}

BuoyancyResult run_column (Real implicit_fac, Real c)
{
    const Box domain(IntVect(0, 0, 0), IntVect(3, 3, 11));
    const BoxArray ba(domain);
    const DistributionMapping dm(ba);
    const Box domain_2d(IntVect(0, 0, 0), IntVect(3, 3, 0));
    const BoxArray ba_2d(domain_2d);

    MultiFab cell_data(ba, dm, NVAR_max, 2);
    MultiFab cell_prim(ba, dm, NPRIMVAR_max, 2);
    MultiFab cell_rhs(ba, dm, NVAR_max, 0);
    MultiFab xvel(BoxArray(surroundingNodes(domain, 0)), dm, 1, 1);
    MultiFab yvel(BoxArray(surroundingNodes(domain, 1)), dm, 1, 1);
    MultiFab xflux(BoxArray(surroundingNodes(domain, 0)), dm, 1, 0);
    MultiFab yflux(BoxArray(surroundingNodes(domain, 1)), dm, 1, 0);
    MultiFab zflux(BoxArray(surroundingNodes(domain, 2)), dm, 1, 0);
    MultiFab hfx_x(BoxArray(surroundingNodes(domain, 0)), dm, 1, 0);
    MultiFab hfx_y(BoxArray(surroundingNodes(domain, 1)), dm, 1, 0);
    MultiFab hfx_z(BoxArray(surroundingNodes(domain, 2)), dm, 1, 1);
    MultiFab qfx1_x(BoxArray(surroundingNodes(domain, 0)), dm, 1, 0);
    MultiFab qfx1_y(BoxArray(surroundingNodes(domain, 1)), dm, 1, 0);
    MultiFab qfx1_z(BoxArray(surroundingNodes(domain, 2)), dm, 1, 1);
    MultiFab qfx2_z(BoxArray(surroundingNodes(domain, 2)), dm, 1, 0);
    MultiFab diss(ba, dm, 1, 1);
    MultiFab smn(ba, dm, 1, 1);
    MultiFab mu_turb(ba, dm, EddyDiff::NumDiffs, 2);
    MultiFab mf_mx(ba_2d, dm, 1, 3);
    MultiFab mf_ux(BoxArray(surroundingNodes(domain_2d, 0)), dm, 1, 3);
    MultiFab mf_vx(BoxArray(surroundingNodes(domain_2d, 1)), dm, 1, 3);
    MultiFab mf_my(ba_2d, dm, 1, 3);
    MultiFab mf_uy(BoxArray(surroundingNodes(domain_2d, 0)), dm, 1, 3);
    MultiFab mf_vy(BoxArray(surroundingNodes(domain_2d, 1)), dm, 1, 3);

    cell_data.setVal(Real(0.0));
    cell_prim.setVal(Real(0.0));
    cell_rhs.setVal(Real(0.0));
    xvel.setVal(Real(0.0));
    yvel.setVal(Real(0.0));
    hfx_z.setVal(Real(0.0));
    qfx1_z.setVal(Real(0.0));
    diss.setVal(Real(0.0));
    smn.setVal(Real(0.0));
    mu_turb.setVal(Real(0.0));
    mf_mx.setVal(Real(1.0));
    mf_ux.setVal(Real(1.0));
    mf_vx.setVal(Real(1.0));
    mf_my.setVal(Real(1.0));
    mf_uy.setVal(Real(1.0));
    mf_vy.setVal(Real(1.0));
    fill_column(cell_data, cell_prim, mu_turb, c);

    Vector<BCRec> bcs(NBCVAR_max);
    for (auto& bc : bcs) {
        for (int dir = 0; dir < AMREX_SPACEDIM; ++dir) {
            bc.setLo(dir, ERFBCType::foextrap);
            bc.setHi(dir, ERFBCType::foextrap);
        }
    }
    Gpu::DeviceVector<BCRec> bcs_d(bcs.size());
    Gpu::copy(Gpu::hostToDevice, bcs.begin(), bcs.end(), bcs_d.begin());

    SolverChoice solver_choice;
    solver_choice.diffChoice.molec_diff_type = MolecDiffType::None;
    solver_choice.turbChoice.resize(1);
    solver_choice.turbChoice[0].use_kturb = true;
    solver_choice.turbChoice[0].use_keqn  = true;
    solver_choice.turbChoice[0].les_type  = LESType::Deardorff;
    solver_choice.turbChoice[0].theta_ref = theta0;
    solver_choice.turbChoice[0].implicit_tke_dissipation = false;

    const GpuArray<Real, AMREX_SPACEDIM> dx_inv = {Real(1.0), Real(1.0), Real(1.0)/dz};
    const GpuArray<Real, AMREX_SPACEDIM> grav = {Real(0.0), Real(0.0), -gabs};
    Array4<const Real> tm_arr{};
    auto hfx_x_arr = hfx_x[0].array();
    auto hfx_y_arr = hfx_y[0].array();
    auto hfx_z_arr = hfx_z[0].array();
    auto qfx1_x_arr = qfx1_x[0].array();
    auto qfx1_y_arr = qfx1_y[0].array();
    auto qfx1_z_arr = qfx1_z[0].array();
    auto qfx2_arr = qfx2_z[0].array();
    auto diss_arr = diss[0].array();
    Vector<std::unique_ptr<SurfaceLayer>> surf_layer(6);

    // theta and k in one call, as in erf_slow_rhs_pre / post; then qv
    DiffusionSrcForState_N(
        domain, domain, RhoTheta_comp, 2,
        xvel[0].const_array(), yvel[0].const_array(), cell_data[0].const_array(), cell_prim[0].const_array(),
        cell_rhs[0].array(), xflux[0].array(), yflux[0].array(), zflux[0].array(), dx_inv,
        smn[0].const_array(), mf_mx[0].const_array(), mf_ux[0].const_array(), mf_vx[0].const_array(),
        mf_my[0].const_array(), mf_uy[0].const_array(), mf_vy[0].const_array(), hfx_x_arr, hfx_y_arr, hfx_z_arr,
        qfx1_x_arr, qfx1_y_arr, qfx1_z_arr, qfx2_arr, diss_arr, mu_turb[0].const_array(), solver_choice, 0,
        tm_arr, grav, bcs_d.data(), false, surf_layer, implicit_fac);
    DiffusionSrcForState_N(
        domain, domain, RhoQ1_comp, 1,
        xvel[0].const_array(), yvel[0].const_array(), cell_data[0].const_array(), cell_prim[0].const_array(),
        cell_rhs[0].array(), xflux[0].array(), yflux[0].array(), zflux[0].array(), dx_inv,
        smn[0].const_array(), mf_mx[0].const_array(), mf_ux[0].const_array(), mf_vx[0].const_array(),
        mf_my[0].const_array(), mf_uy[0].const_array(), mf_vy[0].const_array(), hfx_x_arr, hfx_y_arr, hfx_z_arr,
        qfx1_x_arr, qfx1_y_arr, qfx1_z_arr, qfx2_arr, diss_arr, mu_turb[0].const_array(), solver_choice, 0,
        tm_arr, grav, bcs_d.data(), false, surf_layer, implicit_fac);
    Gpu::streamSynchronize();

    // Interior faces k = 1..11 and interior cells k = 1..10 (away from the one-sided boundary faces)
    const Box faces(IntVect(0, 0, 1), IntVect(3, 3, 11), IntVect(0, 0, 1));
    const Box cells(IntVect(0, 0, 1), IntVect(3, 3, 10));

    BuoyancyResult r;
    r.max_hfx_error = max_abs_error(hfx_z,  0, faces, 0, c);
    r.max_qfx_error = max_abs_error(qfx1_z, 0, faces, 1, c);
    r.max_src_error = max_abs_error(cell_rhs, RhoKE_comp, cells, 2, c);
    const Real zmax = Real(12.0) * dz;
    r.flux_scale = K * (std::abs(b) + Real(2.0) * std::abs(c) * zmax + std::abs(qb)) + K * a / dz;
    r.src_scale  = gabs / theta0 * r.flux_scale;
    return r;
}

Real tol (Real scale)
{
    return Real(128.0) * std::numeric_limits<Real>::epsilon() * std::max(Real(1.0), scale);
}

// The buoyancy production of k reads the flux of the theta diffusion at the two faces of each cell.
// With the implicit vertical solve (implicit_fac = 1) the stored face fluxes must still be the full
// fluxes, not the explicit fraction (which is zero), and the source must not change.
TEST(TKEBuoyancySource, FaceFluxesAreFullForExplicitAndImplicitDiffusion)
{
    for (Real fac : {Real(0.0), Real(0.5), Real(1.0)}) {
        const auto r = run_column(fac, Real(0.0));
        EXPECT_LE(r.max_hfx_error, tol(r.flux_scale)) << "implicit_fac = " << fac;
        EXPECT_LE(r.max_qfx_error, tol(r.flux_scale)) << "implicit_fac = " << fac;
        EXPECT_LE(r.max_src_error, tol(r.src_scale))  << "implicit_fac = " << fac;
    }
}

// With theta quadratic in z the flux varies linearly, so the average of the two face fluxes is the
// exact cell-centred flux; reading the lower face alone is off by K c dz.
TEST(TKEBuoyancySource, CellCentredFluxIsTheAverageOfTheTwoFaces)
{
    const Real c = Real(0.1);
    // The lower-face error the check has to be able to see, also in single precision
    ASSERT_GT(gabs / theta0 * K * c * dz, Real(4.0) * tol(gabs / theta0 * K * a / dz));
    for (Real fac : {Real(0.0), Real(1.0)}) {
        const auto r = run_column(fac, c);
        EXPECT_LE(r.max_hfx_error, tol(r.flux_scale)) << "implicit_fac = " << fac;
        EXPECT_LE(r.max_src_error, tol(r.src_scale))  << "implicit_fac = " << fac;
    }
}

} // namespace
