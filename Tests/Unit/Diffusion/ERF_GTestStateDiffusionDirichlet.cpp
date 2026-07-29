#include <AMReX_MultiFab.H>

#include <ERF_Diffusion.H>

#include <gtest/gtest.h>

#include <array>

using namespace amrex;

namespace {

struct DiffusionCase {
    Real max_rhs = 0.0;
    Real low_flux = 0.0;
    Real high_flux = 0.0;
};

DiffusionCase run_case (int axis, bool low_active, bool high_active,
                        bool primitive_bc, bool quadratic)
{
    const Box domain(IntVect(2, 3, 4), IntVect(9, 14, 20));
    const BoxArray ba(domain);
    const DistributionMapping dm(ba);
    const Box domain_2d(IntVect(domain.smallEnd(0), domain.smallEnd(1), 0),
                        IntVect(domain.bigEnd(0), domain.bigEnd(1), 0));
    const BoxArray ba_2d(domain_2d);
    const std::array<Real,3> dx = {Real(0.3), Real(0.4), Real(0.7)};
    const Real alpha = Real(0.37);
    const Real a = Real(1.2);
    const Real b = Real(0.41);
    const Real c = quadratic ? Real(-0.083) : Real(0.0);
    const int ncell = domain.length(axis);

    MultiFab cell_data(ba, dm, NVAR_max, 2);
    MultiFab cell_prim(ba, dm, NPRIMVAR_max, 2);
    MultiFab cell_rhs(ba, dm, NVAR_max, 0);
    MultiFab xvel(BoxArray(surroundingNodes(domain, 0)), dm, 1, 1);
    MultiFab yvel(BoxArray(surroundingNodes(domain, 1)), dm, 1, 1);
    MultiFab xflux(BoxArray(surroundingNodes(domain, 0)), dm, 1, 0);
    MultiFab yflux(BoxArray(surroundingNodes(domain, 1)), dm, 1, 0);
    MultiFab zflux(BoxArray(surroundingNodes(domain, 2)), dm, 1, 0);
    MultiFab hfx_z(BoxArray(surroundingNodes(domain, 2)), dm, 1, 0);
    MultiFab qfx1_z(BoxArray(surroundingNodes(domain, 2)), dm, 1, 0);
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

    cell_data.setVal(Real(1.0));
    cell_prim.setVal(Real(0.0));
    cell_rhs.setVal(Real(0.0));
    xvel.setVal(Real(0.0));
    yvel.setVal(Real(0.0));
    mu_turb.setVal(Real(0.0));
    smn.setVal(Real(0.0));
    mf_mx.setVal(Real(1.0));
    mf_ux.setVal(Real(1.0));
    mf_vx.setVal(Real(1.0));
    mf_my.setVal(Real(1.0));
    mf_uy.setVal(Real(1.0));
    mf_vy.setVal(Real(1.0));

    for (MFIter mfi(cell_prim); mfi.isValid(); ++mfi) {
        auto prim = cell_prim.array(mfi);
        for (BoxIterator bit(cell_prim[mfi].box()); bit.ok(); ++bit) {
            const IntVect iv = bit();
            const int relative = iv[axis] - domain.smallEnd(axis);
            const Real xi = (Real(relative) + Real(0.5)) * dx[axis];
            prim(iv, PrimTheta_comp) = a + b*xi + c*xi*xi;
        }
        for (int side : {0, 1}) {
            const bool active = side == 0 ? low_active : high_active;
            const int face_index = side == 0 ? domain.smallEnd(axis)-1 : domain.bigEnd(axis)+1;
            if (!active) continue;
            for (BoxIterator bit(cell_prim[mfi].box()); bit.ok(); ++bit) {
                const IntVect iv = bit();
                if (iv[axis] == face_index) {
                    const Real xi = (side == 0 ? Real(0.0) : Real(ncell)*dx[axis]);
                    prim(iv, PrimTheta_comp) = a + b*xi + c*xi*xi;
                }
            }
        }
    }

    Vector<BCRec> bcs(NBCVAR_max);
    for (auto& bc : bcs) {
        bc.setLo(0, ERFBCType::foextrap); bc.setHi(0, ERFBCType::foextrap);
        bc.setLo(1, ERFBCType::foextrap); bc.setHi(1, ERFBCType::foextrap);
        bc.setLo(2, ERFBCType::foextrap); bc.setHi(2, ERFBCType::foextrap);
    }
    const auto active_type = primitive_bc ? ERFBCType::ext_dir_prim : ERFBCType::ext_dir;
    if (axis == 0) {
        if (low_active) bcs[BCVars::RhoTheta_bc_comp].setLo(0, active_type);
        if (high_active) bcs[BCVars::RhoTheta_bc_comp].setHi(0, active_type);
    } else if (axis == 1) {
        if (low_active) bcs[BCVars::RhoTheta_bc_comp].setLo(1, active_type);
        if (high_active) bcs[BCVars::RhoTheta_bc_comp].setHi(1, active_type);
    } else {
        if (low_active) bcs[BCVars::RhoTheta_bc_comp].setLo(2, active_type);
        if (high_active) bcs[BCVars::RhoTheta_bc_comp].setHi(2, active_type);
    }

    SolverChoice solver_choice;
    solver_choice.diffChoice.molec_diff_type = MolecDiffType::ConstantAlpha;
    solver_choice.diffChoice.alpha_T = alpha;
    solver_choice.turbChoice.resize(1);
    solver_choice.turbChoice[0].use_kturb = true;
    solver_choice.turbChoice[0].les_type = LESType::Smagorinsky;

    const GpuArray<Real, AMREX_SPACEDIM> dx_inv = {
        Real(1.0)/dx[0], Real(1.0)/dx[1], Real(1.0)/dx[2]};
    const GpuArray<Real, AMREX_SPACEDIM> grav = {Real(0.0), Real(0.0), Real(0.0)};
    Array4<const Real> tm_arr{};
    auto rhs_arr = cell_rhs[0].array();
    auto xflux_arr = xflux[0].array();
    auto yflux_arr = yflux[0].array();
    auto zflux_arr = zflux[0].array();
    auto hfx_arr = hfx_z[0].array();
    auto qfx1_arr = qfx1_z[0].array();
    auto qfx2_arr = qfx2_z[0].array();
    auto diss_arr = diss[0].array();
    DiffusionSrcForState_N(
        domain, domain, RhoTheta_comp, 1,
        xvel[0].const_array(), yvel[0].const_array(), cell_data[0].const_array(), cell_prim[0].const_array(),
        rhs_arr, xflux_arr, yflux_arr, zflux_arr, dx_inv,
        smn[0].const_array(), mf_mx[0].const_array(), mf_ux[0].const_array(), mf_vx[0].const_array(),
        mf_my[0].const_array(), mf_uy[0].const_array(), mf_vy[0].const_array(), hfx_arr,
        qfx1_arr, qfx2_arr, diss_arr, mu_turb[0].const_array(), solver_choice, 0,
        tm_arr, grav, bcs.data(), false, Real(0.0));
    Gpu::streamSynchronize();

    DiffusionCase result;
    for (MFIter mfi(cell_rhs); mfi.isValid(); ++mfi) {
        const auto rhs = cell_rhs[0].const_array();
        for (BoxIterator bit(domain); bit.ok(); ++bit) {
            const IntVect iv = bit();
            result.max_rhs = std::max(result.max_rhs, std::abs(rhs(iv, RhoTheta_comp)));
        }
        const auto xf = xflux[0].const_array();
        const auto yf = yflux[0].const_array();
        const auto zf = zflux[0].const_array();
        const int i = domain.smallEnd(0) + 1;
        const int j = domain.smallEnd(1) + 1;
        const int k = domain.smallEnd(2) + 1;
        if (axis == 0) {
            result.low_flux = xf(domain.smallEnd(0), j, k);
            result.high_flux = xf(domain.bigEnd(0)+1, j, k);
        } else if (axis == 1) {
            result.low_flux = yf(i, domain.smallEnd(1), k);
            result.high_flux = yf(i, domain.bigEnd(1)+1, k);
        } else {
            result.low_flux = zf(i, j, domain.smallEnd(2));
            result.high_flux = zf(i, j, domain.bigEnd(2)+1);
        }
    }
    return result;
}

// Regression motivation:
// Before Stage 0, the ConstantAlpha + turbulence x-high predicate gated the
// low-wall predicate a second time. These cases call DiffusionSrcForState_N,
// the production operator, on unequal dimensions/spacings and use the exact
// polynomial diffusion identity as the independent oracle. A failure means a
// fixed-value face stencil or the integrated flux divergence is wrong.
TEST(StateDiffusionDirichlet, ConstantAlphaTurbulenceLinearAndQuadraticRotated)
{
    for (int axis = 0; axis < 3; ++axis) {
        for (bool primitive : {false, true}) {
            for (bool quadratic : {false, true}) {
                for (bool low_active : {false, true}) {
                    for (bool high_active : {false, true}) {
                        const auto result = run_case(axis, low_active, high_active,
                                                     primitive, quadratic);
                        const Real scale = Real(1.0) + (quadratic ? Real(0.17) : Real(0.0));
                        const Real tol = Real(2.0e-11) * scale;
                        const Real expected_flux_low = -Real(0.37) * Real(0.41);
                        const Real coord_length = axis == 0 ? Real(8)*Real(0.3) :
                            (axis == 1 ? Real(12)*Real(0.4) : Real(17)*Real(0.7));
                        const Real expected_flux_high = quadratic ?
                            -Real(0.37) * (Real(0.41) + Real(2.0)*Real(-0.083)*coord_length) :
                            expected_flux_low;
                        if (!quadratic) {
                            EXPECT_NEAR(result.max_rhs, Real(0.0), tol);
                        } else {
                            EXPECT_NEAR(result.max_rhs, std::abs(Real(2.0)*Real(0.37)*Real(-0.083)), tol);
                        }
                        if (low_active) EXPECT_NEAR(result.low_flux, expected_flux_low, tol);
                        if (high_active) EXPECT_NEAR(result.high_flux, expected_flux_high, tol);
                    }
                }
            }
        }
    }
}

} // namespace
