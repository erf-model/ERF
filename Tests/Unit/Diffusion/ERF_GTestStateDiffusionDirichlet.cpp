#include <AMReX_MultiFab.H>
#include <AMReX_GpuContainers.H>
#include <AMReX_Reduce.H>

#include <ERF_Diffusion.H>

#include <gtest/gtest.h>

#include <algorithm>
#include <array>
#include <cmath>
#include <limits>

using namespace amrex;

namespace {

struct DiffusionCase {
    Real max_state_error = 0.0;
    Real max_rhs_error = 0.0;
    Real low_flux_error = 0.0;
    Real high_flux_error = 0.0;
    Real integrated_rhs = 0.0;
    Real boundary_flux_rhs = 0.0;
};

Real scaled_tolerance (Real scale, Real factor)
{
    return factor * std::numeric_limits<Real>::epsilon() *
           std::max(Real(1.0), std::abs(scale));
}

Real face_error (const Box& face_box,
                 const Array4<const Real>& flux,
                 Real expected)
{
    ReduceOps<ReduceOpMax> reduce_op;
    ReduceData<Real> reduce_data(reduce_op);
    reduce_op.eval(face_box, reduce_data,
        [=] AMREX_GPU_DEVICE (int i, int j, int k) -> GpuTuple<Real> {
            return std::abs(flux(i,j,k) - expected);
        });
    return get<0>(reduce_data.value());
}

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
    const int active_small = domain.smallEnd(axis);
    const Real active_dx = dx[axis];
    const Real active_length = Real(ncell) * active_dx;
    const Real expected_rhs = Real(2.0) * alpha * c;
    const Real expected_low_flux = -alpha * b;
    const Real expected_high_flux = -alpha * (b + Real(2.0)*c*active_length);

    MultiFab cell_data(ba, dm, NVAR_max, 2);
    MultiFab cell_prim(ba, dm, NPRIMVAR_max, 2);
    MultiFab cell_rhs(ba, dm, NVAR_max, 0);
    MultiFab state_error(ba, dm, 1, 0);
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

    for (MFIter mfi(cell_prim, TilingIfNotGPU()); mfi.isValid(); ++mfi) {
        const Box gbx = mfi.growntilebox(cell_prim.nGrowVect());
        auto prim = cell_prim.array(mfi);
        ParallelFor(gbx, [=] AMREX_GPU_DEVICE (int i, int j, int k) noexcept {
            const int index = (axis == 0) ? i : ((axis == 1) ? j : k);
            Real xi = (Real(index - active_small) + Real(0.5)) * active_dx;
            if (axis == 0) {
                if (low_active && i == domain.smallEnd(0)-1) xi = Real(0.0);
                if (high_active && i == domain.bigEnd(0)+1) xi = active_length;
            } else if (axis == 1) {
                if (low_active && j == domain.smallEnd(1)-1) xi = Real(0.0);
                if (high_active && j == domain.bigEnd(1)+1) xi = active_length;
            } else {
                if (low_active && k == domain.smallEnd(2)-1) xi = Real(0.0);
                if (high_active && k == domain.bigEnd(2)+1) xi = active_length;
            }
            prim(i,j,k,PrimTheta_comp) = a + b*xi + c*xi*xi;
        });
    }
    for (MFIter mfi(state_error, TilingIfNotGPU()); mfi.isValid(); ++mfi) {
        const Box bx = mfi.tilebox();
        const auto prim = cell_prim.const_array(mfi);
        const auto error = state_error.array(mfi);
        ParallelFor(bx, [=] AMREX_GPU_DEVICE (int i, int j, int k) noexcept {
            const int index = (axis == 0) ? i : ((axis == 1) ? j : k);
            const Real xi = (Real(index - active_small) + Real(0.5)) * active_dx;
            error(i,j,k) = std::abs(prim(i,j,k,PrimTheta_comp) -
                                    (a + b*xi + c*xi*xi));
        });
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

    Gpu::DeviceVector<BCRec> bcs_d(bcs.size());
    Gpu::copy(Gpu::hostToDevice, bcs.begin(), bcs.end(), bcs_d.begin());

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
    auto hfx_arr = hfx_z[0].array();
    auto qfx1_arr = qfx1_z[0].array();
    auto qfx2_arr = qfx2_z[0].array();
    auto diss_arr = diss[0].array();
    DiffusionSrcForState_N(
        domain, domain, RhoTheta_comp, 1,
        xvel[0].const_array(), yvel[0].const_array(), cell_data[0].const_array(), cell_prim[0].const_array(),
        cell_rhs[0].array(), xflux[0].array(), yflux[0].array(), zflux[0].array(), dx_inv,
        smn[0].const_array(), mf_mx[0].const_array(), mf_ux[0].const_array(), mf_vx[0].const_array(),
        mf_my[0].const_array(), mf_uy[0].const_array(), mf_vy[0].const_array(), hfx_arr,
        qfx1_arr, qfx2_arr, diss_arr, mu_turb[0].const_array(), solver_choice, 0,
        tm_arr, grav, bcs_d.data(), false, Real(0.0));
    Gpu::streamSynchronize();

    MultiFab rhs_error(ba, dm, 1, 0);
    MultiFab weighted_rhs(ba, dm, 1, 0);
    MultiFab boundary_terms(ba, dm, 1, 0);
    const Real cell_volume = dx[0] * dx[1] * dx[2];
    const Real face_area = (axis == 0) ? dx[1]*dx[2] :
                           ((axis == 1) ? dx[0]*dx[2] : dx[0]*dx[1]);
    for (MFIter mfi(rhs_error, TilingIfNotGPU()); mfi.isValid(); ++mfi) {
        const Box bx = mfi.tilebox();
        const auto rhs = cell_rhs.const_array(mfi);
        const auto rhs_err = rhs_error.array(mfi);
        const auto weighted = weighted_rhs.array(mfi);
        const auto boundary = boundary_terms.array(mfi);
        const auto xf = xflux.const_array(mfi);
        const auto yf = yflux.const_array(mfi);
        const auto zf = zflux.const_array(mfi);
        ParallelFor(bx, [=] AMREX_GPU_DEVICE (int i, int j, int k) noexcept {
            rhs_err(i,j,k) = std::abs(rhs(i,j,k,RhoTheta_comp) - expected_rhs);
            weighted(i,j,k) = cell_volume * rhs(i,j,k,RhoTheta_comp);
            Real contribution = Real(0.0);
            if (axis == 0) {
                if (i == domain.smallEnd(0)) {
                    contribution += face_area * xf(i,j,k);
                }
                if (i == domain.bigEnd(0)) {
                    contribution -= face_area * xf(i+1,j,k);
                }
            } else if (axis == 1) {
                if (j == domain.smallEnd(1)) {
                    contribution += face_area * yf(i,j,k);
                }
                if (j == domain.bigEnd(1)) {
                    contribution -= face_area * yf(i,j+1,k);
                }
            } else {
                if (k == domain.smallEnd(2)) {
                    contribution += face_area * zf(i,j,k);
                }
                if (k == domain.bigEnd(2)) {
                    contribution -= face_area * zf(i,j,k+1);
                }
            }
            boundary(i,j,k) = contribution;
        });
    }
    Gpu::streamSynchronize();

    DiffusionCase result;
    result.max_rhs_error = rhs_error.norm0(0, 0, false);
    result.integrated_rhs = weighted_rhs.sum(0);
    result.boundary_flux_rhs = boundary_terms.sum(0);

    const Box lo_face(IntVect(domain.smallEnd(0), domain.smallEnd(1), domain.smallEnd(2)),
                      IntVect(domain.smallEnd(0), domain.bigEnd(1), domain.bigEnd(2)));
    const Box hi_face(IntVect(domain.bigEnd(0)+1, domain.smallEnd(1), domain.smallEnd(2)),
                      IntVect(domain.bigEnd(0)+1, domain.bigEnd(1), domain.bigEnd(2)));
    const Box lo_yface(IntVect(domain.smallEnd(0), domain.smallEnd(1), domain.smallEnd(2)),
                       IntVect(domain.bigEnd(0), domain.smallEnd(1), domain.bigEnd(2)));
    const Box hi_yface(IntVect(domain.smallEnd(0), domain.bigEnd(1)+1, domain.smallEnd(2)),
                       IntVect(domain.bigEnd(0), domain.bigEnd(1)+1, domain.bigEnd(2)));
    const Box lo_zface(IntVect(domain.smallEnd(0), domain.smallEnd(1), domain.smallEnd(2)),
                       IntVect(domain.bigEnd(0), domain.bigEnd(1), domain.smallEnd(2)));
    const Box hi_zface(IntVect(domain.smallEnd(0), domain.smallEnd(1), domain.bigEnd(2)+1),
                       IntVect(domain.bigEnd(0), domain.bigEnd(1), domain.bigEnd(2)+1));
    if (axis == 0) {
        if (low_active) result.low_flux_error = face_error(lo_face, xflux[0].const_array(), expected_low_flux);
        if (high_active) result.high_flux_error = face_error(hi_face, xflux[0].const_array(), expected_high_flux);
    } else if (axis == 1) {
        if (low_active) result.low_flux_error = face_error(lo_yface, yflux[0].const_array(), expected_low_flux);
        if (high_active) result.high_flux_error = face_error(hi_yface, yflux[0].const_array(), expected_high_flux);
    } else {
        if (low_active) result.low_flux_error = face_error(lo_zface, zflux[0].const_array(), expected_low_flux);
        if (high_active) result.high_flux_error = face_error(hi_zface, zflux[0].const_array(), expected_high_flux);
    }
    result.max_state_error = state_error.norm0(0, 0, false);
    return result;
}

// Motivation: The production scalar-diffusion operator must remain device-safe
// and preserve the signed face-flux/RHS balance for rotated linear and
// quadratic profiles.
TEST(StateDiffusionDirichlet, ConstantAlphaTurbulenceLinearAndQuadraticRotated)
{
    for (int axis = 0; axis < 3; ++axis) {
        for (bool primitive : {false, true}) {
            for (bool quadratic : {false, true}) {
                for (bool low_active : {false, true}) {
                    for (bool high_active : {false, true}) {
                        const auto result = run_case(axis, low_active, high_active,
                                                     primitive, quadratic);
                        const Real alpha = Real(0.37);
                        const Real c = quadratic ? Real(-0.083) : Real(0.0);
                        const Real b = Real(0.41);
                        const Real length = axis == 0 ? Real(8)*Real(0.3) :
                            (axis == 1 ? Real(12)*Real(0.4) : Real(17)*Real(0.7));
                        const Real expected_rhs = Real(2.0)*alpha*c;
                        const Real expected_low_flux = -alpha*b;
                        const Real expected_high_flux = -alpha*(b + Real(2.0)*c*length);
                        const Real cell_volume = Real(0.3)*Real(0.4)*Real(0.7);
                        const Real ntransverse = axis == 0 ? Real(12*17) :
                            (axis == 1 ? Real(8*17) : Real(8*12));
                        const Real face_area = axis == 0 ? Real(0.4)*Real(0.7) :
                            (axis == 1 ? Real(0.3)*Real(0.7) : Real(0.3)*Real(0.4));
                        const Real expected_integrated = expected_rhs * cell_volume *
                            Real(8*12*17);
                        const Real expected_boundary = face_area*ntransverse *
                            (expected_low_flux - expected_high_flux);
                        const Real state_scale = std::max(Real(1.0),
                            std::max(std::abs(Real(1.2)), std::abs(Real(1.2) + b*length + c*length*length)));
                        const Real flux_scale = std::max(Real(1.0),
                            alpha * (std::abs(b) + Real(2.0)*std::abs(c)*length));
                        const Real rhs_scale = std::max(Real(1.0),
                            std::max(std::abs(expected_rhs), alpha*state_scale/(Real(0.3)*Real(0.3))));
                        const Real budget_scale = std::max(Real(1.0),
                            std::max(std::abs(expected_integrated),
                                     std::abs(face_area*ntransverse*expected_low_flux) +
                                     std::abs(face_area*ntransverse*expected_high_flux)));
                        const Real state_tol = scaled_tolerance(state_scale, Real(64.0));
                        const Real flux_tol = scaled_tolerance(flux_scale, Real(128.0));
                        const Real rhs_tol = scaled_tolerance(rhs_scale, Real(256.0));
                        const Real budget_tol = scaled_tolerance(budget_scale, Real(512.0));
                        EXPECT_LE(result.max_state_error, state_tol);
                        EXPECT_LE(result.max_rhs_error, rhs_tol);
                        if (low_active) {
                            EXPECT_LE(result.low_flux_error, flux_tol);
                        }
                        if (high_active) {
                            EXPECT_LE(result.high_flux_error, flux_tol);
                        }
                        EXPECT_NEAR(result.integrated_rhs, expected_integrated, budget_tol);
                        EXPECT_NEAR(result.boundary_flux_rhs, expected_boundary, budget_tol);
                        EXPECT_NEAR(result.integrated_rhs, result.boundary_flux_rhs, budget_tol);
                    }
                }
            }
        }
    }
}

} // namespace
