#include "ERF_IndexDefines.H"
#include "ERF_Constants.H"
#include "ERF_EOS.H"
#include "ERF_ShocColumnData.H"
#include "ERF_ShocDriver.H"
#include "ERF_ShocPreprocess.H"
#include "ERF_ShocStructure.H"
#include "ERF_ShocTestUtils.H"
#include "ERF_ShocTKE.H"

#include <gtest/gtest.h>

#include <AMReX.H>
#include <AMReX_BoxArray.H>
#include <AMReX_DistributionMapping.H>
#include <AMReX_Geometry.H>
#include <AMReX_GpuContainers.H>
#include <AMReX_MultiFab.H>
#include <AMReX_ParmParse.H>
#include <AMReX_RealBox.H>
#include <AMReX_Vector.H>

#include <cmath>
#include <string>

namespace
{
using namespace amrex::literals;

using amrex::Box;
using amrex::DistributionMapping;
using amrex::Geometry;
using amrex::GpuArray;
using amrex::IntVect;
using amrex::MFIter;
using amrex::MultiFab;
using amrex::ParmParse;
using amrex::Real;

constexpr int nx = 5;
constexpr int ny = 1;
constexpr int nz = 5;
constexpr Real tol = 1.0e-12_rt;

using FixtureMap = std::map<std::string, amrex::Vector<Real>>;

amrex::BoxArray
make_nodal_zphys_boxarray (const amrex::BoxArray& ba)
{
    amrex::BoxArray zphys_ba(ba);
    zphys_ba.surroundingNodes();
    return zphys_ba;
}

template <int N>
GpuArray<Real, N>
fixture_gpu_array (const FixtureMap& fixture,
                   const std::string& key)
{
    const auto& values = shoc_test::fixture_values(fixture, key);
    AMREX_ALWAYS_ASSERT_WITH_MESSAGE(static_cast<int>(values.size()) == N,
                                     "SHOC driver fixture size mismatch");

    GpuArray<Real, N> result{};
    for (int i = 0; i < N; ++i) {
        result[i] = values[i];
    }
    return result;
}

class ScopedParmParseString
{
public:
    ScopedParmParseString (const char* prefix,
                           const char* name,
                           const std::string& value)
        : m_pp(prefix),
          m_name(name)
    {
        m_had_previous = m_pp.query(m_name, m_previous);
        m_pp.remove(m_name);
        m_pp.add(m_name, value);
    }

    ~ScopedParmParseString ()
    {
        m_pp.remove(m_name);
        if (m_had_previous) {
            m_pp.add(m_name, m_previous);
        }
    }

private:
    ParmParse m_pp;
    std::string m_name;
    std::string m_previous;
    bool m_had_previous = false;
};

void
expect_multifab_finite (const MultiFab& mf,
                        int start_comp,
                        int ncomp,
                        const char* name)
{
    EXPECT_FALSE(mf.contains_nan(start_comp, ncomp, 0)) << name << " contains NaN";
    for (int comp = start_comp; comp < start_comp + ncomp; ++comp) {
        const Real comp_min = mf.min(comp);
        const Real comp_max = mf.max(comp);
        EXPECT_TRUE(std::isfinite(comp_min)) << name << " min is not finite for comp " << comp;
        EXPECT_TRUE(std::isfinite(comp_max)) << name << " max is not finite for comp " << comp;
    }
}

void
expect_component_nonnegative (const MultiFab& mf,
                              int comp,
                              const char* name)
{
    expect_multifab_finite(mf, comp, 1, name);
    EXPECT_GE(mf.min(comp), -tol) << name << " has negative values";
}

void
expect_component_between (const MultiFab& mf,
                          int comp,
                          Real lo,
                          Real hi,
                          const char* name)
{
    expect_multifab_finite(mf, comp, 1, name);
    EXPECT_GE(mf.min(comp), lo - tol) << name << " falls below lower bound";
    EXPECT_LE(mf.max(comp), hi + tol) << name << " exceeds upper bound";
}

Real
component_max_abs (const MultiFab& mf, int comp)
{
    const Real comp_min = mf.min(comp);
    const Real comp_max = mf.max(comp);
    return amrex::max(amrex::Math::abs(comp_min), amrex::Math::abs(comp_max));
}

Real
component_max_abs_difference (const MultiFab& after,
                              int after_comp,
                              const MultiFab& before,
                              int before_comp)
{
    MultiFab delta(after.boxArray(), after.DistributionMap(), 1, 0);
    MultiFab::Copy(delta, after, after_comp, 0, 1, 0);
    MultiFab::Subtract(delta, before, before_comp, 0, 1, 0);
    return delta.norm0(0, 0, true);
}

void
expect_component_minmax_match (const MultiFab& lhs,
                               const MultiFab& rhs,
                               int comp,
                               const char* name)
{
    EXPECT_NEAR(lhs.min(comp), rhs.min(comp), tol) << name << " min mismatch";
    EXPECT_NEAR(lhs.max(comp), rhs.max(comp), tol) << name << " max mismatch";
}

void
expect_shoc_column_finite (const amrex::FArrayBox& fab,
                           int nlev,
                           const char* name)
{
    const auto arr = fab.const_array();
    for (int ic = 0; ic < fab.box().length(0); ++ic) {
        for (int k = 0; k < nlev; ++k) {
            EXPECT_TRUE(std::isfinite(arr(ic,k,0))) << name << " is not finite at (" << ic << "," << k << ")";
        }
    }
}

void
seed_first_step_turbulence (ShocColumnData& col,
                            const MultiFab& cons,
                            const MultiFab& eddy_diffs,
                            const MFIter& mfi)
{
    const auto rho = cons.const_array(mfi);
    const auto host_diff = eddy_diffs.const_array(mfi);
    auto tk = col.tk.array();
    auto tkh = col.tkh.array();
    const auto layout = col.layout;
    const Box xy_box = amrex::makeSlab(mfi.validbox(), 2, layout.kmin);

    amrex::ParallelFor(xy_box, [=] AMREX_GPU_DEVICE (int i, int j, int) noexcept
    {
        const int ic = shoc_column_index(layout, i, j);
        for (int kk = 0; kk < layout.nlev; ++kk) {
            const int k = layout.kmin + kk;
            const Real rho_safe = amrex::max(rho(i,j,k,Rho_comp), 1.0e-12_rt);
            tk(ic,kk,0) = amrex::max(0.0_rt, host_diff(i,j,k,EddyDiff::Mom_v) / rho_safe);
            tkh(ic,kk,0) = amrex::max(0.0_rt, host_diff(i,j,k,EddyDiff::Theta_v) / rho_safe);
        }
    });
}

void
fill_multifab (MultiFab& mf, Real value)
{
    for (MFIter mfi(mf, false); mfi.isValid(); ++mfi) {
        const Box& bx = mfi.validbox();
        auto arr = mf.array(mfi);
        const int ncomp = mf.nComp();
        amrex::ParallelFor(bx, ncomp, [=] AMREX_GPU_DEVICE (int i, int j, int k, int n) noexcept
        {
            arr(i,j,k,n) = value;
        });
    }
}

void
set_eddy_diff_sentinel (MultiFab& eddy_diffs,
                        Real km,
                        Real kh,
                        Real tke,
                        Real scalar,
                        Real q,
                        Real lengthscale)
{
    fill_multifab(eddy_diffs, 0.0_rt);

    for (MFIter mfi(eddy_diffs, false); mfi.isValid(); ++mfi) {
        const Box& bx = mfi.validbox();
        auto arr = eddy_diffs.array(mfi);
        amrex::ParallelFor(bx, [=] AMREX_GPU_DEVICE (int i, int j, int k) noexcept
        {
            arr(i,j,k,EddyDiff::Mom_v) = km;
            arr(i,j,k,EddyDiff::Theta_v) = kh;
            arr(i,j,k,EddyDiff::KE_v) = tke;
            arr(i,j,k,EddyDiff::Scalar_v) = scalar;
            arr(i,j,k,EddyDiff::Q_v) = q;
            arr(i,j,k,EddyDiff::Turb_lengthscale) = lengthscale;
        });
    }
}

void
initialize_basic_column_state (MultiFab& cons,
                               MultiFab& zvel,
                               MultiFab& z_phys_nd,
                               MultiFab& hfx3,
                               MultiFab& qfx3,
                               MultiFab& eddy_diffs)
{
    fill_multifab(cons, 0.0_rt);
    fill_multifab(zvel, 0.0_rt);
    fill_multifab(z_phys_nd, 0.0_rt);
    fill_multifab(hfx3, 0.0_rt);
    fill_multifab(qfx3, 0.0_rt);
    fill_multifab(eddy_diffs, 0.0_rt);

    for (MFIter mfi(cons, false); mfi.isValid(); ++mfi) {
        const Box& bx = mfi.validbox();
        auto arr = cons.array(mfi);
        amrex::ParallelFor(bx, [=] AMREX_GPU_DEVICE (int i, int j, int k) noexcept
        {
            const Real rho = 1.0_rt + 0.05_rt * k;
            const Real theta = 300.0_rt + 0.5_rt * k;
            const Real qv = 0.010_rt - 0.001_rt * k;
            const Real qc = 0.0_rt;

            arr(i,j,k,Rho_comp) = rho;
            arr(i,j,k,RhoTheta_comp) = rho * theta;
            arr(i,j,k,RhoKE_comp) = rho * (0.2_rt + 0.01_rt * k);
            arr(i,j,k,RhoQ1_comp) = rho * qv;
            arr(i,j,k,RhoQ2_comp) = rho * qc;
        });
    }

    for (MFIter mfi(zvel, false); mfi.isValid(); ++mfi) {
        const Box& bx = mfi.validbox();
        auto arr = zvel.array(mfi);
        amrex::ParallelFor(bx, [=] AMREX_GPU_DEVICE (int i, int j, int k) noexcept
        {
            arr(i,j,k) = 0.0_rt;
        });
    }

    for (MFIter mfi(z_phys_nd, false); mfi.isValid(); ++mfi) {
        const Box& bx = mfi.validbox();
        auto arr = z_phys_nd.array(mfi);
        amrex::ParallelFor(bx, [=] AMREX_GPU_DEVICE (int i, int j, int k) noexcept
        {
            arr(i,j,k) = 100.0_rt * k;
        });
    }

    for (MFIter mfi(hfx3, false); mfi.isValid(); ++mfi) {
        const Box& bx = mfi.validbox();
        auto arr = hfx3.array(mfi);
        amrex::ParallelFor(bx, [=] AMREX_GPU_DEVICE (int i, int j, int k) noexcept
        {
            arr(i,j,k) = (k == 0) ? 0.02_rt : 0.0_rt;
        });
    }

    for (MFIter mfi(qfx3, false); mfi.isValid(); ++mfi) {
        const Box& bx = mfi.validbox();
        auto arr = qfx3.array(mfi);
        amrex::ParallelFor(bx, [=] AMREX_GPU_DEVICE (int i, int j, int k) noexcept
        {
            arr(i,j,k) = (k == 0) ? 1.0e-4_rt : 0.0_rt;
        });
    }

    for (MFIter mfi(eddy_diffs, false); mfi.isValid(); ++mfi) {
        const Box& bx = mfi.validbox();
        auto arr = eddy_diffs.array(mfi);
        amrex::ParallelFor(bx, [=] AMREX_GPU_DEVICE (int i, int j, int k) noexcept
        {
            const Real coeff = 0.5_rt + 0.05_rt * k;
            arr(i,j,k,EddyDiff::Mom_v) = coeff;
            arr(i,j,k,EddyDiff::Theta_v) = coeff;
            arr(i,j,k,EddyDiff::KE_v) = coeff;
            arr(i,j,k,EddyDiff::Scalar_v) = coeff;
            arr(i,j,k,EddyDiff::Q_v) = coeff;
            arr(i,j,k,EddyDiff::Turb_lengthscale) = 10.0_rt;
        });
    }
}

void
initialize_terrain_geometry (MultiFab& z_phys_nd)
{
    for (MFIter mfi(z_phys_nd, false); mfi.isValid(); ++mfi) {
        const Box& bx = mfi.validbox();
        auto arr = z_phys_nd.array(mfi);
        amrex::ParallelFor(bx, [=] AMREX_GPU_DEVICE (int i, int j, int k) noexcept
        {
            arr(i,j,k) = 100.0_rt * k + 10.0_rt * static_cast<Real>(i) + 2.0_rt * static_cast<Real>(j);
        });
    }
}

void
initialize_face_velocities_with_seam_offsets (MultiFab& xvel,
                                              MultiFab& yvel,
                                              const Box& domain,
                                              Real seam_delta_x,
                                              Real seam_delta_y)
{
    for (MFIter mfi(xvel, false); mfi.isValid(); ++mfi) {
        const Box& bx = mfi.validbox();
        const bool is_left = bx.smallEnd(0) == domain.smallEnd(0);
        const int seam_x = is_left ? bx.bigEnd(0) : bx.smallEnd(0);
        const Real seam_sign = is_left ? 1.0_rt : -1.0_rt;
        auto arr = xvel.array(mfi);
        amrex::ParallelFor(bx, [=] AMREX_GPU_DEVICE (int i, int j, int k) noexcept
        {
            Real value = 2.0_rt + 0.1_rt * k;
            if (i == seam_x) {
                value += seam_sign * seam_delta_x;
            }
            arr(i,j,k) = value;
        });
    }

    for (MFIter mfi(yvel, false); mfi.isValid(); ++mfi) {
        const Box& bx = mfi.validbox();
        const bool is_bottom = bx.smallEnd(1) == domain.smallEnd(1);
        const int seam_y = is_bottom ? bx.bigEnd(1) : bx.smallEnd(1);
        const Real seam_sign = is_bottom ? 1.0_rt : -1.0_rt;
        auto arr = yvel.array(mfi);
        amrex::ParallelFor(bx, [=] AMREX_GPU_DEVICE (int i, int j, int k) noexcept
        {
            Real value = 1.0_rt + 0.05_rt * k;
            if (j == seam_y) {
                value += seam_sign * seam_delta_y;
            }
            arr(i,j,k) = value;
        });
    }
}

void
initialize_face_stresses_with_seam_offsets (MultiFab& tau13,
                                            MultiFab& tau23,
                                            const Box& domain,
                                            Real seam_delta_x,
                                            Real seam_delta_y)
{
    for (MFIter mfi(tau13, false); mfi.isValid(); ++mfi) {
        const Box& bx = mfi.validbox();
        const bool is_left = bx.smallEnd(0) == domain.smallEnd(0);
        const int seam_x = is_left ? bx.bigEnd(0) : bx.smallEnd(0);
        const Real seam_sign = is_left ? 1.0_rt : -1.0_rt;
        auto arr = tau13.array(mfi);
        amrex::ParallelFor(bx, [=] AMREX_GPU_DEVICE (int i, int j, int k) noexcept
        {
            Real value = 0.04_rt;
            if (i == seam_x) {
                value += seam_sign * seam_delta_x;
            }
            arr(i,j,k) = (k == 0) ? value : 0.0_rt;
        });
    }

    for (MFIter mfi(tau23, false); mfi.isValid(); ++mfi) {
        const Box& bx = mfi.validbox();
        const bool is_bottom = bx.smallEnd(1) == domain.smallEnd(1);
        const int seam_y = is_bottom ? bx.bigEnd(1) : bx.smallEnd(1);
        const Real seam_sign = is_bottom ? 1.0_rt : -1.0_rt;
        auto arr = tau23.array(mfi);
        amrex::ParallelFor(bx, [=] AMREX_GPU_DEVICE (int i, int j, int k) noexcept
        {
            Real value = 0.03_rt;
            if (j == seam_y) {
                value += seam_sign * seam_delta_y;
            }
            arr(i,j,k) = (k == 0) ? value : 0.0_rt;
        });
    }
}

void
initialize_face_stress_gradients (MultiFab& tau13,
                                  MultiFab& tau23)
{
    for (MFIter mfi(tau13, false); mfi.isValid(); ++mfi) {
        const Box& bx = mfi.validbox();
        auto arr = tau13.array(mfi);
        amrex::ParallelFor(bx, [=] AMREX_GPU_DEVICE (int i, int j, int k) noexcept
        {
            arr(i,j,k) = (k == 0) ? (10.0_rt + static_cast<Real>(i)) : 0.0_rt;
        });
    }

    for (MFIter mfi(tau23, false); mfi.isValid(); ++mfi) {
        const Box& bx = mfi.validbox();
        auto arr = tau23.array(mfi);
        amrex::ParallelFor(bx, [=] AMREX_GPU_DEVICE (int i, int j, int k) noexcept
        {
            arr(i,j,k) = (k == 0) ? (20.0_rt + static_cast<Real>(j)) : 0.0_rt;
        });
    }
}

void
expect_face_overlap_matches (const MultiFab& mf,
                             const Box& domain,
                             int seam_dir,
                             int seam_index,
                             int fixed_i,
                             int fixed_j,
                             int fixed_k,
                             const char* name)
{
    Real a = 0.0_rt;
    Real b = 0.0_rt;
    bool have_a = false;
    bool have_b = false;

    for (MFIter mfi(mf, false); mfi.isValid(); ++mfi) {
        const Box& bx = mfi.validbox();
        const auto arr = mf.const_array(mfi);

        if (seam_dir == 0) {
            if (bx.smallEnd(0) == domain.smallEnd(0) &&
                bx.smallEnd(1) == domain.smallEnd(1)) {
                a = arr(seam_index, fixed_j, fixed_k);
                have_a = true;
            } else if (bx.smallEnd(0) == seam_index &&
                       bx.smallEnd(1) == domain.smallEnd(1)) {
                b = arr(seam_index, fixed_j, fixed_k);
                have_b = true;
            }
        } else {
            if (bx.smallEnd(0) == domain.smallEnd(0) &&
                bx.smallEnd(1) == domain.smallEnd(1)) {
                a = arr(fixed_i, seam_index, fixed_k);
                have_a = true;
            } else if (bx.smallEnd(0) == domain.smallEnd(0) &&
                       bx.smallEnd(1) == seam_index) {
                b = arr(fixed_i, seam_index, fixed_k);
                have_b = true;
            }
        }
    }

    EXPECT_TRUE(have_a) << name << " missing first seam copy";
    EXPECT_TRUE(have_b) << name << " missing second seam copy";
    EXPECT_NEAR(a, b, tol) << name << " seam copies diverged";
}

FixtureMap
load_driver_fixture ()
{
    return shoc_test::read_named_fixture_vectors(
        "implicit_energy/e3sm_update_prognostics_implicit_multicolumn.txt");
}

void
initialize_state (MultiFab& cons,
                  const FixtureMap& fixture)
{
    const auto rho_zt = fixture_gpu_array<nz>(fixture, "rho_zt");
    const auto thetal = fixture_gpu_array<nz>(fixture, "thetal");
    const auto qw = fixture_gpu_array<nz>(fixture, "qw");
    const auto tke = fixture_gpu_array<nz>(fixture, "tke");

    fill_multifab(cons, 0.0_rt);

    for (MFIter mfi(cons, false); mfi.isValid(); ++mfi) {
        const Box& bx = mfi.validbox();
        auto arr = cons.array(mfi);

        amrex::ParallelFor(bx, [=] AMREX_GPU_DEVICE (int i, int j, int k) noexcept
        {
            const Real rho = rho_zt[k];
            const Real theta = thetal[k];
            const Real qv = qw[k];
            const Real qc = 0.0_rt;

            arr(i,j,k,Rho_comp) = rho;
            arr(i,j,k,RhoTheta_comp) = rho * theta;
            arr(i,j,k,RhoKE_comp) = rho * tke[k];
            arr(i,j,k,RhoQ1_comp) = rho * qv;
            arr(i,j,k,RhoQ2_comp) = rho * qc;
        });
    }
}

void
initialize_velocity (MultiFab& xvel,
                     MultiFab& yvel,
                     MultiFab& zvel,
                     const FixtureMap& fixture)
{
    const auto u = fixture_gpu_array<nz>(fixture, "u");
    const auto v = fixture_gpu_array<nz>(fixture, "v");

    fill_multifab(zvel, 0.0_rt);

    for (MFIter mfi(xvel, false); mfi.isValid(); ++mfi) {
        const Box& bx = mfi.validbox();
        auto arr = xvel.array(mfi);
        amrex::ParallelFor(bx, [=] AMREX_GPU_DEVICE (int i, int j, int k) noexcept
        {
            arr(i,j,k) = u[k];
        });
    }

    for (MFIter mfi(yvel, false); mfi.isValid(); ++mfi) {
        const Box& bx = mfi.validbox();
        auto arr = yvel.array(mfi);
        amrex::ParallelFor(bx, [=] AMREX_GPU_DEVICE (int i, int j, int k) noexcept
        {
            arr(i,j,k) = v[k];
        });
    }
}

void
initialize_geometry (MultiFab& z_phys_nd,
                     const FixtureMap& fixture)
{
    const auto zi_grid = fixture_gpu_array<nz + 1>(fixture, "zi_grid");

    for (MFIter mfi(z_phys_nd, false); mfi.isValid(); ++mfi) {
        const Box& bx = mfi.validbox();
        auto arr = z_phys_nd.array(mfi);
        amrex::ParallelFor(bx, [=] AMREX_GPU_DEVICE (int i, int j, int k) noexcept
        {
            arr(i,j,k) = zi_grid[k];
        });
    }
}

void
initialize_surface_fluxes (MultiFab& hfx3,
                           MultiFab& qfx3,
                           MultiFab& tau13,
                           MultiFab& tau23,
                           const FixtureMap& fixture)
{
    const auto rho_zt = fixture_gpu_array<nz>(fixture, "rho_zt");
    const auto wthl_sfc = fixture_gpu_array<nx * ny>(fixture, "wthl_sfc");
    const auto wqw_sfc = fixture_gpu_array<nx * ny>(fixture, "wqw_sfc");
    const auto uw_sfc = fixture_gpu_array<nx * ny>(fixture, "uw_sfc");
    const auto vw_sfc = fixture_gpu_array<nx * ny>(fixture, "vw_sfc");

    const Real rho_sfc = rho_zt[0];

    fill_multifab(hfx3, 0.0_rt);
    fill_multifab(qfx3, 0.0_rt);
    fill_multifab(tau13, 0.0_rt);
    fill_multifab(tau23, 0.0_rt);

    for (MFIter mfi(hfx3, false); mfi.isValid(); ++mfi) {
        const Box& bx = mfi.validbox();
        auto hfx = hfx3.array(mfi);
        auto qfx = qfx3.array(mfi);
        amrex::ParallelFor(bx, [=] AMREX_GPU_DEVICE (int i, int j, int k) noexcept
        {
            const int ic = i + nx * j;
            hfx(i,j,k) = (k == 0) ? rho_sfc * wthl_sfc[ic] : 0.0_rt;
            qfx(i,j,k) = (k == 0) ? rho_sfc * wqw_sfc[ic] : 0.0_rt;
        });
    }

    for (MFIter mfi(tau13, false); mfi.isValid(); ++mfi) {
        const Box& bx = mfi.validbox();
        auto arr = tau13.array(mfi);
        amrex::ParallelFor(bx, [=] AMREX_GPU_DEVICE (int i, int j, int k) noexcept
        {
            const int ic = amrex::min(i, nx - 1) + nx * j;
            arr(i,j,k) = (k == 0) ? rho_sfc * uw_sfc[ic] : 0.0_rt;
        });
    }

    for (MFIter mfi(tau23, false); mfi.isValid(); ++mfi) {
        const Box& bx = mfi.validbox();
        auto arr = tau23.array(mfi);
        amrex::ParallelFor(bx, [=] AMREX_GPU_DEVICE (int i, int j, int k) noexcept
        {
            const int jc = amrex::min(j, ny - 1);
            const int ic = i + nx * jc;
            arr(i,j,k) = (k == 0) ? rho_sfc * vw_sfc[ic] : 0.0_rt;
        });
    }
}

void
initialize_eddy_diffs (MultiFab& eddy_diffs,
                       const FixtureMap& fixture)
{
    const auto rho_zt = fixture_gpu_array<nz>(fixture, "rho_zt");
    const auto tkh = fixture_gpu_array<nz>(fixture, "tkh");

    fill_multifab(eddy_diffs, 0.0_rt);

    for (MFIter mfi(eddy_diffs, false); mfi.isValid(); ++mfi) {
        const Box& bx = mfi.validbox();
        auto arr = eddy_diffs.array(mfi);
        amrex::ParallelFor(bx, [=] AMREX_GPU_DEVICE (int i, int j, int k) noexcept
        {
            const Real coeff = rho_zt[k] * tkh[k];
            arr(i,j,k,EddyDiff::Mom_v) = coeff;
            arr(i,j,k,EddyDiff::Theta_v) = coeff;
            arr(i,j,k,EddyDiff::KE_v) = coeff;
            arr(i,j,k,EddyDiff::Scalar_v) = coeff;
            arr(i,j,k,EddyDiff::Q_v) = coeff;
            arr(i,j,k,EddyDiff::Turb_lengthscale) = 100.0_rt;
        });
    }
}

void
initialize_passive_scalar (MultiFab& cons)
{
    for (MFIter mfi(cons, false); mfi.isValid(); ++mfi) {
        const Box& bx = mfi.validbox();
        auto arr = cons.array(mfi);
        amrex::ParallelFor(bx, [=] AMREX_GPU_DEVICE (int i, int j, int k) noexcept
        {
            const Real rho = arr(i,j,k,Rho_comp);
            arr(i,j,k,RhoScalar_comp) = rho * (0.1_rt + 0.01_rt * k + 0.001_rt * i);
        });
    }
}

void
reset_rhs (amrex::Vector<MultiFab>& rhs)
{
    for (auto& mf : rhs) {
        fill_multifab(mf, 0.0_rt);
    }
}

void
expect_driver_diagnostics_stable (const ShocDriver& driver,
                                  Real max_mix_len)
{
    expect_multifab_finite(driver.native_diagnostics(), 0, EddyDiff::NumDiffs, "multistep native_diagnostics");
    expect_component_nonnegative(driver.native_diagnostics(), EddyDiff::Mom_v, "multistep Kmv");
    expect_component_nonnegative(driver.native_diagnostics(), EddyDiff::Theta_v, "multistep Khv");
    expect_component_nonnegative(driver.native_diagnostics(), EddyDiff::KE_v, "multistep KE_v");
    expect_component_nonnegative(driver.native_diagnostics(), EddyDiff::Scalar_v, "multistep Scalar_v");
    expect_component_nonnegative(driver.native_diagnostics(), EddyDiff::Q_v, "multistep Q_v");
    expect_component_nonnegative(driver.native_diagnostics(), EddyDiff::Turb_lengthscale, "multistep lengthscale");
    EXPECT_LE(driver.native_diagnostics().max(EddyDiff::Turb_lengthscale), max_mix_len + tol)
        << "multistep native lengthscale exceeds the geometry cap enforced by structure diagnostics";

    expect_multifab_finite(driver.pblh_diagnostics(), 0, 1, "multistep pblh");
    expect_component_nonnegative(driver.pblh_diagnostics(), 0, "multistep pblh");

    expect_multifab_finite(driver.shoc_cldfrac_diagnostics(), 0, 1, "multistep cloud fraction");
    expect_component_between(driver.shoc_cldfrac_diagnostics(), 0, 0.0_rt, 1.0_rt, "multistep cloud fraction");

    expect_multifab_finite(driver.shoc_ql_diagnostics(), 0, 1, "multistep ql");
    expect_component_nonnegative(driver.shoc_ql_diagnostics(), 0, "multistep ql");

    expect_multifab_finite(driver.w_sec_diagnostics(), 0, 1, "multistep w_sec");
    expect_component_nonnegative(driver.w_sec_diagnostics(), 0, "multistep w_sec");

    expect_multifab_finite(driver.wthv_sec_diagnostics(), 0, 1, "multistep wthv_sec");
    EXPECT_LT(component_max_abs(driver.wthv_sec_diagnostics(), 0), 10.0_rt)
        << "multistep wthv_sec leaves the bounded range already enforced in SHOC PDF property tests";
}
}

TEST(ShocPreprocess, SurfaceStressAveragesFaceValuesToColumnCenters)
{
    Box domain(IntVect(0,0,0), IntVect(1,1,1));
    amrex::RealBox real_box(0.0_rt, 0.0_rt, 0.0_rt,
                            200.0_rt, 200.0_rt, 200.0_rt);
    int is_periodic[AMREX_SPACEDIM] = {1, 1, 0};
    Geometry geom(domain, &real_box, amrex::CoordSys::cartesian, is_periodic);

    amrex::BoxArray ba(domain);
    ba.maxSize(IntVect(2, 2, 2));
    DistributionMapping dm(ba);

    MultiFab cons(ba, dm, NVAR_max, 0);
    amrex::BoxArray xba = amrex::convert(ba, IntVect(1,0,0));
    amrex::BoxArray yba = amrex::convert(ba, IntVect(0,1,0));
    amrex::BoxArray zba = amrex::convert(ba, IntVect(0,0,1));
    amrex::BoxArray xzba = amrex::convert(ba, IntVect(1,0,1));
    amrex::BoxArray yzba = amrex::convert(ba, IntVect(0,1,1));

    MultiFab xvel(xba, dm, 1, 0);
    MultiFab yvel(yba, dm, 1, 0);
    MultiFab zvel(zba, dm, 1, 0);
    MultiFab z_phys_nd(make_nodal_zphys_boxarray(ba), dm, 1, 0);
    MultiFab hfx3(ba, dm, 1, 0);
    MultiFab qfx3(ba, dm, 1, 0);
    MultiFab tau13(xzba, dm, 1, 0);
    MultiFab tau23(yzba, dm, 1, 0);
    MultiFab eddy_diffs(ba, dm, EddyDiff::NumDiffs, 0);

    initialize_basic_column_state(cons, zvel, z_phys_nd, hfx3, qfx3, eddy_diffs);
    initialize_face_velocities_with_seam_offsets(xvel, yvel, domain, 0.0_rt, 0.0_rt);
    initialize_face_stress_gradients(tau13, tau23);
    shoc_test::sync();

    ShocColumnWorkspace workspace;
    const MoistureComponentIndices moisture_indices(RhoQ1_comp, RhoQ2_comp);

    for (MFIter mfi(cons, false); mfi.isValid(); ++mfi) {
        const ShocColumnLayout active_layout = make_shoc_layout(mfi.validbox(), geom);
        workspace.ensure_capacity(active_layout, shoc_test::test_arena(), shoc::InitRunOn::Host);

        shoc_test::run_and_sync([&] {
            ShocPreprocess::fill_columns(workspace.col, mfi, cons, xvel, yvel, zvel,
                                         &hfx3, &qfx3, &tau13, &tau23, z_phys_nd, geom,
                                         moisture_indices);
        });

        const auto tauu = workspace.col.surf_tau_u.const_array();
        const auto tauv = workspace.col.surf_tau_v.const_array();
        const auto layout = workspace.col.layout;

        for (int j = 0; j < layout.ny; ++j) {
            for (int i = 0; i < layout.nx; ++i) {
                const int ic = shoc_column_index(layout, layout.imin + i, layout.jmin + j);
                EXPECT_NEAR(tauu(ic,0,0), 10.5_rt + static_cast<Real>(i), tol);
                EXPECT_NEAR(tauv(ic,0,0), 20.5_rt + static_cast<Real>(j), tol);
            }
        }
    }
}

TEST(ShocPreprocess, UsesFourNodeAveragedTerrainGeometry)
{
    Box domain(IntVect(0,0,0), IntVect(1,1,2));
    amrex::RealBox real_box(0.0_rt, 0.0_rt, 0.0_rt,
                            200.0_rt, 200.0_rt, 300.0_rt);
    int is_periodic[AMREX_SPACEDIM] = {0, 0, 0};
    Geometry geom(domain, &real_box, amrex::CoordSys::cartesian, is_periodic);

    amrex::BoxArray ba(domain);
    ba.maxSize(IntVect(2, 2, 3));
    DistributionMapping dm(ba);

    MultiFab cons(ba, dm, NVAR_max, 0);
    amrex::BoxArray xba = amrex::convert(ba, IntVect(1,0,0));
    amrex::BoxArray yba = amrex::convert(ba, IntVect(0,1,0));
    amrex::BoxArray zba = amrex::convert(ba, IntVect(0,0,1));
    amrex::BoxArray xzba = amrex::convert(ba, IntVect(1,0,1));
    amrex::BoxArray yzba = amrex::convert(ba, IntVect(0,1,1));
    amrex::BoxArray zphys_ba = make_nodal_zphys_boxarray(ba);

    MultiFab xvel(xba, dm, 1, 0);
    MultiFab yvel(yba, dm, 1, 0);
    MultiFab zvel(zba, dm, 1, 0);
    MultiFab z_phys_nd(zphys_ba, dm, 1, 0);
    MultiFab hfx3(ba, dm, 1, 0);
    MultiFab qfx3(ba, dm, 1, 0);
    MultiFab tau13(xzba, dm, 1, 0);
    MultiFab tau23(yzba, dm, 1, 0);
    MultiFab eddy_diffs(ba, dm, EddyDiff::NumDiffs, 0);

    initialize_basic_column_state(cons, zvel, z_phys_nd, hfx3, qfx3, eddy_diffs);
    initialize_face_velocities_with_seam_offsets(xvel, yvel, domain, 0.0_rt, 0.0_rt);
    initialize_face_stress_gradients(tau13, tau23);

    initialize_terrain_geometry(z_phys_nd);
    shoc_test::sync();

    ShocColumnWorkspace workspace;
    const MoistureComponentIndices moisture_indices(RhoQ1_comp, RhoQ2_comp);

    for (MFIter mfi(cons, false); mfi.isValid(); ++mfi) {
        const ShocColumnLayout active_layout = make_shoc_layout(mfi.validbox(), geom);
        workspace.ensure_capacity(active_layout, shoc_test::test_arena(), shoc::InitRunOn::Host);

        shoc_test::run_and_sync([&] {
            ShocPreprocess::fill_columns(workspace.col, mfi, cons, xvel, yvel, zvel,
                                         &hfx3, &qfx3, &tau13, &tau23, z_phys_nd, geom,
                                         moisture_indices);
        });

        const auto zi = workspace.col.zi.const_array();
        const auto zt = workspace.col.zt.const_array();
        const auto dz = workspace.col.dz.const_array();
        const auto dse = workspace.col.host_dse.const_array();
        const auto tabs = workspace.col.tabs.const_array();
        const auto zarr = z_phys_nd.const_array(mfi);
        const auto layout = workspace.col.layout;

        for (int j = 0; j < layout.ny; ++j) {
            for (int i = 0; i < layout.nx; ++i) {
                const int ic = shoc_column_index(layout, layout.imin + i, layout.jmin + j);
                for (int kk = 0; kk < layout.nlev; ++kk) {
                    const int k = layout.kmin + kk;
                    const Real zlo = 0.25_rt * (zarr(i,  j,  k) +
                                                zarr(i+1,j,  k) +
                                                zarr(i,  j+1,k) +
                                                zarr(i+1,j+1,k));
                    const Real zhi = 0.25_rt * (zarr(i,  j,  k+1) +
                                                zarr(i+1,j,  k+1) +
                                                zarr(i,  j+1,k+1) +
                                                zarr(i+1,j+1,k+1));
                    const Real zc = 0.5_rt * (zlo + zhi);
                    const Real expected_dse = Cp_d * tabs(ic,kk,0) + CONST_GRAV * zc;

                    EXPECT_NEAR(zi(ic,kk,0), zlo, tol);
                    EXPECT_NEAR(zt(ic,kk,0), zc, tol);
                    EXPECT_NEAR(dz(ic,kk,0), zhi - zlo, tol);
                    EXPECT_NEAR(dse(ic,kk,0), expected_dse, tol);
                }
            }
        }
    }
}

TEST(ShocDriver, FaceSyncKeepsSeamCopiesAlignedAfterStateUpdate)
{
    Box domain(IntVect(0,0,0), IntVect(3,3,4));
    amrex::RealBox real_box(0.0_rt, 0.0_rt, 0.0_rt,
                            400.0_rt, 400.0_rt, 500.0_rt);
    int is_periodic[AMREX_SPACEDIM] = {1, 1, 0};
    Geometry geom(domain, &real_box, amrex::CoordSys::cartesian, is_periodic);

    amrex::BoxArray ba(domain);
    ba.maxSize(IntVect(2, 2, nz));
    DistributionMapping dm(ba);

    MultiFab cons(ba, dm, NVAR_max, 0);
    amrex::BoxArray xba = amrex::convert(ba, IntVect(1,0,0));
    amrex::BoxArray yba = amrex::convert(ba, IntVect(0,1,0));
    amrex::BoxArray zba = amrex::convert(ba, IntVect(0,0,1));
    amrex::BoxArray xzba = amrex::convert(ba, IntVect(1,0,1));
    amrex::BoxArray yzba = amrex::convert(ba, IntVect(0,1,1));

    MultiFab xvel(xba, dm, 1, 0);
    MultiFab yvel(yba, dm, 1, 0);
    MultiFab zvel(zba, dm, 1, 0);
    MultiFab z_phys_nd(make_nodal_zphys_boxarray(ba), dm, 1, 0);
    MultiFab hfx3(ba, dm, 1, 0);
    MultiFab qfx3(ba, dm, 1, 0);
    MultiFab tau13(xzba, dm, 1, 0);
    MultiFab tau23(yzba, dm, 1, 0);
    MultiFab eddy_diffs(ba, dm, EddyDiff::NumDiffs, 0);

    initialize_basic_column_state(cons, zvel, z_phys_nd, hfx3, qfx3, eddy_diffs);
    initialize_face_velocities_with_seam_offsets(xvel, yvel, domain, 0.25_rt, 0.15_rt);
    initialize_face_stresses_with_seam_offsets(tau13, tau23, domain, 0.02_rt, 0.01_rt);
    shoc_test::sync();

    SolverChoice solver_choice;
    solver_choice.moisture_type = MoistureType::Morrison;
    solver_choice.moisture_indices = MoistureComponentIndices(RhoQ1_comp, RhoQ2_comp);

    ShocDriver driver(0, solver_choice);
    EXPECT_TRUE(driver.uses_state_update());
    EXPECT_FALSE(driver.uses_host_diffusion());

    shoc_test::run_and_sync([&] {
        driver.advance(cons, xvel, yvel, zvel,
                       &tau13, &tau23, &hfx3, &qfx3, &eddy_diffs,
                       z_phys_nd, geom, 5.0_rt);
    });

    const int seam_x = domain.smallEnd(0) + 2;
    const int seam_y = domain.smallEnd(1) + 2;

    expect_face_overlap_matches(xvel, domain, 0, seam_x, 0, domain.smallEnd(1), 0,
                                "xvel");
    expect_face_overlap_matches(yvel, domain, 1, seam_y, domain.smallEnd(0), 0, 0,
                                "yvel");
    expect_face_overlap_matches(tau13, domain, 0, seam_x, 0, domain.smallEnd(1), 0,
                                "tau13");
    expect_face_overlap_matches(tau23, domain, 1, seam_y, domain.smallEnd(0), 0, 0,
                                "tau23");
}

TEST(ShocDriver, AdvanceAppliesStateUpdateAndProducesFiniteDiagnostics)
{
    Box domain(IntVect(0,0,0), IntVect(nx-1, ny-1, nz-1));
    amrex::RealBox real_box(0.0_rt, 0.0_rt, 0.0_rt,
                            500.0_rt, 100.0_rt, 900.0_rt);
    int is_periodic[AMREX_SPACEDIM] = {1, 1, 0};
    Geometry geom(domain, &real_box, amrex::CoordSys::cartesian, is_periodic);
    const FixtureMap fixture = load_driver_fixture();

    amrex::BoxArray ba(domain);
    // The current SHOC driver path maps each MFIter tile to a full-height
    // column workspace, so keep tiling in x only here.
    ba.maxSize(IntVect(3, ny, nz));
    DistributionMapping dm(ba);

    MultiFab cons(ba, dm, NVAR_max, 0);
    amrex::BoxArray xba = amrex::convert(ba, IntVect(1,0,0));
    amrex::BoxArray yba = amrex::convert(ba, IntVect(0,1,0));
    amrex::BoxArray zba = amrex::convert(ba, IntVect(0,0,1));
    amrex::BoxArray xzba = amrex::convert(ba, IntVect(1,0,1));
    amrex::BoxArray yzba = amrex::convert(ba, IntVect(0,1,1));

    MultiFab xvel(xba, dm, 1, 0);
    MultiFab yvel(yba, dm, 1, 0);
    MultiFab zvel(zba, dm, 1, 0);
    MultiFab z_phys_nd(make_nodal_zphys_boxarray(ba), dm, 1, 0);
    MultiFab hfx3(ba, dm, 1, 0);
    MultiFab qfx3(ba, dm, 1, 0);
    MultiFab tau13(xzba, dm, 1, 0);
    MultiFab tau23(yzba, dm, 1, 0);
    MultiFab eddy_diffs(ba, dm, EddyDiff::NumDiffs, 0);

    initialize_state(cons, fixture);
    initialize_velocity(xvel, yvel, zvel, fixture);
    initialize_geometry(z_phys_nd, fixture);
    initialize_surface_fluxes(hfx3, qfx3, tau13, tau23, fixture);
    initialize_eddy_diffs(eddy_diffs, fixture);
    shoc_test::sync();

    expect_multifab_finite(cons, Rho_comp, 1, "cons rho input");
    expect_multifab_finite(cons, RhoTheta_comp, 1, "cons rho-theta input");
    expect_multifab_finite(cons, RhoQ1_comp, 1, "cons rho-qv input");

    SolverChoice solver_choice;
    solver_choice.moisture_type = MoistureType::None;

    ShocDriver driver(0, solver_choice);
    EXPECT_TRUE(driver.uses_state_update());
    EXPECT_FALSE(driver.uses_host_diffusion());
    EXPECT_TRUE(driver.owns_scalar_surface_fluxes());
    EXPECT_TRUE(driver.uses_momentum_host_diffusion());
    const Real dt = 10.0_rt;
    const Real max_mix_len = std::sqrt(geom.CellSizeArray()[0] * geom.CellSizeArray()[1]);

    MultiFab cons_before(ba, dm, NVAR_max, 0);
    MultiFab xvel_before(xba, dm, 1, 0);
    MultiFab yvel_before(yba, dm, 1, 0);
    MultiFab zvel_before(zba, dm, 1, 0);
    MultiFab scalar_before(ba, dm, 1, 0);
    MultiFab hfx3_before(ba, dm, 1, 0);
    MultiFab qfx3_before(ba, dm, 1, 0);
    MultiFab tau13_before(xzba, dm, 1, 0);
    MultiFab tau23_before(yzba, dm, 1, 0);
    MultiFab::Copy(cons_before, cons, 0, 0, NVAR_max, 0);
    MultiFab::Copy(xvel_before, xvel, 0, 0, 1, 0);
    MultiFab::Copy(yvel_before, yvel, 0, 0, 1, 0);
    MultiFab::Copy(zvel_before, zvel, 0, 0, 1, 0);
    MultiFab::Copy(scalar_before, cons, RhoScalar_comp, 0, 1, 0);
    MultiFab::Copy(hfx3_before, hfx3, 0, 0, 1, 0);
    MultiFab::Copy(qfx3_before, qfx3, 0, 0, 1, 0);
    MultiFab::Copy(tau13_before, tau13, 0, 0, 1, 0);
    MultiFab::Copy(tau23_before, tau23, 0, 0, 1, 0);

    shoc_test::run_and_sync([&] {
        driver.advance(cons, xvel, yvel, zvel,
                       &tau13, &tau23, &hfx3, &qfx3, &eddy_diffs,
                       z_phys_nd, geom, dt);
    });
    EXPECT_TRUE(driver.has_native_diagnostics());
    EXPECT_TRUE(driver.uses_state_update());

    expect_multifab_finite(driver.native_diagnostics(), 0, EddyDiff::NumDiffs, "native_diagnostics");
    expect_component_nonnegative(driver.native_diagnostics(), EddyDiff::Mom_v, "native Kmv");
    expect_component_nonnegative(driver.native_diagnostics(), EddyDiff::Theta_v, "native Khv");
    expect_component_nonnegative(driver.native_diagnostics(), EddyDiff::KE_v, "native KE_v");
    expect_component_nonnegative(driver.native_diagnostics(), EddyDiff::Scalar_v, "native Scalar_v");
    expect_component_nonnegative(driver.native_diagnostics(), EddyDiff::Q_v, "native Q_v");
    expect_component_nonnegative(driver.native_diagnostics(), EddyDiff::Turb_lengthscale, "native lengthscale");
    EXPECT_LE(driver.native_diagnostics().max(EddyDiff::Turb_lengthscale), max_mix_len + tol)
        << "native lengthscale exceeds the geometry cap enforced by structure diagnostics";

    expect_multifab_finite(driver.pblh_diagnostics(), 0, 1, "pblh");
    expect_component_nonnegative(driver.pblh_diagnostics(), 0, "pblh");

    expect_multifab_finite(driver.shoc_cldfrac_diagnostics(), 0, 1, "cloud fraction");
    expect_component_between(driver.shoc_cldfrac_diagnostics(), 0, 0.0_rt, 1.0_rt, "cloud fraction");

    expect_multifab_finite(driver.shoc_ql_diagnostics(), 0, 1, "ql");
    expect_component_nonnegative(driver.shoc_ql_diagnostics(), 0, "ql");

    expect_multifab_finite(driver.w_sec_diagnostics(), 0, 1, "w_sec");
    expect_component_nonnegative(driver.w_sec_diagnostics(), 0, "w_sec");
    expect_multifab_finite(driver.wthv_sec_diagnostics(), 0, 1, "wthv_sec");
    EXPECT_LT(component_max_abs(driver.wthv_sec_diagnostics(), 0), 10.0_rt)
        << "wthv_sec leaves the bounded range already enforced in SHOC PDF property tests";

    EXPECT_EQ(component_max_abs_difference(cons, Rho_comp, cons_before, Rho_comp), 0.0_rt);
    EXPECT_GT(component_max_abs_difference(cons, RhoTheta_comp, cons_before, RhoTheta_comp), 0.0_rt);
    expect_component_nonnegative(cons, RhoQ1_comp, "state-update qv");
    expect_component_nonnegative(cons, RhoKE_comp, "state-update tke");
    EXPECT_EQ(component_max_abs_difference(xvel, 0, xvel_before, 0), 0.0_rt);
    EXPECT_EQ(component_max_abs_difference(yvel, 0, yvel_before, 0), 0.0_rt);
    EXPECT_EQ(component_max_abs_difference(zvel, 0, zvel_before, 0), 0.0_rt);
    EXPECT_EQ(component_max_abs_difference(cons, RhoScalar_comp, scalar_before, 0), 0.0_rt);

    shoc_test::run_and_sync([&] {
        driver.set_eddy_diffs();
    });

    expect_multifab_finite(eddy_diffs, 0, EddyDiff::NumDiffs, "eddy_diffs");
    expect_component_minmax_match(eddy_diffs, driver.native_diagnostics(),
                                  EddyDiff::Mom_v, "eddy_diffs Mom_v writeback");
    expect_component_between(eddy_diffs, EddyDiff::Theta_v, 0.0_rt, 0.0_rt, "eddy_diffs Theta_v cleared");
    expect_component_between(eddy_diffs, EddyDiff::KE_v, 0.0_rt, 0.0_rt, "eddy_diffs KE_v cleared");
    expect_component_between(eddy_diffs, EddyDiff::Scalar_v, 0.0_rt, 0.0_rt, "eddy_diffs Scalar_v cleared");
    expect_component_between(eddy_diffs, EddyDiff::Q_v, 0.0_rt, 0.0_rt, "eddy_diffs Q_v cleared");
    expect_component_nonnegative(eddy_diffs, EddyDiff::Turb_lengthscale, "eddy_diffs lengthscale");
    expect_component_minmax_match(eddy_diffs, driver.native_diagnostics(),
                                  EddyDiff::Turb_lengthscale, "lengthscale writeback");

    shoc_test::run_and_sync([&] {
        driver.set_diff_stresses();
    });

    expect_component_between(hfx3, 0, 0.0_rt, 0.0_rt, "hfx3 cleared");
    expect_component_between(qfx3, 0, 0.0_rt, 0.0_rt, "qfx3 cleared");
    expect_component_minmax_match(tau13, tau13_before, 0, "tau13 preserved for host momentum diffusion");
    expect_component_minmax_match(tau23, tau23_before, 0, "tau23 preserved for host momentum diffusion");

    amrex::Vector<MultiFab> rhs(IntVars::NumTypes);
    rhs[IntVars::cons].define(ba, dm, NVAR_max, 0);
    rhs[IntVars::xmom].define(xba, dm, 1, 0);
    rhs[IntVars::ymom].define(yba, dm, 1, 0);
    rhs[IntVars::zmom].define(zba, dm, 1, 0);
    for (auto& mf : rhs) {
        fill_multifab(mf, 0.0_rt);
    }

    shoc_test::run_and_sync([&] {
        driver.add_fast_tend(rhs);
    });

    shoc_test::run_and_sync([&] {
        for (MFIter mfi(rhs[IntVars::cons], false); mfi.isValid(); ++mfi) {
            driver.add_slow_tend(mfi, mfi.validbox(), rhs[IntVars::cons].array(mfi));
        }
    });

    expect_component_between(rhs[IntVars::cons], RhoTheta_comp, 0.0_rt, 0.0_rt, "theta rhs");
    expect_component_between(rhs[IntVars::cons], RhoQ1_comp, 0.0_rt, 0.0_rt, "qv rhs");
    expect_component_between(rhs[IntVars::cons], RhoKE_comp, 0.0_rt, 0.0_rt, "tke rhs");
    expect_component_between(rhs[IntVars::xmom], 0, 0.0_rt, 0.0_rt, "xmom rhs");
    expect_component_between(rhs[IntVars::ymom], 0, 0.0_rt, 0.0_rt, "ymom rhs");
    expect_component_between(rhs[IntVars::zmom], 0, 0.0_rt, 0.0_rt, "zmom rhs");
}

TEST(ShocDriver, MomentumTransportNoneLeavesVelocitiesUntouchedAndSkipsMomentumExports)
{
    ScopedParmParseString momentum_transport("erf.shoc", "momentum_transport", "none");

    Box domain(IntVect(0,0,0), IntVect(nx-1, ny-1, nz-1));
    amrex::RealBox real_box(0.0_rt, 0.0_rt, 0.0_rt,
                            500.0_rt, 100.0_rt, 900.0_rt);
    int is_periodic[AMREX_SPACEDIM] = {1, 1, 0};
    Geometry geom(domain, &real_box, amrex::CoordSys::cartesian, is_periodic);
    const FixtureMap fixture = load_driver_fixture();

    amrex::BoxArray ba(domain);
    ba.maxSize(IntVect(3, ny, nz));
    DistributionMapping dm(ba);

    MultiFab cons(ba, dm, NVAR_max, 0);
    amrex::BoxArray xba = amrex::convert(ba, IntVect(1,0,0));
    amrex::BoxArray yba = amrex::convert(ba, IntVect(0,1,0));
    amrex::BoxArray zba = amrex::convert(ba, IntVect(0,0,1));
    amrex::BoxArray xzba = amrex::convert(ba, IntVect(1,0,1));
    amrex::BoxArray yzba = amrex::convert(ba, IntVect(0,1,1));

    MultiFab xvel(xba, dm, 1, 0);
    MultiFab yvel(yba, dm, 1, 0);
    MultiFab zvel(zba, dm, 1, 0);
    MultiFab z_phys_nd(make_nodal_zphys_boxarray(ba), dm, 1, 0);
    MultiFab hfx3(ba, dm, 1, 0);
    MultiFab qfx3(ba, dm, 1, 0);
    MultiFab tau13(xzba, dm, 1, 0);
    MultiFab tau23(yzba, dm, 1, 0);
    MultiFab eddy_diffs(ba, dm, EddyDiff::NumDiffs, 0);

    initialize_state(cons, fixture);
    initialize_velocity(xvel, yvel, zvel, fixture);
    initialize_geometry(z_phys_nd, fixture);
    initialize_surface_fluxes(hfx3, qfx3, tau13, tau23, fixture);
    initialize_eddy_diffs(eddy_diffs, fixture);
    shoc_test::sync();

    MultiFab xvel_before(xba, dm, 1, 0);
    MultiFab yvel_before(yba, dm, 1, 0);
    MultiFab::Copy(xvel_before, xvel, 0, 0, 1, 0);
    MultiFab::Copy(yvel_before, yvel, 0, 0, 1, 0);

    SolverChoice solver_choice;
    solver_choice.moisture_type = MoistureType::None;

    ShocDriver driver(0, solver_choice);
    EXPECT_TRUE(driver.uses_state_update());
    EXPECT_TRUE(driver.owns_scalar_surface_fluxes());
    EXPECT_TRUE(driver.disables_momentum_transport());
    EXPECT_FALSE(driver.uses_momentum_host_diffusion());
    EXPECT_FALSE(driver.uses_momentum_state_update());

    shoc_test::run_and_sync([&] {
        driver.advance(cons, xvel, yvel, zvel,
                       &tau13, &tau23, &hfx3, &qfx3, &eddy_diffs,
                       z_phys_nd, geom, 10.0_rt);
        driver.set_eddy_diffs();
        driver.set_diff_stresses();
    });

    EXPECT_EQ(component_max_abs_difference(xvel, 0, xvel_before, 0), 0.0_rt);
    EXPECT_EQ(component_max_abs_difference(yvel, 0, yvel_before, 0), 0.0_rt);
    expect_component_between(eddy_diffs, EddyDiff::Mom_v, 0.0_rt, 0.0_rt, "none momentum export");
    expect_component_minmax_match(eddy_diffs, driver.native_diagnostics(),
                                  EddyDiff::Turb_lengthscale, "none lengthscale writeback");
    expect_component_between(hfx3, 0, 0.0_rt, 0.0_rt, "none hfx3 cleared");
    expect_component_between(qfx3, 0, 0.0_rt, 0.0_rt, "none qfx3 cleared");
    expect_component_between(tau13, 0, 0.0_rt, 0.0_rt, "none tau13 cleared");
    expect_component_between(tau23, 0, 0.0_rt, 0.0_rt, "none tau23 cleared");
}

TEST(ShocDriver, MomentumTransportStateUpdateUpdatesFaceVelocitiesWhenRequested)
{
    ScopedParmParseString momentum_transport("erf.shoc", "momentum_transport", "state_update");

    Box domain(IntVect(0,0,0), IntVect(nx-1, ny-1, nz-1));
    amrex::RealBox real_box(0.0_rt, 0.0_rt, 0.0_rt,
                            500.0_rt, 100.0_rt, 900.0_rt);
    int is_periodic[AMREX_SPACEDIM] = {1, 1, 0};
    Geometry geom(domain, &real_box, amrex::CoordSys::cartesian, is_periodic);
    const FixtureMap fixture = load_driver_fixture();

    amrex::BoxArray ba(domain);
    ba.maxSize(IntVect(3, ny, nz));
    DistributionMapping dm(ba);

    MultiFab cons(ba, dm, NVAR_max, 0);
    amrex::BoxArray xba = amrex::convert(ba, IntVect(1,0,0));
    amrex::BoxArray yba = amrex::convert(ba, IntVect(0,1,0));
    amrex::BoxArray zba = amrex::convert(ba, IntVect(0,0,1));
    amrex::BoxArray xzba = amrex::convert(ba, IntVect(1,0,1));
    amrex::BoxArray yzba = amrex::convert(ba, IntVect(0,1,1));

    MultiFab xvel(xba, dm, 1, 0);
    MultiFab yvel(yba, dm, 1, 0);
    MultiFab zvel(zba, dm, 1, 0);
    MultiFab z_phys_nd(make_nodal_zphys_boxarray(ba), dm, 1, 0);
    MultiFab hfx3(ba, dm, 1, 0);
    MultiFab qfx3(ba, dm, 1, 0);
    MultiFab tau13(xzba, dm, 1, 0);
    MultiFab tau23(yzba, dm, 1, 0);
    MultiFab eddy_diffs(ba, dm, EddyDiff::NumDiffs, 0);

    initialize_state(cons, fixture);
    initialize_velocity(xvel, yvel, zvel, fixture);
    initialize_geometry(z_phys_nd, fixture);
    initialize_surface_fluxes(hfx3, qfx3, tau13, tau23, fixture);
    initialize_eddy_diffs(eddy_diffs, fixture);
    shoc_test::sync();

    MultiFab xvel_before(xba, dm, 1, 0);
    MultiFab yvel_before(yba, dm, 1, 0);
    MultiFab tau13_before(xzba, dm, 1, 0);
    MultiFab tau23_before(yzba, dm, 1, 0);
    MultiFab::Copy(xvel_before, xvel, 0, 0, 1, 0);
    MultiFab::Copy(yvel_before, yvel, 0, 0, 1, 0);
    MultiFab::Copy(tau13_before, tau13, 0, 0, 1, 0);
    MultiFab::Copy(tau23_before, tau23, 0, 0, 1, 0);

    SolverChoice solver_choice;
    solver_choice.moisture_type = MoistureType::None;

    ShocDriver driver(0, solver_choice);
    EXPECT_TRUE(driver.uses_state_update());
    EXPECT_TRUE(driver.owns_scalar_surface_fluxes());
    EXPECT_TRUE(driver.uses_momentum_state_update());
    EXPECT_FALSE(driver.uses_momentum_host_diffusion());

    shoc_test::run_and_sync([&] {
        driver.advance(cons, xvel, yvel, zvel,
                       &tau13, &tau23, &hfx3, &qfx3, &eddy_diffs,
                       z_phys_nd, geom, 10.0_rt);
        driver.set_diff_stresses();
    });

    EXPECT_GT(component_max_abs_difference(xvel, 0, xvel_before, 0), 0.0_rt);
    EXPECT_GT(component_max_abs_difference(yvel, 0, yvel_before, 0), 0.0_rt);
    expect_component_between(hfx3, 0, 0.0_rt, 0.0_rt, "state_update hfx3 cleared");
    expect_component_between(qfx3, 0, 0.0_rt, 0.0_rt, "state_update qfx3 cleared");
    expect_component_between(tau13, 0, 0.0_rt, 0.0_rt, "state_update tau13 cleared");
    expect_component_between(tau23, 0, 0.0_rt, 0.0_rt, "state_update tau23 cleared");
}

TEST(ShocDriver, FullHeightBoxPredicateRejectsVerticallySplitBoxes)
{
    Box domain(IntVect(0,0,0), IntVect(nx-1, ny-1, nz-1));
    amrex::BoxArray ba(domain);
    ba.maxSize(3);

    EXPECT_FALSE(shoc_boxarray_spans_full_height(ba, domain));

    amrex::BoxArray full_height(domain);
    full_height.maxSize(IntVect(3, ny, nz));
    EXPECT_TRUE(shoc_boxarray_spans_full_height(full_height, domain));
}

TEST(ShocDriver, PassiveScalarIsNotModifiedByNativeShoc)
{
    // Native SHOC applies its coupled heat, moisture, momentum, and TKE
    // increment directly to the state. It should not modify passive scalar
    // state components.
    ASSERT_LT(RhoScalar_comp, NVAR_max);

    Box domain(IntVect(0,0,0), IntVect(nx-1, ny-1, nz-1));
    amrex::RealBox real_box(0.0_rt, 0.0_rt, 0.0_rt,
                            500.0_rt, 100.0_rt, 900.0_rt);
    int is_periodic[AMREX_SPACEDIM] = {1, 1, 0};
    Geometry geom(domain, &real_box, amrex::CoordSys::cartesian, is_periodic);
    const FixtureMap fixture = load_driver_fixture();

    amrex::BoxArray ba(domain);
    ba.maxSize(IntVect(3, ny, nz));
    DistributionMapping dm(ba);

    MultiFab cons(ba, dm, NVAR_max, 0);
    amrex::BoxArray xba = amrex::convert(ba, IntVect(1,0,0));
    amrex::BoxArray yba = amrex::convert(ba, IntVect(0,1,0));
    amrex::BoxArray zba = amrex::convert(ba, IntVect(0,0,1));
    amrex::BoxArray xzba = amrex::convert(ba, IntVect(1,0,1));
    amrex::BoxArray yzba = amrex::convert(ba, IntVect(0,1,1));

    MultiFab xvel(xba, dm, 1, 0);
    MultiFab yvel(yba, dm, 1, 0);
    MultiFab zvel(zba, dm, 1, 0);
    MultiFab z_phys_nd(make_nodal_zphys_boxarray(ba), dm, 1, 0);
    MultiFab hfx3(ba, dm, 1, 0);
    MultiFab qfx3(ba, dm, 1, 0);
    MultiFab tau13(xzba, dm, 1, 0);
    MultiFab tau23(yzba, dm, 1, 0);
    MultiFab eddy_diffs(ba, dm, EddyDiff::NumDiffs, 0);

    initialize_state(cons, fixture);
    initialize_velocity(xvel, yvel, zvel, fixture);
    initialize_geometry(z_phys_nd, fixture);
    initialize_surface_fluxes(hfx3, qfx3, tau13, tau23, fixture);
    initialize_eddy_diffs(eddy_diffs, fixture);
    initialize_passive_scalar(cons);
    shoc_test::sync();

    MultiFab scalar_before(ba, dm, 1, 0);
    MultiFab::Copy(scalar_before, cons, RhoScalar_comp, 0, 1, 0);

    SolverChoice solver_choice;
    solver_choice.moisture_type = MoistureType::Morrison;

    ShocDriver driver(0, solver_choice);
    const Real dt = 10.0_rt;

    shoc_test::run_and_sync([&] {
        driver.advance(cons, xvel, yvel, zvel,
                       &tau13, &tau23, &hfx3, &qfx3, &eddy_diffs,
                       z_phys_nd, geom, dt);
    });

    EXPECT_EQ(component_max_abs_difference(cons, RhoScalar_comp, scalar_before, 0), 0.0_rt);

    amrex::Vector<MultiFab> rhs(IntVars::NumTypes);
    rhs[IntVars::cons].define(ba, dm, NVAR_max, 0);
    rhs[IntVars::xmom].define(xba, dm, 1, 0);
    rhs[IntVars::ymom].define(yba, dm, 1, 0);
    rhs[IntVars::zmom].define(zba, dm, 1, 0);
    reset_rhs(rhs);

    shoc_test::run_and_sync([&] {
        driver.add_fast_tend(rhs);
        for (MFIter mfi(rhs[IntVars::cons], false); mfi.isValid(); ++mfi) {
            driver.add_slow_tend(mfi, mfi.validbox(), rhs[IntVars::cons].array(mfi));
        }
    });

    EXPECT_EQ(component_max_abs(rhs[IntVars::cons], RhoScalar_comp), 0.0_rt);
}

TEST(ShocDriver, SecondAdvanceIgnoresHostDiffSeedsAfterCarryStateIsEstablished)
{
    Box domain(IntVect(0,0,0), IntVect(nx-1, ny-1, nz-1));
    amrex::RealBox real_box(0.0_rt, 0.0_rt, 0.0_rt,
                            500.0_rt, 100.0_rt, 900.0_rt);
    int is_periodic[AMREX_SPACEDIM] = {1, 1, 0};
    Geometry geom(domain, &real_box, amrex::CoordSys::cartesian, is_periodic);
    const FixtureMap fixture = load_driver_fixture();

    amrex::BoxArray ba(domain);
    // The current SHOC driver path maps each MFIter tile to a full-height
    // column workspace, so keep tiling in x only here.
    ba.maxSize(IntVect(3, ny, nz));
    DistributionMapping dm(ba);

    MultiFab cons_a(ba, dm, NVAR_max, 0);
    MultiFab cons_b(ba, dm, NVAR_max, 0);
    amrex::BoxArray xba = amrex::convert(ba, IntVect(1,0,0));
    amrex::BoxArray yba = amrex::convert(ba, IntVect(0,1,0));
    amrex::BoxArray zba = amrex::convert(ba, IntVect(0,0,1));
    amrex::BoxArray xzba = amrex::convert(ba, IntVect(1,0,1));
    amrex::BoxArray yzba = amrex::convert(ba, IntVect(0,1,1));

    MultiFab xvel_a(xba, dm, 1, 0);
    MultiFab xvel_b(xba, dm, 1, 0);
    MultiFab yvel_a(yba, dm, 1, 0);
    MultiFab yvel_b(yba, dm, 1, 0);
    MultiFab zvel_a(zba, dm, 1, 0);
    MultiFab zvel_b(zba, dm, 1, 0);
    MultiFab z_phys_nd_a(make_nodal_zphys_boxarray(ba), dm, 1, 0);
    MultiFab z_phys_nd_b(make_nodal_zphys_boxarray(ba), dm, 1, 0);
    MultiFab hfx3_a(ba, dm, 1, 0);
    MultiFab hfx3_b(ba, dm, 1, 0);
    MultiFab qfx3_a(ba, dm, 1, 0);
    MultiFab qfx3_b(ba, dm, 1, 0);
    MultiFab tau13_a(xzba, dm, 1, 0);
    MultiFab tau13_b(xzba, dm, 1, 0);
    MultiFab tau23_a(yzba, dm, 1, 0);
    MultiFab tau23_b(yzba, dm, 1, 0);
    MultiFab eddy_diffs_a(ba, dm, EddyDiff::NumDiffs, 0);
    MultiFab eddy_diffs_b(ba, dm, EddyDiff::NumDiffs, 0);

    initialize_state(cons_a, fixture);
    initialize_state(cons_b, fixture);
    initialize_velocity(xvel_a, yvel_a, zvel_a, fixture);
    initialize_velocity(xvel_b, yvel_b, zvel_b, fixture);
    initialize_geometry(z_phys_nd_a, fixture);
    initialize_geometry(z_phys_nd_b, fixture);
    initialize_surface_fluxes(hfx3_a, qfx3_a, tau13_a, tau23_a, fixture);
    initialize_surface_fluxes(hfx3_b, qfx3_b, tau13_b, tau23_b, fixture);
    initialize_eddy_diffs(eddy_diffs_a, fixture);
    initialize_eddy_diffs(eddy_diffs_b, fixture);
    shoc_test::sync();

    SolverChoice solver_choice;
    solver_choice.moisture_type = MoistureType::None;

    ShocDriver driver_a(0, solver_choice);
    ShocDriver driver_b(0, solver_choice);

    shoc_test::run_and_sync([&] {
        driver_a.advance(cons_a, xvel_a, yvel_a, zvel_a,
                         &tau13_a, &tau23_a, &hfx3_a, &qfx3_a, &eddy_diffs_a,
                         z_phys_nd_a, geom, 10.0_rt);
        driver_b.advance(cons_b, xvel_b, yvel_b, zvel_b,
                         &tau13_b, &tau23_b, &hfx3_b, &qfx3_b, &eddy_diffs_b,
                         z_phys_nd_b, geom, 10.0_rt);
    });

    fill_multifab(eddy_diffs_a, 0.0_rt);
    set_eddy_diff_sentinel(eddy_diffs_b,
                           1.0e6_rt,
                           2.0e6_rt,
                           3.0e6_rt,
                           4.0e6_rt,
                           5.0e6_rt,
                           9.0e3_rt);
    shoc_test::sync();

    shoc_test::run_and_sync([&] {
        driver_a.advance(cons_a, xvel_a, yvel_a, zvel_a,
                         &tau13_a, &tau23_a, &hfx3_a, &qfx3_a, &eddy_diffs_a,
                         z_phys_nd_a, geom, 10.0_rt);
        driver_b.advance(cons_b, xvel_b, yvel_b, zvel_b,
                         &tau13_b, &tau23_b, &hfx3_b, &qfx3_b, &eddy_diffs_b,
                         z_phys_nd_b, geom, 10.0_rt);
    });

    for (int comp = 0; comp < EddyDiff::NumDiffs; ++comp) {
        expect_component_minmax_match(driver_a.native_diagnostics(), driver_b.native_diagnostics(),
                                      comp, "second-advance native diagnostics");
    }
    expect_component_minmax_match(driver_a.wthv_sec_diagnostics(), driver_b.wthv_sec_diagnostics(),
                                  0, "second-advance wthv_sec diagnostics");
    expect_component_minmax_match(driver_a.pblh_diagnostics(), driver_b.pblh_diagnostics(),
                                  0, "second-advance pblh diagnostics");
}

TEST(ShocDriver, HostDiffusionModeExportsDiffusivitiesWithoutInternalTransport)
{
    ScopedParmParseString transport_mode("erf.shoc", "transport_mode", "host_diffusion");

    Box domain(IntVect(0,0,0), IntVect(nx-1, ny-1, nz-1));
    amrex::RealBox real_box(0.0_rt, 0.0_rt, 0.0_rt,
                            500.0_rt, 100.0_rt, 900.0_rt);
    int is_periodic[AMREX_SPACEDIM] = {1, 1, 0};
    Geometry geom(domain, &real_box, amrex::CoordSys::cartesian, is_periodic);
    const FixtureMap fixture = load_driver_fixture();

    amrex::BoxArray ba(domain);
    // The current SHOC driver path maps each MFIter tile to a full-height
    // column workspace, so keep tiling in x only here.
    ba.maxSize(IntVect(3, ny, nz));
    DistributionMapping dm(ba);

    MultiFab cons(ba, dm, NVAR_max, 0);
    amrex::BoxArray xba = amrex::convert(ba, IntVect(1,0,0));
    amrex::BoxArray yba = amrex::convert(ba, IntVect(0,1,0));
    amrex::BoxArray zba = amrex::convert(ba, IntVect(0,0,1));
    amrex::BoxArray xzba = amrex::convert(ba, IntVect(1,0,1));
    amrex::BoxArray yzba = amrex::convert(ba, IntVect(0,1,1));

    MultiFab xvel(xba, dm, 1, 0);
    MultiFab yvel(yba, dm, 1, 0);
    MultiFab zvel(zba, dm, 1, 0);
    MultiFab z_phys_nd(make_nodal_zphys_boxarray(ba), dm, 1, 0);
    MultiFab hfx3(ba, dm, 1, 0);
    MultiFab qfx3(ba, dm, 1, 0);
    MultiFab tau13(xzba, dm, 1, 0);
    MultiFab tau23(yzba, dm, 1, 0);
    MultiFab eddy_diffs(ba, dm, EddyDiff::NumDiffs, 0);

    initialize_state(cons, fixture);
    initialize_velocity(xvel, yvel, zvel, fixture);
    initialize_geometry(z_phys_nd, fixture);
    initialize_surface_fluxes(hfx3, qfx3, tau13, tau23, fixture);
    initialize_eddy_diffs(eddy_diffs, fixture);
    shoc_test::sync();

    MultiFab hfx3_before(ba, dm, 1, 0);
    MultiFab qfx3_before(ba, dm, 1, 0);
    MultiFab tau13_before(xzba, dm, 1, 0);
    MultiFab tau23_before(yzba, dm, 1, 0);
    MultiFab cons_before(ba, dm, NVAR_max, 0);
    MultiFab xvel_before(xba, dm, 1, 0);
    MultiFab yvel_before(yba, dm, 1, 0);
    MultiFab zvel_before(zba, dm, 1, 0);
    MultiFab::Copy(hfx3_before, hfx3, 0, 0, 1, 0);
    MultiFab::Copy(qfx3_before, qfx3, 0, 0, 1, 0);
    MultiFab::Copy(tau13_before, tau13, 0, 0, 1, 0);
    MultiFab::Copy(tau23_before, tau23, 0, 0, 1, 0);
    MultiFab::Copy(cons_before, cons, 0, 0, NVAR_max, 0);
    MultiFab::Copy(xvel_before, xvel, 0, 0, 1, 0);
    MultiFab::Copy(yvel_before, yvel, 0, 0, 1, 0);
    MultiFab::Copy(zvel_before, zvel, 0, 0, 1, 0);

    SolverChoice solver_choice;
    solver_choice.moisture_type = MoistureType::None;

    ShocDriver driver(0, solver_choice);
    EXPECT_FALSE(driver.uses_state_update());
    EXPECT_TRUE(driver.uses_host_diffusion());
    EXPECT_TRUE(driver.uses_momentum_host_diffusion());

    shoc_test::run_and_sync([&] {
        driver.advance(cons, xvel, yvel, zvel,
                       &tau13, &tau23, &hfx3, &qfx3, &eddy_diffs,
                       z_phys_nd, geom, 10.0_rt);
    });

    shoc_test::run_and_sync([&] {
        driver.set_eddy_diffs();
        driver.set_diff_stresses();
    });

    for (int comp = 0; comp < EddyDiff::NumDiffs; ++comp) {
        expect_component_minmax_match(eddy_diffs, driver.native_diagnostics(), comp,
                                      "host-diffusion eddy_diffs writeback");
    }
    expect_component_minmax_match(hfx3, hfx3_before, 0, "host-diffusion hfx3 preservation");
    expect_component_minmax_match(qfx3, qfx3_before, 0, "host-diffusion qfx3 preservation");
    expect_component_minmax_match(tau13, tau13_before, 0, "host-diffusion tau13 preservation");
    expect_component_minmax_match(tau23, tau23_before, 0, "host-diffusion tau23 preservation");

    amrex::Vector<MultiFab> rhs(IntVars::NumTypes);
    rhs[IntVars::cons].define(ba, dm, NVAR_max, 0);
    rhs[IntVars::xmom].define(xba, dm, 1, 0);
    rhs[IntVars::ymom].define(yba, dm, 1, 0);
    rhs[IntVars::zmom].define(zba, dm, 1, 0);
    for (auto& mf : rhs) {
        fill_multifab(mf, 0.0_rt);
    }

    shoc_test::run_and_sync([&] {
        driver.add_fast_tend(rhs);
        for (MFIter mfi(rhs[IntVars::cons], false); mfi.isValid(); ++mfi) {
            driver.add_slow_tend(mfi, mfi.validbox(), rhs[IntVars::cons].array(mfi));
        }
    });

    EXPECT_EQ(component_max_abs_difference(cons, Rho_comp, cons_before, Rho_comp), 0.0_rt);
    EXPECT_EQ(component_max_abs_difference(cons, RhoTheta_comp, cons_before, RhoTheta_comp), 0.0_rt);
    EXPECT_EQ(component_max_abs_difference(cons, RhoQ1_comp, cons_before, RhoQ1_comp), 0.0_rt);
    EXPECT_EQ(component_max_abs_difference(xvel, 0, xvel_before, 0), 0.0_rt);
    EXPECT_EQ(component_max_abs_difference(yvel, 0, yvel_before, 0), 0.0_rt);
    EXPECT_EQ(component_max_abs_difference(zvel, 0, zvel_before, 0), 0.0_rt);

    expect_component_between(rhs[IntVars::cons], RhoTheta_comp, 0.0_rt, 0.0_rt, "host-diffusion theta rhs");
    expect_component_between(rhs[IntVars::cons], RhoQ1_comp, 0.0_rt, 0.0_rt, "host-diffusion qv rhs");
    expect_component_between(rhs[IntVars::cons], RhoKE_comp, 0.0_rt, 0.0_rt, "host-diffusion tke rhs");
    expect_component_between(rhs[IntVars::xmom], 0, 0.0_rt, 0.0_rt, "host-diffusion xmom rhs");
    expect_component_between(rhs[IntVars::ymom], 0, 0.0_rt, 0.0_rt, "host-diffusion ymom rhs");
}

TEST(ShocDriver, NumberAwareMoistureLayoutHelperRecognizesNcAndNi)
{
    EXPECT_TRUE(shoc_layout_requires_number_closure(
        MoistureComponentIndices(
            RhoQ1_comp, RhoQ2_comp, RhoQ3_comp, RhoQ4_comp, RhoQ5_comp, RhoQ6_comp,
            RhoQ7_comp, RhoQ8_comp, RhoQ9_comp, RhoQ10_comp, RhoQ11_comp),
        NVAR_max));

    EXPECT_FALSE(shoc_layout_requires_number_closure(
        MoistureComponentIndices(
            RhoQ1_comp, RhoQ2_comp, RhoQ3_comp, RhoQ4_comp, RhoQ5_comp, RhoQ6_comp),
        NVAR_max));

    EXPECT_FALSE(shoc_layout_requires_number_closure(
        MoistureComponentIndices(RhoQ1_comp, RhoQ2_comp),
        NVAR_max));

    EXPECT_FALSE(shoc_layout_requires_number_closure(
        MoistureComponentIndices(
            RhoQ1_comp, RhoQ2_comp, -1, RhoQ3_comp),
        NVAR_max));

    EXPECT_FALSE(shoc_layout_requires_number_closure(
        MoistureComponentIndices(
            RhoQ1_comp, RhoQ2_comp, RhoQ3_comp, RhoQ4_comp, RhoQ5_comp, RhoQ6_comp,
            -1, -1, -1, -1, -1),
        NVAR_max));
}

TEST(ShocDriver, StateUpdateModeRejectsNumberAwareMoistureLayoutsWithoutNumberClosure)
{
    Box domain(IntVect(0,0,0), IntVect(nx-1, ny-1, nz-1));
    amrex::RealBox real_box(0.0_rt, 0.0_rt, 0.0_rt,
                            500.0_rt, 100.0_rt, 900.0_rt);
    int is_periodic[AMREX_SPACEDIM] = {1, 1, 0};
    Geometry geom(domain, &real_box, amrex::CoordSys::cartesian, is_periodic);
    const FixtureMap fixture = load_driver_fixture();

    amrex::BoxArray ba(domain);
    ba.maxSize(IntVect(3, ny, nz));
    DistributionMapping dm(ba);

    MultiFab cons(ba, dm, NVAR_max, 0);
    amrex::BoxArray xba = amrex::convert(ba, IntVect(1,0,0));
    amrex::BoxArray yba = amrex::convert(ba, IntVect(0,1,0));
    amrex::BoxArray zba = amrex::convert(ba, IntVect(0,0,1));
    amrex::BoxArray xzba = amrex::convert(ba, IntVect(1,0,1));
    amrex::BoxArray yzba = amrex::convert(ba, IntVect(0,1,1));

    MultiFab xvel(xba, dm, 1, 0);
    MultiFab yvel(yba, dm, 1, 0);
    MultiFab zvel(zba, dm, 1, 0);
    MultiFab z_phys_nd(make_nodal_zphys_boxarray(ba), dm, 1, 0);
    MultiFab hfx3(ba, dm, 1, 0);
    MultiFab qfx3(ba, dm, 1, 0);
    MultiFab tau13(xzba, dm, 1, 0);
    MultiFab tau23(yzba, dm, 1, 0);
    MultiFab eddy_diffs(ba, dm, EddyDiff::NumDiffs, 0);

    initialize_state(cons, fixture);
    initialize_velocity(xvel, yvel, zvel, fixture);
    initialize_geometry(z_phys_nd, fixture);
    initialize_surface_fluxes(hfx3, qfx3, tau13, tau23, fixture);
    initialize_eddy_diffs(eddy_diffs, fixture);
    shoc_test::sync();

    SolverChoice solver_choice;
    solver_choice.moisture_type = MoistureType::Morrison;
    solver_choice.moisture_indices = MoistureComponentIndices(
        RhoQ1_comp, RhoQ2_comp, RhoQ3_comp, RhoQ4_comp, RhoQ5_comp, RhoQ6_comp,
        RhoQ7_comp, RhoQ8_comp, RhoQ9_comp, RhoQ10_comp, RhoQ11_comp);

    ShocRuntimeOptions opts;
    opts.transport_mode = ShocTransportMode::StateUpdate;
    std::string error_message;
    EXPECT_FALSE(shoc_driver_state_update_layout_supported(
        opts, solver_choice.moisture_indices, cons.nComp(), error_message));
    EXPECT_NE(error_message.find("number closure"), std::string::npos);
}

TEST(ShocDriver, HostDiffusionWithMoistureIsRejected)
{
    ShocRuntimeOptions opts;
    opts.transport_mode = ShocTransportMode::HostDiffusion;

    std::string error_message;
    EXPECT_FALSE(shoc_driver_host_diffusion_moisture_supported(
        opts, MoistureType::Morrison, error_message));
    EXPECT_NE(error_message.find("host_diffusion with moisture"), std::string::npos);
}

TEST(ShocDriver, StateUpdateModeSupportsIceCapableMoistureLayouts)
{
    Box domain(IntVect(0,0,0), IntVect(nx-1, ny-1, nz-1));
    amrex::RealBox real_box(0.0_rt, 0.0_rt, 0.0_rt,
                            500.0_rt, 100.0_rt, 900.0_rt);
    int is_periodic[AMREX_SPACEDIM] = {1, 1, 0};
    Geometry geom(domain, &real_box, amrex::CoordSys::cartesian, is_periodic);
    const FixtureMap fixture = load_driver_fixture();

    amrex::BoxArray ba(domain);
    // The current SHOC driver path maps each MFIter tile to a full-height
    // column workspace, so keep tiling in x only here.
    ba.maxSize(IntVect(3, ny, nz));
    DistributionMapping dm(ba);

    MultiFab cons(ba, dm, NVAR_max, 0);
    amrex::BoxArray xba = amrex::convert(ba, IntVect(1,0,0));
    amrex::BoxArray yba = amrex::convert(ba, IntVect(0,1,0));
    amrex::BoxArray zba = amrex::convert(ba, IntVect(0,0,1));
    amrex::BoxArray xzba = amrex::convert(ba, IntVect(1,0,1));
    amrex::BoxArray yzba = amrex::convert(ba, IntVect(0,1,1));

    MultiFab xvel(xba, dm, 1, 0);
    MultiFab yvel(yba, dm, 1, 0);
    MultiFab zvel(zba, dm, 1, 0);
    MultiFab z_phys_nd(make_nodal_zphys_boxarray(ba), dm, 1, 0);
    MultiFab hfx3(ba, dm, 1, 0);
    MultiFab qfx3(ba, dm, 1, 0);
    MultiFab tau13(xzba, dm, 1, 0);
    MultiFab tau23(yzba, dm, 1, 0);
    MultiFab eddy_diffs(ba, dm, EddyDiff::NumDiffs, 0);

    initialize_state(cons, fixture);
    initialize_velocity(xvel, yvel, zvel, fixture);
    initialize_geometry(z_phys_nd, fixture);
    initialize_surface_fluxes(hfx3, qfx3, tau13, tau23, fixture);
    initialize_eddy_diffs(eddy_diffs, fixture);
    shoc_test::sync();

    SolverChoice solver_choice;
    solver_choice.moisture_type = MoistureType::WSM6;
    solver_choice.moisture_indices = MoistureComponentIndices(
        RhoQ1_comp, RhoQ2_comp, RhoQ3_comp, RhoQ4_comp, RhoQ5_comp, RhoQ6_comp);

    ShocDriver driver(0, solver_choice);
    EXPECT_TRUE(driver.uses_state_update());
    EXPECT_FALSE(driver.uses_host_diffusion());

    shoc_test::run_and_sync([&] {
        driver.advance(cons, xvel, yvel, zvel,
                       &tau13, &tau23, &hfx3, &qfx3, &eddy_diffs,
                       z_phys_nd, geom, 10.0_rt);
    });

    amrex::Vector<MultiFab> rhs(IntVars::NumTypes);
    rhs[IntVars::cons].define(ba, dm, NVAR_max, 0);
    rhs[IntVars::xmom].define(xba, dm, 1, 0);
    rhs[IntVars::ymom].define(yba, dm, 1, 0);
    rhs[IntVars::zmom].define(zba, dm, 1, 0);
    reset_rhs(rhs);

    shoc_test::run_and_sync([&] {
        driver.add_fast_tend(rhs);
        for (MFIter mfi(rhs[IntVars::cons], false); mfi.isValid(); ++mfi) {
            driver.add_slow_tend(mfi, mfi.validbox(), rhs[IntVars::cons].array(mfi));
        }
    });

    expect_component_between(rhs[IntVars::cons], RhoTheta_comp, 0.0_rt, 0.0_rt, "state-update theta rhs");
    expect_component_between(rhs[IntVars::cons], RhoQ1_comp, 0.0_rt, 0.0_rt, "state-update qv rhs");
    expect_component_between(rhs[IntVars::cons], RhoQ7_comp, 0.0_rt, 0.0_rt,
                             "state-update nc rhs stays zero");
    expect_component_between(rhs[IntVars::cons], RhoQ8_comp, 0.0_rt, 0.0_rt,
                             "state-update ni rhs stays zero");
    expect_component_between(rhs[IntVars::cons], RhoQ9_comp, 0.0_rt, 0.0_rt,
                             "state-update nr rhs stays zero");
    expect_component_between(rhs[IntVars::cons], RhoQ10_comp, 0.0_rt, 0.0_rt,
                             "state-update ns rhs stays zero");
    expect_component_between(rhs[IntVars::cons], RhoQ11_comp, 0.0_rt, 0.0_rt,
                             "state-update ng rhs stays zero");
}

TEST(ShocDriver, MultiStepStateUpdateModeKeepsColumnDiagnosticsStable)
{
    Box domain(IntVect(0,0,0), IntVect(nx-1, ny-1, nz-1));
    amrex::RealBox real_box(0.0_rt, 0.0_rt, 0.0_rt,
                            500.0_rt, 100.0_rt, 900.0_rt);
    int is_periodic[AMREX_SPACEDIM] = {1, 1, 0};
    Geometry geom(domain, &real_box, amrex::CoordSys::cartesian, is_periodic);
    const FixtureMap fixture = load_driver_fixture();
    const Real dt = 10.0_rt;
    const int nsteps = 5;
    const Real max_mix_len = std::sqrt(geom.CellSizeArray()[0] * geom.CellSizeArray()[1]);

    amrex::BoxArray ba(domain);
    // The current SHOC driver path maps each MFIter tile to a full-height
    // column workspace, so keep tiling in x only here.
    ba.maxSize(IntVect(3, ny, nz));
    DistributionMapping dm(ba);

    MultiFab cons(ba, dm, NVAR_max, 0);
    amrex::BoxArray xba = amrex::convert(ba, IntVect(1,0,0));
    amrex::BoxArray yba = amrex::convert(ba, IntVect(0,1,0));
    amrex::BoxArray zba = amrex::convert(ba, IntVect(0,0,1));
    amrex::BoxArray xzba = amrex::convert(ba, IntVect(1,0,1));
    amrex::BoxArray yzba = amrex::convert(ba, IntVect(0,1,1));

    MultiFab xvel(xba, dm, 1, 0);
    MultiFab yvel(yba, dm, 1, 0);
    MultiFab zvel(zba, dm, 1, 0);
    MultiFab z_phys_nd(make_nodal_zphys_boxarray(ba), dm, 1, 0);
    MultiFab hfx3(ba, dm, 1, 0);
    MultiFab qfx3(ba, dm, 1, 0);
    MultiFab tau13(xzba, dm, 1, 0);
    MultiFab tau23(yzba, dm, 1, 0);
    MultiFab eddy_diffs(ba, dm, EddyDiff::NumDiffs, 0);

    initialize_state(cons, fixture);
    initialize_velocity(xvel, yvel, zvel, fixture);
    initialize_geometry(z_phys_nd, fixture);
    initialize_surface_fluxes(hfx3, qfx3, tau13, tau23, fixture);
    initialize_eddy_diffs(eddy_diffs, fixture);
    shoc_test::sync();

    SolverChoice solver_choice;
    solver_choice.moisture_type = MoistureType::None;

    ShocDriver driver(0, solver_choice);
    EXPECT_TRUE(driver.uses_state_update());

    amrex::Vector<MultiFab> rhs(IntVars::NumTypes);
    rhs[IntVars::cons].define(ba, dm, NVAR_max, 0);
    rhs[IntVars::xmom].define(xba, dm, 1, 0);
    rhs[IntVars::ymom].define(yba, dm, 1, 0);
    rhs[IntVars::zmom].define(zba, dm, 1, 0);

    for (int step = 0; step < nsteps; ++step) {
        shoc_test::run_and_sync([&] {
            driver.advance(cons, xvel, yvel, zvel,
                           &tau13, &tau23, &hfx3, &qfx3, &eddy_diffs,
                           z_phys_nd, geom, dt);
        });

        expect_driver_diagnostics_stable(driver, max_mix_len);

        shoc_test::run_and_sync([&] {
            driver.set_eddy_diffs();
            driver.set_diff_stresses();
            reset_rhs(rhs);
            driver.add_fast_tend(rhs);
            for (MFIter mfi(rhs[IntVars::cons], false); mfi.isValid(); ++mfi) {
                driver.add_slow_tend(mfi, mfi.validbox(), rhs[IntVars::cons].array(mfi));
            }
        });

        expect_multifab_finite(cons, Rho_comp, 1, "multistep cons rho");
        expect_multifab_finite(cons, RhoTheta_comp, 1, "multistep cons rho-theta");
        expect_multifab_finite(cons, RhoQ1_comp, 1, "multistep cons rho-qv");
        expect_multifab_finite(cons, RhoKE_comp, 1, "multistep cons rho-ke");
        expect_multifab_finite(xvel, 0, 1, "multistep xvel");
        expect_multifab_finite(yvel, 0, 1, "multistep yvel");
    }
}

TEST(ShocDriver, MultiStepHostDiffusionModeKeepsExportsStable)
{
    ScopedParmParseString transport_mode("erf.shoc", "transport_mode", "host_diffusion");

    Box domain(IntVect(0,0,0), IntVect(nx-1, ny-1, nz-1));
    amrex::RealBox real_box(0.0_rt, 0.0_rt, 0.0_rt,
                            500.0_rt, 100.0_rt, 900.0_rt);
    int is_periodic[AMREX_SPACEDIM] = {1, 1, 0};
    Geometry geom(domain, &real_box, amrex::CoordSys::cartesian, is_periodic);
    const FixtureMap fixture = load_driver_fixture();
    const Real dt = 10.0_rt;
    const int nsteps = 5;
    const Real max_mix_len = std::sqrt(geom.CellSizeArray()[0] * geom.CellSizeArray()[1]);

    amrex::BoxArray ba(domain);
    // The current SHOC driver path maps each MFIter tile to a full-height
    // column workspace, so keep tiling in x only here.
    ba.maxSize(IntVect(3, ny, nz));
    DistributionMapping dm(ba);

    MultiFab cons(ba, dm, NVAR_max, 0);
    amrex::BoxArray xba = amrex::convert(ba, IntVect(1,0,0));
    amrex::BoxArray yba = amrex::convert(ba, IntVect(0,1,0));
    amrex::BoxArray zba = amrex::convert(ba, IntVect(0,0,1));
    amrex::BoxArray xzba = amrex::convert(ba, IntVect(1,0,1));
    amrex::BoxArray yzba = amrex::convert(ba, IntVect(0,1,1));

    MultiFab xvel(xba, dm, 1, 0);
    MultiFab yvel(yba, dm, 1, 0);
    MultiFab zvel(zba, dm, 1, 0);
    MultiFab z_phys_nd(make_nodal_zphys_boxarray(ba), dm, 1, 0);
    MultiFab hfx3(ba, dm, 1, 0);
    MultiFab qfx3(ba, dm, 1, 0);
    MultiFab tau13(xzba, dm, 1, 0);
    MultiFab tau23(yzba, dm, 1, 0);
    MultiFab eddy_diffs(ba, dm, EddyDiff::NumDiffs, 0);

    initialize_state(cons, fixture);
    initialize_velocity(xvel, yvel, zvel, fixture);
    initialize_geometry(z_phys_nd, fixture);
    initialize_surface_fluxes(hfx3, qfx3, tau13, tau23, fixture);
    initialize_eddy_diffs(eddy_diffs, fixture);
    shoc_test::sync();

    SolverChoice solver_choice;
    solver_choice.moisture_type = MoistureType::None;

    ShocDriver driver(0, solver_choice);
    EXPECT_FALSE(driver.uses_state_update());
    EXPECT_TRUE(driver.uses_host_diffusion());
    EXPECT_TRUE(driver.uses_momentum_host_diffusion());

    amrex::Vector<MultiFab> rhs(IntVars::NumTypes);
    rhs[IntVars::cons].define(ba, dm, NVAR_max, 0);
    rhs[IntVars::xmom].define(xba, dm, 1, 0);
    rhs[IntVars::ymom].define(yba, dm, 1, 0);
    rhs[IntVars::zmom].define(zba, dm, 1, 0);

    MultiFab hfx3_before(ba, dm, 1, 0);
    MultiFab qfx3_before(ba, dm, 1, 0);
    MultiFab tau13_before(xzba, dm, 1, 0);
    MultiFab tau23_before(yzba, dm, 1, 0);
    MultiFab cons_before(ba, dm, NVAR_max, 0);
    MultiFab xvel_before(xba, dm, 1, 0);
    MultiFab yvel_before(yba, dm, 1, 0);
    MultiFab zvel_before(zba, dm, 1, 0);
    MultiFab::Copy(hfx3_before, hfx3, 0, 0, 1, 0);
    MultiFab::Copy(qfx3_before, qfx3, 0, 0, 1, 0);
    MultiFab::Copy(tau13_before, tau13, 0, 0, 1, 0);
    MultiFab::Copy(tau23_before, tau23, 0, 0, 1, 0);
    MultiFab::Copy(cons_before, cons, 0, 0, NVAR_max, 0);
    MultiFab::Copy(xvel_before, xvel, 0, 0, 1, 0);
    MultiFab::Copy(yvel_before, yvel, 0, 0, 1, 0);
    MultiFab::Copy(zvel_before, zvel, 0, 0, 1, 0);

    for (int step = 0; step < nsteps; ++step) {
        shoc_test::run_and_sync([&] {
            driver.advance(cons, xvel, yvel, zvel,
                           &tau13, &tau23, &hfx3, &qfx3, &eddy_diffs,
                           z_phys_nd, geom, dt);
            driver.set_eddy_diffs();
            driver.set_diff_stresses();
            reset_rhs(rhs);
            driver.add_fast_tend(rhs);
            for (MFIter mfi(rhs[IntVars::cons], false); mfi.isValid(); ++mfi) {
                driver.add_slow_tend(mfi, mfi.validbox(), rhs[IntVars::cons].array(mfi));
            }
        });

        EXPECT_EQ(component_max_abs_difference(cons, Rho_comp, cons_before, Rho_comp), 0.0_rt);
        EXPECT_EQ(component_max_abs_difference(cons, RhoTheta_comp, cons_before, RhoTheta_comp), 0.0_rt);
        EXPECT_EQ(component_max_abs_difference(cons, RhoQ1_comp, cons_before, RhoQ1_comp), 0.0_rt);
        EXPECT_EQ(component_max_abs_difference(xvel, 0, xvel_before, 0), 0.0_rt);
        EXPECT_EQ(component_max_abs_difference(yvel, 0, yvel_before, 0), 0.0_rt);
        EXPECT_EQ(component_max_abs_difference(zvel, 0, zvel_before, 0), 0.0_rt);
        expect_driver_diagnostics_stable(driver, max_mix_len);
        for (int comp = 0; comp < EddyDiff::NumDiffs; ++comp) {
            expect_component_minmax_match(eddy_diffs, driver.native_diagnostics(), comp,
                                          "multistep host-diffusion eddy_diffs writeback");
        }
        expect_component_minmax_match(hfx3, hfx3_before, 0, "multistep host-diffusion hfx3 preservation");
        expect_component_minmax_match(qfx3, qfx3_before, 0, "multistep host-diffusion qfx3 preservation");
        expect_component_minmax_match(tau13, tau13_before, 0, "multistep host-diffusion tau13 preservation");
        expect_component_minmax_match(tau23, tau23_before, 0, "multistep host-diffusion tau23 preservation");

        expect_component_between(rhs[IntVars::cons], RhoTheta_comp, 0.0_rt, 0.0_rt, "multistep host-diffusion theta rhs");
        expect_component_between(rhs[IntVars::cons], RhoQ1_comp, 0.0_rt, 0.0_rt, "multistep host-diffusion qv rhs");
        expect_component_between(rhs[IntVars::cons], RhoKE_comp, 0.0_rt, 0.0_rt, "multistep host-diffusion tke rhs");
        expect_component_between(rhs[IntVars::xmom], 0, 0.0_rt, 0.0_rt, "multistep host-diffusion xmom rhs");
        expect_component_between(rhs[IntVars::ymom], 0, 0.0_rt, 0.0_rt, "multistep host-diffusion ymom rhs");
    }
}

TEST(ShocDriver, FixturePreprocessStructureAndTkeStayBoundedBeforePdf)
{
    Box domain(IntVect(0,0,0), IntVect(nx-1, ny-1, nz-1));
    amrex::RealBox real_box(0.0_rt, 0.0_rt, 0.0_rt,
                            500.0_rt, 100.0_rt, 900.0_rt);
    int is_periodic[AMREX_SPACEDIM] = {1, 1, 0};
    Geometry geom(domain, &real_box, amrex::CoordSys::cartesian, is_periodic);
    const FixtureMap fixture = load_driver_fixture();
    const Real max_mix_len = std::sqrt(geom.CellSizeArray()[0] * geom.CellSizeArray()[1]);

    amrex::BoxArray ba(domain);
    // The current SHOC driver path maps each MFIter tile to a full-height
    // column workspace, so keep tiling in x only here.
    ba.maxSize(IntVect(3, ny, nz));
    DistributionMapping dm(ba);

    MultiFab cons(ba, dm, NVAR_max, 0);
    amrex::BoxArray xba = amrex::convert(ba, IntVect(1,0,0));
    amrex::BoxArray yba = amrex::convert(ba, IntVect(0,1,0));
    amrex::BoxArray zba = amrex::convert(ba, IntVect(0,0,1));
    amrex::BoxArray xzba = amrex::convert(ba, IntVect(1,0,1));
    amrex::BoxArray yzba = amrex::convert(ba, IntVect(0,1,1));

    MultiFab xvel(xba, dm, 1, 0);
    MultiFab yvel(yba, dm, 1, 0);
    MultiFab zvel(zba, dm, 1, 0);
    MultiFab z_phys_nd(make_nodal_zphys_boxarray(ba), dm, 1, 0);
    MultiFab hfx3(ba, dm, 1, 0);
    MultiFab qfx3(ba, dm, 1, 0);
    MultiFab tau13(xzba, dm, 1, 0);
    MultiFab tau23(yzba, dm, 1, 0);
    MultiFab eddy_diffs(ba, dm, EddyDiff::NumDiffs, 0);

    initialize_state(cons, fixture);
    initialize_velocity(xvel, yvel, zvel, fixture);
    initialize_geometry(z_phys_nd, fixture);
    initialize_surface_fluxes(hfx3, qfx3, tau13, tau23, fixture);
    initialize_eddy_diffs(eddy_diffs, fixture);
    shoc_test::sync();

    expect_multifab_finite(cons, Rho_comp, 1, "cons rho input");
    expect_multifab_finite(cons, RhoTheta_comp, 1, "cons rho-theta input");
    expect_multifab_finite(cons, RhoQ1_comp, 1, "cons rho-qv input");

    ShocRuntimeOptions opts;
    ShocColumnWorkspace workspace;

    for (MFIter mfi(cons, false); mfi.isValid(); ++mfi) {
        const ShocColumnLayout active_layout = make_shoc_layout(mfi.validbox(), geom);
        workspace.ensure_capacity(active_layout, shoc_test::test_arena(), shoc::InitRunOn::Host);
        ShocColumnData& col = workspace.col;

        shoc_test::run_and_sync([&] {
            ShocPreprocess::fill_columns(col, mfi, cons, xvel, yvel, zvel,
                                         &hfx3, &qfx3, &tau13, &tau23, z_phys_nd, geom,
                                         MoistureComponentIndices(RhoQ1_comp, RhoQ2_comp));
            seed_first_step_turbulence(col, cons, eddy_diffs, mfi);
        });

        expect_shoc_column_finite(col.rho, col.layout.nlev, "rho after preprocess");
        expect_shoc_column_finite(col.theta, col.layout.nlev, "theta after preprocess");
        expect_shoc_column_finite(col.qv, col.layout.nlev, "qv after preprocess");
        expect_shoc_column_finite(col.p_mid, col.layout.nlev, "p_mid after preprocess");
        expect_shoc_column_finite(col.tabs, col.layout.nlev, "tabs after preprocess");
        expect_shoc_column_finite(col.exner, col.layout.nlev, "exner after preprocess");
        expect_shoc_column_finite(col.thetal, col.layout.nlev, "thetal after preprocess");

        amrex::Vector<Real> expected_kbfs(col.layout.ncell, 0.0_rt);
        const auto thetal = col.thetal.const_array();
        const auto sens = col.surf_sens_flux.const_array();
        const auto lat = col.surf_lat_flux.const_array();
        for (int ic = 0; ic < col.layout.ncell; ++ic) {
            ASSERT_TRUE(std::isfinite(thetal(ic,0,0)));
            ASSERT_TRUE(std::isfinite(sens(ic,0,0)));
            ASSERT_TRUE(std::isfinite(lat(ic,0,0)));
            expected_kbfs[ic] = sens(ic,0,0) + 0.61_rt * thetal(ic,0,0) * lat(ic,0,0);
            ASSERT_TRUE(std::isfinite(expected_kbfs[ic]));
        }

        shoc_test::run_and_sync([&] {
            ShocStructure::diagnose_surface_layer(col);
            ShocStructure::diagnose_pblh(col);
            ShocStructure::diagnose_length_and_brunt(col, opts,
                                                    geom.CellSizeArray()[0],
                                                    geom.CellSizeArray()[1]);
            ShocTKE::diagnose_tke_and_diffusivities(col, opts, 10.0_rt);
        });

        expect_shoc_column_finite(col.shoc_mix, col.layout.nlev, "shoc_mix");
        expect_shoc_column_finite(col.tk, col.layout.nlev, "tk");
        expect_shoc_column_finite(col.tkh, col.layout.nlev, "tkh");
        expect_shoc_column_finite(col.wthv_sec, col.layout.nlev, "wthv_sec");

        const auto mix = col.shoc_mix.const_array();
        const auto tk = col.tk.const_array();
        const auto tkh = col.tkh.const_array();
        const auto wthv = col.wthv_sec.const_array();
        for (int ic = 0; ic < col.layout.ncell; ++ic) {
            EXPECT_NEAR(wthv(ic,0,0), expected_kbfs[ic], 1.0e-12_rt)
                << "surface buoyancy flux diverged from the preprocess-derived kbfs at column " << ic;
            for (int k = 0; k < col.layout.nlev; ++k) {
                EXPECT_GE(mix(ic,k,0), 20.0_rt);
                EXPECT_LE(mix(ic,k,0), max_mix_len + tol)
                    << "pre-PDF lengthscale exceeded geometry cap at (" << ic << "," << k << ")";
                EXPECT_NEAR(tk(ic,k,0), tkh(ic,k,0), 1.0e-12_rt)
                    << "tk and tkh diverged before implicit/PDF at (" << ic << "," << k << ")";
            }
        }
    }
}
