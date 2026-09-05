#ifndef ERF_GTEST_MOST_AVERAGE_COMMON_H_
#define ERF_GTEST_MOST_AVERAGE_COMMON_H_

#include <AMReX_Geometry.H>
#include <AMReX_Gpu.H>
#include <AMReX_MultiFab.H>
#include <AMReX_ParmParse.H>

#include <ERF_Constants.H>
#include <ERF_MOSTAverage.H>

#include <algorithm>
#include <array>
#include <cmath>
#include <limits>
#include <memory>

namespace erf_most_average_test {

using namespace amrex;

constexpr Real test_u = Real(3.0);
constexpr Real test_v = Real(4.0);
constexpr Real test_w = Real(12.0);
constexpr Real test_theta = Real(300.0);
constexpr Real test_qv = Real(0.01);

class ScopedMOSTParams
{
public:
    explicit ScopedMOSTParams (const char* prefix,
                               const int policy = 0,
                               const int radius = 0)
        : m_pp(prefix)
    {
        m_pp.add("most.average_policy", policy);
        m_pp.add("most.radius", radius);
    }

    void add_zref (const Real zref)
    {
        m_pp.add("most.zref", zref);
    }

    void add_k_index (const int k_index)
    {
        m_pp.addarr("most.k_arr_in", Vector<int>{k_index});
    }

    ~ScopedMOSTParams ()
    {
        for (const auto* name : m_names) {
            m_pp.remove(name);
        }
    }

private:
    ParmParse m_pp;
    const std::array<const char*, 4> m_names{{
        "most.radius",
        "most.average_policy",
        "most.zref",
        "most.k_arr_in"
    }};
};

inline std::array<Real, 9>
expected_average_values ()
{
    return {{
        test_u,
        test_v,
        test_w,
        test_theta,
        test_qv,
        test_theta * (one + epsv * test_qv),
        std::sqrt(test_u*test_u + test_v*test_v),
        std::sqrt(test_u*test_u + test_w*test_w),
        std::sqrt(test_v*test_v + test_w*test_w)
    }};
}

inline Real
tolerance (const Real expected)
{
    return Real(128.0) * std::numeric_limits<Real>::epsilon() *
           std::max(Real(1.0), std::abs(expected));
}

inline void
initialize_uniform_terrain_height (MultiFab& z_phys_nd, const Geometry& geom)
{
    const Real zlo = geom.ProbLo(2);
    const Real dz = geom.CellSize(2);
    for (MFIter mfi(z_phys_nd, false); mfi.isValid(); ++mfi) {
        Box box = z_phys_nd.boxArray()[mfi.index()];
        box.grow(z_phys_nd.nGrowVect());
        auto z_arr = z_phys_nd.array(mfi);
        ParallelFor(box, [=] AMREX_GPU_DEVICE (int i, int j, int k) noexcept
        {
            z_arr(i,j,k) = zlo + static_cast<Real>(k) * dz;
        });
    }
    Gpu::streamSynchronize();
}

} // namespace erf_most_average_test

#endif
