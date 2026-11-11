#include <algorithm>
#include <cmath>
#include <limits>
#include <memory>
#include <string>
#include <vector>

#include <AMReX_Geometry.H>
#include <AMReX_Gpu.H>
#include <AMReX_MultiFab.H>
#include <AMReX_ParmParse.H>

#include <gtest/gtest.h>

#include "../../../Exec/ERF_Prob.H"

namespace rico_moisture_source_test {

// The RICO large-scale moisture tendency (van Zanten et al. 2011): -1 g/kg/day at the
// surface, rising linearly by 1.3456 g/kg/day over 2980 m and constant above
constexpr amrex::Real kSlope = amrex::Real(5.226199353716133e-12);   // 1/s per m
constexpr amrex::Real kRate  = amrex::Real(-1.157e-8);              // 1/s at z = 0
constexpr int kNz = 40;                                              // dz = 100 m
constexpr amrex::Real kDz = amrex::Real(100.0);

// The source on a 1 x 1 x 40 column with cell centres at 50, 150, ..., 3950 m
std::vector<amrex::Real> source_column ()
{
    const amrex::Box domain(amrex::IntVect(0,0,0), amrex::IntVect(0,0,kNz-1));
    const amrex::RealBox rb({AMREX_D_DECL(amrex::Real(0.0), amrex::Real(0.0), amrex::Real(0.0))},
                            {AMREX_D_DECL(amrex::Real(100.0), amrex::Real(100.0), kNz * kDz)});
    const amrex::Geometry geom(domain, rb, 0, {AMREX_D_DECL(0,0,0)});
    const amrex::BoxArray ba(domain);
    const amrex::DistributionMapping dm(ba);
    amrex::MFInfo info;
#ifdef AMREX_USE_GPU
    info.SetArena(amrex::The_Pinned_Arena());
#endif
    amrex::MultiFab qsrc(ba, dm, 1, 0, info);
    qsrc.setVal(amrex::Real(0.0));
    auto z_cc = std::make_unique<amrex::MultiFab>(ba, dm, 1, 0, info);
    amrex::Vector<amrex::Real> zlevels_stag(kNz + 1);
    for (int k = 0; k <= kNz; ++k) { zlevels_stag[k] = k * kDz; }
    for (amrex::MFIter mfi(*z_cc); mfi.isValid(); ++mfi) {
        const auto& a = z_cc->array(mfi);
        for (int k = 0; k < kNz; ++k) { a(0,0,k) = (k + amrex::Real(0.5)) * kDz; }
    }

    Problem prob(geom.ProbLo(), geom.ProbHi());
    prob.update_rhoqt_sources(0.0, &qsrc, geom, z_cc, zlevels_stag);
    amrex::Gpu::streamSynchronize();

    std::vector<amrex::Real> out(kNz, std::numeric_limits<amrex::Real>::quiet_NaN());
    for (amrex::MFIter mfi(qsrc); mfi.isValid(); ++mfi) {
        const auto& a = qsrc.const_array(mfi);
        for (int k = 0; k < kNz; ++k) { out[k] = a(0,0,k); }
    }
    return out;
}

// Sets erf.prob_name = RICO and the moisture inputs for one test, and removes them after
struct RICOInputs
{
    explicit RICOInputs (bool set_cutoff, amrex::Real cutoff = 0)
    {
        amrex::ParmParse pp_erf("erf");
        pp_erf.add("prob_name", std::string("RICO"));
        amrex::ParmParse pp_prob("prob");
        pp_prob.add("advection_moisture_rate", kRate);
        if (set_cutoff) { pp_prob.add("moisture_source_cutoff", cutoff); }
    }
    ~RICOInputs ()
    {
        amrex::ParmParse pp_erf("erf");
        pp_erf.remove("prob_name");
        amrex::ParmParse pp_prob("prob");
        pp_prob.remove("advection_moisture_rate");
        pp_prob.remove("moisture_source_cutoff");
    }
};

amrex::Real tol () { return amrex::Real(64.0) * std::numeric_limits<amrex::Real>::epsilon() * std::abs(kRate); }

} // namespace rico_moisture_source_test

using namespace rico_moisture_source_test;

// Motivation: above prob.moisture_source_cutoff the source held the linear profile's value
// at 2980 m whatever the cutoff, so a cutoff of 500 m made it jump by 1.12 g/kg/day there
// (erf-model/ERF#4626). It now holds the value at the cutoff.
TEST(RICOMoistureSource, ContinuousAtAConfiguredCutoff)
{
    RICOInputs inputs(true, amrex::Real(500.0));
    const auto q = source_column();
    for (int k = 0; k < kNz; ++k) {
        const amrex::Real z = (k + amrex::Real(0.5)) * kDz;
        const amrex::Real expected = kSlope * std::min(z, amrex::Real(500.0)) + kRate;
        EXPECT_NEAR(q[k], expected, tol()) << "z = " << z;
    }
    // the cells either side of the cutoff differ by at most the slope over one cell
    EXPECT_LE(std::abs(q[5] - q[4]), kSlope * kDz + tol());
}

// Motivation: RICO's profile changes at 2980 m, so that is the default cutoff (it was
// 500 m, where the source jumped).
TEST(RICOMoistureSource, DefaultCutoffIsTheRICOHeight)
{
    RICOInputs inputs(false);
    const auto q = source_column();
    for (int k = 0; k < kNz; ++k) {
        const amrex::Real z = (k + amrex::Real(0.5)) * kDz;
        const amrex::Real expected = kSlope * std::min(z, amrex::Real(2980.0)) + kRate;
        EXPECT_NEAR(q[k], expected, tol()) << "z = " << z;
    }
}
