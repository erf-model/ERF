#include <cmath>

#include <AMReX_Geometry.H>
#include <AMReX_MultiFab.H>

#include <gtest/gtest.h>

#include <ERF_Constants.H>
#include <ERF_DataStruct.H>
#include <ERF_EOS.H>
#include <ERF_IndexDefines.H>
#include <ERF_RadiationSimple.H>

namespace {

// One column of nz cells, moist (qv > 8 g/kg) all the way to the top, with cloud
// water only in the top cell, so the inversion is the domain top (itop = nz) and
// the column's whole optical depth sits in that cell.
struct SimpleColumn
{
    static constexpr int nz = 4;
    static constexpr amrex::Real dz  = 100.0;
    static constexpr amrex::Real rho = 1.0;
    static constexpr amrex::Real theta = 300.0;
    static constexpr amrex::Real qv  = 0.01;
    static constexpr amrex::Real qc_top = 1.0e-4;

    amrex::Geometry geom;
    amrex::BoxArray ba;
    amrex::DistributionMapping dm;
    amrex::MultiFab cons;
    amrex::MultiFab qheating;

    SimpleColumn ()
    {
        const amrex::Box domain(amrex::IntVect(0, 0, 0), amrex::IntVect(0, 0, nz - 1));
        const amrex::RealBox rb({0.0, 0.0, 0.0}, {dz, dz, nz * dz});
        geom.define(domain, rb, amrex::CoordSys::cartesian, {1, 1, 0});
        ba = amrex::BoxArray(domain);
        dm = amrex::DistributionMapping(ba);

        // One ghost cell: the top-face flux reads the density at k = nz
        cons.define(ba, dm, RhoQ2_comp + 1, 1);
        cons.setVal(0.0);
        cons.setVal(rho, Rho_comp, 1, 1);
        cons.setVal(rho * theta, RhoTheta_comp, 1, 1);
        cons.setVal(rho * qv, RhoQ1_comp, 1, 1);

        // Cloud water in the top cell only.  (A ParallelFor here would put an extended
        // __device__ lambda in a constructor, which nvcc does not allow.)
        const amrex::Box top(amrex::IntVect(0, 0, nz - 1), amrex::IntVect(0, 0, nz - 1));
        cons.setVal(rho * qc_top, top, RhoQ2_comp, 1);

        qheating.define(ba, dm, 2, 0);
        qheating.setVal(0.0);
    }

    void run ()
    {
        SolverChoice sc;
        sc.moisture_indices.qv = RhoQ1_comp;
        sc.moisture_indices.qc = RhoQ2_comp;

        RadiationSimple rad(0, sc);
        rad.Init(geom, ba, &cons);

        int level = 0;
        int step = 0;
        double time = 0.0;
        const double dt = 1.0;
        amrex::Real declin = 0.0;
        amrex::Real calday = 0.0;
        const amrex::Vector<const amrex::MultiFab*> lsm_in{};
        const amrex::Vector<amrex::MultiFab*> lsm_out{};
        rad.Run(level, step, time, dt, ba, geom, &cons, nullptr, nullptr, lsm_in, lsm_out,
                &qheating, nullptr, nullptr, nullptr, nullptr, false, declin, calday);
        amrex::Gpu::streamSynchronize();
    }

    [[nodiscard]] amrex::Real heating_at (int k) const
    {
        amrex::Real value = 0.0;
        for (amrex::MFIter mfi(qheating); mfi.isValid(); ++mfi) {
            const amrex::Box cell(amrex::IntVect(0, 0, k), amrex::IntVect(0, 0, k));
            if ((cell & mfi.validbox()).isEmpty()) { continue; }
            amrex::FArrayBox host(cell, 1, amrex::The_Pinned_Arena());
            host.copy<amrex::RunOn::Device>(qheating[mfi], cell, 1, cell, 0, 1);
            amrex::Gpu::streamSynchronize();
            value = host(amrex::IntVect(0, 0, k), 0);
        }
        amrex::ParallelDescriptor::ReduceRealSum(value);
        return value;
    }
};

} // namespace

// erf-model/ERF#4484: with the moist layer reaching the top cell the inversion is the
// domain top, and the flux sweep used to run to k = nz and read deltaq there, one cell
// past its box (deltaq has no ghost cells). A bounds-checked build stopped; a release
// build took whatever lay past the box into the top cell's heating. The top cell now
// gets the closed-form answer: below the cloud the flux is 70 e^-tau + 22, at the
// top 70 + 22 e^-tau, with no above-inversion term since zi = z_itop there.
TEST(RadiationSimple, MoistTopCellStaysInsideTheColumn)
{
    SimpleColumn col;
    col.run();

    const amrex::Real tau = 85.0 * SimpleColumn::rho * SimpleColumn::qc_top * SimpleColumn::dz;
    const amrex::Real flux_below = 70.0 * std::exp(-tau) + 22.0;
    const amrex::Real flux_top   = 70.0 + 22.0 * std::exp(-tau);

    const amrex::Real pres  = getPgivenRTh(SimpleColumn::rho * SimpleColumn::theta, SimpleColumn::qv);
    const amrex::Real exner = getExnergivenP(pres, R_d / Cp_d);
    const amrex::Real expected_top =
        -(flux_top - flux_below) / (1015.0 * SimpleColumn::rho * SimpleColumn::dz) / exner;

    for (int k = 0; k < SimpleColumn::nz - 1; ++k) {
        EXPECT_NEAR(col.heating_at(k), 0.0, 1.0e-14) << "k = " << k;
    }
    EXPECT_NEAR(col.heating_at(SimpleColumn::nz - 1), expected_top, 1.0e-12 * std::abs(expected_top));
}
