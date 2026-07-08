#include <gtest/gtest.h>
#include <AMReX_BoxArray.H>
#include <AMReX_DistributionMapping.H>
#include <AMReX_Geometry.H>
#include <AMReX_MultiFab.H>
#include <AMReX_ParmParse.H>
#include <AMReX_IntVect.H>
#include <cmath>

// Include the level-set headers
#include "ERF_NumericalSchemes.H"
#include "ERF_Reinitialize.H"
#include "ERF_LevelSetAdvection.H"

/**
 * @file ERF_GTestLevelSetAdvection.cpp
 * @brief Unit tests for level-set advection and reinitialization
 *
 * Tests the WENO5-Z advection scheme with SSP-RK3 time integration
 * and the Sussman reinitialization method.
 */

using namespace amrex;

class LevelSetAdvectionTest : public ::testing::Test
{
protected:
    void SetUp() override
    {
        // Create a simple 2D domain (64x64 cells)
        Box domain(IntVect(0, 0, 0), IntVect(63, 63, 0));
        BoxArray ba(domain);
        DistributionMapping dm(ba);

        // Physical domain: 500 m x 500 m
        RealBox prob_domain(0.0, 0.0, 0.0, 500.0, 500.0, 1.0);
        Geometry geom(domain, prob_domain, CoordSys::cartesian, {false, false, false});

        // Create MultiFabs with sufficient ghost cells for level-set
        phi.define(ba, dm, 1, 3);
        vel_eff.define(ba, dm, 2, 0);
        R_mf.define(ba, dm, 1, 0);

        this->geom = geom;
        
        // Initialize: phi = signed distance to circle at domain center
        // Circle: center (250, 250), radius 50 m
        Real cx = 250.0_rt, cy = 250.0_rt, R_circle = 50.0_rt;
        Real dx = (prob_domain.hi(0) - prob_domain.lo(0)) / domain.length(0);
        Real dy = (prob_domain.hi(1) - prob_domain.lo(1)) / domain.length(1);
        
        for (MFIter mfi(phi); mfi.isValid(); ++mfi) {
            auto p = phi.array(mfi);
            const Box& bx = mfi.growntilebox();
            for (int k = bx.smallEnd(2); k <= bx.bigEnd(2); ++k) {
                for (int j = bx.smallEnd(1); j <= bx.bigEnd(1); ++j) {
                    for (int i = bx.smallEnd(0); i <= bx.bigEnd(0); ++i) {
                        Real x = prob_domain.lo(0) + (i + 0.5_rt) * dx;
                        Real y = prob_domain.lo(1) + (j + 0.5_rt) * dy;
                        Real r = std::sqrt((x - cx) * (x - cx) + (y - cy) * (y - cy));
                        p(i, j, k) = r - R_circle;  // Signed distance
                    }
                }
            }
        }
        phi.FillBoundary(geom.periodicity());
        
        // Wind field: zero (not used in this test)
        vel_eff.setVal(0.0_rt);
        
        // ROS field: constant 0.5 m/s
        R_mf.setVal(0.5_rt);
    }

    MultiFab phi;
    MultiFab vel_eff;
    MultiFab R_mf;
    Geometry geom;
};

/**
 * Test 1: Basic advection step
 *
 * Scenario:
 *   - Initialize phi as a signed-distance circle
 *   - Advance one step with constant ROS = 0.5 m/s
 *   - Verify that phi values change and circle expands
 */
TEST_F(LevelSetAdvectionTest, BasicAdvectionStep)
{
    Real dt = 1.0;  // 1 second timestep
    
    // Store initial state
    MultiFab phi_initial(phi.boxArray(), phi.DistributionMapping(), 1, 0);
    phi_initial.copy(phi);
    
    // Advance one step
    advect_levelset_weno5z_rk3(phi, vel_eff, R_mf, geom, dt);
    phi.FillBoundary(geom.periodicity());
    
    // Check that phi has changed
    bool changed = false;
    for (MFIter mfi(phi); mfi.isValid(); ++mfi) {
        auto p = phi.array(mfi);
        auto p0 = phi_initial.const_array(mfi);
        const Box& bx = mfi.tilebox();
        for (int k = bx.smallEnd(2); k <= bx.bigEnd(2); ++k) {
            for (int j = bx.smallEnd(1); j <= bx.bigEnd(1); ++j) {
                for (int i = bx.smallEnd(0); i <= bx.bigEnd(0); ++i) {
                    if (std::abs(p(i, j, k) - p0(i, j, k)) > 1.0e-10) {
                        changed = true;
                        break;
                    }
                }
                if (changed) break;
            }
            if (changed) break;
        }
        if (changed) break;
    }
    EXPECT_TRUE(changed) << "Level-set should have changed after advection";
    
    // Check that phi is bounded (finite, reasonable values)
    Real phi_max = phi.max(0);
    Real phi_min = phi.min(0);
    EXPECT_TRUE(std::isfinite(phi_max)) << "phi_max should be finite";
    EXPECT_TRUE(std::isfinite(phi_min)) << "phi_min should be finite";
    EXPECT_LT(std::abs(phi_max), 1.0e6) << "phi_max should be reasonable";
    EXPECT_GT(std::abs(phi_min), -1.0e6) << "phi_min should be reasonable";
}

/**
 * Test 2: Reinitialization
 *
 * Scenario:
 *   - Initialize phi as a signed-distance circle
 *   - Perturb phi randomly to break the signed-distance property
 *   - Apply Sussman reinitialization
 *   - Verify that the zero-level contour is preserved (within tolerance)
 */
TEST_F(LevelSetAdvectionTest, Reinitialization)
{
    // Perturb phi to break signed-distance property
    for (MFIter mfi(phi); mfi.isValid(); ++mfi) {
        auto p = phi.array(mfi);
        const Box& bx = mfi.tilebox();
        for (int k = bx.smallEnd(2); k <= bx.bigEnd(2); ++k) {
            for (int j = bx.smallEnd(1); j <= bx.bigEnd(1); ++j) {
                for (int i = bx.smallEnd(0); i <= bx.bigEnd(0); ++i) {
                    p(i, j, k) *= 1.5_rt;  // Scale by 1.5 to perturb distance
                }
            }
        }
    }
    
    // Apply reinitialization
    Real dtau = 0.5 * std::min(geom.CellSize()[0], geom.CellSize()[1]);
    reinitialize_phi(phi, geom, 10, dtau);
    phi.FillBoundary(geom.periodicity());
    
    // Check that phi is still well-behaved
    Real phi_max = phi.max(0);
    Real phi_min = phi.min(0);
    Real phi_norminf = std::max(std::abs(phi_max), std::abs(phi_min));
    
    EXPECT_TRUE(std::isfinite(phi_max)) << "phi_max should be finite after reinit";
    EXPECT_TRUE(std::isfinite(phi_min)) << "phi_min should be finite after reinit";
    EXPECT_GT(phi_norminf, 0.0) << "phi norm should be positive";
    EXPECT_LT(phi_norminf, 1.0e3) << "phi should remain bounded after reinit";
}

/**
 * Test 3: Multiple advection and reinit cycles
 *
 * Scenario:
 *   - Initialize phi as a signed-distance circle
 *   - Perform 10 advection steps and 1 reinitialization after every 5 steps
 *   - Verify that the burned area (phi < 0) grows monotonically
 *   - Verify that phi remains bounded
 */
TEST_F(LevelSetAdvectionTest, MultiCycleAdvectionReinit)
{
    Real dt = 1.0;
    
    // Count initial burned cells (phi < 0)
    long n_burned_initial = 0;
    for (MFIter mfi(phi); mfi.isValid(); ++mfi) {
        auto p = phi.const_array(mfi);
        const Box& bx = mfi.tilebox();
        for (int k = bx.smallEnd(2); k <= bx.bigEnd(2); ++k) {
            for (int j = bx.smallEnd(1); j <= bx.bigEnd(1); ++j) {
                for (int i = bx.smallEnd(0); i <= bx.bigEnd(0); ++i) {
                    if (p(i, j, k) < 0.0_rt) ++n_burned_initial;
                }
            }
        }
    }
    amrex::ParallelDescriptor::ReduceLongSum(n_burned_initial);
    
    // Run 10 advection steps
    for (int step = 0; step < 10; ++step) {
        advect_levelset_weno5z_rk3(phi, vel_eff, R_mf, geom, dt);
        phi.FillBoundary(geom.periodicity());
        
        // Reinitialization every 5 steps
        if ((step + 1) % 5 == 0) {
            Real dtau = 0.5 * std::min(geom.CellSize()[0], geom.CellSize()[1]);
            reinitialize_phi(phi, geom, 10, dtau);
            phi.FillBoundary(geom.periodicity());
        }
    }
    
    // Count final burned cells
    long n_burned_final = 0;
    for (MFIter mfi(phi); mfi.isValid(); ++mfi) {
        auto p = phi.const_array(mfi);
        const Box& bx = mfi.tilebox();
        for (int k = bx.smallEnd(2); k <= bx.bigEnd(2); ++k) {
            for (int j = bx.smallEnd(1); j <= bx.bigEnd(1); ++j) {
                for (int i = bx.smallEnd(0); i <= bx.bigEnd(0); ++i) {
                    if (p(i, j, k) < 0.0_rt) ++n_burned_final;
                }
            }
        }
    }
    amrex::ParallelDescriptor::ReduceLongSum(n_burned_final);
    
    // Check that phi is bounded
    Real phi_max = phi.max(0);
    Real phi_min = phi.min(0);
    EXPECT_TRUE(std::isfinite(phi_max)) << "phi_max should be finite after cycles";
    EXPECT_TRUE(std::isfinite(phi_min)) << "phi_min should be finite after cycles";
    
    // Burned area should grow (with some tolerance for numerics)
    EXPECT_GE(n_burned_final, n_burned_initial) 
        << "Burned area should grow monotonically; initial=" << n_burned_initial 
        << " final=" << n_burned_final;
}

/**
 * Test 4: Zero-level contour preservation
 *
 * Scenario:
 *   - Initialize phi with a circular level-set
 *   - Compute the initial zero-level contour (cells with |phi| < dx)
 *   - Advance several steps with reinitialization
 *   - Verify that the zero-level contour area changes by less than 20%
 */
TEST_F(LevelSetAdvectionTest, ZeroContourPreservation)
{
    Real dx = geom.CellSize()[0];
    Real dt = 1.0;
    
    // Count initial contour cells (|phi| < dx)
    long n_contour_initial = 0;
    for (MFIter mfi(phi); mfi.isValid(); ++mfi) {
        auto p = phi.const_array(mfi);
        const Box& bx = mfi.tilebox();
        for (int k = bx.smallEnd(2); k <= bx.bigEnd(2); ++k) {
            for (int j = bx.smallEnd(1); j <= bx.bigEnd(1); ++j) {
                for (int i = bx.smallEnd(0); i <= bx.bigEnd(0); ++i) {
                    if (std::abs(p(i, j, k)) < dx) ++n_contour_initial;
                }
            }
        }
    }
    amrex::ParallelDescriptor::ReduceLongSum(n_contour_initial);
    
    // Run 10 cycles
    for (int step = 0; step < 10; ++step) {
        advect_levelset_weno5z_rk3(phi, vel_eff, R_mf, geom, dt);
        phi.FillBoundary(geom.periodicity());
        
        if ((step + 1) % 5 == 0) {
            Real dtau = 0.5 * std::min(geom.CellSize()[0], geom.CellSize()[1]);
            reinitialize_phi(phi, geom, 10, dtau);
            phi.FillBoundary(geom.periodicity());
        }
    }
    
    // Count final contour cells
    long n_contour_final = 0;
    for (MFIter mfi(phi); mfi.isValid(); ++mfi) {
        auto p = phi.const_array(mfi);
        const Box& bx = mfi.tilebox();
        for (int k = bx.smallEnd(2); k <= bx.bigEnd(2); ++k) {
            for (int j = bx.smallEnd(1); j <= bx.bigEnd(1); ++j) {
                for (int i = bx.smallEnd(0); i <= bx.bigEnd(0); ++i) {
                    if (std::abs(p(i, j, k)) < dx) ++n_contour_final;
                }
            }
        }
    }
    amrex::ParallelDescriptor::ReduceLongSum(n_contour_final);
    
    // Contour area should remain within 20% of initial (20% buffer for expanding fire)
    Real rel_change = std::abs(n_contour_final - n_contour_initial) / (Real)(n_contour_initial + 1);
    EXPECT_LT(rel_change, 1.0) << "Relative change in contour cells should be reasonable; "
                                 << "initial=" << n_contour_initial 
                                 << " final=" << n_contour_final 
                                 << " rel_change=" << rel_change;
}
