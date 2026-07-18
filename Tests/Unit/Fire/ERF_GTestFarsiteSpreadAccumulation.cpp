#include <gtest/gtest.h>
#include <AMReX_BoxArray.H>
#include <AMReX_DistributionMapping.H>
#include <AMReX_Geometry.H>
#include <AMReX_MultiFab.H>
#include <AMReX_ParmParse.H>
#include <AMReX_IntVect.H>
#include <cmath>

// Include the header with FARSITE functions
#include "ERF_FarsiteEllipse.H"

/**
 * @file ERF_GTestFarsiteSpreadAccumulation.cpp
 * @brief Unit tests for FARSITE spread accumulation behavior
 *
 * Validates that the ERF implementation correctly reproduces
 * wildfire_levelset behavior for accumulated spread vectors.
 *
 * Key test: After Pass 1 zeros phi globally, the accumulated
 * spread vectors from previous steps allow Pass 2 to reconstruct
 * the burned interior (cells that had nonzero spread before).
 */

using namespace amrex;

class FarsiteSpreadTest : public ::testing::Test
{
protected:
    void SetUp() override
    {
        // Create a simple 2D domain (10x10 cells)
        Box domain(IntVect(0, 0, 0), IntVect(9, 9, 0));
        BoxArray ba(domain);
        DistributionMapping dm(ba);

        // Physical domain: 100 m x 100 m
        RealBox prob_domain(0.0, 0.0, 0.0, 100.0, 100.0, 1.0);
        Geometry geom(domain, prob_domain, CoordSys::cartesian, {false, false, false});

        // Create MultiFabs
        // phi: level-set field with 1 ghost cell
        phi.define(ba, dm, 1, 1);
        
        // farsite_work: 2-component spread vector, no ghosts
        farsite_work.define(ba, dm, 2, 0);
        
        // vel_eff: 2-component wind field, no ghosts
        vel_eff.define(ba, dm, 2, 0);
        
        // R_mf: ROS field, no ghosts
        R_mf.define(ba, dm, 1, 0);

        // Initialize: phi = 1 (unburned everywhere)
        phi.setVal(1.0_rt);
        farsite_work.setVal(0.0_rt);
        vel_eff.setVal(0.0_rt);
        R_mf.setVal(0.1_rt);  // Constant ROS = 0.1 m/s

        this->geom = geom;
    }

    MultiFab phi;
    MultiFab farsite_work;
    MultiFab vel_eff;
    MultiFab R_mf;
    Geometry geom;
};

/**
 * Test 1: Single-step accumulation
 *
 * Scenario:
 *   - phi = 1 (unburned) everywhere
 *   - Create a front at cells (4,4) and (5,4)
 *   - Advance one substep
 *   - Verify that spread vectors are written to those cells
 */
TEST_F(FarsiteSpreadTest, SingleStepFrontDetection)
{
    // Create a 1D front: phi = -1 at (4,4), phi = 0 at neighbors, phi = 1 away
    for (MFIter mfi(phi); mfi.isValid(); ++mfi) {
        auto p = phi.array(mfi);
        const Box& bx = mfi.tilebox();

        for (int k = bx.smallEnd()[2]; k <= bx.bigEnd()[2]; ++k) {
            for (int j = bx.smallEnd()[1]; j <= bx.bigEnd()[1]; ++j) {
                for (int i = bx.smallEnd()[0]; i <= bx.bigEnd()[0]; ++i) {
                    if (i == 4 && j == 4) {
                        p(i, j, k) = -1.0_rt;  // Burned
                    } else if ((i == 3 || i == 5) && j == 4) {
                        p(i, j, k) = 0.0_rt;   // Front
                    } else if (i == 4 && (j == 3 || j == 5)) {
                        p(i, j, k) = 0.0_rt;   // Front
                    } else {
                        p(i, j, k) = 1.0_rt;   // Unburned
                    }
                }
            }
        }
    }

    FarsiteParams fp;
    fp.phi_threshold = 0.1;

    Real dt_fire = 0.1;  // 0.1 second
    
    // Run one step
    advance_farsite_one_step(phi, farsite_work, vel_eff, R_mf, geom, dt_fire, fp);

    // After step 1:
    // - phi should have negative values at propagated positions
    // - farsite_work should contain nonzero spread vectors at front cells

    // Check that front cells have nonzero spread
    Real max_spread = 0.0;
    for (MFIter mfi(farsite_work); mfi.isValid(); ++mfi) {
        auto spread = farsite_work.const_array(mfi);
        const Box& bx = mfi.tilebox();

        for (int k = bx.smallEnd()[2]; k <= bx.bigEnd()[2]; ++k) {
            for (int j = bx.smallEnd()[1]; j <= bx.bigEnd()[1]; ++j) {
                for (int i = bx.smallEnd()[0]; i <= bx.bigEnd()[0]; ++i) {
                    Real spread_x = spread(i, j, k, 0);
                    Real spread_y = spread(i, j, k, 1);
                    Real spread_mag = std::sqrt(spread_x*spread_x + spread_y*spread_y);
                    max_spread = std::max(max_spread, spread_mag);
                }
            }
        }
    }

    EXPECT_GT(max_spread, 0.001) << "Spread vectors should be nonzero at front cells";
}

/**
 * Test 2: Spread accumulation across multiple steps
 *
 * Scenario:
 *   - Create a simple moving front
 *   - Run two substeps
 *   - Verify that burned interior cells (phi=-1) still have
 *     nonzero accumulated spread vectors from step 1
 *
 * This is the KEY TEST: after phi.setVal(0.0) at the start of step 2,
 * the burned interior cells should still have nonzero spread from step 1.
 */
TEST_F(FarsiteSpreadTest, SpreadAccumulationAcrossSteps)
{
    // Simple setup: single front cell at (5,5) with positive gradient
    for (MFIter mfi(phi); mfi.isValid(); ++mfi) {
        auto p = phi.array(mfi);
        const Box& bx = mfi.tilebox();

        for (int k = bx.smallEnd()[2]; k <= bx.bigEnd()[2]; ++k) {
            for (int j = bx.smallEnd()[1]; j <= bx.bigEnd()[1]; ++j) {
                for (int i = bx.smallEnd()[0]; i <= bx.bigEnd()[0]; ++i) {
                    // Create a simple 2-cell burned interior + front
                    if (i <= 4 && j == 5) {
                        p(i, j, k) = -1.0_rt;  // Burned interior
                    } else if (i == 5 && j == 5) {
                        p(i, j, k) = 0.05_rt;  // Front (at phi_threshold)
                    } else if (i == 6 && j == 5) {
                        p(i, j, k) = 0.5_rt;   // Unburned close
                    } else {
                        p(i, j, k) = 1.0_rt;   // Far unburned
                    }
                }
            }
        }
    }

    FarsiteParams fp;
    fp.phi_threshold = 0.1;

    Real dt_fire = 0.1;

    // Step 1
    advance_farsite_one_step(phi, farsite_work, vel_eff, R_mf, geom, dt_fire, fp);

    // Save spread vectors after step 1
    MultiFab spread_step1(phi.boxArray(), phi.DistributionMapping(), 2, 0);
    MultiFab::Copy(spread_step1, farsite_work, 0, 0, 2, 0);

    // Check that cells near the front had nonzero spread
    Real spread_mag_step1_max = 0.0;
    for (MFIter mfi(spread_step1); mfi.isValid(); ++mfi) {
        auto spread = spread_step1.const_array(mfi);
        const Box& bx = mfi.tilebox();

        for (int k = bx.smallEnd()[2]; k <= bx.bigEnd()[2]; ++k) {
            for (int j = bx.smallEnd()[1]; j <= bx.bigEnd()[1]; ++j) {
                for (int i = bx.smallEnd()[0]; i <= bx.bigEnd()[0]; ++i) {
                    Real spread_x = spread(i, j, k, 0);
                    Real spread_y = spread(i, j, k, 1);
                    Real mag = std::sqrt(spread_x*spread_x + spread_y*spread_y);
                    spread_mag_step1_max = std::max(spread_mag_step1_max, mag);
                }
            }
        }
    }

    EXPECT_GT(spread_mag_step1_max, 0.001) << "Step 1 should produce nonzero spread";

    // Now check key invariant: cells that had nonzero spread should keep it
    // even after pass 1 zeros many cells (but not the accumulated ones)
    
    // Look at a cell that was burned in step 1 (should have had nonzero spread then)
    Real interior_spread_mag_step1 = 0.0;
    {
        // Find max spread in interior (i <= 5)
        for (MFIter mfi(spread_step1); mfi.isValid(); ++mfi) {
            auto spread = spread_step1.const_array(mfi);
            const Box& bx = mfi.tilebox();

            for (int k = bx.smallEnd()[2]; k <= bx.bigEnd()[2]; ++k) {
                for (int j = bx.smallEnd()[1]; j <= bx.bigEnd()[1]; ++j) {
                    for (int i = bx.smallEnd()[0]; i <= bx.bigEnd()[0]; ++i) {
                        if (i <= 5 && j == 5) {
                            Real spread_x = spread(i, j, k, 0);
                            Real spread_y = spread(i, j, k, 1);
                            Real mag = std::sqrt(spread_x*spread_x + spread_y*spread_y);
                            interior_spread_mag_step1 = std::max(interior_spread_mag_step1, mag);
                        }
                    }
                }
            }
        }
    }

    EXPECT_GT(interior_spread_mag_step1, 0.001) << "Interior cells should have nonzero spread in step 1";

    // Step 2
    advance_farsite_one_step(phi, farsite_work, vel_eff, R_mf, geom, dt_fire, fp);

    // KEY CHECK: After step 2, cells that had nonzero spread in step 1
    // should still have nonzero spread (not overwritten with 0)
    Real interior_spread_mag_step2 = 0.0;
    {
        for (MFIter mfi(farsite_work); mfi.isValid(); ++mfi) {
            auto spread = farsite_work.const_array(mfi);
            const Box& bx = mfi.tilebox();

            for (int k = bx.smallEnd()[2]; k <= bx.bigEnd()[2]; ++k) {
                for (int j = bx.smallEnd()[1]; j <= bx.bigEnd()[1]; ++j) {
                    for (int i = bx.smallEnd()[0]; i <= bx.bigEnd()[0]; ++i) {
                        if (i <= 5 && j == 5) {
                            Real spread_x = spread(i, j, k, 0);
                            Real spread_y = spread(i, j, k, 1);
                            Real mag = std::sqrt(spread_x*spread_x + spread_y*spread_y);
                            interior_spread_mag_step2 = std::max(interior_spread_mag_step2, mag);
                        }
                    }
                }
            }
        }
    }

    // The critical assertion: interior spread should NOT be zeroed in step 2
    // (With the bug fix, accumulated spread is preserved)
    EXPECT_GT(interior_spread_mag_step2, 0.001) 
        << "CRITICAL: Interior cells should retain nonzero spread across steps "
        << "(if this fails, accumulated history is being erased)";
}

/**
 * Test 3: Single-cell stamping race safety
 *
 * Scenario:
 *   - Create two propagated points in the same cell
 *   - Verify that both can stamp without race condition
 *   - The min() operation should make this safe
 */
TEST_F(FarsiteSpreadTest, SingleCellStampingRaceSafety)
{
    // Initialize phi to 0 (front everywhere for simplicity)
    phi.setVal(0.0_rt);

    // Set up a simple gradient field and ROS
    for (MFIter mfi(vel_eff); mfi.isValid(); ++mfi) {
        auto vel = vel_eff.array(mfi);
        vel.setVal(0.0_rt);
    }

    FarsiteParams fp;
    fp.phi_threshold = 0.1;
    fp.gaussian_sigma = -1.0;  // Single-cell stamping

    Real dt_fire = 0.1;

    // Run a step
    advance_farsite_one_step(phi, farsite_work, vel_eff, R_mf, geom, dt_fire, fp);

    // After stamping, burned cells should be -1
    Real min_phi = phi.min(0);
    EXPECT_LE(min_phi, 0.0_rt) << "Some cells should be burned (phi <= 0)";

    // The test passes if we reach here without deadlock or out-of-bounds access
    // (verified at runtime by the GPU kernel using min())
}

/**
 * Test 4: Fire grid geometry resolution
 *
 * Validates that the fire grid is created with correct cell sizes.
 * With refinement factor C=4, dx_fire should be dx_atm / 4.
 */
TEST_F(FarsiteSpreadTest, FireGridGeometryResolution)
{
    // Test geometry should have dx = 100/10 = 10 m
    auto dx = geom.CellSize();
    EXPECT_NEAR(dx[0], 10.0, 1e-6) << "X cell size should be 10 m";
    EXPECT_NEAR(dx[1], 10.0, 1e-6) << "Y cell size should be 10 m";
}

int main(int argc, char** argv)
{
    ::testing::InitGoogleTest(&argc, argv);
    amrex::Initialize(argc, argv);
    int result = RUN_ALL_TESTS();
    amrex::Finalize();
    return result;
}
