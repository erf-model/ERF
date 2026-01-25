/**
 * \file ERF_InitCustomTerrain.cpp
 */

#include "AMReX_ParmParse.H"
#include "ERF_Constants.H"
#include "ERF_TerrainMetrics.H"

using namespace amrex;

void
init_my_custom_terrain ( const Geometry& geom,
                         FArrayBox& terrain_fab,
                         const Real& time )
{
    //
    // We put this here as a convenience for testing the map factor implementation
    // Note that these factors must match those in Source/ERF_MakeNewArrays.cpp
    //
    ParmParse pp("erf");
    bool test_mapfactor = false;
    pp.query("test_mapfactor",test_mapfactor);

    Real mf_m;
    if (test_mapfactor) {
        mf_m = 0.5;
    } else {
        mf_m = 1.;
    }

    // Domain cell size and real bounds
    auto dx = geom.CellSizeArray();
    auto ProbLoArr = geom.ProbLoArray();
    auto ProbHiArr = geom.ProbHiArray();

    const amrex::Box& domain = geom.Domain();
    int domlo_x = domain.smallEnd(0); int domhi_x = domain.bigEnd(0) + 1;
    int domlo_y = domain.smallEnd(1); int domhi_y = domain.bigEnd(1) + 1;
    int domlo_z = domain.smallEnd(2);

    // User function parameters
    Real a    = 0.5;
    Real num  = 8. * a * a * a;
    Real xcen = 0.5 * (ProbLoArr[0] + ProbHiArr[0]) / mf_m;
    Real ycen = 0.5 * (ProbLoArr[1] + ProbHiArr[1]) / mf_m;

    // Populate bottom plane
    int k0 = domlo_z;

    std::string custom_terrain_type = "None";
    ParmParse pp_prob("prob"); pp_prob.query("custom_terrain_type", custom_terrain_type);
    amrex::Print() << "IN CUSTOM TERRAIN WITH TYPE = " << custom_terrain_type << std::endl;

    amrex::Box zbx = terrain_fab.box();
    if (zbx.smallEnd(2) <= k0)
    {
        amrex::Array4<Real> const& z_arr = terrain_fab.array();

        if (custom_terrain_type == "WoA") {

            // Default to x-direction
            int dir = 0;  pp_prob.query("dir", dir);

            Real  L        = 100.0; pp_prob.query("L"        , L);
            Real  z_offset =   0.0; pp_prob.query("z_offset" , z_offset);

            // If hm is nonzero, then use alternate hill definition
            Real  hm       =   0.0; pp_prob.query("hmax" , hm);

            // This is a 2D hill with variation in only the x-direction
            if (dir == 0) {
                ParallelFor(zbx, [=] AMREX_GPU_DEVICE (int i, int j, int)
                {
                    // Clip indices for ghost-cells
                    int ii = amrex::min(amrex::max(i,domlo_x),domhi_x);

                    // Location of nodes
                    Real x = (ProbLoArr[0] + ii * dx[0] - xcen) * mf_m;

                    // WoA Hill in x-direction
                    if (hm==0) {
                        z_arr(i,j,k0) = num / (x*x + 4 * a * a);
                    } else {
                        Real x_L = x / L;
                        z_arr(i,j,k0) = hm / (1 + x_L*x_L) + z_offset;
                    }
                });
            } else if (dir == 1) {
                ParallelFor(zbx, [=] AMREX_GPU_DEVICE (int i, int j, int)
                {
                    // Clip indices for ghost-cells
                    int jj = amrex::min(amrex::max(j,domlo_y),domhi_y);

                    // Location of nodes
                    Real y = (ProbLoArr[1] + jj * dx[1] - ycen) * mf_m;

                    // WoA Hill in y-direction
                    if (hm==0) {
                        z_arr(i,j,k0) = num / (y*y + 4.0 * a * a);
                    } else {
                        Real y_L = y / L;
                        z_arr(i,j,k0) = hm / (1.0 + y_L*y_L) + z_offset;
                    }
                });
            } else if (dir == 2) {
                ParallelFor(zbx, [=] AMREX_GPU_DEVICE (int i, int j, int)
                {
                    // Clip indices for ghost-cells
                    int ii = amrex::min(amrex::max(i,domlo_x),domhi_x);
                    int jj = amrex::min(amrex::max(j,domlo_y),domhi_y);

                    // Location of nodes
                    Real x = (ProbLoArr[0] + ii * dx[0] - xcen) * mf_m;
                    Real y = (ProbLoArr[1] + jj * dx[1] - ycen) * mf_m;
                    Real r = std::sqrt(x*x + y*y);

                    // WoA Hill in radial direction
                    if (hm==0) {
                        z_arr(i,j,k0) = num / (r*r + 4.0 * a * a);
                    } else {
                        Real r_L = r / L;
                        z_arr(i,j,k0) = hm / (1.0 + r_L*r_L) + z_offset;
                    }
                });
            } else {
                amrex::Abort("Unknown dir in ERF_Prob.cpp");
            }

        } else if (custom_terrain_type == "ScharMountain") {

            Real asq    = 5000.0 * 5000.0;
            Real Hm     =  250.0;
            Real lambda = 4000.0;

            ParallelFor(zbx, [=] AMREX_GPU_DEVICE (int i, int j, int)
            {
                // Clip indices for ghost-cells
                int ii = amrex::min(amrex::max(i,domlo_x),domhi_x);

                // Location of nodes
                Real x = (ProbLoArr[0] + ii * dx[0] - xcen);

                Real cosx = cos(PI * x / lambda);

                z_arr(i,j,k0) = Hm * std::exp(-x*x/asq) * cosx * cosx;
            });

        } else if (custom_terrain_type == "HalfCylinder") {

            Real asq = 0.5 * 0.5;

            ParallelFor(zbx, [=] AMREX_GPU_DEVICE (int i, int j, int)
            {
                // Clip indices for ghost-cells
                int ii = amrex::min(amrex::max(i,domlo_x),domhi_x);

                // Location of nodes
                Real x = (ProbLoArr[0] + ii * dx[0] - xcen);

                Real rsq = x*x;

                if (rsq < asq) {
                    z_arr(i,j,k0) = pow(asq - rsq, 0.5);
                } else {
                    z_arr(i,j,k0) = 0.0;
                }
            });

        } else if (custom_terrain_type == "Hemisphere") {

            Real asq = 0.5 * 0.5;

            ParallelFor(zbx, [=] AMREX_GPU_DEVICE (int i, int j, int)
            {
                // Clip indices for ghost-cells
                int ii = amrex::min(amrex::max(i,domlo_x),domhi_x);
                int jj = amrex::min(amrex::max(j,domlo_y),domhi_y);

                // Location of nodes
                Real x = (ProbLoArr[0] + ii * dx[0] - xcen);
                Real y = (ProbLoArr[1] + jj * dx[1] - ycen);

                Real rsq = x*x + y*y;

                if (rsq < asq) {
                    z_arr(i,j,k0) = std::pow(asq-rsq, 0.5);
                } else {
                    z_arr(i,j,k0) = 0.0;
                }
            });

        } else if (custom_terrain_type == "MovingSineWave") {

            Real Ampl        = 0.0;  pp_prob.query("Ampl", Ampl);
            Real wavelength  = 100.; pp_prob.query("wavelength", wavelength);

            Real kp          = 2.0 * PI / wavelength;
            Real g           = CONST_GRAV;
            Real omega       = std::sqrt(g * kp);

            ParallelFor(zbx, [=] AMREX_GPU_DEVICE (int i, int j, int)
            {
                // Clip indices for ghost-cells
                int ii = amrex::min(amrex::max(i,domlo_x),domhi_x);

                // Location of nodes
                Real x = ii  * dx[0];

                // Wave height
                Real height = Ampl * std::sin(kp * x - omega * time);

                // Populate terrain height
                z_arr(i,j,0) = height;
            });

        } else if (custom_terrain_type == "WindFarmTest") {

            ParallelFor(zbx, [=] AMREX_GPU_DEVICE (int i, int j, int)
            {
                // Clip indices for ghost-cells
                int ii = amrex::min(amrex::max(i,domlo_x),domhi_x);
                int jj = amrex::min(amrex::max(j,domlo_y),domhi_y);

                // Location of nodes
                Real x = (ProbLoArr[0] + ii * dx[0] - xcen);
                Real y = (ProbLoArr[1] + jj * dx[1] - ycen);

                Real x_L = x/100.0;
                Real y_L = y/100.0;

                z_arr(i,j,k0) = 100.0 / (1.0 + x_L*x_L + y_L*y_L);
            });

        } else if (custom_terrain_type == "RaisedFlat") {
            Real  z_offset = 0.0; pp_prob.query("z_offset" , z_offset);
            ParallelFor(zbx, [=] AMREX_GPU_DEVICE (int i, int j, int)
            {
                z_arr(i,j,k0) = z_offset;
            });

        } else if (custom_terrain_type == "None") {
            ParallelFor(zbx, [=] AMREX_GPU_DEVICE (int i, int j, int)
            {
                z_arr(i,j,k0) = 0.0;
            });
        } else {
            Abort("Don't know this custom_terrain_type");
        }
    }
}
