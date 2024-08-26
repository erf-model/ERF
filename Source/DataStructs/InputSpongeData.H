#ifndef _INPUT_SPONGE_DATA_H_
#define _INPUT_SPONGE_DATA_H_

#include <string>
#include <iostream>

#include <AMReX_ParmParse.H>
#include <AMReX_Print.H>
#include <AMReX_Gpu.H>
#include <AMReX_Geometry.H>

#include <ERF_Constants.H>
#include <Interpolation_1D.H>

/**
 * Data structure storing input sponge data. Also
 * handles reading the input file for sponge data
 */
struct InputSpongeData {
public:
    InputSpongeData ()
    {
        // Text input_sounding file
        amrex::ParmParse pp("erf");
        pp.query("input_sponge_file", input_sponge_file);
    }

    void read_from_file (const amrex::Geometry &geom,
                         const amrex::Vector<amrex::Real>& zlevels_stag)
    {
        const int klo = 0;
        const int khi = geom.Domain().bigEnd()[AMREX_SPACEDIM-1];
        const int Nz = geom.Domain().size()[AMREX_SPACEDIM-1];
        const amrex::Real dz = geom.CellSize()[AMREX_SPACEDIM-1];

        const bool use_terrain = (zlevels_stag.size() > 0);
        const amrex::Real zbot = (use_terrain) ? zlevels_stag[klo]   : geom.ProbLo(AMREX_SPACEDIM-1);
        const amrex::Real ztop = (use_terrain) ? zlevels_stag[khi+1] : geom.ProbHi(AMREX_SPACEDIM-1);

        z_inp_sponge.resize(Nz+2);
        U_inp_sponge.resize(Nz+2);
        V_inp_sponge.resize(Nz+2);

        // Read the input_sponge file
        amrex::Print() << "input_sponge file location : " << input_sponge_file << std::endl;
        std::ifstream input_sponge_reader(input_sponge_file);
        if(!input_sponge_reader.is_open()) {
            amrex::Error("Error opening the input_sponge file.\n");
        }
        else {
            // Read the contents of the input_sponge file
            amrex::Print() << "Successfully opened the input_sponge file. Now reading... " << std::endl;
            std::string line;

            // First, read the input data into temp vectors; then, interpolate vectors to the
            // domain lo/hi and cell centers (from level 0)
            amrex::Vector<amrex::Real> z_inp_sponge_tmp, U_inp_sponge_tmp, V_inp_sponge_tmp;

            // Add surface
            z_inp_sponge_tmp.push_back(zbot); // height above sea level [m]
            U_inp_sponge_tmp.push_back(0);
            V_inp_sponge_tmp.push_back(0);

            // Read the vertical profile at each given height
            amrex::Real z, U, V;
            while(std::getline(input_sponge_reader, line)) {
                std::istringstream iss_z(line);
                iss_z >> z >> U >> V;
                if (z == zbot) {
                    U_inp_sponge_tmp[0] = U;
                    V_inp_sponge_tmp[0] = V;
                } else {
                    AMREX_ALWAYS_ASSERT(z > z_inp_sponge_tmp[z_inp_sponge_tmp.size()-1]); // sounding is increasing in height
                    z_inp_sponge_tmp.push_back(z);
                    U_inp_sponge_tmp.push_back(U);
                    V_inp_sponge_tmp.push_back(V);
                    if (z >= ztop) break;
                }
            }

            // At this point, we have an input_sponge from zbot up to
            // z_inp_sponge_tmp[N-1] >= ztop. Now, interpolate to grid level 0 heights
            const int Ninp = z_inp_sponge_tmp.size();
            z_inp_sponge[0]     = zbot;
            U_inp_sponge[0]     = U_inp_sponge_tmp[0];
            V_inp_sponge[0]     = V_inp_sponge_tmp[0];
            for (int k=0; k < Nz; ++k) {
                z_inp_sponge[k+1] = (use_terrain) ? 0.5 * (zlevels_stag[k] + zlevels_stag[k+1])
                                                 : zbot + (k + 0.5) * dz;
                U_inp_sponge[k+1] = interpolate_1d(z_inp_sponge_tmp.dataPtr(), U_inp_sponge_tmp.dataPtr(), z_inp_sponge[k+1], Ninp);
                V_inp_sponge[k+1] = interpolate_1d(z_inp_sponge_tmp.dataPtr(), V_inp_sponge_tmp.dataPtr(), z_inp_sponge[k+1], Ninp);
            }
            z_inp_sponge[Nz+1]     = ztop;
            U_inp_sponge[Nz+1] = interpolate_1d(z_inp_sponge_tmp.dataPtr(), U_inp_sponge_tmp.dataPtr(), ztop, Ninp);
            V_inp_sponge[Nz+1] = interpolate_1d(z_inp_sponge_tmp.dataPtr(), V_inp_sponge_tmp.dataPtr(), ztop, Ninp);
        }

        amrex::Print() << "Successfully read the input_sponge file..." << std::endl;
        input_sponge_reader.close();
    }

    int size () const
    {
        AMREX_ALWAYS_ASSERT(z_inp_sponge.size() == U_inp_sponge.size());
        AMREX_ALWAYS_ASSERT(z_inp_sponge.size() == V_inp_sponge.size());
        return z_inp_sponge.size();
    }

    // Members

    std::string input_sponge_file = "input_sponge_file.txt";

    // - read from file
    amrex::Vector<amrex::Real> z_inp_sponge, U_inp_sponge, V_inp_sponge;
};
#endif
