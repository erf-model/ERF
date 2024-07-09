#ifndef _INPUT_SOUNDING_DATA_H_
#define _INPUT_SOUNDING_DATA_H_

#include <string>
#include <iostream>

#include <AMReX_ParmParse.H>
#include <AMReX_Print.H>
#include <AMReX_Gpu.H>
#include <AMReX_Geometry.H>

#include <ERF_Constants.H>
#include <Interpolation_1D.H>

/**
 * Data structure storing input sounding data. Also
 * handles reading the input file for sounding data and
 * hydrostatic column integration.
 */
struct InputSoundingData {
public:
    InputSoundingData () {}

    void read_from_file (const std::string input_sounding_file,
                         const amrex::Geometry &geom,
                         const amrex::Vector<amrex::Real>& zlevels_stag)
    {
        const int klo = 0;
        const int khi = geom.Domain().bigEnd()[AMREX_SPACEDIM-1];
        const int Nz = geom.Domain().size()[AMREX_SPACEDIM-1];
        const amrex::Real dz = geom.CellSize()[AMREX_SPACEDIM-1];

        const bool use_terrain = (zlevels_stag.size() > 0);
        const amrex::Real zbot = (use_terrain) ? zlevels_stag[klo]   : geom.ProbLo(AMREX_SPACEDIM-1);
        const amrex::Real ztop = (use_terrain) ? zlevels_stag[khi+1] : geom.ProbHi(AMREX_SPACEDIM-1);

        z_inp_sound.resize(Nz+2);
        theta_inp_sound.resize(Nz+2);
        qv_inp_sound.resize(Nz+2);
        U_inp_sound.resize(Nz+2);
        V_inp_sound.resize(Nz+2);

        // Read the input_sounding file
        amrex::Print() << "input_sounding file location : " << input_sounding_file << std::endl;
        std::ifstream input_sounding_reader(input_sounding_file);
        if(!input_sounding_reader.is_open()) {
            amrex::Error("Error opening the input_sounding file\n");
        }
        else {
            // Read the contents of the input_sounding file
            amrex::Print() << "Successfully opened the input_sounding file. Now reading... " << std::endl;
            std::string line;

            // First, read the input data into temp vectors; then, interpolate vectors to the
            // domain lo/hi and cell centers (from level 0)
            amrex::Vector<amrex::Real> z_inp_sound_tmp, theta_inp_sound_tmp, qv_inp_sound_tmp,
                                       U_inp_sound_tmp, V_inp_sound_tmp;

            // Read surface quantities from the first line
            std::getline(input_sounding_reader, line);
            std::istringstream iss(line);
            iss >> press_ref_inp_sound >> theta_ref_inp_sound >> qv_ref_inp_sound;
            press_ref_inp_sound *= 100; // convert from hPa to Pa
            qv_ref_inp_sound *= 0.001; // convert from g/kg to kg/kg

            // Add surface
            z_inp_sound_tmp.push_back(zbot); // height above sea level [m]
            theta_inp_sound_tmp.push_back(theta_ref_inp_sound);
            qv_inp_sound_tmp.push_back(qv_ref_inp_sound);
            U_inp_sound_tmp.push_back(0);
            V_inp_sound_tmp.push_back(0);

            // Read the vertical profile at each given height
            amrex::Real z, theta, qv, U, V;
            while(std::getline(input_sounding_reader, line)) {
                std::istringstream iss_z(line);
                iss_z >> z >> theta >> qv >> U >> V;
                if (z == zbot) {
                    AMREX_ALWAYS_ASSERT(theta == theta_inp_sound_tmp[0]);
                    AMREX_ALWAYS_ASSERT(qv*0.001 == qv_inp_sound_tmp[0]); // convert from g/kg to kg/kg
                    U_inp_sound_tmp[0] = U;
                    V_inp_sound_tmp[0] = V;
                } else {
                    AMREX_ALWAYS_ASSERT(z > z_inp_sound_tmp[z_inp_sound_tmp.size()-1]); // sounding is increasing in height
                    z_inp_sound_tmp.push_back(z);
                    theta_inp_sound_tmp.push_back(theta);
                    qv_inp_sound_tmp.push_back(qv*0.001); // convert from g/kg to kg/kg
                    U_inp_sound_tmp.push_back(U);
                    V_inp_sound_tmp.push_back(V);
                    if (z >= ztop) break;
                }
            }

            // At this point, we have an input_sounding from zbot up to
            // z_inp_sound_tmp[N-1] >= ztop. Now, interpolate to grid level 0 heights
            const int Ninp = z_inp_sound_tmp.size();
            z_inp_sound[0]     = zbot;
            theta_inp_sound[0] = theta_inp_sound_tmp[0];
            qv_inp_sound[0]    = qv_inp_sound_tmp[0];
            U_inp_sound[0]     = U_inp_sound_tmp[0];
            V_inp_sound[0]     = V_inp_sound_tmp[0];
            for (int k=0; k < Nz; ++k) {
                z_inp_sound[k+1] = (use_terrain) ? 0.5 * (zlevels_stag[k] + zlevels_stag[k+1])
                                                 : zbot + (k + 0.5) * dz;
                theta_inp_sound[k+1] = interpolate_1d(z_inp_sound_tmp.dataPtr(), theta_inp_sound_tmp.dataPtr(), z_inp_sound[k+1], Ninp);
                   qv_inp_sound[k+1] = interpolate_1d(z_inp_sound_tmp.dataPtr(),    qv_inp_sound_tmp.dataPtr(), z_inp_sound[k+1], Ninp);
                    U_inp_sound[k+1] = interpolate_1d(z_inp_sound_tmp.dataPtr(),     U_inp_sound_tmp.dataPtr(), z_inp_sound[k+1], Ninp);
                    V_inp_sound[k+1] = interpolate_1d(z_inp_sound_tmp.dataPtr(),     V_inp_sound_tmp.dataPtr(), z_inp_sound[k+1], Ninp);
            }
            z_inp_sound[Nz+1]     = ztop;
            theta_inp_sound[Nz+1] = interpolate_1d(z_inp_sound_tmp.dataPtr(), theta_inp_sound_tmp.dataPtr(), ztop, Ninp);
               qv_inp_sound[Nz+1] = interpolate_1d(z_inp_sound_tmp.dataPtr(),    qv_inp_sound_tmp.dataPtr(), ztop, Ninp);
                U_inp_sound[Nz+1] = interpolate_1d(z_inp_sound_tmp.dataPtr(),     U_inp_sound_tmp.dataPtr(), ztop, Ninp);
                V_inp_sound[Nz+1] = interpolate_1d(z_inp_sound_tmp.dataPtr(),     V_inp_sound_tmp.dataPtr(), ztop, Ninp);
        }

        amrex::Print() << "Successfully read the input_sounding file..." << std::endl;
        input_sounding_reader.close();

        host_to_device();
    }

    void calc_rho_p ()
    {
        /* Calculate density and pressure, roughly following the procedure in
         * WRF dyn_em/module_initialize_ideal.F. We integrate hydrostatically
         * from the surface up through the air column to get the dry density
         * and moist pressure.
         */
        const int maxiter = 10;
        const int Ninp = size();
        pm_integ.resize(Ninp);
        rhod_integ.resize(Ninp);

        // evaluate surface quantities (k=0): total pressure and dry air
          pm_integ[0] = press_ref_inp_sound;
        rhod_integ[0] = getRhogivenThetaPress(theta_ref_inp_sound,
                                              press_ref_inp_sound,
                                              R_d/Cp_d,
                                              qv_ref_inp_sound);

        amrex::Print() << "ideal sounding init: surface density of moist air = "
                       << rhod_integ[0] << " kg/m^3" << std::endl;

        // Note:
        //   p_dry = rho_d R_d T
        //   p_tot = rho_m R_d T_v
        //         = rho_d(1 + q_v) R_d T_v

        // integrate from surface to domain top
        amrex::Real qv, dz;
#if 1   // Printing
        // In this absence of moisture, this moist profile will match the
        // following dry profile
        amrex::Print() << "z  p_m  rho_d  theta  qv  U  V" << std::endl;
        amrex::Print() << z_inp_sound[0]
                       << " " << pm_integ[0]
                       << " " << rhod_integ[0]
                       << " " << theta_inp_sound[0]
                       << " " << qv_inp_sound[0]
                       << " " << U_inp_sound[0]
                       << " " << V_inp_sound[0]
                       << std::endl;
#endif
        for (int k=1; k < size(); ++k)
        {
            dz = z_inp_sound[k] - z_inp_sound[k-1];
            rhod_integ[k] = rhod_integ[k-1]; // guess
            for (int it=0; it < maxiter; ++it)
            {
                amrex::Real rho_tot_hi = rhod_integ[k  ] * (1. + qv_inp_sound[k  ]);
                amrex::Real rho_tot_lo = rhod_integ[k-1] * (1. + qv_inp_sound[k-1]);
                pm_integ[k] = pm_integ[k-1] - 0.5*dz*(rho_tot_hi + rho_tot_lo)*CONST_GRAV;
                AMREX_ALWAYS_ASSERT(pm_integ[k] > 0);
                rhod_integ[k] = getRhogivenThetaPress(theta_inp_sound[k],
                                                      pm_integ[k],
                                                      R_d/Cp_d,
                                                      qv_inp_sound[k]);
            }
#if 1       // Printing
            amrex::Print() << z_inp_sound[k]
                           << " " << pm_integ[k]
                           << " " << rhod_integ[k]
                           << " " << theta_inp_sound[k]
                           << " " << qv_inp_sound[k]
                           << " " << U_inp_sound[k]
                           << " " << V_inp_sound[k]
                           << std::endl;
#endif
        }
        // Note: at this point, the surface pressure, density of the dry air
        // column is stored in pm_integ[0], rhod_integ[0]

        // update
        host_to_device();
    }

    void host_to_device ()
    {
        const int Ninp = size();
        z_inp_sound_d.resize(Ninp);
        theta_inp_sound_d.resize(Ninp);
        qv_inp_sound_d.resize(Ninp);
        U_inp_sound_d.resize(Ninp);
        V_inp_sound_d.resize(Ninp);

        amrex::Gpu::copy(amrex::Gpu::hostToDevice,
                         z_inp_sound.begin(), z_inp_sound.end(),
                         z_inp_sound_d.begin());
        amrex::Gpu::copy(amrex::Gpu::hostToDevice,
                         theta_inp_sound.begin(), theta_inp_sound.end(),
                         theta_inp_sound_d.begin());
        amrex::Gpu::copy(amrex::Gpu::hostToDevice,
                         qv_inp_sound.begin(), qv_inp_sound.end(),
                         qv_inp_sound_d.begin());
        amrex::Gpu::copy(amrex::Gpu::hostToDevice,
                         U_inp_sound.begin(), U_inp_sound.end(),
                         U_inp_sound_d.begin());
        amrex::Gpu::copy(amrex::Gpu::hostToDevice,
                         V_inp_sound.begin(), V_inp_sound.end(),
                         V_inp_sound_d.begin());

        if (rhod_integ.size() > 0)
        {
            //amrex::Print() << "Copying rho_d, p_d to device" << std::endl;
            rho_inp_sound_d.resize(size()+2);
            p_inp_sound_d.resize(size()+2);
            amrex::Gpu::copy(amrex::Gpu::hostToDevice,
                             rhod_integ.begin(), rhod_integ.end(),
                             rho_inp_sound_d.begin());
            amrex::Gpu::copy(amrex::Gpu::hostToDevice,
                             pm_integ.begin(), pm_integ.end(),
                             p_inp_sound_d.begin());
        }
    }

    int size () const
    {
        AMREX_ALWAYS_ASSERT(z_inp_sound.size() == theta_inp_sound.size());
        AMREX_ALWAYS_ASSERT(z_inp_sound.size() == qv_inp_sound.size());
        AMREX_ALWAYS_ASSERT(z_inp_sound.size() == U_inp_sound.size());
        AMREX_ALWAYS_ASSERT(z_inp_sound.size() == V_inp_sound.size());
        return z_inp_sound.size();
    }

    // Members
    // - read from file
    amrex::Real press_ref_inp_sound, theta_ref_inp_sound, qv_ref_inp_sound; // input
    amrex::Vector<amrex::Real> z_inp_sound, theta_inp_sound, qv_inp_sound, U_inp_sound, V_inp_sound;
    amrex::Gpu::DeviceVector<amrex::Real> z_inp_sound_d, theta_inp_sound_d, qv_inp_sound_d, U_inp_sound_d, V_inp_sound_d;
    // - moist profiles
    amrex::Vector<amrex::Real> pm_integ; // from integrating up air column
    // - dry profiles
    amrex::Vector<amrex::Real> rhod_integ; // from integrating down air column
    // - to set solution fields
    amrex::Gpu::DeviceVector<amrex::Real> p_inp_sound_d, rho_inp_sound_d;
};
#endif
