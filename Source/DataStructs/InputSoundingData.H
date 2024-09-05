#ifndef _INPUT_SOUNDING_DATA_H_
#define _INPUT_SOUNDING_DATA_H_

#include <string>
#include <iostream>

#include <AMReX_ParmParse.H>
#include <AMReX_Print.H>
#include <AMReX_Gpu.H>
#include <AMReX_Geometry.H>

#include <EOS.H>
#include <ERF_Constants.H>
#include <Interpolation_1D.H>
#include <HSE_utils.H>

/**
 * Data structure storing input sounding data. Also
 * handles reading the input file for sounding data and
 * hydrostatic column integration.
 */
struct InputSoundingData {
public:
    InputSoundingData ()
    {
        amrex::ParmParse pp("erf");
        pp.query("tau_nudging", tau_nudging);

        // Read in input_sounding filename
        n_sounding_files = pp.countval("input_sounding_file");
        if (n_sounding_files > 0) {
            input_sounding_file.resize(n_sounding_files);
            pp.queryarr("input_sounding_file", input_sounding_file, 0, n_sounding_files);
        } else {
            n_sounding_files = 1;
            input_sounding_file.resize(n_sounding_files);
            input_sounding_file[0] = "input_sounding";
        }

        // Read in input_sounding times
        n_sounding_times = pp.countval("input_sounding_time");

        if (n_sounding_times > 0) {
            input_sounding_time.resize(n_sounding_times);
            pp.queryarr("input_sounding_time", input_sounding_time, 0, n_sounding_times);
        } else {
            n_sounding_times = 1;
            input_sounding_time.resize(n_sounding_times);
            input_sounding_time[0] = 0.0;
        }

        // If we have more files than times or times than files we just use the minimum
        int n = std::min(n_sounding_times, n_sounding_files);
        n_sounding_files = n;
        n_sounding_times = n;
        input_sounding_file.resize(n);
        input_sounding_time.resize(n);
    }

    void resize_arrays ()
    {
        ntimes = n_sounding_files;

        z_inp_sound.resize(ntimes);
        theta_inp_sound.resize(ntimes);
        qv_inp_sound.resize(ntimes);
        U_inp_sound.resize(ntimes);
        V_inp_sound.resize(ntimes);

        z_inp_sound_d.resize(ntimes);
        theta_inp_sound_d.resize(ntimes);
        qv_inp_sound_d.resize(ntimes);
        U_inp_sound_d.resize(ntimes);
        V_inp_sound_d.resize(ntimes);
    }

    void read_from_file (const amrex::Geometry &geom,
                         const amrex::Vector<amrex::Real>& zlevels_stag,
                         int itime)
    {
        const int klo = 0;
        const int khi = geom.Domain().bigEnd()[AMREX_SPACEDIM-1];
        const int Nz = geom.Domain().size()[AMREX_SPACEDIM-1];
        const amrex::Real dz = geom.CellSize()[AMREX_SPACEDIM-1];

        const bool use_terrain = (zlevels_stag.size() > 0);
        const amrex::Real zbot = (use_terrain) ? zlevels_stag[klo]   : geom.ProbLo(AMREX_SPACEDIM-1);
        const amrex::Real ztop = (use_terrain) ? zlevels_stag[khi+1] : geom.ProbHi(AMREX_SPACEDIM-1);

        z_inp_sound[itime].resize(Nz+2);
        theta_inp_sound[itime].resize(Nz+2);
        qv_inp_sound[itime].resize(Nz+2);
        U_inp_sound[itime].resize(Nz+2);
        V_inp_sound[itime].resize(Nz+2);

        // Read the input_sounding file
        amrex::Print() << "input_sounding file location : " << input_sounding_file[itime] << std::endl;
        std::ifstream input_sounding_reader(input_sounding_file[itime]);
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
            z_inp_sound[itime][0]     = zbot;
            theta_inp_sound[itime][0] = theta_inp_sound_tmp[0];
            qv_inp_sound[itime][0]    = qv_inp_sound_tmp[0];
            U_inp_sound[itime][0]     = U_inp_sound_tmp[0];
            V_inp_sound[itime][0]     = V_inp_sound_tmp[0];
            for (int k=0; k < Nz; ++k) {
                z_inp_sound[itime][k+1] = (use_terrain) ? 0.5 * (zlevels_stag[k] + zlevels_stag[k+1])
                                                 : zbot + (k + 0.5) * dz;
                theta_inp_sound[itime][k+1] = interpolate_1d(z_inp_sound_tmp.dataPtr(), theta_inp_sound_tmp.dataPtr(), z_inp_sound[itime][k+1], Ninp);
                   qv_inp_sound[itime][k+1] = interpolate_1d(z_inp_sound_tmp.dataPtr(),    qv_inp_sound_tmp.dataPtr(), z_inp_sound[itime][k+1], Ninp);
                    U_inp_sound[itime][k+1] = interpolate_1d(z_inp_sound_tmp.dataPtr(),     U_inp_sound_tmp.dataPtr(), z_inp_sound[itime][k+1], Ninp);
                    V_inp_sound[itime][k+1] = interpolate_1d(z_inp_sound_tmp.dataPtr(),     V_inp_sound_tmp.dataPtr(), z_inp_sound[itime][k+1], Ninp);
            }
            z_inp_sound[itime][Nz+1]     = ztop;
            theta_inp_sound[itime][Nz+1] = interpolate_1d(z_inp_sound_tmp.dataPtr(), theta_inp_sound_tmp.dataPtr(), ztop, Ninp);
               qv_inp_sound[itime][Nz+1] = interpolate_1d(z_inp_sound_tmp.dataPtr(),    qv_inp_sound_tmp.dataPtr(), ztop, Ninp);
                U_inp_sound[itime][Nz+1] = interpolate_1d(z_inp_sound_tmp.dataPtr(),     U_inp_sound_tmp.dataPtr(), ztop, Ninp);
                V_inp_sound[itime][Nz+1] = interpolate_1d(z_inp_sound_tmp.dataPtr(),     V_inp_sound_tmp.dataPtr(), ztop, Ninp);
        }

        amrex::Print() << "Successfully read the " << itime << "th input_sounding file..." << std::endl;
        input_sounding_reader.close();

        host_to_device(itime);
    }

    void calc_rho_p (int itime)
    {
        /* Calculate density and pressure, roughly following the procedure in
         * WRF dyn_em/module_initialize_ideal.F. We integrate hydrostatically
         * from the surface up through the air column to get the dry density
         * and moist pressure.
         */
        const amrex::Real tol = 1.0e-12;
        const int Ninp = size(itime);
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

#if 0   // Printing
        // In this absence of moisture, this moist profile will match the
        // following dry profile
        amrex::Print() << "z  p_m  rho_d  theta  qv  U  V" << std::endl;
        amrex::Print() << z_inp_sound[itime][0]
                       << " " << pm_integ[0]
                       << " " << rhod_integ[0]
                       << " " << theta_inp_sound[itime][0]
                       << " " << qv_inp_sound[itime][0]
                       << " " << U_inp_sound[itime][0]
                       << " " << V_inp_sound[itime][0]
                       << std::endl;
#endif

        // integrate from surface to domain top
        amrex::Real dz, F, C;
        amrex::Real rho_tot_hi, rho_tot_lo;
        for (int k=1; k < size(itime); ++k)
        {
            // Vertical grid spacing
            dz = z_inp_sound[itime][k] - z_inp_sound[itime][k-1];

            // Establish known constant
            rho_tot_lo = rhod_integ[k-1] * (1. + qv_inp_sound[itime][k-1]);
            C  = -pm_integ[k-1] + 0.5*rho_tot_lo*CONST_GRAV*dz;

            // Initial guess and residual
            pm_integ[k]   = pm_integ[k-1];
            rhod_integ[k] = getRhogivenThetaPress(theta_inp_sound[itime][k],
                                                  pm_integ[k],
                                                  R_d/Cp_d,
                                                  qv_inp_sound[itime][k]);
            rho_tot_hi = rhod_integ[k] * (1. + qv_inp_sound[itime][k]);
            F = pm_integ[k] + 0.5*rho_tot_hi*CONST_GRAV*dz + C;

            // Do iterations
            if (std::abs(F)>tol) HSEutils::Newton_Raphson_hse(tol, R_d/Cp_d, dz,
                                                              CONST_GRAV, C, theta_inp_sound[itime][k],
                                                              qv_inp_sound[itime][k], qv_inp_sound[itime][k],
                                                              pm_integ[k], rhod_integ[k], F);

#if 0       // Printing
            amrex::Print() << z_inp_sound[itime][k]
                           << " " << pm_integ[k]
                           << " " << rhod_integ[k]
                           << " " << theta_inp_sound[itime][k]
                           << " " << qv_inp_sound[itime][k]
                           << " " << U_inp_sound[itime][k]
                           << " " << V_inp_sound[itime][k]
                           << std::endl;
#endif
        }
        // Note: at this point, the surface pressure, density of the dry air
        // column is stored in pm_integ[0], rhod_integ[0]

        // update
        host_to_device(itime);
    }

    void host_to_device (int itime)
    {
        const int Ninp = size(itime);
        z_inp_sound_d[itime].resize(Ninp);
        theta_inp_sound_d[itime].resize(Ninp);
        qv_inp_sound_d[itime].resize(Ninp);
        U_inp_sound_d[itime].resize(Ninp);
        V_inp_sound_d[itime].resize(Ninp);

        amrex::Gpu::copy(amrex::Gpu::hostToDevice,
                         z_inp_sound[itime].begin(), z_inp_sound[itime].end(),
                         z_inp_sound_d[itime].begin());
        amrex::Gpu::copy(amrex::Gpu::hostToDevice,
                         theta_inp_sound[itime].begin(), theta_inp_sound[itime].end(),
                         theta_inp_sound_d[itime].begin());
        amrex::Gpu::copy(amrex::Gpu::hostToDevice,
                         qv_inp_sound[itime].begin(), qv_inp_sound[itime].end(),
                         qv_inp_sound_d[itime].begin());
        amrex::Gpu::copy(amrex::Gpu::hostToDevice,
                         U_inp_sound[itime].begin(), U_inp_sound[itime].end(),
                         U_inp_sound_d[itime].begin());
        amrex::Gpu::copy(amrex::Gpu::hostToDevice,
                         V_inp_sound[itime].begin(), V_inp_sound[itime].end(),
                         V_inp_sound_d[itime].begin());

        if (rhod_integ.size() > 0)
        {
            //amrex::Print() << "Copying rho_d, p_d to device" << std::endl;
            rho_inp_sound_d.resize(size(itime)+2);
            p_inp_sound_d.resize(size(itime)+2);
            amrex::Gpu::copy(amrex::Gpu::hostToDevice,
                             rhod_integ.begin(), rhod_integ.end(),
                             rho_inp_sound_d.begin());
            amrex::Gpu::copy(amrex::Gpu::hostToDevice,
                             pm_integ.begin(), pm_integ.end(),
                             p_inp_sound_d.begin());
        }
    }

    int size (int itime) const
    {
        AMREX_ALWAYS_ASSERT(z_inp_sound[itime].size() == theta_inp_sound[itime].size());
        AMREX_ALWAYS_ASSERT(z_inp_sound[itime].size() == qv_inp_sound[itime].size());
        AMREX_ALWAYS_ASSERT(z_inp_sound[itime].size() == U_inp_sound[itime].size());
        AMREX_ALWAYS_ASSERT(z_inp_sound[itime].size() == V_inp_sound[itime].size());
        return z_inp_sound[itime].size();
    }

    // Members
    int ntimes;

    amrex::Real tau_nudging = 5.0; // time scale used for nudging

    amrex::Vector<std::string> input_sounding_file = {};
    amrex::Vector<amrex::Real> input_sounding_time = {};
    int n_sounding_files = 0;
    int n_sounding_times = 0;

    // - read from file
    amrex::Real press_ref_inp_sound, theta_ref_inp_sound, qv_ref_inp_sound; // input

    // This is a vector (over time) of Vectors
    amrex::Vector<amrex::Vector<amrex::Real>> z_inp_sound, theta_inp_sound, qv_inp_sound, U_inp_sound, V_inp_sound;

    // This is a vector (over time) of DeviceVectors
    amrex::Vector<amrex::Gpu::DeviceVector<amrex::Real>> z_inp_sound_d, theta_inp_sound_d, qv_inp_sound_d, U_inp_sound_d, V_inp_sound_d;

    // - moist profiles
    amrex::Vector<amrex::Real> pm_integ; // from integrating up air column
    // - dry profiles
    amrex::Vector<amrex::Real> rhod_integ; // from integrating down air column
    // - to set solution fields
    amrex::Gpu::DeviceVector<amrex::Real> p_inp_sound_d, rho_inp_sound_d;
};
#endif
