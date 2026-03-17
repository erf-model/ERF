#ifndef ERF_LARGE_SCALE_FORCING_DATA_H_
#define ERF_LARGE_SCALE_FORCING_DATA_H_

#include <string>
#include <iostream>

#include <AMReX_ParmParse.H>
#include <AMReX_Print.H>
#include <AMReX_Gpu.H>
#include <AMReX_Geometry.H>

#include <ERF_EOS.H>
#include <ERF_Constants.H>
#include <ERF_Interpolation_1D.H>
#include <ERF_HSEUtils.H>
#include <ERF_InputSoundingData.H>

/**
 * Data structure storing large scale forcing data.
 */
struct LargeScaleForcingData {
public:
    LargeScaleForcingData ()
    {
        amrex::ParmParse pp("erf");
        pp.query("forcing_timescale", tau_lsf);
        pp.query("large_scale_forcing_file", lsf_file);


    };

    int get_columns_in_str(std::string &str)
    {
        std::istringstream iss(str);
        std::string tmp;
        // Get the number of columns in the file
        int ncols = 0;
        while (iss >> tmp) {
            ncols+= 1;
        }
        return ncols;
    }


    void read_forcing_file()
    {
        amrex::Print() << "  Opening LSF file '" << lsf_file << "'" << std::endl;
        std::ifstream ifs(lsf_file);
        if (!ifs.is_open())
        {
            amrex::Error("Error opening input forcing file " + lsf_file);
        }

        std::string line;

        // temporary vectors for storing the values at each time
        amrex::Vector<amrex::Real> z_in, p_in, t_in, q_in, u_in, v_in, w_in;

        //   skip first header line
        std::getline(ifs, line);
        // header should have 7 columns: z, p, tls, qls, uls, vls, wls
        AMREX_ALWAYS_ASSERT(get_columns_in_str(line) == 7);

        num_times = 0;
        nz_lsf = 0;

        amrex::Real start_time = 0.0; // first time in file, in days, other inputs are relative to this time

        amrex::Real pres0;
        while(std::getline(ifs, line)) {
            std::istringstream iss(line);

            // The start of each day has 6 columns: day, levels, pres0, "day, levels, pres0"
            int ncol = get_columns_in_str(line);

            if (ncol == 6)
            {
                // this is a new set of inputs at a given time

                // have to read these as strings first to split on commas
                std::string s_time, s_lev, s_pres0;
                iss >> s_time >> s_lev >> s_pres0;
                amrex::Real time = std::stod(s_time);
                int nlev = std::stoi(s_lev);
                pres0 = std::stod(s_pres0);

                // convert time to seconds (relative to first time), and pres from mb to Pa
                if (num_times == 0)
                {
                    start_time = time;
                }
                time = (time - start_time) * 86400.0;

                if (times_lsf.size() > 0 && time != times_lsf.back())
                {
                    // if we started a new time, store data for previous time and reset tmp arrays
                    AMREX_ALWAYS_ASSERT(z_in.size() == nz_lsf + 1);
                    z_lsf.push_back(z_in);
                    p_lsf.push_back(p_in);
                    t_lsf.push_back(t_in);
                    q_lsf.push_back(q_in);
                    u_lsf.push_back(u_in);
                    v_lsf.push_back(v_in);
                    w_lsf.push_back(w_in);

                    z_in.clear();
                    p_in.clear();
                    t_in.clear();
                    q_in.clear();
                    u_in.clear();
                    v_in.clear();
                    w_in.clear();
                }

                times_lsf.push_back(time);
                num_times++;

                pres0 *= 100.0;
                if (nz_lsf == 0)
                {
                    nz_lsf = nlev;
                } else {
                    // make sure this has the same number of inputs.
                    // TODO: this probably doesn't have to be true since we are interpolating onto the grid?
                    AMREX_ALWAYS_ASSERT(nlev == nz_lsf);
                }

                z_in.push_back(0.0); // assuming zbot = 0.0 for now, this will get interpolated later once we have grid info
                p_in.push_back(pres0);
                t_in.push_back(0.0); // assuming temperature change at surf is 0
                q_in.push_back(0.0);
                u_in.push_back(0.0);
                v_in.push_back(0.0);
                w_in.push_back(0.0);
            }
            else if (ncol == 7)
            {
                // this is defining a level at the current time
                amrex::Real z, p, tls, qls, uls, vls, wls;
                iss >> z >> p >> tls >> qls >> uls >> vls >> wls;

                /*
                if (z_in.size() == 1)
                {
                    u_in[0] = uls;
                    v_in[0] = vls;
                    w_in[0] = wls;
                    //t_in[0] = getThgivenPandT(tls, p*100.0, R_d/Cp_d);
                    //t_in[0] = tls*std::pow(pres0/(p*100.0), R_d/Cp_d);
                    t_in[0] = tls;
                    q_in[0] = qls;

                    z_in[0] = 0.0;
                }
                */
                if (z < 0.0) {
                    // if z < 0, then LSF is defined on pressure levels
                    use_z_grid = false;
                }

                z_in.push_back(z);
                p_in.push_back(p * 100.0);
                //t_in.push_back(getThgivenPandT(tls, p*100.0, R_d/Cp_d));
                //t_in.push_back(tls*std::pow(pres0/(p*100.0), R_d/Cp_d));
                t_in.push_back(tls);
                q_in.push_back(qls);
                u_in.push_back(uls);
                v_in.push_back(vls);
                w_in.push_back(wls);
            }
            else {
                amrex::Error("unexpected line in input file: " + line);
            }
        }

        ifs.close();

        // store the last set from the EOF
        AMREX_ALWAYS_ASSERT(z_in.size() == nz_lsf + 1);
        z_lsf.push_back(z_in);
        p_lsf.push_back(p_in);
        t_lsf.push_back(t_in);
        q_lsf.push_back(q_in);
        u_lsf.push_back(u_in);
        v_lsf.push_back(v_in);
        w_lsf.push_back(w_in);

        num_times = times_lsf.size();
        amrex::Print() << "  Read " << nz_lsf << " levels at " << num_times << " times" << std::endl;

        // debugging print
        for (int i = 0; i < times_lsf.size(); i++)
        {
            amrex::Print() << " LSF at time = " << times_lsf[i] << std::endl;
            for (int lev = 0; lev < nz_lsf + 1; lev++)
            {
                amrex::Print() << "   level " << lev << ": z = " << z_lsf[i][lev] << " p = " << p_lsf[i][lev] << " tls = " << t_lsf[i][lev] << " qls = " << q_lsf[i][lev] << " uls = " << u_lsf[i][lev] << " vls = " << v_lsf[i][lev] << " wls = " << w_lsf[i][lev] << std::endl;
            }
            break;
        }
    }


    void
    interp_forcing(const amrex::GeometryData &geom,
                   const amrex::Vector<amrex::Real>& zlevels_stag,
                   const InputSoundingData &sounding)
    {
        const int klo = 0;
        const int khi = geom.Domain().bigEnd()[AMREX_SPACEDIM-1];
        const int Nz = geom.Domain().size()[AMREX_SPACEDIM-1];
        const amrex::Real dz = geom.CellSize()[AMREX_SPACEDIM-1];

        const bool use_terrain = (zlevels_stag.size() > 0);
        const amrex::Real zbot = (use_terrain) ? zlevels_stag[klo]   : geom.ProbLo(AMREX_SPACEDIM-1);
        const amrex::Real ztop = (use_terrain) ? zlevels_stag[khi+1] : geom.ProbHi(AMREX_SPACEDIM-1);

        int num_times = times_lsf.size();
        AMREX_ALWAYS_ASSERT(num_times == num_times);

        z_int_lsf.resize(num_times);
        t_int_lsf.resize(num_times);
        q_int_lsf.resize(num_times);
        u_int_lsf.resize(num_times);
        v_int_lsf.resize(num_times);
        w_int_lsf.resize(num_times);

        z_int_lsf_d.resize(num_times);
        t_int_lsf_d.resize(num_times);
        q_int_lsf_d.resize(num_times);
        u_int_lsf_d.resize(num_times);
        v_int_lsf_d.resize(num_times);
        w_int_lsf_d.resize(num_times);
        for (int itime = 0; itime < num_times; itime++)
        {
            z_int_lsf[itime].resize(Nz+1);
            t_int_lsf[itime].resize(Nz+1);
            q_int_lsf[itime].resize(Nz+1);
            u_int_lsf[itime].resize(Nz+1);
            v_int_lsf[itime].resize(Nz+1);
            w_int_lsf[itime].resize(Nz+1);

            for (int k = 0; k < Nz; k++)
            {
                bool set = false;
                if (use_z_grid) {
                    for (int lk = 2; lk < nz_lsf + 1; lk++)
                    {
                        z_int_lsf[itime][k] = (use_terrain) ? 0.5 * (zlevels_stag[k] + zlevels_stag[k+1])
                                                            : zbot + (k + 0.5) * dz;
                        if (z_int_lsf[itime][k] <= z_lsf[itime][lk])
                        {
                            amrex::Real coef = (z_int_lsf[itime][k] - z_lsf[itime][lk - 1]) / (z_lsf[itime][lk] - z_lsf[itime][lk - 1]);
                            t_int_lsf[itime][k] = t_lsf[itime][lk - 1] + (t_lsf[itime][lk] - t_lsf[itime][lk-1])*coef;
                            q_int_lsf[itime][k] = q_lsf[itime][lk - 1] + (q_lsf[itime][lk] - q_lsf[itime][lk-1])*coef;
                            u_int_lsf[itime][k] = u_lsf[itime][lk - 1] + (u_lsf[itime][lk] - u_lsf[itime][lk-1])*coef;
                            v_int_lsf[itime][k] = v_lsf[itime][lk - 1] + (v_lsf[itime][lk] - v_lsf[itime][lk-1])*coef;
                            w_int_lsf[itime][k] = w_lsf[itime][lk - 1] + (w_lsf[itime][lk] - w_lsf[itime][lk-1])*coef;
                            set = true;
                            break;
                        }
                    }
                } else {
                    for (int lk = 2; lk < nz_lsf + 1; lk++)
                    {
                        z_int_lsf[itime][k] = (use_terrain) ? 0.5 * (zlevels_stag[k] + zlevels_stag[k+1])
                                                            : zbot + (k + 0.5) * dz;
                        if (sounding.pm_integ[k] >= p_lsf[itime][lk])
                        {
                            amrex::Real coef = (sounding.pm_integ[k] - p_lsf[itime][lk - 1]) / (p_lsf[itime][lk] - p_lsf[itime][lk - 1]);
                            t_int_lsf[itime][k] = t_lsf[itime][lk - 1] + (t_lsf[itime][lk] - t_lsf[itime][lk-1])*coef;
                            q_int_lsf[itime][k] = q_lsf[itime][lk - 1] + (q_lsf[itime][lk] - q_lsf[itime][lk-1])*coef;
                            u_int_lsf[itime][k] = u_lsf[itime][lk - 1] + (u_lsf[itime][lk] - u_lsf[itime][lk-1])*coef;
                            v_int_lsf[itime][k] = v_lsf[itime][lk - 1] + (v_lsf[itime][lk] - v_lsf[itime][lk-1])*coef;
                            w_int_lsf[itime][k] = w_lsf[itime][lk - 1] + (w_lsf[itime][lk] - w_lsf[itime][lk-1])*coef;
                            set = true;
                            break;
                        }
                    }
                }

                if (set) {
                    set = false;
                    continue;
                }

                t_int_lsf[itime][k] = 0.0;
                q_int_lsf[itime][k] = 0.0;
                u_int_lsf[itime][k] = u_int_lsf[itime][k-1];
                v_int_lsf[itime][k] = v_int_lsf[itime][k-1];
                w_int_lsf[itime][k] = 0.0;
            }

            z_int_lsf[itime][Nz] = ztop;
            t_int_lsf[itime][Nz] = 0.0;
            q_int_lsf[itime][Nz] = 0.0;
            u_int_lsf[itime][Nz] = u_int_lsf[itime][Nz-1];
            v_int_lsf[itime][Nz] = v_int_lsf[itime][Nz-1];
            w_int_lsf[itime][Nz] = 0.0;

            host_to_device(itime);
        }

        // debugging print
        for (int i = 0; i < times_lsf.size(); i++)
        {
            amrex::Print() << " INTERPOLATED LSF at time = " << times_lsf[i] << std::endl;
            for (int k = 0; k <= Nz; k++)
            {
                amrex::Print() << "   level " << k << ": z = " << z_int_lsf[i][k] << " tls = " << t_int_lsf[i][k] << " qls = " << q_int_lsf[i][k] << " uls = " << u_int_lsf[i][k] << " vls = " << v_int_lsf[i][k] << " wls = " << w_int_lsf[i][k] << std::endl;
            }
            break;
        }
    }

    void
    get_forcing_time_coeffs(const amrex::Real& time,
                            int& curr,
                            int& next,
                            amrex::Real& coeff_curr,
                            amrex::Real& coeff_next)
    {
        // helper function to get the current time indices of forcing data to use and
        // coefficients based on the current simulation time
        curr = 0;
        next = 0;
        coeff_curr = 1.0;
        coeff_next = 0.0;

        // find the time index of forcing data to use
        for (int nt = 1; nt < num_times; nt++)
        {
            if (time - start_time > times_lsf[nt])
            {
                curr = nt;
            }
        }

        if (curr == num_times - 1)
        {
            // use last set if time > last forcing time
            next = num_times;
        } else {
            next = curr + 1;
            coeff_next = ((time - start_time) - times_lsf[curr]) / (times_lsf[next] - times_lsf[curr]);
            coeff_curr = (1.0 - coeff_next);
        }
    }

    void host_to_device (int itime)
    {
        const int nz = z_int_lsf[itime].size();
        z_int_lsf_d[itime].resize(nz);
        t_int_lsf_d[itime].resize(nz);
        q_int_lsf_d[itime].resize(nz);
        u_int_lsf_d[itime].resize(nz);
        v_int_lsf_d[itime].resize(nz);
        w_int_lsf_d[itime].resize(nz);

        amrex::Gpu::copy(amrex::Gpu::hostToDevice,
                         z_int_lsf[itime].begin(), z_int_lsf[itime].end(),
                         z_int_lsf_d[itime].begin());
        amrex::Gpu::copy(amrex::Gpu::hostToDevice,
                         t_int_lsf[itime].begin(), t_int_lsf[itime].end(),
                         t_int_lsf_d[itime].begin());
        amrex::Gpu::copy(amrex::Gpu::hostToDevice,
                         q_int_lsf[itime].begin(), q_int_lsf[itime].end(),
                         q_int_lsf_d[itime].begin());
        amrex::Gpu::copy(amrex::Gpu::hostToDevice,
                         u_int_lsf[itime].begin(), u_int_lsf[itime].end(),
                         u_int_lsf_d[itime].begin());
        amrex::Gpu::copy(amrex::Gpu::hostToDevice,
                         v_int_lsf[itime].begin(), v_int_lsf[itime].end(),
                         v_int_lsf_d[itime].begin());
        amrex::Gpu::copy(amrex::Gpu::hostToDevice,
                         w_int_lsf[itime].begin(), w_int_lsf[itime].end(),
                         w_int_lsf_d[itime].begin());
    }

    amrex::Real start_time = 0.0; // simulation start time
    amrex::Real tau_lsf = 0.0; // timescale to apply lsf
    std::string lsf_file = "";

    bool use_z_grid = true;

    // inputs from large scale forcing file
    int nz_lsf = -1;    // number of z levels in the LSF file
    int num_times = -1; // number of times found in the LSF file
    amrex::Vector<amrex::Real> times_lsf;
    // original data from the file
    amrex::Vector<amrex::Vector<amrex::Real>> z_lsf;
    amrex::Vector<amrex::Vector<amrex::Real>> p_lsf;
    amrex::Vector<amrex::Vector<amrex::Real>> t_lsf;
    amrex::Vector<amrex::Vector<amrex::Real>> q_lsf;
    amrex::Vector<amrex::Vector<amrex::Real>> u_lsf;
    amrex::Vector<amrex::Vector<amrex::Real>> v_lsf;
    amrex::Vector<amrex::Vector<amrex::Real>> w_lsf;

    // interpolated data from file
    amrex::Vector<amrex::Vector<amrex::Real>> z_int_lsf;
    amrex::Vector<amrex::Vector<amrex::Real>> p_int_lsf;
    amrex::Vector<amrex::Vector<amrex::Real>> t_int_lsf;
    amrex::Vector<amrex::Vector<amrex::Real>> q_int_lsf;
    amrex::Vector<amrex::Vector<amrex::Real>> u_int_lsf;
    amrex::Vector<amrex::Vector<amrex::Real>> v_int_lsf;
    amrex::Vector<amrex::Vector<amrex::Real>> w_int_lsf;

    amrex::Vector<amrex::Gpu::DeviceVector<amrex::Real>> z_int_lsf_d, t_int_lsf_d, q_int_lsf_d, u_int_lsf_d, v_int_lsf_d, w_int_lsf_d;

};
#endif
