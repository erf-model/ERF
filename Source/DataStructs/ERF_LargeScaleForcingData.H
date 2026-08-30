#ifndef ERF_LARGE_SCALE_FORCING_DATA_H_
#define ERF_LARGE_SCALE_FORCING_DATA_H_

#include <string>
#include <iostream>
#include <sstream>
#include <fstream>

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
 * Data structure storing large-scale forcing data.
 *
 * This provides time-varying vertical profiles for large scale temperature,
 * and water vapor tendencies, as well as observed large scale zonal and meridonal wind,
 * and vertical subsidence.
 * Profiles are read from a text file, interpolated to the model vertical grid, and then copied
 * to device memory for use during time integration.
 *
 * This mirrors the lsf functionality as implemented in SAM.
 */
struct LargeScaleForcingData {
public:
    /**
     * @brief Construct large-scale forcing metadata from the ERF input namespace.
     */
    LargeScaleForcingData ()
    {
        amrex::ParmParse pp("erf");
        pp.queryAdd("forcing_timescale", tau_lsf);
        pp.queryAdd("large_scale_forcing_file", lsf_file);
    }

    /**
     * @brief Count the number of whitespace-separated columns in a line.
     * @param str Input line to inspect.
     * @return Number of columns in the line.
     */
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


    /**
     * @brief Read the forcing file and cache the raw time slices on host.
     */
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

        amrex::Real file_start_time = zero; // first time in file, in days, other inputs are relative to this time

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
                    file_start_time = time;
                }
                time = (time - file_start_time) * amrex::Real(86400.0);

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

                pres0 *= amrex::Real(100.0);
                if (nz_lsf == 0)
                {
                    nz_lsf = nlev;
                } else {
                    // make sure this has the same number of inputs.
                    // TODO: this probably doesn't have to be true since we are interpolating onto the grid?
                    AMREX_ALWAYS_ASSERT(nlev == nz_lsf);
                }

                z_in.push_back(zero); // assuming zbot = 0.0 for now, this will get interpolated later once we have grid info
                p_in.push_back(pres0);
                t_in.push_back(zero); // assuming temperature change at surf is 0
                q_in.push_back(zero);
                u_in.push_back(zero);
                v_in.push_back(zero);
                w_in.push_back(zero);
            }
            else if (ncol == 7)
            {
                // this is defining a level at the current time
                amrex::Real z, p, tls, qls, uls, vls, wls;
                iss >> z >> p >> tls >> qls >> uls >> vls >> wls;

                if (z < zero) {
                    // if z < 0, then LSF is defined on pressure levels
                    use_z_grid = false;
                }

                z_in.push_back(z);
                p_in.push_back(p * amrex::Real(100.0));
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
        AMREX_ALWAYS_ASSERT_WITH_MESSAGE(z_lsf.size() == num_times,
                                         "repeated timestamp in the large-scale forcing file");
        amrex::Print() << "  Read " << nz_lsf << " levels at " << num_times << " times" << std::endl;

        // debugging print
        if (verbose_print && !times_lsf.empty()) {
            const int i = 0; // print only the first time block
            amrex::Print() << " LSF at time = " << times_lsf[i] << std::endl;
            for (int lev = 0; lev < nz_lsf + 1; lev++)
            {
                amrex::Print() << "   level " << lev << ": z = " << z_lsf[i][lev] << " p = " << p_lsf[i][lev] << " tls = "
                               << t_lsf[i][lev] << " qls = " << q_lsf[i][lev] << " uls = " << u_lsf[i][lev] << " vls = "
                               << v_lsf[i][lev] << " wls = " << w_lsf[i][lev] << std::endl;
            }
        }
    }


    /**
     * @brief Interpolate the forcing profiles onto the model grid.
     * @param geom Geometry defining the vertical grid.
     * @param zlevels_stag Staggered vertical coordinates used for interpolation.
     * @param sounding Input sounding data used when the forcing file is on pressure levels.
     */
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

        AMREX_ALWAYS_ASSERT(num_times == static_cast<int>(times_lsf.size()));

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
                z_int_lsf[itime][k] = (use_terrain) ? myhalf * (zlevels_stag[k] + zlevels_stag[k+1])
                                                    : zbot + (k + myhalf) * dz;

                bool set = false;
                for (int lk = 2; lk < nz_lsf + 1; lk++)
                {
                    // Determine if grid cell k falls in bracket [lk-1, lk]:
                    //   z-grid mode: cell z <= forcing z (z increases with lk)
                    //   p-grid mode: cell p >= forcing p (p decreases with lk)
                    bool in_bracket = use_z_grid
                        ? (z_int_lsf[itime][k] <= z_lsf[itime][lk])
                        : (sounding.pm_integ[k] >= p_lsf[itime][lk]);

                    if (in_bracket)
                    {
                        // Interpolation coordinate: z or p depending on mode
                        const auto& coord_lsf = use_z_grid ? z_lsf[itime] : p_lsf[itime];
                        amrex::Real coord_k = use_z_grid ? z_int_lsf[itime][k] : sounding.pm_integ[k];
                        amrex::Real coef = (coord_k - coord_lsf[lk - 1]) / (coord_lsf[lk] - coord_lsf[lk - 1]);

                        t_int_lsf[itime][k] = t_lsf[itime][lk - 1] + (t_lsf[itime][lk] - t_lsf[itime][lk-1])*coef;
                        q_int_lsf[itime][k] = q_lsf[itime][lk - 1] + (q_lsf[itime][lk] - q_lsf[itime][lk-1])*coef;
                        u_int_lsf[itime][k] = u_lsf[itime][lk - 1] + (u_lsf[itime][lk] - u_lsf[itime][lk-1])*coef;
                        v_int_lsf[itime][k] = v_lsf[itime][lk - 1] + (v_lsf[itime][lk] - v_lsf[itime][lk-1])*coef;
                        w_int_lsf[itime][k] = w_lsf[itime][lk - 1] + (w_lsf[itime][lk] - w_lsf[itime][lk-1])*coef;
                        set = true;
                        break;
                    }
                }

                if (set) {
                    continue;
                }

                t_int_lsf[itime][k] = zero;
                q_int_lsf[itime][k] = zero;
                u_int_lsf[itime][k] = (k > 0) ? u_int_lsf[itime][k-1] : zero;
                v_int_lsf[itime][k] = (k > 0) ? v_int_lsf[itime][k-1] : zero;
                w_int_lsf[itime][k] = zero;
            }

            z_int_lsf[itime][Nz] = ztop;
            t_int_lsf[itime][Nz] = zero;
            q_int_lsf[itime][Nz] = zero;
            u_int_lsf[itime][Nz] = u_int_lsf[itime][Nz-1];
            v_int_lsf[itime][Nz] = v_int_lsf[itime][Nz-1];
            w_int_lsf[itime][Nz] = zero;

            host_to_device(itime);
        }

        // debugging print
        if (verbose_print && !times_lsf.empty()) {
            const int i = 0; // print only the first time block
            amrex::Print() << " INTERPOLATED LSF at time = " << times_lsf[i] << std::endl;
            for (int k = 0; k <= Nz; k++)
            {
                amrex::Print() << "   level " << k << ": z = " << z_int_lsf[i][k] << " tls = " << t_int_lsf[i][k]
                               << " qls = " << q_int_lsf[i][k] << " uls = " << u_int_lsf[i][k] << " vls = "
                               << v_int_lsf[i][k] << " wls = " << w_int_lsf[i][k] << std::endl;
            }
        }
    }

    /**
     * @brief Determine interpolation coefficients to apply tendencies for the given time
     * @param time Elapsed simulation time in seconds.
     * @param curr Index of the lower bracketing forcing time.
     * @param next Index of the upper bracketing forcing time.
     * @param coeff_curr Interpolation weight for the lower time slice.
     * @param coeff_next Interpolation weight for the upper time slice.
     */
    void
    get_forcing_time_coeffs(const amrex::Real& time,
                            int& curr,
                            int& next,
                            amrex::Real& coeff_curr,
                            amrex::Real& coeff_next)
    {
        // helper function to get the current time indices of forcing data to use and
        // coefficients based on the current simulation time.
        curr = 0;
        next = 0;
        coeff_curr = one;
        coeff_next = zero;

        // find the time index of forcing data to use
        for (int nt = 1; nt < num_times; nt++)
        {
            if (time > times_lsf[nt])
            {
                curr = nt;
            }
        }

        if (curr == num_times - 1)
        {
            // use last set if time > last forcing time
            next = num_times - 1;
            coeff_curr = one;
            coeff_next = zero;
        } else {
            next = curr + 1;
            coeff_next = (time - times_lsf[curr]) / (times_lsf[next] - times_lsf[curr]);
            coeff_curr = (one - coeff_next);
        }
    }

    /**
     * @brief Copy one interpolated forcing slice from host memory to device memory.
     * @param itime Time-slice index to copy.
     */
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

    /**
     * Relaxation time scale used by the large-scale forcing and nudging terms.
     */
    amrex::Real tau_lsf = zero;
    /**
     * Path to the large-scale forcing input file.
     */
    std::string lsf_file = "";

    /**
     * True when the forcing file is tabulated on height levels, false when it is
     * tabulated on pressure levels.
     */
    bool use_z_grid = true;

    /**
     * Whether to display the lsf at initialization (before and after interpolation)
     */
    bool verbose_print = false;

    /**
     * Number of vertical levels and times found in the forcing file.
     */
    int nz_lsf = -1;
    int num_times = -1;
    /**
     * Forcing times read from the file, converted to elapsed seconds.
     */
    amrex::Vector<amrex::Real> times_lsf;
    /**
     * Raw forcing profiles read from the input file, indexed by time then level.
     */
    amrex::Vector<amrex::Vector<amrex::Real>> z_lsf;
    amrex::Vector<amrex::Vector<amrex::Real>> p_lsf;
    amrex::Vector<amrex::Vector<amrex::Real>> t_lsf;
    amrex::Vector<amrex::Vector<amrex::Real>> q_lsf;
    amrex::Vector<amrex::Vector<amrex::Real>> u_lsf;
    amrex::Vector<amrex::Vector<amrex::Real>> v_lsf;
    amrex::Vector<amrex::Vector<amrex::Real>> w_lsf;

    /**
     * Forcing profiles interpolated onto the model grid, indexed by time then level.
     */
    amrex::Vector<amrex::Vector<amrex::Real>> z_int_lsf;
    amrex::Vector<amrex::Vector<amrex::Real>> t_int_lsf;
    amrex::Vector<amrex::Vector<amrex::Real>> q_int_lsf;
    amrex::Vector<amrex::Vector<amrex::Real>> u_int_lsf;
    amrex::Vector<amrex::Vector<amrex::Real>> v_int_lsf;
    amrex::Vector<amrex::Vector<amrex::Real>> w_int_lsf;

    /**
     * Device copies of the interpolated forcing profiles.
     */
    amrex::Vector<amrex::Gpu::DeviceVector<amrex::Real>> z_int_lsf_d, t_int_lsf_d, q_int_lsf_d, u_int_lsf_d, v_int_lsf_d, w_int_lsf_d;

};
#endif
