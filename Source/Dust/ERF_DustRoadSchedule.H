/**
 * @file ERF_DustRoadSchedule.H
 * @brief Haul road vehicle schedule emission for mine dust.
 *
 * Reads a CSV road schedule file and applies per-road vehicle
 * resuspension emission to dust_emission_flux each timestep.
 *
 * Emission formula (EPA AP-42 Chapter 13.2.2, unpaved roads):
 *   E_road [g/VKT] = k * (s/12)^a * (W/3)^b
 *   k=2.6, a=0.8, b=0.4 (TSP, AP-42 Table 13.2.2-2)
 *   F_vehicle [kg/m^2/s] = E_road * vmt_per_h / (road_width * road_length * 1e6 * 3600)
 *
 * CSV format:
 *   road_name, x_lo_m, y_lo_m, x_hi_m, y_hi_m,
 *   road_width_m, vehicle_weight_t, silt_pct, vmt_per_h,
 *   start_time_s, end_time_s
 *   Lines beginning with # or ! are comments.
 *   end_time_s = -1 means active for full simulation.
 *
 * Pattern: ERF_DustBlastSchedule.H (Phase 7) and ERF_IgnitionSchedule.H:
 *   https://github.com/hgopalan/ERF/blob/f264ce7dd0f5e727c71c939659a6a5caf001d21b/
 *   Source/Fire/ERF_IgnitionSchedule.H
 *
 * References:
 *   EPA AP-42 Ch. 13.2.2 (Unpaved Roads):
 *     https://www.epa.gov/air-emissions-factors-and-quantification/
 *     ap-42-compilation-air-emissions-factors
 *   Abu-Allaban et al. (2003). Atmos. Environ., 37, 1417.
 *     https://doi.org/10.1016/S1352-2310(02)01015-8
 *   Watson et al. (1994). Atmos. Environ., 28, 2493.
 *     https://doi.org/10.1016/1352-2310(94)90400-6
 */
#ifdef ERF_USE_DUST
#ifndef ERF_DUST_ROAD_SCHEDULE_H_
#define ERF_DUST_ROAD_SCHEDULE_H_

#include <AMReX_MultiFab.H>
#include <AMReX_Geometry.H>
#include <AMReX_ParallelDescriptor.H>
#include <AMReX_Vector.H>
#include <fstream>
#include <sstream>
#include <string>
#include <algorithm>
#include <cmath>

namespace DustRoadConst {
    constexpr amrex::Real K_AP42   = 2.6;   ///< AP-42 emission factor k [g/VKT]
    constexpr amrex::Real A_AP42   = 0.8;   ///< AP-42 silt exponent a
    constexpr amrex::Real B_AP42   = 0.4;   ///< AP-42 weight exponent b
    constexpr amrex::Real FLUX_MAX = 1.0e-3;///< Maximum road flux [kg/m^2/s]
    constexpr amrex::Real G_PER_KG = 1.0e-3;///< g -> kg
    constexpr amrex::Real S_PER_H  = 3600.0;///< s/h
}

struct RoadEvent {
    std::string name;
    amrex::Real x_lo, y_lo, x_hi, y_hi; ///< Road bounding box [m]
    amrex::Real road_width_m;            ///< Road width [m]
    amrex::Real vehicle_weight_t;        ///< Mean vehicle weight [tons]
    amrex::Real silt_pct;                ///< Silt content [%]
    amrex::Real vmt_per_h;              ///< Vehicle miles/km travelled per hour
    amrex::Real start_time_s;           ///< Activation start time [s]
    amrex::Real end_time_s;             ///< Activation end time [s]. -1 = whole sim.
};

struct RoadSchedule {
    amrex::Vector<RoadEvent> roads;
    bool loaded = false;
};

/**
 * @brief Load road schedule from CSV file.
 *
 * Rank-0 reads file and broadcasts to all ranks.
 * Pattern: load_blast_schedule in ERF_DustBlastSchedule.H (Phase 7).
 */
inline void load_road_schedule(const std::string& filename,
                                RoadSchedule& schedule)
{
    schedule.roads.clear();
    schedule.loaded = false;
    if (filename.empty()) return;

    amrex::Vector<amrex::Real> buf;
    amrex::Vector<std::string> names;
    int n_roads = 0;

    if (amrex::ParallelDescriptor::IOProcessor()) {
        std::ifstream in(filename);
        if (!in.is_open()) {
            amrex::Print() << "[DUST] WARNING: Cannot open road schedule: "
                           << filename << "\n";
            return;
        }
        std::string line;
        while (std::getline(in, line)) {
            if (line.empty() || line[0] == '#' || line[0] == '!') continue;
            std::replace(line.begin(), line.end(), ',', ' ');
            std::istringstream iss(line);
            RoadEvent ev;
            iss >> ev.name >> ev.x_lo >> ev.y_lo >> ev.x_hi >> ev.y_hi
                >> ev.road_width_m >> ev.vehicle_weight_t >> ev.silt_pct
                >> ev.vmt_per_h >> ev.start_time_s >> ev.end_time_s;
            if (!iss.fail()) {
                names.push_back(ev.name);
                buf.push_back(ev.x_lo); buf.push_back(ev.y_lo);
                buf.push_back(ev.x_hi); buf.push_back(ev.y_hi);
                buf.push_back(ev.road_width_m); buf.push_back(ev.vehicle_weight_t);
                buf.push_back(ev.silt_pct); buf.push_back(ev.vmt_per_h);
                buf.push_back(ev.start_time_s); buf.push_back(ev.end_time_s);
                ++n_roads;
            }
        }
    }
    amrex::ParallelDescriptor::Bcast(&n_roads, 1,
        amrex::ParallelDescriptor::IOProcessorNumber());
    if (n_roads == 0) { schedule.loaded = true; return; }

    buf.resize(n_roads * 10);
    names.resize(n_roads, "road");
    amrex::ParallelDescriptor::Bcast(buf.data(), (int)buf.size(),
        amrex::ParallelDescriptor::IOProcessorNumber());

    schedule.roads.resize(n_roads);
    for (int r = 0; r < n_roads; ++r) {
        auto& ev = schedule.roads[r];
        ev.name            = names[r];
        ev.x_lo            = buf[r*10+0]; ev.y_lo = buf[r*10+1];
        ev.x_hi            = buf[r*10+2]; ev.y_hi = buf[r*10+3];
        ev.road_width_m    = buf[r*10+4]; ev.vehicle_weight_t = buf[r*10+5];
        ev.silt_pct        = buf[r*10+6]; ev.vmt_per_h = buf[r*10+7];
        ev.start_time_s    = buf[r*10+8]; ev.end_time_s = buf[r*10+9];
    }
    schedule.loaded = true;

    if (amrex::ParallelDescriptor::IOProcessor())
        amrex::Print() << "[DUST] Road schedule loaded: " << n_roads
                       << " roads from " << filename << "\n";
}

/**
 * @brief Apply active road emissions to dust_emission_flux (additive).
 *
 * For each active road (start <= cur_time <= end), identifies dust grid
 * cells in the road bounding box and adds F_vehicle to dust_emission_flux.
 * Pattern: apply_blast_schedule in ERF_DustBlastSchedule.H (Phase 7).
 *
 * @param flux       dust_emission_flux MultiFab (component 0, additive).
 * @param geom_dust  Dust grid geometry.
 * @param schedule   Loaded road schedule.
 * @param cur_time   Current simulation time [s].
 * @param dt         Timestep [s] (not used in flux; kept for consistency).
 * @param debug      Print per-road activation if true.
 * @param diag_file  CSV file for road diagnostics (appended each active step).
 * @param nstep      Current step number.
 */
inline void apply_road_schedule(
    amrex::MultiFab&       flux,
    const amrex::Geometry& geom_dust,
    const RoadSchedule&    schedule,
    amrex::Real            cur_time,
    amrex::Real            dt,
    bool                   debug,
    const std::string&     diag_file,
    int                    nstep)
{
    // Documented above as kept for signature consistency with the other schedule
    // appliers; the road flux is a rate and does not integrate over the step.
    amrex::ignore_unused(dt);

    if (!schedule.loaded || schedule.roads.empty()) return;

    // GpuArray by value: CellSize()/ProbLo() return pointers into host Geometry
    // members, and the road-stamping kernel below dereferences these.
    const auto dxd = geom_dust.CellSizeArray();
    const auto lod = geom_dust.ProbLoArray();
    bool any_active = false;

    for (const auto& ev : schedule.roads) {
        bool active = (cur_time >= ev.start_time_s) &&
                      (ev.end_time_s < 0.0 || cur_time <= ev.end_time_s);
        if (!active) continue;
        any_active = true;

        // EPA AP-42 Ch. 13.2.2 emission factor [g/VKT]
        amrex::Real E_road = DustRoadConst::K_AP42
            * std::pow(ev.silt_pct / 12.0, DustRoadConst::A_AP42)
            * std::pow(ev.vehicle_weight_t / 3.0, DustRoadConst::B_AP42);

        // Road length and flux conversion.
        amrex::Real road_len = std::max(ev.x_hi - ev.x_lo, ev.y_hi - ev.y_lo);
        if (road_len <= 0.0) continue;
        // F_vehicle [kg/m^2/s]
        amrex::Real F_v = E_road * DustRoadConst::G_PER_KG
                        * ev.vmt_per_h
                        / (ev.road_width_m * road_len * DustRoadConst::S_PER_H);
        F_v = amrex::min(F_v, DustRoadConst::FLUX_MAX);

        const amrex::Real x_lo = ev.x_lo, y_lo = ev.y_lo;
        const amrex::Real x_hi = ev.x_hi, y_hi = ev.y_hi;

        for (amrex::MFIter mfi(flux, amrex::TilingIfNotGPU()); mfi.isValid(); ++mfi) {
            const amrex::Box& bx = mfi.tilebox();
            auto fl = flux.array(mfi);
            amrex::ParallelFor(bx, [=] AMREX_GPU_DEVICE(int i, int j, int k) noexcept {
                amrex::Real xc = lod[0] + (i + 0.5) * dxd[0];
                amrex::Real yc = lod[1] + (j + 0.5) * dxd[1];
                if (xc >= x_lo && xc <= x_hi && yc >= y_lo && yc <= y_hi)
                    fl(i,j,k) += F_v;
            });
        }

        if (debug) {
            amrex::Print() << "[DUST DEBUG] Phase 22: road=" << ev.name
                           << " t=" << cur_time << " s"
                           << " F_v=" << F_v << " kg/m^2/s\n";
        }
    }

    // Write road diagnostics CSV.
    if (any_active && !diag_file.empty()) {
        if (!amrex::ParallelDescriptor::IOProcessor()) return;
        bool file_exists = std::ifstream(diag_file).good();
        std::ofstream out(diag_file, std::ios::app);
        if (!file_exists)
            out << "# EPA AP-42 Chapter 13.2.2 road emission diagnostics\n"
                << "step,time_s,road_name,flux_add_kg_m2_s\n";
        for (const auto& ev : schedule.roads) {
            bool active = (cur_time >= ev.start_time_s) &&
                          (ev.end_time_s < 0.0 || cur_time <= ev.end_time_s);
            if (!active) continue;
            amrex::Real E_road = DustRoadConst::K_AP42
                * std::pow(ev.silt_pct/12.0, DustRoadConst::A_AP42)
                * std::pow(ev.vehicle_weight_t/3.0, DustRoadConst::B_AP42);
            amrex::Real road_len = std::max(ev.x_hi-ev.x_lo, ev.y_hi-ev.y_lo);
            amrex::Real F_v = (road_len > 0.0)
                ? amrex::min(E_road * DustRoadConst::G_PER_KG * ev.vmt_per_h
                             / (ev.road_width_m * road_len * DustRoadConst::S_PER_H),
                             DustRoadConst::FLUX_MAX)
                : 0.0;
            out << nstep << "," << cur_time << "," << ev.name << ","
                << F_v << "\n";
        }
    }
}
#endif
#endif
