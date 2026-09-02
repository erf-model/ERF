/**
 * @file ERF_DustBlastSchedule.H
 * @brief Time-scheduled blasting dust injection for the ERF-Dust emission layer.
 *
 * Implements timed mass injection of freshly blasted mineral particles into
 * the dust_emission_flux MultiFab on the 2D dust grid. Adapted from
 * ERF_IgnitionSchedule.H in Source/Fire/ which uses the same rank-0 CSV
 * read, MPI broadcast, and GPU ParallelFor stamp pattern.
 *
 * CSV file format (space or comma delimited, # or ! = comment):
 *   time_s  cx  cy  radius  mass_kg_m2  [mineral_type]  [priority]
 *
 *   time_s       — simulation time at which blast occurs [s]
 *   cx, cy       — blast center coordinates [m] in dust grid physical domain
 *   radius       — blast affected radius [m]
 *   mass_kg_m2   — dust mass per unit area injected per blast [kg/m^2]
 *   mineral_type — integer: 0=quartz tailings (default), 1=lithium brine,
 *                  2=REE tailings, 3=copper tailings (matches soil type codes)
 *   priority     — integer; higher value fires first within the same timestep
 *
 * The blast injects mass_kg_m2 uniformly into all cells within radius of
 * (cx, cy) by adding to dust_emission_flux. The total injection is divided
 * equally among all n_size_bins components.
 *
 * Blasting at a mine site generates freshly exposed mineral surfaces with
 * higher reactivity than weathered tailings. In Phase 4, PHREEQC provides
 * the composition for aged surfaces. For blast events, a simple fresh-surface
 * multiplier (blast_reactivity_factor) is applied to scale the injected flux.
 *
 * References:
 *   Finney, M.A. (2004). FARSITE: Fire Area Simulator. USDA Forest Service
 *     RMRS-RP-4. (ERF_IgnitionSchedule.H source pattern)
 *   U.S. EPA AP-42, Chapter 13.2.3 (Heavy Construction Operations).
 *     https://www.epa.gov/air-emissions-factors-and-quantification/
 *     ap-42-compilation-air-emission-factors
 *   MSHA 30 CFR Part 70 (Coal Dust Exposure Limits).
 *     https://www.ecfr.gov/current/title-30/chapter-I/subchapter-O/part-70
 */

#ifndef ERF_DUST_BLAST_SCHEDULE_H
#define ERF_DUST_BLAST_SCHEDULE_H

#ifdef ERF_USE_DUST

#include <AMReX_MultiFab.H>
#include <AMReX_MFIter.H>
#include <AMReX_Print.H>
#include <AMReX_ParallelDescriptor.H>
#include <vector>
#include <string>
#include <fstream>
#include <sstream>
#include <algorithm>
#include <cmath>

// Forward declaration of DustGrid
struct DustGrid;

namespace DustBlastConst {
    // Default blast reactivity multiplier for fresh surfaces.
    // Fresh mineral surfaces from blasting are more erodible than
    // weathered tailings. This factor scales the injected mass flux.
    constexpr amrex::Real BLAST_REACTIVITY_DEFAULT = 2.0;

    // Maximum injected mass per blast per cell per timestep [kg/m^2/s].
    // Prevents unreasonable values from a single blast event.
    constexpr amrex::Real BLAST_FLUX_MAX = 1.0e-2;
}

/**
 * @struct BlastEvent
 * @brief A single timed blasting event on the dust grid.
 *
 * Adapted from IgnitionEvent in ERF_IgnitionSchedule.H.
 */
struct BlastEvent
{
    amrex::Real time_s;          ///< Time at which blast occurs [s]
    amrex::Real cx;              ///< Blast center x [m]
    amrex::Real cy;              ///< Blast center y [m]
    amrex::Real radius;          ///< Affected radius [m]
    amrex::Real mass_kg_m2;      ///< Dust mass injected per unit area [kg/m^2]
    int         mineral_type;    ///< 0=quartz, 1=Li brine, 2=REE, 3=Cu
    int         priority;        ///< Higher = fires first in same timestep
    mutable bool fired;          ///< Set to true after event is applied

    BlastEvent()
        : time_s(0.0), cx(0.0), cy(0.0), radius(50.0),
          mass_kg_m2(1.0e-4), mineral_type(0), priority(5), fired(false)
    {}
};

/**
 * @struct BlastSchedule
 * @brief Container for a list of blast events loaded from a CSV file.
 *
 * Adapted from IgnitionSchedule in ERF_IgnitionSchedule.H.
 */
struct BlastSchedule
{
    std::vector<BlastEvent> events;
    bool empty() const { return events.empty(); }
};

/**
 * @brief Load blast schedule from a CSV file.
 *
 * Rank-0 reads CSV file; parses required fields time_s cx cy radius mass_kg_m2
 * and optional mineral_type priority. Comment lines start with # or !; 
 * commas converted to spaces before parsing. Sorts by time_s ascending,
 * then priority descending. Broadcasts all fields via MPI_Bcast using the 
 * same flat-array serialization pattern as ERF_IgnitionSchedule.H.
 *
 * @param[in]  filename  Path to CSV file.
 * @param[out] schedule  BlastSchedule to populate.
 */
inline
void load_blast_schedule(const std::string& filename,
                         BlastSchedule& schedule)
{
    schedule.events.clear();

    if (amrex::ParallelDescriptor::IOProcessor()) {
        std::ifstream file(filename);
        if (!file.is_open()) {
            amrex::Abort(std::string("Could not open blast schedule file: ") + filename);
        }

        std::string line;
        while (std::getline(file, line)) {
            // Skip comments and empty lines
            if (line.empty() || line[0] == '#' || line[0] == '!') {
                continue;
            }

            // Replace commas with spaces
            for (auto& c : line) {
                if (c == ',') c = ' ';
            }

            std::istringstream iss(line);
            amrex::Real time_s, cx, cy, radius, mass_kg_m2;
            int mineral_type = 0;
            int priority = 5;

            // Parse required fields
            if (!(iss >> time_s >> cx >> cy >> radius >> mass_kg_m2)) {
                continue;  // Skip malformed lines
            }

            // Parse optional fields
            iss >> mineral_type;
            iss >> priority;

            BlastEvent event;
            event.time_s = time_s;
            event.cx = cx;
            event.cy = cy;
            event.radius = radius;
            event.mass_kg_m2 = mass_kg_m2;
            event.mineral_type = mineral_type;
            event.priority = priority;
            event.fired = false;

            schedule.events.push_back(event);
        }

        file.close();

        // Sort by time ascending, then priority descending
        std::sort(schedule.events.begin(), schedule.events.end(),
                  [](const BlastEvent& a, const BlastEvent& b) {
                      if (std::abs(a.time_s - b.time_s) > 1.0e-14) {
                          return a.time_s < b.time_s;
                      }
                      return a.priority > b.priority;
                  });
    }

    // Broadcast event count
    int n_events = schedule.events.size();
    amrex::ParallelDescriptor::Bcast(&n_events, 1, amrex::ParallelDescriptor::IOProcessorNumber());

    // Resize on non-IO ranks
    if (!amrex::ParallelDescriptor::IOProcessor()) {
        schedule.events.resize(n_events);
    }

    // Broadcast event data via flat array
    for (auto& event : schedule.events) {
        amrex::ParallelDescriptor::Bcast(&event.time_s, 1, amrex::ParallelDescriptor::IOProcessorNumber());
        amrex::ParallelDescriptor::Bcast(&event.cx, 1, amrex::ParallelDescriptor::IOProcessorNumber());
        amrex::ParallelDescriptor::Bcast(&event.cy, 1, amrex::ParallelDescriptor::IOProcessorNumber());
        amrex::ParallelDescriptor::Bcast(&event.radius, 1, amrex::ParallelDescriptor::IOProcessorNumber());
        amrex::ParallelDescriptor::Bcast(&event.mass_kg_m2, 1, amrex::ParallelDescriptor::IOProcessorNumber());
        amrex::ParallelDescriptor::Bcast(&event.mineral_type, 1, amrex::ParallelDescriptor::IOProcessorNumber());
        amrex::ParallelDescriptor::Bcast(&event.priority, 1, amrex::ParallelDescriptor::IOProcessorNumber());
    }

    if (amrex::ParallelDescriptor::IOProcessor()) {
        amrex::Print() << "[DUST] Blast schedule loaded: " << n_events << " events\n";
    }
}

/**
 * @brief Apply blast events due in (prev_time, current_time] to dust_emission_flux.
 *
 * For each event in (prev_time, current_time] (sorted by priority descending,
 * time ascending), all dust grid cells within radius of (cx, cy) receive an
 * additional emission flux contribution:
 *
 *   dust_emission_flux(i, j, k, b) += blast_flux_per_bin
 *
 * where:
 *   blast_flux_per_bin = mass_kg_m2 * blast_reactivity_factor / (dt_s * n_size_bins)
 *   clamped to DustBlastConst::BLAST_FLUX_MAX per bin.
 *
 * This converts a total injected mass per unit area [kg/m^2] to a flux rate
 * [kg/m^2/s] over the current timestep duration dt_s.
 *
 * Adapted from apply_scheduled_ignitions in ERF_IgnitionSchedule.H.
 * Uses host-side event loop + GPU ParallelFor per event, same pattern.
 *
 * @param[in,out] dust_emission_flux  Emission flux MultiFab (ncomp=n_size_bins).
 * @param[in]     dg                  DustGrid geometry.
 * @param[in,out] schedule            BlastSchedule; fired flag set after apply.
 * @param[in]     current_time        Time at end of current step [s].
 * @param[in]     prev_time           Time at end of previous step [s].
 * @param[in]     dt_s                Current timestep duration [s].
 * @param[in]     n_size_bins         Number of size bins.
 * @param[in]     blast_reactivity    Reactivity multiplier for fresh surface [>= 1].
 */
inline
void apply_blast_schedule(
    amrex::MultiFab& dust_emission_flux,
    const DustGrid&  dg,
    BlastSchedule&   schedule,
    amrex::Real      current_time,
    amrex::Real      prev_time,
    amrex::Real      dt_s,
    int              n_size_bins,
    amrex::Real      blast_reactivity)
{
    // Collect events in (prev_time, current_time]
    std::vector<BlastEvent*> events_due;
    for (auto& evt : schedule.events) {
        if (evt.time_s > prev_time && evt.time_s <= current_time && !evt.fired) {
            events_due.push_back(&evt);
        }
    }

    // Sort by priority descending, time ascending (stable for same-time events)
    std::sort(events_due.begin(), events_due.end(),
              [](const BlastEvent* a, const BlastEvent* b) {
                  if (a->priority != b->priority) {
                      return a->priority > b->priority;
                  }
                  return a->time_s < b->time_s;
              });

    // Prevent division by zero and very small timesteps
    amrex::Real dt_safe = (dt_s > 1.0e-16) ? dt_s : 1.0e-16;

    // Get physical domain from dg
    amrex::Geometry const& dust_geom = dg.geom;
    auto const& plo = dust_geom.ProbLo();
    auto const& dx = dust_geom.CellSize();

    // Apply each event via GPU kernel
    for (auto evt : events_due) {
        amrex::Real blast_flux_per_bin = evt->mass_kg_m2 * blast_reactivity 
                                        / (dt_safe * n_size_bins);
        blast_flux_per_bin = std::min(blast_flux_per_bin, DustBlastConst::BLAST_FLUX_MAX);

        // GPU ParallelFor over dust grid valid region
#ifdef AMREX_USE_OMP
#pragma omp parallel if (amrex::Gpu::notInLaunchRegion())
#endif
        for (amrex::MFIter mfi(dust_emission_flux, true); mfi.isValid(); ++mfi) {
            const auto& box = mfi.tilebox();
            auto const& dest = dust_emission_flux.array(mfi);

            amrex::ParallelFor(box, n_size_bins,
                [=] AMREX_GPU_DEVICE(int i, int j, int k, int n) noexcept
                {
                    // Compute physical position of cell center
                    amrex::Real x_cell = plo[0] + (i + 0.5) * dx[0];
                    amrex::Real y_cell = plo[1] + (j + 0.5) * dx[1];

                    // Compute distance from blast center
                    amrex::Real dist_x = x_cell - evt->cx;
                    amrex::Real dist_y = y_cell - evt->cy;
                    amrex::Real dist = std::sqrt(dist_x * dist_x + dist_y * dist_y);

                    // If within radius, add blast flux
                    if (dist <= evt->radius) {
                        dest(i, j, k, n) += blast_flux_per_bin;
                    }
                });
        }

        // Mark event as fired and print status
        evt->fired = true;
        if (amrex::ParallelDescriptor::IOProcessor()) {
            amrex::Print() << "[DUST] Blast event at t=" << evt->time_s << " s, "
                          << "center=(" << evt->cx << ", " << evt->cy << "), "
                          << "radius=" << evt->radius << " m, "
                          << "mass=" << evt->mass_kg_m2 << " kg/m^2\n";
        }
    }
}

#endif // ERF_USE_DUST

#endif // ERF_DUST_BLAST_SCHEDULE_H
