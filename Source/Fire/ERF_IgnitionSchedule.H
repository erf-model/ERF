#ifndef ERF_IGNITION_SCHEDULE_H
#define ERF_IGNITION_SCHEDULE_H

/**
 * @file ERF_IgnitionSchedule.H
 * @brief Time-scheduled ignition events for multi-start fire scenarios
 *
 * Implements time-scheduled ignition events: events at specified times and
 * locations are stamped onto the level-set field phi using the same SDF
 * convention as the primary ignition. Multiple fire fronts merge implicitly
 * via the min(phi, new_val) stamping rule.
 *
 * CSV file format:
 *   time_s  cx  cy  radius  [source_type]  [priority]
 *   # or ! lines are comments
 *   source_type: "sphere" (default, integer 0) or "polyline" (integer 1)
 *   priority: integer; higher priority events fire first within the same timestep
 *   suppress_if_burning: if 1, skip cells already on fire (phi < 0)
 *
 * Merging of multiple fire fronts from different ignitions is implicit — the
 * FARSITE phi stamping uses min(phi, new_val) so fronts that overlap simply
 * share the phi field.
 *
 * Reference: Finney (2004) RMRS-RP-4
 */

#include <AMReX_MultiFab.H>
#include <AMReX_Geometry.H>
#include <AMReX_ParallelDescriptor.H>
#include <string>
#include <vector>
#include <cmath>
#include <algorithm>
#include <sstream>
#include <fstream>

/**
 * @struct IgnitionEvent
 * @brief A single ignition event at a specified time and location
 */
struct IgnitionEvent
{
    amrex::Real time_s;                 ///< Time at which event fires [seconds]
    amrex::Real cx;                     ///< Ignition center x [m]
    amrex::Real cy;                     ///< Ignition center y [m]
    amrex::Real radius;                 ///< Ignition sphere radius [m]
    int source_type;                    ///< 0 = sphere, 1 = polyline (unused in Phase 11)
    int priority;                       ///< Priority; higher = fires first in same timestep
    bool suppress_if_burning;           ///< If true, skip cells with phi < 0
    mutable bool fired;                 ///< True after event has been applied

    IgnitionEvent()
        : time_s(0.0), cx(0.0), cy(0.0), radius(30.0),
          source_type(0), priority(5), suppress_if_burning(false), fired(false)
    {}
};

/**
 * @struct IgnitionSchedule
 * @brief Container for a list of ignition events, loaded from CSV file
 */
struct IgnitionSchedule
{
    std::vector<IgnitionEvent> events;

    bool empty() const { return events.empty(); }
};

/**
 * @brief Load and parse a CSV schedule file.
 *
 * File is read only on the IO rank (amrex::ParallelDescriptor::IOProcessor()).
 * Events are parsed, sorted by time (ascending) then priority (descending),
 * and broadcast to all ranks via MPI_Bcast.
 *
 * CSV format:
 *   time_s  cx  cy  radius  [source_type]  [priority]  [suppress_if_burning]
 *   # or ! lines are comments; blank lines are ignored
 *   commas and spaces are treated as delimiters; all commas converted to spaces
 *
 * If the file is empty at startup (which is not an error if primary ignition
 * is via ignition_x/y/r), a warning is printed and the schedule remains empty.
 *
 * @param[in]  filename     Path to CSV schedule file
 * @param[out] schedule     IgnitionSchedule to populate
 */
inline
void load_ignition_schedule(const std::string& filename,
                            IgnitionSchedule& schedule)
{
    using namespace amrex;

    schedule.events.clear();

    std::vector<IgnitionEvent> host_events;

    if (ParallelDescriptor::IOProcessor()) {
        std::ifstream file(filename);
        if (!file.is_open()) {
            amrex::Abort("Cannot open ignition schedule file: " + filename);
        }

        std::string line;
        int line_num = 0;

        while (std::getline(file, line)) {
            ++line_num;

            // Skip empty lines
            if (line.empty()) continue;

            // Skip comment lines
            if (line[0] == '#' || line[0] == '!') continue;

            // Convert commas to spaces
            for (auto& c : line) {
                if (c == ',') c = ' ';
            }

            // Parse line
            std::istringstream iss(line);
            IgnitionEvent event;

            if (!(iss >> event.time_s >> event.cx >> event.cy >> event.radius)) {
                amrex::Print() << "[WARNING] Line " << line_num
                               << " in " << filename
                               << " could not be parsed; skipping\n";
                continue;
            }

            // Optional fields
            event.source_type = 0;           // default: sphere
            event.priority = 5;              // default: medium priority
            event.suppress_if_burning = false;  // default: apply everywhere

            iss >> event.source_type;        // may or may not be present
            iss >> event.priority;           // may or may not be present
            int suppress_int = 0;
            iss >> suppress_int;
            event.suppress_if_burning = (suppress_int != 0);

            event.fired = false;             // not yet fired

            host_events.push_back(event);
        }

        file.close();

        if (host_events.empty()) {
            amrex::Print() << "[WARNING] Ignition schedule file '" << filename
                           << "' is empty or contains only comments.\n";
        }

        // Sort by time (ascending), then by priority (descending)
        std::sort(host_events.begin(), host_events.end(),
                  [](const IgnitionEvent& a, const IgnitionEvent& b) {
                      if (std::abs(a.time_s - b.time_s) > 1.0e-10_rt) {
                          return a.time_s < b.time_s;
                      }
                      return a.priority > b.priority;  // higher priority first
                  });

        if (amrex::ParallelDescriptor::NProcs() > 1) {
            // Broadcast schedule size
            int n_events = static_cast<int>(host_events.size());
            MPI_Bcast(&n_events, 1, MPI_INT, amrex::ParallelDescriptor::IOProcessorNumber(),
                      MPI_COMM_WORLD);

            if (n_events > 0) {
                // Serialize events into flat arrays for broadcast
                std::vector<amrex::Real> d_time_s(n_events);
                std::vector<amrex::Real> d_cx(n_events);
                std::vector<amrex::Real> d_cy(n_events);
                std::vector<amrex::Real> d_radius(n_events);
                std::vector<int> d_source_type(n_events);
                std::vector<int> d_priority(n_events);
                std::vector<int> d_suppress(n_events);

                for (int i = 0; i < n_events; ++i) {
                    d_time_s[i] = host_events[i].time_s;
                    d_cx[i] = host_events[i].cx;
                    d_cy[i] = host_events[i].cy;
                    d_radius[i] = host_events[i].radius;
                    d_source_type[i] = host_events[i].source_type;
                    d_priority[i] = host_events[i].priority;
                    d_suppress[i] = host_events[i].suppress_if_burning ? 1 : 0;
                }

                MPI_Bcast(d_time_s.data(), n_events, MPI_DOUBLE,
                          amrex::ParallelDescriptor::IOProcessorNumber(), MPI_COMM_WORLD);
                MPI_Bcast(d_cx.data(), n_events, MPI_DOUBLE,
                          amrex::ParallelDescriptor::IOProcessorNumber(), MPI_COMM_WORLD);
                MPI_Bcast(d_cy.data(), n_events, MPI_DOUBLE,
                          amrex::ParallelDescriptor::IOProcessorNumber(), MPI_COMM_WORLD);
                MPI_Bcast(d_radius.data(), n_events, MPI_DOUBLE,
                          amrex::ParallelDescriptor::IOProcessorNumber(), MPI_COMM_WORLD);
                MPI_Bcast(d_source_type.data(), n_events, MPI_INT,
                          amrex::ParallelDescriptor::IOProcessorNumber(), MPI_COMM_WORLD);
                MPI_Bcast(d_priority.data(), n_events, MPI_INT,
                          amrex::ParallelDescriptor::IOProcessorNumber(), MPI_COMM_WORLD);
                MPI_Bcast(d_suppress.data(), n_events, MPI_INT,
                          amrex::ParallelDescriptor::IOProcessorNumber(), MPI_COMM_WORLD);

                host_events.resize(n_events);
                for (int i = 0; i < n_events; ++i) {
                    host_events[i].time_s = d_time_s[i];
                    host_events[i].cx = d_cx[i];
                    host_events[i].cy = d_cy[i];
                    host_events[i].radius = d_radius[i];
                    host_events[i].source_type = d_source_type[i];
                    host_events[i].priority = d_priority[i];
                    host_events[i].suppress_if_burning = (d_suppress[i] != 0);
                    host_events[i].fired = false;
                }
            }
        }
    } else {
        // Non-IO ranks participate in broadcast
        if (amrex::ParallelDescriptor::NProcs() > 1) {
            int n_events = 0;
            MPI_Bcast(&n_events, 1, MPI_INT, amrex::ParallelDescriptor::IOProcessorNumber(),
                      MPI_COMM_WORLD);

            if (n_events > 0) {
                std::vector<amrex::Real> d_time_s(n_events);
                std::vector<amrex::Real> d_cx(n_events);
                std::vector<amrex::Real> d_cy(n_events);
                std::vector<amrex::Real> d_radius(n_events);
                std::vector<int> d_source_type(n_events);
                std::vector<int> d_priority(n_events);
                std::vector<int> d_suppress(n_events);

                MPI_Bcast(d_time_s.data(), n_events, MPI_DOUBLE,
                          amrex::ParallelDescriptor::IOProcessorNumber(), MPI_COMM_WORLD);
                MPI_Bcast(d_cx.data(), n_events, MPI_DOUBLE,
                          amrex::ParallelDescriptor::IOProcessorNumber(), MPI_COMM_WORLD);
                MPI_Bcast(d_cy.data(), n_events, MPI_DOUBLE,
                          amrex::ParallelDescriptor::IOProcessorNumber(), MPI_COMM_WORLD);
                MPI_Bcast(d_radius.data(), n_events, MPI_DOUBLE,
                          amrex::ParallelDescriptor::IOProcessorNumber(), MPI_COMM_WORLD);
                MPI_Bcast(d_source_type.data(), n_events, MPI_INT,
                          amrex::ParallelDescriptor::IOProcessorNumber(), MPI_COMM_WORLD);
                MPI_Bcast(d_priority.data(), n_events, MPI_INT,
                          amrex::ParallelDescriptor::IOProcessorNumber(), MPI_COMM_WORLD);
                MPI_Bcast(d_suppress.data(), n_events, MPI_INT,
                          amrex::ParallelDescriptor::IOProcessorNumber(), MPI_COMM_WORLD);

                host_events.resize(n_events);
                for (int i = 0; i < n_events; ++i) {
                    host_events[i].time_s = d_time_s[i];
                    host_events[i].cx = d_cx[i];
                    host_events[i].cy = d_cy[i];
                    host_events[i].radius = d_radius[i];
                    host_events[i].source_type = d_source_type[i];
                    host_events[i].priority = d_priority[i];
                    host_events[i].suppress_if_burning = (d_suppress[i] != 0);
                    host_events[i].fired = false;
                }
            }
        }
    }

    schedule.events = host_events;
}

/**
 * @brief Apply scheduled ignitions that are due within the given time window.
 *
 * Events are applied if their time_s falls in the interval (prev_time, current_time].
 * Within a single timestep, events are ordered by priority (higher first), then by
 * time (earlier first). Each event is stamped onto phi using:
 *   phi(i,j,k) = min(phi(i,j,k), -(radius - dist))
 *
 * This SDF stamp convention (same as the primary ignition) ensures:
 * - Cells within radius of the ignition center become negative (burned)
 * - Cells outside the radius remain unaffected if they are positive (unburned)
 * - Multiple ignitions merge implicitly via min() accumulation
 *
 * With suppress_if_burning=true, cells where phi < 0 (already burning) are not
 * modified; only positive cells (unburned) receive the stamp.
 *
 * The fired flag is set to true after successful application to prevent re-firing.
 *
 * @param[in,out] phi         Level-set field (on fire grid)
 * @param[in]     geom_fire   Fire grid geometry
 * @param[in,out] schedule    IgnitionSchedule with events to apply
 * @param[in]     current_time Time at end of current step [seconds]
 * @param[in]     prev_time   Time at end of previous step [seconds]
 */
inline
void apply_scheduled_ignitions(
    amrex::MultiFab& phi,
    const amrex::Geometry& geom_fire,
    IgnitionSchedule& schedule,
    amrex::Real current_time,
    amrex::Real prev_time)
{
    using namespace amrex;

    const auto* prob_lo = geom_fire.ProbLo();
    const auto* dx      = geom_fire.CellSize();

    // Collect events that should fire in this timestep, in priority order
    std::vector<IgnitionEvent*> events_to_fire;

    for (auto& event : schedule.events) {
        // Skip if already fired
        if (event.fired) continue;

        // Check if time falls in (prev_time, current_time]
        if (event.time_s > prev_time && event.time_s <= current_time) {
            events_to_fire.push_back(&event);
        }
    }

    // Sort by priority (higher first), then by time (earlier first)
    std::sort(events_to_fire.begin(), events_to_fire.end(),
              [](const IgnitionEvent* a, const IgnitionEvent* b) {
                  if (a->priority != b->priority) {
                      return a->priority > b->priority;
                  }
                  return a->time_s < b->time_s;
              });

    // Apply each event
    for (auto* evt : events_to_fire) {
        Real cx = evt->cx;
        Real cy = evt->cy;
        Real radius = evt->radius;
        bool suppress = evt->suppress_if_burning;

        for (MFIter mfi(phi); mfi.isValid(); ++mfi) {
            const Box& bx = mfi.validbox();
            auto p = phi.array(mfi);

            ParallelFor(bx, [=] AMREX_GPU_DEVICE (int i, int j, int k) noexcept {
                Real x = prob_lo[0] + (i + 0.5_rt) * dx[0];
                Real y = prob_lo[1] + (j + 0.5_rt) * dx[1];
                Real dist = std::sqrt((x - cx) * (x - cx) + (y - cy) * (y - cy));

                Real new_phi = -(radius - dist);
                // Clamp to [-1, 1] to match primary ignition convention
                new_phi = amrex::max(amrex::min(new_phi, 1.0_rt), -1.0_rt);

                // Apply suppress_if_burning logic
                if (suppress && p(i, j, k) < 0.0_rt) {
                    // Cell is already on fire; skip it
                    return;
                }

                // Stamp using min() semantics: existing fire is never overwritten
                p(i, j, k) = amrex::min(p(i, j, k), new_phi);
            });
        }

        evt->fired = true;

        if (ParallelDescriptor::IOProcessor()) {
            amrex::Print() << "[FIRE] Ignition event at t=" << evt->time_s
                           << " s, center=(" << cx << ", " << cy
                           << "), radius=" << radius << " m\n";
        }
    }
}

#endif // ERF_IGNITION_SCHEDULE_H
