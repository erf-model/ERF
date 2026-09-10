#ifndef ERF_CLOUD_CHAMBER_BUDGET_H_
#define ERF_CLOUD_CHAMBER_BUDGET_H_

#include <AMReX_Geometry.H>
#include <AMReX_MultiFab.H>
#include <AMReX_ParallelDescriptor.H>
#include <ERF_IndexDefines.H>

#include <array>
#include <algorithm>
#include <cmath>
#include <fstream>
#include <iomanip>
#include <limits>
#include <string>

/**
 * Low-overhead Stage 1 conserved-scalar budget accumulator.
 *
 * Fluxes are positive in the positive coordinate direction.  The reported
 * boundary balance is therefore low-minus-high for each coordinate pair.
 * Local face sums are held until a configured report interval, then one
 * packed MPI reduction produces the six-face values.  For the anelastic RK2
 * path, nrk=0 is retained and nrk=1 applies the trapezoidal stage weight;
 * this avoids counting both RHS evaluations as two physical timesteps.
 */
class CloudChamberBudget
{
public:
    enum Scalar : int { RhoTheta = 0, RhoQv = 1, RhoQc = 2, NumScalars = 3 };
    static constexpr int NumFaces = 2 * AMREX_SPACEDIM;

    explicit CloudChamberBudget (int interval) : m_interval(interval) {}

    bool enabled () const noexcept { return m_interval > 0; }
    /**
     * Determine if a budget report is due for the given timestep.
     *
     * @param step Current timestep
     * @return True if a report is due, false otherwise
     */
    bool due (int step) const noexcept {
        return enabled() && step > 0 && (step % m_interval == 0);
    }

    /**
     * Initialize the base state for the budget accumulation period.
     *
     * @param cons Conserved variable MultiFab
     * @param geom Geometry used for volume calculation
     * @param step Timestep at the start of the interval
     * @param time Simulation time at the start of the interval
     */
    void set_initial_state (const amrex::MultiFab& cons,
                            const amrex::Geometry& geom,
                            int step = 0,
                            double time = 0.0)
    {
        if (m_have_initial) { return; }
        const amrex::Real volume = geom.CellSize(0) * geom.CellSize(1) * geom.CellSize(2);
        for (int s = 0; s < NumScalars; ++s) {
            const int comp = (s == RhoTheta) ? RhoTheta_comp : RhoQ1_comp + (s-1);
            m_initial_state[s] = (comp < cons.nComp()) ? cons.sum(comp, true) * volume : amrex::Real(0.0);
        }
        m_interval_start_step = step;
        m_interval_start_time = time;
        m_have_initial = true;
    }

    /**
     * Integrate boundary fluxes for a specific Runge-Kutta stage.
     *
     * @param scalar Index of the scalar being tracked
     * @param nrk Runge-Kutta stage index
     * @param dt Timestep size
     * @param xflux Fluxes in the x-direction
     * @param yflux Fluxes in the y-direction
     * @param zflux Fluxes in the z-direction
     * @param geom Geometry used for face area calculation
     * @param flux_comp Component index of the flux
     */
    void capture_stage (int scalar, int nrk, amrex::Real dt,
                        const amrex::MultiFab& xflux,
                        const amrex::MultiFab& yflux,
                        const amrex::MultiFab& zflux,
                        const amrex::Geometry& geom,
                        int flux_comp = 0)
    {
        if (!enabled() || scalar < 0 || scalar >= NumScalars) { return; }
        std::array<amrex::Real, NumFaces> current{};
        for (int d = 0; d < AMREX_SPACEDIM; ++d) {
            amrex::Box low = geom.Domain();
            amrex::Box high = low;
            low.surroundingNodes(d);
            high.surroundingNodes(d);
            low.setSmall(d, geom.Domain().smallEnd(d));
            low.setBig(d, geom.Domain().smallEnd(d));
            high.setSmall(d, geom.Domain().bigEnd(d)+1);
            high.setBig(d, geom.Domain().bigEnd(d)+1);
            const amrex::MultiFab* flux = (d == 0) ? &xflux : ((d == 1) ? &yflux : &zflux);
            const amrex::Real area = (d == 0) ? geom.CellSize(1)*geom.CellSize(2) :
                                     ((d == 1) ? geom.CellSize(0)*geom.CellSize(2) :
                                                 geom.CellSize(0)*geom.CellSize(1));
            current[2*d]   = flux->sum(low, flux_comp, true) * area * dt;
            current[2*d+1] = flux->sum(high, flux_comp, true) * area * dt;
        }

        if (nrk == 0) {
            m_stage0[scalar] = current;
            m_have_stage0[scalar] = true;
        } else if (nrk == 1 && m_have_stage0[scalar]) {
            for (int face = 0; face < NumFaces; ++face) {
                m_cumulative[scalar][face] +=
                    amrex::Real(0.5) * (m_stage0[scalar][face] + current[face]);
            }
            m_have_stage0[scalar] = false;
        }
    }

    std::array<amrex::Real, NumScalars>
    state_integrals (const amrex::MultiFab& cons,
                     const amrex::Geometry& geom) const
    {
        const amrex::Real volume = geom.CellSize(0) * geom.CellSize(1) * geom.CellSize(2);
        std::array<amrex::Real, NumScalars> result{};
        result[RhoTheta] = cons.sum(RhoTheta_comp, true) * volume;
        result[RhoQv] = (RhoQ1_comp < cons.nComp()) ? cons.sum(RhoQ1_comp, true) * volume : amrex::Real(0.0);
        result[RhoQc] = (RhoQ2_comp < cons.nComp()) ? cons.sum(RhoQ2_comp, true) * volume : amrex::Real(0.0);
        return result;
    }

    /**
     * Accumulate internal source terms for the tracked scalars.
     *
     * @param before State of scalars before source application
     * @param after State of scalars after source application
     */
    void record_internal_source (const std::array<amrex::Real, NumScalars>& before,
                                 const std::array<amrex::Real, NumScalars>& after)
    {
        for (int s = 0; s < NumScalars; ++s) {
            m_internal_source[s] += after[s] - before[s];
        }
    }

    /**
     * Perform global reduction and write the budget report to file.
     *
     * @param step Current timestep
     * @param time Current simulation time
     * @param cons Current conserved variable MultiFab
     * @param geom Geometry used for state integration
     * @param cloudy Flag indicating if the state is cloudy
     */
    void report (int step, double time, const amrex::MultiFab& cons,
                 const amrex::Geometry& geom, bool cloudy)
    {
        if (!due(step) || !m_have_initial) { return; }
        const auto current = state_integrals(cons, geom);
        std::array<amrex::Real, NumScalars*NumFaces + NumScalars*3> packed{};
        int index = 0;
        for (int s = 0; s < NumScalars; ++s) {
            for (int f = 0; f < NumFaces; ++f) { packed[index++] = m_cumulative[s][f]; }
        }
        for (int s = 0; s < NumScalars; ++s) { packed[index++] = m_initial_state[s]; }
        for (int s = 0; s < NumScalars; ++s) { packed[index++] = current[s]; }
        for (int s = 0; s < NumScalars; ++s) { packed[index++] = m_internal_source[s]; }
        amrex::ParallelDescriptor::ReduceRealSum(packed.data(), static_cast<int>(packed.size()));

        std::array<std::array<amrex::Real, NumFaces>, NumScalars> global_faces{};
        index = 0;
        for (int s = 0; s < NumScalars; ++s) {
            for (int f = 0; f < NumFaces; ++f) { global_faces[s][f] = packed[index++]; }
        }
        std::array<amrex::Real, NumScalars> initial{};
        std::array<amrex::Real, NumScalars> global_current{};
        index = NumScalars*NumFaces;
        for (int s = 0; s < NumScalars; ++s) { initial[s] = packed[index++]; }
        for (int s = 0; s < NumScalars; ++s) { global_current[s] = packed[index++]; }
        std::array<amrex::Real, NumScalars> global_internal{};
        for (int s = 0; s < NumScalars; ++s) { global_internal[s] = packed[index++]; }

        if (amrex::ParallelDescriptor::IOProcessor()) {
            std::ofstream out("cloud_chamber_budget.dat", std::ios::app);
            if (out.tellp() == std::streampos(0)) {
                out << "# start_step start_time end_step end_time scalar units "
                       "xlo xhi ylo yhi zlo zhi net_boundary volume_change "
                       "internal_source residual tolerance status\n";
            }
            out << std::setprecision(17);
            const char* names[NumScalars] = {"rhoTheta", "water_vapor", "cloud_water"};
            const char* units[NumScalars] = {"conserved_scalar", "mixing_ratio_density", "mixing_ratio_density"};
            for (int s = 0; s < NumScalars; ++s) {
                std::array<amrex::Real, NumFaces> faces = global_faces[s];
                amrex::Real net = faces[0] - faces[1] + faces[2] - faces[3] + faces[4] - faces[5];
                const amrex::Real change = global_current[s] - initial[s];
                const amrex::Real internal = global_internal[s];
            const amrex::Real residual = change - net - internal;
            const amrex::Real scale = std::max({amrex::Real(1.0),
                    std::abs(initial[s]), std::abs(global_current[s]),
                    std::abs(net), std::abs(internal)});
                const amrex::Real effective_operations =
                    static_cast<amrex::Real>(geom.Domain().numPts()) +
                    static_cast<amrex::Real>(NumFaces);
                amrex::Real boundary_scale = amrex::Real(0.0);
                for (const auto value : faces) { boundary_scale += std::abs(value); }
                const amrex::Real tolerance = amrex::Real(128.0) *
                    std::numeric_limits<amrex::Real>::epsilon() *
                    (effective_operations * scale + boundary_scale + std::abs(internal));
                const bool status = std::abs(residual) <= tolerance;
                const char* status_text = status ? "PASS" :
                    (cloudy && s == RhoTheta ? "UNSUPPORTED_SOURCE" : "FAIL");
                out << m_interval_start_step << ' ' << m_interval_start_time << ' '
                    << step << ' ' << time << ' ' << names[s] << ' ' << units[s];
                for (const auto value : faces) { out << ' ' << value; }
                out << ' ' << net << ' ' << change << ' ' << internal << ' '
                    << residual << ' ' << tolerance << ' ' << status_text << '\n';
            }
            std::array<amrex::Real, NumFaces> total_faces{};
            for (int f = 0; f < NumFaces; ++f) {
                total_faces[f] = global_faces[RhoQv][f] + global_faces[RhoQc][f];
            }
            const amrex::Real total_net = total_faces[0] - total_faces[1] +
                total_faces[2] - total_faces[3] + total_faces[4] - total_faces[5];
            const amrex::Real total_change = global_current[RhoQv] + global_current[RhoQc] -
                initial[RhoQv] - initial[RhoQc];
            const amrex::Real total_internal = global_internal[RhoQv] + global_internal[RhoQc];
            const amrex::Real total_residual = total_change - total_net - total_internal;
            const amrex::Real total_scale = std::max({amrex::Real(1.0), std::abs(total_change),
                std::abs(total_net), std::abs(total_internal)});
            amrex::Real total_boundary_scale = amrex::Real(0.0);
            for (const auto value : total_faces) { total_boundary_scale += std::abs(value); }
            const amrex::Real total_tolerance = amrex::Real(128.0) *
                std::numeric_limits<amrex::Real>::epsilon() *
                (static_cast<amrex::Real>(geom.Domain().numPts() + NumFaces) * total_scale +
                 total_boundary_scale + std::abs(total_internal));
            const bool total_status = std::abs(total_residual) <= total_tolerance;
            out << m_interval_start_step << ' ' << m_interval_start_time << ' '
                << step << ' ' << time << " total_nonprecipitating_water mixing_ratio_density";
            for (const auto value : total_faces) { out << ' ' << value; }
            out << ' ' << total_net << ' ' << total_change << ' ' << total_internal << ' '
                << total_residual << ' ' << total_tolerance << ' '
                << (total_status ? "PASS" : "FAIL") << '\n';
        }
        // The next report is an interval-local balance beginning at the
        // state just reported.  Keep this reference rank-local: assigning the
        // MPI-global value here would make every rank contribute the same
        // global reference at the next reduction.
        m_initial_state = current;
        m_cumulative = {};
        m_internal_source = {};
        m_stage0 = {};
        m_have_stage0 = {};
        m_interval_start_step = step;
        m_interval_start_time = time;
    }

private:
    int m_interval = 0;
    int m_interval_start_step = 0;
    double m_interval_start_time = 0.0;
    bool m_have_initial = false;
    std::array<bool, NumScalars> m_have_stage0{};
    std::array<amrex::Real, NumScalars> m_initial_state{};
    std::array<amrex::Real, NumScalars> m_internal_source{};
    std::array<std::array<amrex::Real, NumFaces>, NumScalars> m_stage0{};
    std::array<std::array<amrex::Real, NumFaces>, NumScalars> m_cumulative{};
};

#endif
