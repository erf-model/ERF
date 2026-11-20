#ifndef ERF_RUNNING_STATS_H
#define ERF_RUNNING_STATS_H

/** Running mean, rms, minimum and maximum of a body's scalar diagnostics, carried across checkpoints. */

#include <cmath>
#include <fstream>
#include <iomanip>
#include <limits>
#include <map>
#include <sstream>
#include <string>
#include <vector>

#include <AMReX.H>
#include <AMReX_BLassert.H>
#include <AMReX_ParallelDescriptor.H>
#include <AMReX_Print.H>
#include <AMReX_REAL.H>
#include <AMReX_Utility.H>
#include <AMReX_Vector.H>

#include "ERF_DiagnosticsLog.H"

namespace erf_actuator {

/**
 * One set of running statistics: the mean, the root mean square, the minimum and the maximum of
 * each named scalar quantity over every sample added since the caller's averaging start.
 * write() rewrites <output_root>_stats.csv with the current
 * values; the sums are part of the checkpoint, so the statistics continue across a restart.
 * A non-finite sample aborts, naming the statistics and the quantity.
 */
class RunningStats
{
public:
    RunningStats (std::string name, std::string output_root, std::vector<std::string> quantities)
        : m_name(std::move(name)), m_output_root(std::move(output_root)), m_names(std::move(quantities)),
          m_sum(m_names.size(), 0.0), m_sum2(m_names.size(), 0.0),
          m_min(m_names.size(), std::numeric_limits<amrex::Real>::max()),
          m_max(m_names.size(), std::numeric_limits<amrex::Real>::lowest()) {}

    const std::string& name () const { return m_name; }
    const std::vector<std::string>& quantities () const { return m_names; }

    /**
     * Accept a checkpoint listing quantity `old_name` where `now` stands: an earlier name of the same
     * quantity (`same` true: its statistics continue), or of another quantity (`same` false: the
     * statistics start afresh at the restart, and say so).
     */
    void renamed_from (const std::string& old_name, const std::string& now, bool same)
    {
        m_renamed[old_name] = std::make_pair(now, same);
    }
    int num_samples () const { return m_count; }
    double t_first () const { return m_t_first; }
    double t_last () const { return m_t_last; }

    /** Add one sample of every quantity, taken at time t (s). */
    void accumulate (double t, const std::vector<amrex::Real>& values)
    {
        AMREX_ALWAYS_ASSERT_WITH_MESSAGE(values.size() == m_names.size(), "RunningStats::accumulate: one value per quantity");
        for (std::size_t i = 0; i < values.size(); ++i) {
            if (!std::isfinite(values[i])) {
                std::ostringstream msg;
                msg << "RunningStats " << m_name << ": quantity " << m_names[i] << " is not finite at t = " << t << " s";
                amrex::Abort(msg.str());
            }
        }
        if (m_count == 0) { m_t_first = t; }
        m_t_last = t;
        for (std::size_t i = 0; i < values.size(); ++i) {
            const double v = static_cast<double>(values[i]);
            m_sum[i] += v;
            m_sum2[i] += v * v;
            m_min[i] = std::min(m_min[i], values[i]);
            m_max[i] = std::max(m_max[i], values[i]);
        }
        ++m_count;
    }
    amrex::Real mean (std::size_t i) const { return (m_count > 0) ? static_cast<amrex::Real>(m_sum[i] / m_count) : amrex::Real(0.0); }
    amrex::Real rms (std::size_t i) const
    {
        return (m_count > 0) ? static_cast<amrex::Real>(std::sqrt(std::max(m_sum2[i] / m_count, 0.0))) : amrex::Real(0.0);
    }
    amrex::Real min (std::size_t i) const { return (m_count > 0) ? m_min[i] : amrex::Real(0.0); }
    amrex::Real max (std::size_t i) const { return (m_count > 0) ? m_max[i] : amrex::Real(0.0); }

    /** Rewrite <output_root>_stats.csv: one row per quantity with the current statistics. */
    void write () const
    {
        if (!amrex::ParallelDescriptor::IOProcessor()) { return; }
        std::ofstream out;
        open_log(out, m_output_root + "_stats.csv", true);
        out << "samples,t_first,t_last,quantity,mean,rms,min,max\n" << std::setprecision(10);
        for (std::size_t i = 0; i < m_names.size(); ++i) {
            out << m_count << "," << m_t_first << "," << m_t_last << "," << m_names[i] << ","
                << mean(i) << "," << rms(i) << "," << min(i) << "," << max(i) << "\n";
        }
    }

    /** The state for a checkpoint: <dir>/<name>_stats.dat. */
    void write_state (const std::string& dir) const
    {
        if (!amrex::ParallelDescriptor::IOProcessor()) { return; }
        std::ofstream out(dir + "/" + m_name + "_stats.dat", std::ios::trunc);
        if (!out) { amrex::Abort("cannot write the statistics checkpoint for " + m_name + " in '" + dir + "'"); }
        out << std::setprecision(17) << "count = " << m_count << "\nsize = " << m_names.size()
            << "\nt_first = " << m_t_first << "\nt_last = " << m_t_last << "\n";
        for (std::size_t i = 0; i < m_names.size(); ++i) {
            out << m_names[i] << " " << m_sum[i] << " " << m_sum2[i] << " " << m_min[i] << " " << m_max[i] << "\n";
        }
    }
    /**
     * Restore the state written by write_state(); false when there is none. Aborts, changing
     * nothing, when the file is malformed: a negative count, t_last before t_first, a non-finite
     * value, or quantities other than this object's (but those accepted by renamed_from).
     */
    bool read_state (const std::string& dir)
    {
        const std::string fname = dir + "/" + m_name + "_stats.dat";
        if (!file_exists_everywhere(fname)) {
            amrex::Print() << "RunningStats " << m_name << ": no " << fname << " in the checkpoint; the statistics start afresh\n";
            return false;
        }
        amrex::Vector<char> chars;
        amrex::ParallelDescriptor::ReadAndBcastFile(fname, chars);
        std::istringstream in(std::string(chars.dataPtr(), chars.size()));
        std::string key, eq, qname;
        std::size_t size = 0;
        int count = 0;
        double t0 = 0.0, t1 = 0.0;
        if (!(in >> key >> eq >> count) || key != "count" || !(in >> key >> eq >> size) || key != "size" ||
            !(in >> key >> eq >> t0) || key != "t_first" || !(in >> key >> eq >> t1) || key != "t_last") {
            amrex::Abort("malformed statistics checkpoint '" + fname + "'");
        }
        if (size != m_names.size()) {
            amrex::Abort("the statistics checkpoint '" + fname + "' holds " + std::to_string(size) + " quantities but " +
                         m_name + " has " + std::to_string(m_names.size()));
        }
        if (count < 0 || !std::isfinite(t0) || !std::isfinite(t1) || (count > 0 && t1 < t0)) {
            amrex::Abort("malformed statistics checkpoint '" + fname + "': the count must be >= 0 and t_first <= t_last, both finite");
        }
        std::vector<double> sum(size), sum2(size);
        std::vector<amrex::Real> mn(size), mx(size);
        std::string afresh;
        for (std::size_t i = 0; i < size; ++i) {
            if (!(in >> qname >> sum[i] >> sum2[i] >> mn[i] >> mx[i])) {
                amrex::Abort("the statistics checkpoint '" + fname + "' is truncated: " + m_names[i] + " is missing or incomplete");
            }
            if (qname != m_names[i]) {
                const auto r = m_renamed.find(qname);
                if (r == m_renamed.end() || r->second.first != m_names[i]) {
                    amrex::Abort("the statistics checkpoint '" + fname + "' does not list " + m_names[i] + " where expected");
                }
                if (r->second.second) {
                    amrex::Print() << "RunningStats " << m_name << ": the checkpoint's " << qname << " statistics continue as "
                                   << m_names[i] << "\n";
                } else if (afresh.empty()) {
                    afresh = qname;
                }
            }
            if (!std::isfinite(sum[i]) || !std::isfinite(sum2[i]) || !std::isfinite(mn[i]) || !std::isfinite(mx[i])) {
                amrex::Abort("malformed statistics checkpoint '" + fname + "': a value of " + m_names[i] + " is not finite");
            }
        }
        if (!afresh.empty()) {
            amrex::Print() << "RunningStats " << m_name << ": the checkpoint lists " << afresh << ", another quantity than "
                           << m_renamed[afresh].first << " here; all statistics of " << m_name << " start afresh at the restart\n";
            return true;
        }
        m_sum = sum; m_sum2 = sum2; m_min = mn; m_max = mx;
        m_count = count; m_t_first = t0; m_t_last = t1;
        return true;
    }

private:
    std::string m_name, m_output_root;
    std::vector<std::string> m_names;
    std::vector<double> m_sum, m_sum2;          //!< in double: a float sum over many steps would lose the samples
    std::vector<amrex::Real> m_min, m_max;
    int m_count{0};
    double m_t_first{0.0}, m_t_last{0.0};
    std::map<std::string, std::pair<std::string, bool>> m_renamed;   //!< old name -> (name now, same quantity)
};

} // namespace erf_actuator

#endif
