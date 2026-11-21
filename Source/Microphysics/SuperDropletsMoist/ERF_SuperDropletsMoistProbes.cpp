#include "ERF_SuperDropletsMoist.H"

#ifdef ERF_USE_PARTICLES
#include "ERF_TerrainConversion.H"
#include <AMReX_Algorithm.H>
#include <AMReX_ParmParse.H>
#include <algorithm>
#include <cmath>
#include <cstdint>
#include <filesystem>
#include <iomanip>
#include <set>
#include <sstream>

using namespace amrex;

namespace {
// Reconcile the series once, at the checkpoint time (even on an unsampled step).
// Times are written with round-trip precision. Never append different bins or volumes.
void open_series (std::ofstream& out, const std::string& path, const std::string& header,
                  bool restart, double time)
{
    const bool append = restart && std::filesystem::exists(path);
    if (append) {
        std::ifstream in(path, std::ios::binary);
        std::istringstream expected(header);
        std::string line, want;
        while (std::getline(expected, want)) {
            if (!std::getline(in, line) || line != want) {
                Abort("SDM DSD: incompatible probe header in " + path);
            }
        }
        for (;;) {
            const auto offset = in.tellg();
            if (!std::getline(in, line)) {
                if (!in.eof()) { Abort("SDM DSD: error reading " + path); }
                break;
            }
            // A last row without a newline was cut short (e.g. by a crash): drop it.
            if (in.eof()) {
                in.close();
                std::filesystem::resize_file(path, static_cast<std::uintmax_t>(offset));
                break;
            }
            double t;
            if (!(std::istringstream(line) >> t)) { Abort("SDM DSD: invalid time in " + path); }
            if (t >= time) {
                in.close();
                std::filesystem::resize_file(path, static_cast<std::uintmax_t>(offset));
                break;
            }
        }
    }
    // Keep line endings and byte offsets consistent with the binary restart scan.
    out.open(path, std::ios::binary | (append ? std::ios::app : std::ios::trunc));
    if (!out) { FileOpenFailed(path); }
    if (!append) { out << header; }
    out << std::setprecision(std::numeric_limits<double>::max_digits10);
}
}

void SuperDropletsMoist::readDSDProbes ()
{
    ParmParse pp(m_name);
    Vector<std::string> names;
    if (!pp.queryarr("dsd_probe_names", names)) { return; }
    pp.query("dsd_probe_interval", m_dsd_probe_interval);
    if (m_dsd_probe_interval <= 0) { Abort(m_name + ".dsd_probe_interval must be positive"); }
    std::set<std::string> seen;
    for (const auto& name : names) {
        if (name.empty() || name.find_first_not_of("abcdefghijklmnopqrstuvwxyzABCDEFGHIJKLMNOPQRSTUVWXYZ0123456789_-")
                != std::string::npos || !seen.insert(name).second) {
            Abort("SDM DSD: invalid or duplicate probe name: " + name);
        }
        auto& p = m_dsd_probes.emplace_back();
        p.name = name;
        const auto prefix = m_name + "." + name;
        ParmParse cfg(prefix);
        auto fail = [&](const std::string& msg) { Abort(prefix + ": " + msg); };
        p.by_cell = cfg.contains("cell");
        const bool corners = cfg.contains("box_lo") || cfg.contains("box_hi");
        const bool centered = cfg.contains("center") || cfg.contains("side_lengths");
        if (int(p.by_cell) + int(corners) + int(centered) != 1) {
            fail("give exactly one of cell, box_lo/box_hi, or center/side_lengths");
        }
        if (p.by_cell) {
            cfg.getarr("cell", p.cell);
            if (!m_geom.Domain().contains(p.cell)) { fail("cell is outside the level-0 grid"); }
        } else {
            Vector<Real> a, b;
            cfg.getarr(corners ? "box_lo" : "center", a);
            cfg.getarr(corners ? "box_hi" : "side_lengths", b);
            if (centered && b.size() == 1) { b.assign(3, Real(b[0])); }
            if (a.size() != 3 || b.size() != 3) { fail("expected three coordinates per corner or center and one or three side lengths"); }
            for (int d = 0; d < 3; ++d) {
                p.lo[d] = corners ? a[d] : a[d] - b[d]/2;
                p.hi[d] = corners ? b[d] : a[d] + b[d]/2;
                if (!std::isfinite(p.lo[d]) || !std::isfinite(p.hi[d]) || p.hi[d] <= p.lo[d]) { fail("invalid box extent"); }
                // On a stretched/terrain mesh ProbHi(2) is not a physical height.
                if ((d < 2 || !m_variable_dz) && (p.lo[d] < m_geom.ProbLo(d) || p.hi[d] > m_geom.ProbHi(d))) {
                    fail("box extends outside the domain");
                }
            }
        }
        constexpr int max_bins = (std::numeric_limits<int>::max()-1)/2;
        if (cfg.contains("r_edges")) {
            if (cfg.contains("rmin") || cfg.contains("rmax") || cfg.contains("nbins") || cfg.contains("spacing")) {
                fail("r_edges cannot be combined with generated-bin parameters");
            }
            cfg.getarr("r_edges", p.edges);
        } else {
            double lo, hi;
            int n;
            std::string spacing = "linear";
            cfg.get("rmin", lo); cfg.get("rmax", hi); cfg.get("nbins", n); cfg.query("spacing", spacing);
            if (n < 1 || n > max_bins || !std::isfinite(lo) || !std::isfinite(hi) || !(0 <= lo && lo < hi)
                || (spacing != "linear" && spacing != "log") || (spacing == "log" && lo == 0)) {
                fail("invalid rmin/rmax/nbins/spacing (log requires rmin > 0)");
            }
            p.edges.resize(n+1);
            p.edges.front() = lo; p.edges.back() = hi;
            for (int i = 1; i < n; ++i) {
                const double f = static_cast<double>(i)/n;
                p.edges[i] = spacing == "log" ? std::exp(std::lerp(std::log(lo), std::log(hi), f)) : std::lerp(lo, hi, f);
            }
        }
        if (p.edges.size() < 2) { fail("expected at least two radius edges"); }
        if (p.edges.size()-1 > max_bins) { fail("too many radius bins"); }
        double prev = -1;
        for (const auto edge : p.edges) {
            if (!std::isfinite(edge) || edge < 0 || edge <= prev) { fail("radius edges must be finite, nonnegative and strictly increasing"); }
            prev = edge;
        }
        cfg.query("write_mass", p.write_mass);
        p.device_edges.resize(p.edges.size());
        Gpu::copy(Gpu::hostToDevice, p.edges.begin(), p.edges.end(), p.device_edges.begin());
        p.sums.resize((p.edges.size()-1)*(p.write_mass ? 2 : 1) + 1); // last slot: volume from one owner
        p.histogram.resize(p.sums.size());
    }
}

void SuperDropletsMoist::sampleDSDProbes (int step, Real time, const Vector<MFPtr>& z_nd)
{
    BL_PROFILE("SuperDropletsMoist::sampleDSDProbes()");
    const bool due = step % m_dsd_probe_interval == 0;
    const auto& pc = *m_super_droplets;
    const int water = ridx_s(m_idx_w, m_num_aerosols, m_num_species);
    const double rho_w = pc.getSpeciesMaterial(Species::Name::H2O).m_density;
    const auto& g0 = pc.Geom(0);
    const auto io = ParallelDescriptor::IOProcessorNumber();
    for (auto& probe : m_dsd_probes) {
        const int n = static_cast<int>(probe.edges.size())-1;
        const bool by_cell = probe.by_cell, mass = probe.write_mass;
        const auto lo = probe.lo, hi = probe.hi;
        if (ParallelDescriptor::IOProcessor() && !probe.output.is_open()) {
            std::ostringstream header;
            header << std::setprecision(std::numeric_limits<double>::max_digits10)
                   << "# SDM DSD v1: liquid-water-equivalent radius; bin-integrated concentrations\n";
            if (by_cell) {
                header << "# level0_cell " << probe.cell << " domain " << g0.Domain() << " bounds_m " << g0.ProbDomain();
            } else {
                header << "# box_lo_m " << lo[0] << ' ' << lo[1] << ' ' << lo[2]
                       << " box_hi_m " << hi[0] << ' ' << hi[1] << ' ' << hi[2];
            }
            header << "\n# radius_edges_m:";
            for (double r : probe.edges) { header << ' ' << r; }
            header << "\n# time_s";
            for (int b = 0; b < n; ++b) { header << " N" << b << "[m^-3]"; }
            if (mass) { for (int b = 0; b < n; ++b) { header << " M" << b << "[kg/m^3]"; } }
            header << '\n';
            open_series(probe.output, probe.name + ".dsd", header.str(), m_restarted_particles, time);
        }
        if (!due) { continue; }

        // One small histogram per rank (not per mesh cell), one reduction per probe.
        std::fill(probe.sums.begin(), probe.sums.end(), 0.0);
        Gpu::copy(Gpu::hostToDevice, probe.sums.begin(), probe.sums.end(), probe.histogram.begin());
        auto* hist = probe.histogram.data();
        const auto* edges = probe.device_edges.data();
        Box cells(probe.cell, probe.cell);
        for (int lev = 0; lev <= pc.finestLevel(); ++lev) {
            if (by_cell && lev > 0) { cells.refine(pc.GetParGDB()->refRatio(lev-1)); }
            const auto& geom = pc.Geom(lev);
            const auto plo = geom.ProbLoArray(), dxi = geom.InvCellSizeArray();
            const auto domain = geom.Domain();
            // Skip distant particle tiles before launching a kernel. The bounding cell
            // box is conservative; the per-particle physical test decides faces exactly.
            Box candidates = by_cell ? cells : domain;
            if (!by_cell) {
                for (int d = 0; d < (m_variable_dz ? 2 : 3); ++d) {
                    candidates.setSmall(d, domain.smallEnd(d) + static_cast<int>(std::floor((lo[d]-plo[d])*dxi[d])));
                    candidates.setBig(d, domain.smallEnd(d) + static_cast<int>(std::floor((hi[d]-plo[d])*dxi[d])));
                }
            }
            for (SuperDropletPC::ParConstIterType pti(pc, lev); pti.isValid(); ++pti) {
                if (!pti.tilebox().intersects(candidates)) { continue; }
                const auto ptd = pti.GetParticleTile().getConstParticleTileData();
                const auto z = !by_cell && m_variable_dz ? z_nd[lev]->const_array(pti) : Array4<const Real>{};
                ParallelFor(pti.numParticles(), [=] AMREX_GPU_DEVICE (int i) noexcept {
                    const auto& p = ptd.m_aos[i];
                    if (p.id() <= 0 || ptd.m_runtime_idata[SDPCDefn::SuperDropletsIntIdxSoA_RT::active][i] <= 0) { return; }
                    if (by_cell) {
                        if (!cells.contains(getParticleCell(p, plo, dxi, domain))) { return; }
                    } else {
                        const Real x = p.pos(0), y = p.pos(1);
                        if (x < lo[0] || x >= hi[0] || y < lo[1] || y >= hi[1]) { return; }
                        const Real h = ERF::ParticlePos::z_from_zeta(x, y, p.pos(2), plo, dxi, z);
                        if (h < lo[2] || h >= hi[2]) { return; }
                    }
                    const double mw = ptd.m_runtime_rdata[water][i];
                    if (!(mw > 0)) { return; } // no liquid droplet: ice/dry aerosol alone does not count
                    const double r = std::cbrt(mw/(four_thirds_pi*rho_w));
                    if (!(r >= edges[0] && r < edges[n])) { return; }
                    const int b = static_cast<int>(amrex::upper_bound(edges, edges+n+1, r) - edges)-1;
                    const double xi = ptd.m_runtime_rdata[SDPCDefn::SuperDropletsRealIdxSoA_RT::multiplicity][i];
                    Gpu::Atomic::AddNoRet(hist+b, xi);
                    if (mass) { Gpu::Atomic::AddNoRet(hist+n+b, xi*mw); }
                });
            }
        }
        Gpu::copy(Gpu::deviceToHost, probe.histogram.begin(), probe.histogram.end(), probe.sums.begin());
        // The level-0 owner supplies the physical volume even if there are no SDs in the cell.
        if (by_cell && m_variable_dz) {
            const auto c = probe.cell;
            for (MFIter mfi(*z_nd[0]); mfi.isValid(); ++mfi) {
                if (enclosedCells(mfi.validbox()).contains(c)) {
                    const SDCellVolume vol{g0.CellSizeArray(), z_nd[0]->const_array(mfi), true};
                    probe.sums.back() = Reduce::Sum<double>(1, [=] AMREX_GPU_DEVICE (int) noexcept -> double {
                        return vol(c[0],c[1],c[2]);
                    });
                }
            }
        } else if (ParallelDescriptor::IOProcessor()) {
            probe.sums.back() = by_cell ? g0.CellSize(0)*g0.CellSize(1)*g0.CellSize(2)
                                        : (hi[0]-lo[0])*(hi[1]-lo[1])*(hi[2]-lo[2]);
        }
        ParallelDescriptor::ReduceRealSum(probe.sums.data(), static_cast<int>(probe.sums.size()), io);
        if (ParallelDescriptor::IOProcessor()) {
            const double volume = probe.sums.back();
            if (!std::isfinite(volume) || volume <= 0) { Abort("SDM DSD: invalid sample volume for " + probe.name); }
            probe.output << time;
            for (int b = 0; b < n*(mass ? 2 : 1); ++b) { probe.output << ' ' << probe.sums[b]/volume; }
            probe.output << std::endl;
            if (!probe.output) { Abort("SDM DSD: error writing " + probe.name + ".dsd"); }
        }
    }
}
#endif
