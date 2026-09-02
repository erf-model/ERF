/**
 * @file ERF_DustCriticalMaterials.H
 * @brief DOE Critical Materials Assessment output for mine dust emissions.
 *
 * Computes critical material (CM) emission flux from dust emission flux
 * and user-specified per-bin mass fractions:
 *   F_cm(i,j) = sum_b( F_bin(i,j,b) * cm_fraction[b] )  [kg_CM/m^2/s]
 *
 * Provides domain and per-site mass budget CSV output.
 *
 * Pattern references:
 *   ERF_DustPM.H (Phase 17): per-bin classification.
 *   ERF_FireDiagnostics.H:
 *     https://github.com/hgopalan/ERF/blob/f264ce7dd0f5e727c71c939659a6a5caf001d21b/
 *     Source/Fire/ERF_FireDiagnostics.H
 *   ERF_FireStatsOutput.H (CSV append):
 *     https://github.com/hgopalan/ERF/blob/f264ce7dd0f5e727c71c939659a6a5caf001d21b/
 *     Source/Fire/ERF_FireStatsOutput.H
 *
 * References:
 *   US DOE (2023). Critical Materials Assessment.
 *     https://www.energy.gov/eere/vehicles/articles/us-doe-critical-materials-assessment
 *   Hofstra et al. (2020). Resources, Conservation and Recycling.
 *     https://doi.org/10.1016/j.resconrec.2020.104802
 *   USGS Mineral Resources Program:
 *     https://www.usgs.gov/programs/mineral-resources-program
 *   Nriagu & Pacyna (1988). Nature, 333, 134.
 *     https://doi.org/10.1038/333134a0
 */
#ifdef ERF_USE_DUST
#ifndef ERF_DUST_CRITICAL_MATERIALS_H_
#define ERF_DUST_CRITICAL_MATERIALS_H_

#include <AMReX_MultiFab.H>
#include <AMReX_MFIter.H>
#include <AMReX_Geometry.H>
#include <AMReX_ParallelDescriptor.H>
#include <AMReX_Vector.H>
#include <fstream>
#include <string>
#include <ERF_HostFabView.H>

/**
 * @brief Compute critical material emission flux from dust emission flux.
 *
 * F_cm(i,j) = sum_b( F_dust(i,j,b) * cm_fraction[b] )
 *
 * When transport_bins_separately=false, emission_flux has 1 component
 * and cm_fraction[0] applies. When true, one component per active bin.
 *
 * @param[out] cm_flux       Critical material flux [kg_CM/m^2/s] on dust grid.
 * @param[in]  emission_flux dust_emission_flux MultiFab (1 or n_bins components).
 * @param[in]  cm_fractions  Per-bin CM mass fractions [kg_CM/kg_dust].
 * @param[in]  n_active      Number of active bins (1 or n_size_bins).
 */
inline void compute_cm_flux(
    amrex::MultiFab&       cm_flux,
    const amrex::MultiFab& emission_flux,
    const amrex::Vector<amrex::Real>& cm_fractions,
    int n_active)
{
    cm_flux.setVal(0.0);
    for (int b = 0; b < n_active; ++b) {
        int fi = (b < (int)cm_fractions.size())
               ? b : (int)cm_fractions.size()-1;
        amrex::Real frac = cm_fractions[fi];
        if (frac <= 0.0) continue;
        for (amrex::MFIter mfi(cm_flux, amrex::TilingIfNotGPU());
             mfi.isValid(); ++mfi) {
            const amrex::Box& bx = mfi.tilebox();
            auto cm  = cm_flux.array(mfi);
            auto fl  = emission_flux.const_array(mfi);
            const int src_comp = b;
            amrex::ParallelFor(bx,
                [=] AMREX_GPU_DEVICE(int i, int j, int k) noexcept {
                    cm(i,j,k) += amrex::max(fl(i,j,k,src_comp), 0.0) * frac;
            });
        }
    }
}

/**
 * @brief Append domain and per-site CM mass budget to CSV.
 *
 * Columns: step, time_s, site_index, site_name, cm_emission_kg_s.
 * One row per site per call (plus one row for domain total, site_index=-1).
 * Pattern: append_fire_stats in ERF_FireStatsOutput.H.
 *
 * @param filename      CSV file path.
 * @param cm_flux       Critical material flux [kg_CM/m^2/s] on dust grid.
 * @param site_id       dust_site_id from Phase 20 (may be null).
 * @param geom_dust     Dust grid geometry.
 * @param cur_time      Current time [s].
 * @param nstep         Current step.
 * @param site_names    Phase 20 site names (may be empty).
 */
inline void append_cm_budget(
    const std::string&                filename,
    const amrex::MultiFab&            cm_flux,
    const amrex::MultiFab*            site_id,
    const amrex::Geometry&            geom_dust,
    amrex::Real                       cur_time,
    int                               nstep,
    const amrex::Vector<std::string>& site_names)
{
    const auto dx = geom_dust.CellSizeArray();
    const amrex::Real cell_area = dx[0] * dx[1];
    int n_sites = (int)site_names.size();

    amrex::Vector<amrex::Real> site_cm(n_sites + 1, 0.0);

    // Scatter-accumulation into named per-site bins, so it stays a host loop --
    // but the MultiFabs are device-only under CUDA, hence the staged views.
    for (amrex::MFIter mfi(cm_flux); mfi.isValid(); ++mfi) {
        const amrex::Box& bx = mfi.validbox();
        const ERFHostFabView cm_view(cm_flux[mfi]);
        auto cm_arr = cm_view.array();
        if (site_id) {
            const ERFHostFabView sid_view((*site_id)[mfi]);
            auto sid = sid_view.array();
            amrex::LoopOnCpu(bx, [&](int i, int j, int k) noexcept {
                int s = amrex::max(0, amrex::min(n_sites, (int)sid(i,j,k)));
                site_cm[s] += cm_arr(i,j,k) * cell_area;
            });
        } else {
            amrex::LoopOnCpu(bx, [&](int i, int j, int k) noexcept {
                site_cm[0] += cm_arr(i,j,k) * cell_area;
            });
        }
    }
    for (int s = 0; s <= n_sites; ++s)
        amrex::ParallelDescriptor::ReduceRealSum(site_cm[s]);

    amrex::Real domain_total = 0.0;
    for (int s = 0; s <= n_sites; ++s) domain_total += site_cm[s];

    if (!amrex::ParallelDescriptor::IOProcessor()) return;

    bool file_exists = std::ifstream(filename).good();
    std::ofstream out(filename, std::ios::app);
    if (!out.is_open()) {
        amrex::Print() << "[DUST] WARNING: Cannot open CM budget file: "
                       << filename << "\n";
        return;
    }
    if (!file_exists) {
        out << "# DOE Critical Materials mass budget\n"
            << "# DOE (2023): https://www.energy.gov/eere/vehicles/articles/"
            << "us-doe-critical-materials-assessment\n"
            << "step,time_s,site_index,site_name,cm_emission_kg_s\n";
    }
    // Domain total (site_index = -1).
    out << nstep << "," << cur_time << ",-1,domain," << domain_total << "\n";
    for (int s = 0; s <= n_sites; ++s) {
        std::string sname = (s == 0) ? "unassigned"
                          : site_names[s-1];
        out << nstep << "," << cur_time << "," << s
            << "," << sname << "," << site_cm[s] << "\n";
    }
}
#endif
#endif
