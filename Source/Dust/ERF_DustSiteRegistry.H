/**
 * @file ERF_DustSiteRegistry.H
 * @brief Multi-mine site registry: per-site bounding boxes and PHREEQC tables.
 *
 * Assigns each dust grid cell to a named mine site based on its (x,y) position
 * relative to per-site bounding boxes. Cells outside all boxes receive site_id=0
 * and use the global Phase 4 PHREEQC table.
 *
 * Site assignment at cell (i,j):
 *   x_c = lo[0] + (i + 0.5) * dx[0]
 *   y_c = lo[1] + (j + 0.5) * dx[1]
 *   site_id(i,j) = last s in [1..N] where x_c in [x_lo_s, x_hi_s]
 *                                       and y_c in [y_lo_s, y_hi_s]
 *   site_id(i,j) = 0 if no site matches.
 *
 * Spatial pattern reference: ERF_FireGrid.H (Source/Fire/):
 *   https://github.com/hgopalan/ERF/blob/ERF-Fire/Source/Fire/ERF_FireGrid.H
 *
 * PHREEQC reference:
 *   Parkhurst & Appelo (2013). PHREEQC Version 3.
 *     https://pubs.usgs.gov/tm/06/a43/
 *
 * Threshold reduction references:
 *   Marticorena & Bergametti (1995). J. Geophys. Res., 100, 16415.
 *     https://doi.org/10.1029/95JD00690
 *   Shao & Lu (2000). J. Geophys. Res., 105, 22437.
 *     https://doi.org/10.1029/2000JD900304
 */
#ifdef ERF_USE_DUST
#ifndef ERF_DUST_SITE_REGISTRY_H_
#define ERF_DUST_SITE_REGISTRY_H_

#include <AMReX_MultiFab.H>
#include <AMReX_Geometry.H>
#include <AMReX_MFIter.H>
#include <AMReX_Vector.H>
#include <string>

/**
 * @brief Populate dust_site_id with per-cell site indices.
 *
 * Iterates over all dust grid cells. For each cell, checks site bounding
 * boxes in order 1..N and assigns the last matching site index.
 * Cells not covered by any site receive 0.
 *
 * @param[out] site_id      2D MultiFab on dust grid [int-valued Real].
 * @param[in]  geom_dust    Dust grid geometry.
 * @param[in]  x_lo, y_lo   Per-site bounding box lower corners [m].
 * @param[in]  x_hi, y_hi   Per-site bounding box upper corners [m].
 */
inline void populate_dust_site_id(
    amrex::MultiFab& site_id,
    const amrex::Geometry& geom_dust,
    const amrex::Vector<amrex::Real>& x_lo,
    const amrex::Vector<amrex::Real>& y_lo,
    const amrex::Vector<amrex::Real>& x_hi,
    const amrex::Vector<amrex::Real>& y_hi)
{
    site_id.setVal(0.0);
    const int n_sites = (int)x_lo.size();
    if (n_sites == 0) return;

    const auto& dx  = geom_dust.CellSize();
    const auto& lo  = geom_dust.ProbLo();

    for (amrex::MFIter mfi(site_id, amrex::TilingIfNotGPU());
         mfi.isValid(); ++mfi) {
        const amrex::Box& bx = mfi.tilebox();
        auto sid = site_id.array(mfi);
        // CPU loop — site assignment is a one-time initialisation.
        amrex::LoopOnCpu(bx, [&](int i, int j, int k) noexcept {
            amrex::Real xc = lo[0] + (i + 0.5) * dx[0];
            amrex::Real yc = lo[1] + (j + 0.5) * dx[1];
            for (int s = 0; s < n_sites; ++s) {
                if (xc >= x_lo[s] && xc <= x_hi[s] &&
                    yc >= y_lo[s] && yc <= y_hi[s]) {
                    sid(i,j,k) = amrex::Real(s + 1); // 1-based index
                }
            }
        });
    }
}

/**
 * @brief Count dust grid cells assigned to each site.
 *
 * Returns a vector of length n_sites+1 where index 0 is unassigned count
 * and indices 1..n_sites are per-site counts.
 *
 * @param site_id  2D MultiFab with site indices.
 * @param n_sites  Number of registered sites.
 */
inline amrex::Vector<long> count_site_cells(
    const amrex::MultiFab& site_id, int n_sites)
{
    amrex::Vector<long> counts(n_sites + 1, 0L);
    for (amrex::MFIter mfi(site_id); mfi.isValid(); ++mfi) {
        const amrex::Box& bx = mfi.validbox();
        auto sid = site_id.const_array(mfi);
        amrex::LoopOnCpu(bx, [&](int i, int j, int k) noexcept {
            int s = (int)sid(i,j,k);
            if (s >= 0 && s <= n_sites) ++counts[s];
        });
    }
    for (int s = 0; s <= n_sites; ++s)
        amrex::ParallelDescriptor::ReduceLongSum(counts[s]);
    return counts;
}
#endif
#endif
