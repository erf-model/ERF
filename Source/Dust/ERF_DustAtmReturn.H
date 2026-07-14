/**
 * @file ERF_DustAtmReturn.H
 * @brief Return fields from the ERF 3D conserved state to the 2D dust grid.
 *
 * Implements the ERF-Atm -> ERF-Dust return path (Phase 13). Two fields
 * are extracted from the 3D solver and mapped to the 2D dust grid:
 *
 *   dust_conc_sfc:   RhoDust at k=klo [kg/m^3] from S_new conserved state.
 *                    Feeds Phase 5 Shao (2001) loading feedback.
 *
 *   dust_surf_moist: Q1fx3 at z-face klo [kg/m^2/s]. Null-safe: zeroed when
 *                    moisture_type == None (SFS_q1fx3_lev is null in that case,
 *                    confirmed in ERF_MakeNewArrays.cpp). Feeds Phase 5 Fecan
 *                    moisture inhibition when moisture is active.
 *
 * Coarsening pattern (same as fill_dust_scalar_from_atm in ERF_DustWindExtract.H):
 *   dust_field(i,j,0) = atm_field(i/C, j/C, klo)
 *
 * Pattern references:
 *   ERF_FireWindExtract.H, Source/Fire/:
 *     https://github.com/hgopalan/ERF/blob/ecb9b7656bbb98fcc5d4f781e486825ee9478a73/
 *     Source/Fire/ERF_FireWindExtract.H
 *   fill_dust_scalar_from_atm in ERF_DustWindExtract.H (Phase 9).
 *
 * References:
 *   Shao, Y. (2001). J. Geophys. Res., 106, 20599.
 *     https://doi.org/10.1029/2001JD900171
 *   Fecan, F., Marticorena, B., & Bergametti, G. (1999). Ann. Geophys., 17, 149.
 *     https://doi.org/10.1007/s00585-999-0149-7
 *   Hong, S.-Y., & Pan, H.-L. (1996). Mon. Wea. Rev., 124, 2322.
 *     https://doi.org/10.1175/1520-0493(1996)124<2322:NBLVDI>2.0.CO;2
 */

#ifndef ERF_DUST_ATM_RETURN_H
#define ERF_DUST_ATM_RETURN_H

#include <AMReX_MultiFab.H>
#include <AMReX_Geometry.H>
#include <AMReX_MFIter.H>

/**
 * @brief Map RhoDust at k=klo from the 3D conserved state to the 2D dust grid.
 *
 * dust_conc_sfc(i,j,0) = max(S_cons(i/C, j/C, klo, dust_comp), 0)
 *
 * @param[out] dust_conc_sfc  2D MultiFab on dust grid [kg/m^3].
 * @param[in]  S_cons         3D conserved state (S_new) on atm grid.
 * @param[in]  dust_comp      Component index for dust scalar.
 * @param[in]  geom_atm       Atmospheric geometry.
 * @param[in]  C              Grid ratio (1 when dust grid == atm 2D grid).
 */
inline void fill_dust_conc_from_atm(
    amrex::MultiFab&       dust_conc_sfc,
    const amrex::MultiFab& S_cons,
    int                    dust_comp,
    const amrex::Geometry& geom_atm,
    int                    C)
{
    const int klo = geom_atm.Domain().smallEnd(2);
    for (amrex::MFIter mfi(dust_conc_sfc, amrex::TilingIfNotGPU());
         mfi.isValid(); ++mfi) {
        const amrex::Box& bx = mfi.tilebox();
        auto conc = dust_conc_sfc.array(mfi);
        auto cons = S_cons.const_array(mfi);
        const int comp = dust_comp;
        const int kl   = klo;
        amrex::ParallelFor(bx, [=] AMREX_GPU_DEVICE(int i, int j, int k) noexcept {
            conc(i,j,k) = amrex::max(cons(i/C, j/C, kl, comp), 0.0);
        });
    }
}

/**
 * @brief Map vertical moisture flux at z-face klo from atm grid to dust grid.
 *
 * dust_surf_moist(i,j,0) = Q1fx3(i/C, j/C, klo)   [kg/m^2/s]
 *
 * Q1fx3 is face-staggered (z-faces, IntVect(0,0,1)). The surface flux is
 * at face index klo (the bottom face of the lowest cell).
 *
 * Null-safe: when Q1fx3 == nullptr (moisture_type == None, which is the
 * case for all current dust tests), dust_surf_moist is zeroed and the
 * function returns immediately. The Phase 5 threshold computation then
 * uses the static moisture flag from the Phase 3 surface map unchanged.
 * The null path costs one pointer comparison per call.
 *
 * When Q1fx3 is non-null (moisture_type != None), the full extraction
 * runs and dust_surf_moist carries the surface latent flux for use in
 * the Fecan et al. (1999) moisture inhibition term.
 *
 * @param[out] dust_surf_moist  2D MultiFab on dust grid [kg/m^2/s].
 * @param[in]  Q1fx3            Vertical moisture flux (face-staggered z).
 *                              SFS_q1fx3_lev[level].get() from ERF.H.
 *                              Null when moisture_type == None.
 * @param[in]  geom_atm         Atmospheric geometry.
 * @param[in]  C                Grid ratio.
 */
inline void fill_dust_moist_from_atm(
    amrex::MultiFab&       dust_surf_moist,
    const amrex::MultiFab* Q1fx3,
    const amrex::Geometry& geom_atm,
    int                    C)
{
    if (!Q1fx3) {
        // Moisture scheme not active. SFS_q1fx3_lev is null when
        // l_use_moist = false in ERF_MakeNewArrays.cpp. Static Phase 3
        // moisture flag is used by Phase 5 threshold in this case.
        dust_surf_moist.setVal(0.0);
        return;
    }
    const int klo = geom_atm.Domain().smallEnd(2);
    for (amrex::MFIter mfi(dust_surf_moist, amrex::TilingIfNotGPU());
         mfi.isValid(); ++mfi) {
        const amrex::Box& bx = mfi.tilebox();
        auto moist = dust_surf_moist.array(mfi);
        auto qfx   = Q1fx3->const_array(mfi);
        const int kl = klo;
        amrex::ParallelFor(bx, [=] AMREX_GPU_DEVICE(int i, int j, int k) noexcept {
            // Q1fx3 is at z-faces; face klo is the domain bottom surface.
            moist(i,j,k) = qfx(i/C, j/C, kl);
        });
    }
}

#endif // ERF_DUST_ATM_RETURN_H
