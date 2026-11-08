#ifndef ERF_UTILS_H_
#define ERF_UTILS_H_

#include <utility>

#include "AMReX.H"
#include "AMReX_MultiFab.H"
#include "AMReX_BCRec.H"
#include "ERF_NumericalConstants.H"
#include "ERF_MicrophysicsConstants.H"
#include "ERF_DataStruct.H"
#include "ERF_IndexDefines.H"
#include "ERF_SurfaceLayer.H"
#include "ERF_FillPatcher.H"
#include "ERF_ImmersedWallCell.H"
#include "ERF_ReadBndryPlanes.H"

/** Build one density-weighted, nonnegative WRF moisture target. */
AMREX_GPU_HOST_DEVICE
AMREX_FORCE_INLINE
amrex::Real
wrf_moisture_target (const amrex::Real rho,
                     const amrex::Real q_n,
                     const amrex::Real q_np1,
                     const amrex::Real alpha) noexcept
{
    const amrex::Real q = (one - alpha) * q_n + alpha * q_np1;
    return rho * amrex::max(q, amrex::Real(0.0));
}

/** True when ERF carries all five hydrometeor mass species separately. */
AMREX_GPU_HOST_DEVICE
AMREX_FORCE_INLINE
bool
wrf_bdy_has_separate_hydrometeors (const MoistureComponentIndices& moisture_indices) noexcept
{
    return moisture_indices.qc >= 0 && moisture_indices.qi >= 0 &&
           moisture_indices.qr >= 0 && moisture_indices.qs >= 0 &&
           moisture_indices.qg >= 0;
}

/** Map an active conserved moisture component to its real-boundary field. */
inline int
wrf_bdy_var_for_moisture_component (const int comp,
                                    const MoistureComponentIndices& moisture_indices,
                                    const bool use_wrf_bdy_qc_qi,
                                    const bool separate_hydrometeors = false)
{
    if (comp == moisture_indices.qv) { return RealBdyVars::QV; }
    if (use_wrf_bdy_qc_qi) {
        if (comp == moisture_indices.qc) { return RealBdyVars::QC; }
        if (moisture_indices.qi >= 0 && comp == moisture_indices.qi) {
            return RealBdyVars::QI;
        }
        if (separate_hydrometeors) {
            if (comp == moisture_indices.qr) { return RealBdyHydrometeorVars::QR; }
            if (comp == moisture_indices.qs) { return RealBdyHydrometeorVars::QS; }
            if (comp == moisture_indices.qg) { return RealBdyHydrometeorVars::QG; }
        }
    }
    return -1;
}

/**
 * The real-boundary maps of the conserved state, one entry per conserved
 * component up to NVAR_max: whether the component is read from the WRF/metgrid
 * boundary data (read), from which boundary field (map), and from which
 * boundary-plane slot when planes are read as well (bnd, cons_bc_comp).
 * Density is read with erf.use_wrf_bdy_density, theta and the water vapour
 * always; passive scalars never. Hydrometeor mappings retain the upstream
 * read flags, which do not enable hydrometeor ghost-cell fills.
 * Built from the component layout, so every passive scalar (NSCALARS of them)
 * sits between RhoKE_comp and RhoQ1_comp.
 */
struct RealBdyConsMaps
{
    amrex::Vector<int> read;
    amrex::Vector<int> map;
    amrex::Vector<int> bnd;
};

inline RealBdyConsMaps
real_bdy_cons_maps (const bool read_density,
                    const MoistureComponentIndices& moisture_indices,
                    const bool use_wrf_bdy_qc_qi,
                    const bool separate_hydrometeors)
{
    RealBdyConsMaps m;
    m.read.assign(NVAR_max, 0);
    m.map.resize(NVAR_max);
    m.bnd.resize(NVAR_max);
    for (int comp = 0; comp < NVAR_max; ++comp) {
        m.map[comp] = comp;
        m.bnd[comp] = cons_bc_comp(comp);
    }
    m.read[Rho_comp]      = read_density ? 1 : 0;
    m.map [Rho_comp]      = read_density ? RealBdyVars::R : Rho_comp;
    m.read[RhoTheta_comp] = 1;
    m.map [RhoTheta_comp] = RealBdyVars::T;
    m.read[RhoQ1_comp]    = 1;
    m.map [RhoQ1_comp]    = RealBdyVars::QV;
    for (int comp = 0; comp < NVAR_max; ++comp) {
        const int bdy_var = wrf_bdy_var_for_moisture_component(
            comp, moisture_indices, use_wrf_bdy_qc_qi, separate_hydrometeors);
        if (bdy_var >= 0) {
            m.map[comp]  = bdy_var;
        }
    }
    return m;
}

/*
 * Rebalance state to satisfy both EOS and HSE one column at a time
 */
void
rebalance_columns (amrex::MultiFab& rho,
                   amrex::MultiFab& theta,
                   const amrex::MultiFab& qt,
                   const amrex::MultiFab& qv,
                   const amrex::MultiFab* z_phys,
                   const amrex::Geometry& geom,
                   const bool& maintain_Th,
                   bool use_sfc = false);

/*
 * Create a new BoxArray in which every grid touches the bottom boundary
 */
void ChopGrids2D (amrex::BoxArray& ba, const amrex::Box& domain, int target_size);

/*
 * Create a new BoxArray with exactly the number of boxes as MPI ranks.
 * This is designed to work only on the base level which covers the entire domain.
 */
amrex::BoxArray
ERFPostProcessBaseGrids (const amrex::Box& domain, bool decompose_in_z);

/**
 * Convert conservative state to primitive variables.
 * @param[in] cons_state Conservative state MultiFab.
 * @param[out] S_prim Primitive state MultiFab.
 * @param[in] ng Number of ghost cells.
 */
void
cons_to_prim (const amrex::MultiFab& cons_state, amrex::MultiFab& S_prim, int ng);

/**
 * Fill the ghost cells of the wall distance field, which is only computed on
 * the valid region, so that operators reading from grown boxes do not see the
 * sentinel value used to initialize it.
 * @param[inout] wdist MultiFab holding the wall distance.
 * @param[in] geom Geometry at this level.
 */
void
fill_wall_dist_ghost_cells (amrex::MultiFab& wdist, const amrex::Geometry& geom);

/**
 * Compute total water mixing ratio from conservative state.
 * @param[in] cons_state Conservative state MultiFab.
 * @param[out] qt MultiFab to store total water.
 * @param[in] n_qstate_into_total Number of components to include in total water.
 */
void
make_qt (const amrex::MultiFab& cons_state, amrex::MultiFab& qt, int n_qstate_into_total);

/**
 * Average down a cell-centred state over immersed-forcing levels (erf.if_average_down =
 * fluid_mass). The fine cells are weighted by their fluid mass,
 * (rho S)_c = rho_c * sum (1 - beta_f) V_f (rho S)_f / sum (1 - beta_f) V_f rho_f,
 * with beta the solid fraction, V the fine cell volume (detJ / (m_x m_y), the weight the plain
 * average carries in the state) and rho_c the coarse density: the fluid-volume-weighted mean of
 * the fine density, sum (1 - beta_f) V_f rho_f / sum (1 - beta_f) V_f, when the density is in the
 * component range, else the coarse cell's own. The coarse value is the fluid mean over the whole
 * coarse cell (intensive), not the sum over the fine cells. Coarse cells that are solid on their
 * own level, or with no fluid fine cell under them, keep their value. The plain average mixes the
 * values of solid fine cells into coarse fluid cells.
 * @param[in] S_fine Fine state (unscaled).
 * @param[in,out] S_crse Coarse state (unscaled); the covered cells are written.
 * @param[in] beta_fine Fine solid fraction (terrain_blanking).
 * @param[in] beta_crse Coarse solid fraction.
 * @param[in] vol_fine Fine cell volume weight (any constant factor cancels).
 * @param[in] scomp First component to average.
 * @param[in] ncomp Number of components.
 * @param[in] ratio Refinement ratio.
 */
void
if_average_down (const amrex::MultiFab& S_fine, amrex::MultiFab& S_crse,
                 const amrex::MultiFab& beta_fine, const amrex::MultiFab& beta_crse,
                 const amrex::MultiFab& vol_fine,
                 int scomp, int ncomp, const amrex::IntVect& ratio);

/**
 * RANS wall distance over terrain represented by immersed forcing with the fraction-stress wall
 * law (erf.if_wall_form = fraction_stress): measured vertically from the surface itself,
 * h = z_bot + beta dz of the wall cell (the partial cell, or the fluid cell on a solid cell;
 * ERF_ImmersedWallCell.H), the wall cell at its fluid centroid (1 - beta) dz / 2 above it and
 * every other cell at |z_cc - h|; a column without solid measures from the bottom of the domain.
 * A box that does not see the surface of a column (a fine box inside the solid, or above the
 * surface without reaching the bottom) measures from the coarser level's wall height.
 * @param[out] wdist Wall distance (valid cells; ghost cells are left untouched).
 * @param[in] terrain_blank Cell-centred solid fraction.
 * @param[in] z_phys_nd Node heights.
 * @param[in] z_phys_cc Cell-centre heights.
 * @param[in] domain Index domain of the level.
 * @param[in] small_volfrac Solid fraction below which a cell counts as fluid.
 * @param[out] wall_height Optional: the height measured from in each column, on the level's
 *            BoxArray flattened to the bottom cell (same DistributionMapping).
 * @param[in] wall_height_crse Optional, fine levels: wall_height of the coarser level brought to
 *            this level's flattened BoxArray. Without it every box must hold whole columns.
 */
void
immersed_wall_dist (amrex::MultiFab& wdist,
                    const amrex::MultiFab& terrain_blank,
                    const amrex::MultiFab& z_phys_nd,
                    const amrex::MultiFab& z_phys_cc,
                    const amrex::Box& domain,
                    amrex::Real small_volfrac,
                    amrex::MultiFab* wall_height = nullptr,
                    const amrex::MultiFab* wall_height_crse = nullptr);

/**
 * Wall-face masks of the fraction-stress wall law (erf.if_wall_form = fraction_stress):
 * zero on the xz (yz) edge that is the bottom face of the wall cell of an x-face (y-face)
 * column, the spacing factor on the face above it, one elsewhere. The vertical flux of x (y)
 * momentum is multiplied by them, so no resolved flux crosses the wall face (see
 * ERF_ImmersedWallCell.H).
 * @param[in] terrain_blank Cell-centred solid fraction; needs one more ghost cell than the masks.
 * @param[out] mask13 On convert(ba, IntVect(1,0,1)), valid and ghost cells.
 * @param[out] mask23 On convert(ba, IntVect(0,1,1)), valid and ghost cells.
 * @param[out] mask33 On convert(ba, IntVect(0,0,1)): the same for the cell columns, applied to
 *             the vertical diffusive fluxes of theta and TKE (zero on the bottom face of the wall cell).
 * @param[in] domain Index domain of the level.
 * @param[in] small_volfrac Solid fraction below which a cell counts as fluid.
 * @param[in] spacing erf.if_wall_face_spacing (0 grid, 1 physical, 2 log): the face ABOVE the wall
 *            cell gets ib_wall::face_above_factor instead of one.
 * @param[in] dz Cell height (constant-dz mesh).
 * @param[in] z0 Roughness length of the wall law.
 */
void
make_ib_wall_face_masks (const amrex::MultiFab& terrain_blank,
                         amrex::MultiFab& mask13,
                         amrex::MultiFab& mask23,
                         amrex::MultiFab& mask33,
                         const amrex::Box& domain,
                         amrex::Real small_volfrac,
                         int spacing = 0,
                         amrex::Real dz = amrex::Real(1.0),
                         amrex::Real z0 = amrex::Real(0.1));

/**
 * Where the coarse and fine levels disagree on an immersed surface along the lateral coarse/fine
 * faces of the fine level (erf.if_cf_mismatch_warning): the surface on a coarse cell's mid-plane
 * leaves the coarse cell half solid over fine cells that are all solid and all fluid. For each
 * coarse cell covered by the fine level and next to an uncovered one, the spread (max - min) of
 * the fine solid fractions under it is formed; cells with a spread of at least warn are counted.
 *
 * @param[in] fine_tb    fine solid fraction (terrain_blanking)
 * @param[in] fine_grids fine BoxArray
 * @param[in] crse_grids coarse BoxArray
 * @param[in] crse_dmap  coarse DistributionMapping
 * @param[in] crse_geom  coarse Geometry
 * @param[in] ratio      refinement ratio
 * @param[in] warn       spread at or above which a cell is counted
 * @return the number of counted cells over all ranks, and the largest spread along those faces
 */
std::pair<int, amrex::Real>
if_cf_surface_mismatch (const amrex::MultiFab& fine_tb,
                        const amrex::BoxArray& fine_grids,
                        const amrex::BoxArray& crse_grids,
                        const amrex::DistributionMapping& crse_dmap,
                        const amrex::Geometry& crse_geom,
                        const amrex::IntVect& ratio,
                        amrex::Real warn);


/*
 * Create the Jacobian for the metric transformation when use_terrain is true
 */
void make_J (const amrex::Geometry& geom,
             amrex::MultiFab& z_phys_nd,
             amrex::MultiFab& detJ_cc);

void make_areas (const amrex::Geometry& geom,
                 amrex::MultiFab& z_phys_nd,
                 amrex::MultiFab& ax,
                 amrex::MultiFab& ay,
                 amrex::MultiFab& az);

/*
 * Average z_phys_nd on nodes to cell centers
 */
void make_zcc (const amrex::Geometry& geom,
               amrex::MultiFab& z_phys_nd,
               amrex::MultiFab& z_phys_cc);

/*
 * Convert momentum to velocity by dividing by density averaged onto faces
 */
void MomentumToVelocity (amrex::MultiFab& xvel_out,
                         amrex::MultiFab& yvel_out,
                         amrex::MultiFab& zvel_out,
                         const amrex::MultiFab& cons_in,
                         const amrex::MultiFab& xmom_in,
                         const amrex::MultiFab& ymom_in,
                         const amrex::MultiFab& zmom_in,
                         const amrex::Box& domain,
                         const amrex::Vector<amrex::BCRec>& domain_bcs_type_h,
                         const amrex::MultiFab* c_vfrac = nullptr // optional
                         );

/*
 * Convert velocity to momentum by multiplying by density averaged onto faces
 */
void VelocityToMomentum (const amrex::MultiFab& xvel_in,
                         const amrex::IntVect & xvel_ngrow,
                         const amrex::MultiFab& yvel_in,
                         const amrex::IntVect & yvel_ngrow,
                         const amrex::MultiFab& zvel_in,
                         const amrex::IntVect & zvel_ngrow,
                         const amrex::MultiFab& cons_in,
                         amrex::MultiFab& xmom_out,
                         amrex::MultiFab& ymom_out,
                         amrex::MultiFab& zmom_out,
                         const amrex::Box& domain,
                         const amrex::Vector<amrex::BCRec>& domain_bcs_type_h,
                         const amrex::MultiFab* c_vfrac = nullptr // optional
                         );

/*
 * Convert (den_div u) to (den_mlt u) by multiplying by (den_mlt/den_div)
 */
void
ConvertForProjection (const amrex::MultiFab& den_div, const amrex::MultiFab& den_mlt,
                      amrex::MultiFab& xmom, amrex::MultiFab& ymom, amrex::MultiFab& zmom,
                      const amrex::Box& domain, const amrex::Vector<amrex::BCRec>& domain_bcs_type_h);

/*
 * Compute lateral boundary mass influx/outflux from base-state-density-weighted
 * boundary normal velocities.
 */
void
compute_influx_outflux_bdy (amrex::FArrayBox& bdy_data_xlo, amrex::FArrayBox& bdy_data_xhi,
                            amrex::FArrayBox& bdy_data_ylo, amrex::FArrayBox& bdy_data_yhi,
                            amrex::Array<amrex::MultiFab*, AMREX_SPACEDIM>& area_vec,
                            const amrex::Geometry& geom,
                            amrex::Real& influx, amrex::Real& outflux,
                            const int n);

/*
 * \brief Enforces solvability on lateral boundary data by scaling outflow to
 * match with inflow.
 */
void enforceInOutSolvability_bdy (const amrex::MultiFab& rho0,
                                  amrex::FArrayBox& bdy_data_xlo,
                                  amrex::FArrayBox& bdy_data_xhi,
                                  amrex::FArrayBox& bdy_data_ylo,
                                  amrex::FArrayBox& bdy_data_yhi,
                                  amrex::Array<amrex::MultiFab*, AMREX_SPACEDIM>& area_vec,
                                  const amrex::Geometry& geom,
                                  const amrex::Vector<amrex::BCRec>& domain_bcs_type_h);

/*
 * \brief Enforces solvability by scaling outflow to match with inflow.
 *
 */
void enforceInOutSolvability (int lev,
                              amrex::Array<amrex::MultiFab*, AMREX_SPACEDIM>& vels_vec,
                              amrex::Array<amrex::MultiFab*, AMREX_SPACEDIM>& area_vec,
                              const amrex::Geometry& geom);

/*
 * Compute boxes for looping over interior/exterior ghost cells
 * for use by fillpatch, erf_slow_rhs_pre, and erf_slow_rhs_post
 */
void realbdy_interior_bxs_xy (const amrex::Box& bx,
                              const amrex::Box& domain,
                              const int& width,
                              amrex::Box& bx_xlo,
                              amrex::Box& bx_xhi,
                              amrex::Box& bx_ylo,
                              amrex::Box& bx_yhi,
                              const amrex::IntVect& ng_vect=amrex::IntVect(0,0,0),
                              const bool get_int_ng=false);

/*
 * Compute boxes for looping over set region cells
 * for use by fillpatch, erf_slow_rhs_pre, and erf_slow_rhs_post
 */
void realbdy_bc_bxs_xy (const amrex::Box& bx,
                        const amrex::Box& domain,
                        const int& set_width,
                        amrex::Box& bx_xlo,
                        amrex::Box& bx_xhi,
                        amrex::Box& bx_ylo,
                        amrex::Box& bx_yhi,
                        const amrex::IntVect& ng_vect=amrex::IntVect(0,0,0));

/*
 * The two boundary-data slices that bracket a time, and the weight of the later one
 */
struct BdyTimeWeights {
    int n_lo;          // slice at or before the time
    int n_hi;          // slice after it (n_lo once the data has run out)
    amrex::Real alpha; // weight of slice n_hi; slice n_lo gets (1 - alpha)
};

/**
 * Time interpolation weights for boundary data at fixed intervals.
 *
 * The time must be measured from the first slice, and formed as
 * elapsed + (start_time - start_bdy_time), not as (elapsed + start_time) - start_bdy_time:
 * start_time and start_bdy_time are epoch seconds (~1.7e9), where doubles are spaced
 * ~2.4e-7 s apart, so the second form rounds the elapsed time to that spacing.  The
 * difference of the two epoch times is exact.
 *
 * @param[in] time_since_start time since the first slice
 * @param[in] data_span        time from the first slice to the last one
 * @param[in] dT               interval between slices
 */
AMREX_FORCE_INLINE
BdyTimeWeights
bdy_time_weights (const double time_since_start, const double data_span, const double dT)
{
    BdyTimeWeights w;
    if (time_since_start >= data_span) {
        // Do not over run the last slice
        w.n_lo  = static_cast<int>(data_span / dT);
        w.n_hi  = w.n_lo;
        w.alpha = zero;
    } else {
        w.n_lo  = static_cast<int>(time_since_start / dT);
        w.n_hi  = w.n_lo + 1;
        w.alpha = static_cast<amrex::Real>((time_since_start - w.n_lo * dT) / dT);
    }
    AMREX_ALWAYS_ASSERT(w.alpha >= zero && w.alpha <= one);
    return w;
}

/*
 * Compute relaxation region RHS with wrfbdy
 */
void realbdy_compute_interior_ghost_rhs (const double& elapsed_time,
                                         const double& start_time,
                                         const double& delta_t,
                                         const double& start_bdy_time,
                                         const double& final_bdy_time,
                                         const double& bdy_time_interval,
                                         const amrex::Real& nudge_factor,
                                         int width,
                                         const amrex::Geometry& geom,
                                         amrex::Vector<amrex::MultiFab>& S_rhs,
                                         amrex::Vector<amrex::MultiFab>& S_cur_data,
                                         amrex::Vector<amrex::Vector<amrex::FArrayBox>>& bdy_data_xlo,
                                         amrex::Vector<amrex::Vector<amrex::FArrayBox>>& bdy_data_xhi,
                                         amrex::Vector<amrex::Vector<amrex::FArrayBox>>& bdy_data_ylo,
                                         amrex::Vector<amrex::Vector<amrex::FArrayBox>>& bdy_data_yhi,
                                         std::unique_ptr<ReadBndryPlanes>& m_r2d,
                                         const amrex::Real& c_p,
                                         const amrex::Real& rdOcp,
                                         const bool use_wrf_bdy_density,
                                         const amrex::Real& bdy_rho_nudge_factor);

/*
 * Compute relaxation region RHS at fine-crse interface
 */
void
fine_compute_interior_ghost_rhs (const double& time,
                                 const double& delta_t,
                                 const int& width,
                                 const int& set_width,
                                 const amrex::Geometry& geom,
                                 ERFFillPatcher* FPr_c,
                                 ERFFillPatcher* FPr_u,
                                 ERFFillPatcher* FPr_v,
                                 ERFFillPatcher* FPr_w,
                                 amrex::Vector<amrex::BCRec>& domain_bcs_type,
                                 amrex::Vector<amrex::MultiFab>& S_rhs_f,
                                 amrex::Vector<amrex::MultiFab>& S_data_f);

/*
 * Accumulate time averaged velocity fields
 */
void
Time_Avg_Vel_atCC (double dt,
                   double& t_avg_cnt,
                   amrex::MultiFab* vel_t_avg,
                   amrex::MultiFab& xvel,
                   amrex::MultiFab& yvel,
                   amrex::MultiFab& zvel);

/**
 * Accumulate time-weighted cell-centered velocity/theta first and second moments.
 * @param[in] dt accumulation interval
 * @param[in,out] t_mean_cnt accumulated time
 * @param[in,out] interval_means ten-component moment storage
 * @param[in] xvel x-face velocity
 * @param[in] yvel y-face velocity
 * @param[in] zvel z-face velocity
 * @param[in] cons conserved cell-centered state
 */
void
Accumulate_Interval_Means (double dt,
                           double& t_mean_cnt,
                           amrex::MultiFab* interval_means,
                           amrex::MultiFab& xvel,
                           amrex::MultiFab& yvel,
                           amrex::MultiFab& zvel,
                           amrex::MultiFab& cons);

/**
 * Compute the nudging relaxation
 *
 * @param[in] delta_t time step
 * @param[in] icomp component offset
 * @param[in] num_var number of variables to loop
 * @param[in] width width of wrf bdy file
 * @param[in] dom_lo low bound of domain
 * @param[in] dom_hi high bound of domain
 * @param[in] F1 drift relaxation parameter
 * @param[in] bx_xlo box for low x relaxation
 * @param[in] bx_xhi box for high x relaxation
 * @param[in] bx_ylo box for low y relaxation
 * @param[in] bx_yhi box for high y relaxation
 * @param[in] arr_xlo array for low x relaxation
 * @param[in] arr_xhi array for high x relaxation
 * @param[in] arr_ylo array for low y relaxation
 * @param[in] arr_yhi array for high y relaxation
 * @param[in] data_arr data array
 * @param[out] rhs_arr RHS array
 */
AMREX_GPU_HOST
AMREX_FORCE_INLINE
void
realbdy_compute_relaxation (const int& icomp,
                            const int& num_var,
                            const int& width,
                            const amrex::GpuArray<amrex::Real,AMREX_SPACEDIM>& dx,
                            const amrex::GpuArray<amrex::Real,AMREX_SPACEDIM>& ProbLo,
                            const amrex::GpuArray<amrex::Real,AMREX_SPACEDIM>& ProbHi,
                            const amrex::Real& F1,
                            const amrex::Box& bx_xlo,
                            const amrex::Box& bx_xhi,
                            const amrex::Box& bx_ylo,
                            const amrex::Box& bx_yhi,
                            const amrex::Array4<const amrex::Real>& arr_xlo,
                            const amrex::Array4<const amrex::Real>& arr_xhi,
                            const amrex::Array4<const amrex::Real>& arr_ylo,
                            const amrex::Array4<const amrex::Real>& arr_yhi,
                            const amrex::Array4<const amrex::Real>& data_arr,
                            const amrex::Array4<amrex::Real>& rhs_arr,
                            const amrex::Real& c_p,
                            const amrex::Real& rdOcp,
                            const int bdy_moist_nudge_type=0)
{
    amrex::IntVect iv = bx_xlo.type();
    amrex::Real ioff  = (iv[0]==1) ? zero : myhalf;
    amrex::Real joff  = (iv[1]==1) ? zero : myhalf;

    // These are defined in ERF_Constants.H
    amrex::Real l_rdOcp  = rdOcp;
    amrex::Real cond_fac = lcond / c_p; // condensation
    amrex::Real  sub_fac = lsub  / c_p;  // sublimation

    int nq = num_var;

    amrex::ParallelFor( bx_xlo, bx_xhi,
    [=] AMREX_GPU_DEVICE (int i, int j, int k) noexcept
    {
        // Corners with x boxes
        amrex::Real x = ProbLo[0] + (i + ioff) * dx[0];
        amrex::Real y = ProbLo[1] + (j + joff) * dx[1];
        amrex::Real x_end  = ProbLo[0] + width * dx[0];
        amrex::Real y_end  = ProbLo[1] + width * dx[1];
        amrex::Real y_strt = ProbHi[1] - width * dx[1];
        amrex::Real xi     = (x_end - x) / (x_end - ProbLo[0]);
        amrex::Real eta_lo = (y < y_end ) ? (y_end  - y) / (y_end  - ProbLo[1]) : zero;
        amrex::Real eta_hi = (y > y_strt) ? (y - y_strt) / (ProbHi[1] - y_strt) : zero;
        amrex::Real eta    = std::max(eta_lo,eta_hi);
        amrex::Real Factor = std::max(xi*xi,eta*eta) * F1;

        if (icomp == RhoQ1_comp)
        {
            amrex::Real exner_pi = getExnergivenRTh(data_arr(i,j,k,RhoTheta_comp),l_rdOcp,
                                                    data_arr(i,j,k,RhoQ1_comp)/data_arr(i,j,k,Rho_comp));

            // qv
            amrex::Real delta = arr_xlo(i,j,k,0) - data_arr(i,j,k,icomp);
            rhs_arr(i,j,k,icomp) += Factor*delta;

            if (bdy_moist_nudge_type==2) {
                // NOTE: we only operate on qv (above), qc and qi here; the
                //       precipitating species (qr/qs/qg) are left untouched.

                // qc and heat source term due to evaporation
                delta = -data_arr(i,j,k,RhoQ2_comp); // This effectively nudges to 0
                rhs_arr(i,j,k,RhoQ2_comp   ) += Factor * delta;
                rhs_arr(i,j,k,RhoTheta_comp) += Factor * delta * cond_fac / exner_pi;

                // qi and heat source term due to sublimation
                if (nq > 3) { // nq > 3 guarantees that RhoQ3 is ice
                    delta = -data_arr(i,j,k,RhoQ3_comp); // This effectively nudges to 0
                    rhs_arr(i,j,k,RhoQ3_comp   ) += Factor * delta;
                    rhs_arr(i,j,k,RhoTheta_comp) += Factor * delta * sub_fac / exner_pi;
                }
            } // moist_nudge_type
        } else {
            for (int n = 0; n < nq; n++) {
                amrex::Real delta = arr_xlo(i,j,k,n) - data_arr(i,j,k,n+icomp);
                rhs_arr(i,j,k,n+icomp) += Factor*delta;
            }
        }
    },
    [=] AMREX_GPU_DEVICE (int i, int j, int k) noexcept
    {
        // Corners with x boxes
        amrex::Real x = ProbLo[0] + (i + ioff) * dx[0];
        amrex::Real y = ProbLo[1] + (j + joff) * dx[1];
        amrex::Real x_strt = ProbHi[0] - width * dx[0];
        amrex::Real y_strt = ProbHi[1] - width * dx[1];
        amrex::Real y_end  = ProbLo[1] + width * dx[1];
        amrex::Real xi     = (x - x_strt) / (ProbHi[0] - x_strt);
        amrex::Real eta_lo = (y < y_end ) ? (y_end  - y) / (y_end  - ProbLo[1]) : zero;
        amrex::Real eta_hi = (y > y_strt) ? (y - y_strt) / (ProbHi[1] - y_strt) : zero;
        amrex::Real eta    = std::max(eta_lo,eta_hi);
        amrex::Real Factor = std::max(xi*xi,eta*eta) * F1;

        if (icomp == RhoQ1_comp)
        {
            amrex::Real exner_pi = getExnergivenRTh(data_arr(i,j,k,RhoTheta_comp),l_rdOcp,
                                                    data_arr(i,j,k,RhoQ1_comp)/data_arr(i,j,k,Rho_comp));

            // qv
            amrex::Real delta = arr_xhi(i,j,k,0) - data_arr(i,j,k,icomp);
            rhs_arr(i,j,k,icomp) += Factor*delta;

            if (bdy_moist_nudge_type==2) {
                // NOTE: we only operate on qv (above), qc and qi here; the
                //       precipitating species (qr/qs/qg) are left untouched.

                // qc and heat source term due to evaporation
                delta = -data_arr(i,j,k,RhoQ2_comp); // This effectively nudges to 0
                rhs_arr(i,j,k,RhoQ2_comp   ) += Factor*delta;
                rhs_arr(i,j,k,RhoTheta_comp) += Factor * delta * cond_fac / exner_pi;

                // qi and heat source term due to sublimation
                if (nq > 3) { // nq > 3 guarantees that RhoQ3 is ice
                    delta = -data_arr(i,j,k,RhoQ3_comp); // This effectively nudges to 0
                    rhs_arr(i,j,k,RhoQ3_comp   ) += Factor * delta;
                    rhs_arr(i,j,k,RhoTheta_comp) += Factor * delta * sub_fac / exner_pi;
                }
            } // moist_nudge_type
        } else {
            for (int n = 0; n < nq; n++) {
                amrex::Real delta = arr_xhi(i,j,k,n) - data_arr(i,j,k,n+icomp);
                rhs_arr(i,j,k,n+icomp) += Factor*delta;
            }
        }
    });

    amrex::ParallelFor( bx_ylo, bx_yhi,
    [=] AMREX_GPU_DEVICE (int i, int j, int k) noexcept
    {
        // No corners for y boxes
        amrex::Real y      = ProbLo[1] + (j + joff) * dx[1];
        amrex::Real y_end  = ProbLo[1] + width * dx[1];
        amrex::Real eta    = (y_end - y) / (y_end - ProbLo[1]);
        amrex::Real Factor = eta*eta * F1;

        if (icomp == RhoQ1_comp)
        {
            amrex::Real exner_pi = getExnergivenRTh(data_arr(i,j,k,RhoTheta_comp),l_rdOcp,
                                                    data_arr(i,j,k,RhoQ1_comp)/data_arr(i,j,k,Rho_comp));
            // qv
            amrex::Real delta = arr_ylo(i,j,k,0) - data_arr(i,j,k,icomp);
            rhs_arr(i,j,k,icomp) += Factor*delta;

            if (bdy_moist_nudge_type==2) {
                // NOTE: we only operate on qv (above), qc and qi here; the
                //       precipitating species (qr/qs/qg) are left untouched.

                // qc and heat source term due to evaporation
                delta = -data_arr(i,j,k,RhoQ2_comp); // This effectively nudges to 0
                rhs_arr(i,j,k,RhoQ2_comp   ) += Factor*delta;
                rhs_arr(i,j,k,RhoTheta_comp) += Factor * delta * cond_fac / exner_pi;

                // qi and heat source term due to sublimation
                if (nq > 3) { // nq > 3 guarantees that RhoQ3 is ice
                    delta = -data_arr(i,j,k,RhoQ3_comp); // This effectively nudges to 0
                    rhs_arr(i,j,k,RhoQ3_comp   ) += Factor * delta;
                    rhs_arr(i,j,k,RhoTheta_comp) += Factor * delta * sub_fac / exner_pi;
                }
            } // moist_nudge_type
        } else {
            for (int n = 0; n < nq; n++) {
                amrex::Real delta = arr_ylo(i,j,k,n) - data_arr(i,j,k,n+icomp);
                rhs_arr(i,j,k,n+icomp) += Factor*delta;
            }
        }
    },
    [=] AMREX_GPU_DEVICE (int i, int j, int k) noexcept
    {
        // No corners for y boxes
        amrex::Real y      = ProbLo[1] + (j + joff) * dx[1];
        amrex::Real y_strt = ProbHi[1] - width * dx[1];
        amrex::Real eta    = (y - y_strt) / (ProbHi[1] - y_strt);
        amrex::Real Factor = eta*eta * F1;

        if (icomp == RhoQ1_comp)
        {
            amrex::Real exner_pi = getExnergivenRTh(data_arr(i,j,k,RhoTheta_comp),l_rdOcp,
                                                    data_arr(i,j,k,RhoQ1_comp)/data_arr(i,j,k,Rho_comp));
            // qv
            amrex::Real delta = arr_yhi(i,j,k,0) - data_arr(i,j,k,icomp);
            rhs_arr(i,j,k,icomp) += Factor*delta;

            if (bdy_moist_nudge_type==2) {
                // NOTE: we only operate on qv (above), qc and qi here; the
                //       precipitating species (qr/qs/qg) are left untouched.

                // qc and heat source term due to evaporation
                delta = -data_arr(i,j,k,RhoQ2_comp); // This effectively nudges to 0
                rhs_arr(i,j,k,RhoQ2_comp   ) += Factor * delta;
                rhs_arr(i,j,k,RhoTheta_comp) += Factor * delta * cond_fac / exner_pi;

                // qi and heat source term due to sublimation
                if (nq > 3) { // nq > 3 guarantees that RhoQ3 is ice
                    delta = -data_arr(i,j,k,RhoQ3_comp); // This effectively nudges to 0
                    rhs_arr(i,j,k,RhoQ3_comp   ) += Factor * delta;
                    rhs_arr(i,j,k,RhoTheta_comp) += Factor * delta * sub_fac / exner_pi;
                }
            } // moist_nudge_type
        } else {
            for (int n = 0; n < nq; n++) {
                amrex::Real delta = arr_yhi(i,j,k,n) - data_arr(i,j,k,n+icomp);
                rhs_arr(i,j,k,n+icomp) += Factor*delta;
            }
        }
    });
}

/**
 * Apply Davies relaxation to explicitly mapped moisture components.  Target
 * component n corresponds to state component moisture_comps[n].
 */
AMREX_GPU_HOST
AMREX_FORCE_INLINE
void
realbdy_compute_moisture_relaxation (
                            const amrex::GpuArray<int,6>& moisture_comps,
                            const int n_targets,
                            const int& width,
                            const amrex::GpuArray<amrex::Real,AMREX_SPACEDIM>& dx,
                            const amrex::GpuArray<amrex::Real,AMREX_SPACEDIM>& ProbLo,
                            const amrex::GpuArray<amrex::Real,AMREX_SPACEDIM>& ProbHi,
                            const amrex::Real& F1,
                            const amrex::Box& bx_xlo,
                            const amrex::Box& bx_xhi,
                            const amrex::Box& bx_ylo,
                            const amrex::Box& bx_yhi,
                            const amrex::Array4<const amrex::Real>& arr_xlo,
                            const amrex::Array4<const amrex::Real>& arr_xhi,
                            const amrex::Array4<const amrex::Real>& arr_ylo,
                            const amrex::Array4<const amrex::Real>& arr_yhi,
                            const amrex::Array4<const amrex::Real>& data_arr,
                            const amrex::Array4<amrex::Real>& rhs_arr)
{
    const amrex::IntVect iv = bx_xlo.type();
    const amrex::Real ioff = (iv[0] == 1) ? zero : myhalf;
    const amrex::Real joff = (iv[1] == 1) ? zero : myhalf;

    amrex::ParallelFor(bx_xlo, bx_xhi,
    [=] AMREX_GPU_DEVICE (int i, int j, int k) noexcept
    {
        const amrex::Real x = ProbLo[0] + (i + ioff) * dx[0];
        const amrex::Real y = ProbLo[1] + (j + joff) * dx[1];
        const amrex::Real x_end = ProbLo[0] + width * dx[0];
        const amrex::Real y_end = ProbLo[1] + width * dx[1];
        const amrex::Real y_strt = ProbHi[1] - width * dx[1];
        const amrex::Real xi = (x_end - x) / (x_end - ProbLo[0]);
        const amrex::Real eta_lo = (y < y_end) ? (y_end - y) / (y_end - ProbLo[1]) : amrex::Real(0.0);
        const amrex::Real eta_hi = (y > y_strt) ? (y - y_strt) / (ProbHi[1] - y_strt) : amrex::Real(0.0);
        const amrex::Real factor = amrex::max(xi*xi, amrex::max(eta_lo,eta_hi)*amrex::max(eta_lo,eta_hi)) * F1;
        for (int n = 0; n < n_targets; ++n) {
            const int comp = moisture_comps[n];
            if (comp >= 0) {
                rhs_arr(i,j,k,comp) += factor * (arr_xlo(i,j,k,n) - data_arr(i,j,k,comp));
            }
        }
    },
    [=] AMREX_GPU_DEVICE (int i, int j, int k) noexcept
    {
        const amrex::Real x = ProbLo[0] + (i + ioff) * dx[0];
        const amrex::Real y = ProbLo[1] + (j + joff) * dx[1];
        const amrex::Real x_strt = ProbHi[0] - width * dx[0];
        const amrex::Real y_end = ProbLo[1] + width * dx[1];
        const amrex::Real y_strt = ProbHi[1] - width * dx[1];
        const amrex::Real xi = (x - x_strt) / (ProbHi[0] - x_strt);
        const amrex::Real eta_lo = (y < y_end) ? (y_end - y) / (y_end - ProbLo[1]) : amrex::Real(0.0);
        const amrex::Real eta_hi = (y > y_strt) ? (y - y_strt) / (ProbHi[1] - y_strt) : amrex::Real(0.0);
        const amrex::Real eta = amrex::max(eta_lo,eta_hi);
        const amrex::Real factor = amrex::max(xi*xi, eta*eta) * F1;
        for (int n = 0; n < n_targets; ++n) {
            const int comp = moisture_comps[n];
            if (comp >= 0) {
                rhs_arr(i,j,k,comp) += factor * (arr_xhi(i,j,k,n) - data_arr(i,j,k,comp));
            }
        }
    });

    amrex::ParallelFor(bx_ylo, bx_yhi,
    [=] AMREX_GPU_DEVICE (int i, int j, int k) noexcept
    {
        const amrex::Real y = ProbLo[1] + (j + joff) * dx[1];
        const amrex::Real y_end = ProbLo[1] + width * dx[1];
        const amrex::Real eta = (y_end - y) / (y_end - ProbLo[1]);
        const amrex::Real factor = eta*eta * F1;
        for (int n = 0; n < n_targets; ++n) {
            const int comp = moisture_comps[n];
            if (comp >= 0) {
                rhs_arr(i,j,k,comp) += factor * (arr_ylo(i,j,k,n) - data_arr(i,j,k,comp));
            }
        }
    },
    [=] AMREX_GPU_DEVICE (int i, int j, int k) noexcept
    {
        const amrex::Real y = ProbLo[1] + (j + joff) * dx[1];
        const amrex::Real y_strt = ProbHi[1] - width * dx[1];
        const amrex::Real eta = (y - y_strt) / (ProbHi[1] - y_strt);
        const amrex::Real factor = eta*eta * F1;
        for (int n = 0; n < n_targets; ++n) {
            const int comp = moisture_comps[n];
            if (comp >= 0) {
                rhs_arr(i,j,k,comp) += factor * (arr_yhi(i,j,k,n) - data_arr(i,j,k,comp));
            }
        }
    });
}

// Backward-compatible three-species overload used by existing unit tests and
// callers that only target qv/qc/qi.
AMREX_GPU_HOST
AMREX_FORCE_INLINE
void
realbdy_compute_moisture_relaxation (
                            const amrex::GpuArray<int,3>& moisture_comps,
                            const int& width,
                            const amrex::GpuArray<amrex::Real,AMREX_SPACEDIM>& dx,
                            const amrex::GpuArray<amrex::Real,AMREX_SPACEDIM>& ProbLo,
                            const amrex::GpuArray<amrex::Real,AMREX_SPACEDIM>& ProbHi,
                            const amrex::Real& F1,
                            const amrex::Box& bx_xlo,
                            const amrex::Box& bx_xhi,
                            const amrex::Box& bx_ylo,
                            const amrex::Box& bx_yhi,
                            const amrex::Array4<const amrex::Real>& arr_xlo,
                            const amrex::Array4<const amrex::Real>& arr_xhi,
                            const amrex::Array4<const amrex::Real>& arr_ylo,
                            const amrex::Array4<const amrex::Real>& arr_yhi,
                            const amrex::Array4<const amrex::Real>& data_arr,
                            const amrex::Array4<amrex::Real>& rhs_arr)
{
    const amrex::GpuArray<int,6> extended_comps = {
        moisture_comps[0], moisture_comps[1], moisture_comps[2], -1, -1, -1};
    realbdy_compute_moisture_relaxation(
        extended_comps, 3, width, dx, ProbLo, ProbHi, F1,
        bx_xlo, bx_xhi, bx_ylo, bx_yhi,
        arr_xlo, arr_xhi, arr_ylo, arr_yhi, data_arr, rhs_arr);
}

/*
 * Effectively a Multiply for a MultiFab and an iMultiFab mask
 */
AMREX_GPU_HOST
AMREX_FORCE_INLINE
void
ApplyMask (amrex::MultiFab& dst,
           const amrex::iMultiFab& imask,
           const int nghost = 0)
{
    for (amrex::MFIter mfi(dst,amrex::TilingIfNotGPU()); mfi.isValid(); ++mfi)
    {
        const amrex::Box& bx = mfi.growntilebox(nghost);
        if (bx.ok())
        {
            auto        dstFab = dst.array(mfi);
            const auto maskFab = imask.const_array(mfi);
            ParallelFor(bx, [=] AMREX_GPU_DEVICE (int i, int j, int k) noexcept
            {
                dstFab(i,j,k) *= maskFab(i,j,k);
            });
        }
    }
}

/**
 * Multiply a MultiFab by the inverse of an integer mask.
 * @param[in,out] dst MultiFab to be masked.
 * @param[in] imask Integer mask MultiFab.
 * @param[in] nghost Number of ghost cells to include.
 */
AMREX_GPU_HOST
AMREX_FORCE_INLINE
void
ApplyInvertedMask (amrex::MultiFab& dst,
                   const amrex::iMultiFab& imask,
                   const int nghost = 0)
{
    for (amrex::MFIter mfi(dst,amrex::TilingIfNotGPU()); mfi.isValid(); ++mfi)
    {
        const amrex::Box& bx = mfi.growntilebox(nghost);
        if (bx.ok())
        {
            auto        dstFab = dst.array(mfi);
            const auto maskFab = imask.const_array(mfi);
            ParallelFor(bx, [=] AMREX_GPU_DEVICE (int i, int j, int k) noexcept
            {
                dstFab(i,j,k) *= (1-maskFab(i,j,k));
            });
        }
    }
}

/**
 * Compute the distance to thin-body walls.
 * @param[out] wdist MultiFab to store wall distances.
 * @param[in] xfaces Indices of x-faces forming walls.
 * @param[in] yfaces Indices of y-faces forming walls.
 * @param[in] zfaces Indices of z-faces forming walls.
 * @param[in] geomdata Geometry used for physical distance.
 * @param[in] z_phys_cc Cell-centered physical height.
 */
void
thinbody_wall_dist (std::unique_ptr<amrex::MultiFab>& wdist,
                    amrex::Vector<amrex::IntVect>& xfaces,
                    amrex::Vector<amrex::IntVect>& yfaces,
                    amrex::Vector<amrex::IntVect>& zfaces,
                    const amrex::Geometry& geomdata,
                    std::unique_ptr<amrex::MultiFab>& z_phys_cc);


void WeatherDataInterpolation (const double time);

/**
 * Perform quadratic interpolation in one dimension.
 * @param[in] z Interpolation coordinate.
 * @param[in] z0 Coordinate of first point.
 * @param[in] p0 Value at first point.
 * @param[in] z1 Coordinate of second point.
 * @param[in] p1 Value at second point.
 * @param[in] z2 Coordinate of third point.
 * @param[in] p2 Value at third point.
 * @return Interpolated value.
 */
AMREX_GPU_HOST_DEVICE
AMREX_FORCE_INLINE
amrex::Real
quad_interp_1d (amrex::Real z,
                amrex::Real z0, amrex::Real p0,
                amrex::Real z1, amrex::Real p1,
                amrex::Real z2, amrex::Real p2)
{
    return p0 * ( (z - z1) * (z - z2) ) / ( (z0 - z1) * (z0 - z2) )
         + p1 * ( (z - z0) * (z - z2) ) / ( (z1 - z0) * (z1 - z2) )
         + p2 * ( (z - z0) * (z - z1) ) / ( (z2 - z0) * (z2 - z1) );
}

/**
 * Which time level(s) a fillpatch at a given time should draw from.
 */
enum struct TimeLevelSel { Old, New, Interp };

/**
 * Decide whether data at time "time" should come from the old-time data, the new-time data,
 * or an interpolation between the two.
 *
 * This is done in double precision, before any cast to Real, and with a single tolerance
 * at both ends.  The tolerance is a small fraction of (t_new - t_old), but never smaller than
 * a few ulps of the times themselves, since those carry roundoff proportional to their
 * magnitude -- which matters once the time is large (e.g. months measured in seconds).
 *
 * Times outside [t_old, t_new] are clamped to the nearer end rather than extrapolated.
 * If t_old == t_new (to within the tolerance), Old is returned.
 *
 * @param[in] time  time at which data are requested
 * @param[in] t_old old time of the level
 * @param[in] t_new new time of the level
 */
inline TimeLevelSel
select_time_level (double time, double t_old, double t_new)
{
    const double t_mag = std::max({std::abs(time), std::abs(t_old), std::abs(t_new)});
    const double tol   = std::max(1.e-8 * std::abs(t_new - t_old),
                                  16.0 * std::numeric_limits<double>::epsilon() * t_mag);

    if (time <= t_old + tol) { return TimeLevelSel::Old; }
    if (time >= t_new - tol) { return TimeLevelSel::New; }
    return TimeLevelSel::Interp;
}

/**
 * The pair of MultiFabs to hand to an AMReX FillPatch routine for a given time-level selection.
 * When only one time level is wanted it is passed twice, so that AMReX's own comparison of the
 * (rounded) times can never bring in the other one.
 *
 * @param[in] sel    which time level(s) to use
 * @param[in] old_mf data at the old time
 * @param[in] new_mf data at the new time
 */
inline amrex::Vector<amrex::MultiFab*>
select_time_data (TimeLevelSel sel, amrex::MultiFab& old_mf, amrex::MultiFab& new_mf)
{
    if (sel == TimeLevelSel::Old) { return {&old_mf, &old_mf}; }
    if (sel == TimeLevelSel::New) { return {&new_mf, &new_mf}; }
    return {&old_mf, &new_mf};
}

#endif
