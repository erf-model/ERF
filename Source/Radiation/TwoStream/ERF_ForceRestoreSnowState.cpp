#include <ERF_ForceRestoreSnowState.H>

#include <limits>

#include <AMReX_MFIter.H>
#include <AMReX_ParallelDescriptor.H>
#include <AMReX_PlotFileUtil.H>
#include <AMReX_Reduce.H>
#include <AMReX_Utility.H>
#include <AMReX_VisMF.H>
#include <AMReX_MultiFabUtil.H>  // amrex::cast, used by ERF_ReadCheckpointMF.H
#include <AMReX_iMultiFab.H>

#include <ERF_Constants.H>               // L_v
#include <ERF_RadStruct.H>
#include <ERF_SEBLandMask.H>
#include <ERF_ReadCheckpointMF.H>

using namespace amrex;

ForceRestoreSnowParams
ForceRestoreSnow::params (const RadChoice& rc)
{
    ForceRestoreSnowParams p;
    p.w_crn  = rc.seb_snow_critical_swe_kg_m2;
    p.albedo = rc.seb_snow_albedo;
    p.c_s    = rc.seb_surface_heat_capacity;
    p.t_max  = rc.seb_prognostic_t_max_k;
    return p;
}

void
ForceRestoreSnow::resize (int nlevs_max)
{
    m_swe.resize(nlevs_max);
    m_totals.resize(nlevs_max);
    m_soil_water.resize(nlevs_max);
}

void
ForceRestoreSnow::define_level (int lev, const RadChoice& rc, const BoxArray& ba2d,
                                const DistributionMapping& dm, const IntVect& ng)
{
    if (!rc.seb_snow_enable) { return; }
    m_swe[lev]         = std::make_unique<MultiFab>(ba2d, dm, 1, ng);
    m_totals[lev]      = std::make_unique<MultiFab>(ba2d, dm, NumTotals, ng);
    m_soil_water[lev]  = std::make_unique<MultiFab>(ba2d, dm, 1, 0);
    m_swe[lev]->setVal(rc.seb_snow_initial_swe_kg_m2);
    m_totals[lev]->setVal(0.0);
    m_soil_water[lev]->setVal(0.0);
}

ForceRestoreSnow::StepSums
ForceRestoreSnow::step (int lev, MultiFab& t_s, const MultiFab& lh, const MultiFab* footprint,
                        Real dt, const ForceRestoreSnowParams& prm,
                        const ForceRestorePrecipInput& snowfall,
                        Real& t_s_sum_local, Real& t_s_max_local,
                        const iMultiFab* land_mask)
{
    StepSums sums;
    t_s_sum_local = Real(0.0);
    t_s_max_local = std::numeric_limits<Real>::lowest();
    if (!defined(lev)) { return sums; }
    seb_assert_land_mask_layout(land_mask, t_s, "ForceRestoreSnow::step");

    // The snowfall of the step: the rate times its scale (ERF_ForceRestorePrecipRate.H).
    const bool has_snowfall = snowfall.rate != nullptr;
    const Real snow_scale = snowfall.scale;
    const int snow_comp = snowfall.comp;

    // The last sum counts the columns the step left in a state the caller refuses.
    ReduceOps<ReduceOpSum, ReduceOpSum, ReduceOpSum, ReduceOpSum, ReduceOpMax, ReduceOpSum,
              ReduceOpSum, ReduceOpMax, ReduceOpSum> reduce_op;
    ReduceData<Real, Real, Real, Real, Real, Long, Real, Real, Long> reduce_data(reduce_op);
    using ReduceTuple = typename decltype(reduce_data)::Type;

    // Untiled, on the valid boxes: the update is in place.
    for (MFIter mfi(*m_swe[lev], false); mfi.isValid(); ++mfi) {
        const Box& bx = mfi.validbox();
        Array4<Real>       ts  = t_s.array(mfi);
        Array4<Real>       wn  = m_swe[lev]->array(mfi);
        Array4<Real>       tot = m_totals[lev]->array(mfi);
        Array4<Real>       out = m_soil_water[lev]->array(mfi);
        Array4<const Real> le  = lh.const_array(mfi);
        Array4<const Real> pn  = has_snowfall ? snowfall.rate->const_array(mfi)
                                              : Array4<const Real>{};
        const bool has_foot = footprint != nullptr;
        Array4<const Real> foot = has_foot ? footprint->const_array(mfi) : Array4<const Real>{};
        const SEBLandMaskView mask = seb_land_mask_view(land_mask, mfi);

        reduce_op.eval(bx, reduce_data,
        [=] AMREX_GPU_DEVICE (int i, int j, int k) noexcept -> ReduceTuple
        {
            if ((has_foot && foot(i, j, k) > Real(0.5)) || mask.water(i, j)) {
                // No snow lies on a building footprint, and the layer is not kept on water
                // (the balance does not run there): the initial cover is taken off.
                wn(i, j, k) = Real(0.0);
                out(i, j, k) = Real(0.0);
                return {Real(0.), Real(0.), Real(0.), Real(0.), Real(0.), Long(0), Real(0.),
                        std::numeric_limits<Real>::lowest(), Long(0)};
            }
            const Real p = has_snowfall ? pn(i, j, k, snow_comp) * snow_scale : Real(0.0);
            const Real e = le(i, j, k) / L_v;

            const Real w_old = wn(i, j, k);
            const Real t_old = ts(i, j, k);
            Real w_c = w_old, t_c = t_old;
            const ForceRestoreSnowFluxes f = force_restore_snow_step(w_c, t_c, p, e, dt, prm);
            wn(i, j, k) = w_c;
            ts(i, j, k) = t_c;
            out(i, j, k) = f.melt + f.sublimation;

            tot(i, j, k, Snowfall)    += f.snowfall;
            tot(i, j, k, Sublimation) += f.sublimation;
            tot(i, j, k, Melt)        += f.melt;

            // The column's snow budget, DRM95 Eq. (6), which closes to round-off.
            const Real residual = amrex::Math::abs((w_c - w_old)
                                                   - (f.snowfall - f.sublimation - f.melt));
            // A column is refused on what the step was given as well as on what it leaves: the
            // bounds of the step (the melt, the upper bound of T_s) can hide a bad input.
            const bool bad = !amrex::Math::isfinite(t_old) || !amrex::Math::isfinite(w_old) ||
                             w_old < Real(0.0) ||
                             !amrex::Math::isfinite(t_c) || !amrex::Math::isfinite(w_c) ||
                             w_c < Real(0.0);
            return {f.snowfall, f.sublimation, f.melt, w_c, residual, Long(1), t_c,
                    amrex::Math::abs(t_c), bad ? Long(1) : Long(0)};
        });
    }

    ReduceTuple h = reduce_data.value(reduce_op);
    Real s[4] = {amrex::get<0>(h), amrex::get<1>(h), amrex::get<2>(h), amrex::get<3>(h)};
    Real rmax = amrex::get<4>(h);
    // The two counts go through one collective.
    Long counts[2] = {amrex::get<5>(h), amrex::get<8>(h)};
    t_s_sum_local = amrex::get<6>(h);
    t_s_max_local = amrex::get<7>(h);
    ParallelDescriptor::ReduceRealSum(s, 4);
    ParallelDescriptor::ReduceRealMax(rmax);
    ParallelDescriptor::ReduceLongSum(counts, 2);
    sums.snowfall = s[0]; sums.sublimation = s[1]; sums.melt = s[2]; sums.swe = s[3];
    sums.max_residual = rmax; sums.columns = counts[0]; sums.bad_columns = counts[1];
    return sums;
}

void
ForceRestoreSnow::apply_albedo (int lev, MultiFab& alb_sw, const ForceRestoreSnowParams& prm,
                                const MultiFab* footprint, const iMultiFab* land_mask) const
{
    if (!defined(lev)) { return; }
    seb_assert_land_mask_layout(land_mask, alb_sw, "ForceRestoreSnow::apply_albedo");
    // Valid cells only: step() does not refresh the ghost cells of W_n, and the sweep reads
    // the albedo on valid cells.
    for (MFIter mfi(alb_sw, TilingIfNotGPU()); mfi.isValid(); ++mfi) {
        const Box& bx = mfi.tilebox();
        Array4<Real> a = alb_sw.array(mfi);
        Array4<const Real> wn = m_swe[lev]->const_array(mfi);
        const bool has_foot = footprint != nullptr;
        Array4<const Real> foot = has_foot ? footprint->const_array(mfi) : Array4<const Real>{};
        const SEBLandMaskView mask = seb_land_mask_view(land_mask, mfi);
        ParallelFor(bx, [=] AMREX_GPU_DEVICE (int i, int j, int k) noexcept
        {
            // step() takes the snow off a footprint, but W_n there still holds the initial
            // value before the first step (at the start, and after a restart from a
            // checkpoint without the snow layer).
            if (has_foot && foot(i, j, k) > Real(0.5)) { return; }
            // Nor on a water column, the same columns step() leaves out.
            if (mask.water(i, j)) { return; }
            a(i, j, k) = force_restore_snow_albedo(a(i, j, k), wn(i, j, k), prm);
        });
    }
}

void
ForceRestoreSnow::write_checkpoint (int lev, const std::string& checkpointname) const
{
    if (!defined(lev)) { return; }
    VisMF::Write(*m_swe[lev],
                 MultiFabFileFullPrefix(lev, checkpointname, "Level_", "TwoStream_SnowSWE"));
    VisMF::Write(*m_totals[lev],
                 MultiFabFileFullPrefix(lev, checkpointname, "Level_", "TwoStream_SnowTotals"));
}

void
ForceRestoreSnow::read_checkpoint (int lev, const std::string& restart_chkfile,
                                   const Periodicity& period)
{
    if (!defined(lev)) { return; }
    const std::string swe_name =
        MultiFabFileFullPrefix(lev, restart_chkfile, "Level_", "TwoStream_SnowSWE");
    const std::string totals_name =
        MultiFabFileFullPrefix(lev, restart_chkfile, "Level_", "TwoStream_SnowTotals");
    const bool has_swe = amrex::FileExists(swe_name + "_H");
    const bool has_totals = amrex::FileExists(totals_name + "_H");
    // A checkpoint written without the snow layer has no snow: it starts from the input.
    if (!has_swe && !has_totals) {
        // Once, not per level.
        if (lev == 0) {
            amrex::Print() << "NOTE: the checkpoint has no force-restore snow; seb_snow_swe "
                           << "starts from erf.radiation.seb_snow_initial_swe_kg_m2 and the "
                           << "snow totals restart from zero.\n";
        }
        return;
    }
    // One file without the other is not a checkpoint without the snow layer: starting the
    // snow from the input here would drop the snow the run had.
    if (has_swe != has_totals) {
        amrex::Abort("The checkpoint " + restart_chkfile + " holds the force-restore snow file " +
                     (has_swe ? swe_name : totals_name) + " but not " +
                     (has_swe ? totals_name : swe_name) + ": the checkpoint is incomplete");
    }
    erf_checkpoint::read_checkpoint_multifab(*m_swe[lev], swe_name, 0, 0, 1,
                                             m_swe[lev]->nGrowVect(), true, period);
    erf_checkpoint::read_checkpoint_multifab(*m_totals[lev], totals_name, 0, 0, NumTotals,
                                             m_totals[lev]->nGrowVect(), true, period);
}
