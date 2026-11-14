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
#include <ERF_ForceRestoreHydrology.H>   // force_restore_rain_rate
#include <ERF_RadStruct.H>
#include <ERF_ReadCheckpointMF.H>

using namespace amrex;

ForceRestoreSnowParams
ForceRestoreSnow::params (const RadChoice& rc)
{
    ForceRestoreSnowParams p;
    p.w_crn  = rc.seb_snow_critical_swe_kg_m2;
    p.albedo = rc.seb_snow_albedo;
    p.c_s    = rc.seb_surface_heat_capacity;
    return p;
}

void
ForceRestoreSnow::resize (int nlevs_max)
{
    m_swe.resize(nlevs_max);
    m_totals.resize(nlevs_max);
    m_accum_start.resize(nlevs_max);
    m_frozen_rate.resize(nlevs_max);
    m_soil_water.resize(nlevs_max);
    m_started.resize(nlevs_max, 0);
    m_unused.resize(nlevs_max, 0);
    m_elapsed.resize(nlevs_max, 0.0);
}

void
ForceRestoreSnow::define_level (int lev, const RadChoice& rc, const BoxArray& ba2d,
                                const DistributionMapping& dm, const IntVect& ng)
{
    if (!rc.seb_snow_enable) { return; }
    m_swe[lev]         = std::make_unique<MultiFab>(ba2d, dm, 1, ng);
    m_totals[lev]      = std::make_unique<MultiFab>(ba2d, dm, NumTotals, ng);
    m_accum_start[lev] = std::make_unique<MultiFab>(ba2d, dm, 1, 0);
    m_frozen_rate[lev] = std::make_unique<MultiFab>(ba2d, dm, 1, 0);
    m_soil_water[lev]  = std::make_unique<MultiFab>(ba2d, dm, 1, 0);
    m_swe[lev]->setVal(rc.seb_snow_initial_swe_kg_m2);
    m_totals[lev]->setVal(0.0);
    m_accum_start[lev]->setVal(0.0);
    m_frozen_rate[lev]->setVal(0.0);
    m_soil_water[lev]->setVal(0.0);
    m_started[lev] = 0;
    m_unused[lev] = 0;
    m_elapsed[lev] = 0.0;
}

void
ForceRestoreSnow::begin_step (int lev, const SurfacePrecipAccumulationSources& sources, int klo)
{
    if (!defined(lev)) { return; }
    if (!m_noted_no_frozen &&
        !(surface_precip_has_source(sources.snow) || surface_precip_has_source(sources.graupel) ||
          surface_precip_has_source(sources.hail))) {
        m_noted_no_frozen = true;
        amrex::Print() << "NOTE: erf.radiation.seb_snow_enable is on, but the microphysics scheme "
                          "has no snow, graupel or hail accumulation at the surface: no snow will "
                          "fall, and the snow is that of erf.radiation.seb_snow_initial_swe_kg_m2\n";
    }
    if (m_unused[lev] == 0) {
        surface_precip_accum_k0(*m_accum_start[lev], sources, klo, SurfacePrecipPhase::Frozen);
    }
    m_started[lev] = 1;
}

void
ForceRestoreSnow::end_step (int lev, const SurfacePrecipAccumulationSources& sources, int klo,
                            Real dt)
{
    if (!defined(lev)) { return; }
    if (m_started[lev] == 0) { m_frozen_rate[lev]->setVal(0.0); return; }
    MultiFab accum_end(m_accum_start[lev]->boxArray(), m_accum_start[lev]->DistributionMap(), 1, 0);
    surface_precip_accum_k0(accum_end, sources, klo, SurfacePrecipPhase::Frozen);
    // The time since the stored start: this step, plus any step whose snowfall step() did
    // not use.
    m_elapsed[lev] += dt;
    const Real elapsed = m_elapsed[lev];
    for (MFIter mfi(*m_frozen_rate[lev], TilingIfNotGPU()); mfi.isValid(); ++mfi) {
        const Box& bx = mfi.tilebox();
        Array4<Real> r = m_frozen_rate[lev]->array(mfi);
        Array4<const Real> a0 = m_accum_start[lev]->const_array(mfi);
        Array4<const Real> a1 = accum_end.const_array(mfi);
        ParallelFor(bx, [=] AMREX_GPU_DEVICE (int i, int j, int k) noexcept
        {
            r(i, j, k) = force_restore_rain_rate(a1(i, j, k), a0(i, j, k), elapsed);
        });
    }
    m_started[lev] = 0;
    m_unused[lev] = 1;
}

ForceRestoreSnow::StepSums
ForceRestoreSnow::step (int lev, MultiFab& t_s, const MultiFab& lh, const MultiFab* footprint,
                        Real dt, const ForceRestoreSnowParams& prm,
                        Real& t_s_sum_local, Real& t_s_max_local)
{
    StepSums sums;
    t_s_sum_local = Real(0.0);
    t_s_max_local = std::numeric_limits<Real>::lowest();
    if (!defined(lev)) { return sums; }

    // The rate is a mean over m_elapsed, which exceeds dt when an earlier step's snowfall
    // was carried: the step takes the whole of it.
    const Real snow_scale = (m_elapsed[lev] > Real(0.0)) ? m_elapsed[lev] / dt : Real(1.0);

    ReduceOps<ReduceOpSum, ReduceOpSum, ReduceOpSum, ReduceOpSum, ReduceOpMax, ReduceOpSum,
              ReduceOpSum, ReduceOpMax> reduce_op;
    ReduceData<Real, Real, Real, Real, Real, Long, Real, Real> reduce_data(reduce_op);
    using ReduceTuple = typename decltype(reduce_data)::Type;

    // Untiled, on the valid boxes: the update is in place.
    for (MFIter mfi(*m_swe[lev], false); mfi.isValid(); ++mfi) {
        const Box& bx = mfi.validbox();
        Array4<Real>       ts  = t_s.array(mfi);
        Array4<Real>       wn  = m_swe[lev]->array(mfi);
        Array4<Real>       tot = m_totals[lev]->array(mfi);
        Array4<Real>       out = m_soil_water[lev]->array(mfi);
        Array4<const Real> le  = lh.const_array(mfi);
        Array4<const Real> pn  = m_frozen_rate[lev]->const_array(mfi);
        const bool has_foot = footprint != nullptr;
        Array4<const Real> foot = has_foot ? footprint->const_array(mfi) : Array4<const Real>{};

        reduce_op.eval(bx, reduce_data,
        [=] AMREX_GPU_DEVICE (int i, int j, int k) noexcept -> ReduceTuple
        {
            if (has_foot && foot(i, j, k) > Real(0.5)) {
                // No snow lies on a building footprint: the initial cover is taken off.
                wn(i, j, k) = Real(0.0);
                out(i, j, k) = Real(0.0);
                return {Real(0.), Real(0.), Real(0.), Real(0.), Real(0.), Long(0), Real(0.),
                        std::numeric_limits<Real>::lowest()};
            }
            const Real p = pn(i, j, k) * snow_scale;
            const Real e = le(i, j, k) / L_v;

            const Real w_old = wn(i, j, k);
            Real w_c = w_old, t_c = ts(i, j, k);
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
            return {f.snowfall, f.sublimation, f.melt, w_c, residual, Long(1), t_c,
                    amrex::Math::abs(t_c)};
        });
    }

    m_unused[lev] = 0;
    m_elapsed[lev] = 0.0;

    ReduceTuple h = reduce_data.value(reduce_op);
    Real s[4] = {amrex::get<0>(h), amrex::get<1>(h), amrex::get<2>(h), amrex::get<3>(h)};
    Real rmax = amrex::get<4>(h);
    Long n = amrex::get<5>(h);
    t_s_sum_local = amrex::get<6>(h);
    t_s_max_local = amrex::get<7>(h);
    ParallelDescriptor::ReduceRealSum(s, 4);
    ParallelDescriptor::ReduceRealMax(rmax);
    ParallelDescriptor::ReduceLongSum(n);
    sums.snowfall = s[0]; sums.sublimation = s[1]; sums.melt = s[2]; sums.swe = s[3];
    sums.max_residual = rmax; sums.columns = n;
    return sums;
}

void
ForceRestoreSnow::apply_albedo (int lev, MultiFab& alb_sw, const ForceRestoreSnowParams& prm) const
{
    if (!defined(lev)) { return; }
    // Valid cells only: step() does not refresh the ghost cells of W_n, and the sweep reads
    // the albedo on valid cells.
    for (MFIter mfi(alb_sw, TilingIfNotGPU()); mfi.isValid(); ++mfi) {
        const Box& bx = mfi.tilebox();
        Array4<Real> a = alb_sw.array(mfi);
        Array4<const Real> wn = m_swe[lev]->const_array(mfi);
        ParallelFor(bx, [=] AMREX_GPU_DEVICE (int i, int j, int k) noexcept
        {
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
    // A checkpoint written without the snow layer has no snow: it starts from the input.
    if (!amrex::FileExists(swe_name + "_H") || !amrex::FileExists(totals_name + "_H")) {
        amrex::Print() << "NOTE: the checkpoint has no force-restore snow at level " << lev
                       << "; seb_snow_swe starts from erf.radiation.seb_snow_initial_swe_kg_m2 "
                       << "and the snow totals restart from zero.\n";
        return;
    }
    erf_checkpoint::read_checkpoint_multifab(*m_swe[lev], swe_name, 0, 0, 1,
                                             m_swe[lev]->nGrowVect(), true, period);
    erf_checkpoint::read_checkpoint_multifab(*m_totals[lev], totals_name, 0, 0, NumTotals,
                                             m_totals[lev]->nGrowVect(), true, period);
}
