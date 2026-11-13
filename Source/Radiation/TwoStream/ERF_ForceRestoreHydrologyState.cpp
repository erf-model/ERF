#include <ERF_ForceRestoreHydrologyState.H>

#include <AMReX_MFIter.H>
#include <AMReX_ParallelDescriptor.H>
#include <AMReX_PlotFileUtil.H>
#include <AMReX_Reduce.H>
#include <AMReX_Utility.H>
#include <AMReX_VisMF.H>
#include <AMReX_MultiFabUtil.H>  // amrex::cast, used by ERF_ReadCheckpointMF.H
#include <AMReX_iMultiFab.H>

#include <ERF_Constants.H>          // L_v
#include <ERF_RadStruct.H>
#include <ERF_ReadCheckpointMF.H>

using namespace amrex;

ForceRestoreHydrologyParams
ForceRestoreHydrology::params (const RadChoice& rc)
{
    ForceRestoreHydrologyParams p;
    p.rain  = rc.seb_precipitation_source != SEBPrecipitationSource::None;
    p.deep  = rc.seb_deep_prognostic_enable;
    p.drain = rc.seb_drainage_enable;
    p.d_s   = rc.seb_moisture_layer_depth_m;
    p.d_2   = rc.seb_deep_layer_depth_m;
    p.tau   = rc.seb_restore_timescale_s;
    p.w_fc  = rc.seb_soil_moisture_fc;
    p.w_sat = rc.seb_soil_porosity;
    // MN96 Eq. (16): C_3 scales as 1/d_2; Table 1 gives C_3 d_2 [m].
    p.c3    = p.drain ? rc.seb_drainage_c3_m / rc.seb_deep_layer_depth_m : Real(0.0);
    p.q_min = rc.seb_prognostic_q_min;
    p.q_max = rc.seb_prognostic_q_max;
    return p;
}

void
ForceRestoreHydrology::resize (int nlevs_max)
{
    m_accum_start.resize(nlevs_max);
    m_rain_rate.resize(nlevs_max);
    m_totals.resize(nlevs_max);
    m_started.resize(nlevs_max, 0);
    m_unused.resize(nlevs_max, 0);
    m_elapsed.resize(nlevs_max, 0.0);
}

void
ForceRestoreHydrology::define_level (int lev, const RadChoice& rc, const BoxArray& ba2d,
                                     const DistributionMapping& dm, const IntVect& ng)
{
    if (!rc.seb_hydrology_enabled()) { return; }
    m_accum_start[lev] = std::make_unique<MultiFab>(ba2d, dm, 1, 0);
    m_rain_rate[lev]   = std::make_unique<MultiFab>(ba2d, dm, 1, 0);
    m_totals[lev]      = std::make_unique<MultiFab>(ba2d, dm, NumTotals, ng);
    m_accum_start[lev]->setVal(0.0);
    m_rain_rate[lev]->setVal(0.0);
    m_totals[lev]->setVal(0.0);
    m_started[lev] = 0;
    m_unused[lev] = 0;
    m_elapsed[lev] = 0.0;
}

void
ForceRestoreHydrology::begin_step (int lev, const SurfacePrecipAccumulationSources& sources, int klo)
{
    if (!defined(lev)) { return; }
    if (m_unused[lev] == 0) { surface_precip_accum_k0(*m_accum_start[lev], sources, klo); }
    m_started[lev] = 1;
}

void
ForceRestoreHydrology::end_step (int lev, const SurfacePrecipAccumulationSources& sources, int klo,
                                 Real dt)
{
    if (!defined(lev)) { return; }
    if (m_started[lev] == 0) { m_rain_rate[lev]->setVal(0.0); return; }
    MultiFab accum_end(m_accum_start[lev]->boxArray(), m_accum_start[lev]->DistributionMap(), 1, 0);
    surface_precip_accum_k0(accum_end, sources, klo);
    // The time since the stored start: this step, plus any step whose rain step() did not use.
    m_elapsed[lev] += dt;
    const Real elapsed = m_elapsed[lev];
    for (MFIter mfi(*m_rain_rate[lev], TilingIfNotGPU()); mfi.isValid(); ++mfi) {
        const Box& bx = mfi.tilebox();
        Array4<Real> r = m_rain_rate[lev]->array(mfi);
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

ForceRestoreHydrology::StepSums
ForceRestoreHydrology::step (int lev, const MultiFab& t_s, MultiFab& q_s,
                             MultiFab& q_deep, MultiFab& t_deep,
                             const MultiFab& lh, const MultiFab* footprint,
                             Real dt, const ForceRestoreHydrologyParams& prm)
{
    StepSums sums;
    if (!defined(lev)) { return sums; }

    // The rain rate is a mean over m_elapsed, which exceeds dt when an earlier step's rain was
    // carried: the step takes the whole of it, P m_elapsed, in its P dt.
    const Real rain_scale = (prm.rain && m_elapsed[lev] > Real(0.0)) ? m_elapsed[lev] / dt : Real(1.0);

    ReduceOps<ReduceOpSum, ReduceOpSum, ReduceOpSum, ReduceOpSum, ReduceOpSum, ReduceOpMax,
              ReduceOpSum> reduce_op;
    ReduceData<Real, Real, Real, Real, Real, Real, Long> reduce_data(reduce_op);
    using ReduceTuple = typename decltype(reduce_data)::Type;

    // Untiled, on the valid boxes: the update is in place.
    for (MFIter mfi(q_s, false); mfi.isValid(); ++mfi) {
        const Box& bx = mfi.validbox();
        Array4<const Real> ts  = t_s.const_array(mfi);
        Array4<Real>       qs  = q_s.array(mfi);
        Array4<Real>       qd  = q_deep.array(mfi);
        Array4<Real>       td  = t_deep.array(mfi);
        Array4<const Real> le  = lh.const_array(mfi);
        Array4<const Real> pr  = m_rain_rate[lev]->const_array(mfi);
        Array4<Real>       tot = m_totals[lev]->array(mfi);
        const bool has_foot = footprint != nullptr;
        Array4<const Real> foot = has_foot ? footprint->const_array(mfi) : Array4<const Real>{};

        reduce_op.eval(bx, reduce_data,
        [=] AMREX_GPU_DEVICE (int i, int j, int k) noexcept -> ReduceTuple
        {
            if (has_foot && foot(i, j, k) > Real(0.5)) {
                return {Real(0.), Real(0.), Real(0.), Real(0.), Real(0.), Real(0.), Long(0)};
            }
            Real p = prm.rain ? pr(i, j, k) * rain_scale : Real(0.0);
            if (!amrex::Math::isfinite(p)) { p = Real(0.0); }
            Real e = le(i, j, k) / L_v;
            if (!amrex::Math::isfinite(e)) { e = Real(0.0); }

            const Real qd_old = qd(i, j, k);
            Real qs_c = qs(i, j, k), qd_c = qd_old, td_c = td(i, j, k);
            const ForceRestoreHydrologyFluxes f =
                force_restore_hydrology_step(qs_c, qd_c, td_c, ts(i, j, k), p, e, dt, prm);
            qs(i, j, k) = qs_c;
            if (prm.deep) { qd(i, j, k) = qd_c; td(i, j, k) = td_c; }

            tot(i, j, k, Rain)     += p * dt;
            tot(i, j, k, Runoff)   += f.runoff;
            tot(i, j, k, Drainage) += f.drainage;
            tot(i, j, k, Deficit)  += f.deficit;

            // The column's water budget, which closes to round-off with the deep reservoir.
            Real residual = Real(0.0);
            if (prm.deep) {
                residual = amrex::Math::abs(rhor * prm.d_2 * (qd_c - qd_old)
                                            - ((p - e) * dt - f.drainage - f.runoff + f.deficit));
            }
            return {p * dt, e * dt, f.runoff, f.drainage, f.deficit, residual, Long(1)};
        });
    }

    m_unused[lev] = 0;
    m_elapsed[lev] = 0.0;

    ReduceTuple h = reduce_data.value(reduce_op);
    Real s[5] = {amrex::get<0>(h), amrex::get<1>(h), amrex::get<2>(h), amrex::get<3>(h),
                 amrex::get<4>(h)};
    Real rmax = amrex::get<5>(h);
    Long n = amrex::get<6>(h);
    ParallelDescriptor::ReduceRealSum(s, 5);
    ParallelDescriptor::ReduceRealMax(rmax);
    ParallelDescriptor::ReduceLongSum(n);
    sums.rain = s[0]; sums.evaporation = s[1]; sums.runoff = s[2];
    sums.drainage = s[3]; sums.deficit = s[4]; sums.max_residual = rmax; sums.columns = n;
    return sums;
}

void
ForceRestoreHydrology::write_checkpoint (int lev, const std::string& checkpointname) const
{
    if (!defined(lev)) { return; }
    VisMF::Write(*m_totals[lev],
                 MultiFabFileFullPrefix(lev, checkpointname, "Level_", "TwoStream_HydrologyTotals"));
}

void
ForceRestoreHydrology::read_checkpoint (int lev, const std::string& restart_chkfile,
                                        const Periodicity& period)
{
    if (!defined(lev)) { return; }
    const std::string name =
        MultiFabFileFullPrefix(lev, restart_chkfile, "Level_", "TwoStream_HydrologyTotals");
    // A checkpoint written without the hydrology has no totals: they restart from zero.
    if (!amrex::FileExists(name + "_H")) {
        amrex::Print() << "NOTE: the checkpoint has no force-restore hydrology totals at level "
                       << lev << "; seb_rain_total, seb_runoff_total, seb_drainage_total and "
                       << "seb_deficit_total restart from zero.\n";
        return;
    }
    erf_checkpoint::read_checkpoint_multifab(*m_totals[lev], name, 0, 0, NumTotals,
                                             m_totals[lev]->nGrowVect(), true, period);
}
