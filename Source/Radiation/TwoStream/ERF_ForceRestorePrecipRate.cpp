#include <ERF_ForceRestorePrecipRate.H>

#include <AMReX_MFIter.H>
#include <AMReX_Print.H>

#include <ERF_ForceRestoreHydrology.H>   // force_restore_rain_rate
#include <ERF_RadStruct.H>

using namespace amrex;

void
ForceRestorePrecipRate::resize (int nlevs_max)
{
    m_accum_start.resize(nlevs_max);
    m_accum_now.resize(nlevs_max);
    m_rate.resize(nlevs_max);
    m_started.resize(nlevs_max, 0);
    m_unused.resize(nlevs_max, 0);
    m_elapsed.resize(nlevs_max, 0.0);
}

void
ForceRestorePrecipRate::define_level (int lev, const RadChoice& rc, const BoxArray& ba2d,
                                      const DistributionMapping& dm)
{
    if (rc.seb_precipitation_source == SEBPrecipitationSource::None) { return; }
    // The one place that decides what the rate holds: with the snow layer the rain and the
    // frozen part apart, otherwise the total.
    m_split = rc.seb_snow_enable;
    const int ncomp = m_split ? 2 : 1;
    m_accum_start[lev] = std::make_unique<MultiFab>(ba2d, dm, ncomp, 0);
    m_accum_now[lev]   = std::make_unique<MultiFab>(ba2d, dm, ncomp, 0);
    m_rate[lev]        = std::make_unique<MultiFab>(ba2d, dm, ncomp, 0);
    m_accum_start[lev]->setVal(0.0);
    m_accum_now[lev]->setVal(0.0);
    m_rate[lev]->setVal(0.0);
    m_started[lev] = 0;
    m_unused[lev] = 0;
    m_elapsed[lev] = 0.0;
}

void
ForceRestorePrecipRate::fill_accumulation (MultiFab& accum,
                                           const SurfacePrecipAccumulationSources& sources,
                                           int klo) const
{
    if (m_split) {
        surface_precip_liquid_frozen_accum_k0(accum, sources, klo);
    } else {
        surface_precip_accum_k0(accum, sources, klo, SurfacePrecipPhase::Total);
    }
}

void
ForceRestorePrecipRate::begin_step (int lev, const SurfacePrecipAccumulationSources& sources,
                                    int klo)
{
    if (!defined(lev)) { return; }
    if (m_split && !m_noted_no_frozen &&
        !(surface_precip_has_source(sources.snow) || surface_precip_has_source(sources.graupel) ||
          surface_precip_has_source(sources.hail))) {
        m_noted_no_frozen = true;
        amrex::Print() << "NOTE: erf.radiation.seb_snow_enable is on, but the microphysics scheme "
                          "has no snow, graupel or hail accumulation at the surface: no snow will "
                          "fall, and the snow is that of erf.radiation.seb_snow_initial_swe_kg_m2\n";
    }
    if (m_unused[lev] == 0) { fill_accumulation(*m_accum_start[lev], sources, klo); }
    m_started[lev] = 1;
}

void
ForceRestorePrecipRate::end_step (int lev, const SurfacePrecipAccumulationSources& sources,
                                  int klo, Real dt)
{
    if (!defined(lev)) { return; }
    if (m_started[lev] == 0) { m_rate[lev]->setVal(0.0); return; }
    fill_accumulation(*m_accum_now[lev], sources, klo);
    // The time since the stored start: this step, plus any step whose rate no step took.
    m_elapsed[lev] += dt;
    const Real elapsed = m_elapsed[lev];
    const int ncomp = m_rate[lev]->nComp();
    for (MFIter mfi(*m_rate[lev], TilingIfNotGPU()); mfi.isValid(); ++mfi) {
        const Box& bx = mfi.tilebox();
        Array4<Real> r = m_rate[lev]->array(mfi);
        Array4<const Real> a0 = m_accum_start[lev]->const_array(mfi);
        Array4<const Real> a1 = m_accum_now[lev]->const_array(mfi);
        ParallelFor(bx, ncomp, [=] AMREX_GPU_DEVICE (int i, int j, int k, int n) noexcept
        {
            // The change of each part over the step, not below zero: the bound is on the
            // change, never on the accumulation itself.
            r(i, j, k, n) = force_restore_rain_rate(a1(i, j, k, n), a0(i, j, k, n), elapsed);
        });
    }
    m_started[lev] = 0;
    m_unused[lev] = 1;
}

ForceRestorePrecipInput
ForceRestorePrecipRate::input (int lev, int comp, Real dt) const noexcept
{
    ForceRestorePrecipInput in;
    if (!defined(lev)) { return in; }
    in.rate = m_rate[lev].get();
    in.comp = comp;
    // The rate is a mean over m_elapsed, which exceeds dt when an earlier step's
    // precipitation was carried: the step takes the whole of it. A rate a step already took
    // (no end_step since) gives nothing.
    in.scale = (m_unused[lev] != 0 && dt > Real(0.0)) ? m_elapsed[lev] / dt : Real(0.0);
    return in;
}

void
ForceRestorePrecipRate::mark_used (int lev) noexcept
{
    if (!defined(lev)) { return; }
    m_unused[lev] = 0;
    m_elapsed[lev] = 0.0;
}
