#include <cmath>
#include <functional>
#include <limits>

#include <AMReX_BoxIterator.H>
#include <AMReX_ParmParse.H>
#include <gtest/gtest.h>

#include <ERF_EOS.H>
#include <ERF_MicrophysicsUtils.H>
#include <ERF_SuperDropletsMoist.H>

namespace {

using namespace amrex;

class ThermoModel : public SuperDropletsMoist
{
public:
    ThermoModel ()
    {
        m_super_droplets = nullptr;
        m_name = "sdm_thermodynamics_test";
        m_species = {Species::Name::H2O, Species::Name::ice, Species::Name::water};
        m_num_species = 3;
        m_idx_w = 0;
        m_idx_i = 1;
        m_istart_sp = 2;
        m_num_nonmoist_sp = 1;
        m_qstate_nonmoist_size = 2;
        m_with_ice = true;
        m_aerosols = {Species::Name::NaCl};
        m_num_aerosols = 1;
        m_init_phase_change = false;
        m_num_substeps_phase_change = 2;
        m_diagnostics_iter = 0;
    }

    ~ThermoModel () override { delete m_super_droplets; }

    MultiFab& field (int n, int lev = 0) { return *m_mic_fab_vars[lev][n]; }
    void restored_particles () { m_restarted_particles = true; }

    std::function<void()> before_substep;
    std::function<void(int, bool)> process_observer;
    void phaseChange_SL_w (const Real& dt, const Vector<std::unique_ptr<MultiFab>>& z,
                           int lev, const iMultiFab& mask) override
    {
        if (before_substep) { before_substep(); }
        if (process_observer) { process_observer(0, true); }
        SuperDropletsMoist::phaseChange_SL_w(dt, z, lev, mask);
        if (process_observer) { process_observer(0, false); }
    }

    void phaseChange_LV_w (const Real& dt, const Vector<std::unique_ptr<MultiFab>>& z,
                           int lev, const iMultiFab& mask) override
    {
        if (process_observer) { process_observer(1, true); }
        SuperDropletsMoist::phaseChange_LV_w(dt, z, lev, mask);
        if (process_observer) { process_observer(1, false); }
    }

    void phaseChange_SV_i (const Real& dt, const Vector<std::unique_ptr<MultiFab>>& z,
                           int lev, const iMultiFab& mask) override
    {
        if (process_observer) { process_observer(2, true); }
        SuperDropletsMoist::phaseChange_SV_i(dt, z, lev, mask);
        if (process_observer) { process_observer(2, false); }
    }

    void phaseChange_LV_s (int species, const Real& dt,
                           const Vector<std::unique_ptr<MultiFab>>& z,
                           int lev, const iMultiFab& mask) override
    {
        if (process_observer) { process_observer(3, true); }
        SuperDropletsMoist::phaseChange_LV_s(species, dt, z, lev, mask);
        if (process_observer) { process_observer(3, false); }
    }
};

struct ThermoFixture {
    Geometry geom;
    BoxArray boxes;
    DistributionMapping dm;
    MultiFab cons;
    MultiFab base;
    SolverChoice sc;
    ThermoModel model;
    Vector<std::unique_ptr<MultiFab>> z;

    explicit ThermoFixture (Vector<int> modes = {1, 0}, Real cp = Cp_d,
                            bool split = false, bool particles = false)
        : geom(Box(IntVect(0), IntVect(3)),
               RealBox(Real(0), Real(0), Real(0), Real(4), Real(4), Real(4)),
               CoordSys::cartesian, {0, 0, 0}),
          boxes(geom.Domain()),
          dm(make_boxes(split)),
          cons(boxes, dm, RhoQ1_comp + 8, 1),
          base(boxes, dm, BaseState::num_comps, 1)
    {
        ParmParse pp("sdm_thermodynamics_test");
        pp.add("num_initializations", particles ? 1 : 0);
        pp.add("num_injections", 0);
        pp.add("initial_number_density", Real(1.e7));
        pp.add("initial_particles_per_cell", 1);
        pp.add("multiplicity_type", "constant");
        pp.add("initial_species_mean_mass_H2O", Real(1.e-12));
        pp.add("initial_species_mean_mass_water", Real(2.e-12));
        pp.add("initial_aerosol_mean_mass_NaCl", Real(1.e-19));
        pp.add("place_randomly_in_cells", false);
        sc.anelastic = std::move(modes);
        sc.c_p = cp;
        sc.rdOcp = R_d / cp;
        model.Define(sc);
        initialize_state(particles);
    }

    // CUDA extended device lambdas need an addressable, non-constructor parent.
    void initialize_state (bool particles)
    {
        cons.setVal(Real(0));
        base.setVal(Real(-12345));
        for (MFIter mfi(cons); mfi.isValid(); ++mfi) {
            const auto s = cons.array(mfi);
            const auto b = base.array(mfi);
            ParallelFor(mfi.fabbox(), [=] AMREX_GPU_DEVICE (int i, int j, int k) noexcept {
                const bool left = i < 2;
                const Real rho = left ? Real(1.1) : Real(0.9);
                s(i,j,k,Rho_comp) = rho;
                s(i,j,k,RhoTheta_comp) = rho * (left ? Real(303) : Real(287));
                s(i,j,k,RhoQ1_comp) = rho * (left ? Real(0.003) : Real(0.007));
                s(i,j,k,RhoQ7_comp) = rho * Real(0.004);
                b(i,j,k,BaseState::p0_comp) = left ? Real(90000) : Real(82000);
            });
        }
        std::unique_ptr<MultiFab> unused;
        model.Init(cons, boxes, geom, Real(0.001), unused, unused);
        if (particles) {
            BoxArray nodes(boxes);
            nodes.surroundingNodes();
            z.push_back(std::make_unique<MultiFab>(nodes, dm, 1, 2));
            for (MFIter mfi(*z[0]); mfi.isValid(); ++mfi) {
                const auto height = z[0]->array(mfi);
                ParallelFor(mfi.fabbox(), [=] AMREX_GPU_DEVICE (int i, int j, int k) noexcept {
                    height(i,j,k) = Real(k);
                });
            }
            model.InitParticles(0, z[0]);
            auto& pc = *static_cast<SuperDropletPC*>(model.getParticleContainer());
            for (SuperDropletPC::ParIterType pti(pc, 0); pti.isValid(); ++pti) {
                const auto data = pti.GetParticleTile().getParticleTileData();
                ParallelFor(pti.numParticles(), [=] AMREX_GPU_DEVICE (int n) noexcept {
                    data.m_runtime_rdata[ridx_ice_Tfz(1, 3)][n] = ParticleReal(273.15);
                });
            }
            update(0, sc.anelastic[0] ? &base : nullptr);
            model.FinishInit(0, cons, z);
        }
    }

    const BoxArray& make_boxes (bool split)
    {
        if (split) { boxes.maxSize(IntVect(2, 4, 4)); }
        return boxes;
    }

    void update (int lev = 0, const MultiFab* pressure_base = nullptr)
    {
        model.SetCurrentLevel(lev);
        NullMoist& interface = model;
        interface.Update_Micro_Vars(cons, pressure_base);
        Gpu::streamSynchronize();
    }
};

// Host mirrors also inspect physical ghosts without host-dereferencing device arrays.
template <typename F>
void check_field (const MultiFab& mf, F&& expected, Real operations = Real(64))
{
    const Real eps = std::numeric_limits<Real>::epsilon();
    Gpu::streamSynchronize();
    for (MFIter mfi(mf); mfi.isValid(); ++mfi) {
        FArrayBox host(mfi.fabbox(), 1, The_Pinned_Arena());
        Gpu::dtoh_memcpy(host.dataPtr(), mf[mfi].dataPtr(), mf[mfi].size() * sizeof(Real));
        const auto values = host.const_array();
        for (BoxIterator cell(mfi.fabbox()); cell.ok(); ++cell) {
            const auto iv = cell();
            const Real want = expected(iv);
            EXPECT_NEAR(values(iv), want, operations * eps * std::max(Real(1), std::abs(want)))
                << iv;
        }
    }
}

void check_copy_in (ThermoFixture& f, int lev, bool anelastic)
{
    check_field(f.model.field(MicVar_SD::pressure, lev), [=](const IntVect& iv) {
        const bool left = iv[0] < 2;
        if (anelastic) { return left ? Real(90000) : Real(82000); }
        const Real rho = left ? Real(1.1) : Real(0.9);
        return getPgivenRTh(rho * (left ? Real(303) : Real(287)),
                           left ? Real(0.003) : Real(0.007));
    });
    const Real kappa = f.sc.rdOcp;
    check_field(f.model.field(MicVar_SD::temperature, lev), [=](const IntVect& iv) {
        const bool left = iv[0] < 2;
        const Real theta = left ? Real(303) : Real(287);
        if (anelastic) {
            return theta * std::pow((left ? Real(90000) : Real(82000)) / Real(100000), kappa);
        }
        const Real rho = left ? Real(1.1) : Real(0.9);
        return getTgivenRandRTh(rho, rho * theta, left ? Real(0.003) : Real(0.007));
    });
}

void copy_in_test ()
{
    for (Real cp : {Real(Cp_d), Real(1200)}) {
        for (bool split : {false, true}) {
            ThermoFixture f({1, 0}, cp, split);
            f.update(0, &f.base);
            check_copy_in(f, 0, true);
        }
    }
}

// Motivation: Anelastic copy-in must use cellwise base pressure, including
// physical ghosts. Varying c_p and box layouts detects EOS fallback, a fixed
// exponent, or incomplete pressure/temperature fills.
TEST(SDMAnelasticThermodynamics, CopyInUsesHydrostaticPressureAndConfiguredCp)
{
    copy_in_test();
}

void compressible_test ()
{
    ThermoFixture f({0, 1}, Real(1200));
    f.update();
    check_copy_in(f, 0, false);
    // An incompatible, zero-ghost dummy must be completely ignored in this mode.
    MultiFab dummy(BoxArray(Box(IntVect(0), IntVect(0))), DistributionMapping{Vector<int>{0}}, 1, 0);
    dummy.setVal(Real(-1));
    f.update(0, &dummy);
    check_copy_in(f, 0, false);
}

// Motivation: Anelastic support must not change compressible EOS arithmetic.
// An unusable base-state argument must be ignored on a compressible level,
// even when a nondefault c_p is configured.
TEST(SDMAnelasticThermodynamics, CompressibleIgnoresBasePointerAndCustomCp)
{
    compressible_test();
}

void mixed_level_test ()
{
    for (Vector<int> modes : {Vector<int>{1, 0}, Vector<int>{0, 1}}) {
        ThermoFixture f(modes);
        for (int lev : {1, 0, 1}) {
            f.model.InitLevel(lev, f.cons);
            f.update(lev, &f.base);
            check_copy_in(f, lev, modes[lev] != 0);
        }
    }
}

// Motivation: One SDM instance serves all AMR levels. Both mode orders and
// repeated level switches must select the current level's pressure source,
// rather than retaining the first or most recently used level's mode.
TEST(SDMAnelasticThermodynamics, SharedModelSelectsEachLevelIndependently)
{
    mixed_level_test();
}

class ThrowOnAbort
{
public:
    ThrowOnAbort () : saved(system::throw_exception) { system::throw_exception = true; }
    ~ThrowOnAbort () { system::throw_exception = saved; }
private:
    bool saved;
};

// Motivation: Missing pressure, incompatible layouts or ghosts, unconfigured
// levels, and invalid c_p must be rejected before source kernels can consume
// inconsistent thermodynamic inputs.
TEST(SDMAnelasticThermodynamics, MissingOrIncompatibleBaseFailsBeforeCopyIn)
{
    ThermoFixture f;
    ThrowOnAbort throwing;
    EXPECT_THROW(f.update(), RuntimeError);
    EXPECT_THROW(f.model.Copy_State_to_Micro(f.cons), RuntimeError);
    MultiFab short_ghosts(f.boxes, f.dm, BaseState::num_comps, 0);
    EXPECT_THROW(f.update(0, &short_ghosts), RuntimeError);
    MultiFab short_components(f.boxes, f.dm, BaseState::p0_comp, 1);
    EXPECT_THROW(f.update(0, &short_components), RuntimeError);
    BoxArray different(f.boxes);
    different.maxSize(2);
    MultiFab wrong_layout(different, DistributionMapping(different), BaseState::num_comps, 1);
    EXPECT_THROW(f.update(0, &wrong_layout), RuntimeError);
    EXPECT_THROW(f.update(2, &f.base), RuntimeError);
    f.sc.c_p = Real(0);
    EXPECT_THROW(f.model.Define(f.sc), RuntimeError);
}

void reallocation_test ()
{
    ThermoFixture f;
    f.update(0, &f.base);
    MultiFab cons(f.boxes, f.dm, f.cons.nComp(), 2);
    MultiFab base(f.boxes, f.dm, BaseState::num_comps, 2);
    cons.setVal(Real(0));
    cons.setVal(Real(1), Rho_comp, 1, 2);
    cons.setVal(Real(300), RhoTheta_comp, 1, 2);
    base.setVal(Real(-1));
    base.setVal(Real(75000), BaseState::p0_comp, 1, 2);
    NullMoist& interface = f.model;
    interface.Update_Micro_Vars(cons, &base);
    EXPECT_EQ(f.model.field(MicVar_SD::pressure).nGrowVect(), IntVect(2));
    check_field(f.model.field(MicVar_SD::pressure), [](const IntVect&) { return Real(75000); });
    const Real temperature = Real(300) * std::pow(Real(0.75), f.sc.rdOcp);
    check_field(f.model.field(MicVar_SD::temperature), [=](const IntVect&) { return temperature; });
}

// Motivation: Reallocated AMR state must receive current base pressure and
// temperature over its new ghost extent, not stale values from an earlier
// allocation or a retained base-state pointer.
TEST(SDMAnelasticThermodynamics, CopyInReseedsPressureAfterArrayReallocation)
{
    reallocation_test();
}

// Motivation: Temperature is lagged within a substep while theta changes.
// Sequential sources with both signs expose the spurious feedback from theta/T
// on anelastic levels while checking the compressible rule independently.
TEST(SDMAnelasticThermodynamics, SequentialLatentSourcesUseFixedExnerWithBothSigns)
{
    const Real cp = Real(1200);
    const Real kappa = Real(287) / cp;
    const Real pb = Real(82000);
    const Real theta0 = Real(310);
    const Real pi = std::pow(pb / Real(100000), kappa);
    const Real stored_temperature = theta0 * pi;
    const Real tolerance = Real(32) * std::numeric_limits<Real>::epsilon() * theta0;
    for (Real sign1 : {Real(-1), Real(1)}) {
        for (Real sign2 : {Real(-1), Real(1)}) {
            Real theta = theta0;
            Real compressible = theta0;
            Real old_rule = theta0;
            Real heat = Real(0);
            const Real latent[] = {Real(2.5e6), Real(3.34e5)};
            const Real dq[] = {sign1 * Real(0.004), sign2 * Real(0.006)};
            for (int process = 0; process < 2; ++process) {
                const Real dT = latent[process] / cp * dq[process];
                theta += sdm_theta_per_temperature(theta, stored_temperature, pb, kappa, true) * dT;
                heat += dT;
                old_rule += (old_rule / stored_temperature) * dT;
                compressible += sdm_theta_per_temperature(
                    compressible, stored_temperature, Real(-1), kappa, false) * dT;
                EXPECT_NEAR(theta, theta0 + heat / pi, tolerance);
                EXPECT_NEAR((theta - theta0) * pi, heat, tolerance);
                EXPECT_EQ(compressible, old_rule);
            }
            EXPECT_GT(std::abs(old_rule - theta), Real(8) * tolerance);
        }
    }
}

void check_current_thermo (ThermoFixture& f)
{
    auto& temperature = f.model.field(MicVar_SD::temperature);
    MultiFab error(f.boxes, f.dm, 1, temperature.nGrowVect());
    const Real kappa = f.sc.rdOcp;
    for (MFIter mfi(error); mfi.isValid(); ++mfi) {
        const auto e = error.array(mfi);
        const auto t = temperature.const_array(mfi);
        const auto theta = f.model.field(MicVar_SD::theta).const_array(mfi);
        const auto base = f.base.const_array(mfi);
        ParallelFor(mfi.fabbox(), [=] AMREX_GPU_DEVICE (int i, int j, int k) noexcept {
            e(i,j,k) = t(i,j,k) -
                theta(i,j,k) * std::pow(base(i,j,k,BaseState::p0_comp) / Real(100000), kappa);
        });
    }
    EXPECT_LE(error.norm0(0, 1, error.nGrowVect()),
              Real(128) * std::numeric_limits<Real>::epsilon() * temperature.norm0());
    check_field(f.model.field(MicVar_SD::pressure), [](const IntVect& iv) {
        return iv[0] < 2 ? Real(90000) : Real(82000);
    });
}

void check_saturation (ThermoFixture& f)
{
    for (int species = 0; species < 3; ++species) {
        const int sr_index = species == 0 ? MicVar_SD::rh_w :
                             species == 1 ? MicVar_SD::rh_i : s_sr_idx(2, 2);
        const int q_index = species == 2 ? s_qv_idx(2, 2) : MicVar_SD::q_v;
        auto& sr = f.model.field(sr_index);
        MultiFab error(f.boxes, f.dm, 1, sr.nGrowVect());
        for (MFIter mfi(error); mfi.isValid(); ++mfi) {
            const auto e = error.array(mfi);
            const auto actual = sr.const_array(mfi);
            const auto qv = f.model.field(q_index).const_array(mfi);
            const auto t = f.model.field(MicVar_SD::temperature).const_array(mfi);
            const auto base = f.base.const_array(mfi);
            ParallelFor(mfi.fabbox(), [=] AMREX_GPU_DEVICE (int i, int j, int k) noexcept {
                const Real es_hpa = species == 1 ? erf_esati(t(i,j,k)) : erf_esatw(t(i,j,k));
                const Real qs = (R_d / R_v) * es_hpa /
                    amrex::max(es_hpa, base(i,j,k,BaseState::p0_comp) / Real(100) - es_hpa);
                const Real expected = qs > Real(0) ? qv(i,j,k) / qs : Real(0);
                e(i,j,k) = (actual(i,j,k) - expected) / amrex::max(Real(1), expected);
            });
        }
        EXPECT_LE(error.norm0(0, 1, error.nGrowVect()),
                  Real(128) * std::numeric_limits<Real>::epsilon());
    }
}

// Motivation: Correct temperature alone can hide an incorrect saturation
// pressure source. Independent water, ice, and additional-species expectations
// detect EOS fallback or a Pa/hPa conversion error.
TEST(SDMAnelasticThermodynamics, SaturationUsesBasePressureInPaForWaterIceAndSpecies)
{
    ThermoFixture f;
    f.update(0, &f.base);
    check_saturation(f);
}

Real reservoir (ThermoModel& model, bool other_species)
{
    const int first = other_species ? s_qv_idx(2, 2) : MicVar_SD::q_v;
    const int last = other_species ? s_qc_idx(2, 2) : MicVar_SD::q_g;
    Real sum = Real(0);
    for (int n = first; n <= last; ++n) {
        sum += MultiFab::Dot(model.field(MicVar_SD::rho), 0, model.field(n), 0, 1, 0);
    }
    return sum;
}

void check_particles (ThermoModel& model)
{
    auto& pc = *static_cast<SuperDropletPC*>(model.getParticleContainer());
    const int offset = SDPCDefn::SuperDropletsRealIdx::ncomps;
    for (SuperDropletPC::ParIterType pti(pc, 0); pti.isValid(); ++pti) {
        const auto& soa = pti.GetStructOfArrays();
        for (int component : {int(SDPCDefn::SuperDropletsRealIdx::mass),
                              offset + SDPCDefn::SuperDropletsRealIdxSoA_RT::radius,
                              offset + SDPCDefn::SuperDropletsRealIdxSoA_RT::multiplicity,
                              idx_s(0, 1, 3), idx_s(1, 1, 3), idx_s(2, 1, 3)}) {
            const auto& device = soa.GetRealData(component);
            Vector<ParticleReal> values(device.size());
            Gpu::copy(Gpu::deviceToHost, device.begin(), device.end(), values.begin());
            for (ParticleReal value : values) {
                EXPECT_TRUE(std::isfinite(value));
                EXPECT_GE(value, ParticleReal(0));
            }
        }
    }
}

void particle_substep_test ()
{
    for (bool split : {false, true}) {
        ThermoFixture f({1}, Real(1200), split, true);
        const Real water = reservoir(f.model, false);
        const Real species = reservoir(f.model, true);
        const Real eps = std::numeric_limits<Real>::epsilon();
        const Real expected_water = Real(32) * (Real(1.1) * Real(0.003) + Real(0.9) * Real(0.007))
                                    + Real(64) * Real(1.e7) * Real(1.e-12);
        const Real expected_species = Real(32) * (Real(1.1) + Real(0.9)) * Real(0.004)
                                      + Real(64) * Real(1.e7) * Real(2.e-12);
        EXPECT_NEAR(water, expected_water, Real(128) * eps * expected_water);
        EXPECT_NEAR(species, expected_species, Real(128) * eps * expected_species);
        const Real initial_theta = f.model.field(MicVar_SD::theta).sum();
        const Real initial_vapor = f.model.field(MicVar_SD::q_v).sum();
        int visited = 0;
        f.model.before_substep = [&]() {
            ++visited;
            check_current_thermo(f);
            check_saturation(f);
            check_particles(f.model);
        };
        f.model.phaseChange(Real(0.02), f.z, 0);
        check_particles(f.model);
        EXPECT_EQ(visited, 2);
        check_current_thermo(f);
        check_saturation(f);
        const Real tol = Real(512) * std::numeric_limits<Real>::epsilon();
        EXPECT_NEAR(reservoir(f.model, false), water, tol * water);
        EXPECT_NEAR(reservoir(f.model, true), species, tol * species);
        EXPECT_GT(std::abs(f.model.field(MicVar_SD::theta).sum() - initial_theta), Real(1.e-4));
        EXPECT_GT(std::abs(f.model.field(MicVar_SD::q_v).sum() - initial_vapor), Real(1.e-8));
        for (int n = MicVar_SD::q_v; n <= MicVar_SD::q_g; ++n) {
            EXPECT_GE(f.model.field(n).min(0), Real(0));
            EXPECT_FALSE(f.model.field(n).contains_nan());
        }
    }
}

// Motivation: Real particle exchange must preserve dry-air-referenced reservoirs
// and held base pressure across active substeps and split boxes. Inspecting the
// next substep also detects an erroneous EOS pressure refresh.
TEST(SDMAnelasticThermodynamics, ActiveParticleSubstepsKeepPressureAndConserveDryAirReferencedMass)
{
    particle_substep_test();
}

// Motivation: Checkpointed condensate can be a fine-level average that differs
// from local particle deposition. Restart initialization must preserve it
// rather than perturb the moisture-dependent resumed dynamics.
TEST(SDMAnelasticThermodynamics, RestartFinishInitPreservesCheckpointedCondensate)
{
    ThermoFixture f({1}, Real(1200), true, true);
    f.cons.setVal(Real(2.e-5), RhoQ2_comp, 1, 0);
    f.cons.setVal(Real(3.e-5), RhoQ8_comp, 1, 0);
    MultiFab expected(f.boxes, f.dm, f.cons.nComp(), 0);
    MultiFab::Copy(expected, f.cons, 0, 0, f.cons.nComp(), 0);
    f.model.restored_particles();
    f.update(0, &f.base);
    f.model.FinishInit(0, f.cons, f.z);
    MultiFab::Subtract(expected, f.cons, 0, 0, f.cons.nComp(), 0);
    EXPECT_EQ(expected.norm0(0, f.cons.nComp(), IntVect(0)), Real(0));
}

void real_latent_sites_test ()
{
    ThermoFixture f({1}, Real(1200), true, true);
    for (MFIter mfi(f.cons); mfi.isValid(); ++mfi) {
        const auto state = f.cons.array(mfi);
        ParallelFor(mfi.fabbox(), [=] AMREX_GPU_DEVICE (int i, int j, int k) noexcept {
            state(i,j,k,RhoTheta_comp) = state(i,j,k,Rho_comp) * (i < 2 ? Real(303) : Real(280));
        });
    }
    f.update(0, &f.base);
    MultiFab old_theta(f.boxes, f.dm, 1, 0);
    MultiFab old_q(f.boxes, f.dm, 1, 0);
    MultiFab delta(f.boxes, f.dm, 1, 0);
    MultiFab error(f.boxes, f.dm, 1, 0);
    Real exchanged[4] = {0, 0, 0, 0};
    f.model.process_observer = [&](int process, bool before) {
        const int qindex = process == 0 ? MicVar_SD::q_c :
                           process == 3 ? s_qv_idx(2, 2) : MicVar_SD::q_v;
        if (before) {
            MultiFab::Copy(old_theta, f.model.field(MicVar_SD::theta), 0, 0, 1, 0);
            MultiFab::Copy(old_q, f.model.field(qindex), 0, 0, 1, 0);
            if (process == 0) { MultiFab::Add(old_q, f.model.field(MicVar_SD::q_r), 0, 0, 1, 0); }
            return;
        }
        const Real latent = process == 0 ? Real(lfus) :
                            process == 2 ? Real(2.8342e6) : Real(L_v);
        const Real cp = f.sc.c_p;
        const Real kappa = Real(287) / cp;
        for (MFIter mfi(error); mfi.isValid(); ++mfi) {
            const auto theta0 = old_theta.const_array(mfi);
            const auto q0 = old_q.const_array(mfi);
            const auto theta = f.model.field(MicVar_SD::theta).const_array(mfi);
            const auto q = f.model.field(qindex).const_array(mfi);
            const auto rain = f.model.field(MicVar_SD::q_r).const_array(mfi);
            const auto base = f.base.const_array(mfi);
            const auto e = error.array(mfi);
            const auto d = delta.array(mfi);
            ParallelFor(mfi.validbox(), [=] AMREX_GPU_DEVICE (int i, int j, int k) noexcept {
                const Real dq = q0(i,j,k) - q(i,j,k) - (process == 0 ? rain(i,j,k) : Real(0));
                const Real inv_pi = std::pow(Real(100000) / base(i,j,k,BaseState::p0_comp), kappa);
                e(i,j,k) = theta(i,j,k) - theta0(i,j,k) - inv_pi * (latent / cp) * dq;
                d(i,j,k) = dq;
            });
        }
        EXPECT_LE(error.norm0(), Real(128) * std::numeric_limits<Real>::epsilon() * old_theta.norm0());
        exchanged[process] += delta.norm0();
    };
    f.model.phaseChange(Real(0.02), f.z, 0);
    check_particles(f.model);
    check_current_thermo(f);
    for (int process = 0; process < 4; ++process) {
        EXPECT_GT(exchanged[process], Real(0)) << "inactive latent site " << process;
    }
}

// Motivation: A correct helper does not prove that every particle process uses
// it. Measuring actual exchange and theta changes at all four active latent
// sites detects a missed call site or reuse of the evolving theta/T factor.
TEST(SDMAnelasticThermodynamics, AllFourParticleLatentSitesUseFixedExner)
{
    real_latent_sites_test();
}

} // namespace
