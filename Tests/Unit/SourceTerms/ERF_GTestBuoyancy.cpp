#include <algorithm>
#include <cmath>
#include <limits>
#include <memory>
#include <vector>

#include <AMReX_BoxArray.H>
#include <AMReX_DistributionMapping.H>
#include <AMReX_Gpu.H>
#include <AMReX_GpuContainers.H>
#include <AMReX_MultiFab.H>

#include <gtest/gtest.h>

#include "ERF_Constants.H"
#include "ERF_IndexDefines.H"
#include "ERF_SrcHeaders.H"
#include "ERF_SlowRhsPreUtils.H"

namespace {

constexpr int kNz = 8;
constexpr amrex::Real kDz = amrex::Real(1000.0);
constexpr amrex::Real kT0 = amrex::Real(300.0);
constexpr amrex::Real kTopSentinel = amrex::Real(12345.0);

#ifdef AMREX_USE_FLOAT
constexpr amrex::Real kAbsTol = amrex::Real(2.0e-4);
constexpr amrex::Real kRelTol = amrex::Real(3.0e-5);
#else
constexpr amrex::Real kAbsTol = amrex::Real(5.0e-13);
constexpr amrex::Real kRelTol = amrex::Real(3.0e-12);
#endif

struct StateParameters {
    amrex::Real qv0 = amrex::Real(0.0);
    amrex::Real theta_factor = amrex::Real(1.0);
    amrex::Real pressure_factor = amrex::Real(1.0);
    amrex::Real vapor_anomaly = amrex::Real(0.0);
    amrex::Real condensate = amrex::Real(0.0);
};

class BuoyancyFixture {
public:
    BuoyancyFixture ()
        : domain(amrex::IntVect(0, 0, 0), amrex::IntVect(0, 0, kNz-1)),
          grids(domain), mapping(grids),
          real_box({amrex::Real(0.0), amrex::Real(0.0), amrex::Real(0.0)},
                   {amrex::Real(1.0), amrex::Real(1.0), amrex::Real(kNz)*kDz}),
          periodic{1, 1, 0},
          geom(domain, &real_box, amrex::CoordSys::cartesian, periodic.data()),
          base_state(grids, mapping, BaseState::num_comps, 1),
          primitive(grids, mapping, NPRIMVAR_max, 1),
          total_water(grids, mapping, 1, 1),
          buoyancy(amrex::convert(grids, amrex::IntVect(0, 0, 1)), mapping, 1, 0),
          conserved(IntVars::NumTypes)
    {
        conserved[IntVars::cons].define(grids, mapping, NVAR_max, 1);

        solver.terrain_type = TerrainType::None;
        solver.gravity = CONST_GRAV;
        solver.rdOcp = RdoCp;
        solver.moisture_type = MoistureType::None;
        solver.buoyancy_type.resize(1);
        eb_factory = std::make_unique<eb_>();
    }

    void fill_state (const StateParameters& parameters, bool fixed_density = false)
    {
        auto& cons = conserved[IntVars::cons];
        for (amrex::MFIter mfi(cons); mfi.isValid(); ++mfi) {
            const auto bx = mfi.growntilebox();
            const auto state = cons.array(mfi);
            const auto prim = primitive.array(mfi);
            const auto qt = total_water.array(mfi);
            const auto base = base_state.array(mfi);
            const auto qv0 = parameters.qv0;
            const auto theta_factor = parameters.theta_factor;
            const auto pressure_factor = parameters.pressure_factor;
            const auto vapor_anomaly = parameters.vapor_anomaly;
            const auto condensate = parameters.condensate;
            const auto use_fixed_density = fixed_density;

            amrex::ParallelFor(bx, [=] AMREX_GPU_DEVICE (int i, int j, int k) noexcept {
                const amrex::Real z = (static_cast<amrex::Real>(k) + amrex::Real(0.5))*kDz;
                const amrex::Real moist_hydrostatic_factor = (amrex::Real(1.0) + qv0) /
                    (amrex::Real(1.0) + (R_v/R_d)*qv0);
                const amrex::Real scale_height = R_d*kT0 /
                    (CONST_GRAV*moist_hydrostatic_factor);
                const amrex::Real p0 = p_0*std::exp(-z/scale_height);
                const amrex::Real rho0 = p0 /
                    (R_d*kT0*(amrex::Real(1.0) + (R_v/R_d)*qv0));
                const amrex::Real theta0 = kT0*std::pow(p_0/p0, RdoCp);
                const amrex::Real qv = qv0 + vapor_anomaly;
                const amrex::Real theta = theta_factor*theta0;
                const amrex::Real pressure = pressure_factor*p0;
                const amrex::Real rho = use_fixed_density
                    ? rho0 : getRhogivenThetaPress(theta, pressure, RdoCp, qv);
                const amrex::Real qt_value = qv + condensate;

                base(i,j,k,BaseState::r0_comp) = rho0;
                base(i,j,k,BaseState::p0_comp) = p0;
                base(i,j,k,BaseState::pi0_comp) = amrex::Real(1.0);
                base(i,j,k,BaseState::th0_comp) = theta0;
                base(i,j,k,BaseState::qv0_comp) = qv0;

                for (int n = 0; n < NVAR_max; ++n) state(i,j,k,n) = amrex::Real(0.0);
                for (int n = 0; n < NPRIMVAR_max; ++n) prim(i,j,k,n) = amrex::Real(0.0);
                state(i,j,k,Rho_comp) = rho;
                state(i,j,k,RhoTheta_comp) = rho*theta;
                state(i,j,k,RhoQ1_comp) = rho*qv;
                state(i,j,k,RhoQ2_comp) = rho*condensate;
                prim(i,j,k,PrimTheta_comp) = theta;
                prim(i,j,k,PrimQ1_comp) = qv;
                prim(i,j,k,PrimQ2_comp) = condensate;
                qt(i,j,k) = qt_value;
            });
        }
        amrex::Gpu::streamSynchronize();
    }

    std::vector<amrex::Real> run (int selector,
                                  int anelastic,
                                  MoistureType moisture,
                                  amrex::Real gravity = CONST_GRAV)
    {
        solver.gravity = gravity;
        solver.moisture_type = moisture;
        solver.buoyancy_type[0] = selector;
        buoyancy.setVal(kTopSentinel);
        make_buoyancy(0, conserved, primitive, total_water, buoyancy, geom,
                      solver, base_state, moisture == MoistureType::None ? 0 : 11,
                      *eb_factory, anelastic);
        amrex::Gpu::streamSynchronize();

        amrex::Gpu::DeviceVector<amrex::Real> device_values(kNz+1);
        amrex::Gpu::HostVector<amrex::Real> host_values(kNz+1);
        auto* output = device_values.data();
        for (amrex::MFIter mfi(buoyancy); mfi.isValid(); ++mfi) {
            const auto values = buoyancy.const_array(mfi);
            const amrex::Box faces(amrex::IntVect(0, 0, 0), amrex::IntVect(0, 0, kNz));
            amrex::ParallelFor(faces, [=] AMREX_GPU_DEVICE (int i, int j, int k) noexcept {
                output[k] = values(i,j,k);
            });
        }
        amrex::Gpu::streamSynchronize();
        amrex::Gpu::copy(amrex::Gpu::deviceToHost, device_values.begin(), device_values.end(),
                         host_values.begin());
        return std::vector<amrex::Real>(host_values.begin(), host_values.end());
    }

    std::vector<amrex::Real> density_and_base_density () const
    {
        const auto& cons = conserved[IntVars::cons];
        amrex::Gpu::DeviceVector<amrex::Real> device_values(2*kNz);
        amrex::Gpu::HostVector<amrex::Real> host_values(2*kNz);
        auto* output = device_values.data();
        for (amrex::MFIter mfi(cons); mfi.isValid(); ++mfi) {
            const auto state = cons.const_array(mfi);
            const auto base = base_state.const_array(mfi);
            amrex::ParallelFor(domain, [=] AMREX_GPU_DEVICE (int i, int j, int k) noexcept {
                output[2*k] = state(i,j,k,Rho_comp);
                output[2*k+1] = base(i,j,k,BaseState::r0_comp);
            });
        }
        amrex::Gpu::streamSynchronize();
        amrex::Gpu::copy(amrex::Gpu::deviceToHost, device_values.begin(), device_values.end(),
                         host_values.begin());
        return std::vector<amrex::Real>(host_values.begin(), host_values.end());
    }

    amrex::Real run_z_pressure_buoyancy_rhs (int face, amrex::Real qt_lo,
                                             amrex::Real qt_hi, amrex::Real gpz,
                                             amrex::Real abl_pressure_grad_z,
                                             amrex::Real buoyancy_force,
                                             bool use_moisture)
    {
        total_water.setVal(amrex::Real(0.0));
        const amrex::Box adjacent_cells(amrex::IntVect(0, 0, face-1),
                                        amrex::IntVect(0, 0, face));
        for (amrex::MFIter mfi(total_water); mfi.isValid(); ++mfi) {
            const auto qt = total_water.array(mfi);
            amrex::ParallelFor(adjacent_cells, [=] AMREX_GPU_DEVICE (int i, int j, int k) noexcept {
                qt(i,j,k) = (k == face-1) ? qt_lo : qt_hi;
            });
        }
        amrex::Gpu::streamSynchronize();

        amrex::Gpu::DeviceVector<amrex::Real> device_value(1);
        amrex::Gpu::HostVector<amrex::Real> host_value(1);
        auto* output = device_value.data();
        const amrex::Box face_box(amrex::IntVect(0, 0, face),
                                  amrex::IntVect(0, 0, face));
        for (amrex::MFIter mfi(total_water); mfi.isValid(); ++mfi) {
            const auto qt = total_water.const_array(mfi);
            amrex::ParallelFor(face_box, [=] AMREX_GPU_DEVICE (int i, int j, int k) noexcept {
                output[0] = slow_rhs_pre_z_pressure_buoyancy(
                    qt, i, j, k, use_moisture, gpz, abl_pressure_grad_z, buoyancy_force);
            });
        }
        amrex::Gpu::streamSynchronize();
        amrex::Gpu::copy(amrex::Gpu::deviceToHost, device_value.begin(), device_value.end(),
                         host_value.begin());
        return host_value[0];
    }

    amrex::Real base_density (int k, amrex::Real qv0) const
    {
        const amrex::Real z = (static_cast<amrex::Real>(k) + amrex::Real(0.5))*kDz;
        const amrex::Real moist_hydrostatic_factor = (amrex::Real(1.0) + qv0) /
            (amrex::Real(1.0) + (R_v/R_d)*qv0);
        const amrex::Real scale_height = R_d*kT0 /
            (CONST_GRAV*moist_hydrostatic_factor);
        const amrex::Real p0 = p_0*std::exp(-z/scale_height);
        return p0/(R_d*kT0*(amrex::Real(1.0) + (R_v/R_d)*qv0));
    }

    amrex::Real base_theta (int k, amrex::Real qv0) const
    {
        const amrex::Real z = (static_cast<amrex::Real>(k) + amrex::Real(0.5))*kDz;
        const amrex::Real moist_hydrostatic_factor = (amrex::Real(1.0) + qv0) /
            (amrex::Real(1.0) + (R_v/R_d)*qv0);
        const amrex::Real scale_height = R_d*kT0 /
            (CONST_GRAV*moist_hydrostatic_factor);
        const amrex::Real p0 = p_0*std::exp(-z/scale_height);
        return kT0*std::pow(p_0/p0, RdoCp);
    }

    amrex::Box domain;
    amrex::BoxArray grids;
    amrex::DistributionMapping mapping;
    amrex::RealBox real_box;
    amrex::Array<int, AMREX_SPACEDIM> periodic;
    amrex::Geometry geom;
    amrex::MultiFab base_state;
    amrex::MultiFab primitive;
    amrex::MultiFab total_water;
    amrex::MultiFab buoyancy;
    amrex::Vector<amrex::MultiFab> conserved;
    SolverChoice solver;
    std::unique_ptr<eb_> eb_factory;
};

void expect_near (amrex::Real actual, amrex::Real expected)
{
    EXPECT_LE(std::abs(actual-expected), kAbsTol + kRelTol*std::abs(expected));
}

void fill_anelastic_state (BuoyancyFixture& fixture, const StateParameters& state)
{
    fixture.fill_state(state, true);
    const auto density = fixture.density_and_base_density();
    for (int k = 0; k < kNz; ++k) {
        SCOPED_TRACE("fixed-density anelastic cell k=" + std::to_string(k));
        expect_near(density[static_cast<std::size_t>(2*k)],
                    density[static_cast<std::size_t>(2*k+1)]);
    }
}

amrex::Real type1_oracle (const BuoyancyFixture& fixture,
                           int face,
                           const StateParameters& state)
{
    const amrex::Real gz = -CONST_GRAV;
    const auto cell_perturbation = [&fixture, &state] (int k) {
        const amrex::Real rho0 = fixture.base_density(k, state.qv0);
        const amrex::Real theta0 = fixture.base_theta(k, state.qv0);
        const amrex::Real qv = state.qv0 + state.vapor_anomaly;
        const amrex::Real theta = state.theta_factor*theta0;
        // Derive rho/rho0 from the EOS closure: p scales as pressure_factor,
        // theta scales as theta_factor, and water changes the virtual-temperature factor.
        const amrex::Real rho = rho0 *
            std::pow(state.pressure_factor, iGamma) / state.theta_factor *
            (amrex::Real(1.0) + RvoRd*state.qv0) /
            (amrex::Real(1.0) + RvoRd*qv);
        const amrex::Real qt = qv + state.condensate;
        return rho*(amrex::Real(1.0)+qt) - rho0*(amrex::Real(1.0)+state.qv0);
    };
    return gz*amrex::Real(0.5)*(cell_perturbation(face-1)+cell_perturbation(face));
}


TEST(ERFBuoyancy, DryCompressibleType4NeutralUsesPrimitiveTheta)
{
    BuoyancyFixture fixture;
    const StateParameters neutral{};
    fixture.fill_state(neutral);
    const auto values = fixture.run(4, 0, MoistureType::None);

    EXPECT_EQ(values.front(), kTopSentinel);
    EXPECT_EQ(values.back(), kTopSentinel);
    for (int k = 1; k < kNz; ++k) {
        SCOPED_TRACE("interior face k=" + std::to_string(k));
        expect_near(values[static_cast<std::size_t>(k)], amrex::Real(0.0));
    }
}

TEST(ERFBuoyancy, DryCompressibleTemperatureSchemesKeepSelectorsAndPressureResponse)
{
    BuoyancyFixture fixture;
    const StateParameters pressure_only{amrex::Real(0.0), amrex::Real(1.0),
                                        amrex::Real(1.01), amrex::Real(0.0), amrex::Real(0.0)};
    fixture.fill_state(pressure_only);
    const auto type1 = fixture.run(1, 0, MoistureType::None);
    const auto type2 = fixture.run(2, 0, MoistureType::None);
    const auto type3 = fixture.run(3, 0, MoistureType::None);
    const auto type4 = fixture.run(4, 0, MoistureType::None);

    for (int k = 1; k < kNz; ++k) {
        SCOPED_TRACE("pressure-only interior face k=" + std::to_string(k));
        const auto idx = static_cast<std::size_t>(k);
        EXPECT_LT(type1[idx], amrex::Real(0.0));
        EXPECT_GT(type2[idx], amrex::Real(0.0));
        expect_near(type2[idx], type3[idx]);
        expect_near(type4[idx], amrex::Real(0.0));
        expect_near(type1[idx], type1_oracle(fixture, k, pressure_only));
    }

    for (const amrex::Real theta_factor : {amrex::Real(1.01), amrex::Real(0.99)}) {
        const StateParameters thermal{amrex::Real(0.0), theta_factor,
                                      amrex::Real(1.0), amrex::Real(0.0), amrex::Real(0.0)};
        fixture.fill_state(thermal);
        const auto theta_result = fixture.run(4, 0, MoistureType::None);
        for (int k = 1; k < kNz; ++k) {
            const amrex::Real rho_face = amrex::Real(0.5)*(
                fixture.base_density(k-1, amrex::Real(0.0)) +
                fixture.base_density(k, amrex::Real(0.0)));
            const amrex::Real expected = -rho_face*(-CONST_GRAV)*(theta_factor-amrex::Real(1.0));
            expect_near(theta_result[static_cast<std::size_t>(k)], expected);
            if (theta_factor > amrex::Real(1.0)) {
                EXPECT_GT(theta_result[static_cast<std::size_t>(k)], amrex::Real(0.0));
            } else {
                EXPECT_LT(theta_result[static_cast<std::size_t>(k)], amrex::Real(0.0));
            }
        }
    }
}

TEST(ERFBuoyancy, CompressibleMoistDensityReferenceAndCondensateLoading)
{
    BuoyancyFixture fixture;
    const StateParameters moist_reference{amrex::Real(0.01)};
    fixture.fill_state(moist_reference);
    const auto neutral = fixture.run(1, 0, MoistureType::Morrison);
    for (int k = 1; k < kNz; ++k) {
        expect_near(neutral[static_cast<std::size_t>(k)], amrex::Real(0.0));
    }

    const StateParameters loaded{amrex::Real(0.01), amrex::Real(1.0), amrex::Real(1.0),
                                 amrex::Real(0.0), amrex::Real(0.002)};
    fixture.fill_state(loaded);
    const auto buoyant_force = fixture.run(1, 0, MoistureType::Morrison);
    for (int k = 1; k < kNz; ++k) {
        const amrex::Real expected = type1_oracle(fixture, k, loaded);
        expect_near(buoyant_force[static_cast<std::size_t>(k)], expected);
        EXPECT_LT(buoyant_force[static_cast<std::size_t>(k)], amrex::Real(0.0));
    }
}

TEST(ERFBuoyancy, MoistTemperatureSelectorsRemainEquivalent)
{
    BuoyancyFixture fixture;
    const StateParameters moist_state{amrex::Real(0.01), amrex::Real(1.002),
                                      amrex::Real(1.0), amrex::Real(0.0005), amrex::Real(0.0002)};
    fixture.fill_state(moist_state);
    const auto type2 = fixture.run(2, 0, MoistureType::Morrison);
    const auto type3 = fixture.run(3, 0, MoistureType::Morrison);
    for (int k = 1; k < kNz; ++k) {
        expect_near(type2[static_cast<std::size_t>(k)], type3[static_cast<std::size_t>(k)]);
    }
}

TEST(ERFBuoyancy, MoistType4RetainsFiniteHumidityApproximation)
{
    BuoyancyFixture fixture;
    // Keep the anomaly small while separating the buoyancy signal from float roundoff.
    const amrex::Real theta_anomaly = amrex::Real(1.0e-2);
    StateParameters dry_state{amrex::Real(0.0), amrex::Real(1.0)+theta_anomaly};
    fixture.fill_state(dry_state);
    const auto dry_type1 = fixture.run(1, 0, MoistureType::Morrison);
    const auto dry_type4 = fixture.run(4, 0, MoistureType::Morrison);

    StateParameters moist_state{amrex::Real(0.02), amrex::Real(1.0)+theta_anomaly};
    fixture.fill_state(moist_state);
    const auto moist_type1 = fixture.run(1, 0, MoistureType::Morrison);
    const auto moist_type4 = fixture.run(4, 0, MoistureType::Morrison);

    for (int k = 1; k < kNz; ++k) {
        const auto idx = static_cast<std::size_t>(k);
        expect_near(dry_type4[idx]/dry_type1[idx], amrex::Real(1.0)+theta_anomaly);
        expect_near(moist_type4[idx]/moist_type1[idx],
                    (amrex::Real(1.0)+theta_anomaly)/(amrex::Real(1.0)+moist_state.qv0));
        EXPECT_LT(std::abs(moist_type4[idx]), std::abs(moist_type1[idx]));
    }
}

TEST(ERFBuoyancy, ActiveDryAnelasticUsesFaceRatioAndNeutralReference)
{
    BuoyancyFixture fixture;
    StateParameters state{};
    fill_anelastic_state(fixture, state);
    const auto neutral = fixture.run(3, 1, MoistureType::None);
    for (int k = 1; k < kNz; ++k) {
        expect_near(neutral[static_cast<std::size_t>(k)], amrex::Real(0.0));
    }

    state.theta_factor = amrex::Real(1.015);
    fill_anelastic_state(fixture, state);
    const auto warm = fixture.run(3, 1, MoistureType::None);
    for (int k = 1; k < kNz; ++k) {
        const amrex::Real theta_lo = state.theta_factor*fixture.base_theta(k-1, amrex::Real(0.0));
        const amrex::Real theta_hi = state.theta_factor*fixture.base_theta(k, amrex::Real(0.0));
        const amrex::Real theta0_lo = fixture.base_theta(k-1, amrex::Real(0.0));
        const amrex::Real theta0_hi = fixture.base_theta(k, amrex::Real(0.0));
        const amrex::Real rho0 = amrex::Real(0.5)*(
            fixture.base_density(k-1, amrex::Real(0.0)) +
            fixture.base_density(k, amrex::Real(0.0)));
        const amrex::Real expected = -rho0*(-CONST_GRAV)*
            (amrex::Real(0.5)*(theta_lo+theta_hi)-amrex::Real(0.5)*(theta0_lo+theta0_hi)) /
            (amrex::Real(0.5)*(theta0_lo+theta0_hi));
        expect_near(warm[static_cast<std::size_t>(k)], expected);
        EXPECT_GT(warm[static_cast<std::size_t>(k)], amrex::Real(0.0));
    }
}

TEST(ERFBuoyancy, ActiveMoistAnelasticIsNeutralAndRespondsToVaporAndLoading)
{
    BuoyancyFixture fixture;
    StateParameters state{amrex::Real(0.01)};
    fill_anelastic_state(fixture, state);
    const auto neutral = fixture.run(3, 1, MoistureType::Morrison);
    for (int k = 1; k < kNz; ++k) {
        expect_near(neutral[static_cast<std::size_t>(k)], amrex::Real(0.0));
    }

    state.vapor_anomaly = amrex::Real(0.001);
    fill_anelastic_state(fixture, state);
    const auto vapor = fixture.run(3, 1, MoistureType::Morrison);
    for (int k = 1; k < kNz; ++k) {
        const amrex::Real rho0 = amrex::Real(0.5)*(
            fixture.base_density(k-1, state.qv0) + fixture.base_density(k, state.qv0));
        const amrex::Real expected = -rho0*(-CONST_GRAV)*epsv*state.vapor_anomaly;
        expect_near(vapor[static_cast<std::size_t>(k)], expected);
        EXPECT_GT(vapor[static_cast<std::size_t>(k)], amrex::Real(0.0));
    }

    state.vapor_anomaly = amrex::Real(0.0);
    state.condensate = amrex::Real(0.001);
    fill_anelastic_state(fixture, state);
    const auto loaded = fixture.run(3, 1, MoistureType::Morrison);
    for (int k = 1; k < kNz; ++k) {
        const amrex::Real rho0 = amrex::Real(0.5)*(
            fixture.base_density(k-1, state.qv0) + fixture.base_density(k, state.qv0));
        const amrex::Real expected = -rho0*(-CONST_GRAV)*(-state.condensate);
        expect_near(loaded[static_cast<std::size_t>(k)], expected);
        EXPECT_LT(loaded[static_cast<std::size_t>(k)], amrex::Real(0.0));
    }
}

TEST(ERFBuoyancy, MoistFaceInertiaWeightsPressureAndBuoyancyTogether)
{
    BuoyancyFixture fixture;
    constexpr int face = 4;
    const amrex::Real qt_lo = amrex::Real(0.004);
    const amrex::Real qt_hi = amrex::Real(0.020);
    const amrex::Real gpz = amrex::Real(3.25);
    const amrex::Real abl_pressure_grad_z = amrex::Real(-1.75);
    const amrex::Real qt_face = amrex::Real(0.5)*(qt_lo+qt_hi);

    const auto check_force = [&] (amrex::Real buoyancy_force) {
        const amrex::Real numerator = -gpz - abl_pressure_grad_z + buoyancy_force;
        const amrex::Real expected_moist = numerator/(amrex::Real(1.0)+qt_face);
        const amrex::Real actual_moist = fixture.run_z_pressure_buoyancy_rhs(
            face, qt_lo, qt_hi, gpz, abl_pressure_grad_z, buoyancy_force, true);
        expect_near(actual_moist, expected_moist);

        const amrex::Real actual_dry = fixture.run_z_pressure_buoyancy_rhs(
            face, qt_lo, qt_hi, gpz, abl_pressure_grad_z, buoyancy_force, false);
        expect_near(actual_dry, numerator);
    };

    check_force(amrex::Real(2.5));
    check_force(amrex::Real(-1.25));

    const amrex::Real dry_limit = -gpz - abl_pressure_grad_z;
    const amrex::Real zero_water = fixture.run_z_pressure_buoyancy_rhs(
        face, amrex::Real(0.0), amrex::Real(0.0), gpz,
        abl_pressure_grad_z, amrex::Real(0.0), true);
    expect_near(zero_water, dry_limit);
}

TEST(ERFBuoyancy, ZeroGravityAndPhysicalVerticalFacesAreHandled)
{
    BuoyancyFixture fixture;
    fixture.fill_state(StateParameters{amrex::Real(0.0), amrex::Real(1.02)});
    const auto values = fixture.run(4, 0, MoistureType::None, amrex::Real(0.0));
    EXPECT_EQ(values.front(), kTopSentinel);
    EXPECT_EQ(values.back(), kTopSentinel);
    for (int k = 1; k < kNz; ++k) {
        EXPECT_EQ(values[static_cast<std::size_t>(k)], amrex::Real(0.0));
    }
}

} // namespace
