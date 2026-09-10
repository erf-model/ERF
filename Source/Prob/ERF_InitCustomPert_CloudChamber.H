// Cloud Chamber Stage 1 initializer.
//
// ERF has already initialized the hydrostatic/base state when this snippet is
// called.  `state` is therefore the current conserved state and state_pert is
// an additive correction.  In anelastic mode density is fixed by that base
// state, so only RhoTheta and (for SatAdj) RhoQ1/RhoQ2 are corrected here.
// The analytic perturbation is evaluated at cell centers and uses no RNG;
// consequently it is independent of MPI decomposition, tiling, and restart
// layout.  theta is potential temperature, qv is vapor mixing ratio, and qc
// is SatAdj cloud-water mixing ratio; none of these prescribed values claim
// physical panel temperature, wet-wall evaporation, or a wall-transfer law.
AMREX_ALWAYS_ASSERT_WITH_MESSAGE(
    m_cloud_chamber_config.active,
    "Cloud Chamber: initializer reached without a parsed configuration");
AMREX_ALWAYS_ASSERT_WITH_MESSAGE(
    sc.anelastic[lev] == 1 && sc.fixed_density[lev] == 1,
    "Cloud Chamber: the selected level must use fixed-density anelastic dynamics");
AMREX_ALWAYS_ASSERT_WITH_MESSAGE(
    sc.terrain_type == TerrainType::None &&
    sc.buildings_type == BuildingsType::None &&
    sc.mesh_type == MeshType::ConstantDz,
    "Cloud Chamber: terrain, buildings, and non-Cartesian mesh modes are unsupported");
if (m_cloud_chamber_config.cloudy) {
    AMREX_ALWAYS_ASSERT_WITH_MESSAGE(
        sc.moisture_type == MoistureType::SatAdj,
        "Cloud Chamber: cloudy mode requires erf.moisture_model = SatAdj");
} else {
    AMREX_ALWAYS_ASSERT_WITH_MESSAGE(
        sc.moisture_type == MoistureType::None,
        "Cloud Chamber: dry mode requires no active moisture model");
}

const GpuArray<Real, AMREX_SPACEDIM> prob_lo = {
    geomdata.ProbLo(0), geomdata.ProbLo(1), geomdata.ProbLo(2)};
const GpuArray<Real, AMREX_SPACEDIM> dx = {
    geomdata.CellSize(0), geomdata.CellSize(1), geomdata.CellSize(2)};
const auto config = m_cloud_chamber_config;
const Real rdOcp = sc.rdOcp;

ParallelFor(bx, [=] AMREX_GPU_DEVICE (int i, int j, int k) noexcept
{
    const Real x = prob_lo[0] + (Real(i) + Real(0.5)) * dx[0];
    const Real y = prob_lo[1] + (Real(j) + Real(0.5)) * dx[1];
    const Real z = prob_lo[2] + (Real(k) + Real(0.5)) * dx[2];
    const Real rho = state(i,j,k,Rho_comp);
    const Real temperature = erf_cloud_chamber::linear_profile(
        config.initial_temperature_bottom, config.initial_temperature_top,
        z, config.prob_lo[2], config.prob_hi[2] - config.prob_lo[2]) +
        erf_cloud_chamber::deterministic_perturbation(
            x, y, z, config.prob_lo,
            GpuArray<Real, AMREX_SPACEDIM>{
                config.prob_hi[0] - config.prob_lo[0],
                config.prob_hi[1] - config.prob_lo[1],
                config.prob_hi[2] - config.prob_lo[2]},
            config.temperature_perturbation_amplitude);
    const Real target_theta = config.physical_initialization ?
        temperature * std::pow(p_0 / p_hse(i,j,k), rdOcp) :
        erf_cloud_chamber::theta_at(config, x, y, z);
    state_pert(i,j,k,RhoTheta_comp) =
        rho * target_theta - state(i,j,k,RhoTheta_comp);

    if (config.cloudy) {
        const Real target_qv = config.physical_initialization ?
            erf_cloud_chamber::vapor_mixing_ratio_from_relative_humidity(
                temperature, p_hse(i,j,k), config.initial_relative_humidity) :
            erf_cloud_chamber::qv_at(config, z);
        state_pert(i,j,k,RhoQ1_comp) =
            rho * target_qv - state(i,j,k,RhoQ1_comp);
        // SatAdj owns the first adjustment of cloud water.  Starting from
        // zero makes that source explicit instead of carrying stale state.
        state_pert(i,j,k,RhoQ2_comp) = -state(i,j,k,RhoQ2_comp);
    }
});
