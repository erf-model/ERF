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

const auto prob_lo = geomdata.ProbLo();
const auto dx = geomdata.CellSize();
const auto config = m_cloud_chamber_config;

ParallelFor(bx, [=] AMREX_GPU_DEVICE (int i, int j, int k) noexcept
{
    const Real x = prob_lo[0] + (Real(i) + Real(0.5)) * dx[0];
    const Real y = prob_lo[1] + (Real(j) + Real(0.5)) * dx[1];
    const Real z = prob_lo[2] + (Real(k) + Real(0.5)) * dx[2];
    const Real rho = state(i,j,k,Rho_comp);
    const Real target_theta = erf_cloud_chamber::theta_at(config, x, y, z);
    state_pert(i,j,k,RhoTheta_comp) =
        rho * target_theta - state(i,j,k,RhoTheta_comp);

    if (config.cloudy) {
        const Real target_qv = erf_cloud_chamber::qv_at(config, z);
        state_pert(i,j,k,RhoQ1_comp) =
            rho * target_qv - state(i,j,k,RhoQ1_comp);
        // SatAdj owns the first adjustment of cloud water.  Starting from
        // zero makes that source explicit instead of carrying stale state.
        state_pert(i,j,k,RhoQ2_comp) = -state(i,j,k,RhoQ2_comp);
    }
});
