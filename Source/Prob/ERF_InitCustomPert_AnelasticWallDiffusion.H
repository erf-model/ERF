// Stage 0 manufactured state for the six-wall anelastic chamber fixtures.
// The perturbation is written as rho*theta so the existing initialization
// wrapper can add it to the conserved state without changing density.
ParmParse pp_prob("prob");
int axis = 0;
Real theta_lo = 300.0;
Real theta_hi = 301.0;
pp_prob.query("axis", axis);
pp_prob.query("theta_lo", theta_lo);
pp_prob.query("theta_hi", theta_hi);

const Real problo_axis = geomdata.ProbLo(axis);
const Real dx_axis     = geomdata.CellSize(axis);
const Real length      = geomdata.ProbHi(axis) - geomdata.ProbLo(axis);

ParallelFor(bx, [=] AMREX_GPU_DEVICE(int i, int j, int k) noexcept
{
    const int index = (axis == 0) ? i : ((axis == 1) ? j : k);
    const Real coordinate = problo_axis + (static_cast<Real>(index) + Real(0.5)) * dx_axis;
    const Real theta = theta_lo + (theta_hi - theta_lo) * (coordinate - problo_axis) / length;
    const Real rho = state(i,j,k,Rho_comp);
    state_pert(i,j,k,RhoTheta_comp) = rho * theta - state(i,j,k,RhoTheta_comp);
});
