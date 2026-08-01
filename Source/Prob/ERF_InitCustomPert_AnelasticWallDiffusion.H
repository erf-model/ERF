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

const auto prob_lo = geomdata.ProbLo();
const auto prob_hi = geomdata.ProbHi();
const auto dx = geomdata.CellSize();
const Real length = prob_hi[axis] - prob_lo[axis];

ParallelFor(bx, [=] AMREX_GPU_DEVICE(int i, int j, int k) noexcept
{
    const int index = (axis == 0) ? i : ((axis == 1) ? j : k);
    const Real coordinate = prob_lo[axis] + (index + Real(0.5)) * dx[axis];
    const Real theta = theta_lo + (theta_hi - theta_lo) *
                       (coordinate - prob_lo[axis]) / length;
    const Real rho = state(i,j,k,Rho_comp);
    state_pert(i,j,k,RhoTheta_comp) = rho * theta -
                                       state(i,j,k,RhoTheta_comp);
});
