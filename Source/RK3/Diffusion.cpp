// Created by Pankaj Jha on 6/15/21.
#include <RK3.H>

using namespace amrex;

// Compute tau_ij (m + 1/2), tau_ij (m - 1/2) where m = {i, j, k} for DNS or Smagorinsky
Real ComputeStressTerm (const int &i, const int &j, const int &k,
                        const Array4<Real>& u, const Array4<Real>& v, const Array4<Real>& w,
                        const enum NextOrPrev &nextOrPrev,
                        const enum MomentumEqn &momentumEqn,
                        const enum DiffusionDir &diffDir,
                        const Geometry &geom,
                        const Array4<Real>& nut,
                        const SolverChoice &solverChoice,
                        const Real &molViscosity) {
    Real stressTerm = 0.0;
    Real turbViscInterpolated = 0.0;

    // Here, we have computed strain rate on the fly.
    // TODO: It may be better to store S11, S12 etc. at all the (m+1/2) and (m-1/2) grid points (edges) and use them here.
    Real strainRate = ComputeStrainRate(i, j, k, u, v, w, nextOrPrev, momentumEqn, diffDir, geom);

    // TODO: Consider passing turbModel to this function instead of computing it here from SolverChoice
    enum TurbulenceModel turbModel;

    //TODO: Update this to account for other turbulence models in future. This would alo require update in SolverChoice
    if (solverChoice.use_smagorinsky)
        turbModel = TurbulenceModel::Smagorinsky;
    else
        turbModel = TurbulenceModel::DNS;

    switch (turbModel) {
        case TurbulenceModel::DNS:
            stressTerm = 2.0*molViscosity * strainRate; // 2*nu*Sij(m+1/2) or 2*nu*Sij(m-1/2)
            break;
        case TurbulenceModel::Smagorinsky:
            turbViscInterpolated = InterpolateTurbulentViscosity(i, j, k, u, v, w, nextOrPrev, momentumEqn, diffDir, geom, nut);
            stressTerm = turbViscInterpolated * strainRate; // // K_interp*Sij(m+1/2) or K_interp*Sij(m-1/2)
            break;
        default:
            amrex::Abort("Error: Turbulence model is unrecognized");
    }

    return stressTerm;
}

Real ComputeDiffusionTermForState(const int &i, const int &j, const int &k,
                                  const Array4<Real>& qty, const int & qty_index,
                                  const enum Coord& coordDir) {
    Real diffusionTerm = 0.0;

    switch (coordDir) {
        case Coord::x:
            diffusionTerm = qty(i+1, j, k, qty_index) -2.0*qty(i, j, k, qty_index) + qty(i-1, j, k, qty_index);
            break;
        case Coord::y:
            diffusionTerm = qty(i, j+1, k, qty_index) -2.0*qty(i, j, k, qty_index) + qty(i, j-1, k, qty_index);
            break;
        case Coord::z:
            diffusionTerm = qty(i, j, k+1, qty_index) -2.0*qty(i, j, k, qty_index) + qty(i, j, k-1, qty_index);
            break;
        default:
            amrex::Abort("Error: Coord direction is unrecognized");
    }

    return diffusionTerm;
}
