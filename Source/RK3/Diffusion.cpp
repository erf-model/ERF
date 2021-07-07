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
                        const SolverChoice &solverChoice) {
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
            stressTerm = 2.0 * solverChoice.kinematicViscosity * strainRate; // 2*nu*Sij(m+1/2) or 2*nu*Sij(m-1/2)
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

Real
DiffusionContributionForMom(const int &i, const int &j, const int &k,
                            const Array4<Real>& u, const Array4<Real>& v, const Array4<Real>& w,
                            const enum MomentumEqn &momentumEqn,
                            const amrex::Geometry &geom,
                            const Array4<Real>& nut,
                            const SolverChoice &solverChoice) {

    const GpuArray<Real, AMREX_SPACEDIM> cellSize = geom.CellSizeArray();
    auto dx = cellSize[0], dy = cellSize[1], dz = cellSize[2];
    Real diffusionContribution = 0.0;

    switch (momentumEqn) {
        case MomentumEqn::x:
            Real tau11Next, tau11Prev, tau12Next, tau12Prev, tau13Next, tau13Prev;
            tau11Next = ComputeStressTerm(i, j, k, u, v, w, NextOrPrev::next, momentumEqn,DiffusionDir::x, geom, nut, solverChoice);
            tau11Prev = ComputeStressTerm(i, j, k, u, v, w, NextOrPrev::prev, momentumEqn,DiffusionDir::x, geom, nut, solverChoice);
            tau12Next = ComputeStressTerm(i, j, k, u, v, w, NextOrPrev::next, momentumEqn,DiffusionDir::y, geom, nut, solverChoice);
            tau12Prev = ComputeStressTerm(i, j, k, u, v, w, NextOrPrev::prev, momentumEqn,DiffusionDir::y, geom, nut, solverChoice);
            tau13Next = ComputeStressTerm(i, j, k, u, v, w, NextOrPrev::next, momentumEqn,DiffusionDir::z, geom, nut, solverChoice);
            tau13Prev = ComputeStressTerm(i, j, k, u, v, w, NextOrPrev::prev, momentumEqn,DiffusionDir::z, geom, nut, solverChoice);

            diffusionContribution = (tau11Next - tau11Prev) / dx  // Contribution to x-mom eqn from diffusive flux in x-dir
                                  + (tau12Next - tau12Prev) / dy  // Contribution to x-mom eqn from diffusive flux in y-dir
                                  + (tau13Next - tau13Prev) / dz; // Contribution to x-mom eqn from diffusive flux in z-dir
            break;
        case MomentumEqn::y:
            Real tau21Next, tau21Prev, tau22Next, tau22Prev, tau23Next, tau23Prev;
            tau21Next = ComputeStressTerm(i, j, k, u, v, w, NextOrPrev::next, momentumEqn, DiffusionDir::x, geom, nut, solverChoice);
            tau21Prev = ComputeStressTerm(i, j, k, u, v, w, NextOrPrev::prev, momentumEqn, DiffusionDir::x, geom, nut, solverChoice);
            tau22Next = ComputeStressTerm(i, j, k, u, v, w, NextOrPrev::next, momentumEqn, DiffusionDir::y, geom, nut, solverChoice);
            tau22Prev = ComputeStressTerm(i, j, k, u, v, w, NextOrPrev::prev, momentumEqn, DiffusionDir::y, geom, nut, solverChoice);
            tau23Next = ComputeStressTerm(i, j, k, u, v, w, NextOrPrev::next, momentumEqn, DiffusionDir::z, geom, nut, solverChoice);
            tau23Prev = ComputeStressTerm(i, j, k, u, v, w, NextOrPrev::prev, momentumEqn, DiffusionDir::z, geom, nut, solverChoice);

            diffusionContribution = (tau21Next - tau21Prev) / dx  // Contribution to y-mom eqn from diffusive flux in x-dir
                                  + (tau22Next - tau22Prev) / dy  // Contribution to y-mom eqn from diffusive flux in y-dir
                                  + (tau23Next - tau23Prev) / dz; // Contribution to y-mom eqn from diffusive flux in z-dir
            break;
        case MomentumEqn::z:
            Real tau31Next, tau31Prev, tau32Next, tau32Prev, tau33Next, tau33Prev;
            tau31Next = ComputeStressTerm(i, j, k, u, v, w, NextOrPrev::next, momentumEqn, DiffusionDir::x, geom, nut, solverChoice);
            tau31Prev = ComputeStressTerm(i, j, k, u, v, w, NextOrPrev::prev, momentumEqn, DiffusionDir::x, geom, nut, solverChoice);
            tau32Next = ComputeStressTerm(i, j, k, u, v, w, NextOrPrev::next, momentumEqn, DiffusionDir::y, geom, nut, solverChoice);
            tau32Prev = ComputeStressTerm(i, j, k, u, v, w, NextOrPrev::prev, momentumEqn, DiffusionDir::y, geom, nut, solverChoice);
            tau33Next = ComputeStressTerm(i, j, k, u, v, w, NextOrPrev::next, momentumEqn, DiffusionDir::z, geom, nut, solverChoice);
            tau33Prev = ComputeStressTerm(i, j, k, u, v, w, NextOrPrev::prev, momentumEqn, DiffusionDir::z, geom, nut, solverChoice);

            diffusionContribution = (tau31Next - tau31Prev) / dx  // Contribution to z-mom eqn from diffusive flux in x-dir
                                  + (tau32Next - tau32Prev) / dy  // Contribution to z-mom eqn from diffusive flux in y-dir
                                  + (tau33Next - tau33Prev) / dz; // Contribution to z-mom eqn from diffusive flux in z-dir
            break;
        default:
            amrex::Abort("Error: Momentum equation is unrecognized");
    }

    return diffusionContribution;
}

Real ComputeDiffusionTermForState(const int &i, const int &j, const int &k,
                                  const Array4<Real>& cell_data, const int & qty_index,
                                  const enum Coord& coordDir) {
    Real diffusionTerm = 0.0;

    switch (coordDir) {
        case Coord::x:
            diffusionTerm = cell_data(i+1, j, k, qty_index) -2.0*cell_data(i, j, k, qty_index) + cell_data(i-1, j, k, qty_index);
            break;
        case Coord::y:
            diffusionTerm = cell_data(i, j+1, k, qty_index) -2.0*cell_data(i, j, k, qty_index) + cell_data(i, j-1, k, qty_index);
            break;
        case Coord::z:
            diffusionTerm = cell_data(i, j, k+1, qty_index) -2.0*cell_data(i, j, k, qty_index) + cell_data(i, j, k-1, qty_index);
            break;
        default:
            amrex::Abort("Error: Coord direction is unrecognized");
    }

    return diffusionTerm;
}

Real
DiffusionContributionForState(const int &i, const int &j, const int &k,
                              const Array4<Real>& cell_data, const int & qty_index,
                              const amrex::Geometry &geom,
                              const SolverChoice &solverChoice) {

    const GpuArray<Real, AMREX_SPACEDIM> cellSize = geom.CellSizeArray();
    auto dx = cellSize[0], dy = cellSize[1], dz = cellSize[2];
    Real xDiffFlux = 0.0, yDiffFlux = 0.0, zDiffFlux = 0.0, diffCoeff = 0.0;

    xDiffFlux = ComputeDiffusionTermForState(i, j, k, cell_data, qty_index, Coord::x);
    yDiffFlux = ComputeDiffusionTermForState(i, j, k, cell_data, qty_index, Coord::y);
    zDiffFlux = ComputeDiffusionTermForState(i, j, k, cell_data, qty_index, Coord::z);

    switch(qty_index) {
        case RhoTheta_comp: // Temperature
            diffCoeff = solverChoice.alpha_T;
            break;
        case RhoScalar_comp: // Scalar
            diffCoeff = solverChoice.alpha_S;
            break;
        default:
            amrex::Abort("Error: Diffusion term for the data index isn't implemented");
    }

    // Assemble diffusion contribution.
    Real diffusionContribution = diffCoeff *
                                 (  xDiffFlux / (dx*dx)  // Diffusive flux in x-dir
                                  + yDiffFlux / (dy*dy)  // Diffusive flux in y-dir
                                  + zDiffFlux / (dz*dz)  // Diffusive flux in z-dir
                                 );

    return diffusionContribution;
}
