#include <TimeIntegration.H>

using namespace amrex;

AMREX_GPU_DEVICE
Real
ComputeAdvectedQuantityForMom(const int &i, const int &j, const int &k,
                              const Array4<Real>& rho_u, const Array4<Real>& rho_v, const Array4<Real>& rho_w,
                              const Array4<Real>& u, const Array4<Real>& v, const Array4<Real>& w,
                              const enum NextOrPrev &nextOrPrev,
                              const enum AdvectedQuantity &advectedQuantity,
                              const enum AdvectingQuantity &advectingQuantity,
                              const int &spatial_order) {
  Real advectingQty = 0.0;
  Real advectedQty = 1.0;
  if (nextOrPrev == NextOrPrev::next) {
    switch(advectedQuantity) {
    case AdvectedQuantity::u: //x-momentum, reference face index is (i, j, k)
      switch (advectingQuantity) {
      case AdvectingQuantity::rho_u:
        advectedQty = InterpolateFromCellOrFace(i, j, k, u, 0, nextOrPrev, Coord::x, spatial_order); // u(i+1/2,    j, k    )
        advectingQty = 0.5*(rho_u(i+1, j, k) + rho_u(i, j, k)); // Effectively rho_u (i+1/2, j, k)
        break;
      case AdvectingQuantity::rho_v:
        advectedQty = InterpolateFromCellOrFace(i, j, k, u, 0, nextOrPrev, Coord::y, spatial_order); // u(i   , j+1/2, k    )
        advectingQty = 0.5*(rho_v(i, j+1, k) + rho_v(i-1, j+1, k)); // Effectively rho_v (i-1/2, j+1, k)
        break;
      case AdvectingQuantity::rho_w:
        advectedQty = InterpolateFromCellOrFace(i, j, k, u, 0, nextOrPrev, Coord::z, spatial_order); // u(i   , j    , k+1/2)
        advectingQty = 0.5*(rho_w(i, j, k+1) + rho_w(i-1, j, k+1)); // Effectively rho_w (i-1/2, j, k+1)
        break;
      default:
        amrex::Abort("Error: Advecting quantity is unrecognized");
      }
      break;
    case AdvectedQuantity::v: //y-momentum, reference face index is (i, j, k)
      switch (advectingQuantity) {
      case AdvectingQuantity::rho_u:
        advectedQty = InterpolateFromCellOrFace(i, j, k, v, 0, nextOrPrev, Coord::x, spatial_order); // v(i+1/2,    j, k    )
        advectingQty = 0.5*(rho_u(i+1, j, k) + rho_u(i+1, j-1, k)); // Effectively rho_u (i+1, j-1/2, k)
        break;
      case AdvectingQuantity::rho_v:
        advectedQty = InterpolateFromCellOrFace(i, j, k, v, 0, nextOrPrev, Coord::y, spatial_order); // v(i   , j+1/2, k    )
        advectingQty = 0.5*(rho_v(i, j+1, k) + rho_v(i, j, k)); // Effectively rho_v (i, j+1/2, k)
        break;
      case AdvectingQuantity::rho_w:
        advectedQty = InterpolateFromCellOrFace(i, j, k, v, 0, nextOrPrev, Coord::z, spatial_order); // v(i   , j    , k+1/2)
        advectingQty = 0.5*(rho_w(i, j, k+1) + rho_w(i, j-1, k+1)); // Effectively rho_w (i, j-1/2, k+1)
        break;
      default:
        amrex::Abort("Error: Advecting quantity is unrecognized");
      }
      break;
    case AdvectedQuantity::w: //z-momentum, reference face index is (i, j, k)
      switch (advectingQuantity) {
      case AdvectingQuantity::rho_u:
        advectedQty = InterpolateFromCellOrFace(i, j, k, w, 0, nextOrPrev, Coord::x, spatial_order); // w(i+1/2,    j, k    )
        advectingQty = 0.5*(rho_u(i+1, j, k) + rho_u(i+1, j, k-1)); // Effectively rho_u (i+1, j, k-1/2)
        break;
      case AdvectingQuantity::rho_v:
        advectedQty = InterpolateFromCellOrFace(i, j, k, w, 0, nextOrPrev, Coord::y, spatial_order); // w(i   , j+1/2, k    )
        advectingQty = 0.5*(rho_v(i, j+1, k) + rho_v(i, j+1, k-1)); // Effectively rho_v (i, j+1, k-1/2)
        break;
      case AdvectingQuantity::rho_w:
        advectedQty = InterpolateFromCellOrFace(i, j, k, w, 0, nextOrPrev, Coord::z, spatial_order); // w(i   , j    , k+1/2)
        advectingQty = 0.5*(rho_w(i, j, k+1) + rho_w(i, j, k)); // Effectively rho_w (i, j, k+1/2)
        break;
      default:
        amrex::Abort("Error: Advecting quantity is unrecognized");
      }
      break;
    default:
      amrex::Abort("Error: Advected quantity is unrecognized");
    }
  }
  else { // nextOrPrev == NextOrPrev::prev
    switch(advectedQuantity) {
    case AdvectedQuantity::u: //x-momentum, reference face index is (i, j, k)
      switch (advectingQuantity) {
      case AdvectingQuantity::rho_u:
        advectedQty = InterpolateFromCellOrFace(i, j, k, u, 0, nextOrPrev, Coord::x, spatial_order); // u(i-1/2,    j, k    )
        advectingQty = 0.5*(rho_u(i-1, j, k) + rho_u(i, j, k)); // Effectively rho_u (i-1/2, j, k)
        break;
      case AdvectingQuantity::rho_v:
        advectedQty = InterpolateFromCellOrFace(i, j, k, u, 0, nextOrPrev, Coord::y, spatial_order); // u(i   , j-1/2, k    )
        advectingQty = 0.5*(rho_v(i, j, k) + rho_v(i-1, j, k)); // Effectively rho_v (i-1/2, j, k)
        break;
      case AdvectingQuantity::rho_w:
        advectedQty = InterpolateFromCellOrFace(i, j, k, u, 0, nextOrPrev, Coord::z, spatial_order); // u(i   , j    , k-1/2)
        advectingQty = 0.5*(rho_w(i, j, k) + rho_w(i-1, j, k)); // Effectively rho_w (i-1/2, j, k)
        break;
      default:
        amrex::Abort("Error: Advecting quantity is unrecognized");
      }
      break;
    case AdvectedQuantity::v: //y-momentum, reference face index is (i, j, k)
      switch (advectingQuantity) {
      case AdvectingQuantity::rho_u:
        advectedQty = InterpolateFromCellOrFace(i, j, k, v, 0, nextOrPrev, Coord::x, spatial_order); // v(i-1/2,    j, k    )
        advectingQty = 0.5*(rho_u(i, j, k) + rho_u(i, j-1, k)); // Effectively rho_u (i, j-1/2, k)
        break;
      case AdvectingQuantity::rho_v:
        advectedQty = InterpolateFromCellOrFace(i, j, k, v, 0, nextOrPrev, Coord::y, spatial_order); // v(i   , j-1/2, k    )
        advectingQty = 0.5*(rho_v(i, j, k) + rho_v(i, j-1, k)); // Effectively rho_v (i, j-1/2, k)
        break;
      case AdvectingQuantity::rho_w:
        advectedQty = InterpolateFromCellOrFace(i, j, k, v, 0, nextOrPrev, Coord::z, spatial_order); // v(i   , j    , k-1/2)
        advectingQty = 0.5*(rho_w(i, j, k) + rho_w(i, j-1, k)); // Effectively rho_w (i, j-1/2, k)
        break;
      default:
        amrex::Abort("Error: Advecting quantity is unrecognized");
      }
      break;
    case AdvectedQuantity::w: //z-momentum, reference face index is (i, j, k)
      switch (advectingQuantity) {
      case AdvectingQuantity::rho_u:
        advectedQty = InterpolateFromCellOrFace(i, j, k, w, 0, nextOrPrev, Coord::x, spatial_order); // w(i-1/2,    j, k    )
        advectingQty = 0.5*(rho_u(i, j, k) + rho_u(i, j, k-1)); // Effectively rho_u (i, j, k-1/2)
        break;
      case AdvectingQuantity::rho_v:
        advectedQty = InterpolateFromCellOrFace(i, j, k, w, 0, nextOrPrev, Coord::y, spatial_order); // w(i   , j-1/2, k    )
        advectingQty = 0.5*(rho_v(i, j, k) + rho_v(i, j, k-1)); // Effectively rho_v (i, j, k-1/2)
        break;
      case AdvectingQuantity::rho_w:
        advectedQty = InterpolateFromCellOrFace(i, j, k, w, 0, nextOrPrev, Coord::z, spatial_order); // w(i   , j    , k-1/2)
        advectingQty = 0.5*(rho_w(i, j, k) + rho_w(i, j, k-1)); // Effectively rho_w (i, j, k-1/2)
        break;
      default:
        amrex::Abort("Error: Advecting quantity is unrecognized");
      }
      break;
    default:
      amrex::Abort("Error: Advected quantity is unrecognized");
    }
  }

  return advectingQty * advectedQty;
}

AMREX_GPU_DEVICE
Real
AdvectionContributionForMom(const int &i, const int &j, const int &k,
                            const Array4<Real>& rho_u, const Array4<Real>& rho_v, const Array4<Real>& rho_w,
                            const Array4<Real>& u, const Array4<Real>& v, const Array4<Real>& w,
                            const enum MomentumEqn &momentumEqn,
                            const GpuArray<Real, AMREX_SPACEDIM>& cellSize,
                            const SolverChoice &solverChoice) {

    auto dx = cellSize[0], dy = cellSize[1], dz = cellSize[2];
    Real advectionContribution = 0.0;

    int l_spatial_order = solverChoice.spatial_order;

    switch (momentumEqn) {
        case MomentumEqn::x: //x-momentum, reference face index is (i, j, k)
            Real centFluxXXNext, centFluxXXPrev, edgeFluxXYNext, edgeFluxXYPrev, edgeFluxXZNext, edgeFluxXZPrev;
            centFluxXXNext = ComputeAdvectedQuantityForMom(i, j, k, rho_u, rho_v, rho_w, u, v, w,                                                                     NextOrPrev::next, AdvectedQuantity::u, AdvectingQuantity::rho_u, l_spatial_order);
            centFluxXXPrev = ComputeAdvectedQuantityForMom(i, j, k, rho_u, rho_v, rho_w, u, v, w,                                                                     NextOrPrev::prev, AdvectedQuantity::u, AdvectingQuantity::rho_u, l_spatial_order);
            edgeFluxXYNext = ComputeAdvectedQuantityForMom(i, j, k, rho_u, rho_v, rho_w, u, v, w,                                                                     NextOrPrev::next, AdvectedQuantity::u, AdvectingQuantity::rho_v, l_spatial_order);
            edgeFluxXYPrev = ComputeAdvectedQuantityForMom(i, j, k, rho_u, rho_v, rho_w, u, v, w,                                                                     NextOrPrev::prev, AdvectedQuantity::u, AdvectingQuantity::rho_v, l_spatial_order);
            edgeFluxXZNext = ComputeAdvectedQuantityForMom(i, j, k, rho_u, rho_v, rho_w, u, v, w,                                                                     NextOrPrev::next, AdvectedQuantity::u, AdvectingQuantity::rho_w, l_spatial_order);
            edgeFluxXZPrev = ComputeAdvectedQuantityForMom(i, j, k, rho_u, rho_v, rho_w, u, v, w,                                                                     NextOrPrev::prev, AdvectedQuantity::u, AdvectingQuantity::rho_w, l_spatial_order);

            advectionContribution = (centFluxXXNext - centFluxXXPrev) / dx   // Contribution to x-mom eqn from advective flux in x-dir
                                  + (edgeFluxXYNext - edgeFluxXYPrev) / dy   // Contribution to x-mom eqn from advective flux in y-dir
                                  + (edgeFluxXZNext - edgeFluxXZPrev) / dz;  // Contribution to x-mom eqn from advective flux in z-dir
            break;
        case MomentumEqn::y: //y-momentum, reference face index is (i, j, k)
            Real centFluxYYNext, centFluxYYPrev, edgeFluxYXNext, edgeFluxYXPrev, edgeFluxYZNext, edgeFluxYZPrev;
            edgeFluxYXNext = ComputeAdvectedQuantityForMom(i, j, k, rho_u, rho_v, rho_w, u, v, w,                                                                     NextOrPrev::next, AdvectedQuantity::v, AdvectingQuantity::rho_u, l_spatial_order);
            edgeFluxYXPrev = ComputeAdvectedQuantityForMom(i, j, k, rho_u, rho_v, rho_w, u, v, w,                                                                     NextOrPrev::prev, AdvectedQuantity::v, AdvectingQuantity::rho_u, l_spatial_order);
            centFluxYYNext = ComputeAdvectedQuantityForMom(i, j, k, rho_u, rho_v, rho_w, u, v, w,                                                                     NextOrPrev::next, AdvectedQuantity::v, AdvectingQuantity::rho_v, l_spatial_order);
            centFluxYYPrev = ComputeAdvectedQuantityForMom(i, j, k, rho_u, rho_v, rho_w, u, v, w,                                                                     NextOrPrev::prev, AdvectedQuantity::v, AdvectingQuantity::rho_v, l_spatial_order);
            edgeFluxYZNext = ComputeAdvectedQuantityForMom(i, j, k, rho_u, rho_v, rho_w, u, v, w,                                                                     NextOrPrev::next, AdvectedQuantity::v, AdvectingQuantity::rho_w, l_spatial_order);
            edgeFluxYZPrev = ComputeAdvectedQuantityForMom(i, j, k, rho_u, rho_v, rho_w, u, v, w,                                                                     NextOrPrev::prev, AdvectedQuantity::v, AdvectingQuantity::rho_w, l_spatial_order);

            advectionContribution = (edgeFluxYXNext - edgeFluxYXPrev) / dx   // Contribution to y-mom eqn from advective flux in x-dir
                                  + (centFluxYYNext - centFluxYYPrev) / dy   // Contribution to y-mom eqn from advective flux in y-dir
                                  + (edgeFluxYZNext - edgeFluxYZPrev) / dz;  // Contribution to y-mom eqn from advective flux in z-dir
            break;
        case MomentumEqn::z: //z-momentum, reference face index is (i, j, k)
            Real centFluxZZNext, centFluxZZPrev, edgeFluxZXNext, edgeFluxZXPrev, edgeFluxZYNext, edgeFluxZYPrev;
            edgeFluxZXNext = ComputeAdvectedQuantityForMom(i, j, k, rho_u, rho_v, rho_w, u, v, w,                                                                     NextOrPrev::next, AdvectedQuantity::w, AdvectingQuantity::rho_u, l_spatial_order);
            edgeFluxZXPrev = ComputeAdvectedQuantityForMom(i, j, k, rho_u, rho_v, rho_w, u, v, w,                                                                     NextOrPrev::prev, AdvectedQuantity::w, AdvectingQuantity::rho_u, l_spatial_order);
            edgeFluxZYNext = ComputeAdvectedQuantityForMom(i, j, k, rho_u, rho_v, rho_w, u, v, w,                                                                     NextOrPrev::next, AdvectedQuantity::w, AdvectingQuantity::rho_v, l_spatial_order);
            edgeFluxZYPrev = ComputeAdvectedQuantityForMom(i, j, k, rho_u, rho_v, rho_w, u, v, w,                                                                     NextOrPrev::prev, AdvectedQuantity::w, AdvectingQuantity::rho_v, l_spatial_order);
            centFluxZZNext = ComputeAdvectedQuantityForMom(i, j, k, rho_u, rho_v, rho_w, u, v, w,                                                                     NextOrPrev::next, AdvectedQuantity::w, AdvectingQuantity::rho_w, l_spatial_order);
            centFluxZZPrev = ComputeAdvectedQuantityForMom(i, j, k, rho_u, rho_v, rho_w, u, v, w,                                                                     NextOrPrev::prev, AdvectedQuantity::w, AdvectingQuantity::rho_w, l_spatial_order);

            advectionContribution = (edgeFluxZXNext - edgeFluxZXPrev) / dx  // Contribution to z-mom eqn from advective flux in x-dir
                                  + (edgeFluxZYNext - edgeFluxZYPrev) / dy  // Contribution to z-mom eqn from advective flux in y-dir
                                  + (centFluxZZNext - centFluxZZPrev) / dz; // Contribution to z-mom eqn from advective flux in z-dir
            break;
        default:
            amrex::Abort("Error: Momentum equation is unrecognized");
    }

    return advectionContribution;
}

AMREX_GPU_DEVICE
Real
ComputeAdvectedQuantityForState(const int &i, const int &j, const int &k,
                                const Array4<Real>& cell_data, const int& qty_index,
                                const enum AdvectingQuantity &advectingQuantity,
                                const int &spatial_order) {
  Real advectingQty = 0.0;
  Real advectedQty = 1.0;

  AdvectedQuantity advectedQuantity;

  switch(qty_index) {
        case RhoTheta_comp: // Temperature
            advectedQuantity = AdvectedQuantity::theta;
            break;
        case RhoScalar_comp: // Scalar
            advectedQuantity = AdvectedQuantity::scalar;
            break;
        default:
            amrex::Abort("Error: Conserved quantity index is unrecognized");
    }

  // HACK HACK HACK 
  NextOrPrev nextOrPrev = NextOrPrev::prev;

  // Compute advected quantity for different choice of AdvectingQuantity
  switch(advectedQuantity) {

  case AdvectedQuantity::theta:
    switch (advectingQuantity) { // reference cell is (i, j, k)
    case AdvectingQuantity::rho_u:
      // Get theta (i+1/2,    j, k    ) = theta on face (i+1, j  , k  ) for x-dir if nextOrPrev = NextOrPrev::next
      // Get theta (i-1/2,    j, k    ) = theta on face (i,   j  , k  ) for x-dir if nextOrPrev = NextOrPrev::prev
      advectedQty = InterpolateRhoThetaFromCellToFace(i, j, k, cell_data, nextOrPrev, Coord::x, spatial_order);
      advectedQty/= InterpolateDensityFromCellToFace(i, j, k, cell_data, nextOrPrev, Coord::x, spatial_order);
      break;
    case AdvectingQuantity::rho_v:
      // Get theta (i   , j+1/2, k    ) = theta on face (i  , j+1, k  ) for y-dir if nextOrPrev = NextOrPrev::next
      // Get theta (i   , j-1/2, k    ) = theta on face (i  , j  , k  ) for y-dir if nextOrPrev = NextOrPrev::prev
      advectedQty = InterpolateRhoThetaFromCellToFace(i, j, k, cell_data, nextOrPrev, Coord::y, spatial_order);
      advectedQty/= InterpolateDensityFromCellToFace(i, j, k, cell_data, nextOrPrev, Coord::y, spatial_order);
      break;
    case AdvectingQuantity::rho_w:
      // Get theta (i   , j    , k+1/2) = theta on face (i  , j  , k+1) for z-dir if nextOrPrev = NextOrPrev::next
      // Get theta (i   , j    , k-1/2) = theta on face (i  , j  , k  ) for z-dir if nextOrPrev = NextOrPrev::prev
      advectedQty = InterpolateRhoThetaFromCellToFace(i, j, k, cell_data, nextOrPrev, Coord::z, spatial_order);
      advectedQty/= InterpolateDensityFromCellToFace(i, j, k, cell_data, nextOrPrev, Coord::z, spatial_order);
      break;
    default:
      amrex::Abort("Error: Advecting quantity is unrecognized");
    }
    break;

  case AdvectedQuantity::scalar:
    switch (advectingQuantity) { // reference cell is (i, j, k)
    case AdvectingQuantity::rho_u:
      // Get scalar (i+1/2,    j, k    ) = scalar on face (i+1, j  , k  ) for x-dir if nextOrPrev = NextOrPrev::next
      // Get scalar (i-1/2,    j, k    ) = scalar on face (i,   j  , k  ) for x-dir if nextOrPrev = NextOrPrev::prev
      advectedQty = InterpolateRhoScalarFromCellToFace(i, j, k, cell_data, nextOrPrev, Coord::x, spatial_order);
      advectedQty/= InterpolateDensityFromCellToFace(i, j, k, cell_data, nextOrPrev, Coord::x, spatial_order);
      break;
    case AdvectingQuantity::rho_v:
      // Get scalar (i   , j+1/2, k    ) = scalar on face (i  , j+1, k  ) for y-dir if nextOrPrev = NextOrPrev::next
      // Get scalar (i   , j-1/2, k    ) = scalar on face (i  , j  , k  ) for y-dir if nextOrPrev = NextOrPrev::prev
      advectedQty = InterpolateRhoScalarFromCellToFace(i, j, k, cell_data, nextOrPrev, Coord::y, spatial_order);
      advectedQty/= InterpolateDensityFromCellToFace(i, j, k, cell_data, nextOrPrev, Coord::y, spatial_order);
      break;
    case AdvectingQuantity::rho_w:
      // Get scalar (i   , j    , k+1/2) = scalar on face (i  , j  , k+1) for z-dir if nextOrPrev = NextOrPrev::next
      // Get scalar (i   , j    , k-1/2) = scalar on face (i  , j  , k  ) for z-dir if nextOrPrev = NextOrPrev::prev
      advectedQty = InterpolateRhoScalarFromCellToFace(i, j, k, cell_data, nextOrPrev, Coord::z, spatial_order);
      advectedQty/= InterpolateDensityFromCellToFace(i, j, k, cell_data, nextOrPrev, Coord::z, spatial_order);
      break;
    default:
      amrex::Abort("Error: Advecting quantity is unrecognized");
    }
    break;
  default:
    amrex::Abort("Error: Advected quantity is unrecognized");
  }

  // Return the advected quantity
  return advectedQty;
}

AMREX_GPU_DEVICE
Real
AdvectionContributionForState(const int &i, const int &j, const int &k,
                              const Array4<Real>& rho_u, const Array4<Real>& rho_v, const Array4<Real>& rho_w,
                              const Array4<Real>& cell_data, const int &qty_index,
                              const Array4<Real>& xflux, const Array4<Real>& yflux, const Array4<Real>& zflux,
                              const GpuArray<Real, AMREX_SPACEDIM>& cellSizeInv,
                              const int &spatial_order) {

    auto dxInv = cellSizeInv[0], dyInv = cellSizeInv[1], dzInv = cellSizeInv[2];

    Real advectionContribution;

    if (qty_index == Rho_comp)
    {
        advectionContribution = (rho_u(i+1,j,k,qty_index) - rho_u(i  ,j,k,qty_index)) * dxInv
                              + (rho_v(i,j+1,k,qty_index) - rho_v(i,j  ,k,qty_index)) * dyInv
                              + (rho_w(i,j,k+1,qty_index) - rho_w(i,j,k  ,qty_index)) * dzInv;
    } else {

        xflux(i+1,j,k,qty_index) = ComputeAdvectedQuantityForState(i+1, j, k, cell_data, qty_index,
                                   AdvectingQuantity::rho_u, spatial_order) * rho_u(i+1,j,k);
        xflux(i  ,j,k,qty_index) = ComputeAdvectedQuantityForState(i  , j, k, cell_data, qty_index,
                                   AdvectingQuantity::rho_u, spatial_order) * rho_u(i  ,j,k);
        yflux(i,j+1,k,qty_index) = ComputeAdvectedQuantityForState(i, j+1, k, cell_data, qty_index,
                                   AdvectingQuantity::rho_v, spatial_order) * rho_v(i,j+1,k);
        yflux(i,j  ,k,qty_index) = ComputeAdvectedQuantityForState(i, j  , k, cell_data, qty_index,
                                   AdvectingQuantity::rho_v, spatial_order) * rho_v(i,j  ,k);
        zflux(i,j,k+1,qty_index) = ComputeAdvectedQuantityForState(i, j, k+1, cell_data, qty_index,
                                   AdvectingQuantity::rho_w, spatial_order) * rho_w(i,j,k+1);
        zflux(i,j,k  ,qty_index) = ComputeAdvectedQuantityForState(i, j, k  , cell_data, qty_index,
                                   AdvectingQuantity::rho_w, spatial_order) * rho_w(i,j,k  );

        advectionContribution = (xflux(i+1,j,k,qty_index) - xflux(i  ,j,k,qty_index)) * dxInv
                              + (yflux(i,j+1,k,qty_index) - yflux(i,j  ,k,qty_index)) * dyInv
                              + (zflux(i,j,k+1,qty_index) - zflux(i,j,k  ,qty_index)) * dzInv;
    }

    return advectionContribution;
}
