#include <TimeIntegration.H>

using namespace amrex;

AMREX_GPU_DEVICE
Real
ComputeAdvectedQuantityForMom(const int &i, const int &j, const int &k,
                              const Array4<const Real>& rho_u, const Array4<const Real>& rho_v, const Array4<const Real>& rho_w,
                              const Array4<Real>& u, const Array4<Real>& v, const Array4<Real>& w,
                              const enum AdvectedQuantity &advectedQuantity,
                              const enum AdvectingQuantity &advectingQuantity,
                              const int &spatial_order)
{
  Real flux = 0.;
  switch(advectedQuantity) {
    case AdvectedQuantity::u: //x-momentum
      switch (advectingQuantity) {
      case AdvectingQuantity::rho_u:
        flux = 0.5*(rho_u(i-1, j, k) + rho_u(i, j, k)) *
                    InterpolateFromCellOrFace(i, j, k, u, 0, Coord::x, spatial_order);
        break;
      case AdvectingQuantity::rho_v:
        flux = 0.5*(rho_v(i, j, k) + rho_v(i-1, j, k)) *
                    InterpolateFromCellOrFace(i, j, k, u, 0, Coord::y, spatial_order);
        break;
      case AdvectingQuantity::rho_w:
        flux = 0.5*(rho_w(i, j, k) + rho_w(i-1, j, k)) *
                    InterpolateFromCellOrFace(i, j, k, u, 0, Coord::z, spatial_order);
        break;
      default:
        amrex::Abort("Error: Advecting quantity is unrecognized");
      }
      break;

    case AdvectedQuantity::v: //y-momentum
      switch (advectingQuantity) {
      case AdvectingQuantity::rho_u:
        flux = 0.5*(rho_u(i, j, k) + rho_u(i, j-1, k)) *
                    InterpolateFromCellOrFace(i, j, k, v, 0, Coord::x, spatial_order);
        break;
      case AdvectingQuantity::rho_v:
        flux = 0.5*(rho_v(i, j, k) + rho_v(i, j-1, k)) *
                    InterpolateFromCellOrFace(i, j, k, v, 0, Coord::y, spatial_order);
        break;
      case AdvectingQuantity::rho_w:
        flux = 0.5*(rho_w(i, j, k) + rho_w(i, j-1, k)) *
               InterpolateFromCellOrFace(i, j, k, v, 0, Coord::z, spatial_order);
        break;
      default:
        amrex::Abort("Error: Advecting quantity is unrecognized");
      }
      break;

    case AdvectedQuantity::w: //z-momentum
      switch (advectingQuantity) {
      case AdvectingQuantity::rho_u:
        flux = 0.5*(rho_u(i, j, k) + rho_u(i, j, k-1)) *
                    InterpolateFromCellOrFace(i, j, k, w, 0, Coord::x, spatial_order);
        break;
      case AdvectingQuantity::rho_v:
        flux = 0.5*(rho_v(i, j, k) + rho_v(i, j, k-1)) *
                    InterpolateFromCellOrFace(i, j, k, w, 0, Coord::y, spatial_order);
        break;
      case AdvectingQuantity::rho_w:
        flux = 0.5*(rho_w(i, j, k) + rho_w(i, j, k-1)) *
                    InterpolateFromCellOrFace(i, j, k, w, 0, Coord::z, spatial_order);
        break;
      default:
        amrex::Abort("Error: Advecting quantity is unrecognized");
      }
      break;
    default:
      amrex::Abort("Error: Advected quantity is unrecognized");
  }

  return flux;
}

AMREX_GPU_DEVICE
Real
AdvectionContributionForMom(const int &i, const int &j, const int &k,
                            const Array4<const Real>& rho_u, const Array4<const Real>& rho_v, const Array4<const Real>& rho_w,
                            const Array4<Real>& u, const Array4<Real>& v, const Array4<Real>& w,
                            const enum MomentumEqn &momentumEqn,
                            const GpuArray<Real, AMREX_SPACEDIM>& cellSizeInv,
                            const SolverChoice &solverChoice) {

    auto dxInv = cellSizeInv[0], dyInv = cellSizeInv[1], dzInv = cellSizeInv[2];
    Real advectionContribution = 0.0;

    int l_spatial_order = solverChoice.spatial_order;

    switch (momentumEqn) {
        case MomentumEqn::x: //x-momentum, reference face index is (i, j, k)
            Real centFluxXXNext, centFluxXXPrev, edgeFluxXYNext, edgeFluxXYPrev, edgeFluxXZNext, edgeFluxXZPrev;
            centFluxXXNext = ComputeAdvectedQuantityForMom(i+1, j, k, rho_u, rho_v, rho_w, u, v, w,
                AdvectedQuantity::u, AdvectingQuantity::rho_u, l_spatial_order);
            centFluxXXPrev = ComputeAdvectedQuantityForMom(i, j, k, rho_u, rho_v, rho_w, u, v, w,
                AdvectedQuantity::u, AdvectingQuantity::rho_u, l_spatial_order);
            edgeFluxXYNext = ComputeAdvectedQuantityForMom(i, j+1, k, rho_u, rho_v, rho_w, u, v, w,
                AdvectedQuantity::u, AdvectingQuantity::rho_v, l_spatial_order);
            edgeFluxXYPrev = ComputeAdvectedQuantityForMom(i, j, k, rho_u, rho_v, rho_w, u, v, w,
                AdvectedQuantity::u, AdvectingQuantity::rho_v, l_spatial_order);
            edgeFluxXZNext = ComputeAdvectedQuantityForMom(i, j, k+1, rho_u, rho_v, rho_w, u, v, w,
                AdvectedQuantity::u, AdvectingQuantity::rho_w, l_spatial_order);
            edgeFluxXZPrev = ComputeAdvectedQuantityForMom(i, j, k, rho_u, rho_v, rho_w, u, v, w,
                AdvectedQuantity::u, AdvectingQuantity::rho_w, l_spatial_order);

            advectionContribution = (centFluxXXNext - centFluxXXPrev) * dxInv   // Contribution to x-mom eqn from advective flux in x-dir
                                  + (edgeFluxXYNext - edgeFluxXYPrev) * dyInv   // Contribution to x-mom eqn from advective flux in y-dir
                                  + (edgeFluxXZNext - edgeFluxXZPrev) * dzInv;  // Contribution to x-mom eqn from advective flux in z-dir
            break;
        case MomentumEqn::y: //y-momentum, reference face index is (i, j, k)
            Real centFluxYYNext, centFluxYYPrev, edgeFluxYXNext, edgeFluxYXPrev, edgeFluxYZNext, edgeFluxYZPrev;
            edgeFluxYXNext = ComputeAdvectedQuantityForMom(i+1, j, k, rho_u, rho_v, rho_w, u, v, w,
                AdvectedQuantity::v, AdvectingQuantity::rho_u, l_spatial_order);
            edgeFluxYXPrev = ComputeAdvectedQuantityForMom(i, j, k, rho_u, rho_v, rho_w, u, v, w,
                AdvectedQuantity::v, AdvectingQuantity::rho_u, l_spatial_order);
            centFluxYYNext = ComputeAdvectedQuantityForMom(i, j+1, k, rho_u, rho_v, rho_w, u, v, w,
                AdvectedQuantity::v, AdvectingQuantity::rho_v, l_spatial_order);
            centFluxYYPrev = ComputeAdvectedQuantityForMom(i, j, k, rho_u, rho_v, rho_w, u, v, w,
                AdvectedQuantity::v, AdvectingQuantity::rho_v, l_spatial_order);
            edgeFluxYZNext = ComputeAdvectedQuantityForMom(i, j, k+1, rho_u, rho_v, rho_w, u, v, w,
                AdvectedQuantity::v, AdvectingQuantity::rho_w, l_spatial_order);
            edgeFluxYZPrev = ComputeAdvectedQuantityForMom(i, j, k, rho_u, rho_v, rho_w, u, v, w,
                AdvectedQuantity::v, AdvectingQuantity::rho_w, l_spatial_order);

            advectionContribution = (edgeFluxYXNext - edgeFluxYXPrev) * dxInv   // Contribution to y-mom eqn from advective flux in x-dir
                                  + (centFluxYYNext - centFluxYYPrev) * dyInv   // Contribution to y-mom eqn from advective flux in y-dir
                                  + (edgeFluxYZNext - edgeFluxYZPrev) * dzInv;  // Contribution to y-mom eqn from advective flux in z-dir
            break;
        case MomentumEqn::z: //z-momentum, reference face index is (i, j, k)
            Real centFluxZZNext, centFluxZZPrev, edgeFluxZXNext, edgeFluxZXPrev, edgeFluxZYNext, edgeFluxZYPrev;
            edgeFluxZXNext = ComputeAdvectedQuantityForMom(i+1, j, k, rho_u, rho_v, rho_w, u, v, w,
                AdvectedQuantity::w, AdvectingQuantity::rho_u, l_spatial_order);
            edgeFluxZXPrev = ComputeAdvectedQuantityForMom(i, j, k, rho_u, rho_v, rho_w, u, v, w,
                AdvectedQuantity::w, AdvectingQuantity::rho_u, l_spatial_order);
            edgeFluxZYNext = ComputeAdvectedQuantityForMom(i, j+1, k, rho_u, rho_v, rho_w, u, v, w,
                AdvectedQuantity::w, AdvectingQuantity::rho_v, l_spatial_order);
            edgeFluxZYPrev = ComputeAdvectedQuantityForMom(i, j, k, rho_u, rho_v, rho_w, u, v, w,
                AdvectedQuantity::w, AdvectingQuantity::rho_v, l_spatial_order);
            centFluxZZNext = ComputeAdvectedQuantityForMom(i, j, k+1, rho_u, rho_v, rho_w, u, v, w,
                AdvectedQuantity::w, AdvectingQuantity::rho_w, l_spatial_order);
            centFluxZZPrev = ComputeAdvectedQuantityForMom(i, j, k, rho_u, rho_v, rho_w, u, v, w,
                AdvectedQuantity::w, AdvectingQuantity::rho_w, l_spatial_order);

            advectionContribution = (edgeFluxZXNext - edgeFluxZXPrev) * dxInv  // Contribution to z-mom eqn from advective flux in x-dir
                                  + (edgeFluxZYNext - edgeFluxZYPrev) * dyInv  // Contribution to z-mom eqn from advective flux in y-dir
                                  + (centFluxZZNext - centFluxZZPrev) * dzInv; // Contribution to z-mom eqn from advective flux in z-dir
            break;
        default:
            amrex::Abort("Error: Momentum equation is unrecognized");
    }

    return advectionContribution;
}

AMREX_GPU_DEVICE
Real
ComputeAdvectedQuantityForState(const int &i, const int &j, const int &k,
                                const Array4<const Real>& cell_data, const int& qty_index,
                                const enum AdvectingQuantity &advectingQuantity,
                                const int &spatial_order) {
  Real advectedQty = 1.0;

  switch (advectingQuantity) {
    case AdvectingQuantity::rho_u:
      // Get theta (i-1/2,    j, k    ) = theta on face (i,   j  , k  ) for x-dir
      advectedQty = InterpolateFromCellOrFace(i, j, k, cell_data, qty_index, Coord::x, spatial_order);
      advectedQty/= InterpolateFromCellOrFace(i, j, k, cell_data, Rho_comp , Coord::x, spatial_order);
      break;
    case AdvectingQuantity::rho_v:
      // Get theta (i   , j-1/2, k    ) = theta on face (i  , j  , k  ) for y-dir
      advectedQty = InterpolateFromCellOrFace(i, j, k, cell_data, qty_index, Coord::y, spatial_order);
      advectedQty/= InterpolateFromCellOrFace(i, j, k, cell_data, Rho_comp , Coord::y, spatial_order);
      break;
    case AdvectingQuantity::rho_w:
      // Get theta (i   , j    , k-1/2) = theta on face (i  , j  , k  ) for z-dir
      advectedQty = InterpolateFromCellOrFace(i, j, k, cell_data, qty_index, Coord::z, spatial_order);
      advectedQty/= InterpolateFromCellOrFace(i, j, k, cell_data, Rho_comp , Coord::z, spatial_order);
      break;
    default:
      amrex::Abort("Error: Advecting quantity is unrecognized");
  }

  // Return the advected quantity
  return advectedQty;
}

AMREX_GPU_DEVICE
Real
AdvectionContributionForState(const int &i, const int &j, const int &k,
                              const Array4<const Real>& rho_u, const Array4<const Real>& rho_v, const Array4<const Real>& rho_w,
                              const Array4<const Real>& cell_data, const int &qty_index,
                              const Array4<Real>& xflux, const Array4<Real>& yflux, const Array4<Real>& zflux,
                              const GpuArray<Real, AMREX_SPACEDIM>& cellSizeInv,
                              const int &spatial_order) {

    auto dxInv = cellSizeInv[0], dyInv = cellSizeInv[1], dzInv = cellSizeInv[2];

    Real advectionContribution;

    if (qty_index == Rho_comp)
    {
        xflux(i+1,j,k,qty_index) = rho_u(i+1,j,k);
        xflux(i  ,j,k,qty_index) = rho_u(i  ,j,k);
        yflux(i,j+1,k,qty_index) = rho_v(i,j+1,k);
        yflux(i,j  ,k,qty_index) = rho_v(i,j  ,k);
        zflux(i,j,k+1,qty_index) = rho_w(i,j,k+1);
        zflux(i,j,k  ,qty_index) = rho_w(i,j,k  );

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
    }

    advectionContribution = (xflux(i+1,j,k,qty_index) - xflux(i  ,j,k,qty_index)) * dxInv
                          + (yflux(i,j+1,k,qty_index) - yflux(i,j  ,k,qty_index)) * dyInv
                          + (zflux(i,j,k+1,qty_index) - zflux(i,j,k  ,qty_index)) * dzInv;

    return advectionContribution;
}
