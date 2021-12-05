#include <TimeIntegration.H>

using namespace amrex;

AMREX_GPU_DEVICE
Real
InterpolateDensityPertFromCellToFace(
  const int& i,
  const int& j,
  const int& k,
  const Array4<const Real>& cons_in,
  const Coord& coordDir,
  const int& spatial_order,
  const amrex::Real* dptr_hse)
{
  return InterpolatePertFromCell(
    i, j, k, cons_in, Rho_comp, coordDir, spatial_order, dptr_hse);
}

AMREX_GPU_DEVICE
Real
InterpolateFromCellOrFace(
  // (i, j, k) is the reference cell index w.r.t. which a face is being considered
  // The same interpolation routine also works for face-based data and should work for
  // edge-based data
  const int& i, const int& j, const int& k,
  const Array4<const Real>& qty,
  const int& qty_index,
  const Coord& coordDir,
  const int& spatial_order)
{
  Real interpolatedVal = 0.0;
    /*
     w.r.t. the cell index (i, j, k), the face previous to it is the face at cell index m-1/2, where m = {i, j, k}
     Face index is (i, j, k). This means:
     Coordinates of face (i, j, k) = Coordinates of cell (i-1/2, j    , k    ) for x-dir.
     Coordinates of face (i, j, k) = Coordinates of cell (i    , j-1/2, k    ) for y-dir.
     Coordinates of face (i, j, k) = Coordinates of cell (i    , j    , k-1/2) for z-dir.
    */
    /*
     The above explanation is for cell- and face-index relation. However, the interpolation is applicable to any multifab data.
    */
    switch (spatial_order) {
    case 2:
      switch (coordDir) {
        // q = qty(i, j, k, qty_index) = {rho, theta, rhoTheta, scalar, pressure, ...}
      case Coord::x: // m = i, q(m-1/2) = q(i-1/2, j    , k    )
        interpolatedVal = 0.5*(qty(i, j, k, qty_index) + qty(i-1, j, k, qty_index));
        break;
      case Coord::y: // m = j, q(m-1/2) = q(i    , j-1/2, k    )
        interpolatedVal = 0.5*(qty(i, j, k, qty_index) + qty(i, j-1, k, qty_index));
        break;
      case Coord::z: // m = k, q(m-1/2) = q(i    , j    , k-1/2)
        interpolatedVal = 0.5*(qty(i, j, k, qty_index) + qty(i, j, k-1, qty_index));
        break;
      default:
        amrex::Abort("Error: Advection direction is unrecognized");
      }
      break;
    case 4:
      switch (coordDir) {
        // q = qty(i, j, k, qty_index) = {rho, theta, rhoTheta, scalar, pressure, ...}
      case Coord::x: // m = i, q(m-1/2) = q(i-1/2, j    , k    )
        interpolatedVal = (7.0/12.0)*(qty(i, j, k, qty_index) + qty(i-1, j, k, qty_index))
                          -(1.0/12.0)*(qty(i+1, j, k, qty_index) + qty(i-2, j, k, qty_index));
        break;
      case Coord::y: // m = j, q(m-1/2) = q(i    , j-1/2, k    )
        interpolatedVal = (7.0/12.0)*(qty(i, j, k, qty_index) + qty(i, j-1, k, qty_index))
                          -(1.0/12.0)*(qty(i, j+1, k, qty_index) + qty(i, j-2, k, qty_index));
        break;
      case Coord::z: // m = k, q(m-1/2) = q(i    , j    , k-1/2)
        interpolatedVal = (7.0/12.0)*(qty(i, j, k, qty_index) + qty(i, j, k-1, qty_index))
                          -(1.0/12.0)*(qty(i, j, k+1, qty_index) + qty(i, j, k-2, qty_index));
        break;
      default:
        amrex::Abort("Error: Advection direction is unrecognized");
      }
      break;
    case 6: // In order to make this work 'qty' must have indices m-3 and m+2 where m = {i, j, k}
      switch (coordDir) {
        // q = qty(i, j, k, qty_index) = {rho, theta, rhoTheta, scalar, pressure, ...}
      case Coord::x: // m = i, q(m-1/2) = q(i-1/2, j    , k    )
        interpolatedVal = (37.0/60.0)*(qty(i, j, k, qty_index) + qty(i-1, j, k, qty_index))
                          -(2.0/15.0)*(qty(i+1, j, k, qty_index) + qty(i-2, j, k, qty_index))
                          +(1.0/60.0)*(qty(i+2, j, k, qty_index) + qty(i-3, j, k, qty_index));
        break;
      case Coord::y: // m = j, q(m-1/2) = q(i    , j-1/2, k    )
        interpolatedVal = (37.0/60.0)*(qty(i, j, k, qty_index) + qty(i, j-1, k, qty_index))
                          -(2.0/15.0)*(qty(i, j+1, k, qty_index) + qty(i, j-2, k, qty_index))
                          +(1.0/60.0)*(qty(i, j+2, k, qty_index) + qty(i, j-3, k, qty_index));
        break;
      case Coord::z: // m = k, q(m-1/2) = q(i    , j    , k-1/2)
        interpolatedVal = (37.0/60.0)*(qty(i, j, k, qty_index) + qty(i, j, k-1, qty_index))
                          -(2.0/15.0)*(qty(i, j, k+1, qty_index) + qty(i, j, k-2, qty_index))
                          +(1.0/60.0)*(qty(i, j, k+2, qty_index) + qty(i, j, k-3, qty_index));
        break;
      default:
        amrex::Abort("Error: Advection direction is unrecognized");
      }
      break;
    default:
      amrex::Abort("Error: Spatial order has not been implemented");
    }

  return interpolatedVal;
}


AMREX_GPU_DEVICE
Real
InterpolatePertFromCell(
  const int& i, const int& j, const int& k,
  const Array4<const Real>& qty,
  const int& qty_index,
  const Coord& coordDir,
  const int& spatial_order,
  const amrex::Real* dptr_hse)
{
  Real interpolatedVal = 0.0;
    /*
     w.r.t. the cell index (i, j, k), the face previous to it is the face at cell index m-1/2, where m = {i, j, k}
     Face index is (i, j, k). This means:
     Coordinates of face (i, j, k) = Coordinates of cell (i-1/2, j    , k    ) for x-dir. Face is previous to the cell.
     Coordinates of face (i, j, k) = Coordinates of cell (i    , j-1/2, k    ) for y-dir. Face is previous to the cell.
     Coordinates of face (i, j, k) = Coordinates of cell (i    , j    , k-1/2) for z-dir. Face is previous to the cell.
    */
    switch (spatial_order) {
    case 2:
      switch (coordDir) {
        // q = qty(i, j, k, qty_index) = {rho, theta, rhoTheta, scalar, pressure, ...}
      case Coord::x: // m = i, q(m-1/2) = q(i-1/2, j    , k    )
        interpolatedVal = 0.5*(qty(i, j, k, qty_index) + qty(i-1, j, k, qty_index));
        break;
      case Coord::y: // m = j, q(m-1/2) = q(i    , j-1/2, k    )
        interpolatedVal = 0.5*(qty(i, j, k, qty_index) + qty(i, j-1, k, qty_index));
        break;
      case Coord::z: // m = k, q(m-1/2) = q(i    , j    , k-1/2)
        interpolatedVal = 0.5*((qty(i, j, k  , qty_index)-dptr_hse[k  ])
                             + (qty(i, j, k-1, qty_index)-dptr_hse[k-1]));
        break;
      default:
        amrex::Abort("Error: Advection direction is unrecognized");
      }
      break;
    case 4:
      switch (coordDir) {
        // q = qty(i, j, k, qty_index) = {rho, theta, rhoTheta, scalar, pressure, ...}
      case Coord::x: // m = i, q(m-1/2) = q(i-1/2, j    , k    )
        interpolatedVal = (7.0/12.0)*(qty(i, j, k, qty_index) + qty(i-1, j, k, qty_index))
                          -(1.0/12.0)*(qty(i+1, j, k, qty_index) + qty(i-2, j, k, qty_index));
        break;
      case Coord::y: // m = j, q(m-1/2) = q(i    , j-1/2, k    )
        interpolatedVal = (7.0/12.0)*(qty(i, j, k, qty_index) + qty(i, j-1, k, qty_index))
                          -(1.0/12.0)*(qty(i, j+1, k, qty_index) + qty(i, j-2, k, qty_index));
        break;
      case Coord::z: // m = k, q(m-1/2) = q(i    , j    , k-1/2)
        interpolatedVal = (7.0/12.0)*((qty(i, j, k  , qty_index)-dptr_hse[k  ])
                                    + (qty(i, j, k-1, qty_index)-dptr_hse[k-1]))
                         -(1.0/12.0)*((qty(i, j, k+1, qty_index)-dptr_hse[k+1])
                                    + (qty(i, j, k-2, qty_index)-dptr_hse[k-2]));
        break;
      default:
        amrex::Abort("Error: Advection direction is unrecognized");
      }
      break;
    case 6: // In order to make this work 'qty' must have indices m-3 and m+2 where m = {i, j, k}
      switch (coordDir) {
        // q = qty(i, j, k, qty_index) = {rho, theta, rhoTheta, scalar, pressure, ...}
      case Coord::x: // m = i, q(m-1/2) = q(i-1/2, j    , k    )
        interpolatedVal = (37.0/60.0)*(qty(i  , j, k, qty_index) + qty(i-1, j, k, qty_index))
                          -(2.0/15.0)*(qty(i+1, j, k, qty_index) + qty(i-2, j, k, qty_index))
                          +(1.0/60.0)*(qty(i+2, j, k, qty_index) + qty(i-3, j, k, qty_index));
        break;
      case Coord::y: // m = j, q(m-1/2) = q(i    , j-1/2, k    )
        interpolatedVal = (37.0/60.0)*(qty(i, j  , k, qty_index) + qty(i, j-1, k, qty_index))
                          -(2.0/15.0)*(qty(i, j+1, k, qty_index) + qty(i, j-2, k, qty_index))
                          +(1.0/60.0)*(qty(i, j+2, k, qty_index) + qty(i, j-3, k, qty_index));
        break;
      case Coord::z: // m = k, q(m-1/2) = q(i    , j    , k-1/2)
        interpolatedVal = (37.0/60.0)*((qty(i, j, k  , qty_index)-dptr_hse[k  ])
                                     + (qty(i, j, k-1, qty_index)-dptr_hse[k-1]))
                          -(2.0/15.0)*((qty(i, j, k+1, qty_index)-dptr_hse[k+1])
                                     + (qty(i, j, k-2, qty_index)-dptr_hse[k-2]))
                          +(1.0/60.0)*((qty(i, j, k+2, qty_index)-dptr_hse[k+2])
                                     + (qty(i, j, k-3, qty_index)-dptr_hse[k-3]));
        break;
      default:
        amrex::Abort("Error: Advection direction is unrecognized");
      }
      break;
    default:
      amrex::Abort("Error: Spatial order has not been implemented");
  }

  return interpolatedVal;
}
