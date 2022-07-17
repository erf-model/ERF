#include <SpatialStencils.H>

using namespace amrex;

AMREX_GPU_DEVICE
Real
interpolatedVal(const Real& avg1, const Real& avg2, const Real& avg3,
                const Real& diff1, const Real& diff2, const Real& diff3,
                const Real& scaled_upw, const int& spatial_order)
{
    Real myInterpolatedVal;
    switch (spatial_order) {
        case 2:
            myInterpolatedVal = 0.5 * avg1;
            break;
        case 3:
            myInterpolatedVal = (7.0/12.0)*avg1 -(1.0/12.0)*avg2 + (scaled_upw/12.0)*(diff2 - 3.0*diff1);
            break;
        case 4:
            myInterpolatedVal = (7.0/12.0)*avg1 -(1.0/12.0)*avg2;
            break;
        case 5:
            myInterpolatedVal = (37.0/60.0)*avg1 -(2.0/15.0)*avg2 +(1.0/60.0)*avg3
                              -(scaled_upw/60.0)*(diff3 - 5.0*diff2 + 10.0*diff1);
            break;
        case 6:
            myInterpolatedVal = (37.0/60.0)*avg1 -(2.0/15.0)*avg2 +(1.0/60.0)*avg3;
            break;
    }
    return myInterpolatedVal;
}

AMREX_GPU_DEVICE
Real
InterpolateDensityPertFromCellToFace(
  const int& i,
  const int& j,
  const int& k,
  const Array4<const Real>& cons_in,
  const Real& upw,
  const Coord& coordDir,
  const int& spatial_order,
  const Array4<const Real>& r0_arr)
{
  return InterpolatePertFromCell(
    i, j, k, cons_in, Rho_comp, upw, coordDir, spatial_order, r0_arr);
}

AMREX_GPU_DEVICE
Real
InterpolateFromCellOrFace(
  const int& i, const int& j, const int& k,
  const Array4<const Real>& qty,
  const int& qty_index,
  const Real& upw,
  const Coord& coordDir,
  const int& spatial_order)
{
    Real avg1 = 0.; Real avg2 = 0.; Real avg3 = 0.;
    Real diff1 = 0.; Real diff2 = 0.; Real diff3 = 0.;
    Real scaled_upw = 0.;

    // The value that comes in has not been normalized so we do that here
    if (upw != 0.)
        scaled_upw = upw / std::abs(upw);

    if (coordDir ==  Coord::x) {
        avg1  = (qty(i, j, k, qty_index) + qty(i-1, j, k, qty_index));
        diff1 = (qty(i, j, k, qty_index) - qty(i-1, j, k, qty_index));
        if (spatial_order > 2)
        {
            avg2  = (qty(i+1, j, k, qty_index) + qty(i-2, j, k, qty_index));
            diff2 = (qty(i+1, j, k, qty_index) - qty(i-2, j, k, qty_index));
        }
        if (spatial_order > 4)
        {
            avg3  = (qty(i+2, j, k, qty_index) + qty(i-3, j, k, qty_index));
            diff3 = (qty(i+2, j, k, qty_index) - qty(i-3, j, k, qty_index));
        }
    } else if (coordDir ==  Coord::y) {
        avg1  = (qty(i, j  , k, qty_index) + qty(i, j-1, k, qty_index));
        diff1 = (qty(i, j  , k, qty_index) - qty(i, j-1, k, qty_index));
        if (spatial_order > 2)
        {
            avg2  = (qty(i, j+1, k, qty_index) + qty(i, j-2, k, qty_index));
            diff2 = (qty(i, j+1, k, qty_index) - qty(i, j-2, k, qty_index));
        }
        if (spatial_order > 4)
        {
            avg3  = (qty(i, j+2, k, qty_index) + qty(i, j-3, k, qty_index));
            diff3 = (qty(i, j+2, k, qty_index) - qty(i, j-3, k, qty_index));
        }
    } else {
        avg1  = (qty(i, j, k  , qty_index) + qty(i, j, k-1, qty_index));
        diff1 = (qty(i, j, k  , qty_index) - qty(i, j, k-1, qty_index));
        if (spatial_order > 2)
        {
            avg2  = (qty(i, j, k+1, qty_index) + qty(i, j, k-2, qty_index));
            diff2 = (qty(i, j, k+1, qty_index) - qty(i, j, k-2, qty_index));
        }
        if (spatial_order > 4)
        {
            avg3  = (qty(i, j, k+2, qty_index) + qty(i, j, k-3, qty_index));
            diff3 = (qty(i, j, k+2, qty_index) - qty(i, j, k-3, qty_index));
        }
    }

    return interpolatedVal(avg1,avg2,avg3,diff1,diff2,diff3,scaled_upw,spatial_order);
}


AMREX_GPU_DEVICE
Real
InterpolatePertFromCell(
  const int& i, const int& j, const int& k,
  const Array4<const Real>& qty,
  const int& qty_index,
  const Real& upw,
  const Coord& coordDir,
  const int& spatial_order,
  const Array4<const Real>& r0_arr)
{
    Real avg1 = 0.; Real avg2 = 0.; Real avg3 = 0.;
    Real diff1 = 0.; Real diff2 = 0.; Real diff3 = 0.;
    Real scaled_upw = 0.;

    // The value that comes in has not been normalized so we do that here
    if (upw != 0.)
        scaled_upw = upw / std::abs(upw);

    if (coordDir ==  Coord::x) {
        avg1  = (qty(i  , j, k, qty_index) + qty(i-1, j, k, qty_index));
        avg1 -= (r0_arr(i,j,k) + r0_arr(i-1,j,k));
        diff1 = (qty(i  , j, k, qty_index) - qty(i-1, j, k, qty_index));
        if (spatial_order > 2)
        {
            avg2  = (qty(i+1, j, k, qty_index) + qty(i-2, j, k, qty_index));
            avg2 -= (r0_arr(i+1,j,k) + r0_arr(i-2,j,k));
            diff2 = (qty(i+1, j, k, qty_index) - qty(i-2, j, k, qty_index));
        }
        if (spatial_order > 4)
        {
            avg3  = (qty(i+2, j, k, qty_index) + qty(i-3, j, k, qty_index));
            avg3 -= (r0_arr(i+2,j,k) + r0_arr(i-3,j,k));
            diff3 = (qty(i+2, j, k, qty_index) - qty(i-3, j, k, qty_index));
        }
    } else if (coordDir ==  Coord::y) {
        avg1  = (qty(i, j  , k, qty_index) + qty(i, j-1, k, qty_index));
        avg1 -= (r0_arr(i,j,k) + r0_arr(i,j-1,k));
        diff1 = (qty(i, j  , k, qty_index) - qty(i, j-1, k, qty_index));
        if (spatial_order > 2)
        {
            avg2  = (qty(i, j+1, k, qty_index) + qty(i, j-2, k, qty_index));
            avg2 -= (r0_arr(i,j+1,k) + r0_arr(i,j-2,k));
            diff2 = (qty(i, j+1, k, qty_index) - qty(i, j-2, k, qty_index));
        }
        if (spatial_order > 4)
        {
            avg3  = (qty(i, j+2, k, qty_index) + qty(i, j-3, k, qty_index));
            avg3 -= (r0_arr(i,j+2,k) + r0_arr(i,j-3,k));
            diff3 = (qty(i, j+2, k, qty_index) - qty(i, j-3, k, qty_index));
        }
    } else {
        avg1  = (qty(i, j, k  , qty_index) + qty(i, j, k-1, qty_index));
        diff1 = (qty(i, j, k  , qty_index) - qty(i, j, k-1, qty_index));
        avg1  -= (r0_arr(i,j,k) + r0_arr(i,j,k-1));
        diff1 -= (r0_arr(i,j,k) - r0_arr(i,j,k-1));
        if (spatial_order > 2)
        {
            avg2  = (qty(i, j, k+1, qty_index) + qty(i, j, k-2, qty_index));
            diff2 = (qty(i, j, k+1, qty_index) - qty(i, j, k-2, qty_index));
            avg2  -= (r0_arr(i,j,k+1) + r0_arr(i,j,k-2));
            diff2 -= (r0_arr(i,j,k+1) - r0_arr(i,j,k-2));
        }
        if (spatial_order > 4)
        {
            avg3  = (qty(i, j, k+2, qty_index) + qty(i, j, k-3, qty_index));
            diff3 = (qty(i, j, k+2, qty_index) - qty(i, j, k-3, qty_index));
            avg3  -= (r0_arr(i,j,k+2) + r0_arr(i,j,k-3));
            diff3 -= (r0_arr(i,j,k+2) - r0_arr(i,j,k-3));
        }
    }

    return interpolatedVal(avg1,avg2,avg3,diff1,diff2,diff3,scaled_upw,spatial_order);
}

