#ifndef ERF_INTERPOLATE_H_
#define ERF_INTERPOLATE_H_

#include "ERF_DataStruct.H"
#include "ERF_Interpolation_UPW.H"
#include "ERF_Interpolation_WENO.H"
#include "ERF_Interpolation_WENO_Z.H"

/**
 * Interpolation operators used in construction of advective fluxes using non-WENO schemes
 */

AMREX_GPU_DEVICE
AMREX_FORCE_INLINE
amrex::Real
interpolatedVal (amrex::Real avg1,  amrex::Real avg2,  amrex::Real avg3,
                 amrex::Real diff1, amrex::Real diff2, amrex::Real diff3,
                 amrex::Real scaled_upw, const AdvType adv_type)
{
    amrex::Real myInterpolatedVal(amrex::Real(0));
    if (adv_type == AdvType::Centered_2nd) {
        myInterpolatedVal = myhalf * avg1;
    } else if (adv_type == AdvType::Upwind_3rd) {
        myInterpolatedVal = (amrex::Real(7.0)/amrex::Real(12.0))*avg1 -(amrex::Real(1)/amrex::Real(12.0))*avg2 + (scaled_upw/amrex::Real(12.0))*(diff2 - amrex::Real(3)*diff1);
    } else if (adv_type == AdvType::Centered_4th) {
        myInterpolatedVal = (amrex::Real(7.0)/amrex::Real(12.0))*avg1 -(amrex::Real(1)/amrex::Real(12.0))*avg2;
    } else if (adv_type == AdvType::Upwind_5th) {
        myInterpolatedVal = (amrex::Real(37.0)/amrex::Real(60.0))*avg1 -(amrex::Real(2)/amrex::Real(15.0))*avg2 +(amrex::Real(1)/amrex::Real(60.0))*avg3
                          -(scaled_upw/amrex::Real(60.0))*(diff3 - amrex::Real(5.0)*diff2 + amrex::Real(10.0)*diff1);
    } else if (adv_type == AdvType::Centered_6th) {
        myInterpolatedVal = (amrex::Real(37.0)/amrex::Real(60.0))*avg1 -(amrex::Real(2)/amrex::Real(15.0))*avg2 +(amrex::Real(1)/amrex::Real(60.0))*avg3;
    }
    return myInterpolatedVal;
}

AMREX_GPU_DEVICE
AMREX_FORCE_INLINE
amrex::Real
InterpolateInX (int i, int j, int k, const amrex::Array4<const amrex::Real>& qty,
                int qty_index, amrex::Real upw, const AdvType adv_type)
{

    if (adv_type == AdvType::Centered_2nd) {
        return myhalf * (qty(i,j,k,qty_index) + qty(i-1,j,k,qty_index));
    } else {

        amrex::Real avg1 = amrex::Real(0); amrex::Real avg2 = amrex::Real(0); amrex::Real avg3 = amrex::Real(0);
        amrex::Real diff1 = amrex::Real(0); amrex::Real diff2 = amrex::Real(0); amrex::Real diff3 = amrex::Real(0);
        amrex::Real scaled_upw = amrex::Real(0);
        //
        // The value that comes in has not been normalized so we do that here
        if (upw != amrex::Real(0)) { scaled_upw = (upw > 0) ? amrex::Real(1) : -amrex::Real(1); }

        avg1  = (qty(i, j, k, qty_index) + qty(i-1, j, k, qty_index));
        diff1 = (qty(i, j, k, qty_index) - qty(i-1, j, k, qty_index));
        avg2  = (qty(i+1, j, k, qty_index) + qty(i-2, j, k, qty_index));
        diff2 = (qty(i+1, j, k, qty_index) - qty(i-2, j, k, qty_index));
        if (adv_type == AdvType::Upwind_5th || adv_type == AdvType::Centered_6th)
        {
            avg3  = (qty(i+2, j, k, qty_index) + qty(i-3, j, k, qty_index));
            diff3 = (qty(i+2, j, k, qty_index) - qty(i-3, j, k, qty_index));
        }
        return interpolatedVal(avg1,avg2,avg3,diff1,diff2,diff3,scaled_upw,adv_type);
    }
}

AMREX_GPU_DEVICE
AMREX_FORCE_INLINE
amrex::Real
InterpolateInY (int i, int j, int k, const amrex::Array4<const amrex::Real>& qty,
                int qty_index, amrex::Real upw, const AdvType adv_type)
{
    if (adv_type == AdvType::Centered_2nd) {
        return myhalf * (qty(i,j,k,qty_index) + qty(i,j-1,k,qty_index));
    } else {

        amrex::Real avg1; amrex::Real avg2; amrex::Real avg3 = amrex::Real(0);
        amrex::Real diff1; amrex::Real diff2; amrex::Real diff3 = amrex::Real(0);
        amrex::Real scaled_upw = amrex::Real(0);

        // The value that comes in has not been normalized so we do that here
        if (upw != amrex::Real(0)) { scaled_upw = (upw > 0) ? amrex::Real(1) : -amrex::Real(1); }

        avg1  = (qty(i, j  , k, qty_index) + qty(i, j-1, k, qty_index));
        diff1 = (qty(i, j  , k, qty_index) - qty(i, j-1, k, qty_index));
        avg2  = (qty(i, j+1, k, qty_index) + qty(i, j-2, k, qty_index));
        diff2 = (qty(i, j+1, k, qty_index) - qty(i, j-2, k, qty_index));
        if (adv_type == AdvType::Upwind_5th || adv_type == AdvType::Centered_6th)
        {
            avg3  = (qty(i, j+2, k, qty_index) + qty(i, j-3, k, qty_index));
            diff3 = (qty(i, j+2, k, qty_index) - qty(i, j-3, k, qty_index));
        }
        return interpolatedVal(avg1,avg2,avg3,diff1,diff2,diff3,scaled_upw,adv_type);
    }
}

AMREX_GPU_DEVICE
AMREX_FORCE_INLINE
amrex::Real
InterpolateInZ (int i, int j, int k, const amrex::Array4<const amrex::Real>& qty,
                int qty_index, amrex::Real upw, const AdvType adv_type)
{
    if (adv_type == AdvType::Centered_2nd) {
        return myhalf * (qty(i,j,k,qty_index) + qty(i,j,k-1,qty_index));
    } else {

        amrex::Real avg1 = amrex::Real(0); amrex::Real avg2 = amrex::Real(0); amrex::Real avg3 = amrex::Real(0);
        amrex::Real diff1 = amrex::Real(0); amrex::Real diff2 = amrex::Real(0); amrex::Real diff3 = amrex::Real(0);
        amrex::Real scaled_upw = amrex::Real(0);
        // The value that comes in has not been normalized so we do that here
        if (upw != amrex::Real(0)) { scaled_upw = (upw > 0) ? amrex::Real(1) : -amrex::Real(1); }

        avg1  = (qty(i, j, k  , qty_index) + qty(i, j, k-1, qty_index));
        diff1 = (qty(i, j, k  , qty_index) - qty(i, j, k-1, qty_index));
        avg2  = (qty(i, j, k+1, qty_index) + qty(i, j, k-2, qty_index));
        diff2 = (qty(i, j, k+1, qty_index) - qty(i, j, k-2, qty_index));
        if (adv_type == AdvType::Upwind_5th || adv_type == AdvType::Centered_6th)
        {
            avg3  = (qty(i, j, k+2, qty_index) + qty(i, j, k-3, qty_index));
            diff3 = (qty(i, j, k+2, qty_index) - qty(i, j, k-3, qty_index));
        }
        return interpolatedVal(avg1,avg2,avg3,diff1,diff2,diff3,scaled_upw,adv_type);
    }
}

AMREX_GPU_DEVICE
AMREX_FORCE_INLINE
amrex::Real
InterpolatePertFromCell (int i, int j, int k,
                         const amrex::Array4<const amrex::Real>& qty,
                         int qty_index, amrex::Real upw, Coord coordDir,
                         const AdvType adv_type, const amrex::Array4<const amrex::Real>& r0_arr)
{
    amrex::Real avg1 = amrex::Real(0); amrex::Real avg2 = amrex::Real(0); amrex::Real avg3 = amrex::Real(0);
    amrex::Real diff1 = amrex::Real(0); amrex::Real diff2 = amrex::Real(0); amrex::Real diff3 = amrex::Real(0);
    amrex::Real scaled_upw = amrex::Real(0);

    // The value that comes in has not been normalized so we do that here
    if (upw != amrex::Real(0)) { scaled_upw = (upw > 0) ? amrex::Real(1) : -amrex::Real(1); }

    if (coordDir == Coord::x) {
        avg1  = (qty(i  , j, k, qty_index) + qty(i-1, j, k, qty_index));
        avg1 -= (r0_arr(i,j,k) + r0_arr(i-1,j,k));
        diff1 = (qty(i  , j, k, qty_index) - qty(i-1, j, k, qty_index));
        if (adv_type != AdvType::Centered_2nd)
        {
            avg2  = (qty(i+1, j, k, qty_index) + qty(i-2, j, k, qty_index));
            avg2 -= (r0_arr(i+1,j,k) + r0_arr(i-2,j,k));
            diff2 = (qty(i+1, j, k, qty_index) - qty(i-2, j, k, qty_index));
        }
        if (adv_type == AdvType::Upwind_5th || adv_type == AdvType::Centered_6th)
        {
            avg3  = (qty(i+2, j, k, qty_index) + qty(i-3, j, k, qty_index));
            avg3 -= (r0_arr(i+2,j,k) + r0_arr(i-3,j,k));
            diff3 = (qty(i+2, j, k, qty_index) - qty(i-3, j, k, qty_index));
        }
    } else if (coordDir == Coord::y) {
        avg1  = (qty(i, j  , k, qty_index) + qty(i, j-1, k, qty_index));
        avg1 -= (r0_arr(i,j,k) + r0_arr(i,j-1,k));
        diff1 = (qty(i, j  , k, qty_index) - qty(i, j-1, k, qty_index));
        if (adv_type != AdvType::Centered_2nd)
        {
            avg2  = (qty(i, j+1, k, qty_index) + qty(i, j-2, k, qty_index));
            avg2 -= (r0_arr(i,j+1,k) + r0_arr(i,j-2,k));
            diff2 = (qty(i, j+1, k, qty_index) - qty(i, j-2, k, qty_index));
        }
        if (adv_type == AdvType::Upwind_5th || adv_type == AdvType::Centered_6th)
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

        if (adv_type != AdvType::Centered_2nd)
        {
            avg2  = (qty(i, j, k+1, qty_index) + qty(i, j, k-2, qty_index));
            diff2 = (qty(i, j, k+1, qty_index) - qty(i, j, k-2, qty_index));
            avg2  -= (r0_arr(i,j,k+1) + r0_arr(i,j,k-2));
            diff2 -= (r0_arr(i,j,k+1) - r0_arr(i,j,k-2));
        }
        if (adv_type == AdvType::Upwind_5th || adv_type == AdvType::Centered_6th)
        {
            avg3  = (qty(i, j, k+2, qty_index) + qty(i, j, k-3, qty_index));
            diff3 = (qty(i, j, k+2, qty_index) - qty(i, j, k-3, qty_index));
            avg3  -= (r0_arr(i,j,k+2) + r0_arr(i,j,k-3));
            diff3 -= (r0_arr(i,j,k+2) - r0_arr(i,j,k-3));
        }
    }

    return interpolatedVal(avg1,avg2,avg3,diff1,diff2,diff3,scaled_upw,adv_type);
}

AMREX_GPU_DEVICE
AMREX_FORCE_INLINE
amrex::Real
InterpolateDensityPertFromCellToFace (int i, int j, int k, const amrex::Array4<const amrex::Real>& cons_in,
                                      amrex::Real upw, Coord coordDir, const AdvType adv_type,
                                      const amrex::Array4<const amrex::Real>& r0_arr)
{
    return InterpolatePertFromCell(i, j, k, cons_in, Rho_comp, upw, coordDir, adv_type, r0_arr);
}
#endif
