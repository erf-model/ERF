#ifndef INTERPOLATE_UPW_H_
#define INTERPOLATE_UPW_H_

#include "DataStruct.H"

/**
 * Interpolation operators used for 3rd order upwind scheme
 */
struct UPWIND3
{
    UPWIND3(const amrex::Array4<const amrex::Real>& phi,
            const amrex::Array4<const amrex::Real>& psi_x,
            const amrex::Array4<const amrex::Real>& psi_y,
            const amrex::Array4<const amrex::Real>& psi_z)
        : m_phi(phi), m_psi_x(psi_x), m_psi_y(psi_y), m_psi_z(psi_z) {}

    AMREX_GPU_DEVICE
    AMREX_FORCE_INLINE
    void
    InterpolateInX(const int& i,
                   const int& j,
                   const int& k,
                   const int& qty_index,
                   amrex::Real val_hi,
                   amrex::Real val_lo,
                   const amrex::Real& upw_hi_in,
                   const amrex::Real& upw_lo_in)
    {
        // Data to interpolate on
        amrex::Real sp2 = m_phi(i+2, j  , k  , qty_index);
        amrex::Real sp1 = m_phi(i+1, j  , k  , qty_index);
        amrex::Real s   = m_phi(i  , j  , k  , qty_index);
        amrex::Real sm1 = m_phi(i-1, j  , k  , qty_index);
        amrex::Real sm2 = m_phi(i-2, j  , k  , qty_index);

        // Upwinding flags
        amrex::Real upw_hi=0.0; amrex::Real upw_lo=0.0;
        if (upw_hi_in != 0.) upw_hi = (upw_hi_in > 0) ? 1. : -1.;
        if (upw_lo_in != 0.) upw_lo = (upw_lo_in > 0) ? 1. : -1.;

        // Interpolate hi
        val_hi = Evaluate(sp2,sp1,s,sm1,upw_hi);

        // Interpolate lo
        val_lo = Evaluate(sp1,s,sm1,sm2,upw_lo);
    }

    AMREX_GPU_DEVICE
    AMREX_FORCE_INLINE
    void
    InterpolateInY(const int& i,
                   const int& j,
                   const int& k,
                   const int& qty_index,
                   amrex::Real val_hi,
                   amrex::Real val_lo,
                   const amrex::Real& upw_hi_in,
                   const amrex::Real& upw_lo_in)
    {
        // Data to interpolate on
        amrex::Real sp2 = m_phi(i  , j+2, k  , qty_index);
        amrex::Real sp1 = m_phi(i  , j+1, k  , qty_index);
        amrex::Real s   = m_phi(i  , j  , k  , qty_index);
        amrex::Real sm1 = m_phi(i  , j-1, k  , qty_index);
        amrex::Real sm2 = m_phi(i  , j-2, k  , qty_index);

        // Upwinding flags
        amrex::Real upw_hi=0.0; amrex::Real upw_lo=0.0;
        if (upw_hi_in != 0.) upw_hi = (upw_hi_in > 0) ? 1. : -1.;
        if (upw_lo_in != 0.) upw_lo = (upw_lo_in > 0) ? 1. : -1.;

        // Interpolate hi
        val_hi = Evaluate(sp2,sp1,s,sm1,upw_hi);

        // Interpolate lo
        val_lo = Evaluate(sp1,s,sm1,sm2,upw_lo);
    }

    AMREX_GPU_DEVICE
    AMREX_FORCE_INLINE
    void
    InterpolateInZ(const int& i,
                   const int& j,
                   const int& k,
                   const int& qty_index,
                   amrex::Real val_hi,
                   amrex::Real val_lo,
                   const amrex::Real& upw_hi_in,
                   const amrex::Real& upw_lo_in)
    {
        // Data to interpolate on
        amrex::Real sp2 = m_phi(i  , j  , k+2, qty_index);
        amrex::Real sp1 = m_phi(i  , j  , k+1, qty_index);
        amrex::Real s   = m_phi(i  , j  , k  , qty_index);
        amrex::Real sm1 = m_phi(i  , j  , k-1, qty_index);
        amrex::Real sm2 = m_phi(i  , j  , k-2, qty_index);

        // Upwinding flags
        amrex::Real upw_hi=0.0; amrex::Real upw_lo=0.0;
        if (upw_hi_in != 0.) upw_hi = (upw_hi_in > 0) ? 1. : -1.;
        if (upw_lo_in != 0.) upw_lo = (upw_lo_in > 0) ? 1. : -1.;

        // Interpolate hi
        val_hi = Evaluate(sp2,sp1,s,sm1,upw_hi);

        // Interpolate lo
        val_lo = Evaluate(sp1,s,sm1,sm2,upw_lo);
    }

    AMREX_GPU_DEVICE
    AMREX_FORCE_INLINE
    amrex::Real
    Evaluate(const amrex::Real& sp1,
             const amrex::Real& s,
             const amrex::Real& sm1,
             const amrex::Real& sm2,
             const amrex::Real& upw)
    {
        // Averages and diffs
        a1 = (s + sm1);
        d1 = (s - sm1);
        a2 = (sp1 + sm2);
        d2 = (sp1 - sm2);

        // Interpolated value
        return ( g1*a1 - g2*a2 + upw*g2 * (d2 - 3.0*d1) );
    }

private:
    const amrex::Array4<const amrex::Real>& m_phi;   // Quantity to interpolate
    const amrex::Array4<const amrex::Real>& m_psi_x; // Upwind indicator
    const amrex::Array4<const amrex::Real>& m_psi_y; // Upwind indicator
    const amrex::Array4<const amrex::Real>& m_psi_z; // Upwind indicator
    amrex::Real a1 = 0.; amrex::Real a2 = 0.;
    amrex::Real d1 = 0.; amrex::Real d2 = 0.;
    static constexpr amrex::Real g1=(7.0/12.0);
    static constexpr amrex::Real g2=(1.0/12.0);
};

#if 0
/**
 * Interpolation operators used for 4th order upwind scheme
 */
struct UPWIND4
{
    UPWIND4(const amrex::Array4<const amrex::Real>& phi,
            const amrex::Array4<const amrex::Real>& psi)
        : m_phi(phi), m_psi(psi) {}

    AMREX_GPU_DEVICE
    AMREX_FORCE_INLINE
    void
    InterpolateInX(const int& i,
                   const int& j,
                   const int& k,
                   const int& qty_index,
                   amrex::Real val_hi,
                   amrex::Real val_lo,
                   amrex::Real /*upw_hi=0.0*/,
                   amrex::Real /*upw_lo=0.0*/)
    {
        // Data to interpolate on
        amrex::Real sp2 = m_phi(i+2, j  , k  , qty_index);
        amrex::Real sp1 = m_phi(i+1, j  , k  , qty_index);
        amrex::Real s   = m_phi(i  , j  , k  , qty_index);
        amrex::Real sm1 = m_phi(i-1, j  , k  , qty_index);
        amrex::Real sm2 = m_phi(i-2, j  , k  , qty_index);

        // Interpolate hi
        val_hi = Evaluate(sp2,sp1,s,sm1,upw_hi);

        // Interpolate lo
        val_lo = Evaluate(sp1,s,sm1,sm2,upw_lo);
    }

    AMREX_GPU_DEVICE
    AMREX_FORCE_INLINE
    void
    InterpolateInY(const int& i,
                   const int& j,
                   const int& k,
                   const int& qty_index,
                   amrex::Real val_hi,
                   amrex::Real val_lo,
                   amrex::Real /*upw_hi=0.0*/,
                   amrex::Real /*upw_lo=0.0*/)
    {
        // Data to interpolate on
        amrex::Real sp2 = m_phi(i  , j+2, k  , qty_index);
        amrex::Real sp1 = m_phi(i  , j+1, k  , qty_index);
        amrex::Real s   = m_phi(i  , j  , k  , qty_index);
        amrex::Real sm1 = m_phi(i  , j-1, k  , qty_index);
        amrex::Real sm2 = m_phi(i  , j-2, k  , qty_index);

        // Interpolate hi
        val_hi = Evaluate(sp2,sp1,s,sm1);

        // Interpolate lo
        val_lo = Evaluate(sp1,s,sm1,sm2);
    }

    AMREX_GPU_DEVICE
    AMREX_FORCE_INLINE
    void
    InterpolateInZ(const int& i,
                   const int& j,
                   const int& k,
                   const int& qty_index,
                   amrex::Real val_hi,
                   amrex::Real val_lo,
                   amrex::Real /*upw_hi=0.0*/,
                   amrex::Real /*upw_lo=0.0*/)
    {
        // Data to interpolate on
        amrex::Real sp2 = m_phi(i  , j  , k+2, qty_index);
        amrex::Real sp1 = m_phi(i  , j  , k+1, qty_index);
        amrex::Real s   = m_phi(i  , j  , k  , qty_index);
        amrex::Real sm1 = m_phi(i  , j  , k-1, qty_index);
        amrex::Real sm2 = m_phi(i  , j  , k-2, qty_index);

        // Interpolate hi
        val_hi = Evaluate(sp2,sp1,s,sm1);

        // Interpolate lo
        val_lo = Evaluate(sp1,s,sm1,sm2);
    }

    AMREX_GPU_DEVICE
    AMREX_FORCE_INLINE
    amrex::Real
    Evaluate(const amrex::Real& sp1,
             const amrex::Real& s,
             const amrex::Real& sm1,
             const amrex::Real& sm2)
    {
        // Averages and diffs
        a1 = (s + sm1);
        d1 = (s - sm1);
        a2 = (sp1 + sm2);
        d2 = (sp1 - sm2);

        // Interpolated value
        return ( g1*a1 - g2*a2 );
    }

private:
    const amrex::Array4<const amrex::Real>& m_phi; // Quantity to interpolate
    const amrex::Array4<const amrex::Real>& m_psi; // Upwind indicator
    amrex::Real a1 = 0.; amrex::Real a2 = 0.;
    amrex::Real d1 = 0.; amrex::Real d2 = 0.;
    constexpr amrex::Real g1=(7.0/12.0);
    constexpr amrex::Real g2=(1.0/12.0);
};

/**
 * Interpolation operators used for 5th order upwind scheme
 */
struct UPWIND5
{
    UPWIND5(const amrex::Array4<const amrex::Real>& phi,
            const amrex::Array4<const amrex::Real>& psi)
        : m_phi(phi), m_psi(psi) {}

    AMREX_GPU_DEVICE
    AMREX_FORCE_INLINE
    void
    InterpolateInX(const int& i,
                   const int& j,
                   const int& k,
                   const int& qty_index,
                   amrex::Real val_hi,
                   amrex::Real val_lo)
    {
        // Data to interpolate on
        amrex::Real sp3 = m_phi(i+3, j  , k  , qty_index);
        amrex::Real sp2 = m_phi(i+2, j  , k  , qty_index);
        amrex::Real sp1 = m_phi(i+1, j  , k  , qty_index);
        amrex::Real s   = m_phi(i  , j  , k  , qty_index);
        amrex::Real sm1 = m_phi(i-1, j  , k  , qty_index);
        amrex::Real sm2 = m_phi(i-2, j  , k  , qty_index);
        amrex::Real sm3 = m_phi(i-3, j  , k  , qty_index);

        // Upwinding flag
        upw = 0.5 * ( m_psi(i+1, j, k) + m_psi(i, j, k) );
        if (upw != 0.) upw = (upw > 0) ? 1. : -1.;

        // Interpolate hi
        val_hi = Evaluate(sp3,sp2,sp1,s,sm1,sm2,upw);

        // Upwinding flag
        upw = 0.5 * ( m_psi(i-1, j, k) + m_psi(i, j, k) );
        if (upw != 0.) upw = (upw > 0) ? 1. : -1.;

        // Interpolate lo
        val_lo = Evaluate(sp2,sp1,s,sm1,sm2,sm3,upw);
    }

    AMREX_GPU_DEVICE
    AMREX_FORCE_INLINE
    void
    InterpolateInY(const int& i,
                   const int& j,
                   const int& k,
                   const int& qty_index,
                   amrex::Real val_hi,
                   amrex::Real val_lo)
    {
        // Data to interpolate on
        amrex::Real sp3 = m_phi(i  , j+3, k  , qty_index);
        amrex::Real sp2 = m_phi(i  , j+2, k  , qty_index);
        amrex::Real sp1 = m_phi(i  , j+1, k  , qty_index);
        amrex::Real s   = m_phi(i  , j  , k  , qty_index);
        amrex::Real sm1 = m_phi(i  , j-1, k  , qty_index);
        amrex::Real sm2 = m_phi(i  , j-2, k  , qty_index);
        amrex::Real sm2 = m_phi(i  , j-3, k  , qty_index);

        // Upwinding flag
        upw = 0.5 * ( m_psi(i, j+1, k) + m_psi(i, j, k) );
        if (upw != 0.) upw = (upw > 0) ? 1. : -1.;

        // Interpolate hi
        val_hi = Evaluate(sp3,sp2,sp1,s,sm1,sm2,upw);

        // Upwinding flag
        upw = 0.5 * ( m_psi(i, j-1, k) + m_psi(i, j, k) );
        if (upw != 0.) upw = (upw > 0) ? 1. : -1.;

        // Interpolate lo
        val_lo = Evaluate(sp2,sp1,s,sm1,sm2,sm3,upw);
    }

    AMREX_GPU_DEVICE
    AMREX_FORCE_INLINE
    void
    InterpolateInZ(const int& i,
                   const int& j,
                   const int& k,
                   const int& qty_index,
                   amrex::Real val_hi,
                   amrex::Real val_lo)
    {
        // Data to interpolate on
        amrex::Real sp3 = m_phi(i  , j  , k+3, qty_index);
        amrex::Real sp2 = m_phi(i  , j  , k+2, qty_index);
        amrex::Real sp1 = m_phi(i  , j  , k+1, qty_index);
        amrex::Real s   = m_phi(i  , j  , k  , qty_index);
        amrex::Real sm1 = m_phi(i  , j  , k-1, qty_index);
        amrex::Real sm2 = m_phi(i  , j  , k-2, qty_index);
        amrex::Real sm3 = m_phi(i  , j  , k-3, qty_index);

        // Upwinding flag
        upw = 0.5 * ( m_psi(i, j, k+1) + m_psi(i, j, k) );
        if (upw != 0.) upw = (upw > 0) ? 1. : -1.;

        // Interpolate hi
        val_hi = Evaluate(sp3,sp2,sp1,s,sm1,sm2,upw);

        // Upwinding flag
        upw = 0.5 * ( m_psi(i, j, k-1) + m_psi(i, j, k) );
        if (upw != 0.) upw = (upw > 0) ? 1. : -1.;

        // Interpolate lo
        val_lo = Evaluate(sp2,sp1,s,sm1,sm2,sm3,upw);
    }

    AMREX_GPU_DEVICE
    AMREX_FORCE_INLINE
    amrex::Real
    Evaluate(const amrex::Real& sp2,
             const amrex::Real& sp1,
             const amrex::Real& s,
             const amrex::Real& sm1,
             const amrex::Real& sm2,
             const amrex::Real& sm3,
             const amrex::Real& upw)
    {
        // Averages and diffs
        a1 = (s + sm1);
        a2 = (sp1 + sm2);
        a3 = (sp2 + sm3);
        d1 = (s - sm1);
        d2 = (sp1 - sm2);
        d3 = (sp2 - sm3);

        // Interpolated value
        return ( g1*a1 - g2*a2 + g3*a3 - upw*g3 * (d3 - 5.0*d2 + 10.0*d1) );
    }

private:
    const amrex::Array4<const amrex::Real>& m_phi; // Quantity to interpolate
    const amrex::Array4<const amrex::Real>& m_psi; // Upwind indicator
    amrex::Real a1 = 0.; amrex::Real a2 = 0.; amrex::Real a3 = 0.;
    amrex::Real d1 = 0.; amrex::Real d2 = 0.; amrex::Real d3 = 0.;
    amrex::Real upw = 0.;
    constexpr amrex::Real g1=(37.0/60.0);
    constexpr amrex::Real g2=(2.0/15.0);
    constexpr amrex::Real g3=(1.0/60.0);
};

/**
 * Interpolation operators used for 5th order upwind scheme
 */
struct UPWIND6
{
    UPWIND6(const amrex::Array4<const amrex::Real>& phi,
            const amrex::Array4<const amrex::Real>& psi)
        : m_phi(phi), m_psi(psi) {}

    AMREX_GPU_DEVICE
    AMREX_FORCE_INLINE
    void
    InterpolateInX(const int& i,
                   const int& j,
                   const int& k,
                   const int& qty_index,
                   amrex::Real val_hi,
                   amrex::Real val_lo)
    {
        // Data to interpolate on
        amrex::Real sp3 = m_phi(i+3, j  , k  , qty_index);
        amrex::Real sp2 = m_phi(i+2, j  , k  , qty_index);
        amrex::Real sp1 = m_phi(i+1, j  , k  , qty_index);
        amrex::Real s   = m_phi(i  , j  , k  , qty_index);
        amrex::Real sm1 = m_phi(i-1, j  , k  , qty_index);
        amrex::Real sm2 = m_phi(i-2, j  , k  , qty_index);
        amrex::Real sm3 = m_phi(i-3, j  , k  , qty_index);

        // Interpolate hi
        val_hi = Evaluate(sp3,sp2,sp1,s,sm1,sm2,upw);

        // Interpolate lo
        val_lo = Evaluate(sp2,sp1,s,sm1,sm2,sm3,upw);
    }

    AMREX_GPU_DEVICE
    AMREX_FORCE_INLINE
    void
    InterpolateInY(const int& i,
                   const int& j,
                   const int& k,
                   const int& qty_index,
                   amrex::Real val_hi,
                   amrex::Real val_lo)
    {
        // Data to interpolate on
        amrex::Real sp3 = m_phi(i  , j+3, k  , qty_index);
        amrex::Real sp2 = m_phi(i  , j+2, k  , qty_index);
        amrex::Real sp1 = m_phi(i  , j+1, k  , qty_index);
        amrex::Real s   = m_phi(i  , j  , k  , qty_index);
        amrex::Real sm1 = m_phi(i  , j-1, k  , qty_index);
        amrex::Real sm2 = m_phi(i  , j-2, k  , qty_index);
        amrex::Real sm2 = m_phi(i  , j-3, k  , qty_index);

        // Interpolate hi
        val_hi = Evaluate(sp3,sp2,sp1,s,sm1,sm2,upw);

        // Interpolate lo
        val_lo = Evaluate(sp2,sp1,s,sm1,sm2,sm3,upw);
    }

    AMREX_GPU_DEVICE
    AMREX_FORCE_INLINE
    void
    InterpolateInZ(const int& i,
                   const int& j,
                   const int& k,
                   const int& qty_index,
                   amrex::Real val_hi,
                   amrex::Real val_lo)
    {
        // Data to interpolate on
        amrex::Real sp3 = m_phi(i  , j  , k+3, qty_index);
        amrex::Real sp2 = m_phi(i  , j  , k+2, qty_index);
        amrex::Real sp1 = m_phi(i  , j  , k+1, qty_index);
        amrex::Real s   = m_phi(i  , j  , k  , qty_index);
        amrex::Real sm1 = m_phi(i  , j  , k-1, qty_index);
        amrex::Real sm2 = m_phi(i  , j  , k-2, qty_index);
        amrex::Real sm3 = m_phi(i  , j  , k-3, qty_index);

        // Interpolate hi
        val_hi = Evaluate(sp3,sp2,sp1,s,sm1,sm2,upw);

        // Interpolate lo
        val_lo = Evaluate(sp2,sp1,s,sm1,sm2,sm3,upw);
    }

    AMREX_GPU_DEVICE
    AMREX_FORCE_INLINE
    amrex::Real
    Evaluate(const amrex::Real& sp2,
             const amrex::Real& sp1,
             const amrex::Real& s,
             const amrex::Real& sm1,
             const amrex::Real& sm2,
             const amrex::Real& sm3)
    {
        // Averages and diffs
        a1 = (s + sm1);
        a2 = (sp1 + sm2);
        a3 = (sp2 + sm3);
        d1 = (s - sm1);
        d2 = (sp1 - sm2);
        d3 = (sp2 - sm3);

        // Interpolated value
        return ( g1*a1 - g2*a2 + g3*a3 );
    }

private:
    const amrex::Array4<const amrex::Real>& m_phi; // Quantity to interpolate
    const amrex::Array4<const amrex::Real>& m_psi; // Upwind indicator
    amrex::Real a1 = 0.; amrex::Real a2 = 0.; amrex::Real a3 = 0.;
    amrex::Real d1 = 0.; amrex::Real d2 = 0.; amrex::Real d3 = 0.;
    constexpr amrex::Real g1=(37.0/60.0);
    constexpr amrex::Real g2=(2.0/15.0);
    constexpr amrex::Real g3=(1.0/60.0);
};
#endif
#endif
