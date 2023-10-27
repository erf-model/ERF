#ifndef INTERPOLATE_UPW_H_
#define INTERPOLATE_UPW_H_

#include "DataStruct.H"

/**
 * Interpolation operators used for 2nd order centered scheme
 */
struct CENTERED2
{
    CENTERED2 (const amrex::Array4<const amrex::Real>& phi)
        : m_phi(phi) {}

    AMREX_GPU_DEVICE
    AMREX_FORCE_INLINE
    void
    InterpolateInX (const int& i,
                    const int& j,
                    const int& k,
                    const int& qty_index,
                    amrex::Real& val_hi,
                    amrex::Real& val_lo,
                    amrex::Real /*upw_hi*/,
                    amrex::Real /*upw_lo*/) const
    {
        // Data to interpolate on
        amrex::Real sp1 = m_phi(i+1, j  , k  , qty_index);
        amrex::Real s   = m_phi(i  , j  , k  , qty_index);
        amrex::Real sm1 = m_phi(i-1, j  , k  , qty_index);

        // Interpolate hi
        val_hi = Evaluate(sp1,s);

        // Interpolate lo
        val_lo = Evaluate(s,sm1);
    }

    AMREX_GPU_DEVICE
    AMREX_FORCE_INLINE
    void
    InterpolateInY (const int& i,
                    const int& j,
                    const int& k,
                    const int& qty_index,
                    amrex::Real& val_hi,
                    amrex::Real& val_lo,
                    amrex::Real /*upw_hi*/,
                    amrex::Real /*upw_lo*/) const
    {
        // Data to interpolate on
        amrex::Real sp1 = m_phi(i  , j+1, k  , qty_index);
        amrex::Real s   = m_phi(i  , j  , k  , qty_index);
        amrex::Real sm1 = m_phi(i  , j-1, k  , qty_index);

        // Interpolate hi
        val_hi = Evaluate(sp1,s);

        // Interpolate lo
        val_lo = Evaluate(s,sm1);
    }

    AMREX_GPU_DEVICE
    AMREX_FORCE_INLINE
    void
    InterpolateInZ_lo (const int& i,
                       const int& j,
                       const int& k,
                       const int& qty_index,
                       amrex::Real& val_lo,
                       amrex::Real /*upw_lo*/) const
    {
        // Data to interpolate on
        amrex::Real s   = m_phi(i  , j  , k  , qty_index);
        amrex::Real sm1 = m_phi(i  , j  , k-1, qty_index);

        // Interpolate lo
        val_lo = Evaluate(s,sm1);
    }

    AMREX_GPU_DEVICE
    AMREX_FORCE_INLINE
    void
    InterpolateInZ_hi (const int& i,
                       const int& j,
                       const int& k,
                       const int& qty_index,
                       amrex::Real& val_hi,
                       amrex::Real /*upw_hi*/) const
    {
        // Data to interpolate on
        amrex::Real sp1 = m_phi(i  , j  , k+1, qty_index);
        amrex::Real s   = m_phi(i  , j  , k  , qty_index);

        // Interpolate hi
        val_hi = Evaluate(sp1,s);
    }

    AMREX_GPU_DEVICE
    AMREX_FORCE_INLINE
    amrex::Real
    Evaluate(const amrex::Real& s,
             const amrex::Real& sm1) const
    {
        // Averages and diffs
        amrex::Real a1 = (s   + sm1);

        // Interpolated value
        return ( g1*a1 );
    }

private:
    amrex::Array4<const amrex::Real> m_phi;   // Quantity to interpolate
    static constexpr amrex::Real g1=(0.5);
};

/**
 * Interpolation operators used for 3rd order upwind scheme
 */
struct UPWIND3
{
    UPWIND3 (const amrex::Array4<const amrex::Real>& phi)
        : m_phi(phi) {}

    AMREX_GPU_DEVICE
    AMREX_FORCE_INLINE
    void
    InterpolateInX (const int& i,
                    const int& j,
                    const int& k,
                    const int& qty_index,
                    amrex::Real& val_hi,
                    amrex::Real& val_lo,
                    amrex::Real upw_hi,
                    amrex::Real upw_lo) const
    {
        // Data to interpolate on
        amrex::Real sp2 = m_phi(i+2, j  , k  , qty_index);
        amrex::Real sp1 = m_phi(i+1, j  , k  , qty_index);
        amrex::Real s   = m_phi(i  , j  , k  , qty_index);
        amrex::Real sm1 = m_phi(i-1, j  , k  , qty_index);
        amrex::Real sm2 = m_phi(i-2, j  , k  , qty_index);

        // Upwinding flags
        if (upw_hi != 0.) upw_hi = (upw_hi > 0) ? 1. : -1.;
        if (upw_lo != 0.) upw_lo = (upw_lo > 0) ? 1. : -1.;

        // Interpolate hi
        val_hi = Evaluate(sp2,sp1,s,sm1,upw_hi);

        // Interpolate lo
        val_lo = Evaluate(sp1,s,sm1,sm2,upw_lo);
    }

    AMREX_GPU_DEVICE
    AMREX_FORCE_INLINE
    void
    InterpolateInY (const int& i,
                    const int& j,
                    const int& k,
                    const int& qty_index,
                    amrex::Real& val_hi,
                    amrex::Real& val_lo,
                    amrex::Real upw_hi,
                    amrex::Real upw_lo) const
    {
        // Data to interpolate on
        amrex::Real sp2 = m_phi(i  , j+2, k  , qty_index);
        amrex::Real sp1 = m_phi(i  , j+1, k  , qty_index);
        amrex::Real s   = m_phi(i  , j  , k  , qty_index);
        amrex::Real sm1 = m_phi(i  , j-1, k  , qty_index);
        amrex::Real sm2 = m_phi(i  , j-2, k  , qty_index);

        // Upwinding flags
        if (upw_hi != 0.) upw_hi = (upw_hi > 0) ? 1. : -1.;
        if (upw_lo != 0.) upw_lo = (upw_lo > 0) ? 1. : -1.;

        // Interpolate hi
        val_hi = Evaluate(sp2,sp1,s,sm1,upw_hi);

        // Interpolate lo
        val_lo = Evaluate(sp1,s,sm1,sm2,upw_lo);
    }

    AMREX_GPU_DEVICE
    AMREX_FORCE_INLINE
    void
    InterpolateInZ_lo (const int& i,
                       const int& j,
                       const int& k,
                       const int& qty_index,
                       amrex::Real& val_lo,
                       amrex::Real upw_lo) const
    {
        // Data to interpolate on
        amrex::Real sp1 = m_phi(i  , j  , k+1, qty_index);
        amrex::Real s   = m_phi(i  , j  , k  , qty_index);
        amrex::Real sm1 = m_phi(i  , j  , k-1, qty_index);
        amrex::Real sm2 = m_phi(i  , j  , k-2, qty_index);

        // Upwinding flags
        if (upw_lo != 0.) upw_lo = (upw_lo > 0) ? 1. : -1.;

        // Interpolate lo
        val_lo = Evaluate(sp1,s,sm1,sm2,upw_lo);
    }

    AMREX_GPU_DEVICE
    AMREX_FORCE_INLINE
    void
    InterpolateInZ_hi (const int& i,
                       const int& j,
                       const int& k,
                       const int& qty_index,
                       amrex::Real& val_hi,
                       amrex::Real upw_hi) const
    {
        // Data to interpolate on
        amrex::Real sp2 = m_phi(i  , j  , k+2, qty_index);
        amrex::Real sp1 = m_phi(i  , j  , k+1, qty_index);
        amrex::Real s   = m_phi(i  , j  , k  , qty_index);
        amrex::Real sm1 = m_phi(i  , j  , k-1, qty_index);

        // Upwinding flags
        if (upw_hi != 0.) upw_hi = (upw_hi > 0) ? 1. : -1.;

        // Interpolate lo
        val_hi = Evaluate(sp2,sp1,s,sm1,upw_hi);
    }

    AMREX_GPU_DEVICE
    AMREX_FORCE_INLINE
    amrex::Real
    Evaluate (const amrex::Real& sp1,
              const amrex::Real& s,
              const amrex::Real& sm1,
              const amrex::Real& sm2,
              const amrex::Real& upw) const
    {
        // Averages and diffs
        amrex::Real a1 = (s   + sm1);
        amrex::Real d1 = (s   - sm1);
        amrex::Real a2 = (sp1 + sm2);
        amrex::Real d2 = (sp1 - sm2);

        // Interpolated value
        return ( g1*a1 - g2*a2 + upw * g2 * (d2 - 3.0*d1) );
    }

private:
    amrex::Array4<const amrex::Real> m_phi;   // Quantity to interpolate
    static constexpr amrex::Real g1=(7.0/12.0);
    static constexpr amrex::Real g2=(1.0/12.0);
};


/**
 * Interpolation operators used for 4th order centered scheme
 */
struct CENTERED4
{
    CENTERED4 (const amrex::Array4<const amrex::Real>& phi)
        : m_phi(phi) {}

    AMREX_GPU_DEVICE
    AMREX_FORCE_INLINE
    void
    InterpolateInX (const int& i,
                    const int& j,
                    const int& k,
                    const int& qty_index,
                    amrex::Real& val_hi,
                    amrex::Real& val_lo,
                    amrex::Real /*upw_hi*/,
                    amrex::Real /*upw_lo*/) const
    {
        // Data to interpolate on
        amrex::Real sp2 = m_phi(i+2, j  , k  , qty_index);
        amrex::Real sp1 = m_phi(i+1, j  , k  , qty_index);
        amrex::Real s   = m_phi(i  , j  , k  , qty_index);
        amrex::Real sm1 = m_phi(i-1, j  , k  , qty_index);
        amrex::Real sm2 = m_phi(i-2, j  , k  , qty_index);

        // Interpolate hi
        val_hi = Evaluate(sp2,sp1,s,sm1);

        // Interpolate lo
        val_lo = Evaluate(sp1,s,sm1,sm2);
    }

    AMREX_GPU_DEVICE
    AMREX_FORCE_INLINE
    void
    InterpolateInY (const int& i,
                    const int& j,
                    const int& k,
                    const int& qty_index,
                    amrex::Real& val_hi,
                    amrex::Real& val_lo,
                    amrex::Real /*upw_hi*/,
                    amrex::Real /*upw_lo*/) const
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
    InterpolateInZ_lo (const int& i,
                       const int& j,
                       const int& k,
                       const int& qty_index,
                       amrex::Real& val_lo,
                       amrex::Real /*upw_lo*/) const
    {
        // Data to interpolate on
        amrex::Real sp1 = m_phi(i  , j  , k+1, qty_index);
        amrex::Real s   = m_phi(i  , j  , k  , qty_index);
        amrex::Real sm1 = m_phi(i  , j  , k-1, qty_index);
        amrex::Real sm2 = m_phi(i  , j  , k-2, qty_index);

        // Interpolate lo
        val_lo = Evaluate(sp1,s,sm1,sm2);
    }

    AMREX_GPU_DEVICE
    AMREX_FORCE_INLINE
    void
    InterpolateInZ_hi (const int& i,
                       const int& j,
                       const int& k,
                       const int& qty_index,
                       amrex::Real& val_hi,
                       amrex::Real /*upw_hi*/) const
    {
        // Data to interpolate on
        amrex::Real sp2 = m_phi(i  , j  , k+2, qty_index);
        amrex::Real sp1 = m_phi(i  , j  , k+1, qty_index);
        amrex::Real s   = m_phi(i  , j  , k  , qty_index);
        amrex::Real sm1 = m_phi(i  , j  , k-1, qty_index);

        // Interpolate hi
        val_hi = Evaluate(sp2,sp1,s,sm1);
    }

    AMREX_GPU_DEVICE
    AMREX_FORCE_INLINE
    amrex::Real
    Evaluate (const amrex::Real& sp1,
              const amrex::Real& s,
              const amrex::Real& sm1,
              const amrex::Real& sm2) const
    {
        // Averages and diffs
        amrex::Real a1 = (s   + sm1);
        amrex::Real a2 = (sp1 + sm2);

        // Interpolated value
        return ( g1*a1 - g2*a2 );
    }

private:
    amrex::Array4<const amrex::Real> m_phi;   // Quantity to interpolate
    static constexpr amrex::Real g1=(7.0/12.0);
    static constexpr amrex::Real g2=(1.0/12.0);
};

/**
 * Interpolation operators used for 5th order upwind scheme
 */
struct UPWIND5
{
    UPWIND5 (const amrex::Array4<const amrex::Real>& phi)
        : m_phi(phi) {}

    AMREX_GPU_DEVICE
    AMREX_FORCE_INLINE
    void
    InterpolateInX (const int& i,
                    const int& j,
                    const int& k,
                    const int& qty_index,
                    amrex::Real& val_hi,
                    amrex::Real& val_lo,
                    amrex::Real upw_hi,
                    amrex::Real upw_lo) const
    {
        // Data to interpolate on
        amrex::Real sp3 = m_phi(i+3, j  , k  , qty_index);
        amrex::Real sp2 = m_phi(i+2, j  , k  , qty_index);
        amrex::Real sp1 = m_phi(i+1, j  , k  , qty_index);
        amrex::Real s   = m_phi(i  , j  , k  , qty_index);
        amrex::Real sm1 = m_phi(i-1, j  , k  , qty_index);
        amrex::Real sm2 = m_phi(i-2, j  , k  , qty_index);
        amrex::Real sm3 = m_phi(i-3, j  , k  , qty_index);

        // Upwinding flags
        if (upw_hi != 0.) upw_hi = (upw_hi > 0) ? 1. : -1.;
        if (upw_lo != 0.) upw_lo = (upw_lo > 0) ? 1. : -1.;

        // Interpolate hi
        val_hi = Evaluate(sp3,sp2,sp1,s,sm1,sm2,upw_hi);

        // Interpolate lo
        val_lo = Evaluate(sp2,sp1,s,sm1,sm2,sm3,upw_lo);
    }

    AMREX_GPU_DEVICE
    AMREX_FORCE_INLINE
    void
    InterpolateInY (const int& i,
                    const int& j,
                    const int& k,
                    const int& qty_index,
                    amrex::Real& val_hi,
                    amrex::Real& val_lo,
                    amrex::Real upw_hi,
                    amrex::Real upw_lo) const
    {
        // Data to interpolate on
        amrex::Real sp3 = m_phi(i  , j+3, k  , qty_index);
        amrex::Real sp2 = m_phi(i  , j+2, k  , qty_index);
        amrex::Real sp1 = m_phi(i  , j+1, k  , qty_index);
        amrex::Real s   = m_phi(i  , j  , k  , qty_index);
        amrex::Real sm1 = m_phi(i  , j-1, k  , qty_index);
        amrex::Real sm2 = m_phi(i  , j-2, k  , qty_index);
        amrex::Real sm3 = m_phi(i  , j-3, k  , qty_index);

        // Upwinding flags
        if (upw_hi != 0.) upw_hi = (upw_hi > 0) ? 1. : -1.;
        if (upw_lo != 0.) upw_lo = (upw_lo > 0) ? 1. : -1.;

        // Interpolate hi
        val_hi = Evaluate(sp3,sp2,sp1,s,sm1,sm2,upw_hi);

        // Interpolate lo
        val_lo = Evaluate(sp2,sp1,s,sm1,sm2,sm3,upw_lo);
    }

    AMREX_GPU_DEVICE
    AMREX_FORCE_INLINE
    void
    InterpolateInZ_lo (const int& i,
                       const int& j,
                       const int& k,
                       const int& qty_index,
                       amrex::Real& val_lo,
                       amrex::Real upw_lo) const
    {
        // Data to interpolate on
        amrex::Real sp2 = m_phi(i  , j  , k+2, qty_index);
        amrex::Real sp1 = m_phi(i  , j  , k+1, qty_index);
        amrex::Real s   = m_phi(i  , j  , k  , qty_index);
        amrex::Real sm1 = m_phi(i  , j  , k-1, qty_index);
        amrex::Real sm2 = m_phi(i  , j  , k-2, qty_index);
        amrex::Real sm3 = m_phi(i  , j  , k-3, qty_index);

        // Upwinding flags
        if (upw_lo != 0.) upw_lo = (upw_lo > 0) ? 1. : -1.;

        // Interpolate lo
        val_lo = Evaluate(sp2,sp1,s,sm1,sm2,sm3,upw_lo);
    }

    AMREX_GPU_DEVICE
    AMREX_FORCE_INLINE
    void
    InterpolateInZ_hi (const int& i,
                       const int& j,
                       const int& k,
                       const int& qty_index,
                       amrex::Real& val_hi,
                       amrex::Real upw_hi) const
    {
        // Data to interpolate on
        amrex::Real sp3 = m_phi(i  , j  , k+3, qty_index);
        amrex::Real sp2 = m_phi(i  , j  , k+2, qty_index);
        amrex::Real sp1 = m_phi(i  , j  , k+1, qty_index);
        amrex::Real s   = m_phi(i  , j  , k  , qty_index);
        amrex::Real sm1 = m_phi(i  , j  , k-1, qty_index);
        amrex::Real sm2 = m_phi(i  , j  , k-2, qty_index);

        // Upwinding flags
        if (upw_hi != 0.) upw_hi = (upw_hi > 0) ? 1. : -1.;

        // Interpolate hi
        val_hi = Evaluate(sp3,sp2,sp1,s,sm1,sm2,upw_hi);
    }

    AMREX_GPU_DEVICE
    AMREX_FORCE_INLINE
    amrex::Real
    Evaluate (const amrex::Real& sp2,
              const amrex::Real& sp1,
              const amrex::Real& s,
              const amrex::Real& sm1,
              const amrex::Real& sm2,
              const amrex::Real& sm3,
              const amrex::Real& upw) const
    {
        // Averages and diffs
        amrex::Real a1 = (s   + sm1);
        amrex::Real a2 = (sp1 + sm2);
        amrex::Real a3 = (sp2 + sm3);
        amrex::Real d1 = (s   - sm1);
        amrex::Real d2 = (sp1 - sm2);
        amrex::Real d3 = (sp2 - sm3);

        // Interpolated value
        return ( g1*a1 - g2*a2 + g3*a3 - upw * g3 * (d3 - 5.0*d2 + 10.0*d1) );
    }

private:
    amrex::Array4<const amrex::Real> m_phi;   // Quantity to interpolate
    static constexpr amrex::Real g1=(37.0/60.0);
    static constexpr amrex::Real g2=(2.0/15.0);
    static constexpr amrex::Real g3=(1.0/60.0);
};

/**
 * Interpolation operators used for 6th order centered scheme
 */
struct CENTERED6
{
    CENTERED6 (const amrex::Array4<const amrex::Real>& phi)
        : m_phi(phi) {}

    AMREX_GPU_DEVICE
    AMREX_FORCE_INLINE
    void
    InterpolateInX (const int& i,
                    const int& j,
                    const int& k,
                    const int& qty_index,
                    amrex::Real& val_hi,
                    amrex::Real& val_lo,
                    amrex::Real /*upw_hi*/,
                    amrex::Real /*upw_lo*/) const
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
        val_hi = Evaluate(sp3,sp2,sp1,s,sm1,sm2);

        // Interpolate lo
        val_lo = Evaluate(sp2,sp1,s,sm1,sm2,sm3);
    }

    AMREX_GPU_DEVICE
    AMREX_FORCE_INLINE
    void
    InterpolateInY (const int& i,
                    const int& j,
                    const int& k,
                    const int& qty_index,
                    amrex::Real& val_hi,
                    amrex::Real& val_lo,
                    amrex::Real /*upw_hi*/,
                    amrex::Real /*upw_lo*/) const
    {
        // Data to interpolate on
        amrex::Real sp3 = m_phi(i  , j+3, k  , qty_index);
        amrex::Real sp2 = m_phi(i  , j+2, k  , qty_index);
        amrex::Real sp1 = m_phi(i  , j+1, k  , qty_index);
        amrex::Real s   = m_phi(i  , j  , k  , qty_index);
        amrex::Real sm1 = m_phi(i  , j-1, k  , qty_index);
        amrex::Real sm2 = m_phi(i  , j-2, k  , qty_index);
        amrex::Real sm3 = m_phi(i  , j-3, k  , qty_index);

        // Interpolate hi
        val_hi = Evaluate(sp3,sp2,sp1,s,sm1,sm2);

        // Interpolate lo
        val_lo = Evaluate(sp2,sp1,s,sm1,sm2,sm3);
    }

    AMREX_GPU_DEVICE
    AMREX_FORCE_INLINE
    void
    InterpolateInZ_lo (const int& i,
                       const int& j,
                       const int& k,
                       const int& qty_index,
                       amrex::Real& val_lo,
                       amrex::Real /*upw_lo*/) const
    {
        // Data to interpolate on
        amrex::Real sp2 = m_phi(i  , j  , k+2, qty_index);
        amrex::Real sp1 = m_phi(i  , j  , k+1, qty_index);
        amrex::Real s   = m_phi(i  , j  , k  , qty_index);
        amrex::Real sm1 = m_phi(i  , j  , k-1, qty_index);
        amrex::Real sm2 = m_phi(i  , j  , k-2, qty_index);
        amrex::Real sm3 = m_phi(i  , j  , k-3, qty_index);

        // Interpolate lo
        val_lo = Evaluate(sp2,sp1,s,sm1,sm2,sm3);
    }

    AMREX_GPU_DEVICE
    AMREX_FORCE_INLINE
    void
    InterpolateInZ_hi (const int& i,
                       const int& j,
                       const int& k,
                       const int& qty_index,
                       amrex::Real& val_hi,
                       amrex::Real /*upw_hi*/) const
    {
        // Data to interpolate on
        amrex::Real sp3 = m_phi(i  , j  , k+3, qty_index);
        amrex::Real sp2 = m_phi(i  , j  , k+2, qty_index);
        amrex::Real sp1 = m_phi(i  , j  , k+1, qty_index);
        amrex::Real s   = m_phi(i  , j  , k  , qty_index);
        amrex::Real sm1 = m_phi(i  , j  , k-1, qty_index);
        amrex::Real sm2 = m_phi(i  , j  , k-2, qty_index);

        // Interpolate hi
        val_hi = Evaluate(sp3,sp2,sp1,s,sm1,sm2);
    }

    AMREX_GPU_DEVICE
    AMREX_FORCE_INLINE
    amrex::Real
    Evaluate (const amrex::Real& sp2,
              const amrex::Real& sp1,
              const amrex::Real& s,
              const amrex::Real& sm1,
              const amrex::Real& sm2,
              const amrex::Real& sm3) const
    {
        // Averages and diffs
        amrex::Real a1 = (s   + sm1);
        amrex::Real a2 = (sp1 + sm2);
        amrex::Real a3 = (sp2 + sm3);

        // Interpolated value
        return ( g1*a1 - g2*a2 + g3*a3 );
    }

private:
    amrex::Array4<const amrex::Real> m_phi;   // Quantity to interpolate
    static constexpr amrex::Real g1=(37.0/60.0);
    static constexpr amrex::Real g2=(2.0/15.0);
    static constexpr amrex::Real g3=(1.0/60.0);
};

/**
 * Interpolation operators used for all central/upwind schemes
 */
struct UPWINDALL
{
    UPWINDALL (const amrex::Array4<const amrex::Real>& phi)
        : m_phi(phi) {}

    AMREX_GPU_DEVICE
    AMREX_FORCE_INLINE
    void
    InterpolateInZ_lo (const int& i,
                       const int& j,
                       const int& k,
                       const int& qty_index,
                       amrex::Real& val_lo,
                       amrex::Real upw_lo,
                       const AdvType adv_type) const
    {
        if (adv_type == AdvType::Centered_2nd) {
            val_lo = 0.5 * ( m_phi(i,j,k,qty_index) + m_phi(i,j,k-1,qty_index) );
            return;
        } else {
            // Data to interpolate on
            amrex::Real sp2 = (adv_type == AdvType::Upwind_5th || adv_type == AdvType::Centered_6th ) ? m_phi(i  , j  , k+2, qty_index): 0.;
            amrex::Real sp1 = m_phi(i  , j  , k+1, qty_index);
            amrex::Real s   = m_phi(i  , j  , k  , qty_index);
            amrex::Real sm1 = m_phi(i  , j  , k-1, qty_index);
            amrex::Real sm2 = m_phi(i  , j  , k-2, qty_index);
            amrex::Real sm3 = (adv_type == AdvType::Upwind_5th || adv_type == AdvType::Centered_6th ) ? m_phi(i  , j  , k-3, qty_index) : 0.;

            // Upwinding flags
            if (upw_lo != 0.) upw_lo = (upw_lo > 0) ? 1. : -1.;

            // Interpolate lo
            val_lo = Evaluate(sp2,sp1,s,sm1,sm2,sm3,upw_lo,adv_type);
        }
    }

    AMREX_GPU_DEVICE
    AMREX_FORCE_INLINE
    void
    InterpolateInZ_hi (const int& i,
                       const int& j,
                       const int& k,
                       const int& qty_index,
                       amrex::Real& val_hi,
                       amrex::Real upw_hi,
                       const AdvType adv_type) const
    {
        if (adv_type == AdvType::Centered_2nd) {
            val_hi = 0.5 * ( m_phi(i,j,k,qty_index) + m_phi(i,j,k+1,qty_index) );
            return;
        } else {
            // Data to interpolate on
            amrex::Real sp3 = (adv_type == AdvType::Upwind_5th || adv_type == AdvType::Centered_6th ) ? m_phi(i  , j  , k+3, qty_index): 0.;
            amrex::Real sp2 = m_phi(i  , j  , k+2, qty_index);
            amrex::Real sp1 = m_phi(i  , j  , k+1, qty_index);
            amrex::Real s   = m_phi(i  , j  , k  , qty_index);
            amrex::Real sm1 = m_phi(i  , j  , k-1, qty_index);
            amrex::Real sm2 = (adv_type == AdvType::Upwind_5th || adv_type == AdvType::Centered_6th ) ? m_phi(i  , j  , k-2, qty_index) : 0.;

            // Upwinding flags
            if (upw_hi != 0.) upw_hi = (upw_hi > 0) ? 1. : -1.;

            // Interpolate hi
            val_hi = Evaluate(sp3,sp2,sp1,s,sm1,sm2,upw_hi,adv_type);
        }
    }

    AMREX_GPU_DEVICE
    AMREX_FORCE_INLINE
    amrex::Real
    Evaluate (const amrex::Real& sp2,
              const amrex::Real& sp1,
              const amrex::Real& s,
              const amrex::Real& sm1,
              const amrex::Real& sm2,
              const amrex::Real& sm3,
              const amrex::Real& upw,
              const AdvType adv_type) const
    {
        // Averages and diffs
        amrex::Real a1 = (s   + sm1);
        amrex::Real a2 = (sp1 + sm2);
        amrex::Real d1 = (s   - sm1);
        amrex::Real d2 = (sp1 - sm2);
        amrex::Real a3 = (sp2 + sm3);
        amrex::Real d3 = (sp2 - sm3);

        // Interpolated value
        amrex::Real iv(0.0);
        if (adv_type == AdvType::Centered_2nd) {
            iv = 0.5 * a1;
        } else if (adv_type == AdvType::Upwind_3rd) {
            iv = g1_3_4*a1 - g2_3_4*a2 + upw * g2_3_4 * (d2 - 3.0*d1);
        } else if (adv_type == AdvType::Centered_4th) {
            iv = g1_3_4*a1 - g2_3_4*a2;
        } else if (adv_type == AdvType::Upwind_5th) {
            iv = g1_5_6*a1 - g2_5_6*a2 + g3_5_6*a3 - upw * g3_5_6 * (d3 - 5.0*d2 + 10.0*d1);
        } else if (adv_type == AdvType::Centered_6th) {
            iv = g1_5_6*a1 - g2_5_6*a2 + g3_5_6*a3;
        }
        return ( iv );
    }

private:
    amrex::Array4<const amrex::Real> m_phi;   // Quantity to interpolate
    static constexpr amrex::Real g1_3_4=( 7.0/12.0);
    static constexpr amrex::Real g2_3_4=( 1.0/12.0);
    static constexpr amrex::Real g1_5_6=(37.0/60.0);
    static constexpr amrex::Real g2_5_6=( 2.0/15.0);
    static constexpr amrex::Real g3_5_6=( 1.0/60.0);
};
#endif
