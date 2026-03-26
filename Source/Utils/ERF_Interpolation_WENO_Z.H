#ifndef ERF_INTERPOLATE_WENO_Z_H_
#define ERF_INTERPOLATE_WENO_Z_H_

#include "ERF_DataStruct.H"

/**
 * Interpolation operators used for WENO_Z-3 scheme
 */
struct WENO_Z3
{
    WENO_Z3 (const amrex::Array4<const amrex::Real>& phi,
             const amrex::Real /*upw_frac*/)
        : m_phi(phi) {}

    AMREX_GPU_DEVICE
    AMREX_FORCE_INLINE
    void
    InterpolateInX (const int& i,
                    const int& j,
                    const int& k,
                    const int& qty_index,
                    amrex::Real& val_lo,
                    amrex::Real upw_lo) const
    {
        // Upwinding flags
        if (upw_lo != zero) upw_lo = (upw_lo > 0) ? one : -one;

        // Data to interpolate on
        amrex::Real sp1 = m_phi(i+1, j  , k  , qty_index);
        amrex::Real s   = m_phi(i  , j  , k  , qty_index);
        amrex::Real sm1 = m_phi(i-1, j  , k  , qty_index);
        amrex::Real sm2 = m_phi(i-2, j  , k  , qty_index);

        // Left and right fluxes
        amrex::Real Fhi = Evaluate(sm2,sm1,s  );
        amrex::Real Flo = Evaluate(sp1,s  ,sm1);

        // Numerical flux
        val_lo = (one + upw_lo)/two * Fhi + (one - upw_lo)/two * Flo;
    }

    AMREX_GPU_DEVICE
    AMREX_FORCE_INLINE
    void
    InterpolateInY (const int& i,
                    const int& j,
                    const int& k,
                    const int& qty_index,
                    amrex::Real& val_lo,
                    amrex::Real upw_lo) const
    {
        // Upwinding flags
        if (upw_lo != zero) upw_lo = (upw_lo > 0) ? one : -one;

        // Data to interpolate on
        amrex::Real sp1 = m_phi(i  , j+1, k  , qty_index);
        amrex::Real s   = m_phi(i  , j  , k  , qty_index);
        amrex::Real sm1 = m_phi(i  , j-1, k  , qty_index);
        amrex::Real sm2 = m_phi(i  , j-2, k  , qty_index);

        // Left and right fluxes
        amrex::Real Fhi = Evaluate(sm2,sm1,s  );
        amrex::Real Flo = Evaluate(sp1,s  ,sm1);

        // Numerical flux
        val_lo = (one + upw_lo)/two * Fhi + (one - upw_lo)/two * Flo;
    }

    AMREX_GPU_DEVICE
    AMREX_FORCE_INLINE
    void
    InterpolateInZ (const int& i,
                    const int& j,
                    const int& k,
                    const int& qty_index,
                    amrex::Real& val_lo,
                    amrex::Real upw_lo) const
    {
        // Upwinding flags
        if (upw_lo != zero) upw_lo = (upw_lo > 0) ? one : -one;

        // Data to interpolate on
        amrex::Real sp1 = m_phi(i  , j  , k+1, qty_index);
        amrex::Real s   = m_phi(i  , j  , k  , qty_index);
        amrex::Real sm1 = m_phi(i  , j  , k-1, qty_index);
        amrex::Real sm2 = m_phi(i  , j  , k-2, qty_index);

        // Left and right fluxes
        amrex::Real Fhi = Evaluate(sm2,sm1,s  );
        amrex::Real Flo = Evaluate(sp1,s  ,sm1);

        // Numerical flux
        val_lo = (one + upw_lo)/two * Fhi + (one - upw_lo)/two * Flo;
    }

    AMREX_GPU_DEVICE
    AMREX_FORCE_INLINE
    amrex::Real
    Evaluate (const amrex::Real& sm1,
              const amrex::Real& s  ,
              const amrex::Real& sp1) const
    {
        // Smoothing factors
        amrex::Real b1 = (s - sm1) * (s - sm1);
        amrex::Real b2 = (sp1 - s) * (sp1 - s);

        // Weight factors
        amrex::Real t5 = std::abs(b2 - b1);
        amrex::Real w1 = g1 * ( one + (t5*t5) / ((eps + b1) * (eps + b1)) );
        amrex::Real w2 = g2 * ( one + (t5*t5) / ((eps + b2) * (eps + b2)) );

        // Weight factor norm
        amrex::Real wsum = w1 + w2;

        // Taylor expansions
        amrex::Real v1 = -sm1 + three * s;
        amrex::Real v2 =  s   + sp1;

        // Interpolated value
        return ( (w1 * v1 + w2 * v2) / (two * wsum) );
    }

private:
    amrex::Array4<const amrex::Real> m_phi;   // Quantity to interpolate
    const amrex::Real eps=amrex::Real(1.0e-40);
    static constexpr amrex::Real g1=(one/three);
    static constexpr amrex::Real g2=(two/three);
};

/**
 * Interpolation operators used for WENO_MZQ3 scheme
 */
struct WENO_MZQ3
{
    WENO_MZQ3 (const amrex::Array4<const amrex::Real>& phi,
               const amrex::Real /*upw_frac*/)
        : m_phi(phi) {}

    AMREX_GPU_DEVICE
    AMREX_FORCE_INLINE
    void
    InterpolateInX (const int& i,
                    const int& j,
                    const int& k,
                    const int& qty_index,
                    amrex::Real& val_lo,
                    amrex::Real upw_lo) const
    {
        // Upwinding flags
        if (upw_lo != zero) upw_lo = (upw_lo > 0) ? one : -one;

        // Data to interpolate on
        amrex::Real sp1 = m_phi(i+1, j  , k  , qty_index);
        amrex::Real s   = m_phi(i  , j  , k  , qty_index);
        amrex::Real sm1 = m_phi(i-1, j  , k  , qty_index);
        amrex::Real sm2 = m_phi(i-2, j  , k  , qty_index);

        // Left and right fluxes
        amrex::Real Fhi = Evaluate(sm2,sm1,s  );
        amrex::Real Flo = Evaluate(sp1,s  ,sm1);

        // Numerical flux
        val_lo = (one + upw_lo)/two * Fhi + (one - upw_lo)/two * Flo;
    }

    AMREX_GPU_DEVICE
    AMREX_FORCE_INLINE
    void
    InterpolateInY (const int& i,
                    const int& j,
                    const int& k,
                    const int& qty_index,
                    amrex::Real& val_lo,
                    amrex::Real upw_lo) const
    {
        // Upwinding flags
        if (upw_lo != zero) upw_lo = (upw_lo > 0) ? one : -one;

        // Data to interpolate on
        amrex::Real sp1 = m_phi(i  , j+1, k  , qty_index);
        amrex::Real s   = m_phi(i  , j  , k  , qty_index);
        amrex::Real sm1 = m_phi(i  , j-1, k  , qty_index);
        amrex::Real sm2 = m_phi(i  , j-2, k  , qty_index);

        // Left and right fluxes
        amrex::Real Fhi = Evaluate(sm2,sm1,s  );
        amrex::Real Flo = Evaluate(sp1,s  ,sm1);

        // Numerical flux
        val_lo = (one + upw_lo)/two * Fhi + (one - upw_lo)/two * Flo;
    }

    AMREX_GPU_DEVICE
    AMREX_FORCE_INLINE
    void
    InterpolateInZ (const int& i,
                    const int& j,
                    const int& k,
                    const int& qty_index,
                    amrex::Real& val_lo,
                    amrex::Real upw_lo) const
    {
        // Upwinding flags
        if (upw_lo != zero) upw_lo = (upw_lo > 0) ? one : -one;

        // Data to interpolate on
        amrex::Real sp1 = m_phi(i  , j  , k+1, qty_index);
        amrex::Real s   = m_phi(i  , j  , k  , qty_index);
        amrex::Real sm1 = m_phi(i  , j  , k-1, qty_index);
        amrex::Real sm2 = m_phi(i  , j  , k-2, qty_index);

        // Left and right fluxes
        amrex::Real Fhi = Evaluate(sm2,sm1,s  );
        amrex::Real Flo = Evaluate(sp1,s  ,sm1);

        // Numerical flux
        val_lo = (one + upw_lo)/two * Fhi + (one - upw_lo)/two * Flo;
    }

    AMREX_GPU_DEVICE
    AMREX_FORCE_INLINE
    amrex::Real
    Evaluate (const amrex::Real& sm1,
              const amrex::Real& s  ,
              const amrex::Real& sp1) const
    {
        // Smoothing factors
        amrex::Real b1 = (s - sm1) * (s - sm1);
        amrex::Real b2 = (sp1 - s) * (sp1 - s);
        amrex::Real b3 = ( (amrex::Real(13.0) / amrex::Real(12.0)) * ((sm1 - two*s + sp1)*(sm1 - two*s + sp1)) ) + ( ((sp1 - sm1)*(sp1 - sm1)) / amrex::Real(4.0) );

        // Weight factors
        amrex::Real t5 = ( std::abs(b3 - b1) + std::abs(b3 - b2) ) / amrex::Real(32.0);
        amrex::Real a1 = g1 * ( one + (t5*t5) / ((eps + b1) * (eps + b1)) );
        amrex::Real a2 = g2 * ( one + (t5*t5) / ((eps + b2) * (eps + b2)) );
        amrex::Real a3 = g3 * ( one + (t5*t5) / ((eps + b3) * (eps + b3)) );
        amrex::Real asum = a1 + a2 + a3;
        amrex::Real w1 = a1 / asum;
        amrex::Real w2 = a2 / asum;
        amrex::Real w3 = a3 / asum;

        // Taylor expansions
        amrex::Real v1 = (-sm1 + three * s) / two;
        amrex::Real v2 = (  s  +  sp1) / two;
        amrex::Real v3 = (-sm1 + amrex::Real(5.0) * s + two * sp1) / amrex::Real(6.0);

        // Interpolated value
        return ( (w3 / g3) * (v3 - g1 * v1 - g2 * v2) + w1 * v1 + w2 * v2 );
    }

private:
    amrex::Array4<const amrex::Real> m_phi;   // Quantity to interpolate
    const amrex::Real eps=amrex::Real(1.0e-40);
    static constexpr amrex::Real g1=(one/three);
    static constexpr amrex::Real g2=(one/three);
    static constexpr amrex::Real g3=(one/three);
};

/**
 * Interpolation operators used for WENO_Z-5 scheme
 */
struct WENO_Z5
{
    WENO_Z5 (const amrex::Array4<const amrex::Real>& phi,
             const amrex::Real /*upw_frac*/)
        : m_phi(phi) {}

    AMREX_GPU_DEVICE
    AMREX_FORCE_INLINE
    void
    InterpolateInX (const int& i,
                    const int& j,
                    const int& k,
                    const int& qty_index,
                    amrex::Real& val_lo,
                    amrex::Real upw_lo) const
    {
        // Upwinding flags
        if (upw_lo != zero) upw_lo = (upw_lo > 0) ? one : -one;

        // Data to interpolate on
        amrex::Real sp2 = m_phi(i+2, j  , k  , qty_index);
        amrex::Real sp1 = m_phi(i+1, j  , k  , qty_index);
        amrex::Real s   = m_phi(i  , j  , k  , qty_index);
        amrex::Real sm1 = m_phi(i-1, j  , k  , qty_index);
        amrex::Real sm2 = m_phi(i-2, j  , k  , qty_index);
        amrex::Real sm3 = m_phi(i-3, j  , k  , qty_index);

        // Left and right fluxes
        amrex::Real Fhi = Evaluate(sm3,sm2,sm1,s  ,sp1);
        amrex::Real Flo = Evaluate(sp2,sp1,s  ,sm1,sm2);

        // Numerical flux
        val_lo = (one + upw_lo)/two * Fhi + (one - upw_lo)/two * Flo;
    }

    AMREX_GPU_DEVICE
    AMREX_FORCE_INLINE
    void
    InterpolateInY (const int& i,
                    const int& j,
                    const int& k,
                    const int& qty_index,
                    amrex::Real& val_lo,
                    amrex::Real upw_lo) const
    {
        // Upwinding flags
        if (upw_lo != zero) upw_lo = (upw_lo > 0) ? one : -one;

        // Data to interpolate on
        amrex::Real sp2 = m_phi(i  , j+2, k  , qty_index);
        amrex::Real sp1 = m_phi(i  , j+1, k  , qty_index);
        amrex::Real s   = m_phi(i  , j  , k  , qty_index);
        amrex::Real sm1 = m_phi(i  , j-1, k  , qty_index);
        amrex::Real sm2 = m_phi(i  , j-2, k  , qty_index);
        amrex::Real sm3 = m_phi(i  , j-3, k  , qty_index);

        // Left and right fluxes
        amrex::Real Fhi = Evaluate(sm3,sm2,sm1,s  ,sp1);
        amrex::Real Flo = Evaluate(sp2,sp1,s  ,sm1,sm2);

        // Numerical flux
        val_lo = (one + upw_lo)/two * Fhi + (one - upw_lo)/two * Flo;
    }

    AMREX_GPU_DEVICE
    AMREX_FORCE_INLINE
    void
    InterpolateInZ (const int& i,
                    const int& j,
                    const int& k,
                    const int& qty_index,
                    amrex::Real& val_lo,
                    amrex::Real upw_lo) const
    {
        // Upwinding flags
        if (upw_lo != zero) upw_lo = (upw_lo > 0) ? one : -one;

        // Data to interpolate on
        amrex::Real sp2 = m_phi(i  , j  , k+2, qty_index);
        amrex::Real sp1 = m_phi(i  , j  , k+1, qty_index);
        amrex::Real s   = m_phi(i  , j  , k  , qty_index);
        amrex::Real sm1 = m_phi(i  , j  , k-1, qty_index);
        amrex::Real sm2 = m_phi(i  , j  , k-2, qty_index);
        amrex::Real sm3 = m_phi(i  , j  , k-3, qty_index);

        // Left and right fluxes
        amrex::Real Fhi = Evaluate(sm3,sm2,sm1,s  ,sp1);
        amrex::Real Flo = Evaluate(sp2,sp1,s  ,sm1,sm2);

        // Numerical flux
        val_lo = (one + upw_lo)/two * Fhi + (one - upw_lo)/two * Flo;
    }

    AMREX_GPU_DEVICE
    AMREX_FORCE_INLINE
    amrex::Real
    Evaluate (const amrex::Real& sm2,
              const amrex::Real& sm1,
              const amrex::Real& s  ,
              const amrex::Real& sp1,
              const amrex::Real& sp2) const
    {
        // Smoothing factors
        amrex::Real b1 = c1 * (sm2 - two * sm1 + s) * (sm2 - two * sm1 + s) +
                       fourth * (sm2 - amrex::Real(4.0) * sm1 + three * s) * (sm2 - amrex::Real(4.0) * sm1 + three * s);
        amrex::Real b2 = c1 * (sm1 - two * s + sp1) * (sm1 - two * s + sp1) +
                       fourth * (sm1 - sp1) * (sm1 - sp1);
        amrex::Real b3 = c1 * (s - two * sp1 + sp2) * (s - two * sp1 + sp2) +
                       fourth * (three * s - amrex::Real(4.0) * sp1 + sp2) * (three * s - amrex::Real(4.0) * sp1 + sp2);

        // Weight factors
        amrex::Real t5 = std::abs(b3 - b1);
        amrex::Real w1 = g1 * ( one + (t5*t5) / ((eps + b1) * (eps + b1)) );
        amrex::Real w2 = g2 * ( one + (t5*t5) / ((eps + b2) * (eps + b2)) );
        amrex::Real w3 = g3 * ( one + (t5*t5) / ((eps + b3) * (eps + b3)) );

        // Weight factor norm
        amrex::Real wsum = w1 + w2 + w3;

        // Taylor expansions
        amrex::Real v1 = two * sm2 - amrex::Real(7.0) * sm1 + amrex::Real(11.0) * s;
        amrex::Real v2 = -sm1 + amrex::Real(5.0) * s + two * sp1;
        amrex::Real v3 = two * s + amrex::Real(5.0) * sp1 - sp2;

        // Interpolated value
        return ( (w1 * v1 + w2 * v2 + w3 * v3) / (amrex::Real(6.0) * wsum) );
    }

private:
    amrex::Array4<const amrex::Real> m_phi;   // Quantity to interpolate
    const amrex::Real eps=amrex::Real(1.0e-40);
    static constexpr amrex::Real c1=(amrex::Real(13.0)/amrex::Real(12.0));
    static constexpr amrex::Real g1=(one/amrex::Real(10.0));
    static constexpr amrex::Real g2=(three/amrex::Real(5.0));
    static constexpr amrex::Real g3=(three/amrex::Real(10.0));
};

/**
 * Interpolation operators used for WENO_Z-7 scheme
 */
struct WENO_Z7
{
    WENO_Z7 (const amrex::Array4<const amrex::Real>& phi,
             const amrex::Real /*upw_frac*/)
        : m_phi(phi) {}

    AMREX_GPU_DEVICE
    AMREX_FORCE_INLINE
    void
    InterpolateInX (const int& i,
                    const int& j,
                    const int& k,
                    const int& qty_index,
                    amrex::Real& val_lo,
                    amrex::Real upw_lo) const
    {
        // Upwinding flags
        if (upw_lo != zero) upw_lo = (upw_lo > 0) ? one : -one;

        // Data to interpolate on
        amrex::Real sp3 = m_phi(i+3, j  , k  , qty_index);
        amrex::Real sp2 = m_phi(i+2, j  , k  , qty_index);
        amrex::Real sp1 = m_phi(i+1, j  , k  , qty_index);
        amrex::Real s   = m_phi(i  , j  , k  , qty_index);
        amrex::Real sm1 = m_phi(i-1, j  , k  , qty_index);
        amrex::Real sm2 = m_phi(i-2, j  , k  , qty_index);
        amrex::Real sm3 = m_phi(i-3, j  , k  , qty_index);
        amrex::Real sm4 = m_phi(i-4, j  , k  , qty_index);

        // Left and right fluxes
        amrex::Real Fhi = Evaluate(sm4,sm3,sm2,sm1,s  ,sp1,sp2);
        amrex::Real Flo = Evaluate(sp3,sp2,sp1,s  ,sm1,sm2,sm3);

        // Numerical flux
        val_lo = (one + upw_lo)/two * Fhi + (one - upw_lo)/two * Flo;
    }

    AMREX_GPU_DEVICE
    AMREX_FORCE_INLINE
    void
    InterpolateInY (const int& i,
                    const int& j,
                    const int& k,
                    const int& qty_index,
                    amrex::Real& val_lo,
                    amrex::Real upw_lo) const
    {
        // Upwinding flags
        if (upw_lo != zero) upw_lo = (upw_lo > 0) ? one : -one;

        // Data to interpolate on
        amrex::Real sp3 = m_phi(i  , j+3, k  , qty_index);
        amrex::Real sp2 = m_phi(i  , j+2, k  , qty_index);
        amrex::Real sp1 = m_phi(i  , j+1, k  , qty_index);
        amrex::Real s   = m_phi(i  , j  , k  , qty_index);
        amrex::Real sm1 = m_phi(i  , j-1, k  , qty_index);
        amrex::Real sm2 = m_phi(i  , j-2, k  , qty_index);
        amrex::Real sm3 = m_phi(i  , j-3, k  , qty_index);
        amrex::Real sm4 = m_phi(i  , j-4, k  , qty_index);

        // Left and right fluxes
        amrex::Real Fhi = Evaluate(sm4,sm3,sm2,sm1,s  ,sp1,sp2);
        amrex::Real Flo = Evaluate(sp3,sp2,sp1,s  ,sm1,sm2,sm3);

        // Numerical flux
        val_lo = (one + upw_lo)/two * Fhi + (one - upw_lo)/two * Flo;
    }

    AMREX_GPU_DEVICE
    AMREX_FORCE_INLINE
    void
    InterpolateInZ (const int& i,
                    const int& j,
                    const int& k,
                    const int& qty_index,
                    amrex::Real& val_lo,
                    amrex::Real upw_lo) const
    {
        // Upwinding flags
        if (upw_lo != zero) upw_lo = (upw_lo > 0) ? one : -one;

        // Data to interpolate on
        amrex::Real sp3 = m_phi(i  , j  , k+3, qty_index);
        amrex::Real sp2 = m_phi(i  , j  , k+2, qty_index);
        amrex::Real sp1 = m_phi(i  , j  , k+1, qty_index);
        amrex::Real s   = m_phi(i  , j  , k  , qty_index);
        amrex::Real sm1 = m_phi(i  , j  , k-1, qty_index);
        amrex::Real sm2 = m_phi(i  , j  , k-2, qty_index);
        amrex::Real sm3 = m_phi(i  , j  , k-3, qty_index);
        amrex::Real sm4 = m_phi(i  , j  , k-4, qty_index);

        // Left and right fluxes
        amrex::Real Fhi = Evaluate(sm4,sm3,sm2,sm1,s  ,sp1,sp2);
        amrex::Real Flo = Evaluate(sp3,sp2,sp1,s  ,sm1,sm2,sm3);

        // Numerical flux
        val_lo = (one + upw_lo)/two * Fhi + (one - upw_lo)/two * Flo;
    }

    AMREX_GPU_DEVICE
    AMREX_FORCE_INLINE
    amrex::Real
    Evaluate (const amrex::Real& sm3,
              const amrex::Real& sm2,
              const amrex::Real& sm1,
              const amrex::Real& s  ,
              const amrex::Real& sp1,
              const amrex::Real& sp2,
              const amrex::Real& sp3) const
    {
        // Smoothing factors
        amrex::Real b1 = ( sm3 * sm3 * amrex::Real(6649.)/amrex::Real(2880.0)
                         - sm3 * sm2 * amrex::Real(2623.)/amrex::Real(160.0)
                         + sm3 * sm1 * amrex::Real(9449.)/amrex::Real(480.0)
                         - sm3 * s   * amrex::Real(11389.)/amrex::Real(1440.0)
                         + sm2 * sm2 * amrex::Real(28547.)/amrex::Real(960.0)
                         - sm2 * sm1 * amrex::Real(35047.)/amrex::Real(480.0)
                         + sm2 * s   * amrex::Real(14369.)/amrex::Real(480.0)
                         + sm1 * sm1 * amrex::Real(44747.)/amrex::Real(960.0)
                         - sm1 * s   * amrex::Real(6383.)/amrex::Real(160.0)
                         + s   * s   * amrex::Real(25729.)/amrex::Real(2880.0) );
        amrex::Real b2 = ( sm2 * sm2 * 3169/amrex::Real(2880.0)
                         - sm2 * sm1 * 3229/amrex::Real(480.0)
                         + sm2 * s   * 3169/amrex::Real(480.0)
                         - sm2 * sp1 * 2989/amrex::Real(1440.0)
                         + sm1 * sm1 * 11147/amrex::Real(960.0)
                         - sm1 * s   * 11767/amrex::Real(480.0)
                         + sm1 * sp1 * 1283/amrex::Real(160.0)
                         + s   * s   * 13667/amrex::Real(960.0)
                         - s   * sp1 * 5069/amrex::Real(480.0)
                         + sp1 * sp1 * 6649/amrex::Real(2880.0) );
        amrex::Real b3 = ( sm1 * sm1 * amrex::Real(6649.)/amrex::Real(2880.0)
                         - sm1 * s   * amrex::Real(5069.)/amrex::Real(480.0)
                         + sm1 * sp1 * amrex::Real(1283.)/amrex::Real(160.0)
                         - sm1 * sp2 * amrex::Real(2989.)/amrex::Real(1440.0)
                         + s   * s   * amrex::Real(13667.)/amrex::Real(960.0)
                         - s   * sp1 * amrex::Real(11767.)/amrex::Real(480.0)
                         + s   * sp2 * amrex::Real(3169.)/amrex::Real(480.0)
                         + sp1 * sp1 * amrex::Real(11147.)/amrex::Real(960.0)
                         - sp1 * sp2 * amrex::Real(3229.)/amrex::Real(480.0)
                         + sp2 * sp2 * amrex::Real(3169.)/amrex::Real(2880.0) );
        amrex::Real b4 = ( s   * s   * amrex::Real(25729.)/amrex::Real(2880.)
                         - s   * sp1 * amrex::Real(6383.)/amrex::Real(160.)
                         + s   * sp2 * amrex::Real(14369.)/amrex::Real(480.)
                         - s   * sp3 * amrex::Real(11389.)/amrex::Real(1440.)
                         + sp1 * sp1 * amrex::Real(44747.)/amrex::Real(960.)
                         - sp1 * sp2 * amrex::Real(35047.)/amrex::Real(480.)
                         + sp1 * sp3 * amrex::Real(9449.)/amrex::Real(480.)
                         + sp2 * sp2 * amrex::Real(28547.)/amrex::Real(960.)
                         - sp2 * sp3 * amrex::Real(2623.)/amrex::Real(160.)
                         + sp3 * sp3 * amrex::Real(6649.)/amrex::Real(2880.) );

        // Weight factors
        amrex::Real t5 = std::abs(b1 - b2 - b3 + b4);
        amrex::Real w1 = g1 * ( one + (t5*t5) / ((eps + b1) * (eps + b1)) );
        amrex::Real w2 = g2 * ( one + (t5*t5) / ((eps + b2) * (eps + b2)) );
        amrex::Real w3 = g3 * ( one + (t5*t5) / ((eps + b3) * (eps + b3)) );
        amrex::Real w4 = g4 * ( one + (t5*t5) / ((eps + b4) * (eps + b4)) );

        // Weight factor norm
        amrex::Real wsum = w1 + w2 + w3 + w4;

        // Taylor expansions
        amrex::Real v1 = (-amrex::Real(0.3125))*sm3 + ( amrex::Real(1.3125))*sm2 + (-amrex::Real(2.1875))*sm1 + ( amrex::Real(2.1875))*s;
        amrex::Real v2 = ( amrex::Real(0.0625))*sm2 + (-amrex::Real(0.3125))*sm1 + ( amrex::Real(0.9375))*s   + ( amrex::Real(0.3125))*sp1;
        amrex::Real v3 = (-amrex::Real(0.0625))*sm1 + ( amrex::Real(0.5625))*s   + ( amrex::Real(0.5625))*sp1 + (-amrex::Real(0.0625))*sp2;
        amrex::Real v4 = ( amrex::Real(0.3125))*s   + ( amrex::Real(0.9375))*sp1 + (-amrex::Real(0.3125))*sp2 + ( amrex::Real(0.0625))*sp3;

        // Interpolated value
        return ( (w1 * v1 + w2 * v2 + w3 * v3 + w4 * v4) / (wsum) );
    }

private:
    amrex::Array4<const amrex::Real> m_phi;   // Quantity to interpolate
    const amrex::Real eps=amrex::Real(1.0e-40);
    static constexpr amrex::Real g1=( one/amrex::Real(64.0));
    static constexpr amrex::Real g2=(amrex::Real(21.0)/amrex::Real(64.0));
    static constexpr amrex::Real g3=(amrex::Real(35.0)/amrex::Real(64.0));
    static constexpr amrex::Real g4=( amrex::Real(7.0)/amrex::Real(64.0));
};
#endif
