#ifndef ERF_INTERPOLATE_WENO_Z_H_
#define ERF_INTERPOLATE_WENO_Z_H_

#include "ERF_DataStruct.H"

/**
 * Interpolation operators used for WENO_Z-3 scheme
 */
struct WENO_Z3
{
    AMREX_GPU_HOST_DEVICE
    AMREX_FORCE_INLINE
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
        if (upw_lo != amrex::Real(0)) upw_lo = (upw_lo > 0) ? amrex::Real(1) : -amrex::Real(1);

        // Data to interpolate on
        amrex::Real sp1 = m_phi(i+1, j  , k  , qty_index);
        amrex::Real s   = m_phi(i  , j  , k  , qty_index);
        amrex::Real sm1 = m_phi(i-1, j  , k  , qty_index);
        amrex::Real sm2 = m_phi(i-2, j  , k  , qty_index);

        // Left and right fluxes
        amrex::Real Fhi = Evaluate(sm2,sm1,s  );
        amrex::Real Flo = Evaluate(sp1,s  ,sm1);

        // Numerical flux
        val_lo = (amrex::Real(1) + upw_lo)/amrex::Real(2) * Fhi + (amrex::Real(1) - upw_lo)/amrex::Real(2) * Flo;
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
        if (upw_lo != amrex::Real(0)) upw_lo = (upw_lo > 0) ? amrex::Real(1) : -amrex::Real(1);

        // Data to interpolate on
        amrex::Real sp1 = m_phi(i  , j+1, k  , qty_index);
        amrex::Real s   = m_phi(i  , j  , k  , qty_index);
        amrex::Real sm1 = m_phi(i  , j-1, k  , qty_index);
        amrex::Real sm2 = m_phi(i  , j-2, k  , qty_index);

        // Left and right fluxes
        amrex::Real Fhi = Evaluate(sm2,sm1,s  );
        amrex::Real Flo = Evaluate(sp1,s  ,sm1);

        // Numerical flux
        val_lo = (amrex::Real(1) + upw_lo)/amrex::Real(2) * Fhi + (amrex::Real(1) - upw_lo)/amrex::Real(2) * Flo;
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
        if (upw_lo != amrex::Real(0)) upw_lo = (upw_lo > 0) ? amrex::Real(1) : -amrex::Real(1);

        // Data to interpolate on
        amrex::Real sp1 = m_phi(i  , j  , k+1, qty_index);
        amrex::Real s   = m_phi(i  , j  , k  , qty_index);
        amrex::Real sm1 = m_phi(i  , j  , k-1, qty_index);
        amrex::Real sm2 = m_phi(i  , j  , k-2, qty_index);

        // Left and right fluxes
        amrex::Real Fhi = Evaluate(sm2,sm1,s  );
        amrex::Real Flo = Evaluate(sp1,s  ,sm1);

        // Numerical flux
        val_lo = (amrex::Real(1) + upw_lo)/amrex::Real(2) * Fhi + (amrex::Real(1) - upw_lo)/amrex::Real(2) * Flo;
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
        amrex::Real w1 = g1 * ( amrex::Real(1) + (t5*t5) / ((eps + b1) * (eps + b1)) );
        amrex::Real w2 = g2 * ( amrex::Real(1) + (t5*t5) / ((eps + b2) * (eps + b2)) );

        // Weight factor norm
        amrex::Real wsum = w1 + w2;

        // Taylor expansions
        amrex::Real v1 = -sm1 + amrex::Real(3) * s;
        amrex::Real v2 =  s   + sp1;

        // Interpolated value
        return ( (w1 * v1 + w2 * v2) / (amrex::Real(2) * wsum) );
    }

private:
    amrex::Array4<const amrex::Real> m_phi;   // Quantity to interpolate
#ifdef AMREX_USE_FLOAT
    const amrex::Real eps=amrex::Real(1.0e-12);
#else
    const amrex::Real eps=amrex::Real(1.0e-40);
#endif
    static constexpr amrex::Real g1=(amrex::Real(1)/amrex::Real(3));
    static constexpr amrex::Real g2=(amrex::Real(2)/amrex::Real(3));
};

/**
 * Interpolation operators used for WENO_MZQ3 scheme
 */
struct WENO_MZQ3
{
    AMREX_GPU_HOST_DEVICE
    AMREX_FORCE_INLINE
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
        if (upw_lo != amrex::Real(0)) upw_lo = (upw_lo > 0) ? amrex::Real(1) : -amrex::Real(1);

        // Data to interpolate on
        amrex::Real sp1 = m_phi(i+1, j  , k  , qty_index);
        amrex::Real s   = m_phi(i  , j  , k  , qty_index);
        amrex::Real sm1 = m_phi(i-1, j  , k  , qty_index);
        amrex::Real sm2 = m_phi(i-2, j  , k  , qty_index);

        // Left and right fluxes
        amrex::Real Fhi = Evaluate(sm2,sm1,s  );
        amrex::Real Flo = Evaluate(sp1,s  ,sm1);

        // Numerical flux
        val_lo = (amrex::Real(1) + upw_lo)/amrex::Real(2) * Fhi + (amrex::Real(1) - upw_lo)/amrex::Real(2) * Flo;
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
        if (upw_lo != amrex::Real(0)) upw_lo = (upw_lo > 0) ? amrex::Real(1) : -amrex::Real(1);

        // Data to interpolate on
        amrex::Real sp1 = m_phi(i  , j+1, k  , qty_index);
        amrex::Real s   = m_phi(i  , j  , k  , qty_index);
        amrex::Real sm1 = m_phi(i  , j-1, k  , qty_index);
        amrex::Real sm2 = m_phi(i  , j-2, k  , qty_index);

        // Left and right fluxes
        amrex::Real Fhi = Evaluate(sm2,sm1,s  );
        amrex::Real Flo = Evaluate(sp1,s  ,sm1);

        // Numerical flux
        val_lo = (amrex::Real(1) + upw_lo)/amrex::Real(2) * Fhi + (amrex::Real(1) - upw_lo)/amrex::Real(2) * Flo;
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
        if (upw_lo != amrex::Real(0)) upw_lo = (upw_lo > 0) ? amrex::Real(1) : -amrex::Real(1);

        // Data to interpolate on
        amrex::Real sp1 = m_phi(i  , j  , k+1, qty_index);
        amrex::Real s   = m_phi(i  , j  , k  , qty_index);
        amrex::Real sm1 = m_phi(i  , j  , k-1, qty_index);
        amrex::Real sm2 = m_phi(i  , j  , k-2, qty_index);

        // Left and right fluxes
        amrex::Real Fhi = Evaluate(sm2,sm1,s  );
        amrex::Real Flo = Evaluate(sp1,s  ,sm1);

        // Numerical flux
        val_lo = (amrex::Real(1) + upw_lo)/amrex::Real(2) * Fhi + (amrex::Real(1) - upw_lo)/amrex::Real(2) * Flo;
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
        amrex::Real b3 = ( (amrex::Real(13.0) / amrex::Real(12.0)) * ((sm1 - amrex::Real(2)*s + sp1)*(sm1 - amrex::Real(2)*s + sp1)) ) + ( ((sp1 - sm1)*(sp1 - sm1)) / amrex::Real(4.0) );

        // Weight factors
        amrex::Real t5 = ( std::abs(b3 - b1) + std::abs(b3 - b2) ) / amrex::Real(32.0);
        amrex::Real a1 = g1 * ( amrex::Real(1) + (t5*t5) / ((eps + b1) * (eps + b1)) );
        amrex::Real a2 = g2 * ( amrex::Real(1) + (t5*t5) / ((eps + b2) * (eps + b2)) );
        amrex::Real a3 = g3 * ( amrex::Real(1) + (t5*t5) / ((eps + b3) * (eps + b3)) );
        amrex::Real asum = a1 + a2 + a3;
        amrex::Real w1 = a1 / asum;
        amrex::Real w2 = a2 / asum;
        amrex::Real w3 = a3 / asum;

        // Taylor expansions
        amrex::Real v1 = (-sm1 + amrex::Real(3) * s) / amrex::Real(2);
        amrex::Real v2 = (  s  +  sp1) / amrex::Real(2);
        amrex::Real v3 = (-sm1 + amrex::Real(5.0) * s + amrex::Real(2) * sp1) / amrex::Real(6.0);

        // Interpolated value
        return ( (w3 / g3) * (v3 - g1 * v1 - g2 * v2) + w1 * v1 + w2 * v2 );
    }

private:
    amrex::Array4<const amrex::Real> m_phi;   // Quantity to interpolate
#ifdef AMREX_USE_FLOAT
    const amrex::Real eps=amrex::Real(1.0e-12);
#else
    const amrex::Real eps=amrex::Real(1.0e-40);
#endif
    static constexpr amrex::Real g1=(amrex::Real(1)/amrex::Real(3));
    static constexpr amrex::Real g2=(amrex::Real(1)/amrex::Real(3));
    static constexpr amrex::Real g3=(amrex::Real(1)/amrex::Real(3));
};

/**
 * Interpolation operators used for WENO_Z-5 scheme
 */
struct WENO_Z5
{
    AMREX_GPU_HOST_DEVICE
    AMREX_FORCE_INLINE
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
        if (upw_lo != amrex::Real(0)) upw_lo = (upw_lo > 0) ? amrex::Real(1) : -amrex::Real(1);

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
        val_lo = (amrex::Real(1) + upw_lo)/amrex::Real(2) * Fhi + (amrex::Real(1) - upw_lo)/amrex::Real(2) * Flo;
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
        if (upw_lo != amrex::Real(0)) upw_lo = (upw_lo > 0) ? amrex::Real(1) : -amrex::Real(1);

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
        val_lo = (amrex::Real(1) + upw_lo)/amrex::Real(2) * Fhi + (amrex::Real(1) - upw_lo)/amrex::Real(2) * Flo;
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
        if (upw_lo != amrex::Real(0)) upw_lo = (upw_lo > 0) ? amrex::Real(1) : -amrex::Real(1);

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
        val_lo = (amrex::Real(1) + upw_lo)/amrex::Real(2) * Fhi + (amrex::Real(1) - upw_lo)/amrex::Real(2) * Flo;
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
        amrex::Real b1 = c1 * (sm2 - amrex::Real(2) * sm1 + s) * (sm2 - amrex::Real(2) * sm1 + s) +
                       amrex::Real(0.25) * (sm2 - amrex::Real(4.0) * sm1 + amrex::Real(3) * s) * (sm2 - amrex::Real(4.0) * sm1 + amrex::Real(3) * s);
        amrex::Real b2 = c1 * (sm1 - amrex::Real(2) * s + sp1) * (sm1 - amrex::Real(2) * s + sp1) +
                       amrex::Real(0.25) * (sm1 - sp1) * (sm1 - sp1);
        amrex::Real b3 = c1 * (s - amrex::Real(2) * sp1 + sp2) * (s - amrex::Real(2) * sp1 + sp2) +
                       amrex::Real(0.25) * (amrex::Real(3) * s - amrex::Real(4.0) * sp1 + sp2) * (amrex::Real(3) * s - amrex::Real(4.0) * sp1 + sp2);

        // Weight factors
        amrex::Real t5 = std::abs(b3 - b1);
        amrex::Real w1 = g1 * ( amrex::Real(1) + (t5*t5) / ((eps + b1) * (eps + b1)) );
        amrex::Real w2 = g2 * ( amrex::Real(1) + (t5*t5) / ((eps + b2) * (eps + b2)) );
        amrex::Real w3 = g3 * ( amrex::Real(1) + (t5*t5) / ((eps + b3) * (eps + b3)) );

        // Weight factor norm
        amrex::Real wsum = w1 + w2 + w3;

        // Taylor expansions
        amrex::Real v1 = amrex::Real(2) * sm2 - amrex::Real(7.0) * sm1 + amrex::Real(11.0) * s;
        amrex::Real v2 = -sm1 + amrex::Real(5.0) * s + amrex::Real(2) * sp1;
        amrex::Real v3 = amrex::Real(2) * s + amrex::Real(5.0) * sp1 - sp2;

        // Interpolated value
        return ( (w1 * v1 + w2 * v2 + w3 * v3) / (amrex::Real(6.0) * wsum) );
    }

private:
    amrex::Array4<const amrex::Real> m_phi;   // Quantity to interpolate
#ifdef AMREX_USE_FLOAT
    const amrex::Real eps=amrex::Real(1.0e-12);
#else
    const amrex::Real eps=amrex::Real(1.0e-40);
#endif
    static constexpr amrex::Real c1=(amrex::Real(13.0)/amrex::Real(12.0));
    static constexpr amrex::Real g1=(amrex::Real(1)/amrex::Real(10.0));
    static constexpr amrex::Real g2=(amrex::Real(3)/amrex::Real(5.0));
    static constexpr amrex::Real g3=(amrex::Real(3)/amrex::Real(10.0));
};

/**
 * Interpolation operators used for WENO_Z-7 scheme
 */
struct WENO_Z7
{
    AMREX_GPU_HOST_DEVICE
    AMREX_FORCE_INLINE
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
        if (upw_lo != amrex::Real(0)) upw_lo = (upw_lo > 0) ? amrex::Real(1) : -amrex::Real(1);

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
        val_lo = (amrex::Real(1) + upw_lo)/amrex::Real(2) * Fhi + (amrex::Real(1) - upw_lo)/amrex::Real(2) * Flo;
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
        if (upw_lo != amrex::Real(0)) upw_lo = (upw_lo > 0) ? amrex::Real(1) : -amrex::Real(1);

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
        val_lo = (amrex::Real(1) + upw_lo)/amrex::Real(2) * Fhi + (amrex::Real(1) - upw_lo)/amrex::Real(2) * Flo;
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
        if (upw_lo != amrex::Real(0)) upw_lo = (upw_lo > 0) ? amrex::Real(1) : -amrex::Real(1);

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
        val_lo = (amrex::Real(1) + upw_lo)/amrex::Real(2) * Fhi + (amrex::Real(1) - upw_lo)/amrex::Real(2) * Flo;
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
        // Finite-volume Jiang-Shu smoothness indicators for cell-average data.
        amrex::Real b1 = ( amrex::Real(547.0) * sm3 * sm3
                 - amrex::Real(3882.0) * sm3 * sm2
                 + amrex::Real(4642.0) * sm3 * sm1
                 - amrex::Real(1854.0) * sm3 * s
                 + amrex::Real(7043.0) * sm2 * sm2
                 - amrex::Real(17246.0) * sm2 * sm1
                 + amrex::Real(7042.0) * sm2 * s
                 + amrex::Real(11003.0) * sm1 * sm1
                 - amrex::Real(9402.0) * sm1 * s
                 + amrex::Real(2107.0) * s * s ) / amrex::Real(240.0);
        amrex::Real b2 = ( amrex::Real(267.0) * sm2 * sm2
                 - amrex::Real(1642.0) * sm2 * sm1
                 + amrex::Real(1602.0) * sm2 * s
                 - amrex::Real(494.0) * sm2 * sp1
                 + amrex::Real(2843.0) * sm1 * sm1
                 - amrex::Real(5966.0) * sm1 * s
                 + amrex::Real(1922.0) * sm1 * sp1
                 + amrex::Real(3443.0) * s * s
                 - amrex::Real(2522.0) * s * sp1
                 + amrex::Real(547.0) * sp1 * sp1 ) / amrex::Real(240.0);
        amrex::Real b3 = ( amrex::Real(547.0) * sm1 * sm1
                 - amrex::Real(2522.0) * sm1 * s
                 + amrex::Real(1922.0) * sm1 * sp1
                 - amrex::Real(494.0) * sm1 * sp2
                 + amrex::Real(3443.0) * s * s
                 - amrex::Real(5966.0) * s * sp1
                 + amrex::Real(1602.0) * s * sp2
                 + amrex::Real(2843.0) * sp1 * sp1
                 - amrex::Real(1642.0) * sp1 * sp2
                 + amrex::Real(267.0) * sp2 * sp2 ) / amrex::Real(240.0);
        amrex::Real b4 = ( amrex::Real(2107.0) * s * s
                 - amrex::Real(9402.0) * s * sp1
                 + amrex::Real(7042.0) * s * sp2
                 - amrex::Real(1854.0) * s * sp3
                 + amrex::Real(11003.0) * sp1 * sp1
                 - amrex::Real(17246.0) * sp1 * sp2
                 + amrex::Real(4642.0) * sp1 * sp3
                 + amrex::Real(7043.0) * sp2 * sp2
                 - amrex::Real(3882.0) * sp2 * sp3
                 + amrex::Real(547.0) * sp3 * sp3 ) / amrex::Real(240.0);

        // Weight factors
        amrex::Real t5 = std::abs(b1 - b2 - b3 + b4);
        amrex::Real w1 = g1 * ( amrex::Real(1) + (t5*t5) / ((eps + b1) * (eps + b1)) );
        amrex::Real w2 = g2 * ( amrex::Real(1) + (t5*t5) / ((eps + b2) * (eps + b2)) );
        amrex::Real w3 = g3 * ( amrex::Real(1) + (t5*t5) / ((eps + b3) * (eps + b3)) );
        amrex::Real w4 = g4 * ( amrex::Real(1) + (t5*t5) / ((eps + b4) * (eps + b4)) );

        // Weight factor norm
        amrex::Real wsum = w1 + w2 + w3 + w4;

        // Finite-volume face-reconstruction coefficients for the left state at i-1/2.
        amrex::Real v1 = (-amrex::Real(1.0)/amrex::Real(4.0))*sm3 + ( amrex::Real(13.0)/amrex::Real(12.0))*sm2 -
                 ( amrex::Real(23.0)/amrex::Real(12.0))*sm1 + ( amrex::Real(25.0)/amrex::Real(12.0))*s;
        amrex::Real v2 = ( amrex::Real(1.0)/amrex::Real(12.0))*sm2 - ( amrex::Real(5.0)/amrex::Real(12.0))*sm1 +
                 ( amrex::Real(13.0)/amrex::Real(12.0))*s   + ( amrex::Real(1.0)/amrex::Real(4.0))*sp1;
        amrex::Real v3 = (-amrex::Real(1.0)/amrex::Real(12.0))*sm1 + ( amrex::Real(7.0)/amrex::Real(12.0))*s +
                 ( amrex::Real(7.0)/amrex::Real(12.0))*sp1 - ( amrex::Real(1.0)/amrex::Real(12.0))*sp2;
        amrex::Real v4 = ( amrex::Real(1.0)/amrex::Real(4.0))*s + ( amrex::Real(13.0)/amrex::Real(12.0))*sp1 -
                 ( amrex::Real(5.0)/amrex::Real(12.0))*sp2 + ( amrex::Real(1.0)/amrex::Real(12.0))*sp3;

        // Interpolated value
        return ( (w1 * v1 + w2 * v2 + w3 * v3 + w4 * v4) / (wsum) );
    }

private:
    amrex::Array4<const amrex::Real> m_phi;   // Quantity to interpolate
#ifdef AMREX_USE_FLOAT
    const amrex::Real eps=amrex::Real(1.0e-12);
#else
    const amrex::Real eps=amrex::Real(1.0e-40);
#endif
    static constexpr amrex::Real g1=( amrex::Real(1.0)/amrex::Real(35.0));
    static constexpr amrex::Real g2=(amrex::Real(12.0)/amrex::Real(35.0));
    static constexpr amrex::Real g3=(amrex::Real(18.0)/amrex::Real(35.0));
    static constexpr amrex::Real g4=( amrex::Real(4.0)/amrex::Real(35.0));
};
#endif
