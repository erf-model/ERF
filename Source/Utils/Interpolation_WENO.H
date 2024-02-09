#ifndef INTERPOLATE_WENO_H_
#define INTERPOLATE_WENO_H_

#include "DataStruct.H"

/**
 * Interpolation operators used for WENO-5 scheme
 */
struct WENO3
{
    WENO3 (const amrex::Array4<const amrex::Real>& phi)
        : m_phi(phi) {}

    AMREX_GPU_DEVICE
    AMREX_FORCE_INLINE
    void
    InterpolateInX (const int& i,
                    const int& j,
                    const int& k,
                    const int& qty_index,
                    amrex::Real& val_lo,
                    amrex::Real upw_lo,
                    const amrex::Real /*upw_frac*/) const
    {
        // Data to interpolate on
        amrex::Real sp1 = m_phi(i+1, j  , k  , qty_index);
        amrex::Real s   = m_phi(i  , j  , k  , qty_index);
        amrex::Real sm1 = m_phi(i-1, j  , k  , qty_index);
        amrex::Real sm2 = m_phi(i-2, j  , k  , qty_index);

        if (upw_lo > tol) {
            val_lo = Evaluate(sm2,sm1,s  );
        } else if (upw_lo < -tol) {
            val_lo = Evaluate(sp1,s  ,sm1);
        } else {
            val_lo = 0.5 * (s + sm1);
        }
    }

    AMREX_GPU_DEVICE
    AMREX_FORCE_INLINE
    void
    InterpolateInY (const int& i,
                    const int& j,
                    const int& k,
                    const int& qty_index,
                    amrex::Real& val_lo,
                    amrex::Real upw_lo,
                    const amrex::Real /*upw_frac*/) const
    {
        // Data to interpolate on
        amrex::Real sp1 = m_phi(i  , j+1, k  , qty_index);
        amrex::Real s   = m_phi(i  , j  , k  , qty_index);
        amrex::Real sm1 = m_phi(i  , j-1, k  , qty_index);
        amrex::Real sm2 = m_phi(i  , j-2, k  , qty_index);

        if (upw_lo > tol) {
            val_lo = Evaluate(sm2,sm1,s  );
        } else if (upw_lo < -tol) {
            val_lo = Evaluate(sp1,s  ,sm1);
        } else {
            val_lo = 0.5 * (s + sm1);
        }
    }

    AMREX_GPU_DEVICE
    AMREX_FORCE_INLINE
    void
    InterpolateInZ (const int& i,
                       const int& j,
                       const int& k,
                       const int& qty_index,
                       amrex::Real& val_lo,
                       amrex::Real upw_lo,
                       const amrex::Real /*upw_frac*/) const
    {
        // Data to interpolate on
        amrex::Real sp1 = m_phi(i  , j  , k+1, qty_index);
        amrex::Real s   = m_phi(i  , j  , k  , qty_index);
        amrex::Real sm1 = m_phi(i  , j  , k-1, qty_index);
        amrex::Real sm2 = m_phi(i  , j  , k-2, qty_index);

        if (upw_lo > tol) {
            val_lo = Evaluate(sm2,sm1,s  );
        } else if (upw_lo < -tol) {
            val_lo = Evaluate(sp1,s  ,sm1);
        } else {
            val_lo = 0.5 * (s + sm1);
        }
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
        amrex::Real w1 = g1 / ( (eps + b1) * (eps + b1) );
        amrex::Real w2 = g2 / ( (eps + b2) * (eps + b2) );

        // Weight factor norm
        amrex::Real wsum = w1 + w2;

        // Taylor expansions
        amrex::Real v1 = -sm1 + 3.0 * s;
        amrex::Real v2 =  s   + sp1;

        // Interpolated value
        return ( (w1 * v1 + w2 * v2) / (2.0 * wsum) );
    }

private:
    amrex::Array4<const amrex::Real> m_phi;   // Quantity to interpolate
    const amrex::Real eps=1.0e-6;
    const amrex::Real tol=1.0e-12;
    static constexpr amrex::Real g1=(1.0/3.0);
    static constexpr amrex::Real g2=(2.0/3.0);
};

/**
 * Interpolation operators used for WENO-5 scheme
 */
struct WENO5
{
    WENO5 (const amrex::Array4<const amrex::Real>& phi)
        : m_phi(phi) {}

    AMREX_GPU_DEVICE
    AMREX_FORCE_INLINE
    void
    InterpolateInX (const int& i,
                    const int& j,
                    const int& k,
                    const int& qty_index,
                    amrex::Real& val_lo,
                    amrex::Real upw_lo,
                    const amrex::Real /*upw_frac*/) const
    {
        // Data to interpolate on
        amrex::Real sp2 = m_phi(i+2, j  , k  , qty_index);
        amrex::Real sp1 = m_phi(i+1, j  , k  , qty_index);
        amrex::Real s   = m_phi(i  , j  , k  , qty_index);
        amrex::Real sm1 = m_phi(i-1, j  , k  , qty_index);
        amrex::Real sm2 = m_phi(i-2, j  , k  , qty_index);
        amrex::Real sm3 = m_phi(i-3, j  , k  , qty_index);

        if (upw_lo > tol) {
            val_lo = Evaluate(sm3,sm2,sm1,s  ,sp1);
        } else if (upw_lo < -tol) {
            val_lo = Evaluate(sp2,sp1,s,sm1,sm2);
        } else {
            val_lo = 0.5 * (s + sm1);
        }
    }

    AMREX_GPU_DEVICE
    AMREX_FORCE_INLINE
    void
    InterpolateInY (const int& i,
                    const int& j,
                    const int& k,
                    const int& qty_index,
                    amrex::Real& val_lo,
                    amrex::Real upw_lo,
                    const amrex::Real /*upw_frac*/) const
    {
        // Data to interpolate on
        amrex::Real sp2 = m_phi(i  , j+2, k  , qty_index);
        amrex::Real sp1 = m_phi(i  , j+1, k  , qty_index);
        amrex::Real s   = m_phi(i  , j  , k  , qty_index);
        amrex::Real sm1 = m_phi(i  , j-1, k  , qty_index);
        amrex::Real sm2 = m_phi(i  , j-2, k  , qty_index);
        amrex::Real sm3 = m_phi(i  , j-3, k  , qty_index);

        if (upw_lo > tol) {
            val_lo = Evaluate(sm3,sm2,sm1,s  ,sp1);
        } else if (upw_lo < -tol) {
            val_lo = Evaluate(sp2,sp1,s,sm1,sm2);
        } else {
            val_lo = 0.5 * (s + sm1);
        }
    }

    AMREX_GPU_DEVICE
    AMREX_FORCE_INLINE
    void
    InterpolateInZ (const int& i,
                       const int& j,
                       const int& k,
                       const int& qty_index,
                       amrex::Real& val_lo,
                       amrex::Real upw_lo,
                       const amrex::Real /*upw_frac*/) const
    {
        // Data to interpolate on
        amrex::Real sp2 = m_phi(i  , j  , k+2, qty_index);
        amrex::Real sp1 = m_phi(i  , j  , k+1, qty_index);
        amrex::Real s   = m_phi(i  , j  , k  , qty_index);
        amrex::Real sm1 = m_phi(i  , j  , k-1, qty_index);
        amrex::Real sm2 = m_phi(i  , j  , k-2, qty_index);
        amrex::Real sm3 = m_phi(i  , j  , k-3, qty_index);

        if (upw_lo > tol) {
            val_lo = Evaluate(sm3,sm2,sm1,s  ,sp1);
        } else if (upw_lo < -tol) {
            val_lo = Evaluate(sp2,sp1,s,sm1,sm2);
        } else {
            val_lo = 0.5 * (s + sm1);
        }
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
        amrex::Real b1 = c1 * (sm2 - 2.0 * sm1 + s) * (sm2 - 2.0 * sm1 + s) +
                       0.25 * (sm2 - 4.0 * sm1 + 3.0 * s) * (sm2 - 4.0 * sm1 + 3.0 * s);
        amrex::Real b2 = c1 * (sm1 - 2.0 * s + sp1) * (sm1 - 2.0 * s + sp1) +
                       0.25 * (sm1 - sp1) * (sm1 - sp1);
        amrex::Real b3 = c1 * (s - 2.0 * sp1 + sp2) * (s - 2.0 * sp1 + sp2) +
                       0.25 * (3.0 * s - 4.0 * sp1 + sp2) * (3.0 * s - 4.0 * sp1 + sp2);

        // Weight factors
        amrex::Real w1 = g1 / ( (eps + b1) * (eps + b1) );
        amrex::Real w2 = g2 / ( (eps + b2) * (eps + b2) );
        amrex::Real w3 = g3 / ( (eps + b3) * (eps + b3) );

        // Weight factor norm
        amrex::Real wsum = w1 + w2 + w3;

        // Taylor expansions
        amrex::Real v1 = 2.0 * sm2 - 7.0 * sm1 + 11.0 * s;
        amrex::Real v2 = -sm1 + 5.0 * s + 2.0 * sp1;
        amrex::Real v3 = 2.0 * s + 5.0 * sp1 - sp2;

        // Interpolated value
        return ( (w1 * v1 + w2 * v2 + w3 * v3) / (6.0 * wsum) );
    }

private:
    amrex::Array4<const amrex::Real> m_phi;   // Quantity to interpolate
    const amrex::Real eps=1.0e-6;
    const amrex::Real tol=1.0e-12;
    static constexpr amrex::Real c1=(13.0/12.0);
    static constexpr amrex::Real g1=(1.0/10.0);
    static constexpr amrex::Real g2=(3.0/5.0);
    static constexpr amrex::Real g3=(3.0/10.0);
};
#endif
