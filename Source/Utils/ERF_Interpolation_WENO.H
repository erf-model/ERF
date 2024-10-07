#ifndef ERF_INTERPOLATE_WENO_H_
#define ERF_INTERPOLATE_WENO_H_

#include "ERF_DataStruct.H"

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
            val_lo = Evaluate(sp2,sp1,s  ,sm1,sm2);
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
            val_lo = Evaluate(sp2,sp1,s  ,sm1,sm2);
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
            val_lo = Evaluate(sp2,sp1,s  ,sm1,sm2);
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

/**
 * Interpolation operators used for WENO-7 scheme
 */
struct WENO7
{
    WENO7 (const amrex::Array4<const amrex::Real>& phi)
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
        amrex::Real sp3 = m_phi(i+3, j  , k  , qty_index);
        amrex::Real sp2 = m_phi(i+2, j  , k  , qty_index);
        amrex::Real sp1 = m_phi(i+1, j  , k  , qty_index);
        amrex::Real s   = m_phi(i  , j  , k  , qty_index);
        amrex::Real sm1 = m_phi(i-1, j  , k  , qty_index);
        amrex::Real sm2 = m_phi(i-2, j  , k  , qty_index);
        amrex::Real sm3 = m_phi(i-3, j  , k  , qty_index);
        amrex::Real sm4 = m_phi(i-4, j  , k  , qty_index);

        if (upw_lo > tol) {
            val_lo = Evaluate(sm4,sm3,sm2,sm1,s  ,sp1,sp2);
        } else if (upw_lo < -tol) {
            val_lo = Evaluate(sp3,sp2,sp1,s  ,sm1,sm2,sm3);
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
        amrex::Real sp3 = m_phi(i  , j+3, k  , qty_index);
        amrex::Real sp2 = m_phi(i  , j+2, k  , qty_index);
        amrex::Real sp1 = m_phi(i  , j+1, k  , qty_index);
        amrex::Real s   = m_phi(i  , j  , k  , qty_index);
        amrex::Real sm1 = m_phi(i  , j-1, k  , qty_index);
        amrex::Real sm2 = m_phi(i  , j-2, k  , qty_index);
        amrex::Real sm3 = m_phi(i  , j-3, k  , qty_index);
        amrex::Real sm4 = m_phi(i  , j-4, k  , qty_index);

        if (upw_lo > tol) {
            val_lo = Evaluate(sm4,sm3,sm2,sm1,s  ,sp1,sp2);
        } else if (upw_lo < -tol) {
            val_lo = Evaluate(sp3,sp2,sp1,s  ,sm1,sm2,sm3);
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
        amrex::Real sp3 = m_phi(i  , j  , k+2, qty_index);
        amrex::Real sp2 = m_phi(i  , j  , k+2, qty_index);
        amrex::Real sp1 = m_phi(i  , j  , k+1, qty_index);
        amrex::Real s   = m_phi(i  , j  , k  , qty_index);
        amrex::Real sm1 = m_phi(i  , j  , k-1, qty_index);
        amrex::Real sm2 = m_phi(i  , j  , k-2, qty_index);
        amrex::Real sm3 = m_phi(i  , j  , k-3, qty_index);
        amrex::Real sm4 = m_phi(i  , j  , k-3, qty_index);

        if (upw_lo > tol) {
            val_lo = Evaluate(sm4,sm3,sm2,sm1,s  ,sp1,sp2);
        } else if (upw_lo < -tol) {
            val_lo = Evaluate(sp3,sp2,sp1,s  ,sm1,sm2,sm3);
        } else {
            val_lo = 0.5 * (s + sm1);
        }
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
        amrex::Real b1 = ( sm3 * sm3 * 6649./2880.0
                         - sm3 * sm2 * 2623./160.0
                         + sm3 * sm1 * 9449./480.0
                         - sm3 * s   * 11389./1440.0
                         + sm2 * sm2 * 28547./960.0
                         - sm2 * sm1 * 35047./480.0
                         + sm2 * s   * 14369./480.0
                         + sm1 * sm1 * 44747./960.0
                         - sm1 * s   * 6383./160.0
                         + s   * s   * 25729./2880.0 );
        amrex::Real b2 = ( sm2 * sm2 * 3169/2880.0
                         - sm2 * sm1 * 3229/480.0
                         + sm2 * s   * 3169/480.0
                         - sm2 * sp1 * 2989/1440.0
                         + sm1 * sm1 * 11147/960.0
                         - sm1 * s   * 11767/480.0
                         + sm1 * sp1 * 1283/160.0
                         + s   * s   * 13667/960.0
                         - s   * sp1 * 5069/480.0
                         + sp1 * sp1 * 6649/2880.0 );
        amrex::Real b3 = ( sm1 * sm1 * 6649./2880.0
                         - sm1 * s   * 5069./480.0
                         + sm1 * sp1 * 1283./160.0
                         - sm1 * sp2 * 2989./1440.0
                         + s   * s   * 13667./960.0
                         - s   * sp1 * 11767./480.0
                         + s   * sp2 * 3169./480.0
                         + sp1 * sp1 * 11147./960.0
                         - sp1 * sp2 * 3229./480.0
                         + sp2 * sp2 * 3169./2880.0 );
        amrex::Real b4 = ( s   * s   * 25729./2880.
                         - s   * sp1 * 6383./160.
                         + s   * sp2 * 14369./480.
                         - s   * sp3 * 11389./1440.
                         + sp1 * sp1 * 44747./960.
                         - sp1 * sp2 * 35047./480.
                         + sp1 * sp3 * 9449./480.
                         + sp2 * sp2 * 28547./960.
                         - sp2 * sp3 * 2623./160.
                         + sp3 * sp3 * 6649./2880. );

        // Weight factors
        amrex::Real w1 = g1 / ( (eps + b1) * (eps + b1) );
        amrex::Real w2 = g2 / ( (eps + b2) * (eps + b2) );
        amrex::Real w3 = g3 / ( (eps + b3) * (eps + b3) );
        amrex::Real w4 = g4 / ( (eps + b4) * (eps + b4) );

        // Weight factor norm
        amrex::Real wsum = w1 + w2 + w3 + w4;

        // Taylor expansions
        amrex::Real v1 = (-0.3125)*sm3 + ( 1.3125)*sm2 + (-2.1875)*sm1 + ( 2.1875)*s;
        amrex::Real v2 = ( 0.0625)*sm2 + (-0.3125)*sm1 + ( 0.9375)*s   + ( 0.3125)*sp1;
        amrex::Real v3 = (-0.0625)*sm1 + ( 0.5625)*s   + ( 0.5625)*sp1 + (-0.0625)*sp2;
        amrex::Real v4 = ( 0.3125)*s   + ( 0.9375)*sp1 + (-0.3125)*sp2 + ( 0.0625)*sp3;

        // Interpolated value
        return ( (w1 * v1 + w2 * v2 + w3 * v3 + w4 * v4) / (wsum) );
    }

private:
    amrex::Array4<const amrex::Real> m_phi;   // Quantity to interpolate
    const amrex::Real eps=1.0e-8;
    const amrex::Real tol=1.0e-12;
    static constexpr amrex::Real g1=( 1.0/64.0);
    static constexpr amrex::Real g2=(21.0/64.0);
    static constexpr amrex::Real g3=(35.0/64.0);
    static constexpr amrex::Real g4=( 7.0/64.0);
};
#endif
