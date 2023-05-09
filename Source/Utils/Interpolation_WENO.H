#ifndef INTERPOLATE_WENO_H_
#define INTERPOLATE_WENO_H_

#include "DataStruct.H"

/**
 * Interpolation operators used for WENO-5 scheme
 */
struct WENO3
{
    WENO3(const amrex::Array4<const amrex::Real>& phi)
        : m_phi(phi) {}

    AMREX_GPU_DEVICE
    AMREX_FORCE_INLINE
    void
    InterpolateInX(const int& i,
                   const int& j,
                   const int& k,
                   const int& qty_index,
                   amrex::Real& val_hi,
                   amrex::Real& val_lo,
                   amrex::Real /*upw_hi*/,
                   amrex::Real /*upw_lo*/)
    {
        // Data to interpolate on
        amrex::Real sp1 = m_phi(i+1, j  , k  , qty_index);
        amrex::Real s   = m_phi(i  , j  , k  , qty_index);
        amrex::Real sm1 = m_phi(i-1, j  , k  , qty_index);

        val_hi = Evaluate(sm1,s,sp1);
        val_lo = Evaluate(sp1,s,sm1);
    }

    AMREX_GPU_DEVICE
    AMREX_FORCE_INLINE
    void
    InterpolateInY(const int& i,
                   const int& j,
                   const int& k,
                   const int& qty_index,
                   amrex::Real& val_hi,
                   amrex::Real& val_lo,
                   amrex::Real /*upw_hi*/,
                   amrex::Real /*upw_lo*/)
    {
        // Data to interpolate on
        amrex::Real sp1 = m_phi(i  , j+1, k  , qty_index);
        amrex::Real s   = m_phi(i  , j  , k  , qty_index);
        amrex::Real sm1 = m_phi(i  , j-1, k  , qty_index);

        val_hi = Evaluate(sm1,s,sp1);
        val_lo = Evaluate(sp1,s,sm1);
    }

    AMREX_GPU_DEVICE
    AMREX_FORCE_INLINE
    void
    InterpolateInZ_lo(const int& i,
                      const int& j,
                      const int& k,
                      const int& qty_index,
                      amrex::Real& val_lo,
                      amrex::Real /*upw_lo*/)
    {
        // Data to interpolate on
        amrex::Real sp1 = m_phi(i  , j  , k+1, qty_index);
        amrex::Real s   = m_phi(i  , j  , k  , qty_index);
        amrex::Real sm1 = m_phi(i  , j  , k-1, qty_index);

        val_lo = Evaluate(sp1,s,sm1);
    }

    AMREX_GPU_DEVICE
    AMREX_FORCE_INLINE
    void
    InterpolateInZ_hi(const int& i,
                      const int& j,
                      const int& k,
                      const int& qty_index,
                      amrex::Real& val_hi,
                      amrex::Real /*upw_hi*/)
    {
        // Data to interpolate on
        amrex::Real sp1 = m_phi(i  , j  , k+1, qty_index);
        amrex::Real s   = m_phi(i  , j  , k  , qty_index);
        amrex::Real sm1 = m_phi(i  , j  , k-1, qty_index);

        val_hi = Evaluate(sm1,s,sp1);
    }

    AMREX_GPU_DEVICE
    AMREX_FORCE_INLINE
    amrex::Real
    Evaluate(const amrex::Real& sm1,
             const amrex::Real& s  ,
             const amrex::Real& sp1)
    {
        // Smoothing factors
        b1 = (s - sm1) * (s - sm1);
        b2 = (sp1 - s) * (sp1 - s);

        // Weight factors
        w1 = g1 / ( (eps + b1) * (eps + b1) );
        w2 = g2 / ( (eps + b2) * (eps + b2) );

        // Weight factor norm
        wsum = w1 + w2;

        // Taylor expansions
        v1 = -sm1 + 3.0 * s;
        v2 =  s   + sp1;

        // Interpolated value
        return ( (w1 * v1 + w2 * v2) / (2.0 * wsum) );
    }

private:
    const amrex::Array4<const amrex::Real>& m_phi;
    amrex::Real eps=1.0e-6;
    amrex::Real b1=0. , b2=0.;
    amrex::Real w1=0. , w2=0., wsum=0.;
    amrex::Real v1=0. , v2=0.;
    static constexpr amrex::Real g1=(1.0/3.0);
    static constexpr amrex::Real g2=(2.0/3.0);
};

/**
 * Interpolation operators used for WENO-5 scheme
 */
struct WENO5
{
    WENO5(const amrex::Array4<const amrex::Real>& phi)
        : m_phi(phi) {}

    AMREX_GPU_DEVICE
    AMREX_FORCE_INLINE
    void
    InterpolateInX(const int& i,
                   const int& j,
                   const int& k,
                   const int& qty_index,
                   amrex::Real& val_hi,
                   amrex::Real& val_lo,
                   amrex::Real /*upw_hi*/,
                   amrex::Real /*upw_lo*/)
    {
        // Data to interpolate on
        amrex::Real sp2 = m_phi(i+2, j  , k  , qty_index);
        amrex::Real sp1 = m_phi(i+1, j  , k  , qty_index);
        amrex::Real s   = m_phi(i  , j  , k  , qty_index);
        amrex::Real sm1 = m_phi(i-1, j  , k  , qty_index);
        amrex::Real sm2 = m_phi(i-2, j  , k  , qty_index);

        val_hi = Evaluate(sm2,sm1,s,sp1,sp2);
        val_lo = Evaluate(sp2,sp1,s,sm1,sm2);
    }

    AMREX_GPU_DEVICE
    AMREX_FORCE_INLINE
    void
    InterpolateInY(const int& i,
                   const int& j,
                   const int& k,
                   const int& qty_index,
                   amrex::Real& val_hi,
                   amrex::Real& val_lo,
                   amrex::Real /*upw_hi*/,
                   amrex::Real /*upw_lo*/)
    {
        // Data to interpolate on
        amrex::Real sp2 = m_phi(i  , j+2, k  , qty_index);
        amrex::Real sp1 = m_phi(i  , j+1, k  , qty_index);
        amrex::Real s   = m_phi(i  , j  , k  , qty_index);
        amrex::Real sm1 = m_phi(i  , j-1, k  , qty_index);
        amrex::Real sm2 = m_phi(i  , j-2, k  , qty_index);

        val_hi = Evaluate(sm2,sm1,s,sp1,sp2);
        val_lo = Evaluate(sp2,sp1,s,sm1,sm2);
    }

    AMREX_GPU_DEVICE
    AMREX_FORCE_INLINE
    void
    InterpolateInZ_lo(const int& i,
                      const int& j,
                      const int& k,
                      const int& qty_index,
                      amrex::Real& val_lo,
                      amrex::Real /*upw_lo*/)
    {
        // Data to interpolate on
        amrex::Real sp2 = m_phi(i  , j  , k+2, qty_index);
        amrex::Real sp1 = m_phi(i  , j  , k+1, qty_index);
        amrex::Real s   = m_phi(i  , j  , k  , qty_index);
        amrex::Real sm1 = m_phi(i  , j  , k-1, qty_index);
        amrex::Real sm2 = m_phi(i  , j  , k-2, qty_index);

        val_lo = Evaluate(sp2,sp1,s,sm1,sm2);
    }

    AMREX_GPU_DEVICE
    AMREX_FORCE_INLINE
    void
    InterpolateInZ_hi(const int& i,
                      const int& j,
                      const int& k,
                      const int& qty_index,
                      amrex::Real& val_hi,
                      amrex::Real /*upw_hi*/)
    {
        // Data to interpolate on
        amrex::Real sp2 = m_phi(i  , j  , k+2, qty_index);
        amrex::Real sp1 = m_phi(i  , j  , k+1, qty_index);
        amrex::Real s   = m_phi(i  , j  , k  , qty_index);
        amrex::Real sm1 = m_phi(i  , j  , k-1, qty_index);
        amrex::Real sm2 = m_phi(i  , j  , k-2, qty_index);

        val_hi = Evaluate(sm2,sm1,s,sp1,sp2);
    }

    AMREX_GPU_DEVICE
    AMREX_FORCE_INLINE
    amrex::Real
    Evaluate(const amrex::Real& sm2,
             const amrex::Real& sm1,
             const amrex::Real& s  ,
             const amrex::Real& sp1,
             const amrex::Real& sp2)
    {
        // Smoothing factors
        b1 = c1 * (sm2 - 2.0 * sm1 + s) * (sm2 - 2.0 * sm1 + s) +
             0.25 * (sm2 - 4.0 * sm1 + 3.0 * s) * (sm2 - 4.0 * sm1 + 3.0 * s);
        b2 = c1 * (sm1 - 2.0 * s + sp1) * (sm1 - 2.0 * s + sp1) +
             0.25 * (sm1 - sp1) * (sm1 - sp1);
        b3 = c1 * (s - 2.0 * sp1 + sp2) * (s - 2.0 * sp1 + sp2) +
             0.25 * (3.0 * s - 4.0 * sp1 + sp2) * (3.0 * s - 4.0 * sp1 + sp2);

        // Weight factors
        w1 = g1 / ( (eps + b1) * (eps + b1) );
        w2 = g2 / ( (eps + b2) * (eps + b2) );
        w3 = g3 / ( (eps + b3) * (eps + b3) );

        // Weight factor norm
        wsum = w1 + w2 + w3;

        // Taylor expansions
        v1 = 2.0 * sm2 - 7.0 * sm1 + 11.0 * s;
        v2 = -sm1 + 5.0 * s + 2.0 * sp1;
        v3 = 2.0 * s + 5.0 * sp1 - sp2;

        // Interpolated value
        return ( (w1 * v1 + w2 * v2 + w3 * v3) / (6.0 * wsum) );
    }

private:
    const amrex::Array4<const amrex::Real>& m_phi;
    amrex::Real eps=1.0e-6;
    amrex::Real b1=0. , b2=0. , b3=0. ;
    amrex::Real w1=0. , w2=0. , w3=0. , wsum=0.;
    amrex::Real v1=0. , v2=0. , v3=0. ;
    static constexpr amrex::Real c1=(13.0/12.0);
    static constexpr amrex::Real g1=(1.0/10.0);
    static constexpr amrex::Real g2=(3.0/5.0);
    static constexpr amrex::Real g3=(3.0/10.0);
};
#endif
