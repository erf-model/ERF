
#ifndef VelPlaneAverage_H
#define VelPlaneAverage_H

#include "ERF.H"
#include "DirectionSelector.H"

using namespace amrex;

class VelPlaneAverage {
public:
    VelPlaneAverage(const Vector<MultiFab*>& mfab_in,
                    Geometry geom_in, int axis_in,
                    bool compute_deriv);
    ~VelPlaneAverage() = default;

    const Vector<MultiFab*> field() const { return m_mfab; };

    void operator()();
private:
    int m_ncomp; /** number of average components */

    Vector<Real> m_line_xcentroid; /** line storage for centroids of each cell along a line*/
    Real m_dx;  /** mesh spacing in axis direction*/
    Real m_xlo; /** bottom of domain in axis direction */
    int m_ncell_plane; /** number of cells in plane */
    int m_ncell_line;  /** number of cells along line */

    const Vector<MultiFab*> m_mfab;
    Geometry m_geom;
    const int m_axis;
    const bool m_comp_deriv;

    Vector<Real> m_line_hvelmag_average; /** line storage for the average horizontal velocity magnitude */
    //! line storage for the derivative of average horizontal velocity magnitude
    Vector<Real> m_line_hvelmag_deriv;

public:
    template <typename IndexSelector>
    void compute_hvelmag_averages(const IndexSelector& idx_op,
                                  const int h1_idx,
                                  const int h2_idx,
                                  const Vector<MultiFab*> mfab);

    /** return vector containing horizontal velocity magnitude average */
    const Vector<Real>& line_hvelmag_average()
    {
        return m_line_hvelmag_average;
    };

    /** evaluate line haverage at specific location for horizontal velocity magnitude */
    Real line_hvelmag_average_interpolated(Real x) const;

    /** evaluate line average at specific cell for horizontal velocity magnitude */
    Real line_hvelmag_average_cell(int cell) const;

    /** compute derivatives of horizontal velocity magnitude */
    void compute_line_hvelmag_derivatives();

    /** evaluate line average derivative at specific location for horizontal velocity magnitude */
    Real line_hvelmag_derivative_interpolated(Real x) const;

    /** evaluate derivative of a line average at specific cell horizontal velocity magnitude */
    Real line_hvelmag_derivative_of_average_cell(int ind) const;

};

inline VelPlaneAverage::VelPlaneAverage(
             const Vector<MultiFab*>& mfab_in,
             Geometry geom_in,
             int axis_in,
             bool compute_deriv)
             : m_mfab(mfab_in)
             , m_geom(geom_in)
             , m_axis(axis_in)
             , m_comp_deriv(compute_deriv) {
    AMREX_ALWAYS_ASSERT(m_axis >= 0 && m_axis < AMREX_SPACEDIM);
    m_xlo   = m_geom.ProbLo(m_axis);
    m_dx    = m_geom.CellSize(m_axis);

    Box domain = m_geom.Domain();
    IntVect dom_lo(domain.loVect());
    IntVect dom_hi(domain.hiVect());

    m_ncell_line = dom_hi[m_axis] - dom_lo[m_axis] + 1;
    m_ncell_plane = 1;
    for (int i = 0; i < AMREX_SPACEDIM; ++i) {
       if (i != m_axis) m_ncell_plane *= (dom_hi[i] - dom_lo[i] + 1);
    }

    m_line_xcentroid.resize(m_ncell_line);
    for (int i = 0; i < m_ncell_line; ++i) {
       m_line_xcentroid[i] = m_xlo + (i + 0.5) * m_dx;
    }

    m_line_hvelmag_average.resize(m_ncell_line, 0.0);
    if (m_comp_deriv) m_line_hvelmag_deriv.resize(m_ncell_line, 0.0);
}

inline void VelPlaneAverage::operator()()
{
    std::fill(m_line_hvelmag_average.begin(), m_line_hvelmag_average.end(), 0.0);
    switch (m_axis) {
    case 0:
        compute_hvelmag_averages(XDir(), 1, 2, m_mfab);
        break;
    case 1:
        compute_hvelmag_averages(YDir(), 0, 2, m_mfab);
        break;
    case 2:
        compute_hvelmag_averages(ZDir(), 0, 1, m_mfab);
        break;
    default:
        Abort("axis must be equal to 0, 1, or 2");
        break;
    }
    if (m_comp_deriv) compute_line_hvelmag_derivatives();
}

template <typename IndexSelector>
void VelPlaneAverage::compute_hvelmag_averages(
                      const IndexSelector& idx_op,
                      const int h1_idx,
                      const int h2_idx,
                      const Vector<MultiFab*> mfab)
{
    const Real denom = 1.0 / (Real)m_ncell_plane;
    AsyncArray<Real> lavg(m_line_hvelmag_average.data(), m_line_hvelmag_average.size());
    Real* line_avg = lavg.data();

#ifdef _OPENMP
#pragma omp parallel if (Gpu::notInLaunchRegion())
#endif
    for (MFIter mfi(*(mfab[h1_idx]), TilingIfNotGPU()); mfi.isValid(); ++mfi) {
        Box bx = mfi.tilebox();

        auto xfab_arr = mfab[h1_idx]->const_array(mfi);
        auto yfab_arr = mfab[h2_idx]->const_array(mfi);

        Box pbx = PerpendicularBox<IndexSelector>(bx, IntVect{0, 0, 0});
        ParallelFor(Gpu::KernelInfo().setReduction(true), pbx, [=] AMREX_GPU_DEVICE(
            int p_i, int p_j, int p_k, Gpu::Handler const& handler) noexcept {
            // Loop over the direction perpendicular to the plane.
            // This reduces the atomic pressure on the destination arrays.
            Box lbx = ParallelBox<IndexSelector>(bx, IntVect{p_i, p_j, p_k});

            for (int k = lbx.smallEnd(2); k <= lbx.bigEnd(2); ++k) {
                for (int j = lbx.smallEnd(1); j <= lbx.bigEnd(1); ++j) {
                    for (int i = lbx.smallEnd(0); i <= lbx.bigEnd(0); ++i) {
                        const int ind = idx_op.getIndx(i,j,k);
                        const Real hvelmag = std::sqrt(xfab_arr(i,j,k)*xfab_arr(i,j,k)+yfab_arr(i,j,k)*yfab_arr(i,j,k));
                        Gpu::deviceReduceSum(&line_avg[ind],hvelmag*denom,handler);
                    }
                }
           }
      });
    }
    lavg.copyToHost(m_line_hvelmag_average.data(), m_line_hvelmag_average.size());
    ParallelDescriptor::ReduceRealSum(m_line_hvelmag_average.data(), m_line_hvelmag_average.size());
}

inline Real
VelPlaneAverage::line_hvelmag_average_interpolated(Real x) const
{
    Real c = 0.0;
    int ind = 0;
    if (x > m_xlo + 0.5 * m_dx) {
        ind = static_cast<int>(floor((x - m_xlo) / m_dx - 0.5));
        const Real x1 = m_xlo + (ind + 0.5) * m_dx;
        c = (x - x1) / m_dx;
    }
    if (ind + 1 >= m_ncell_line) {
        ind = m_ncell_line - 2;
        c = 1.0;
    }
    AMREX_ALWAYS_ASSERT(ind >= 0 && ind + 1 < m_ncell_line);
    return m_line_hvelmag_average[ind] * (1.0 - c) +
           m_line_hvelmag_average[ind + 1] * c;
}

inline void VelPlaneAverage::compute_line_hvelmag_derivatives()
{
    for (int i = 0; i < m_ncell_line; ++i) {
        m_line_hvelmag_deriv[i] = line_hvelmag_derivative_of_average_cell(i);
    }
}

inline Real
VelPlaneAverage::line_hvelmag_derivative_of_average_cell(int ind) const
{
    AMREX_ALWAYS_ASSERT(ind >= 0 && ind < m_ncell_line);
    Real dudx;
    if (ind == 0)
        dudx =
            (m_line_hvelmag_average[(ind + 1)] - m_line_hvelmag_average[ind]) /
            m_dx;
    else if (ind == m_ncell_line - 1)
        dudx = (m_line_hvelmag_average[ind] - m_line_hvelmag_average[ind - 1]) /
               m_dx;
    else
        dudx = 0.5 *
               (m_line_hvelmag_average[ind + 1] -
                m_line_hvelmag_average[ind - 1]) /
               m_dx;

    return dudx;
}


#endif /* VelPlaneAverage_H */
