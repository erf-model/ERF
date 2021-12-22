
#ifndef PlaneAverage_H
#define PlaneAverage_H

#include "ERF.H"
#include "DirectionSelector.H"

using namespace amrex;

class PlaneAverage {
public:
    explicit PlaneAverage(const MultiFab* field_in,
                          const Geometry geom_in,
                          int axis_in,
                          bool compute_deriv = false);
    PlaneAverage() = default;
    ~PlaneAverage() = default;

    void operator()();

    /** evaluate line average at specific location for any average component */
    Real line_average_interpolated(Real x, int comp) const;
    /** evaluate line average at specific cell for any average component */
    Real line_average_cell(int cell, int comp) const;

    /** evaluate line average derivative at specific location for any average component */
    Real line_derivative_interpolated(Real x, int comp) const;
    /** evaluate derivative of a line average at specific cell for any component
     */
    Real line_derivative_of_average_cell(int ind, int comp) const;

    /** change precision of text file output */
    void set_precision(int p) { m_precision = p; };

    Real dx() const { return m_dx; };
    Real xlo() const { return m_xlo; };

    int axis() const { return m_axis; };
    int level() const { return m_level; };
    int ncomp() const { return m_ncomp; };
    int ncell_plane() const { return m_ncell_plane; };
    int ncell_line() const { return m_ncell_line; };

    const Vector<Real>& line_average() const
    {
        return m_line_average;
    };
    void line_average(int comp, Vector<Real>& l_vec);
    const Vector<Real>& line_centroids() const
    {
        return m_line_xcentroid;
    };

    const MultiFab& field() const { return *m_field; };

protected:
    int m_ncomp; /** number of average components */

    Vector<Real> m_line_average; /** line storage for the average velocity and tracer variables */

    //! line storage for the derivative of average field
    Vector<Real> m_line_deriv;

    Vector<Real> m_line_xcentroid; /** line storage for centroids of each cell along a line*/

    Real m_dx;  /** mesh spacing in axis direction*/
    Real m_xlo; /** bottom of domain in axis direction */

    int m_ncell_plane; /** number of cells in plane */
    int m_ncell_line;  /** number of cells along line */

    int m_precision = 4; /** precision for line plot text file */
    const int m_level = 0; /** level for plane averaging for now fixed at level=0 */

    const MultiFab* m_field;
    Geometry m_geom;
    const int m_axis;
    const bool m_comp_deriv;

public:
    /** fill line storage with averages */
    template <typename IndexSelector>
    void compute_averages(const IndexSelector& idxOp, const MultiFab& mfab);

    /** fill derivatives of line storage */
    void compute_line_derivatives();
};

inline PlaneAverage::PlaneAverage(
             const MultiFab* field_in,
             Geometry geom_in,
             int axis_in,
             bool compute_deriv)
             : m_field(field_in)
             , m_geom(geom_in)
             , m_axis(axis_in)
             , m_comp_deriv(compute_deriv) {

    AMREX_ALWAYS_ASSERT(m_axis >= 0 && m_axis < AMREX_SPACEDIM);

    m_xlo = m_geom.ProbLo(m_axis);
    m_dx  = m_geom.CellSize(m_axis);
    m_ncomp = m_field->nComp();

    Box domain = m_geom.Domain();
    IntVect dom_lo(domain.loVect());
    IntVect dom_hi(domain.hiVect());

    m_ncell_line = dom_hi[m_axis] - dom_lo[m_axis] + 1;

    m_ncell_plane = 1;
    for (int i = 0; i < AMREX_SPACEDIM; ++i) {
       if (i != m_axis) m_ncell_plane *= (dom_hi[i] - dom_lo[i] + 1);
    }

    m_line_average.resize(static_cast<size_t>(m_ncell_line) * m_ncomp, 0.0);
    if (m_comp_deriv)
       m_line_deriv.resize(static_cast<size_t>(m_ncell_line) * m_ncomp, 0.0);
    m_line_xcentroid.resize(m_ncell_line);

    for (int i = 0; i < m_ncell_line; ++i) {
       m_line_xcentroid[i] = m_xlo + (i + 0.5) * m_dx;
    }
}

inline Real
PlaneAverage::line_average_interpolated(Real x, int comp) const
{
    AMREX_ALWAYS_ASSERT(comp >= 0 && comp < m_ncomp);

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

    return m_line_average[m_ncomp * ind + comp] * (1.0 - c) +
           m_line_average[m_ncomp * (ind + 1) + comp] * c;
}

inline void PlaneAverage::line_average(int comp, Vector<Real>& l_vec)
{
    AMREX_ALWAYS_ASSERT(comp >= 0 && comp < m_ncomp);

    for (int i = 0; i < m_ncell_line; i++)
        l_vec[i] = m_line_average[m_ncomp * i + comp];
}

inline Real PlaneAverage::line_average_cell(int ind, int comp) const
{
    AMREX_ALWAYS_ASSERT(comp >= 0 && comp < m_ncomp);
    AMREX_ALWAYS_ASSERT(ind >= 0 && ind < m_ncell_line);

    return m_line_average[m_ncomp * ind + comp];
}

inline void PlaneAverage::compute_line_derivatives()
{
    for (int i = 0; i < m_ncell_line; ++i) {
        for (int n = 0; n < m_ncomp; ++n)
            m_line_deriv[m_ncomp * i + n] =
                line_derivative_of_average_cell(i, n);
    }
}

inline Real PlaneAverage::line_derivative_interpolated(
            Real x, int comp) const
{
    AMREX_ALWAYS_ASSERT(comp >= 0 && comp < m_ncomp);

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

    return m_line_deriv[m_ncomp * ind + comp] * (1.0 - c) +
           m_line_deriv[m_ncomp * (ind + 1) + comp] * c;
}

inline void PlaneAverage::operator()()
{
    std::fill(m_line_average.begin(), m_line_average.end(), 0.0);
    switch (m_axis) {
    case 0:
        compute_averages(XDir(), *m_field);
        break;
    case 1:
        compute_averages(YDir(), *m_field);
        break;
    case 2:
        compute_averages(ZDir(), *m_field);
        break;
    default:
        Abort("axis must be equal to 0, 1, or 2");
        break;
    }
    if (m_comp_deriv) compute_line_derivatives();
}

template <typename IndexSelector>
inline void PlaneAverage::compute_averages(
            const IndexSelector& idxOp, const MultiFab& mfab)
{
    const Real denom = 1.0 / (Real)m_ncell_plane;
    AsyncArray<Real> lavg(m_line_average.data(), m_line_average.size());
    Real* line_avg = lavg.data();
    const int ncomp = m_ncomp;

#ifdef _OPENMP
#pragma omp parallel if (Gpu::notInLaunchRegion())
#endif
    for (MFIter mfi(mfab, TilingIfNotGPU()); mfi.isValid(); ++mfi) {
        Box bx = mfi.tilebox();
        auto fab_arr = mfab.const_array(mfi);
        Box pbx = PerpendicularBox<IndexSelector>(bx, IntVect{0, 0, 0});

        ParallelFor(Gpu::KernelInfo().setReduction(true), pbx, [=] AMREX_GPU_DEVICE( int p_i, int p_j, int p_k,
                   Gpu::Handler const& handler) noexcept {
            // Loop over the direction perpendicular to the plane.
            // This reduces the atomic pressure on the destination arrays.

            Box lbx = ParallelBox<IndexSelector>(bx, IntVect{p_i, p_j, p_k});

            for (int k = lbx.smallEnd(2); k <= lbx.bigEnd(2); ++k) {
                for (int j = lbx.smallEnd(1); j <= lbx.bigEnd(1); ++j) {
                    for (int i = lbx.smallEnd(0); i <= lbx.bigEnd(0); ++i) {
                        int ind = idxOp.getIndx(i, j, k);
                        for (int n = 0; n < ncomp; ++n) {
                            Gpu::deviceReduceSum(&line_avg[ncomp * ind + n], fab_arr(i, j, k, n) * denom, handler);
                        }
                     }
                 }
             }
       });
    }

    lavg.copyToHost(m_line_average.data(), m_line_average.size());
    ParallelDescriptor::ReduceRealSum(m_line_average.data(), m_line_average.size());
}
#endif /* PlaneAverage_H */
