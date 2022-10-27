#ifndef MOSTAverage_H
#define MOSTAverage_H

#include <AMReX_Gpu.H>
#include <AMReX_MultiFab.H>
#include <DirectionSelector.H>

class MOSTAverage {
public:
    explicit MOSTAverage(amrex::Vector<const amrex::MultiFab*> fields,
                         const amrex::MultiFab* z_nd,
                         const amrex::Geometry geom,
                         int axis);
     MOSTAverage() = default;
    ~MOSTAverage() = default;
  
    // Evaluate line average at specific location for any average component
    amrex::Real line_average_interpolated(amrex::Real x, int icomp, int imf) const;

    // Set local vector to m_line_average
    void line_average(int icomp, int imf, amrex::Vector<amrex::Real>& l_vec);

    // Fill line storage with averages
    template <typename IndexSelector>
    void compute_averages(const IndexSelector& idxOp);

    // Set precision of text file output
    void set_precision(int p) { m_precision = p; }

    // Write ASCII file with averages
    void write_most_averages();

    amrex::Real dx()  const { return m_dx; }
    amrex::Real xlo() const { return m_xlo; }
    int axis()        const { return m_axis; }
    int level()       const { return m_level; }    
    int ncell_plane() const { return m_ncell_plane; }
    int ncell_line()  const { return m_ncell_line; }

    int ncomp()      const  { return m_ncomps[0]; }
    int ncomp(int i) const  { return m_ncomps[i]; }
    
    const amrex::MultiFab& field()      const { return *(m_fields[0]); }
    const amrex::MultiFab& field(int i) const { return *(m_fields[i]); }
    
    const amrex::Vector<amrex::Real>& line_average()      const { return m_line_average[0]; }
    const amrex::Vector<amrex::Real>& line_average(int i) const { return m_line_average[i]; }
    
    const amrex::Vector<amrex::Real>& line_centroids() const { return m_line_xcentroid;  }

protected:
    int m_ncell_plane; // number of cells in plane
    int m_ncell_line;  // number of cells along line

    amrex::Real m_dx;  // mesh spacing in axis direction
    amrex::Real m_xlo; // bottom of domain in axis direction

    int m_precision   = 4; // precision for line plot text file
    const int m_level = 0; // level for plane averaging for now fixed at level=0
    
    amrex::Vector<int> m_ncomps; // number of components

    amrex::Vector<amrex::Real> m_line_xcentroid; 
    amrex::Vector<amrex::Vector<amrex::Real>> m_line_average;

    amrex::Vector<const amrex::MultiFab*> m_fields;
    const amrex::MultiFab* m_z_nd;
    amrex::Geometry m_geom;
    const int m_axis;
};

// Evaluate line average at specific location for any average component
AMREX_FORCE_INLINE
amrex::Real
MOSTAverage::line_average_interpolated(amrex::Real x, int icomp, int imf) const
{
    int ncomp = m_ncomps[imf];
    AMREX_ALWAYS_ASSERT(icomp >= 0 && icomp < ncomp);
    
    amrex::Real c = 0.0;
    int ind = 0;

    if (x > m_xlo + 0.5 * m_dx) {
        ind = static_cast<int>(floor((x - m_xlo) / m_dx - 0.5));
        const amrex::Real x1 = m_xlo + (ind + 0.5) * m_dx;
        c = (x - x1) / m_dx;
    }

    if (ind + 1 >= m_ncell_line) {
        ind = m_ncell_line - 2;
        c = 1.0;
    }

    AMREX_ALWAYS_ASSERT(ind >= 0 && ind + 1 < m_ncell_line);

    return m_line_average[imf][ncomp * ind + icomp] * (1.0 - c) +
           m_line_average[imf][ncomp * (ind + 1) + icomp] * c;
}

// Set local vector to m_line_average
AMREX_FORCE_INLINE
void
MOSTAverage::line_average(int icomp, int imf, amrex::Vector<amrex::Real>& l_vec)
{
    int ncomp = m_ncomps[imf];
    AMREX_ALWAYS_ASSERT(icomp >= 0 && icomp < ncomp);

    for (int i = 0; i < m_ncell_line; i++)
        l_vec[i] = m_line_average[imf][ncomp * i + icomp];
}

// Fill line storage with averages
template <typename IndexSelector>
void
MOSTAverage::compute_averages(const IndexSelector& idxOp)
{
    const amrex::Real denom = 1.0 / (amrex::Real)m_ncell_plane;

    for (int imf(0); imf<m_fields.size(); ++imf) {
        const int ncomp = m_ncomps[imf];
        
        amrex::AsyncArray<amrex::Real> lavg(m_line_average[imf].data(), m_line_average[imf].size());
        amrex::Real* line_avg = lavg.data();
       
#ifdef _OPENMP
#pragma omp parallel if (amrex::Gpu::notInLaunchRegion())
#endif
        for (amrex::MFIter mfi(*m_fields[imf], amrex::TilingIfNotGPU()); mfi.isValid(); ++mfi) {
            amrex::Box bx  = mfi.tilebox();
            amrex::Box pbx = PerpendicularBox<IndexSelector>(bx, amrex::IntVect{0, 0, 0});

            auto mf_arr = m_fields[imf]->const_array(mfi);
            
            ParallelFor(amrex::Gpu::KernelInfo().setReduction(true), pbx, [=]
            AMREX_GPU_DEVICE(int p_i, int p_j, int p_k, amrex::Gpu::Handler const& handler) noexcept
            {
                // Loop over the direction perpendicular to the plane.
                // This reduces the atomic pressure on the destination arrays.
                amrex::Box lbx = ParallelBox<IndexSelector>(bx, amrex::IntVect{p_i, p_j, p_k});
                
                for (int k = lbx.smallEnd(2); k <= lbx.bigEnd(2); ++k) {
                    for (int j = lbx.smallEnd(1); j <= lbx.bigEnd(1); ++j) {
                        for (int i = lbx.smallEnd(0); i <= lbx.bigEnd(0); ++i) {
                            int ind = idxOp.getIndx(i, j, k);
                            for (int n = 0; n < ncomp; ++n) {
                                amrex::Gpu::deviceReduceSum(&line_avg[ncomp * ind + n],
                                                            mf_arr(i, j, k, n) * denom, handler);
                            }
                        }
                    }
                }
            });
        }

        lavg.copyToHost(m_line_average[imf].data(), m_line_average[imf].size());
        amrex::ParallelDescriptor::ReduceRealSum(m_line_average[imf].data(), m_line_average[imf].size());
    }
}
#endif
