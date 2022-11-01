#ifndef MOSTAverage_H
#define MOSTAverage_H

#include <AMReX_Gpu.H>
#include <AMReX_FArrayBox.H>
#include <AMReX_MultiFab.H>
#include <DirectionSelector.H>

class MOSTAverage {
public:
    explicit MOSTAverage(amrex::Vector<const amrex::MultiFab*> fields,
                         amrex::Vector<amrex::MultiFab*> averages,
                         const amrex::MultiFab* z_nd,
                         const amrex::Geometry geom,
                         int axis);
     MOSTAverage() = default;
    ~MOSTAverage() = default;

    // Evaluate line average at specific location for any average component
    void line_average_interpolated(amrex::Real x, int icomp, int imf) const;

    // Set local vector to m_line_average
    void line_average(int icomp, int imf, amrex::Vector<amrex::Real>& l_vec);

    // Fill line storage with averages
    template <typename IndexSelector>
    void compute_plane_averages(const IndexSelector& idxOp);

    // Get number of components for a mf in m_fields
    int ncomp()      const  { return m_ncomps[0]; }
    int ncomp(int i) const  { return m_ncomps[i]; }

    // Get pointer to mf in m_fields
    const amrex::MultiFab* get_field()      const { return m_fields[0]; }
    const amrex::MultiFab* get_field(int i) const { return m_fields[i]; }

    // Get pointer to the 2D mf of averages
    const amrex::MultiFab* get_policy_average()      const { return m_averages[0]; };
    const amrex::MultiFab* get_policy_average(int i) const { return m_averages[i]; };

    // Get plane averages long dir()
    const amrex::Vector<amrex::Real>& get_line_average()      const { return m_line_average[0]; }
    const amrex::Vector<amrex::Real>& get_line_average(int i) const { return m_line_average[i]; }

    // Set precision of text file output
    void set_precision(int p) { m_precision = p; }

    // Write text file with m_line_averages
    void write_most_averages();

    // Get geometry values
    amrex::Real dx()  const { return m_dx; }
    amrex::Real xlo() const { return m_xlo; }
    int axis()        const { return m_axis; }
    int level()       const { return m_level; }
    int ncell_plane() const { return m_ncell_plane; }
    int ncell_line()  const { return m_ncell_line; }

protected:
    int m_ncell_plane; // Number of cells in plane
    int m_ncell_line;  // Number of cells along line

    amrex::Real m_dx;  // Mesh spacing in axis direction
    amrex::Real m_xlo; // Bottom of domain in axis direction

    int m_precision   = 4; // Precision for line plot text file
    const int m_level = 0; // Level for plane averaging for now fixed at level=0

    amrex::Vector<int> m_ncomps; // Number of components

    amrex::Vector<amrex::Vector<amrex::Real>> m_line_average; // Planar average along dir()

    amrex::Vector<const amrex::MultiFab*> m_fields;  // Ptr to fields to be averaged
    amrex::Vector<amrex::MultiFab*> m_averages;      // Ptr to 2D mf to hold averages
    const amrex::MultiFab* m_z_nd;                   // Ptr to terrain heigh coords
    amrex::Geometry m_geom;                          // Geometry
    const int m_axis;                                // Averaging axis
};

// Evaluate line average at specific location for any average component
AMREX_FORCE_INLINE
void
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

    // Interpolated value between planar averages
    amrex::Real avg = m_line_average[imf][ncomp * ind       + icomp] * (1.0 - c) +
                      m_line_average[imf][ncomp * (ind + 1) + icomp] * c;

    // Set the 2D multifab to a constant value
    m_averages[imf]->setVal(avg);
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
MOSTAverage::compute_plane_averages(const IndexSelector& idxOp)
{
    const amrex::Real denom = 1.0 / (amrex::Real)m_ncell_plane;

    // Averages over all the fields
    //----------------------------------------------------------
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

    // Averages for the tangential velocity magnitude
    //----------------------------------------------------------
    int imf  = 0;
    int iavg = m_fields.size();
    const int ncomp = 1;

    amrex::AsyncArray<amrex::Real> lavg(m_line_average[iavg].data(), m_line_average[iavg].size());
    amrex::Real* line_avg = lavg.data();

#ifdef _OPENMP
#pragma omp parallel if (amrex::Gpu::notInLaunchRegion())
#endif
    for (amrex::MFIter mfi(*m_fields[imf], amrex::TilingIfNotGPU()); mfi.isValid(); ++mfi) {
        amrex::Box bx  = amrex::enclosedCells(mfi.tilebox());
        amrex::Box pbx = PerpendicularBox<IndexSelector>(bx, amrex::IntVect{0, 0, 0});

        auto u_mf_arr = m_fields[imf]->const_array(mfi);
        auto v_mf_arr = m_fields[imf+1]->const_array(mfi);

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
                            const amrex::Real u_val = 0.5 * (u_mf_arr(i,j,k) + u_mf_arr(i+1,j  ,k));
                            const amrex::Real v_val = 0.5 * (v_mf_arr(i,j,k) + v_mf_arr(i  ,j+1,k));
                            const amrex::Real mag   = std::sqrt(u_val*u_val + v_val*v_val);
                            amrex::Gpu::deviceReduceSum(&line_avg[ncomp * ind + n],
                                                        mag * denom, handler);
                        }
                    }
                }
            }
        });
    }

    lavg.copyToHost(m_line_average[iavg].data(), m_line_average[iavg].size());
    amrex::ParallelDescriptor::ReduceRealSum(m_line_average[iavg].data(), m_line_average[iavg].size());
}
#endif
