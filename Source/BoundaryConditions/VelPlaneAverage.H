
#ifndef VelPlaneAverage_H
#define VelPlaneAverage_H

#include "AMReX_Gpu.H"
#include "AMReX_MultiFab.H"
#include "DirectionSelector.H"

class VelPlaneAverage {
public:
    VelPlaneAverage(amrex::Geometry geom_in);
    ~VelPlaneAverage() = default;
private:

    amrex::Vector<amrex::Real> m_line_xcentroid; /** line storage for centroids of each cell along a line*/
    amrex::Real m_dz;  /** mesh spacing in axis direction*/
    amrex::Real m_xlo; /** bottom of domain in axis direction */
    int m_ncell_plane; /** number of cells in plane */
    int m_ncell_line;  /** number of cells along line */
    amrex::Geometry m_geom;

    amrex::Vector<amrex::Real> m_line_hvelmag_average; /** line storage for the average horizontal velocity magnitude */

public:
    void compute_hvelmag_averages(amrex::MultiFab& u_mf,
                                  amrex::MultiFab& v_mf);

    /** return vector containing horizontal velocity magnitude average */
    const amrex::Vector<amrex::Real>& line_hvelmag_average()
    {
        return m_line_hvelmag_average;
    }

    /** evaluate line haverage at specific location for horizontal velocity magnitude */
    amrex::Real line_hvelmag_average_interpolated(amrex::Real x) const;
};

VelPlaneAverage::VelPlaneAverage(amrex::Geometry geom_in)
     : m_geom(geom_in)
{
    m_xlo   = m_geom.ProbLo(2);
    m_dz    = m_geom.CellSize(2);

    amrex::Box domain = m_geom.Domain();
    amrex::IntVect dom_lo(domain.loVect());
    amrex::IntVect dom_hi(domain.hiVect());

    m_ncell_line  =  dom_hi[2] - dom_lo[2] + 1;
    m_ncell_plane = (dom_hi[0] - dom_lo[0] + 1) * (dom_hi[1] - dom_lo[1] + 1);

    m_line_xcentroid.resize(m_ncell_line);
    for (int i = 0; i < m_ncell_line; ++i) {
       m_line_xcentroid[i] = m_xlo + (i + 0.5) * m_dz;
    }

    m_line_hvelmag_average.resize(m_ncell_line, 0.0);
}

void
VelPlaneAverage::compute_hvelmag_averages(amrex::MultiFab& u_mf,
                                          amrex::MultiFab& v_mf)
{
    const amrex::Real denom = 1.0 / (amrex::Real)m_ncell_plane;

    amrex::Box domain = m_geom.Domain();
    const int dom_lo_z = domain.smallEnd(2);
    const int dom_hi_z = domain.bigEnd(2);

    for (int k = dom_lo_z; k < dom_hi_z; k++) {
        m_line_hvelmag_average[k] = 0.;
    }
    amrex::AsyncArray<amrex::Real> lavg(m_line_hvelmag_average.data(), m_line_hvelmag_average.size());
    amrex::Real* line_avg = lavg.data();

#ifdef _OPENMP
#pragma omp parallel if (amrex::Gpu::notInLaunchRegion())
#endif
    for (amrex::MFIter mfi(u_mf, amrex::TilingIfNotGPU()); mfi.isValid(); ++mfi)
    {
        amrex::Box bx = amrex::enclosedCells(mfi.tilebox());

        auto u_arr = u_mf.const_array(mfi);
        auto v_arr = v_mf.const_array(mfi);

        bx.setSmall(2, dom_lo_z);
        bx.setBig(2, dom_hi_z-1); // the -1 is because in the k loop we start from k = dom_lo_z
                                  // and continue while k < dom_hi_z

        ParallelFor(bx, [=] AMREX_GPU_DEVICE (int i, int j, int k) {
            const amrex::Real xvel = 0.5 * (u_arr(i,j,k) + u_arr(i+1,j,k));
            const amrex::Real yvel = 0.5 * (v_arr(i,j,k) + v_arr(i,j+1,k));
            const amrex::Real hvelmag = std::sqrt(xvel*xvel + yvel*yvel);
            amrex::Gpu::Atomic::AddNoRet(&line_avg[k], hvelmag);
        });
    }


    lavg.copyToHost(m_line_hvelmag_average.data(), m_line_hvelmag_average.size());
    for (int k = dom_lo_z; k < dom_hi_z; k++) {
        m_line_hvelmag_average[k] *= denom;
    }
    amrex::ParallelDescriptor::ReduceRealSum(m_line_hvelmag_average.data(), m_line_hvelmag_average.size());
}

AMREX_FORCE_INLINE
amrex::Real
VelPlaneAverage::line_hvelmag_average_interpolated(amrex::Real x) const
{
    amrex::Real c = 0.0;
    int ind = 0;
    if (x > m_xlo + 0.5 * m_dz) {
        ind = static_cast<int>(floor((x - m_xlo) / m_dz - 0.5));
        const amrex::Real x1 = m_xlo + (ind + 0.5) * m_dz;
        c = (x - x1) / m_dz;
    }
    if (ind + 1 >= m_ncell_line) {
        ind = m_ncell_line - 2;
        c = 1.0;
    }
    AMREX_ALWAYS_ASSERT(ind >= 0 && ind + 1 < m_ncell_line);
    return m_line_hvelmag_average[ind] * (1.0 - c) +
           m_line_hvelmag_average[ind + 1] * c;
}

#endif /* VelPlaneAverage_H */
