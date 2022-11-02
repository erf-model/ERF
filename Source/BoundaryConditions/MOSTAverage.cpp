#include <MOSTAverage.H>

// Constructor
MOSTAverage::MOSTAverage (amrex::Vector<const amrex::MultiFab*> fields,
                          amrex::Vector<amrex::MultiFab*> averages,
                          const amrex::MultiFab* z_nd,
                          amrex::Geometry geom,
                          int axis)
    : m_fields(fields), m_averages(averages), m_z_nd(z_nd), m_geom(geom), m_axis(axis)
{
    AMREX_ALWAYS_ASSERT(m_axis >= 0 && m_axis < AMREX_SPACEDIM);
    AMREX_ALWAYS_ASSERT(fields.size() >= 2);

    // Domain bottom and cell-size
    m_xlo   = m_geom.ProbLo  (m_axis);
    m_dx    = m_geom.CellSize(m_axis);

    // Cells per line
    amrex::Box domain = m_geom.Domain();
    amrex::IntVect dom_lo(domain.loVect());
    amrex::IntVect dom_hi(domain.hiVect());
    m_ncell_line = dom_hi[m_axis] - dom_lo[m_axis] + 1;

    // Num components, line avg, cells per plane
    int asize = m_averages.size();   
    m_ncomps.resize( asize );
    m_line_average.resize( asize );
    m_ncell_plane.resize( asize );
    for (int i(0); i<asize; ++i) {
        m_ncomps[i] = m_averages[i]->nComp();
        m_line_average[i].resize(static_cast<size_t>(m_ncell_line) * m_ncomps[i], 0.0);
        
        m_ncell_plane[i] = 1;
        amrex::IndexType ixt = m_averages[i]->boxArray().ixType();
        for (int j = 0; j < AMREX_SPACEDIM; ++j) {
            if (j != m_axis) {
                if (ixt.nodeCentered(j)) {
                    m_ncell_plane[i] *= (dom_hi[j] - dom_lo[j] + 2);
                } else {
                    m_ncell_plane[i] *= (dom_hi[j] - dom_lo[j] + 1);
                }
            }
        }
    }
}

// Write data
void
MOSTAverage::write_most_averages()
{
    int nmf = m_fields.size();
    std::ofstream ofile;
    ofile.open ("MOSTAverages.txt");
    ofile << "Averages compute via the MOSTAverages class:\n";
    for (int k(0); k<m_ncell_line; ++k) {
        ofile << "Index: " << k << ' ';
        for (int imf(0); imf<=nmf; ++imf) {
            int ncomp = m_ncomps[imf];
            for (int n(0); n<ncomp; ++n) {
                ofile << m_line_average[imf][ncomp * k + n] << ' ';
            }
        }
        ofile << "\n";
    }
    ofile.close();
}
