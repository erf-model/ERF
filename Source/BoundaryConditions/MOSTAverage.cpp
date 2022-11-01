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

    // Cells per plane
    m_ncell_plane = 1;
    for (int i = 0; i < AMREX_SPACEDIM; ++i) {
       if (i != m_axis) m_ncell_plane *= (dom_hi[i] - dom_lo[i] + 1);
    }

    // Line averages and number of components
    int fsize = fields.size();
    m_ncomps.resize( fsize + 1 );
    m_line_average.resize( fsize + 1 );
    for (int i(0); i<fsize; ++i) {
        m_ncomps[i] = m_fields[i]->nComp();
        m_line_average[i].resize(static_cast<size_t>(m_ncell_line) * m_ncomps[i], 0.0);
    }
    m_ncomps[fsize] = m_fields[0]->nComp();
    m_line_average[fsize].resize(static_cast<size_t>(m_ncell_line) * m_ncomps[0], 0.0);
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
