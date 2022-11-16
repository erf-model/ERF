#include <MOSTAverage.H>

// Constructor
MOSTAverage::MOSTAverage (amrex::Real& zref,
                          const amrex::Vector<amrex::Geometry>& geom,
                          const amrex::Vector<amrex::Vector<amrex::MultiFab>>& vars_old,
                          const amrex::Vector<std::unique_ptr<amrex::MultiFab>>& Theta_prim,
                          const amrex::Vector<std::unique_ptr<amrex::MultiFab>>& z_phys_nd)
  : m_zref(zref), m_geom(geom)
{
    // Get PTR to fields and create 2D MFs for averages
    //--------------------------------------------------------
    int nlevs = m_geom.size();
    m_fields.resize(nlevs);
    m_k_indx.resize(nlevs);
    m_averages.resize(nlevs);
    m_z_phys_nd.resize(nlevs);
    for (int lev = 0; lev < nlevs; lev++)
    {
      m_fields[lev].resize(nvar);
      m_averages[lev].resize(navg);
      m_z_phys_nd[lev] = z_phys_nd[lev].get();
      { // Nodal in x
        auto& mf = vars_old[lev][Vars::xvel];
        // Create a 2D ba, dm, & ghost cells
        amrex::BoxArray ba  = mf.boxArray();
        amrex::BoxList bl2d = ba.boxList();
        for (auto& b : bl2d) {
          b.setRange(2,0);
        }
        amrex::BoxArray ba2d(std::move(bl2d));
        const amrex::DistributionMapping& dm = mf.DistributionMap();
        const int ncomp   = 1;
        amrex::IntVect ng = mf.nGrowVect(); ng[2]=0;

          m_fields[lev][0] = &mf;
        m_averages[lev][0] = new amrex::MultiFab(ba2d,dm,ncomp,ng);
        m_averages[lev][0]->setVal(1.E34);
      }
      { // Nodal in y
        auto& mf = vars_old[lev][Vars::yvel];
        // Create a 2D ba, dm, & ghost cells
        amrex::BoxArray ba  = mf.boxArray();
        amrex::BoxList bl2d = ba.boxList();
        for (auto& b : bl2d) {
          b.setRange(2,0);
        }
        amrex::BoxArray ba2d(std::move(bl2d));
        const amrex::DistributionMapping& dm = mf.DistributionMap();
        const int ncomp   = 1;
        amrex::IntVect ng = mf.nGrowVect(); ng[2]=0;

          m_fields[lev][1] = &mf;
        m_averages[lev][1] = new amrex::MultiFab(ba2d,dm,ncomp,ng);
        m_averages[lev][1]->setVal(1.E34);
      }
      { // CC vars
        auto& mf = *Theta_prim[lev];
        // Create a 2D ba, dm, & ghost cells
        amrex::BoxArray ba  = mf.boxArray();
        amrex::BoxList bl2d = ba.boxList();
        for (auto& b : bl2d) {
          b.setRange(2,0);
        }
        amrex::BoxArray ba2d(std::move(bl2d));
        const amrex::DistributionMapping& dm = mf.DistributionMap();
        const int ncomp   = 1;
        const int incomp  = 1;
        amrex::IntVect ng = mf.nGrowVect(); ng[2]=0;

          m_fields[lev][2] = &mf;
        m_averages[lev][2] = new amrex::MultiFab(ba2d,dm,ncomp,ng);
        m_averages[lev][2]->setVal(1.E34);

        m_averages[lev][3] = new amrex::MultiFab(ba2d,dm,ncomp,ng);
        m_averages[lev][3]->setVal(1.E34);

        m_k_indx[lev] = new amrex::iMultiFab(ba2d,dm,incomp,ng);
      }
    } // lev
}



// Driver to call appropriate average member function
void
MOSTAverage::compute_averages()
{
  amrex::Print() << "INSIDE MAC->CA\n";
}
