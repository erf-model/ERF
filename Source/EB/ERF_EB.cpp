#include <AMReX_ParmParse.H>
#include <AMReX_Print.H>
#include <AMReX_WriteEBSurface.H>

#include <ERF_EBIFTerrain.H>
#include <ERF_ProbCommon.H>

#include <ERF_EB.H>

using namespace amrex;

eb_::~eb_()
{
  if (m_factory) { m_factory.reset(nullptr); }
}

eb_:: eb_ ( Geometry const& a_geom, amrex::FArrayBox const& terrain_fab,
            amrex::Gpu::DeviceVector<amrex::Real>& a_dz_stretched,
            bool is_anelastic)
    : m_has_eb(0),
      m_support_level(EBSupport::full),
      m_write_eb_surface(0)
{
    m_type = "terrain";

    int max_coarsening_level;
    if (is_anelastic) {
        max_coarsening_level = 100;
    } else {
        max_coarsening_level = 0;
    }

    int max_level_here = 0;

    if (m_type == "terrain")
    {
        Print() << "\nBuilding EB geometry based on terrain.\n";

        TerrainIF ebterrain(terrain_fab, a_geom, a_dz_stretched);

        auto gshop = EB2::makeShop(ebterrain);

        EB2::Build(gshop, a_geom, max_level_here, max_level_here+max_coarsening_level);

        m_has_eb = 1;

    } else if (m_type == "box") {

        Print() << "\nBuilding box geometry.\n";
        make_box(a_geom);
        m_has_eb = 1;

    } else {

        EB2::AllRegularIF regular;
        auto gshop = EB2::makeShop(regular);
        build_level(a_geom, gshop);
    }
}

void
eb_::
make_factory ( Geometry            const& a_geom,
               DistributionMapping const& a_dmap,
               BoxArray            const& a_grids)
{
  Print() << "making EB factory\n";
  m_factory = std::make_unique<EBFArrayBoxFactory>(*m_eb_level, a_geom, a_grids, a_dmap,
    Vector<int>{nghost_basic(), nghost_volume(), nghost_full()}, m_support_level);

  if (m_write_eb_surface) {
    WriteEBSurface(a_grids, a_dmap, a_geom, m_factory.get());
  }

  { int const idim(0);

    Print() << "making EB staggered u-factory\n";
    //m_u_factory.set_verbose();
    m_u_factory.define(idim, a_geom, a_grids, a_dmap,
      Vector<int>{nghost_basic(), nghost_volume(), nghost_full()},
      m_factory.get());
  }

  { int const idim(1);
    Print() << "making EB staggered v-factory\n";
    //m_v_factory.set_verbose();
    m_v_factory.define(idim, a_geom, a_grids, a_dmap,
      Vector<int>{nghost_basic(), nghost_volume(), nghost_full()},
      m_factory.get());
  }

  { int const idim(2);
    Print() << "making EB staggered w-factory\n";
    //m_w_factory.set_verbose();
    m_w_factory.define(idim, a_geom, a_grids, a_dmap,
      Vector<int>{nghost_basic(), nghost_volume(), nghost_full()},
      m_factory.get());
  }

  Print() << "\nDone making EB factory.\n\n";
}
