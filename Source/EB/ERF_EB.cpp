#include <AMReX_ParmParse.H>
#include <AMReX_Print.H>
#include <AMReX_WriteEBSurface.H>

#include <ERF_EBIFTerrain.H>
#include <ERF_ProbCommon.H>

#include <ERF_EB.H>
#include <ERF_EBToPVD.H>

#include <AMReX.H>
#include <AMReX_Geometry.H>
#include <AMReX_MultiFab.H>
#include <AMReX_MultiCutFab.H>
#include <AMReX_EBFArrayBox.H>
#include <AMReX_EB2.H>

using namespace amrex;

eb_::~eb_()
{
  if (m_factory) { m_factory.reset(nullptr); }
}

eb_:: eb_ ( Geometry const& a_geom, amrex::FArrayBox const& terrain_fab,
            amrex::Gpu::DeviceVector<amrex::Real>& a_dz_stretched,
            bool /*is_anelastic*/)
    : m_has_eb(0),
      m_support_level(EBSupport::full),
      m_write_eb_surface(0)
{
    m_type = "terrain";

    // int max_coarsening_level;
    // if (is_anelastic) {
    //     max_coarsening_level = 100;
    // } else {
    //     max_coarsening_level = 0;
    // }

    // int max_level_here = 0;

    if (m_type == "terrain")
    {
        Print() << "\nBuilding EB geometry based on terrain.\n";

        TerrainIF ebterrain(terrain_fab, a_geom, a_dz_stretched);

        auto gshop = EB2::makeShop(ebterrain);

        // EB2::Build(gshop, a_geom, max_level_here, max_level_here+max_coarsening_level);
        build_level(a_geom, gshop); // This calls EB2::Build

        m_write_eb_surface = 1;

        make_factory(a_geom, m_eb_level->DistributionMap(), m_eb_level->boxArray());

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
    eb_::WriteEBSurface(a_grids, a_dmap, a_geom, m_factory.get());
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


void
eb_::
WriteEBSurface (const BoxArray & ba, const DistributionMapping & dmap, const Geometry & geom,
                     const EBFArrayBoxFactory * ebf)
{

    EBToPVD eb_to_pvd;

    const Real* dx     = geom.CellSize();
    const Real* problo = geom.ProbLo();

    MultiFab mf_ba(ba, dmap, 1, 0, MFInfo(), *ebf);

    for (MFIter mfi(mf_ba); mfi.isValid(); ++mfi) {

        const auto & sfab    = static_cast<EBFArrayBox const &>(mf_ba[mfi]);
        const auto & my_flag = sfab.getEBCellFlagFab();
        const auto * my_flag_ptr = &my_flag;

        const Box & bx = mfi.validbox();

        if (my_flag.getType(bx) == FabType::covered ||
            my_flag.getType(bx) == FabType::regular) { continue; }

        std::array<const CutFab *, AMREX_SPACEDIM> areafrac;
        const CutFab * bndrycent;

        for (int d = 0; d < AMREX_SPACEDIM; ++d) {
            areafrac[d]  =   &(*ebf->getAreaFrac()[d])[mfi];
        }
        bndrycent = &(ebf->getBndryCent()[mfi]);

#ifdef AMREX_USE_GPU
        std::unique_ptr<EBCellFlagFab> host_flag;
        if (my_flag.arena()->isManaged() || my_flag.arena()->isDevice()) {
            host_flag = std::make_unique<EBCellFlagFab>(my_flag.box(), my_flag.nComp(),
                                                  The_Pinned_Arena());
            Gpu::dtoh_memcpy_async(host_flag->dataPtr(), my_flag.dataPtr(),
                                   host_flag->nBytes());
            Gpu::streamSynchronize();
            my_flag_ptr = host_flag.get();
        }

        std::array<std::unique_ptr<CutFab>, AMREX_SPACEDIM> areafrac_h;
        for (int d = 0; d < AMREX_SPACEDIM; ++d) {
            if (areafrac[d]->arena()->isManaged() || areafrac[d]->arena()->isDevice()) {
                areafrac_h[d] = std::make_unique<CutFab>(areafrac[d]->box(), areafrac[d]->nComp(),
                                                           The_Pinned_Arena());
                Gpu::dtoh_memcpy_async(areafrac_h[d]->dataPtr(), areafrac[d]->dataPtr(),
                                       areafrac[d]->size()*sizeof(Real));
                Gpu::streamSynchronize();
                areafrac[d] = areafrac_h[d].get();
            }
        }

        std::unique_ptr<CutFab> bndrycent_h;
        if (bndrycent->arena()->isManaged() || bndrycent->arena()->isDevice()) {
            bndrycent_h = std::make_unique<CutFab>(bndrycent->box(), bndrycent->nComp(),
                                                        The_Pinned_Arena());
            Gpu::dtoh_memcpy_async(bndrycent_h->dataPtr(), bndrycent->dataPtr(),
                                   bndrycent->size()*sizeof(Real));
            Gpu::streamSynchronize();
            bndrycent = bndrycent_h.get();
        }
#endif

        eb_to_pvd.EBToPolygon(
                problo, dx,
                bx, my_flag_ptr->const_array(),
                bndrycent->const_array(),
                areafrac[0]->const_array(),
                areafrac[1]->const_array(),
                areafrac[2]->const_array());
    }

    int cpu = ParallelDescriptor::MyProc();
    int nProcs = ParallelDescriptor::NProcs();

    eb_to_pvd.WriteEBVTP(cpu);

    if(ParallelDescriptor::IOProcessor()) {
        EBToPVD::WritePVTP(nProcs);
    }

    for (MFIter mfi(mf_ba); mfi.isValid(); ++mfi) {

        const auto & sfab    = static_cast<EBFArrayBox const &>(mf_ba[mfi]);
        const auto & my_flag = sfab.getEBCellFlagFab();
        const auto * my_flag_ptr = &my_flag;

        const Box & bx = mfi.validbox();

        if (my_flag.getType(bx) == FabType::covered ||
            my_flag.getType(bx) == FabType::regular) { continue; }

#ifdef AMREX_USE_GPU
        std::unique_ptr<EBCellFlagFab> host_flag;
        if (my_flag.arena()->isManaged() || my_flag.arena()->isDevice()) {
            host_flag = std::make_unique<EBCellFlagFab>(my_flag.box(), my_flag.nComp(),
                                                  The_Pinned_Arena());
            Gpu::dtoh_memcpy_async(host_flag->dataPtr(), my_flag.dataPtr(),
                                   host_flag->nBytes());
            Gpu::streamSynchronize();
            my_flag_ptr = host_flag.get();
        }
#endif

        eb_to_pvd.EBGridCoverage(cpu, problo, dx, bx, my_flag_ptr->const_array());
    }
}