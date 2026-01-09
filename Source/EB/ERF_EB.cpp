#include <AMReX_ParmParse.H>
#include <AMReX_Print.H>
#include <AMReX_WriteEBSurface.H>

#include <ERF_EBIFTerrain.H>
#include <ERF_ProbCommon.H>

#include <ERF.H>
#include <ERF_EB.H>

#include <AMReX.H>
#include <AMReX_Geometry.H>
#include <AMReX_MultiFab.H>
#include <AMReX_MultiCutFab.H>
#include <AMReX_EBFArrayBox.H>
#include <AMReX_EB2.H>
#include <AMReX_EBToPVD.H>

using namespace amrex;

eb_::~eb_()
{
    // if (m_factory) { m_factory.reset(nullptr); }
}

eb_::eb_ ( )
    : m_has_eb(0),
      m_support_level(EBSupport::full),
      m_write_eb_surface(0)
{ }

void
eb_::make_all_factories ([[maybe_unused]] int level,
                        Geometry            const& a_geom,
                        BoxArray            const& ba,
                        DistributionMapping const& dm,
                        EB2::Level const& a_eb_level)
{
    Print() << "making EB factory\n";
    m_factory = std::make_unique<EBFArrayBoxFactory>(a_eb_level, a_geom, ba, dm,
        Vector<int>{nghost_basic(), nghost_volume(), nghost_full()}, m_support_level);

    // Correct cell connectivity
    eb_::set_connection_flags();

    { int const idim(0);
        Print() << "making EB staggered u-factory\n";
        //m_u_factory.set_verbose();
        m_u_factory.define(level, idim, a_geom, ba, dm,
            Vector<int>{nghost_basic(), nghost_volume(), nghost_full()},
            m_factory.get());
    }

    { int const idim(1);
        Print() << "making EB staggered v-factory\n";
        //m_v_factory.set_verbose();
        m_v_factory.define(level, idim, a_geom, ba, dm,
            Vector<int>{nghost_basic(), nghost_volume(), nghost_full()},
            m_factory.get());
    }

    { int const idim(2);
        Print() << "making EB staggered w-factory\n";
        //m_w_factory.set_verbose();
        m_w_factory.define(level, idim, a_geom, ba, dm,
            Vector<int>{nghost_basic(), nghost_volume(), nghost_full()},
            m_factory.get());
    }
    Print() << "\nDone making EB factory at level = " << level << ".\n\n";
}

void
eb_::make_cc_factory ([[maybe_unused]] int level,
                        Geometry            const& a_geom,
                        BoxArray            const& ba,
                        DistributionMapping const& dm,
                        EB2::Level const& a_eb_level)
{
    Print() << "making EB factory\n";
    m_factory = std::make_unique<EBFArrayBoxFactory>(a_eb_level, a_geom, ba, dm,
        Vector<int>{nghost_basic(), nghost_volume(), nghost_full()}, m_support_level);

    Print() << "\nDone making EB factory at level " << level << ".\n\n";
}

/*
Reset cell flags to disconnect cells with zero volume fraction,
via non-const reference from EBFArrayBoxFactory.
*/
void
eb_::set_connection_flags ()
{
    // Get non-const reference to EBCellFlagFab FabArray
    FabArray<EBCellFlagFab>& cellflag = getNonConstEBCellFlags(*m_factory);

    const MultiFab& volfrac = m_factory->getVolFrac();

    for (MFIter mfi(cellflag, false); mfi.isValid(); ++mfi) {
        const Box& bx = mfi.validbox();
        const Box gbx = amrex::grow(bx, cellflag.nGrow()-1); // Leave one cell layer

        Array4<EBCellFlag> const& flag = cellflag.array(mfi);
        Array4<Real const> const& vfrac = volfrac.const_array(mfi);

        ParallelFor(gbx, [=] AMREX_GPU_DEVICE (int i, int j, int k) noexcept
        {
            for(int kk(-1); kk<=1; kk++) {
            for(int jj(-1); jj<=1; jj++) {
            for(int ii(-1); ii<=1; ii++)
            {
                if (vfrac(i+ii,j+jj,k+kk) == 0.0) {
                    flag(i,j,k).setDisconnected(ii,jj,kk);
                }
            }}}
        });

        ParallelFor(gbx, [=] AMREX_GPU_DEVICE (int i, int j, int k) noexcept
        {
            if (vfrac(i,j,k)==0.0) {
                flag(i,j,k).setCovered();
            }
        });

    }
}
