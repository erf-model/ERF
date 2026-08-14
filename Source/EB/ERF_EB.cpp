/**
 * \file ERF_EB.cpp
 * \brief Implements EB factory construction and EB cell-flag connectivity fixes.
 */
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
    eb_::set_connection_flags(m_factory.get());

#if USE_FC_FACTORY
    // New: Native AMReX FC factories (EB2::BuildFC called in ERF.cpp)
    { int const idim(0);
        Print() << "making EB staggered u-factory\n";
        m_u_factory_fc = std::make_unique<EBFArrayBoxFactory>(
            a_eb_level, a_geom, ba, dm,
            Vector<int>{nghost_basic(), nghost_volume(), nghost_full()},
            m_support_level, idim);
        eb_::mask_fc_from_cc(m_u_factory_fc.get(), m_factory.get(), idim);
        eb_::set_connection_flags(m_u_factory_fc.get());
    }

    { int const idim(1);
        Print() << "making EB staggered v-factory\n";
        m_v_factory_fc = std::make_unique<EBFArrayBoxFactory>(
            a_eb_level, a_geom, ba, dm,
            Vector<int>{nghost_basic(), nghost_volume(), nghost_full()},
            m_support_level, idim);
        eb_::mask_fc_from_cc(m_v_factory_fc.get(), m_factory.get(), idim);
        eb_::set_connection_flags(m_v_factory_fc.get());
    }

    { int const idim(2);
        Print() << "making EB staggered w-factory\n";
        m_w_factory_fc = std::make_unique<EBFArrayBoxFactory>(
            a_eb_level, a_geom, ba, dm,
            Vector<int>{nghost_basic(), nghost_volume(), nghost_full()},
            m_support_level, idim);
        eb_::mask_fc_from_cc(m_w_factory_fc.get(), m_factory.get(), idim);
        eb_::set_connection_flags(m_w_factory_fc.get());
    }
#else
    // Original: eb_aux_ factories
    { int const idim(0);
        Print() << "making EB staggered u-factory\n";
        m_u_factory.define(level, idim, a_geom, ba, dm,
            Vector<int>{nghost_basic(), nghost_volume(), nghost_full()},
            m_factory.get());
    }

    { int const idim(1);
        Print() << "making EB staggered v-factory\n";
        m_v_factory.define(level, idim, a_geom, ba, dm,
            Vector<int>{nghost_basic(), nghost_volume(), nghost_full()},
            m_factory.get());
    }

    { int const idim(2);
        Print() << "making EB staggered w-factory\n";
        m_w_factory.define(level, idim, a_geom, ba, dm,
            Vector<int>{nghost_basic(), nghost_volume(), nghost_full()},
            m_factory.get());
    }
#endif
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

/**
 * \brief Reset cell flags to disconnect cells with zero volume fraction.
 *
 * The factory EBCellFlagFab data are updated through a non-const reference.
 */
void
eb_::set_connection_flags (EBFArrayBoxFactory* factory)
{
    // Get non-const reference to EBCellFlagFab FabArray
    FabArray<EBCellFlagFab>& cellflag = getNonConstEBCellFlags(*factory);

    const MultiFab& volfrac = factory->getVolFrac();

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
                if (vfrac(i+ii,j+jj,k+kk) == zero) {
                    flag(i,j,k).setDisconnected(ii,jj,kk);
                }
            }}}
        });

        ParallelFor(gbx, [=] AMREX_GPU_DEVICE (int i, int j, int k) noexcept
        {
            if (vfrac(i,j,k)==zero) {
                flag(i,j,k).setCovered();
            }
        });

    }
}

/**
 * \brief Mask FC faces as covered when both neighboring CC cells are covered.
 *
 * For each face in the FC factory, if BOTH neighboring cell-centered cells
 * are marked as covered, this function sets the FC face's volume fraction to zero
 * and marks it as covered. This maintains consistency between CC and FC EB representations.
 */
void
eb_::mask_fc_from_cc (EBFArrayBoxFactory* fc_factory,
                      const EBFArrayBoxFactory* cc_factory,
                      int idim)
{

    FabArray<EBCellFlagFab>& fc_cellflag = getNonConstEBCellFlags(*fc_factory);
    MultiFab& fc_volfrac = getNonConstVolFrac(*fc_factory);
    const FabArray<EBCellFlagFab>& cc_cellflag = cc_factory->getMultiEBCellFlagFab();

    for (MFIter mfi(fc_cellflag, false); mfi.isValid(); ++mfi) {
        const Box face_box = mfi.nodaltilebox(idim);
        const Box gbx = amrex::grow(face_box, fc_cellflag.nGrow()-1);

        Array4<EBCellFlag> const& fc_flag = fc_cellflag.array(mfi);
        Array4<Real> const& fc_vfrac = fc_volfrac.array(mfi);
        Array4<EBCellFlag const> const& cc_flag = cc_cellflag.const_array(mfi);

        if (idim == 0) {
            ParallelFor(gbx, [=] AMREX_GPU_DEVICE (int i, int j, int k) noexcept {
                if (cc_flag(i-1,j,k).isCovered() && cc_flag(i,j,k).isCovered()) {
                    fc_vfrac(i,j,k) = Real(0.0);
                    fc_flag(i,j,k).setCovered();
                }
            });
        } else if (idim == 1) {
            ParallelFor(gbx, [=] AMREX_GPU_DEVICE (int i, int j, int k) noexcept {
                if (cc_flag(i,j-1,k).isCovered() && cc_flag(i,j,k).isCovered()) {
                    fc_vfrac(i,j,k) = Real(0.0);
                    fc_flag(i,j,k).setCovered();
                }
            });
        } else if (idim == 2) {
            ParallelFor(gbx, [=] AMREX_GPU_DEVICE (int i, int j, int k) noexcept {
                if (cc_flag(i,j,k-1).isCovered() && cc_flag(i,j,k).isCovered()) {
                    fc_vfrac(i,j,k) = Real(0.0);
                    fc_flag(i,j,k).setCovered();
                }
            });
        }
    }
}
