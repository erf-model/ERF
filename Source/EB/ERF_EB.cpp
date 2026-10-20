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
#include <AMReX_Scan.H>

using namespace amrex;

namespace {

#ifdef AMREX_USE_GPU
int
fill_cut_cell_list (Box const& bx, Array4<EBCellFlag const> const& flag,
                    Dim3* grid_cut_cells)
{
    return Scan::PrefixSum<int>(
        bx.numPts(),
        [=] AMREX_GPU_DEVICE (Long offset) noexcept -> int
        {
            Dim3 const cell = bx.atOffset(offset).dim3();
            return flag(cell.x,cell.y,cell.z).isSingleValued();
        },
        [=] AMREX_GPU_DEVICE (Long offset, int cut_offset) noexcept
        {
            Dim3 const cell = bx.atOffset(offset).dim3();
            if (flag(cell.x,cell.y,cell.z).isSingleValued()) {
                grid_cut_cells[cut_offset] = cell;
            }
        },
        Scan::Type::exclusive, Scan::retSum);
}
#endif

} // namespace

eb_::~eb_ ()
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

#if USE_FC_FACTORY
    build_cut_cell_lists(m_u_factory_fc->getMultiEBCellFlagFab(), m_u_cut_cells);
    build_cut_cell_lists(m_v_factory_fc->getMultiEBCellFlagFab(), m_v_cut_cells);
    build_cut_cell_lists(m_w_factory_fc->getMultiEBCellFlagFab(), m_w_cut_cells);
#else
    build_cut_cell_lists(m_u_factory.getMultiEBCellFlagFab(), m_u_cut_cells);
    build_cut_cell_lists(m_v_factory.getMultiEBCellFlagFab(), m_v_cut_cells);
    build_cut_cell_lists(m_w_factory.getMultiEBCellFlagFab(), m_w_cut_cells);
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

    m_u_cut_cells.clear();
    m_v_cut_cells.clear();
    m_w_cut_cells.clear();

    Print() << "\nDone making EB factory at level " << level << ".\n\n";
}

void
eb_::build_cut_cell_lists (FabArray<EBCellFlagFab> const& flags,
                           Vector<CutCellList>& cut_cells)
{
    BL_PROFILE("eb_::build_cut_cell_lists()");

    cut_cells.clear();
    cut_cells.resize(flags.boxArray().size());

#ifdef AMREX_USE_GPU
    for (MFIter mfi(flags, false); mfi.isValid(); ++mfi) {
        Box const& bx = mfi.validbox();
        if (flags[mfi].getType(bx) != FabType::singlevalued) {
            continue;
        }

        Array4<EBCellFlag const> const& flag = flags.const_array(mfi);
        Long const num_candidates = bx.numPts();
        CutCellList grid_cut_cells(num_candidates);
        Dim3* const grid_cut_cells_ptr = grid_cut_cells.data();

        int const num_cut_cells = fill_cut_cell_list(bx, flag, grid_cut_cells_ptr);

        grid_cut_cells.resize(num_cut_cells);
        grid_cut_cells.shrink_to_fit();
        cut_cells[mfi.index()] = std::move(grid_cut_cells);
    }
#else
    amrex::ignore_unused(flags);
#endif
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
