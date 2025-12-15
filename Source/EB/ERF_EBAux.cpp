#include <AMReX_BoxList.H>
#include <AMReX_ParmParse.H>
#include <ERF_EBAux.H>
#include <ERF_EBCutCell.H>
#include <AMReX_MultiFabUtil.H>
#if 0
#include <AMReX_VisMF.H>
#endif

using namespace amrex;

eb_aux_::
~eb_aux_ ()
{
}

eb_aux_::
eb_aux_ ()
  : m_verbose(0)
// ,m_defined(0)
{}

void
eb_aux_::
define( [[maybe_unused]] int const& a_level,
        int const& a_idim,
        Geometry            const& a_geom,
        BoxArray            const& a_grids,
        DistributionMapping const& a_dmap,
        Vector<int>         const& a_ngrow,
        EBFArrayBoxFactory  const* a_factory)
{
  // Box dbox(a_geom.Domain());

  // small_volfrac
  Real small_volfrac = 1.e-14;
  ParmParse pp("eb2");
  pp.queryAdd("small_volfrac", small_volfrac);
  const Real small_value = 1.e-15;

  const IntVect vdim(IntVect::TheDimensionVector(a_idim));

  const BoxArray& my_grids = amrex::convert(a_grids, vdim);

  m_cellflags = new FabArray<EBCellFlagFab>(my_grids, a_dmap, 1, a_ngrow[0], MFInfo(),
                                            DefaultFabFactory<EBCellFlagFab>());

  // Set m_cellflags type to singlevalued
  m_cellflags->setVal(EBCellFlag::TheDefaultCell());
  for (MFIter mfi(*m_cellflags, false); mfi.isValid(); ++mfi) {
    auto& fab = (*m_cellflags)[mfi];
    fab.setType(FabType::singlevalued);
  }

  m_volfrac = new MultiFab(my_grids, a_dmap, 1, a_ngrow[1], MFInfo(), FArrayBoxFactory());
  m_volcent = new MultiFab(my_grids, a_dmap, AMREX_SPACEDIM, a_ngrow[2], MFInfo(), FArrayBoxFactory());

  for (int idim = 0; idim < AMREX_SPACEDIM; ++idim) {
    m_areafrac[idim] = new MultiFab(a_grids, a_dmap,                1, a_ngrow[1]+1, MFInfo(), FArrayBoxFactory());
    m_facecent[idim] = new MultiFab(a_grids, a_dmap, AMREX_SPACEDIM-1, a_ngrow[2], MFInfo(), FArrayBoxFactory());
  }

  m_bndryarea = new MultiFab(my_grids, a_dmap, 1, a_ngrow[2], MFInfo(), FArrayBoxFactory());
  m_bndrycent = new MultiFab(my_grids, a_dmap, AMREX_SPACEDIM, a_ngrow[2], MFInfo(), FArrayBoxFactory());
  m_bndrynorm = new MultiFab(my_grids, a_dmap, AMREX_SPACEDIM, a_ngrow[2], MFInfo(), FArrayBoxFactory());

  // Initialize with zeros
  m_volfrac->setVal(0.0);
  m_volcent->setVal(0.0);

  for (int idim = 0; idim < AMREX_SPACEDIM; ++idim) {
    m_areafrac[idim]->setVal(0.0);
    m_facecent[idim]->setVal(0.0);
  }

  m_bndryarea->setVal(0.0);
  m_bndrycent->setVal(0.0);
  m_bndrynorm->setVal(0.0);

  const auto& FlagFab = a_factory->getMultiEBCellFlagFab(); // EBFArrayBoxFactory, EBDataCollection

  for (MFIter mfi(*m_cellflags, false); mfi.isValid(); ++mfi) {

    const Box& bx = mfi.validbox();
    const Box& bx_grown = mfi.growntilebox();
    const Box tbx = mfi.nodaltilebox(a_idim);
    const Box domain = surroundingNodes(a_geom.Domain(), a_idim);

    GpuArray<Real, AMREX_SPACEDIM> dx = a_geom.CellSizeArray();
    bool l_periodic   = a_geom.isPeriodic(a_idim);

    Array4<EBCellFlag> const& aux_flag  = m_cellflags->array(mfi);
    Array4<Real>       const& aux_vfrac = m_volfrac->array(mfi);
    Array4<Real>       const& aux_afrac_x = m_areafrac[0]->array(mfi);
    Array4<Real>       const& aux_afrac_y = m_areafrac[1]->array(mfi);
    Array4<Real>       const& aux_afrac_z = m_areafrac[2]->array(mfi);

    if (FlagFab[mfi].getType(bx) == FabType::covered ) {

      ParallelFor(tbx, [=] AMREX_GPU_DEVICE (int i, int j, int k) noexcept
      {
        aux_flag(i,j,k).setCovered();
        aux_flag(i,j,k).setDisconnected();
        if (i==bx.bigEnd(0)) {
          aux_flag(i+1,j,k).setCovered();
        }
        if (j==bx.bigEnd(1)) {
          aux_flag(i,j+1,k).setCovered();
        }
        if (k==bx.bigEnd(2)) {
          aux_flag(i,j,k+1).setCovered();
        }
      });

    } else if (FlagFab[mfi].getType(bx) == FabType::regular ) {

      ParallelFor(tbx, [=] AMREX_GPU_DEVICE (int i, int j, int k) noexcept
      {
        aux_flag(i,j,k).setRegular();
        aux_flag(i,j,k).setDisconnected();
        aux_vfrac(i,j,k) = 1.0;
        aux_afrac_x(i,j,k) = 1.0;
        aux_afrac_y(i,j,k) = 1.0;
        aux_afrac_z(i,j,k) = 1.0;
        if (i==bx.bigEnd(0)) {
          aux_flag(i+1,j,k).setRegular();
          aux_vfrac(i+1,j,k) = 1.0;
          aux_afrac_x(i+1,j,k) = 1.0;
        }
        if (j==bx.bigEnd(1)) {
          aux_flag(i,j+1,k).setRegular();
          aux_vfrac(i,j+1,k) = 1.0;
          aux_afrac_y(i,j+1,k) = 1.0;
        }
        if (k==bx.bigEnd(2)) {
          aux_flag(i,j,k+1).setRegular();
          aux_vfrac(i,j,k+1) = 1.0;
          aux_afrac_z(i,j,k+1) = 1.0;
        }
      });

    } else if (FlagFab[mfi].getType(bx) == FabType::singlevalued ) {

      // Initialization

      // CC cell quantities
      Array4<EBCellFlag const> const& flag = FlagFab.const_array(mfi);
      Array4<Real const> const& afrac = (a_factory->getAreaFrac()[a_idim])->const_array(mfi);
      Array4<Real const> const& bnorm = a_factory->getBndryNormal()[mfi].const_array();
      Array4<Real const> const& bcent = a_factory->getBndryCent()[mfi].const_array();

      // aux quantities
      Array4<Real>       const& aux_vcent = m_volcent->array(mfi);
      Array4<Real>       const& aux_fcent_x = m_facecent[0]->array(mfi);
      Array4<Real>       const& aux_fcent_y = m_facecent[1]->array(mfi);
      Array4<Real>       const& aux_fcent_z = m_facecent[2]->array(mfi);
      Array4<Real>       const& aux_barea = m_bndryarea->array(mfi);
      Array4<Real>       const& aux_bcent = m_bndrycent->array(mfi);
      Array4<Real>       const& aux_bnorm = m_bndrynorm->array(mfi);

      // Extended domain in the direction of periodicity
      Box dom_grown = domain;
      for (int idim = 0; idim < AMREX_SPACEDIM; ++idim) {
        if (a_geom.isPeriodic(idim)) {
          dom_grown.grow(idim, a_ngrow[0]);
        }
      }

      const IntVect dom_grown_lo = dom_grown.smallEnd();
      const IntVect dom_grown_hi = dom_grown.bigEnd();

      BoxList diffList = boxDiff(bx_grown, bx);
      for (const Box& b : diffList) {
        ParallelFor(b, [=] AMREX_GPU_DEVICE (int i, int j, int k) noexcept
        {
          if ( i < dom_grown_lo[0] || i > dom_grown_hi[0] ||
               j < dom_grown_lo[1] || j > dom_grown_hi[1] ||
               k < dom_grown_lo[2] || k > dom_grown_hi[2] ) {
            aux_flag(i,j,k).setCovered();
            aux_flag(i,j,k).setDisconnected();
          }
        });
      }

#ifndef AMREX_USE_GPU
      int const verbose=m_verbose;
#endif

      ParallelFor(bx, [=] AMREX_GPU_DEVICE (int i, int j, int k) noexcept
      {
        // defaults to covered and disconnected.
        aux_flag(i,j,k).setCovered();
        aux_flag(i,j,k).setDisconnected();

        if (i==bx.bigEnd(0)) {
          aux_flag(i+1,j,k).setCovered();
        }
        if (j==bx.bigEnd(1)) {
          aux_flag(i,j+1,k).setCovered();
        }
        if (k==bx.bigEnd(2)) {
          aux_flag(i,j,k+1).setCovered();
        }

        // Index for low and hi cells
        IntVect iv_hi(i,j,k);
        IntVect iv_lo(iv_hi - vdim);

        bool lo_isCovered = flag(iv_lo).isCovered();
        bool hi_isCovered = flag(iv_hi).isCovered();
        bool lo_isRegular = flag(iv_lo).isRegular();
        bool hi_isRegular = flag(iv_hi).isRegular();
        bool lo_isSingleValued = flag(iv_lo).isSingleValued();
        bool hi_isSingleValued = flag(iv_hi).isSingleValued();

        const bool at_lo_boundary = (!l_periodic && iv_hi[a_idim]==domain.smallEnd(a_idim));
        const bool at_hi_boundary = (!l_periodic && iv_hi[a_idim]==domain.bigEnd(a_idim));

        // Treatment of lower boundary

        if (at_lo_boundary) {
          if (hi_isCovered) {
            lo_isCovered = true;
            lo_isRegular = false;
            lo_isSingleValued = false;
          } else if (hi_isRegular) {
            lo_isCovered = false;
            lo_isRegular = true;
            lo_isSingleValued = false;
          } else if (hi_isSingleValued) {
            if (almostEqual(afrac(i,j,k),0.0)) {
              lo_isCovered = true;
              lo_isRegular = false;
              lo_isSingleValued = false;
            } else if (almostEqual(afrac(i,j,k),1.0)) {
              lo_isCovered = false;
              lo_isRegular = true;
              lo_isSingleValued = false;
            } else {
              lo_isCovered = false;
              lo_isRegular = false;
              lo_isSingleValued = true;
              iv_lo = iv_hi; // At the lower boundary, low cell takes the values of the high cell.
            }
          }
        }

        // Treatment of upper boundary

        if (at_hi_boundary) {
          if (lo_isCovered) { // Covered
            hi_isCovered = true;
            hi_isRegular = false;
            hi_isSingleValued = false;
          } else if (lo_isRegular) { // Regular
            hi_isCovered = false;
            hi_isRegular = true;
            hi_isSingleValued = false;
          } else if (lo_isSingleValued) { // SingleValued
            if (almostEqual(afrac(i,j,k),0.0)) { //Covered
              hi_isCovered = true;
              hi_isRegular = false;
              hi_isSingleValued = false;
            } else if (almostEqual(afrac(i,j,k),1.0)) { //Regular
              hi_isCovered = false;
              hi_isRegular = true;
              hi_isSingleValued = false;
            } else { // SingleValued
              hi_isCovered = false;
              hi_isRegular = false;
              hi_isSingleValued = true;
              iv_hi = iv_lo; // At the upper boundary, hi cell takes the values of the low cell.
            }
          }
        }

        if ( lo_isCovered && hi_isCovered) {

          // defaults to covered and disconnected.

        } else if ( lo_isRegular && hi_isRegular) {

          aux_flag(i,j,k).setRegular();
          aux_flag(i,j,k).setConnected();

          aux_vfrac(i,j,k) = 1.0;

          aux_afrac_x(i,j,k) = 1.0;
          aux_afrac_y(i,j,k) = 1.0;
          aux_afrac_z(i,j,k) = 1.0;

          if (i==bx.bigEnd(0)) {
            aux_afrac_x(i+1,j,k) = 1.0;
          }
          if (j==bx.bigEnd(1)) {
            aux_afrac_y(i,j+1,k) = 1.0;
          }
          if (k==bx.bigEnd(2)) {
            aux_afrac_z(i,j,k+1) = 1.0;
          }

        } else {

#ifndef AMREX_USE_GPU
          if (verbose) { Print() << "\ncell: " << amrex::IntVect(i,j,k) << "\n"; }
#endif
          Array<Real,AMREX_SPACEDIM> lo_arr = {-0.5,-0.5,-0.5};
          Array<Real,AMREX_SPACEDIM> hi_arr = { 0.5, 0.5, 0.5};

          //-----------------------
          // Low EB cut cell
          //-----------------------

          // Map bcent and bnorm to the isoparametric space for anisotropic grids.
          // (This step is needed because bcent in AMReX is isotropically normalized.)

          RealVect lo_point (bcent(iv_lo,0), bcent(iv_lo,1), bcent(iv_lo,2));
          RealVect lo_normal(bnorm(iv_lo,0), bnorm(iv_lo,1), bnorm(iv_lo,2));

          if (at_lo_boundary) { // At lower boundary
            lo_point[a_idim] += 1.0; // Move the boundary centroid upward in the a_idim direction.
          }

          if (lo_isSingleValued ) {
            Real bnorm_x = bnorm(iv_lo,0) * dx[0];
            Real bnorm_y = bnorm(iv_lo,1) * dx[1];
            Real bnorm_z = bnorm(iv_lo,2) * dx[2];

            Real norm = sqrt( bnorm_x*bnorm_x + bnorm_y*bnorm_y + bnorm_z*bnorm_z);

            RealVect bnorm_isoparam ( bnorm_x / norm, bnorm_y / norm, bnorm_z / norm);

            lo_normal = bnorm_isoparam;
          }

          // High side of low cell
          lo_arr[a_idim] = 0.0;
          hi_arr[a_idim] = 0.5;
          RealBox lo_rbx(lo_arr.data(), hi_arr.data());

          eb_cut_cell_ lo_eb_cc(flag(iv_lo), lo_rbx, lo_point, lo_normal);

          // cell iv_lo covered (regular) imples lo_eb_cc is covered (regular)
          // The inverse is not always true.
          AMREX_ASSERT( !lo_isCovered || lo_eb_cc.isCovered() );
          AMREX_ASSERT( !lo_isRegular || lo_eb_cc.isRegular() );

          //-----------------------
          // High EB cut cell
          //-----------------------

          RealVect hi_point (bcent(iv_hi,0), bcent(iv_hi,1), bcent(iv_hi,2));
          RealVect hi_normal(bnorm(iv_hi,0), bnorm(iv_hi,1), bnorm(iv_hi,2));

          if (at_hi_boundary) {
            hi_point[a_idim] += -1.0; // Move the boundary centroid downward in the a_idim direction.
          }

          if (hi_isSingleValued ) {
            Real bnorm_x = bnorm(iv_hi,0) * dx[0];
            Real bnorm_y = bnorm(iv_hi,1) * dx[1];
            Real bnorm_z = bnorm(iv_hi,2) * dx[2];

            Real norm = sqrt( bnorm_x*bnorm_x + bnorm_y*bnorm_y + bnorm_z*bnorm_z);

            RealVect bnorm_isoparam ( bnorm_x / norm, bnorm_y / norm, bnorm_z / norm);

            hi_normal = bnorm_isoparam;
          }

          // Low side of high cell
          lo_arr[a_idim] = -0.5;
          hi_arr[a_idim] =  0.0;
          RealBox hi_rbx(lo_arr.data(), hi_arr.data());

          eb_cut_cell_ hi_eb_cc(flag(iv_hi), hi_rbx, hi_point, hi_normal);

          // cell iv_hi covered (regular) imples hi_eb_cc is covered (regular)
          // The inverse is not always true.
          AMREX_ASSERT( !hi_isCovered || hi_eb_cc.isCovered() );
          AMREX_ASSERT( !hi_isRegular || hi_eb_cc.isRegular() );

#if 0
#if defined(AMREX_DEBUG) || defined(AMREX_TESTING) || 1

          { /***************************** SANITY CHECK ***********************\
            * Perform some basic sanity checks to verify that what we computed *
            * for cell (i,j,k) compares to what we know to be true.           *
            \******************************************************************/

            // Compute the cut-cell for the high side of the high cell. This is
            // only needed for sanity checks.

            eb_cut_cell_ hi_hi_eb_cc(flag(iv_hi), lo_rbx, hi_point, hi_normal);

            // cell iv_hi covered (regular) imples hi_hi_eb_cc is covered (regular)
            // The inverse is not always true.
#ifndef AMREX_USE_GPU
            if ( !(!hi_isRegular || hi_hi_eb_cc.isRegular()) ||
                 !(!hi_isCovered || hi_hi_eb_cc.isCovered()) ) {
              Print() << "flag(iv_hi) and hi_hi_eb_cc flags do not agree\n"
                      << "\n  isRegular() " << hi_isRegular << "  " << hi_hi_eb_cc.isRegular()
                      << "\n  isCovered() " << hi_isCovered << "  " << hi_hi_eb_cc.isCovered()
                      << "\n";
            }
#endif
            // If cell iv_hi is regular or covered, then hi_hi_eb_cc must also
            // be regular or covered. The inverse is not true.
            AMREX_ALWAYS_ASSERT( !hi_isRegular || hi_hi_eb_cc.isRegular() );
            AMREX_ALWAYS_ASSERT( !hi_isCovered || hi_hi_eb_cc.isCovered() );

            // The area and volume fractions that are computed for the scalar grid
            // are slightly different than those we compute from the geometric
            // reconstruction using the EB point and normal. However, we expect
            // that the area fractions computed here will give back the same
            // normal we used to compute them.
            if ( hi_isSingleValued ) {

              Real const adx = (a_idim == 0)
                             ? (hi_eb_cc.areaLo(0) - hi_hi_eb_cc.areaHi(0)) * dx[1] * dx[2]
                             : (hi_eb_cc.areaLo(0) + hi_hi_eb_cc.areaLo(0)) * dx[1] * dx[2]
                             - (hi_eb_cc.areaHi(0) + hi_hi_eb_cc.areaHi(0)) * dx[1] * dx[2];

              Real const ady = (a_idim == 1)
                             ? (hi_eb_cc.areaLo(1) - hi_hi_eb_cc.areaHi(1)) * dx[0] * dx[2]
                             : (hi_eb_cc.areaLo(1) + hi_hi_eb_cc.areaLo(1)) * dx[0] * dx[2]
                             - (hi_eb_cc.areaHi(1) + hi_hi_eb_cc.areaHi(1)) * dx[0] * dx[2];

              Real const adz = (a_idim == 2)
                             ? (hi_eb_cc.areaLo(2) - hi_hi_eb_cc.areaHi(2)) * dx[0] * dx[1]
                             : (hi_eb_cc.areaLo(2) + hi_hi_eb_cc.areaLo(2)) * dx[0] * dx[1]
                             - (hi_eb_cc.areaHi(2) + hi_hi_eb_cc.areaHi(2)) * dx[0] * dx[1];

              Real const apnorm = std::sqrt(adx*adx + ady*ady + adz*adz);

              // EB normal
              Real const apnorminv = 1. / apnorm;
              RealVect const normal(adx*apnorminv, ady*apnorminv, adz*apnorminv);
              Real const dot_normals = normal.dotProduct(hi_normal);

#ifndef AMREX_USE_GPU
              if ( !amrex::almostEqual(dot_normals, 1.0) ) {
                Print() << "\nFail: check-1 dot_normals " << dot_normals
                        << '\n';

                hi_eb_cc.debug();
                hi_hi_eb_cc.debug();

              } else if (verbose) {
                Print() << "Pass: dot_normals = 1.0\n";

              }
#endif
              AMREX_ALWAYS_ASSERT( amrex::almostEqual(dot_normals, 1.0) );
            }

            // The a_idim area of hi_eb_cc.areaHi() should equal hi_hi_eb_cc.areaLo()
            {
#ifndef AMREX_USE_GPU
            Real const abs_err = std::abs( hi_eb_cc.areaHi(a_idim) - hi_hi_eb_cc.areaLo(a_idim) );
            Real machine_tol = 10.0*std::numeric_limits<amrex::Real>::epsilon();
            if ( abs_err >= machine_tol ) {
                Print() << "\nFail: check-2 area abs_err: " << abs_err
                        << "\n  hi_eb_cc.areaHi " << hi_eb_cc.areaHi(a_idim)
                        << "\n  hi_hi_eb_cc.areaLo " << hi_hi_eb_cc.areaLo(a_idim)
                        << '\n';
            } else if (verbose) {
                Print() << "Pass: hi_eb_cc.areaHi = hi_hi_eb_cc.areaLo"
                        << "  abs_err: " << abs_err << "\n";
            }
            AMREX_ALWAYS_ASSERT( abs_err < machine_tol );
#endif
            }

            // The low-side area of hi_eb_cc should equal a_idim afrac.
            { Real const abs_err = amrex::max(std::abs(lo_eb_cc.areaHi(a_idim) - afrac(iv_hi)),
                                              std::abs(hi_eb_cc.areaLo(a_idim) - afrac(iv_hi)));
              Real compare_tol = 5.0e-6;
#ifndef AMREX_USE_GPU
              if ( abs_err >= compare_tol ) {
                //hi_eb_cc.debug();
                Print() << "\nFail: check-3 area abs_err " << abs_err
                        << "\n  hi_eb_cc.areaLo(" << a_idim << ") = " << hi_eb_cc.areaLo(a_idim)
                        << "\n  lo_eb_cc.areaHi(" << a_idim << ") = " << lo_eb_cc.areaHi(a_idim)
                        << "\n  afrac" << iv_hi << " =  " << afrac(iv_hi)
                        << '\n';
              } else if (verbose) {
                Print() << "Pass: hi_eb_cc.areaLo = afrac = " << afrac(iv_hi)
                        << "  abs_err: " << abs_err << "\n";
              }
#endif
              AMREX_ALWAYS_ASSERT( abs_err < compare_tol );
            }

            // The combined volumes of hi_eb_cc.areaHi() and hi_hi_eb_cc should
            // equal vfrac(iv_hi).
            { Real const vol = hi_eb_cc.volume() + hi_hi_eb_cc.volume();
              Real const abs_err = amrex::Math::abs(vfrac(iv_hi) - vol);
              Real compare_tol = 5.0e-6;
#ifndef AMREX_USE_GPU
              if ( abs_err >= compare_tol ) {
                hi_eb_cc.debug();
                hi_hi_eb_cc.debug();
                amrex::Print() << "\nFail: check-4 volume abs_err: " << abs_err
                               << "\n  point:  " << hi_point
                               << "\n  normal: " << hi_normal
                               << "\n     hi_eb_cc.volume() " <<    hi_eb_cc.volume()
                               << "\n  hi_hi_eb_cc.volume() " << hi_hi_eb_cc.volume()
                               << "\n  vfrac:   " << vfrac(iv_hi)
                               << '\n';
              } else if (verbose) {
                Print() << "Pass: hi_eb_cc + hi_hi_eb_cc = vfrac = " << vfrac(iv_hi)
                        << "  abs_err: " << abs_err << "\n";
              }
#endif
              AMREX_ALWAYS_ASSERT( abs_err < compare_tol );
            }
          } //
#endif
#endif // 0

          //-----------------------
          // Fill out aux_ arrays
          //-----------------------

          if (lo_eb_cc.isCovered() && hi_eb_cc.isCovered()) {

            // defaults to covered and disconnected.

          } else if (lo_eb_cc.isRegular() && hi_eb_cc.isRegular()) {

            aux_flag(i,j,k).setRegular();
            aux_flag(i,j,k).setConnected();

            aux_vfrac(i,j,k) = 1.0;

            aux_afrac_x(i,j,k) = 1.0;
            aux_afrac_y(i,j,k) = 1.0;
            aux_afrac_z(i,j,k) = 1.0;

            aux_fcent_x(i,j,k,0) = 0.0; aux_fcent_x(i,j,k,1) = 0.0;
            aux_fcent_y(i,j,k,0) = 0.0; aux_fcent_y(i,j,k,1) = 0.0;
            aux_fcent_z(i,j,k,0) = 0.0; aux_fcent_z(i,j,k,1) = 0.0;

            if (i==bx.bigEnd(0)) {
              aux_afrac_x(i+1,j,k) = 1.0;
              aux_fcent_x(i+1,j,k,0) = 0.0; aux_fcent_x(i+1,j,k,1) = 0.0;
            }
            if (j==bx.bigEnd(1)) {
              aux_afrac_y(i,j+1,k) = 1.0;
              aux_fcent_y(i,j+1,k,0) = 0.0; aux_fcent_y(i,j+1,k,1) = 0.0;
            }
            if (k==bx.bigEnd(2)) {
              aux_afrac_z(i,j,k+1) = 1.0;
              aux_fcent_z(i,j,k+1,0) = 0.0; aux_fcent_z(i,j,k+1,1) = 0.0;
            }

          } else if ( (lo_eb_cc.isRegular() && hi_eb_cc.isCovered())
                   || (lo_eb_cc.isCovered() && hi_eb_cc.isRegular()) ) {

            // This is a problematic situation.
#ifndef AMREX_USE_GPU
            Print()<< "eb_aux_ / Check: Regular and Covered cut cells are facing each other." << std::endl;
#endif

          } else {

            // 0. Cell Flag

            aux_flag(i,j,k).setSingleValued();

            // 1. Volume Fraction

            Real lo_vol {lo_eb_cc.volume()}; AMREX_ASSERT(lo_vol >= 0.0 && lo_vol <= 0.5);
            Real hi_vol {hi_eb_cc.volume()}; AMREX_ASSERT(hi_vol >= 0.0 && hi_vol <= 0.5);

            aux_vfrac(i,j,k) = lo_vol + hi_vol;

            // 2. Volume Centroid

            /* centVol() returns the coordinates based on m_rbx.
              The coordinates in the a_idim direction are in [0.0,0.5] for the low cell and in [-0.5,0.0] for the hi cell.
              Therefore, they need to be mapped to the eb_aux space, by shifting:
              x' = x - 0.5 (low cell), x + 0.5 (hi cell) if a_idim = 0
              y' = y - 0.5 (low cell), y + 0.5 (hi cell) if a_idim = 1
              z' = z - 0.5 (low cell), z + 0.5 (hi cell) if a_idim = 2
            */

            RealVect lo_vcent {lo_eb_cc.centVol()};
            RealVect hi_vcent {hi_eb_cc.centVol()};

            lo_vcent[a_idim] = lo_vcent[a_idim] - 0.5;
            hi_vcent[a_idim] = hi_vcent[a_idim] + 0.5;

            aux_vcent(i,j,k,0) = ( lo_vol * lo_vcent[0] + hi_vol * hi_vcent[0] ) / aux_vfrac(i,j,k);
            aux_vcent(i,j,k,1) = ( lo_vol * lo_vcent[1] + hi_vol * hi_vcent[1] ) / aux_vfrac(i,j,k);
            aux_vcent(i,j,k,2) = ( lo_vol * lo_vcent[2] + hi_vol * hi_vcent[2] ) / aux_vfrac(i,j,k);

            // 3. Area Fraction

            Real lo_areaLo_x {lo_eb_cc.areaLo(0)};
            Real lo_areaLo_y {lo_eb_cc.areaLo(1)};
            Real lo_areaLo_z {lo_eb_cc.areaLo(2)};

            Real hi_areaLo_x {hi_eb_cc.areaLo(0)};
            Real hi_areaLo_y {hi_eb_cc.areaLo(1)};
            Real hi_areaLo_z {hi_eb_cc.areaLo(2)};

            aux_afrac_x(i,j,k) = (a_idim == 0) ? lo_areaLo_x : lo_areaLo_x + hi_areaLo_x;
            aux_afrac_y(i,j,k) = (a_idim == 1) ? lo_areaLo_y : lo_areaLo_y + hi_areaLo_y;
            aux_afrac_z(i,j,k) = (a_idim == 2) ? lo_areaLo_z : lo_areaLo_z + hi_areaLo_z;

            if (i==bx.bigEnd(0)) {
              Real lo_areaHi_x {lo_eb_cc.areaHi(0)};
              Real hi_areaHi_x {hi_eb_cc.areaHi(0)};
              aux_afrac_x(i+1,j,k) = (a_idim == 0) ? hi_areaHi_x : lo_areaHi_x + hi_areaHi_x;
            }
            if (j==bx.bigEnd(1)) {
              Real lo_areaHi_y {lo_eb_cc.areaHi(1)};
              Real hi_areaHi_y {hi_eb_cc.areaHi(1)};
              aux_afrac_y(i,j+1,k) = (a_idim == 1) ? hi_areaHi_y : lo_areaHi_y + hi_areaHi_y;
            }
            if (k==bx.bigEnd(2)) {
              Real lo_areaHi_z {lo_eb_cc.areaHi(2)};
              Real hi_areaHi_z {hi_eb_cc.areaHi(2)};
              aux_afrac_z(i,j,k+1) = (a_idim == 2) ? hi_areaHi_z : lo_areaHi_z + hi_areaHi_z;
            }

            // 4. Face Centroid

            /* fcentLo returns the coordinates based on m_rbx.
              The coordinates in the a_idim direction are in [0.0,0.5] for the low cell and in [-0.5,0.0] for the hi cell.
              Therefore, they need to be mapped to the eb_aux space, by shifting:
              x' = x - 0.5 (low cell), x + 0.5 (hi cell) if a_idim = 0
              y' = y - 0.5 (low cell), y + 0.5 (hi cell) if a_idim = 1
              z' = z - 0.5 (low cell), z + 0.5 (hi cell) if a_idim = 2
            */

            RealVect lo_centLo_x {lo_eb_cc.centLo(0)};
            RealVect lo_centLo_y {lo_eb_cc.centLo(1)};
            RealVect lo_centLo_z {lo_eb_cc.centLo(2)};

            RealVect hi_centLo_x {hi_eb_cc.centLo(0)};
            RealVect hi_centLo_y {hi_eb_cc.centLo(1)};
            RealVect hi_centLo_z {hi_eb_cc.centLo(2)};

            if (a_idim == 0) {
              aux_fcent_x(i,j,k,0) = lo_centLo_x[1];      // y
              aux_fcent_x(i,j,k,1) = lo_centLo_x[2];      // z
              aux_fcent_y(i,j,k,0) = (aux_afrac_y(i,j,k) > 0.0)   // x (mapped)
                                    ? ( lo_areaLo_y * (lo_centLo_y[0] - 0.5)
                                      + hi_areaLo_y * (hi_centLo_y[0] + 0.5) ) / aux_afrac_y(i,j,k)
                                    : 0.0;
              aux_fcent_y(i,j,k,1) = (aux_afrac_y(i,j,k) > 0.0)   // z
                                    ? ( lo_areaLo_y * lo_centLo_y[2]
                                      + hi_areaLo_y * hi_centLo_y[2] ) / aux_afrac_y(i,j,k)
                                    : 0.0;
              aux_fcent_z(i,j,k,0) = (aux_afrac_z(i,j,k) > 0.0)   // x (mapped)
                                    ? ( lo_areaLo_z * (lo_centLo_z[0] - 0.5)
                                      + hi_areaLo_z * (hi_centLo_z[0] + 0.5) ) / aux_afrac_z(i,j,k)
                                    : 0.0;
              aux_fcent_z(i,j,k,1) = (aux_afrac_z(i,j,k) > 0.0)   // y
                                    ? ( lo_areaLo_z * lo_centLo_z[1]
                                      + hi_areaLo_z * hi_centLo_z[1] ) / aux_afrac_z(i,j,k)
                                    : 0.0;
            } else if (a_idim == 1) {
              aux_fcent_x(i,j,k,0) = (aux_afrac_x(i,j,k) > 0.0)   // y (mapped)
                                    ? ( lo_areaLo_x * (lo_centLo_x[1] - 0.5)
                                      + hi_areaLo_x * (hi_centLo_x[1] + 0.5) ) / aux_afrac_x(i,j,k)
                                    : 0.0;
              aux_fcent_x(i,j,k,1) = (aux_afrac_x(i,j,k) > 0.0)   // z
                                    ? ( lo_areaLo_x * lo_centLo_x[2]
                                      + hi_areaLo_x * hi_centLo_x[2] ) / aux_afrac_x(i,j,k)
                                    : 0.0;
              aux_fcent_y(i,j,k,0) = lo_centLo_y[0];      // x
              aux_fcent_y(i,j,k,1) = lo_centLo_y[2];      // z
              aux_fcent_z(i,j,k,0) = (aux_afrac_z(i,j,k) > 0.0)   // x
                                    ? ( lo_areaLo_z * lo_centLo_z[0]
                                      + hi_areaLo_z * hi_centLo_z[0] ) / aux_afrac_z(i,j,k)
                                    : 0.0;
              aux_fcent_z(i,j,k,1) = (aux_afrac_z(i,j,k) > 0.0)   // y (mapped)
                                    ? ( lo_areaLo_z * (lo_centLo_z[1] - 0.5)
                                      + hi_areaLo_z * (hi_centLo_z[1] + 0.5) ) / aux_afrac_z(i,j,k)
                                    : 0.0;
            } else if (a_idim == 2) {
              aux_fcent_x(i,j,k,0) = (aux_afrac_x(i,j,k) > 0.0)   // y
                                    ? ( lo_areaLo_x * lo_centLo_x[1]
                                      + hi_areaLo_x * hi_centLo_x[1] ) / aux_afrac_x(i,j,k)
                                    : 0.0;
              aux_fcent_x(i,j,k,1) = (aux_afrac_x(i,j,k) > 0.0)   // z (mapped)
                                    ? ( lo_areaLo_x * (lo_centLo_x[2] - 0.5)
                                      + hi_areaLo_x * (hi_centLo_x[2] + 0.5) ) / aux_afrac_x(i,j,k)
                                    : 0.0;
              aux_fcent_y(i,j,k,0) = (aux_afrac_y(i,j,k) > 0.0)   // x
                                    ? ( lo_areaLo_y * lo_centLo_y[0]
                                      + hi_areaLo_y * hi_centLo_y[0] ) / aux_afrac_y(i,j,k)
                                    : 0.0;
              aux_fcent_y(i,j,k,1) = (aux_afrac_y(i,j,k) > 0.0)   // z (mapped)
                                    ? ( lo_areaLo_y * (lo_centLo_y[2] - 0.5)
                                      + hi_areaLo_y * (hi_centLo_y[2] + 0.5) ) / aux_afrac_y(i,j,k)
                                    : 0.0;
              aux_fcent_z(i,j,k,0) = lo_centLo_z[0];      // x
              aux_fcent_z(i,j,k,1) = lo_centLo_z[1];      // y
            }

            if (i==bx.bigEnd(0)) {
              Real lo_areaHi_x {lo_eb_cc.areaHi(0)};
              Real hi_areaHi_x {hi_eb_cc.areaHi(0)};
              RealVect lo_centHi_x {lo_eb_cc.centHi(0)};
              RealVect hi_centHi_x {hi_eb_cc.centHi(0)};
              if (a_idim == 0) {
                aux_fcent_x(i+1,j,k,0) = hi_centHi_x[1];      // y
                aux_fcent_x(i+1,j,k,1) = hi_centHi_x[2];      // z
              } else if (a_idim == 1) {
                aux_fcent_x(i+1,j,k,0) = (aux_afrac_x(i+1,j,k) > 0.0)   // y (mapped)
                                      ? ( lo_areaHi_x * (lo_centHi_x[1] - 0.5)
                                        + hi_areaHi_x * (hi_centHi_x[1] + 0.5) ) / aux_afrac_x(i+1,j,k)
                                      : 0.0;
                aux_fcent_x(i+1,j,k,1) = (aux_afrac_x(i+1,j,k) > 0.0)   // z
                                      ? ( lo_areaHi_x * lo_centHi_x[2]
                                        + hi_areaHi_x * hi_centHi_x[2] ) / aux_afrac_x(i+1,j,k)
                                      : 0.0;
              } else if (a_idim == 2) {
                aux_fcent_x(i+1,j,k,0) = (aux_afrac_x(i+1,j,k) > 0.0)   // y
                                      ? ( lo_areaHi_x * lo_centHi_x[1]
                                        + hi_areaHi_x * hi_centHi_x[1] ) / aux_afrac_x(i+1,j,k)
                                      : 0.0;
                aux_fcent_x(i+1,j,k,1) = (aux_afrac_x(i+1,j,k) > 0.0)   // z (mapped)
                                      ? ( lo_areaHi_x * (lo_centHi_x[2] - 0.5)
                                        + hi_areaHi_x * (hi_centHi_x[2] + 0.5) ) / aux_afrac_x(i+1,j,k)
                                      : 0.0;
              }
            }
            if (j==bx.bigEnd(1)) {
              Real lo_areaHi_y {lo_eb_cc.areaHi(1)};
              Real hi_areaHi_y {hi_eb_cc.areaHi(1)};
              RealVect lo_centHi_y {lo_eb_cc.centHi(1)};
              RealVect hi_centHi_y {hi_eb_cc.centHi(1)};
              if (a_idim == 0) {
                aux_fcent_y(i,j+1,k,0) = (aux_afrac_y(i,j+1,k) > 0.0)   // x (mapped)
                                      ? ( lo_areaHi_y * (lo_centHi_y[0] - 0.5)
                                        + hi_areaHi_y * (hi_centHi_y[0] + 0.5) ) / aux_afrac_y(i,j+1,k)
                                      : 0.0;
                aux_fcent_y(i,j+1,k,1) = (aux_afrac_y(i,j+1,k) > 0.0)   // z
                                      ? ( lo_areaHi_y * lo_centHi_y[2]
                                        + hi_areaHi_y * hi_centHi_y[2] ) / aux_afrac_y(i,j+1,k)
                                      : 0.0;
              } else if (a_idim == 1) {
                aux_fcent_y(i,j+1,k,0) = lo_centHi_y[0];      // x
                aux_fcent_y(i,j+1,k,1) = lo_centHi_y[2];      // z
              } else if (a_idim == 2) {
                aux_fcent_y(i,j+1,k,0) = (aux_afrac_y(i,j+1,k) > 0.0)   // x
                                      ? ( lo_areaHi_y * lo_centHi_y[0]
                                        + hi_areaHi_y * hi_centHi_y[0] ) / aux_afrac_y(i,j+1,k)
                                      : 0.0;
                aux_fcent_y(i,j+1,k,1) = (aux_afrac_y(i,j+1,k) > 0.0)   // z (mapped)
                                      ? ( lo_areaHi_y * (lo_centHi_y[2] - 0.5)
                                        + hi_areaHi_y * (hi_centHi_y[2] + 0.5) ) / aux_afrac_y(i,j+1,k)
                                      : 0.0;
              }
            }
            if (k==bx.bigEnd(2)) {
              Real lo_areaHi_z {lo_eb_cc.areaHi(2)};
              Real hi_areaHi_z {hi_eb_cc.areaHi(2)};
              RealVect lo_centHi_z {lo_eb_cc.centHi(2)};
              RealVect hi_centHi_z {hi_eb_cc.centHi(2)};
              if (a_idim == 0) {
                aux_fcent_z(i,j,k+1,0) = (aux_afrac_z(i,j,k+1) > 0.0)   // x (mapped)
                                      ? ( lo_areaHi_z * (lo_centHi_z[0] - 0.5)
                                        + hi_areaHi_z * (hi_centHi_z[0] + 0.5) ) / aux_afrac_z(i,j,k+1)
                                      : 0.0;
                aux_fcent_z(i,j,k+1,1) = (aux_afrac_z(i,j,k+1) > 0.0)   // y
                                      ? ( lo_areaHi_z * lo_centHi_z[1]
                                        + hi_areaHi_z * hi_centHi_z[1] ) / aux_afrac_z(i,j,k+1)
                                      : 0.0;
              } else if (a_idim == 1) {
                aux_fcent_z(i,j,k+1,0) = (aux_afrac_z(i,j,k+1) > 0.0)   // x
                                      ? ( lo_areaHi_z * lo_centHi_z[0]
                                        + hi_areaHi_z * hi_centHi_z[0] ) / aux_afrac_z(i,j,k+1)
                                      : 0.0;
                aux_fcent_z(i,j,k+1,1) = (aux_afrac_z(i,j,k+1) > 0.0)   // y (mapped)
                                      ? ( lo_areaHi_z * (lo_centHi_z[1] - 0.5)
                                        + hi_areaHi_z * (hi_centHi_z[1] + 0.5) ) / aux_afrac_z(i,j,k+1)
                                      : 0.0;
              } else if (a_idim == 2) {
                aux_fcent_z(i,j,k+1,0) = lo_centHi_z[0];      // x
                aux_fcent_z(i,j,k+1,1) = lo_centHi_z[1];      // y
              }
            }

            // 5. Boundary Area

            Real lo_areaBoun {lo_eb_cc.areaBoun()};
            Real hi_areaBoun {hi_eb_cc.areaBoun()};

            aux_barea(i,j,k) = lo_areaBoun + hi_areaBoun;

            // 6. Boundary Centroid

            RealVect lo_centBoun {lo_eb_cc.centBoun()};
            RealVect hi_centBoun {hi_eb_cc.centBoun()};

            if (a_idim == 0) {
              aux_bcent(i,j,k,0) = ( lo_areaBoun * (lo_centBoun[0]-0.5) + hi_areaBoun * (hi_centBoun[0]+0.5) ) / aux_barea(i,j,k);  // x (mapped)
              aux_bcent(i,j,k,1) = ( lo_areaBoun * lo_centBoun[1] + hi_areaBoun * hi_centBoun[1] ) / aux_barea(i,j,k);              // y
              aux_bcent(i,j,k,2) = ( lo_areaBoun * lo_centBoun[2] + hi_areaBoun * hi_centBoun[2] ) / aux_barea(i,j,k);              // z
            } else if (a_idim == 1) {
              aux_bcent(i,j,k,0) = ( lo_areaBoun * lo_centBoun[0] + hi_areaBoun * hi_centBoun[0] ) / aux_barea(i,j,k);              // x
              aux_bcent(i,j,k,1) = ( lo_areaBoun * (lo_centBoun[1]-0.5) + hi_areaBoun * (hi_centBoun[1]+0.5) ) / aux_barea(i,j,k);  // y (mapped)
              aux_bcent(i,j,k,2) = ( lo_areaBoun * lo_centBoun[2] + hi_areaBoun * hi_centBoun[2] ) / aux_barea(i,j,k);              // z
            } else if (a_idim == 2) {
              aux_bcent(i,j,k,0) = ( lo_areaBoun * lo_centBoun[0] + hi_areaBoun * hi_centBoun[0] ) / aux_barea(i,j,k);              // x
              aux_bcent(i,j,k,1) = ( lo_areaBoun * lo_centBoun[1] + hi_areaBoun * hi_centBoun[1] ) / aux_barea(i,j,k);              // y
              aux_bcent(i,j,k,2) = ( lo_areaBoun * (lo_centBoun[2]-0.5) + hi_areaBoun * (hi_centBoun[2]+0.5) ) / aux_barea(i,j,k);  // z (mapped)
            }

            // 7. Boundary Normal

            RealVect eb_normal = ( lo_areaBoun * lo_normal + hi_areaBoun * hi_normal )/ aux_barea(i,j,k);

            aux_bnorm(i,j,k,0) = eb_normal[0];
            aux_bnorm(i,j,k,1) = eb_normal[1];
            aux_bnorm(i,j,k,2) = eb_normal[2];

          }

        } // flag(iv_lo) and flag(iv_hi)

      });

      ParallelFor(bx, [=] AMREX_GPU_DEVICE (int i, int j, int k) noexcept
      {
        if (aux_vfrac(i,j,k) < small_volfrac)
        {
          aux_vfrac(i,j,k)   = 0.0;
        }
      });

    } // if (FlagFab[mfi].getType(bx) == FabType::singlevalued )

  } // MFIter

  // We FillBoundary volfrac here so that we can use tests on volfrac in ghost cells below
  m_volfrac->FillBoundary(a_geom.periodicity());

  for (MFIter mfi(*m_cellflags, false); mfi.isValid(); ++mfi) {

    const Box& bx = mfi.validbox();
    const Box& bx_grown = mfi.growntilebox();
    const Box domain = surroundingNodes(a_geom.Domain(), a_idim);
    const int dom_lo_i = domain.smallEnd(0);
    const int dom_hi_i = domain.bigEnd(0);
    const int dom_lo_j = domain.smallEnd(1);
    const int dom_hi_j = domain.bigEnd(1);
    const int dom_lo_k = domain.smallEnd(2);
    const int dom_hi_k = domain.bigEnd(2);

    Array4<EBCellFlag> const& aux_flag  = m_cellflags->array(mfi);
    Array4<Real>       const& aux_vfrac = m_volfrac->array(mfi);
    Array4<Real>       const& aux_afrac_x = m_areafrac[0]->array(mfi);
    Array4<Real>       const& aux_afrac_y = m_areafrac[1]->array(mfi);
    Array4<Real>       const& aux_afrac_z = m_areafrac[2]->array(mfi);

    Array4<Real>       const& aux_vcent = m_volcent->array(mfi);
    Array4<Real>       const& aux_fcent_x = m_facecent[0]->array(mfi);
    Array4<Real>       const& aux_fcent_y = m_facecent[1]->array(mfi);
    Array4<Real>       const& aux_fcent_z = m_facecent[2]->array(mfi);
    Array4<Real>       const& aux_barea = m_bndryarea->array(mfi);
    Array4<Real>       const& aux_bcent = m_bndrycent->array(mfi);
    Array4<Real>       const& aux_bnorm = m_bndrynorm->array(mfi);

    if (FlagFab[mfi].getType(bx) == FabType::singlevalued ) {

      // Corrections for small cells
      Box my_xbx(bx); my_xbx.growHi(0,1);
      ParallelFor(my_xbx, [=] AMREX_GPU_DEVICE (int i, int j, int k) noexcept
      {
        if (aux_vfrac(i,j,k) < small_volfrac || aux_vfrac(i-1,j,k) < small_volfrac) {
          // At domain boundary, keep area fraction as is unless inside cell is small
          if ((i == dom_lo_i && aux_vfrac(i,j,k) < small_volfrac) ||
              (i == dom_hi_i+1 && aux_vfrac(i-1,j,k) < small_volfrac) ||
              (i != dom_lo_i && i != dom_hi_i+1)) {
              aux_afrac_x(i,j,k) = 0.0;
          }
        }
      });

      Box my_ybx(bx); my_ybx.growHi(1,1);
      ParallelFor(my_ybx, [=] AMREX_GPU_DEVICE (int i, int j, int k) noexcept
      {
        if (aux_vfrac(i,j,k) < small_volfrac || aux_vfrac(i,j-1,k) < small_volfrac) {
          // At domain boundary, keep area fraction as is unless inside cell is small
          if ((j == dom_lo_j && aux_vfrac(i,j,k) < small_volfrac) ||
              (j == dom_hi_j+1 && aux_vfrac(i,j-1,k) < small_volfrac) ||
              (j != dom_lo_j && j != dom_hi_j+1)) {
              aux_afrac_y(i,j,k) = 0.0;
          }
        }
      });

      Box my_zbx(bx); my_zbx.growHi(2,1);
      ParallelFor(my_zbx, [=] AMREX_GPU_DEVICE (int i, int j, int k) noexcept
      {
        if (aux_vfrac(i,j,k) < small_volfrac || aux_vfrac(i,j,k-1) < small_volfrac) {
          // At domain boundary, keep area fraction as is unless inside cell is small
          if ((k == dom_lo_k && aux_vfrac(i,j,k) < small_volfrac) ||
              (k == dom_hi_k+1 && aux_vfrac(i,j,k-1) < small_volfrac) ||
              (k != dom_lo_k && k != dom_hi_k+1)) {
              aux_afrac_z(i,j,k) = 0.0;
          }
        }
      });

      ParallelFor(bx, [=] AMREX_GPU_DEVICE (int i, int j, int k) noexcept
      {
        if (aux_vfrac(i,j,k) < small_volfrac)
        {
          aux_vcent(i,j,k,0) = 0.0;
          aux_vcent(i,j,k,1) = 0.0;
          aux_vcent(i,j,k,2) = 0.0;

          aux_fcent_x(i  ,j  ,k  ,0) = 0.0;
          aux_fcent_x(i  ,j  ,k  ,1) = 0.0;
          aux_fcent_x(i+1,j  ,k  ,0) = 0.0;
          aux_fcent_x(i+1,j  ,k  ,1) = 0.0;

          aux_fcent_y(i  ,j  ,k  ,0) = 0.0;
          aux_fcent_y(i  ,j  ,k  ,1) = 0.0;
          aux_fcent_y(i  ,j+1,k  ,0) = 0.0;
          aux_fcent_y(i  ,j+1,k  ,1) = 0.0;

          aux_fcent_z(i  ,j  ,k  ,0) = 0.0;
          aux_fcent_z(i  ,j  ,k  ,1) = 0.0;
          aux_fcent_z(i  ,j  ,k+1,0) = 0.0;
          aux_fcent_z(i  ,j  ,k+1,1) = 0.0;

          aux_barea(i,j,k) = 0.0;

          aux_bcent(i,j,k,0) = 0.0;
          aux_bcent(i,j,k,1) = 0.0;
          aux_bcent(i,j,k,2) = 0.0;

          aux_bnorm(i,j,k,0) = 0.0;
          aux_bnorm(i,j,k,1) = 0.0;
          aux_bnorm(i,j,k,2) = 0.0;

          aux_flag(i,j,k).setCovered();
        }

        if (std::abs(aux_vcent(i,j,k,0)) < small_value) aux_vcent(i,j,k,0) = 0.0;
        if (std::abs(aux_vcent(i,j,k,1)) < small_value) aux_vcent(i,j,k,1) = 0.0;
        if (std::abs(aux_vcent(i,j,k,2)) < small_value) aux_vcent(i,j,k,2) = 0.0;
        if (std::abs(aux_bcent(i,j,k,0)) < small_value) aux_bcent(i,j,k,0) = 0.0;
        if (std::abs(aux_bcent(i,j,k,1)) < small_value) aux_bcent(i,j,k,1) = 0.0;
        if (std::abs(aux_bcent(i,j,k,2)) < small_value) aux_bcent(i,j,k,2) = 0.0;
      });

      // Area fraction MultiFab has one more slice at bigEnd(idim),
      // and this slice is not filled by fillBoundary(), for higher levels.
      // (Lower level might be filled by fillBoundary().)
      // Fill the ghost region for the last slice at bigEnd(idim)
      // by the value of the nearst point. And let fillBoundary() overwrite it.

      Box upper_slab = makeSlab(bx_grown, a_idim, bx.bigEnd(a_idim)+1);
      Box bx_grown_1 = bx; bx_grown_1.grow(a_idim,1);
      BoxList slab_diffList = boxDiff(upper_slab, bx_grown_1);

      for (const Box& b : slab_diffList) {
        ParallelFor(b, [=] AMREX_GPU_DEVICE (int i, int j, int k) noexcept
        {
          IntVect iv(AMREX_D_DECL(i,j,k));
          IntVect iv_nearest = iv;
          for (int d=0; d<AMREX_SPACEDIM; ++d) {
              iv_nearest[d] = Clamp(iv[d], bx_grown_1.smallEnd(d), bx_grown_1.bigEnd(d));
          }
          aux_afrac_x(iv) = aux_afrac_x(iv_nearest);
        });
      }

    } // if (FlagFab[mfi].getType(bx) == FabType::singlevalued )

  } // MFIter

  // Fill Boundary

  // The FB call for volfrac is done above
  // m_volfrac->FillBoundary(a_geom.periodicity());

  m_volcent->FillBoundary(a_geom.periodicity());
  for (int idim = 0; idim < AMREX_SPACEDIM; ++idim) {
    m_areafrac[idim]->FillBoundary(a_geom.periodicity());
    m_facecent[idim]->FillBoundary(a_geom.periodicity());
  }
  m_bndryarea->FillBoundary(a_geom.periodicity());
  m_bndrycent->FillBoundary(a_geom.periodicity());
  m_bndrynorm->FillBoundary(a_geom.periodicity());

  // Set Connectivities
  for (MFIter mfi(*m_cellflags, false); mfi.isValid(); ++mfi) {

    const Box& bx = mfi.validbox();
    const Box domain = surroundingNodes(a_geom.Domain(), a_idim);

    if (FlagFab[mfi].getType(bx) == FabType::singlevalued ) {

      Array4<EBCellFlag> const& aux_flag  = m_cellflags->array(mfi);
      Array4<Real>       const& aux_afrac_x = m_areafrac[0]->array(mfi);
      Array4<Real>       const& aux_afrac_y = m_areafrac[1]->array(mfi);
      Array4<Real>       const& aux_afrac_z = m_areafrac[2]->array(mfi);

      ParallelFor(bx, [=] AMREX_GPU_DEVICE (int i, int j, int k) noexcept
      {
        EB2::build_cellflag_from_ap (i, j, k, aux_flag, aux_afrac_x, aux_afrac_y, aux_afrac_z);
      });

      // Set disconnected non-periodicfaces

      bool l_periodic_x = a_geom.isPeriodic(0);
      bool l_periodic_y = a_geom.isPeriodic(1);
      bool l_periodic_z = a_geom.isPeriodic(2);

      if (!l_periodic_x) {
        const Box dom_grown = grow(grow(domain,1,1),2,1);
        const Box bx_grown  = grow(grow(    bx,1,1),2,1);
        const Box bx_face_x_lo = bx_grown & makeSlab(dom_grown,0,domain.smallEnd(0));
        const Box bx_face_x_hi = bx_grown & makeSlab(dom_grown,0,domain.bigEnd(0));

        ParallelFor(bx_face_x_lo, [=] AMREX_GPU_DEVICE (int i, int j, int k) noexcept
        {
          for(int kk(-1); kk<=1; kk++) {
          for(int jj(-1); jj<=1; jj++) {
            aux_flag(i,j,k).setDisconnected(-1,jj,kk);
          }}
        });
        ParallelFor(bx_face_x_hi, [=] AMREX_GPU_DEVICE (int i, int j, int k) noexcept
        {
          for(int kk(-1); kk<=1; kk++) {
          for(int jj(-1); jj<=1; jj++) {
            aux_flag(i,j,k).setDisconnected( 1,jj,kk);
          }}
        });
      }

      if (!l_periodic_y) {
        const Box dom_grown = grow(grow(domain,0,1),2,1);
        const Box bx_grown  = grow(grow(    bx,0,1),2,1);
        const Box bx_face_y_lo = bx_grown & makeSlab(dom_grown,1,domain.smallEnd(1));
        const Box bx_face_y_hi = bx_grown & makeSlab(dom_grown,1,domain.bigEnd(1));

        ParallelFor(bx_face_y_lo, [=] AMREX_GPU_DEVICE (int i, int j, int k) noexcept
        {
          for(int kk(-1); kk<=1; kk++) {
          for(int ii(-1); ii<=1; ii++) {
            aux_flag(i,j,k).setDisconnected(ii,-1,kk);
          }}
        });
        ParallelFor(bx_face_y_hi, [=] AMREX_GPU_DEVICE (int i, int j, int k) noexcept
        {
          for(int kk(-1); kk<=1; kk++) {
          for(int ii(-1); ii<=1; ii++) {
            aux_flag(i,j,k).setDisconnected(ii, 1,kk);
          }}
        });
      }

      if (!l_periodic_z) {
        const Box dom_grown = grow(grow(domain,0,1),1,1);
        const Box bx_grown  = grow(grow(    bx,0,1),1,1);
        const Box bx_face_z_lo = bx_grown & makeSlab(dom_grown,2,domain.smallEnd(2));
        const Box bx_face_z_hi = bx_grown & makeSlab(dom_grown,2,domain.bigEnd(2));

        ParallelFor(bx_face_z_lo, [=] AMREX_GPU_DEVICE (int i, int j, int k) noexcept
        {
          for(int jj(-1); jj<=1; jj++) {
          for(int ii(-1); ii<=1; ii++) {
            aux_flag(i,j,k).setDisconnected(ii,jj,-1);
          }}
        });
        ParallelFor(bx_face_z_hi, [=] AMREX_GPU_DEVICE (int i, int j, int k) noexcept
        {
          for(int jj(-1); jj<=1; jj++) {
          for(int ii(-1); ii<=1; ii++) {
            aux_flag(i,j,k).setDisconnected(ii,jj, 1);
          }}
        });
      }

    } // FabType::singlevalued

  } // MFIter

  // Set disconnected zero-volume-fraction cells
  // (equivalent to eb_::set_connection_flags for CC grids)

  for (MFIter mfi(*m_cellflags, false); mfi.isValid(); ++mfi) {

    const Box& bx = mfi.validbox();
    const Box gbx = amrex::grow(bx, m_cellflags->nGrow()-1); // Leave one cell layer

    Array4<EBCellFlag> const& aux_flag  = m_cellflags->array(mfi);
    Array4<Real>       const& aux_vfrac = m_volfrac->array(mfi);

    ParallelFor(gbx, [=] AMREX_GPU_DEVICE (int i, int j, int k) noexcept
    {
      for(int kk(-1); kk<=1; kk++) {
      for(int jj(-1); jj<=1; jj++) {
      for(int ii(-1); ii<=1; ii++)
      {
        if (aux_vfrac(i+ii,j+jj,k+kk) == 0.0) {
            aux_flag(i,j,k).setDisconnected(ii,jj,kk);
        }
      }}}
    });

    ParallelFor(gbx, [=] AMREX_GPU_DEVICE (int i, int j, int k) noexcept
    {
        if (aux_vfrac(i,j,k)==0.0) {
            aux_flag(i,j,k).setCovered();
        }
    });

  } // MFIter

  // Fill Boundary

  m_cellflags->FillBoundary(a_geom.periodicity());

#if 0
  // We leave these here for debugging if necessary.
  // If you uncomment these, make sure to uncomment AMReX_VisMF include above
  if (a_idim == 0) {
      amrex::VisMF::Write(*m_volfrac,"UVOL");
      amrex::VisMF::Write(*m_areafrac[0],"UAREAX");
      amrex::VisMF::Write(*m_areafrac[1],"UAREAY");
      amrex::VisMF::Write(*m_areafrac[2],"UAREAZ");
  } else if (a_idim == 1) {
      amrex::VisMF::Write(*m_volfrac,"VVOL");
      amrex::VisMF::Write(*m_areafrac[0],"VAREAX");
      amrex::VisMF::Write(*m_areafrac[1],"VAREAY");
      amrex::VisMF::Write(*m_areafrac[2],"VAREAZ");
  } else {
      amrex::VisMF::Write(*m_volfrac,"WVOL");
      amrex::VisMF::Write(*m_areafrac[0],"WAREAX");
      amrex::VisMF::Write(*m_areafrac[1],"WAREAY");
      amrex::VisMF::Write(*m_areafrac[2],"WAREAZ");
  }
#endif
}

const FabArray<EBCellFlagFab>&
eb_aux_::getMultiEBCellFlagFab () const
{
    AMREX_ASSERT(m_cellflags != nullptr);
    return *m_cellflags;
}

const MultiFab&
eb_aux_::getVolFrac () const
{
    AMREX_ASSERT(m_volfrac != nullptr);
    return *m_volfrac;
}

const MultiFab&
eb_aux_::getCentroid () const
{
    AMREX_ASSERT(m_volcent != nullptr);
    return *m_volcent;
}

const MultiFab&
eb_aux_::getBndryArea () const
{
    AMREX_ASSERT(m_bndryarea != nullptr);
    return *m_bndryarea;
}

const MultiFab&
eb_aux_::getBndryCent () const
{
    AMREX_ASSERT(m_bndrycent != nullptr);
    return *m_bndrycent;
}

const MultiFab&
eb_aux_::getBndryNorm () const
{
    AMREX_ASSERT(m_bndrynorm != nullptr);
    return *m_bndrynorm;
}

Array<const MultiFab*, AMREX_SPACEDIM>
eb_aux_::getAreaFrac () const
{
    AMREX_ASSERT(m_areafrac[0] != nullptr);
    return {AMREX_D_DECL(m_areafrac[0], m_areafrac[1], m_areafrac[2])};
}

Array<const MultiFab*, AMREX_SPACEDIM>
eb_aux_::getFaceCent () const
{
    AMREX_ASSERT(m_facecent[0] != nullptr);
    return {AMREX_D_DECL(m_facecent[0], m_facecent[1], m_facecent[2])};
}
