#include <ERF_EBAux.H>
#include <ERF_EBCutCell.H>

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
define( int const& a_idim,
        Geometry            const& a_geom,
        BoxArray            const& a_grids,
        DistributionMapping const& a_dmap,
        Vector<int>         const& a_ngrow,
        EBFArrayBoxFactory  const* a_factory)
{
  // Box dbox(a_geom.Domain());

  const IntVect vdim(IntVect::TheDimensionVector(a_idim));

  const BoxArray& grids = amrex::convert(a_grids, vdim);

  m_cellflags = new FabArray<EBCellFlagFab>(grids, a_dmap, 1, a_ngrow[0], MFInfo(),
                                            DefaultFabFactory<EBCellFlagFab>());

  // Set m_cellflags type to singlevalued
  m_cellflags->setVal(EBCellFlag::TheDefaultCell());
  for (MFIter mfi(*m_cellflags, false); mfi.isValid(); ++mfi) {
    auto& fab = (*m_cellflags)[mfi];
    fab.setType(FabType::singlevalued);
  }

  m_volfrac = new MultiFab(grids, a_dmap, 1, a_ngrow[1], MFInfo(), FArrayBoxFactory());
  m_volcent = new MultiCutFab(grids, a_dmap, AMREX_SPACEDIM, a_ngrow[2], *m_cellflags);

  for (int idim = 0; idim < AMREX_SPACEDIM; ++idim) {
      const BoxArray& faceba = amrex::convert(a_grids, IntVect::TheDimensionVector(idim));
      m_areafrac[idim] = new MultiCutFab(faceba, a_dmap, 1, a_ngrow[2], *m_cellflags);
      m_facecent[idim] = new MultiCutFab(faceba, a_dmap, AMREX_SPACEDIM-1, a_ngrow[2], *m_cellflags);
  }

  m_bndryarea = new MultiCutFab(grids, a_dmap, 1, a_ngrow[2], *m_cellflags);
  m_bndrycent = new MultiCutFab(grids, a_dmap, AMREX_SPACEDIM, a_ngrow[2], *m_cellflags);
  m_bndrynorm = new MultiCutFab(grids, a_dmap, AMREX_SPACEDIM, a_ngrow[2], *m_cellflags);

  const auto& FlagFab = a_factory->getMultiEBCellFlagFab(); // EBFArrayBoxFactory, EBDataCollection

  for (MFIter mfi(*m_cellflags, false); mfi.isValid(); ++mfi) {

    const Box& bx = mfi.validbox();

    if (FlagFab[mfi].getType(bx) == FabType::singlevalued ) {

      GpuArray<Real, AMREX_SPACEDIM> dx = a_geom.CellSizeArray();

      Array4<EBCellFlag const> const& flag = FlagFab.const_array(mfi);

      // Array4<Real const> const& vfrac = (a_factory->getVolFrac()).const_array(mfi);
      // Array4<Real const> const& ccent = (a_factory->getCentroid()).const_array(mfi);
      // Array4<Real const> const& afrac = (a_factory->getAreaFrac()[a_idim])->const_array(mfi);

      // EB normal and face centroid
      Array4<Real const> const& bnorm = a_factory->getBndryNormal()[mfi].const_array();
      Array4<Real const> const& bcent = a_factory->getBndryCent()[mfi].const_array();

      // aux quantities
      Array4<EBCellFlag> const& aux_flag  = m_cellflags->array(mfi);
      Array4<Real>       const& aux_vfrac = m_volfrac->array(mfi);
      Array4<Real>       const& aux_vcent = m_volcent->array(mfi);

      Array4<Real>       const& aux_afrac_x = m_areafrac[0]->array(mfi);
      Array4<Real>       const& aux_afrac_y = m_areafrac[1]->array(mfi);
      Array4<Real>       const& aux_afrac_z = m_areafrac[2]->array(mfi);

      Array4<Real>       const& aux_fcent_x = m_facecent[0]->array(mfi);
      Array4<Real>       const& aux_fcent_y = m_facecent[1]->array(mfi);
      Array4<Real>       const& aux_fcent_z = m_facecent[2]->array(mfi);

      Array4<Real>       const& aux_barea = m_bndryarea->array(mfi);
      Array4<Real>       const& aux_bcent = m_bndrycent->array(mfi);
      Array4<Real>       const& aux_bnorm = m_bndrynorm->array(mfi);

      bool is_per = a_geom.isPeriodic(a_idim);

      ParallelFor(bx, [
#ifndef AMREX_USE_GPU
                  verbose=m_verbose,
#endif
                  dx, bx, bnorm, bcent, flag,
                  aux_flag, aux_vfrac, aux_vcent,
                  aux_afrac_x, aux_afrac_y, aux_afrac_z,
                  aux_fcent_x, aux_fcent_y, aux_fcent_z,
                  aux_barea, aux_bcent, aux_bnorm,
                  vdim, idim=a_idim, is_per ]
      AMREX_GPU_DEVICE (int i, int j, int k) noexcept
      {

        // defaults to covered and disconnected.

        aux_flag(i,j,k).setCovered();
        aux_flag(i,j,k).setDisconnected();

        aux_vfrac(i,j,k) = 0.0;
        aux_vcent(i,j,k,0) = 0.0;
        aux_vcent(i,j,k,1) = 0.0;
        aux_vcent(i,j,k,2) = 0.0;

        aux_afrac_x(i,j,k) = 0.0;
        aux_afrac_y(i,j,k) = 0.0;
        aux_afrac_z(i,j,k) = 0.0;

        aux_fcent_x(i,j,k,0) = 0.0; aux_fcent_x(i,j,k,1) = 0.0;
        aux_fcent_y(i,j,k,0) = 0.0; aux_fcent_y(i,j,k,1) = 0.0;
        aux_fcent_z(i,j,k,0) = 0.0; aux_fcent_z(i,j,k,1) = 0.0;

        if (i==bx.bigEnd(0)) {
          aux_afrac_x(i+1,j,k) = 0.0;
          aux_fcent_x(i+1,j,k,0) = 0.0; aux_fcent_x(i+1,j,k,1) = 0.0;
        }
        if (j==bx.bigEnd(1)) {
          aux_afrac_y(i,j+1,k) = 0.0;
          aux_fcent_y(i,j+1,k,0) = 0.0; aux_fcent_y(i,j+1,k,1) = 0.0;
        }
        if (k==bx.bigEnd(2)) {
          aux_afrac_z(i,j,k+1) = 0.0;
          aux_fcent_z(i,j,k+1,0) = 0.0; aux_fcent_z(i,j,k+1,1) = 0.0;
        }

        aux_barea(i,j,k) = 0.0;

        aux_bcent(i,j,k,0) = 0.0;
        aux_bcent(i,j,k,1) = 0.0;
        aux_bcent(i,j,k,2) = 0.0;

        aux_bnorm(i,j,k,0) = 0.0;
        aux_bnorm(i,j,k,1) = 0.0;
        aux_bnorm(i,j,k,2) = 0.0;

        // Index for low and hi cells
        IntVect iv_hi(i,j,k);
        IntVect iv_lo(iv_hi - vdim);
        if (!is_per && iv_hi[idim]==bx.bigEnd(idim)){
          iv_hi = iv_lo; // At the upper boundary, hi cell takes the values of the low cell.
        }
        if (!is_per && iv_hi[idim]==bx.smallEnd(idim)){
          iv_lo = iv_hi; // At the lower boundary, low cell takes the values of the high cell.
        }

        //

        if ( flag(iv_lo).isCovered() && flag(iv_hi).isCovered()) {

          // defaults to covered and disconnected.

        } else if ( flag(iv_lo).isRegular() && flag(iv_hi).isRegular()) {

          aux_flag(i,j,k).setRegular();
          aux_flag(i,j,k).setConnected(vdim);

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

          if (!is_per && iv_hi[idim]==bx.smallEnd(idim)){
            lo_point[idim] += 1.0; // Move the boundary centroid upward in the idim direction.
          }

          if (flag(iv_lo).isSingleValued() ) {

            Real bnorm_x = bnorm(iv_lo,0) * dx[0];
            Real bnorm_y = bnorm(iv_lo,1) * dx[1];
            Real bnorm_z = bnorm(iv_lo,2) * dx[2];

            Real norm = sqrt( bnorm_x*bnorm_x + bnorm_y*bnorm_y + bnorm_z*bnorm_z);

            RealVect bnorm_isoparam ( bnorm_x / norm, bnorm_y / norm, bnorm_z / norm);

            // plane point and normal
            // lo_point  = bcent_isoparam;
            lo_normal = bnorm_isoparam;

          }

          // High side of low cell
          lo_arr[idim] = 0.0;
          hi_arr[idim] = 0.5;
          RealBox lo_rbx(lo_arr.data(), hi_arr.data());

          eb_cut_cell_ lo_eb_cc(flag(iv_lo), lo_rbx, lo_point, lo_normal);

          // cell iv_lo covered (regular) imples lo_eb_cc is covered (regular)
          // The inverse is not always true.
          AMREX_ASSERT( !flag(iv_lo).isCovered() || lo_eb_cc.isCovered() );
          AMREX_ASSERT( !flag(iv_lo).isRegular() || lo_eb_cc.isRegular() );

          //-----------------------
          // High EB cut cell
          //-----------------------

          RealVect hi_point (bcent(iv_hi,0), bcent(iv_hi,1), bcent(iv_hi,2));
          RealVect hi_normal(bnorm(iv_hi,0), bnorm(iv_hi,1), bnorm(iv_hi,2));

          if (!is_per && iv_hi[idim]==bx.bigEnd(idim)){
            lo_point[idim] += -1.0; // Move the boundary centroid downward in the idim direction.
          }

          if (flag(iv_hi).isSingleValued() ) {

            Real bnorm_x = bnorm(iv_hi,0) * dx[0];
            Real bnorm_y = bnorm(iv_hi,1) * dx[1];
            Real bnorm_z = bnorm(iv_hi,2) * dx[2];

            Real norm = sqrt( bnorm_x*bnorm_x + bnorm_y*bnorm_y + bnorm_z*bnorm_z);

            RealVect bnorm_isoparam ( bnorm_x / norm, bnorm_y / norm, bnorm_z / norm);

            // plane point and normal
            // hi_point  = bcent_isoparam;
            hi_normal = bnorm_isoparam;

          }

          // Low side of high cell
          lo_arr[idim] = -0.5;
          hi_arr[idim] =  0.0;
          RealBox hi_rbx(lo_arr.data(), hi_arr.data());

          eb_cut_cell_ hi_eb_cc(flag(iv_hi), hi_rbx, hi_point, hi_normal);

          // cell iv_hi covered (regular) imples hi_eb_cc is covered (regular)
          // The inverse is not always true.
          AMREX_ASSERT( !flag(iv_hi).isCovered() || hi_eb_cc.isCovered() );
          AMREX_ASSERT( !flag(iv_hi).isRegular() || hi_eb_cc.isRegular() );

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
            if ( !(!flag(iv_hi).isRegular() || hi_hi_eb_cc.isRegular()) ||
                 !(!flag(iv_hi).isCovered() || hi_hi_eb_cc.isCovered()) ) {
              Print() << "flag(iv_hi) and hi_hi_eb_cc flags do not agree\n"
                      << "\n  isRegular() " << flag(iv_hi).isRegular() << "  " << hi_hi_eb_cc.isRegular()
                      << "\n  isCovered() " << flag(iv_hi).isCovered() << "  " << hi_hi_eb_cc.isCovered()
                      << "\n";
            }
#endif
            // If cell iv_hi is regular or covered, then hi_hi_eb_cc must also
            // be regular or covered. The inverse is not true.
            AMREX_ALWAYS_ASSERT( !flag(iv_hi).isRegular() || hi_hi_eb_cc.isRegular() );
            AMREX_ALWAYS_ASSERT( !flag(iv_hi).isCovered() || hi_hi_eb_cc.isCovered() );

            // The area and volume fractions that are computed for the scalar grid
            // are slightly different than those we compute from the geometric
            // reconstruction using the EB point and normal. However, we expect
            // that the area fractions computed here will give back the same
            // normal we used to compute them.
            if ( flag(iv_hi).isSingleValued() ) {

              Real const adx = (idim == 0)
                             ? (hi_eb_cc.areaLo(0) - hi_hi_eb_cc.areaHi(0)) * dx[1] * dx[2]
                             : (hi_eb_cc.areaLo(0) + hi_hi_eb_cc.areaLo(0)) * dx[1] * dx[2]
                             - (hi_eb_cc.areaHi(0) + hi_hi_eb_cc.areaHi(0)) * dx[1] * dx[2];

              Real const ady = (idim == 1)
                             ? (hi_eb_cc.areaLo(1) - hi_hi_eb_cc.areaHi(1)) * dx[0] * dx[2]
                             : (hi_eb_cc.areaLo(1) + hi_hi_eb_cc.areaLo(1)) * dx[0] * dx[2]
                             - (hi_eb_cc.areaHi(1) + hi_hi_eb_cc.areaHi(1)) * dx[0] * dx[2];

              Real const adz = (idim == 2)
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

            // The idim area of hi_eb_cc.areaHi() should equal hi_hi_eb_cc.areaLo()
            {
#ifndef AMREX_USE_GPU
            Real const abs_err = std::abs( hi_eb_cc.areaHi(idim) - hi_hi_eb_cc.areaLo(idim) );
            Real machine_tol = 10.0*std::numeric_limits<amrex::Real>::epsilon();
            if ( abs_err >= machine_tol ) {
                Print() << "\nFail: check-2 area abs_err: " << abs_err
                        << "\n  hi_eb_cc.areaHi " << hi_eb_cc.areaHi(idim)
                        << "\n  hi_hi_eb_cc.areaLo " << hi_hi_eb_cc.areaLo(idim)
                        << '\n';
            } else if (verbose) {
                Print() << "Pass: hi_eb_cc.areaHi = hi_hi_eb_cc.areaLo"
                        << "  abs_err: " << abs_err << "\n";
            }
            AMREX_ALWAYS_ASSERT( abs_err < machine_tol );
#endif
            }

            // The low-side area of hi_eb_cc should equal idim afrac.
            { Real const abs_err = amrex::max(std::abs(lo_eb_cc.areaHi(idim) - afrac(iv_hi)),
                                              std::abs(hi_eb_cc.areaLo(idim) - afrac(iv_hi)));
              Real compare_tol = 5.0e-6;
#ifndef AMREX_USE_GPU
              if ( abs_err >= compare_tol ) {
                //hi_eb_cc.debug();
                Print() << "\nFail: check-3 area abs_err " << abs_err
                        << "\n  hi_eb_cc.areaLo(" << idim << ") = " << hi_eb_cc.areaLo(idim)
                        << "\n  lo_eb_cc.areaHi(" << idim << ") = " << lo_eb_cc.areaHi(idim)
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
            aux_flag(i,j,k).setConnected(vdim);

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

            aux_flag(i,j,k).setSingleValued();
            aux_flag(i,j,k).setConnected(vdim);

            Real lo_vol {lo_eb_cc.volume()};
            Real hi_vol {hi_eb_cc.volume()};

            aux_vfrac(i,j,k) = lo_vol + hi_vol;

            /* centVol() returns the coordinates based on m_rbx.
              The coordinates in the idim direction are in [0.0,0.5] for the low cell and in [-0.5,0.0] for the hi cell.
              Therefore, they need to be mapped to the eb_aux space, by shifting:
              x' = x - 0.5 (low cell), x + 0.5 (hi cell) if idim = 0
              y' = y - 0.5 (low cell), y + 0.5 (hi cell) if idim = 1
              z' = z - 0.5 (low cell), z + 0.5 (hi cell) if idim = 2
            */

            RealVect lo_vcent {lo_eb_cc.centVol()};
            RealVect hi_vcent {hi_eb_cc.centVol()};

            lo_vcent[idim] = lo_vcent[idim] - 0.5;
            hi_vcent[idim] = hi_vcent[idim] + 0.5;

            aux_vcent(i,j,k,0) = ( lo_vol * lo_vcent[0] + hi_vol * hi_vcent[0] ) / aux_vfrac(i,j,k);
            aux_vcent(i,j,k,1) = ( lo_vol * lo_vcent[1] + hi_vol * hi_vcent[1] ) / aux_vfrac(i,j,k);
            aux_vcent(i,j,k,2) = ( lo_vol * lo_vcent[2] + hi_vol * hi_vcent[2] ) / aux_vfrac(i,j,k);

            Real lo_areaLo_x {lo_eb_cc.areaLo(0)};
            Real lo_areaLo_y {lo_eb_cc.areaLo(1)};
            Real lo_areaLo_z {lo_eb_cc.areaLo(2)};

            Real hi_areaLo_x {hi_eb_cc.areaLo(0)};
            Real hi_areaLo_y {hi_eb_cc.areaLo(1)};
            Real hi_areaLo_z {hi_eb_cc.areaLo(2)};

            aux_afrac_x(i,j,k) = (idim == 0) ? lo_areaLo_x : lo_areaLo_x + hi_areaLo_x;
            aux_afrac_y(i,j,k) = (idim == 1) ? lo_areaLo_y : lo_areaLo_y + hi_areaLo_y;
            aux_afrac_z(i,j,k) = (idim == 2) ? lo_areaLo_z : lo_areaLo_z + hi_areaLo_z;

            /* fcentLo returns the coordinates based on m_rbx.
              The coordinates in the idim direction are in [0.0,0.5] for the low cell and in [-0.5,0.0] for the hi cell.
              Therefore, they need to be mapped to the eb_aux space, by shifting:
              x' = x - 0.5 (low cell), x + 0.5 (hi cell) if idim = 0
              y' = y - 0.5 (low cell), y + 0.5 (hi cell) if idim = 1
              z' = z - 0.5 (low cell), z + 0.5 (hi cell) if idim = 2
            */

            RealVect lo_centLo_x {lo_eb_cc.centLo(0)};
            RealVect lo_centLo_y {lo_eb_cc.centLo(1)};
            RealVect lo_centLo_z {lo_eb_cc.centLo(2)};

            RealVect hi_centLo_x {hi_eb_cc.centLo(0)};
            RealVect hi_centLo_y {hi_eb_cc.centLo(1)};
            RealVect hi_centLo_z {hi_eb_cc.centLo(2)};

            if (idim == 0) {
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
            } else if (idim == 1) {
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
            } else if (idim == 2) {
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

            // Need to fill the nodes the big ends?

            Real lo_areaBoun {lo_eb_cc.areaBoun()};
            Real hi_areaBoun {hi_eb_cc.areaBoun()};

            aux_barea(i,j,k) = lo_areaBoun + hi_areaBoun;

            RealVect lo_centBoun {lo_eb_cc.centBoun()};
            RealVect hi_centBoun {hi_eb_cc.centBoun()};

            if (idim == 0) {
              aux_bcent(i,j,k,0) = ( lo_areaBoun * (lo_centBoun[0]-0.5) + hi_areaBoun * (hi_centBoun[0]+0.5) ) / aux_barea(i,j,k);  // x (mapped)
              aux_bcent(i,j,k,1) = ( lo_areaBoun * lo_centBoun[1] + hi_areaBoun * hi_centBoun[1] ) / aux_barea(i,j,k);              // y
              aux_bcent(i,j,k,2) = ( lo_areaBoun * lo_centBoun[2] + hi_areaBoun * hi_centBoun[2] ) / aux_barea(i,j,k);              // z
            } else if (idim == 1) {
              aux_bcent(i,j,k,0) = ( lo_areaBoun * lo_centBoun[0] + hi_areaBoun * hi_centBoun[0] ) / aux_barea(i,j,k);              // x
              aux_bcent(i,j,k,1) = ( lo_areaBoun * (lo_centBoun[1]-0.5) + hi_areaBoun * (hi_centBoun[1]+0.5) ) / aux_barea(i,j,k);  // y (mapped)
              aux_bcent(i,j,k,2) = ( lo_areaBoun * lo_centBoun[2] + hi_areaBoun * hi_centBoun[2] ) / aux_barea(i,j,k);              // z
            } else if (idim == 2) {
              aux_bcent(i,j,k,0) = ( lo_areaBoun * lo_centBoun[0] + hi_areaBoun * hi_centBoun[0] ) / aux_barea(i,j,k);              // x
              aux_bcent(i,j,k,1) = ( lo_areaBoun * lo_centBoun[1] + hi_areaBoun * hi_centBoun[1] ) / aux_barea(i,j,k);              // y
              aux_bcent(i,j,k,2) = ( lo_areaBoun * (lo_centBoun[2]-0.5) + hi_areaBoun * (hi_centBoun[2]+0.5) ) / aux_barea(i,j,k);  // z (mapped)
            }

            RealVect eb_normal = ( lo_areaBoun * lo_normal + hi_areaBoun * hi_normal )/ aux_barea(i,j,k);

            aux_bnorm(i,j,k,0) = eb_normal[0];
            aux_bnorm(i,j,k,1) = eb_normal[1];
            aux_bnorm(i,j,k,2) = eb_normal[2];
          }

        } // flag(iv_lo) and flag(iv_hi)

      });

    }

  }
}

const MultiFab&
eb_aux_::getVolFrac () const
{
    AMREX_ASSERT(m_volfrac != nullptr);
    return *m_volfrac;
}

const MultiCutFab&
eb_aux_::getVolCent () const
{
    AMREX_ASSERT(m_volcent != nullptr);
    return *m_volcent;
}

const MultiCutFab&
eb_aux_::getBndryArea () const
{
    AMREX_ASSERT(m_bndryarea != nullptr);
    return *m_bndryarea;
}

const MultiCutFab&
eb_aux_::getBndryCent () const
{
    AMREX_ASSERT(m_bndrycent != nullptr);
    return *m_bndrycent;
}

const MultiCutFab&
eb_aux_::getBndryNorm () const
{
    AMREX_ASSERT(m_bndrynorm != nullptr);
    return *m_bndrynorm;
}

Array<const MultiCutFab*, AMREX_SPACEDIM>
eb_aux_::getAreaFrac () const
{
    AMREX_ASSERT(m_areafrac[0] != nullptr);
    return {AMREX_D_DECL(m_areafrac[0], m_areafrac[1], m_areafrac[2])};
}

Array<const MultiCutFab*, AMREX_SPACEDIM>
eb_aux_::getFaceCent () const
{
    AMREX_ASSERT(m_facecent[0] != nullptr);
    return {AMREX_D_DECL(m_facecent[0], m_facecent[1], m_facecent[2])};
}
