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

  m_cellflags->setVal(EBCellFlag::TheDefaultCell());

  m_volfrac = new MultiFab(grids, a_dmap, 1, a_ngrow[1], MFInfo(), FArrayBoxFactory());


#if 0
  m_centroid = new MultiCutFab(a_ba, a_dm, AMREX_SPACEDIM, m_ngrow[1], *m_cellflags);
  m_bndrycent = new MultiCutFab(a_ba, a_dm, AMREX_SPACEDIM, m_grow[2], *m_cellflags);

  m_bndryarea = new MultiCutFab(a_ba, a_dm, 1, m_grow[2], *m_cellflags);
  m_bndrynorm = new MultiCutFab(a_ba, a_dm, AMREX_SPACEDIM, m_grow[2], *m_cellflags);

  for (int idim = 0; idim < AMREX_SPACEDIM; ++idim) {
      const BoxArray& faceba = amrex::convert(a_ba, IntVect::TheDimensionVector(idim));
      m_areafrac[idim] = new MultiCutFab(faceba, a_dm, 1, m_grow[2], *m_cellflags);
      m_facecent[idim] = new MultiCutFab(faceba, a_dm, AMREX_SPACEDIM-1, m_grow[2], *m_cellflags);
  }
#endif

  const auto& FlagFab = a_factory->getMultiEBCellFlagFab();

  for (MFIter mfi(*m_cellflags, false); mfi.isValid(); ++mfi) {

    const Box& bx = mfi.validbox();

    if (FlagFab[mfi].getType(bx) == FabType::singlevalued ) {

      GpuArray<Real, AMREX_SPACEDIM> dx = a_geom.CellSizeArray();

      Array4<EBCellFlag const> const& flag = FlagFab.const_array(mfi);

      Array4<Real const> const& vfrac = (a_factory->getVolFrac()).const_array(mfi);
      // Array4<Real const> const& ccent = (a_factory->getCentroid()).const_array(mfi);

      Array4<Real const> const& afrac = (a_factory->getAreaFrac()[a_idim])->const_array(mfi);

      // EB normal and face centroid
      Array4<Real const> const& bnorm = a_factory->getBndryNormal()[mfi].const_array();
      Array4<Real const> const& bcent = a_factory->getBndryCent()[mfi].const_array();

      Array4<EBCellFlag> const& aux_flag  = m_cellflags->array(mfi);
      Array4<Real>       const& aux_vfrac = m_volfrac->array(mfi);

      bool is_per = a_geom.isPeriodic(a_idim);

      ParallelFor(bx, [
#ifndef AMREX_USE_GPU
                  verbose=m_verbose,
#endif
                  dx, bx, vfrac, afrac, bnorm, bcent, flag,
                  aux_flag, aux_vfrac, vdim, idim=a_idim, is_per ]
      AMREX_GPU_DEVICE (int i, int j, int k) noexcept
      {
        aux_flag(i,j,k).setCovered();
        aux_flag(i,j,k).setDisconnected();

        aux_vfrac(i,j,k) = 0.0;

        IntVect iv_hi(i,j,k);
        IntVect iv_lo(iv_hi - vdim);
        if (!is_per && iv_hi[idim]==bx.bigEnd(idim)){
          iv_hi = iv_lo; // At the upper boundary, hi cell takes the values of the low cell.
        }
        if (!is_per && iv_hi[idim]==bx.smallEnd(idim)){
          iv_lo = iv_hi; // At the lower boundary, low cell takes the values of the high cell.
        }

        if ( flag(iv_lo).isRegular() && flag(iv_hi).isRegular()) {

          aux_flag(i,j,k).setRegular();
          aux_flag(i,j,k).setConnected(vdim);

          aux_vfrac(i,j,k) = 1.0;

        } else if ( flag(iv_lo).isCovered() && flag(iv_hi).isCovered()) {

          // defaults to covered and disconnected.

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
          Real norm = ( (bnorm(iv_lo,0)*dx[0])*(bnorm(iv_lo,0)*dx[0])
                      + (bnorm(iv_lo,1)*dx[1])*(bnorm(iv_lo,1)*dx[1])
                      + (bnorm(iv_lo,2)*dx[2])*(bnorm(iv_lo,2)*dx[2]) );
          Real bcent_isoparam_x = bcent(iv_lo,0) / norm * dx[1] * dx[2];
          Real bcent_isoparam_y = bcent(iv_lo,1) / norm * dx[0] * dx[2];
          Real bcent_isoparam_z = bcent(iv_lo,2) / norm * dx[0] * dx[1];

          Real bnorm_x = bnorm(iv_lo,0) / dx[0];
          Real bnorm_y = bnorm(iv_lo,1) / dx[1];
          Real bnorm_z = bnorm(iv_lo,2) / dx[2];

          norm = sqrt( bnorm_x*bnorm_x + bnorm_y*bnorm_y + bnorm_z*bnorm_z);

          Real bnorm_isoparam_x = bnorm_x / norm;
          Real bnorm_isoparam_y = bnorm_y / norm;
          Real bnorm_isoparam_z = bnorm_z / norm;

          // plane point and normal
          RealVect lo_point (bcent_isoparam_x, bcent_isoparam_y, bcent_isoparam_z);
          RealVect lo_normal(bnorm_isoparam_x, bnorm_isoparam_y, bnorm_isoparam_z);

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

          norm = ( (bnorm(iv_hi,0)*dx[0])*(bnorm(iv_hi,0)*dx[0])
                 + (bnorm(iv_hi,1)*dx[1])*(bnorm(iv_hi,1)*dx[1])
                 + (bnorm(iv_hi,2)*dx[2])*(bnorm(iv_hi,2)*dx[2]) );
          bcent_isoparam_x = bcent(iv_hi,0) / norm * dx[1] * dx[2];
          bcent_isoparam_y = bcent(iv_hi,1) / norm * dx[0] * dx[2];
          bcent_isoparam_z = bcent(iv_hi,2) / norm * dx[0] * dx[1];

          bnorm_x = bnorm(iv_hi,0) / dx[0];
          bnorm_y = bnorm(iv_hi,1) / dx[1];
          bnorm_z = bnorm(iv_hi,2) / dx[2];

          norm = sqrt( bnorm_x*bnorm_x + bnorm_y*bnorm_y + bnorm_z*bnorm_z);

          bnorm_isoparam_x = bnorm_x / norm;
          bnorm_isoparam_y = bnorm_y / norm;
          bnorm_isoparam_z = bnorm_z / norm;

          // plane point and normal
          RealVect hi_point (bcent_isoparam_x, bcent_isoparam_y, bcent_isoparam_z);
          RealVect hi_normal(bnorm_isoparam_x, bnorm_isoparam_y, bnorm_isoparam_z);

          // Low side of high cell
          lo_arr[idim] = -0.5;
          hi_arr[idim] =  0.0;
          RealBox hi_rbx(lo_arr.data(), hi_arr.data());

          eb_cut_cell_ hi_eb_cc(flag(iv_hi), hi_rbx, hi_point, hi_normal);

          // cell iv_hi covered (regular) imples hi_eb_cc is covered (regular)
          // The inverse is not always true.
          AMREX_ASSERT( !flag(iv_hi).isCovered() || hi_eb_cc.isCovered() );
          AMREX_ASSERT( !flag(iv_hi).isRegular() || hi_eb_cc.isRegular() );

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
            { Real const abs_err = amrex::min(std::abs(lo_eb_cc.areaHi(idim) - afrac(iv_hi)),
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
          }
#endif

          if (lo_eb_cc.isCovered() && hi_eb_cc.isCovered()) {

            // defaults to covered and disconnected.

          } else if (lo_eb_cc.isRegular() && hi_eb_cc.isRegular()) {

            aux_flag(i,j,k).setRegular();
            aux_flag(i,j,k).setConnected(vdim);

            aux_vfrac(i,j,k) = 1.0;

          } else {

            if (lo_eb_cc.isCovered()) { }
            if (hi_eb_cc.isCovered()) { }

          }
        }
      });

    }

  }
}
