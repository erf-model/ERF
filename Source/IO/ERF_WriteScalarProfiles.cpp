#include <iomanip>

#include "ERF.H"

using namespace amrex;

/**
 * Computes the integrated quantities on the grid such as the
 * total scalar and total mass quantities. Prints and writes to output file.
 *
 * @param time Current time
 */
void
ERF::sum_integrated_quantities (Real time)
{
    BL_PROFILE("ERF::sum_integrated_quantities()");

    if (verbose <= 0)
      return;

    int datwidth = 14;
    int datprecision = 6;

    // Single level sum
    Real mass_sl;

    // Multilevel sums
    Real mass_ml = 0.0;
    Real rhth_ml = 0.0;
    Real scal_ml = 0.0;

#if 1
    mass_sl = volWgtSumMF(0,vars_new[0][Vars::cons],Rho_comp,*mapfac_m[0],false,false);
    for (int lev = 0; lev <= finest_level; lev++) {
        mass_ml += volWgtSumMF(lev,vars_new[lev][Vars::cons],Rho_comp,*mapfac_m[lev],false,true);
    }
#else
    for (int lev = 0; lev <= finest_level; lev++) {
        MultiFab pert_dens(vars_new[lev][Vars::cons].boxArray(),
                           vars_new[lev][Vars::cons].DistributionMap(),
                           1,0);
        MultiFab r_hse (base_state[lev], make_alias, 0, 1); // r_0 is first  component
        for ( MFIter mfi(pert_dens,TilingIfNotGPU()); mfi.isValid(); ++mfi)
        {
            const Box& bx = mfi.tilebox();
            const Array4<Real      >& pert_dens_arr = pert_dens.array(mfi);
            const Array4<Real const>&         S_arr = vars_new[lev][Vars::cons].const_array(mfi);
            const Array4<Real const>&        r0_arr = r_hse.const_array(mfi);
            ParallelFor(bx, [=] AMREX_GPU_DEVICE(int i, int j, int k) noexcept {
                pert_dens_arr(i, j, k, 0) = S_arr(i,j,k,Rho_comp) - r0_arr(i,j,k);
            });
        }
        if (lev == 0) {
            mass_sl = volWgtSumMF(0,pert_dens,0,*mapfac_m[0],false,false);
        }
        mass_ml += volWgtSumMF(lev,pert_dens,0,*mapfac_m[lev],false,true);
    } // lev
#endif

    Real rhth_sl = volWgtSumMF(0,vars_new[0][Vars::cons], RhoTheta_comp,*mapfac_m[0],false,false);
    Real scal_sl = volWgtSumMF(0,vars_new[0][Vars::cons],RhoScalar_comp,*mapfac_m[0],false,false);

    for (int lev = 0; lev <= finest_level; lev++) {
        rhth_ml += volWgtSumMF(lev,vars_new[lev][Vars::cons], RhoTheta_comp,*mapfac_m[lev],false,true);
        scal_ml += volWgtSumMF(lev,vars_new[lev][Vars::cons],RhoScalar_comp,*mapfac_m[lev],false,true);
    }

    if (verbose > 0) {

        Gpu::HostVector<Real> h_avg_ustar; h_avg_ustar.resize(1);
        Gpu::HostVector<Real> h_avg_tstar; h_avg_tstar.resize(1);
        Gpu::HostVector<Real> h_avg_olen; h_avg_olen.resize(1);
        if ((m_most != nullptr) && (NumDataLogs() > 0)) {
            Box domain = geom[0].Domain();
            int zdir = 2;
            h_avg_ustar = sumToLine(*m_most->get_u_star(0),0,1,domain,zdir);
            h_avg_tstar = sumToLine(*m_most->get_t_star(0),0,1,domain,zdir);
            h_avg_olen  = sumToLine(*m_most->get_olen(0),0,1,domain,zdir);

            // Divide by the total number of cells we are averaging over
            Real area_z = static_cast<Real>(domain.length(0)*domain.length(1));
            h_avg_ustar[0] /= area_z;
            h_avg_tstar[0] /= area_z;
            h_avg_olen[0]  /= area_z;

        } else {
            h_avg_ustar[0] = 0.;
            h_avg_tstar[0] = 0.;
            h_avg_olen[0]  = 0.;
        }

        const int nfoo = 6;
        Real foo[nfoo] = {mass_sl,rhth_sl,scal_sl,mass_ml,rhth_ml,scal_ml};
#ifdef AMREX_LAZY
        Lazy::QueueReduction([=]() mutable {
#endif
        ParallelDescriptor::ReduceRealSum(
            foo, nfoo, ParallelDescriptor::IOProcessorNumber());

          if (ParallelDescriptor::IOProcessor()) {
            int i = 0;
            mass_sl = foo[i++];
            rhth_sl = foo[i++];
            scal_sl = foo[i++];
            mass_ml = foo[i++];
            rhth_ml = foo[i++];
            scal_ml = foo[i++];

            Print() << '\n';
            if (finest_level ==  0) {
#if 1
               Print() << "TIME= " << time << "      MASS         = " << mass_sl << '\n';
#else
               Print() << "TIME= " << time << " PERT MASS         = " << mass_sl << '\n';
#endif
               Print() << "TIME= " << time << " RHO THETA         = " << rhth_sl << '\n';
               Print() << "TIME= " << time << " RHO SCALAR        = " << scal_sl << '\n';
            } else {
#if 1
               Print() << "TIME= " << time << "      MASS   SL/ML = " << mass_sl << " " << mass_ml << '\n';
#else
               Print() << "TIME= " << time << " PERT MASS   SL/ML = " << mass_sl << " " << mass_ml << '\n';
#endif
               Print() << "TIME= " << time << " RHO THETA   SL/ML = " << rhth_sl << " " << rhth_ml << '\n';
               Print() << "TIME= " << time << " RHO SCALAR  SL/ML = " << scal_sl << " " << scal_ml << '\n';
            }

            // The first data log only holds scalars
            if (NumDataLogs() > 0)
            {
                int nd = 0;
                std::ostream& data_log1 = DataLog(nd);
                if (data_log1.good()) {
                    if (time == 0.0) {
                        data_log1 << std::setw(datwidth) << "          time";
                        data_log1 << std::setw(datwidth) << "          u_star";
                        data_log1 << std::setw(datwidth) << "          t_star";
                        data_log1 << std::setw(datwidth) << "          olen";
                        data_log1 << std::endl;
                    } // time = 0

                  // Write the quantities at this time
                  data_log1 << std::setw(datwidth) << time;
                  data_log1 << std::setw(datwidth) << std::setprecision(datprecision)
                            << h_avg_ustar[0];
                  data_log1 << std::setw(datwidth) << std::setprecision(datprecision)
                            << h_avg_tstar[0];
                  data_log1 << std::setw(datwidth) << std::setprecision(datprecision)
                            << h_avg_olen[0];
                  data_log1 << std::endl;
                } // if good
            } // loop over i
          } // if IOProcessor
#ifdef AMREX_LAZY
        });
#endif
    } // if verbose

    // This is just an alias for convenience
    int lev = 0;
    if (NumSamplePointLogs() > 0 && NumSamplePoints() > 0) {
        for (int i = 0; i < NumSamplePoints(); ++i)
        {
            sample_points(lev, time, SamplePoint(i), vars_new[lev][Vars::cons]);
        }
    }
    if (NumSampleLineLogs() > 0 && NumSampleLines() > 0) {
        for (int i = 0; i < NumSampleLines(); ++i)
        {
            sample_lines(lev, time, SampleLine(i), vars_new[lev][Vars::cons]);
        }
    }
}

/**
 * Utility function for sampling MultiFab data at a specified cell index.
 *
 * @param lev Level for the associated MultiFab data
 * @param time Current time
 * @param cell IntVect containing the indexes for the cell where we want to sample
 * @param mf MultiFab from which we wish to sample data
 */
void
ERF::sample_points (int /*lev*/, Real time, IntVect cell, MultiFab& mf)
{
    int datwidth = 14;

    int ifile = 0;

    //
    // Sample the data at a single point in space
    //
    int ncomp = mf.nComp();
    Vector<Real> my_point = get_cell_data(mf, cell);

    if (!my_point.empty()) {

        // HERE DO WHATEVER YOU WANT TO THE DATA BEFORE WRITING

        std::ostream& sample_log = SamplePointLog(ifile);
        if (sample_log.good()) {
          sample_log << std::setw(datwidth) << time;
          for (int i = 0; i < ncomp; ++i)
          {
              sample_log << std::setw(datwidth) << my_point[i];
          }
          sample_log << std::endl;
        } // if good
    } // only write from processor that holds the cell
}

/**
 * Utility function for sampling data along a line along the z-dimension
 * at the (x,y) indices specified and writes it to an output file.
 *
 * @param lev Current level
 * @param time Current time
 * @param cell IntVect containing the x,y-dimension indices to sample along z
 * @param mf MultiFab from which we sample the data
 */
void
ERF::sample_lines (int lev, Real time, IntVect cell, MultiFab& mf)
{
    int datwidth = 14;
    int datprecision = 6;

    int ifile = 0;

    const int ncomp = mf.nComp(); // cell-centered state vars

    MultiFab mf_vels(grids[lev], dmap[lev], AMREX_SPACEDIM, 0);
    average_face_to_cellcenter(mf_vels, 0,
                               Array<const MultiFab*,3>{&vars_new[lev][Vars::xvel],&vars_new[lev][Vars::yvel],&vars_new[lev][Vars::zvel]});

    //
    // Sample the data at a line (in direction "dir") in space
    // In this case we sample in the vertical direction so dir = 2
    // The "k" value of "cell" is ignored
    //
    int dir = 2;
    MultiFab my_line       = get_line_data(mf,              dir, cell);
    MultiFab my_line_vels  = get_line_data(mf_vels,         dir, cell);
    MultiFab my_line_tau11 = get_line_data(*Tau11_lev[lev], dir, cell);
    MultiFab my_line_tau12 = get_line_data(*Tau12_lev[lev], dir, cell);
    MultiFab my_line_tau13 = get_line_data(*Tau13_lev[lev], dir, cell);
    MultiFab my_line_tau22 = get_line_data(*Tau22_lev[lev], dir, cell);
    MultiFab my_line_tau23 = get_line_data(*Tau23_lev[lev], dir, cell);
    MultiFab my_line_tau33 = get_line_data(*Tau33_lev[lev], dir, cell);

    for (MFIter mfi(my_line, false); mfi.isValid(); ++mfi)
    {
        // HERE DO WHATEVER YOU WANT TO THE DATA BEFORE WRITING

        std::ostream& sample_log = SampleLineLog(ifile);
        if (sample_log.good()) {
          sample_log << std::setw(datwidth) << std::setprecision(datprecision) << time;
          const auto& my_line_arr = my_line[0].const_array();
          const auto& my_line_vels_arr = my_line_vels[0].const_array();
          const auto& my_line_tau11_arr = my_line_tau11[0].const_array();
          const auto& my_line_tau12_arr = my_line_tau12[0].const_array();
          const auto& my_line_tau13_arr = my_line_tau13[0].const_array();
          const auto& my_line_tau22_arr = my_line_tau22[0].const_array();
          const auto& my_line_tau23_arr = my_line_tau23[0].const_array();
          const auto& my_line_tau33_arr = my_line_tau33[0].const_array();
          const Box&  my_box = my_line[0].box();
          const int klo = my_box.smallEnd(2);
          const int khi = my_box.bigEnd(2);
          int i = cell[0];
          int j = cell[1];
          for (int n = 0; n < ncomp; n++) {
              for (int k = klo; k <= khi; k++) {
                  sample_log << std::setw(datwidth) << std::setprecision(datprecision) << my_line_arr(i,j,k,n);
              }
          }
          for (int n = 0; n < AMREX_SPACEDIM; n++) {
              for (int k = klo; k <= khi; k++) {
                  sample_log << std::setw(datwidth) << std::setprecision(datprecision) << my_line_vels_arr(i,j,k,n);
              }
          }
          for (int k = klo; k <= khi; k++) {
              sample_log << std::setw(datwidth) << std::setprecision(datprecision) << my_line_tau11_arr(i,j,k);
          }
          for (int k = klo; k <= khi; k++) {
              sample_log << std::setw(datwidth) << std::setprecision(datprecision) << my_line_tau12_arr(i,j,k);
          }
          for (int k = klo; k <= khi; k++) {
              sample_log << std::setw(datwidth) << std::setprecision(datprecision) << my_line_tau13_arr(i,j,k);
          }
          for (int k = klo; k <= khi; k++) {
              sample_log << std::setw(datwidth) << std::setprecision(datprecision) << my_line_tau22_arr(i,j,k);
          }
          for (int k = klo; k <= khi; k++) {
              sample_log << std::setw(datwidth) << std::setprecision(datprecision) << my_line_tau23_arr(i,j,k);
          }
          for (int k = klo; k <= khi; k++) {
              sample_log << std::setw(datwidth) << std::setprecision(datprecision) << my_line_tau33_arr(i,j,k);
          }
          sample_log << std::endl;
        } // if good
    } // mfi
}

/**
 * Utility function for computing a volume weighted sum of MultiFab data for a single component
 *
 * @param lev Current level
 * @param mf MultiFab on which we do the volume weighted sum
 * @param comp Index of the component we want to sum
 * @param local Boolean sets whether or not to reduce the sum over the domain (false) or compute sums local to each MPI rank (true)
 * @param finemask If a finer level is available, determines whether we mask fine data
 */
Real
ERF::volWgtSumMF (int lev,
                  const MultiFab& mf, int comp,
                  const MultiFab& mapfac, bool local, bool finemask)
{
    BL_PROFILE("ERF::volWgtSumMF()");

    Real sum = 0.0;
    MultiFab tmp(grids[lev], dmap[lev], 1, 0);
    MultiFab::Copy(tmp, mf, comp, 0, 1, 0);

    // The quantity that is conserved is not (rho S), but rather (rho S / m^2) where
    // m is the map scale factor at cell centers
    for (MFIter mfi(tmp, TilingIfNotGPU()); mfi.isValid(); ++mfi) {
        const Box& bx = mfi.tilebox();
        const Array4<      Real>    tmp_arr =    tmp.array(mfi);
        const Array4<const Real> mapfac_arr = mapfac.const_array(mfi);
        ParallelFor(bx, [=] AMREX_GPU_DEVICE (int i, int j, int k) noexcept
        {
            tmp_arr(i,j,k) /= (mapfac_arr(i,j,0)*mapfac_arr(i,j,0));
        });
    } // mfi

    if (lev < finest_level && finemask) {
        const MultiFab& mask = build_fine_mask(lev+1);
        MultiFab::Multiply(tmp, mask, 0, 0, 1, 0);
    }

    MultiFab volume(grids[lev], dmap[lev], 1, 0);
    auto const& dx = geom[lev].CellSizeArray();
    Real cell_vol = dx[0]*dx[1]*dx[2];
    volume.setVal(cell_vol);
    if (solverChoice.use_terrain)
        MultiFab::Multiply(volume, *detJ_cc[lev], 0, 0, 1, 0);
    sum = MultiFab::Dot(tmp, 0, volume, 0, 1, 0, local);

    if (!local)
      ParallelDescriptor::ReduceRealSum(sum);

    return sum;
}

/**
 * Helper function for constructing a fine mask, that is, a MultiFab
 * masking coarser data at a lower level by zeroing out covered cells
 * in the fine mask MultiFab we compute.
 *
 * @param level Fine level index which masks underlying coarser data
 */
MultiFab&
ERF::build_fine_mask (int level)
{
    // Mask for zeroing covered cells
    AMREX_ASSERT(level > 0);

    const BoxArray& cba = grids[level-1];
    const DistributionMapping& cdm = dmap[level-1];

    // TODO -- we should make a vector of these a member of ERF class
    fine_mask.define(cba, cdm, 1, 0, MFInfo());
    fine_mask.setVal(1.0);

    BoxArray fba = grids[level];
    iMultiFab ifine_mask = makeFineMask(cba, cdm, fba, ref_ratio[level-1], 1, 0);

    const auto  fma =  fine_mask.arrays();
    const auto ifma = ifine_mask.arrays();
    ParallelFor(fine_mask, [=] AMREX_GPU_DEVICE(int bno, int i, int j, int k) noexcept
    {
       fma[bno](i,j,k) = ifma[bno](i,j,k);
    });

    return fine_mask;
}

/**
 * Helper function which uses the current step number, time, and timestep to
 * determine whether it is time to take an action specified at every interval
 * of timesteps.
 *
 * @param nstep Timestep number
 * @param time Current time
 * @param dtlev Timestep for the current level
 * @param action_interval Interval in number of timesteps for taking action
 * @param action_per Interval in simulation time for taking action
 */
bool
ERF::is_it_time_for_action (int nstep, Real time, Real dtlev, int action_interval, Real action_per)
{
  bool int_test = (action_interval > 0 && nstep % action_interval == 0);

  bool per_test = false;
  if (action_per > 0.0) {
    const int num_per_old = static_cast<int>(amrex::Math::floor((time - dtlev) / action_per));
    const int num_per_new = static_cast<int>(amrex::Math::floor((time) / action_per));

    if (num_per_old != num_per_new) {
      per_test = true;
    }
  }

  return int_test || per_test;
}
