#include "AMReX_PhysBCFunct.H"
#include <ERF_PhysBCFunct.H>

//
// mf is the multifab to be filled
// icomp is the index into the MultiFab -- if cell-centered this can be any value
//       from 0 to NVAR-1, if face-centered this must be 0
// ncomp is the number of components -- if cell-centered (var_idx = 0) this can be any value
//       from 1 to NVAR as long as icomp+ncomp <= NVAR-1.  If face-centered this
//       must be 1
// nghost is how many ghost cells to be filled
// time is the time at which the data should be filled
// bccomp is the index into both domain_bcs_type_bcr and bc_extdir_vals for icomp = 0  --
//     so this follows the BCVars enum
//
void ERFPhysBCFunct::operator() (MultiFab& mf, int icomp, int ncomp, IntVect const& nghost,
                                 Real time, int bccomp)
{
        if (m_geom.isAllPeriodic()) return;

        BL_PROFILE("ERFPhysBCFunct::()");

        const auto& domain = m_geom.Domain();
        const auto& dom_lo = amrex::lbound(domain);
        const auto& dom_hi = amrex::ubound(domain);

#ifdef ERF_USE_TERRAIN
        // Private data from constructor
        const GpuArray<Real, AMREX_SPACEDIM> dxInv = m_geom.InvCellSizeArray();
#endif

        // Create a grown domain box containing valid + periodic cells
        Box gdomain = amrex::convert(domain, mf.boxArray().ixType());
        for (int i = 0; i < AMREX_SPACEDIM; ++i) {
            if (m_geom.isPeriodic(i)) {
                gdomain.grow(i, nghost[i]);
            }
        }

#ifdef AMREX_USE_OMP
#pragma omp parallel if (Gpu::notInLaunchRegion())
#endif
        {
            Vector<BCRec> bcrs(ncomp);

#ifdef ERF_USE_TERRAIN
            // We only need the metric terms if we are imposing no-flow on the vertical velocity "w" at the bottom boundary
            // Make new temporary multifabs and copy to match the mf that we are trying to fill
            MultiFab temp_z_phys_nd;
            MultiFab temp_xvel;
            MultiFab temp_yvel;
            if (bccomp == BCVars::zvel_bc) {

               // Copy with all ghost cells
               BoxArray ba_nd(mf.boxArray());
               ba_nd.convert(IntVect(1,1,1));
               amrex::Print() << "ND NGROW " << m_z_phys_nd.nGrow() << std::endl;
               temp_z_phys_nd.define(ba_nd,mf.DistributionMap(),1,m_z_phys_nd.nGrow());
               temp_z_phys_nd.ParallelCopy(m_z_phys_nd,0,0,1,m_z_phys_nd.nGrow(),m_z_phys_nd.nGrow());

               // Copy with all ghost cells
               BoxArray ba_x(mf.boxArray());
               ba_x.convert(IntVect(1,0,0));
               MultiFab xvel_mf(m_data.get_var(Vars::xvel), amrex::make_alias, 0, 1);
               temp_xvel.define(ba_x, mf.DistributionMap(),1,xvel_mf.nGrow());
               amrex::Print() << "XVEL NGROW " << xvel_mf.nGrow() << std::endl;

               temp_xvel.ParallelCopy(xvel_mf,0,0,1,xvel_mf.nGrow(),xvel_mf.nGrow());

               // Copy with all ghost cells
               BoxArray ba_y(mf.boxArray());
               ba_y.convert(IntVect(0,1,0));
               MultiFab yvel_mf(m_data.get_var(Vars::yvel), amrex::make_alias, 0, 1);
               amrex::Print() << "YVEL NGROW " << yvel_mf.nGrow() << std::endl;
            amrex::Print() << " " << std::endl;

               temp_yvel.define(ba_y, mf.DistributionMap(),1,yvel_mf.nGrow());
               temp_yvel.ParallelCopy(yvel_mf,0,0,1,yvel_mf.nGrow(),yvel_mf.nGrow());
            }
#endif

            // Do all BCs except MOST
            for (MFIter mfi(mf); mfi.isValid(); ++mfi)
            {
                FArrayBox& dest = mf[mfi];
                const Array4<Real>& dest_array = mf.array(mfi);
                const Box& bx = mfi.fabbox();

#ifdef ERF_USE_TERRAIN
                auto z_nd     = (bccomp == BCVars::zvel_bc) ? temp_z_phys_nd.const_array(mfi) : Array4<const Real>{};
                auto velx_arr = (bccomp == BCVars::zvel_bc) ?      temp_xvel.const_array(mfi) : Array4<const Real>{};
                auto vely_arr = (bccomp == BCVars::zvel_bc) ?      temp_yvel.const_array(mfi) : Array4<const Real>{};
#endif
                //! if there are cells not in the valid + periodic grown box
                //! we need to fill them here
                //!
                if (!gdomain.contains(bx))
                {
                    //! Based on BCRec for the domain, we need to make BCRec for this Box
                    // bccomp is used as starting index for m_domain_bcs_type
                    //      0 is used as starting index for bcrs
                    amrex::setBC(bx, domain, bccomp, 0, ncomp, m_domain_bcs_type, bcrs);

                    // Call the default fill functions
                    //! Note that we pass 0 as starting component of bcrs.
                    GpuBndryFuncFab<NullFill> bndry_fill_cc_fc_nd(NullFill{});

                    // Calls routines to fill all the foextrap, hoextrap, etc types of bc's
                    bndry_fill_cc_fc_nd(bx, dest, icomp, ncomp, m_geom, time, bcrs, 0, bccomp);

                    // xlo: ori = 0
                    // ylo: ori = 1
                    // zlo: ori = 2
                    // xhi: ori = 3
                    // yhi: ori = 4
                    // zhi: ori = 5

                    amrex::Gpu::DeviceVector<BCRec> bcrs_d(ncomp);
#ifdef AMREX_USE_GPU
                    Gpu::htod_memcpy
                        (bcrs_d.data(), bcrs.data(), sizeof(BCRec)*ncomp);
#else
                    std::memcpy
                        (bcrs_d.data(), bcrs.data(), sizeof(BCRec)*ncomp);
#endif

                    if (m_var_idx == Vars::xvel || m_var_idx == Vars::xmom ||
                        m_var_idx == Vars::yvel || m_var_idx == Vars::ymom ||
                        m_var_idx == Vars::zvel || m_var_idx == Vars::zmom) {
                        AMREX_ALWAYS_ASSERT(ncomp == 1 && icomp == 0);
                    } else {
                        AMREX_ALWAYS_ASSERT(icomp+ncomp <= NVAR);
                    }

                    amrex::GpuArray<amrex::GpuArray<amrex::Real, AMREX_SPACEDIM*2>,
                                                                 AMREX_SPACEDIM+NVAR> l_bc_extdir_vals_d;

                    for (int i = 0; i < ncomp; i++)
                        for (int ori = 0; ori < 2*AMREX_SPACEDIM; ori++)
                            l_bc_extdir_vals_d[i][ori] = m_bc_extdir_vals[bccomp+i][ori];

                    const amrex::BCRec* bc_ptr = bcrs_d.data();

                    // Fill here all the "generic" ext_dir bc's
                    ParallelFor(bx, ncomp, [=] AMREX_GPU_DEVICE (int i, int j, int k, int n)
                    {
#ifdef ERF_USE_TERRAIN
                      // BC for W with terrain (U & V filled before W)
                      if (bccomp == BCVars::zvel_bc)
                      {
                        // Populate W ghost cells on lateral boundaries
                        if (i < dom_lo.x && bc_ptr[n].lo(0) == ERFBCType::ext_dir)
                          dest_array(i,j,k,icomp+n) = WFromOmegaBC(i,j,k,l_bc_extdir_vals_d[n][0],
                                                                   velx_arr,vely_arr,z_nd,dxInv);
                        if (j < dom_lo.y && bc_ptr[n].lo(1) == ERFBCType::ext_dir)
                          dest_array(i,j,k,icomp+n) = WFromOmegaBC(i,j,k,l_bc_extdir_vals_d[n][1],
                                                                   velx_arr,vely_arr,z_nd,dxInv);
                        if (i > dom_hi.x && bc_ptr[n].hi(0) == ERFBCType::ext_dir)
                          dest_array(i,j,k,icomp+n) = WFromOmegaBC(i,j,k,l_bc_extdir_vals_d[n][3],
                                                                   velx_arr,vely_arr,z_nd,dxInv);
                        if (j > dom_hi.y && bc_ptr[n].hi(1) == ERFBCType::ext_dir)
                          dest_array(i,j,k,icomp+n) = WFromOmegaBC(i,j,k,l_bc_extdir_vals_d[n][4],
                                                                   velx_arr,vely_arr,z_nd,dxInv);

                        // Populate W ghost values on bottom boundary
                        // Note we use the same values for all z <= domlo.z
                        if (k <= dom_lo.z && bc_ptr[n].lo(2) == ERFBCType::ext_dir)
                        {
                          dest_array(i,j,k,icomp+n) = WFromOmegaBC(i,j,dom_lo.z,l_bc_extdir_vals_d[n][2],
                                                                   velx_arr,vely_arr,z_nd,dxInv);
                        }

                        // Populate W ghost values on top boundary
                        // Note we use the same values for all z >= domhi.z+1
                        if (k > dom_hi.z && bc_ptr[n].hi(2) == ERFBCType::ext_dir)
                        {
                          dest_array(i,j,k,icomp+n) = WFromOmegaBC(i,j,dom_hi.z+1,l_bc_extdir_vals_d[n][5],
                                                                   velx_arr,vely_arr, z_nd,dxInv);
                        }

                        } else {
                          // Populate ghost cells & upper face on top boundary
                          if (i < dom_lo.x && bc_ptr[n].lo(0) == ERFBCType::ext_dir)
                            dest_array(i,j,k,icomp+n) = l_bc_extdir_vals_d[n][0];
                          if (j < dom_lo.y && bc_ptr[n].lo(1) == ERFBCType::ext_dir)
                            dest_array(i,j,k,icomp+n) = l_bc_extdir_vals_d[n][1];
                          if (k < dom_lo.z && bc_ptr[n].lo(2) == ERFBCType::ext_dir)
                            dest_array(i,j,k,icomp+n) = l_bc_extdir_vals_d[n][2];
                          if (i > dom_hi.x && bc_ptr[n].hi(0) == ERFBCType::ext_dir)
                            dest_array(i,j,k,icomp+n) = l_bc_extdir_vals_d[n][3];
                          if (j > dom_hi.y && bc_ptr[n].hi(1) == ERFBCType::ext_dir)
                            dest_array(i,j,k,icomp+n) = l_bc_extdir_vals_d[n][4];
                          if (k > dom_hi.z && bc_ptr[n].hi(2) == ERFBCType::ext_dir)
                            dest_array(i,j,k,icomp+n) = l_bc_extdir_vals_d[n][5];

                          // Populate face values on bottom boundary (except for W)
                          if (bccomp == BCVars::xvel_bc)
                          {
                            if (i == dom_lo.x && bc_ptr[n].lo(0) == ERFBCType::ext_dir)
                              dest_array(i,j,k,icomp+n) = l_bc_extdir_vals_d[n][0];
                          }
                          if (bccomp == BCVars::yvel_bc)
                          {
                            if (j == dom_lo.y && bc_ptr[n].lo(1) == ERFBCType::ext_dir)
                              dest_array(i,j,k,icomp+n) = l_bc_extdir_vals_d[n][1];
                          }
                      }
#else
                        // Populate ghost cells & upper face on top boundary
                        if (i < dom_lo.x && bc_ptr[n].lo(0) == ERFBCType::ext_dir)
                            dest_array(i,j,k,icomp+n) = l_bc_extdir_vals_d[n][0];
                        if (j < dom_lo.y && bc_ptr[n].lo(1) == ERFBCType::ext_dir)
                            dest_array(i,j,k,icomp+n) = l_bc_extdir_vals_d[n][1];
                        if (k < dom_lo.z && bc_ptr[n].lo(2) == ERFBCType::ext_dir)
                            dest_array(i,j,k,icomp+n) = l_bc_extdir_vals_d[n][2];
                        if (i > dom_hi.x && bc_ptr[n].hi(0) == ERFBCType::ext_dir)
                            dest_array(i,j,k,icomp+n) = l_bc_extdir_vals_d[n][3];
                        if (j > dom_hi.y && bc_ptr[n].hi(1) == ERFBCType::ext_dir)
                            dest_array(i,j,k,icomp+n) = l_bc_extdir_vals_d[n][4];
                        if (k > dom_hi.z && bc_ptr[n].hi(2) == ERFBCType::ext_dir)
                            dest_array(i,j,k,icomp+n) = l_bc_extdir_vals_d[n][5];

                        // Populate face values on bottom boundary
                        if (bccomp == BCVars::xvel_bc)
                        {
                          if (i == dom_lo.x && bc_ptr[n].lo(0) == ERFBCType::ext_dir)
                            dest_array(i,j,k,icomp+n) = l_bc_extdir_vals_d[n][0];
                        }
                        if (bccomp == BCVars::yvel_bc)
                        {
                          if (j == dom_lo.y && bc_ptr[n].lo(1) == ERFBCType::ext_dir)
                            dest_array(i,j,k,icomp+n) = l_bc_extdir_vals_d[n][1];
                        }
                        if (bccomp == BCVars::zvel_bc)
                        {
                          if (k == dom_lo.z && bc_ptr[n].lo(2) == ERFBCType::ext_dir)
                            dest_array(i,j,k,icomp+n) = l_bc_extdir_vals_d[n][2];
                        }
#endif
                    });

                    // This uses data read in as BndryRegisters from a previous ERF run
                    if (m_r2d) {
                        amrex::Vector<std::unique_ptr<PlaneVector>>& bndry_data = m_r2d->interp_in_time(time);
                        const auto& bdatxlo = (*bndry_data[0])[m_lev].const_array();
                        const auto& bdatylo = (*bndry_data[1])[m_lev].const_array();
                        // const auto& bdatzlo = (*bndry_data[2])[m_lev].const_array();
                        const auto& bdatxhi = (*bndry_data[3])[m_lev].const_array();
                        const auto& bdatyhi = (*bndry_data[4])[m_lev].const_array();
                        // const auto& bdatzhi = (*bndry_data[5])[m_lev].const_array();

                        // Fill here all the boundary conditions which are supplied by
                        // planes we have read in and are interpolating in time
                        ParallelFor(bx, ncomp, [=] AMREX_GPU_DEVICE (int i, int j, int k, int n)
                        {
                            if (i < dom_lo.x && bc_ptr[n].lo(0) == ERFBCType::ext_dir_ingested) {
                                int jb = std::min(std::max(j,dom_lo.y),dom_hi.y);
                                int kb = std::min(std::max(k,dom_lo.z),dom_hi.z);
                                dest_array(i,j,k,icomp+n) = bdatxlo(dom_lo.x-1,jb,kb,bccomp+n);
                            }
                            if (j < dom_lo.y && bc_ptr[n].lo(1) == ERFBCType::ext_dir_ingested) {
                                int ib = std::min(std::max(i,dom_lo.x),dom_hi.x);
                                int kb = std::min(std::max(k,dom_lo.z),dom_hi.z);
                                dest_array(i,j,k,icomp+n) = bdatylo(ib,dom_lo.y-1,kb,bccomp+n);
                            }
                            if (k < dom_lo.z && bc_ptr[n].lo(2) == ERFBCType::ext_dir_ingested) {
                                // int ib = std::min(std::max(i,dom_lo.x),dom_hi.x);
                                // int jb = std::min(std::max(j,dom_lo.y),dom_hi.y);
                                // dest_array(i,j,k,icomp+n) = bdatzlo(ib,jb,dom_lo.z-1,bccomp+n);
                            }
                            if (i > dom_hi.x && bc_ptr[n].hi(0) == ERFBCType::ext_dir_ingested) {
                                int jb = std::min(std::max(j,dom_lo.y),dom_hi.y);
                                int kb = std::min(std::max(k,dom_lo.z),dom_hi.z);
                                dest_array(i,j,k,icomp+n) = bdatxhi(dom_hi.x+1,jb,kb,bccomp+n);
                            }
                            if (j > dom_hi.y && bc_ptr[n].hi(1) == ERFBCType::ext_dir_ingested) {
                                int ib = std::min(std::max(i,dom_lo.x),dom_hi.x);
                                int kb = std::min(std::max(k,dom_lo.z),dom_hi.z);
                                dest_array(i,j,k,icomp+n) = bdatyhi(ib,dom_hi.y+1,kb,bccomp+n);
                            }
                            if (k > dom_hi.z && bc_ptr[n].hi(2) == ERFBCType::ext_dir_ingested) {
                                // int ib = std::min(std::max(i,dom_lo.x),dom_hi.x);
                                // int jb = std::min(std::max(j,dom_lo.y),dom_hi.y);
                                // dest_array(i,j,k,icomp+n) = bdatzhi(ib,jb,dom_hi.z+1,bccomp+n);
                            }

                            if (bccomp == BCVars::xvel_bc)
                            {
                                if (i == dom_lo.x && bc_ptr[n].lo(0) == ERFBCType::ext_dir_ingested) {
                                    int jb = std::min(std::max(j,dom_lo.y),dom_hi.y);
                                    int kb = std::min(std::max(k,dom_lo.z),dom_hi.z);
                                    dest_array(i,j,k,icomp+n) = bdatxlo(dom_lo.x-1,jb,kb,bccomp+n);
                                }
                            }
                            if (bccomp == BCVars::yvel_bc)
                            {
                                if (j == dom_lo.y && bc_ptr[n].lo(1) == ERFBCType::ext_dir_ingested) {
                                    int ib = std::min(std::max(i,dom_lo.x),dom_hi.x);
                                    int kb = std::min(std::max(k,dom_lo.z),dom_hi.z);
                                    dest_array(i,j,k,icomp+n) = bdatylo(ib,dom_lo.y-1,kb,bccomp+n);
                                }
                            }
                            if (bccomp == BCVars::zvel_bc)
                            {
                                if (k == dom_lo.z && bc_ptr[n].lo(2) == ERFBCType::ext_dir_ingested) {
                                    // int ib = std::min(std::max(i,dom_lo.x),dom_hi.x);
                                    // int jb = std::min(std::max(j,dom_lo.y),dom_hi.y);
                                    // dest_array(i,j,k,icomp+n) = bdatzlo(ib,jb,dom_lo.z-1,bccomp+n);
                                }
                            }
                        });
                    }

#ifdef ERF_USE_NETCDF
                    //
                    // This uses data read in from WRF real.exe output
                    //
                    if (0) {
                        //
                        // We assume the data has been read in at the level 0 resolution
                        //

                        // This is a total HACK
                        amrex::Real dT = 1.0;

                        int n = time / dT;
                        amrex::Real alpha = (time - n * dT) / dT;
                        amrex::Real oma   = 1.0 - alpha;

                        //
                        // We have data at fixed time intervals we will call dT
                        // Then to interpolate, given time, we can define n = (time/dT)
                        // and alpha = (time - n*dT) / dT, then we define the data at time
                        // as  alpha * (data at time n+1) + (1 - alpha) * (data at time n)
                        //
                        const auto& bdatxlo_n   = m_bdy_data_xlo[n  ].const_array();
                        const auto& bdatxlo_np1 = m_bdy_data_xlo[n+1].const_array();
                        const auto& bdatxhi_n   = m_bdy_data_xhi[n  ].const_array();
                        const auto& bdatxhi_np1 = m_bdy_data_xhi[n+1].const_array();
                        const auto& bdatylo_n   = m_bdy_data_ylo[n  ].const_array();
                        const auto& bdatylo_np1 = m_bdy_data_ylo[n+1].const_array();
                        const auto& bdatyhi_n   = m_bdy_data_yhi[n  ].const_array();
                        const auto& bdatyhi_np1 = m_bdy_data_yhi[n+1].const_array();

                        if (bccomp == BCVars::xvel_bc)
                           amrex::Print() << "XVEL BC " << icomp << " " << ncomp << std::endl;
                        else if (bccomp == BCVars::yvel_bc)
                           amrex::Print() << "YVEL BC " << icomp << " " << ncomp << std::endl;
                        else if (bccomp == BCVars::zvel_bc)
                           amrex::Print() << "ZVEL BC " << icomp << " " << ncomp << std::endl;
                        else
                           amrex::Print() << "BC IC " << bccomp << " " << icomp << " " << ncomp << std::endl;

                        // Fill here all the boundary conditions which are supplied by
                        // planes we have read in and are interpolating in time
                        ParallelFor(bx, ncomp, [=] AMREX_GPU_DEVICE (int i, int j, int k, int n)
                        {
                            if (i < dom_lo.x && bc_ptr[n].lo(0) == ERFBCType::ext_dir_ingested) {
                                int jb = std::min(std::max(j,dom_lo.y),dom_hi.y);
                                int kb = std::min(std::max(k,dom_lo.z),dom_hi.z);
                                dest_array(i,j,k,icomp+n) = oma * bdatxlo_n  (dom_lo.x-1,jb,kb,bccomp+n)
                                                        + alpha * bdatxlo_np1(dom_lo.x-1,jb,kb,bccomp+n);
                            }
                            if (j < dom_lo.y && bc_ptr[n].lo(1) == ERFBCType::ext_dir_ingested) {
                                int ib = std::min(std::max(i,dom_lo.x),dom_hi.x);
                                int kb = std::min(std::max(k,dom_lo.z),dom_hi.z);
                                dest_array(i,j,k,icomp+n) = oma * bdatylo_n  (ib,dom_lo.y-1,kb,bccomp+n)
                                                        + alpha * bdatylo_np1(ib,dom_lo.y-1,kb,bccomp+n);
                            }
                            if (k < dom_lo.z && bc_ptr[n].lo(2) == ERFBCType::ext_dir_ingested) {
                                // int ib = std::min(std::max(i,dom_lo.x),dom_hi.x);
                                // int jb = std::min(std::max(j,dom_lo.y),dom_hi.y);
                                // dest_array(i,j,k,icomp+n) = oma * bdatzlo_n  (ib,jb,dom_lo.z-1,bccomp+n)
                                //                         + alpha * bdatzlo_np1(ib,jb,dom_lo.z-1,bccomp+n);
                            }
                            if (i > dom_hi.x && bc_ptr[n].hi(0) == ERFBCType::ext_dir_ingested) {
                                int jb = std::min(std::max(j,dom_lo.y),dom_hi.y);
                                int kb = std::min(std::max(k,dom_lo.z),dom_hi.z);
                                dest_array(i,j,k,icomp+n) = oma * bdatxhi_n  (dom_hi.x+1,jb,kb,bccomp+n)
                                                        + alpha * bdatxhi_np1(dom_hi.x+1,jb,kb,bccomp+n);
                            }
                            if (j > dom_hi.y && bc_ptr[n].hi(1) == ERFBCType::ext_dir_ingested) {
                                int ib = std::min(std::max(i,dom_lo.x),dom_hi.x);
                                int kb = std::min(std::max(k,dom_lo.z),dom_hi.z);
                                dest_array(i,j,k,icomp+n) = oma * bdatyhi_n  (ib,dom_hi.y+1,kb,bccomp+n)
                                                        + alpha * bdatyhi_np1(ib,dom_hi.y+1,kb,bccomp+n);
                            }
                            if (k > dom_hi.z && bc_ptr[n].hi(2) == ERFBCType::ext_dir_ingested) {
                                // int ib = std::min(std::max(i,dom_lo.x),dom_hi.x);
                                // int jb = std::min(std::max(j,dom_lo.y),dom_hi.y);
                                // dest_array(i,j,k,icomp+n) = oma * bdatzhi_n  (ib,jb,dom_hi.z+1,bccomp+n)
                                //                         + alpha * bdatzhi_np1(ib,jb,dom_hi.z+1,bccomp+n);
                            }

                            if (bccomp == BCVars::xvel_bc)
                            {
                                if (i == dom_lo.x && bc_ptr[n].lo(0) == ERFBCType::ext_dir_ingested) {
                                    int jb = std::min(std::max(j,dom_lo.y),dom_hi.y);
                                    int kb = std::min(std::max(k,dom_lo.z),dom_hi.z);
                                    dest_array(i,j,k,icomp+n) = oma * bdatxlo_n  (dom_lo.x-1,jb,kb,bccomp+n)
                                                            + alpha * bdatxlo_np1(dom_lo.x-1,jb,kb,bccomp+n);
                                }
                            }
                            if (bccomp == BCVars::yvel_bc)
                            {
                                if (j == dom_lo.y && bc_ptr[n].lo(1) == ERFBCType::ext_dir_ingested) {
                                    int ib = std::min(std::max(i,dom_lo.x),dom_hi.x);
                                    int kb = std::min(std::max(k,dom_lo.z),dom_hi.z);
                                    dest_array(i,j,k,icomp+n) = oma * bdatylo_n  (ib,dom_lo.y-1,kb,bccomp+n)
                                                            + alpha * bdatylo_np1(ib,dom_lo.y-1,kb,bccomp+n);
                                }
                            }
                            if (bccomp == BCVars::zvel_bc)
                            {
                                if (k == dom_lo.z && bc_ptr[n].lo(2) == ERFBCType::ext_dir_ingested) {
                                    // int ib = std::min(std::max(i,dom_lo.x),dom_hi.x);
                                    // int jb = std::min(std::max(j,dom_lo.y),dom_hi.y);
                                    // dest_array(i,j,k,icomp+n) = oma * bdatzlo_n  (ib,jb,dom_lo.z-1,bccomp+n)
                                    //                         + alpha * bdatzlo_np1(ib,jb,dom_lo.z-1,bccomp+n);
                                }
                            }
                        });
                    }
#endif
                } // !gdomain.contains(bx)
            } // MFIter
        } // OpenMP
    } // operator()
