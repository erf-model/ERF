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

            for (MFIter mfi(mf); mfi.isValid(); ++mfi)
            {
                FArrayBox& dest = mf[mfi];
                const Array4<Real>& dest_array = mf.array(mfi);
                const Box& bx = mfi.fabbox();

#ifdef ERF_USE_TERRAIN
                // Private data from constructor
                const Array4<const Real>& z_nd = m_z_phys_nd.const_array(mfi);
                const Array4<const Real>& detJ = m_detJ_cc.const_array(mfi);
                const auto velx_arr = m_data.get_var(Vars::xvel)[mfi].array();
                const auto vely_arr = m_data.get_var(Vars::yvel)[mfi].array();
                const auto cons_arr = m_has_most_bcs ? m_data.get_var(Vars::cons)[mfi].array() : Array4<Real>();
                const auto eta_arr  = m_has_most_bcs ? m_viscosity[mfi].array() : Array4<Real>();
#else
                // make Array4's for our data
                const auto cons_arr = m_has_most_bcs ? m_data.get_var(Vars::cons)[mfi].array() : Array4<Real>();
                const auto velx_arr = m_has_most_bcs ? m_data.get_var(Vars::xvel)[mfi].array() : Array4<Real>();
                const auto vely_arr = m_has_most_bcs ? m_data.get_var(Vars::yvel)[mfi].array() : Array4<Real>();
                const auto eta_arr  = m_has_most_bcs ? m_viscosity[mfi].array() : Array4<Real>();
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

                    // Set bc-type at low z to reflect_odd for zvel/zmom when using MOST BCs
                    for (int i = 0; i < bcrs.size(); i++) {
                        // ori = 2 is the zlo side
                        if (bcrs[i].data()[2] == static_cast<int>(ERFBCType::MOST)) {
                            if (m_var_idx == Vars::zvel || m_var_idx == Vars::zmom) {
                                bcrs[i].setLo(2, ERFBCType::reflect_odd);
                            }
                        }
                    }

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
                      if (bccomp == BCVars::zvel_bc) {
                        // Populate W ghost cells & upper face on top boundary
                        if (i < dom_lo.x && bc_ptr[n].lo(0) == ERFBCType::ext_dir)
                          dest_array(i,j,k,icomp+n) = WFromOmegaBC(i,j,k,l_bc_extdir_vals_d[n][0],
                                                                   velx_arr,vely_arr,
                                                                   z_nd,dxInv);
                        if (j < dom_lo.y && bc_ptr[n].lo(1) == ERFBCType::ext_dir)
                          dest_array(i,j,k,icomp+n) = WFromOmegaBC(i,j,k,l_bc_extdir_vals_d[n][1],
                                                                   velx_arr,vely_arr,
                                                                   z_nd,dxInv);
                        if (k < dom_lo.z && bc_ptr[n].lo(2) == ERFBCType::ext_dir)
                          dest_array(i,j,k,icomp+n) = WFromOmegaBC(i,j,k,l_bc_extdir_vals_d[n][2],
                                                                   velx_arr,vely_arr,
                                                                   z_nd,dxInv);
                        if (i > dom_hi.x && bc_ptr[n].hi(0) == ERFBCType::ext_dir)
                          dest_array(i,j,k,icomp+n) = WFromOmegaBC(i,j,k,l_bc_extdir_vals_d[n][3],
                                                                   velx_arr,vely_arr,
                                                                   z_nd,dxInv);
                        if (j > dom_hi.y && bc_ptr[n].hi(1) == ERFBCType::ext_dir)
                          dest_array(i,j,k,icomp+n) = WFromOmegaBC(i,j,k,l_bc_extdir_vals_d[n][4],
                                                                   velx_arr,vely_arr,
                                                                   z_nd,dxInv);
                        if (k > dom_hi.z && bc_ptr[n].hi(2) == ERFBCType::ext_dir)
                          dest_array(i,j,k,icomp+n) = WFromOmegaBC(i,j,k,l_bc_extdir_vals_d[n][5],
                                                                   velx_arr,vely_arr,
                                                                   z_nd,dxInv);
                        // Populate W face value on bottom boundary
                        if (k == dom_lo.z && bc_ptr[n].lo(2) == ERFBCType::ext_dir)
                          dest_array(i,j,k,icomp+n) = WFromOmegaBC(i,j,k,l_bc_extdir_vals_d[n][2],
                                                                   velx_arr,vely_arr,
                                                                   z_nd,dxInv);
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
                        }

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
                        //
                        // TODO: We need to know which level we're actually at -- this is a hack
                        //
                        int lev = 0;
                        amrex::Vector<std::unique_ptr<PlaneVector>>& bndry_data = m_r2d->interp_in_time(time);
                        const auto& bdatxlo = (*bndry_data[0])[lev].const_array();
                        const auto& bdatylo = (*bndry_data[1])[lev].const_array();
                        // const auto& bdatzlo = (*bndry_data[2])[lev].const_array();
                        const auto& bdatxhi = (*bndry_data[3])[lev].const_array();
                        const auto& bdatyhi = (*bndry_data[4])[lev].const_array();
                        // const auto& bdatzhi = (*bndry_data[5])[lev].const_array();

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
                        int lev = 0;

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

                    int zlo = m_geom.Domain().smallEnd(2);

                    // Now handle the MOST bc if we have any
                    // Note: if we use MOST for one variable we use it for all, so we can just test on one component

                    if ( bx.smallEnd()[2] < zlo &&
                         ( (m_var_idx == Vars::cons && bcrs[Cons::RhoTheta].lo(2) == ERFBCType::MOST) ||
                           (m_var_idx == Vars::xvel && bcrs[0].lo(2) == ERFBCType::MOST) ||
                           (m_var_idx == Vars::yvel && bcrs[0].lo(2) == ERFBCType::MOST) ) )
                    {
                        // check m_var_idx to distinguish between Vars::Dvel and Vars::Dmom
                        // (in Legacy this was controlled by the is_derived flag)
                        bool var_is_derived = false;
                        if (m_var_idx == Vars::xvel || m_var_idx == Vars::yvel || m_var_idx == Vars::zvel) {
                            var_is_derived = true;
                        }

                        amrex::Box b2d = bx; // Copy constructor
                        b2d.setBig(2,zlo-1);

                        /**
                        * NOTE: the number of ghost zone for state variables are different from face centered
                        *       variables in the new version.
                        */

                        if (m_var_idx == Vars::cons) {
                            int n = Cons::RhoTheta;
                            ParallelFor(b2d, [=,m_most=m_most] AMREX_GPU_DEVICE (int i, int j, int k)
                            {
                                Real velx, vely, rho, theta, eta;
                                int ix, jx, iy, jy, ie, je;

                                ix = i < lbound(velx_arr).x ? lbound(velx_arr).x : i;
                                jx = j < lbound(velx_arr).y ? lbound(velx_arr).y : j;
                                ix = ix > ubound(velx_arr).x-1 ? ubound(velx_arr).x-1 : ix;
                                jx = jx > ubound(velx_arr).y ? ubound(velx_arr).y : jx;

                                iy = i < lbound(vely_arr).x ? lbound(vely_arr).x : i;
                                jy = j < lbound(vely_arr).y ? lbound(vely_arr).y : j;
                                iy = iy > ubound(vely_arr).x ? ubound(vely_arr).x : iy;
                                jy = jy > ubound(vely_arr).y-1 ? ubound(vely_arr).y-1 : jy;

                                ie = i < lbound(eta_arr).x ? lbound(eta_arr).x : i;
                                je = j < lbound(eta_arr).y ? lbound(eta_arr).y : j;
                                ie = ie > ubound(eta_arr).x ? ubound(eta_arr).x : ie;
                                je = je > ubound(eta_arr).y ? ubound(eta_arr).y : je;

                                velx  = 0.5*(velx_arr(ix,jx,zlo)+velx_arr(ix+1,jx,zlo));
                                vely  = 0.5*(vely_arr(iy,jy,zlo)+vely_arr(iy,jy+1,zlo));
                                rho   = cons_arr(ie,je,zlo,Rho_comp);
                                theta = cons_arr(ie,je,zlo,RhoTheta_comp)/rho;

                                // TODO: Verify this is the correct Diff component
                                eta   = eta_arr(ie,je,zlo,EddyDiff::Mom_h);

                                Real vmag    = sqrt(velx*velx+vely*vely);
                                Real num1    = (theta-m_most.theta_mean)*m_most.vmag_mean;
                                Real num2    = (m_most.theta_mean-m_most.surf_temp)*vmag;
                                Real motheta = (num1+num2)*m_most.utau*m_most.kappa/m_most.phi_h();

                                if (!var_is_derived) {
                                    dest_array(i,j,k,icomp+n) = rho*(m_most.surf_temp + motheta*rho/eta);
                                } else {
                                    dest_array(i,j,k,icomp+n) = m_most.surf_temp + motheta/eta;
                                }
                            });

                        } else if (m_var_idx == Vars::xvel || m_var_idx == Vars::xmom) { //for velx

                            ParallelFor(b2d, [=,m_most=m_most] AMREX_GPU_DEVICE (int i, int j, int k)
                            {
                                Real velx, vely, rho, eta;
                                int jy, ie, je;

                                int iylo = i <= lbound(vely_arr).x ? lbound(vely_arr).x : i-1;
                                int iyhi = i >  ubound(vely_arr).x ? ubound(vely_arr).x : i;

                                jy = j < lbound(vely_arr).y ? lbound(vely_arr).y : j;
                                jy = jy > ubound(vely_arr).y-1 ? ubound(vely_arr).y-1 : jy;

                                ie = i < lbound(eta_arr).x ? lbound(eta_arr).x : i;
                                je = j < lbound(eta_arr).y ? lbound(eta_arr).y : j;
                                ie = ie > ubound(eta_arr).x-1 ? ubound(eta_arr).x-1 : ie;
                                je = je > ubound(eta_arr).y ? ubound(eta_arr).y : je;

                                velx  = velx_arr(i,j,zlo);
                                vely  = 0.25*( vely_arr(iyhi,jy,zlo)+vely_arr(iyhi,jy+1,zlo)
                                              +vely_arr(iylo,jy,zlo)+vely_arr(iylo,jy+1,zlo));
                                rho   = 0.5*(cons_arr(ie-1,je,zlo,Rho_comp)+
                                             cons_arr(ie  ,je,zlo,Rho_comp));
                                eta   = 0.5*( eta_arr(ie-1,je,zlo,EddyDiff::Mom_h)+
                                              eta_arr(ie  ,je,zlo,EddyDiff::Mom_h));

                                Real vmag  = sqrt(velx*velx+vely*vely);
                                Real vgx   = ((velx-m_most.vel_mean[0])*m_most.vmag_mean + vmag*m_most.vel_mean[0])/
                                              (m_most.vmag_mean*m_most.vmag_mean) * m_most.utau*m_most.utau;

                                if (!var_is_derived) {
                                    dest_array(i,j,k,icomp) = dest_array(i,j,zlo,icomp) - vgx*rho/eta;
                                } else {
                                    dest_array(i,j,k,icomp) = dest_array(i,j,zlo,icomp) - vgx/eta;
                                }
                            });

                        } else if (m_var_idx == Vars::yvel || m_var_idx == Vars::ymom) { //for vely

                            ParallelFor(b2d, [=,m_most=m_most] AMREX_GPU_DEVICE (int i, int j, int k)
                            {
                                Real velx, vely, rho, eta;
                                int ix, ie, je;

                                ix = i < lbound(velx_arr).x ? lbound(velx_arr).x : i;
                                ix = ix > ubound(velx_arr).x ? ubound(velx_arr).x : ix;

                                int jxlo = j <= lbound(velx_arr).y ? lbound(velx_arr).y : j-1;
                                int jxhi = j >  ubound(velx_arr).y ? ubound(velx_arr).y : j;

                                ie = i < lbound(eta_arr).x ? lbound(eta_arr).x : i;
                                je = j < lbound(eta_arr).y ? lbound(eta_arr).y : j;
                                ie = ie > ubound(eta_arr).x ? ubound(eta_arr).x : ie;
                                je = je > ubound(eta_arr).y-1 ? ubound(eta_arr).y-1 : je;

                                velx  = 0.25*( velx_arr(ix,jxhi,zlo)+velx_arr(ix+1,jxhi,zlo)
                                              +velx_arr(ix,jxlo,zlo)+velx_arr(ix+1,jxlo,zlo));
                                vely  = vely_arr(i,j,zlo);
                                rho   = 0.5*(cons_arr(ie,je-1,zlo,Rho_comp)+
                                             cons_arr(ie,je  ,zlo,Rho_comp));
                                eta   = 0.5*(eta_arr(ie,je-1,zlo,EddyDiff::Mom_h)+
                                             eta_arr(ie,je  ,zlo,EddyDiff::Mom_h));
                                Real vmag  = sqrt(velx*velx+vely*vely);
                                Real vgy   = ((vely-m_most.vel_mean[1])*m_most.vmag_mean + vmag*m_most.vel_mean[1]) /
                                             (m_most.vmag_mean*m_most.vmag_mean)*m_most.utau*m_most.utau;

                                if (!var_is_derived) {
                                    dest_array(i,j,k,icomp) = dest_array(i,j,zlo,icomp) - vgy*rho/eta;
                                } else {
                                    dest_array(i,j,k,icomp) = dest_array(i,j,zlo,icomp) - vgy/eta;
                                }
                            });

                        } else if (m_var_idx == Vars::zvel || m_var_idx == Vars::zmom) { //for velz
                            ParallelFor(b2d, [=,m_most=m_most] AMREX_GPU_DEVICE (int i, int j, int k) {
                                dest_array(i,j,k,icomp) = (zlo-k+1)*dest_array(i,j,zlo  ,icomp) -
                                                          (zlo-k  )*dest_array(i,j,zlo+1,icomp);
                            });
                        }
                    }

/* EXAMPLE CUSTOM PHYSICAL BC FUNCTION USING ALL VARIABLES
                    // Call our custom BC fill functions here using ccbx and vars_arrays to pass
                    // the cell-centered box and all the state variables. Note that var_idx tells us which
                    // of Vars::cons, Vars::xvel, Vars::yvel, or Vars::zvel we are filling in case we need logic.
                    ParallelFor(bx, [=] AMREX_GPU_DEVICE (int i, int j, int k) {
                        // we probably want some logic to compare (i,j,k) to the domain bounds here ...
                        if (var_idx == Vars::cons) { // depending on which variable type we are filling
                            dest_array(i,j,k,0) = vars_arrays_p[Vars::xvel](i,j,k); // (use other variable types, for example)
                        }
                    });
*/
                }
            }
        }
    }
