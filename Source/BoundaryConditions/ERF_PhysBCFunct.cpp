#include "AMReX_PhysBCFunct.H"
#include "IndexDefines.H"
#include <ERF_PhysBCFunct.H>

using namespace amrex;

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
    const auto dx      = m_geom.CellSizeArray();
    const auto dxInv   = m_geom.InvCellSizeArray();

    // Create a grown domain box containing valid + periodic cells
    Box gdomain = amrex::convert(domain, mf.boxArray().ixType());
    for (int i = 0; i < AMREX_SPACEDIM; ++i) {
        if (m_geom.isPeriodic(i)) {
            gdomain.grow(i, nghost[i]);
        }
    }

    MultiFab* z_phys_ptr = nullptr;
    MultiFab* xvel_ptr   = nullptr;
    MultiFab* yvel_ptr   = nullptr;

    if (m_z_phys_nd) {
            // We must make copies of these MultiFabs onto mf's boxArray for when this operator is
            // called for a MultiFab mf that doesn't have the same boxArray
            BoxArray mf_nodal_grids = mf.boxArray();
            mf_nodal_grids.convert(IntVect(1,1,1));
            bool OnSameGrids = ( (mf_nodal_grids       == m_z_phys_nd->boxArray()        ) &&
                                 (mf.DistributionMap() == m_z_phys_nd->DistributionMap() ) );
            if (!OnSameGrids) {
                IntVect ng_z = m_z_phys_nd->nGrowVect();
                z_phys_ptr = new MultiFab(mf_nodal_grids,mf.DistributionMap(),1,ng_z);
                z_phys_ptr->ParallelCopy(*m_z_phys_nd,0,0,1,ng_z,ng_z);

                IntVect ng_u = m_data.get_var(Vars::xvel).nGrowVect();
                BoxArray ba_u(mf.boxArray());
                ba_u.convert(IntVect(1,0,0));
                xvel_ptr = new MultiFab(ba_u,mf.DistributionMap(),1,ng_u);
                xvel_ptr->ParallelCopy(m_data.get_var(Vars::xvel),0,0,1,ng_u,ng_u);

                IntVect ng_v = m_data.get_var(Vars::yvel).nGrowVect();
                BoxArray ba_v(mf.boxArray());
                ba_v.convert(IntVect(0,1,0));
                yvel_ptr = new MultiFab(ba_v,mf.DistributionMap(),1,ng_v);
                yvel_ptr->ParallelCopy(m_data.get_var(Vars::yvel),0,0,1,ng_v,ng_v);
            }
    }

#ifdef AMREX_USE_OMP
#pragma omp parallel if (Gpu::notInLaunchRegion())
#endif
        {
            Vector<BCRec> bcrs(ncomp);

            // Do all BCs except MOST
            for (MFIter mfi(mf); mfi.isValid(); ++mfi)
            {
                const Array4<Real>& dest_arr = mf.array(mfi);
                Box bx = mfi.validbox(); bx.grow(nghost);

                Array4<const Real> z_nd_arr;
                Array4<const Real> velx_arr;
                Array4<const Real> vely_arr;

                if (m_z_phys_nd)
                {
                    BoxArray mf_nodal_grids = mf.boxArray();
                    mf_nodal_grids.convert(IntVect(1,1,1));
                    bool OnSameGrids = ( (mf_nodal_grids       == m_z_phys_nd->boxArray()        ) &&
                                         (mf.DistributionMap() == m_z_phys_nd->DistributionMap() ) );
                    if (OnSameGrids) {
                        z_nd_arr = m_z_phys_nd->const_array(mfi);
                        velx_arr = m_data.get_var(Vars::xvel).const_array(mfi);
                        vely_arr = m_data.get_var(Vars::yvel).const_array(mfi);
                    } else {
                        z_nd_arr = z_phys_ptr->const_array(mfi);
                        velx_arr = xvel_ptr->const_array(mfi);
                        vely_arr = yvel_ptr->const_array(mfi);
                    }
                }

                //! if there are cells not in the valid + periodic grown box
                //! we need to fill them here
                //!
                if (!gdomain.contains(bx) || (m_var_idx == Vars::zvel))
                {
                    if (m_var_idx == Vars::xvel) {
                        AMREX_ALWAYS_ASSERT(ncomp == 1 && icomp == 0);
                        impose_xvel_bcs(dest_arr,bx,domain,
                                        z_nd_arr,dxInv,
                                        time,bccomp);

                    } else if (m_var_idx == Vars::yvel) {
                        AMREX_ALWAYS_ASSERT(ncomp == 1 && icomp == 0);
                        impose_yvel_bcs(dest_arr,bx,domain,
                                        z_nd_arr,dxInv,
                                        time,bccomp);

                    } else if (m_var_idx == Vars::zvel) {
                        AMREX_ALWAYS_ASSERT(ncomp == 1 && icomp == 0);
                        impose_zvel_bcs(dest_arr,bx,domain,
                                        velx_arr,vely_arr,z_nd_arr,dx,dxInv,
                                        time,m_time_mt,m_delta_t,
                                        bccomp,m_terrain_type);

                    } else if (m_var_idx == Vars::cons) {
                        AMREX_ALWAYS_ASSERT(icomp == 0 && icomp+ncomp <= NVAR);
                        impose_cons_bcs(dest_arr,bx,domain,
                                        z_nd_arr,dxInv,
                                        icomp,ncomp,time,bccomp);
                    } else {
                        amrex::Abort("Dont know this var_idx in ERF_PhysBC");
                    }

                    // ****************************************************************************
                    // Based on BCRec for the domain, we need to make BCRec for this Box
                    // bccomp is used as starting index for m_domain_bcs_type
                    //      0 is used as starting index for bcrs
                    // ****************************************************************************
                    amrex::setBC(bx, domain, bccomp, 0, ncomp, m_domain_bcs_type, bcrs);

                    // xlo: ori = 0
                    // ylo: ori = 1
                    // zlo: ori = 2
                    // xhi: ori = 3
                    // yhi: ori = 4
                    // zhi: ori = 5

                    amrex::Gpu::DeviceVector<BCRec> bcrs_d(ncomp);
#ifdef AMREX_USE_GPU
                    Gpu::htod_memcpy_async
                        (bcrs_d.data(), bcrs.data(), sizeof(BCRec)*ncomp);
#else
                    std::memcpy
                        (bcrs_d.data(), bcrs.data(), sizeof(BCRec)*ncomp);
#endif
                    const amrex::BCRec* bc_ptr = bcrs_d.data();

                    if (m_r2d) fill_from_bndryregs(m_lev, bx, dest_arr, icomp, bccomp, ncomp, domain, bc_ptr, time);

#ifdef ERF_USE_NETCDF
                    if (m_init_type == "real") {
                        int icomp_for_wrfbdy, ncomp_for_wrfbdy, bccomp_for_wrfbdy;
                        if (m_var_idx == Vars::cons) {
                            icomp_for_wrfbdy = RhoTheta_comp;
                            bccomp_for_wrfbdy = BCVars::RhoTheta_bc_comp;
                            ncomp_for_wrfbdy = 1; // (Because we are currently only filling U, V, W, T)
                        } else {
                            icomp_for_wrfbdy = icomp;
                            bccomp_for_wrfbdy = bccomp;
                            ncomp_for_wrfbdy = 1; // (Because we are currently only filling U, V, W, T)
                        }
                        fill_from_wrfbdy(m_lev, bx, dest_arr, icomp_for_wrfbdy, bccomp_for_wrfbdy, ncomp_for_wrfbdy,
                                         domain, bc_ptr,
                                         time, m_bdy_time_interval);
                    }
#endif
                        Gpu::streamSynchronize(); // because of bcrs_d
                } // !gdomain.contains(bx)
            } // MFIter
        } // OpenMP
    } // operator()
