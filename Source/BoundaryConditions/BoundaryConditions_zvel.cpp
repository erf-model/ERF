#include "AMReX_PhysBCFunct.H"
#include <ERF_PhysBCFunct.H>
#include <prob_common.H>

//
// dest_arr is the Array4 to be filled
// time is the time at which the data should be filled
// bccomp is the index into both domain_bcs_type_bcr and bc_extdir_vals
//     so this follows the BCVars enum
//
void ERFPhysBCFunct::impose_zvel_bcs (const Array4<Real>& dest_arr, const Box& bx, const Box& domain,
                                      const Array4<Real const>& velx_arr,
                                      const Array4<Real const>& vely_arr,
                                      const Array4<Real const>& z_nd_arr,
                                      const GpuArray<Real,AMREX_SPACEDIM> dx,
                                      const GpuArray<Real,AMREX_SPACEDIM> dxInv,
                                      Real time, Real time_mt, Real delta_t,
                                      int bccomp, int terrain_type)
{
    const auto& dom_lo = amrex::lbound(domain);
    const auto& dom_hi = amrex::ubound(domain);

    // Based on BCRec for the domain, we need to make BCRec for this Box
    // bccomp is used as starting index for m_domain_bcs_type
    //      0 is used as starting index for bcrs
    int ncomp = 1;
    Vector<BCRec> bcrs(ncomp);
    amrex::setBC(bx, domain, bccomp, 0, ncomp, m_domain_bcs_type, bcrs);

    // xlo: ori = 0
    // ylo: ori = 1
    // zlo: ori = 2
    // xhi: ori = 3
    // yhi: ori = 4
    // zhi: ori = 5

    //! Based on BCRec for the domain, we need to make BCRec for this Box
    // bccomp is used as starting index for m_domain_bcs_type
    //      0 is used as starting index for bcrs
    amrex::setBC(bx, domain, bccomp, 0, ncomp, m_domain_bcs_type, bcrs);

    amrex::Gpu::DeviceVector<BCRec> bcrs_d(ncomp);
#ifdef AMREX_USE_GPU
    Gpu::htod_memcpy_async(bcrs_d.data(), bcrs.data(), sizeof(BCRec)*ncomp);
#else
    std::memcpy(bcrs_d.data(), bcrs.data(), sizeof(BCRec)*ncomp);
#endif
    const amrex::BCRec* bc_ptr = bcrs_d.data();

    amrex::GpuArray<amrex::GpuArray<amrex::Real, AMREX_SPACEDIM*2>,
                                                 AMREX_SPACEDIM+NVAR> l_bc_extdir_vals_d;

    bool l_use_terrain = (m_z_phys_nd != nullptr);
    bool l_moving_terrain = (terrain_type == 1);

    for (int i = 0; i < ncomp; i++)
        for (int ori = 0; ori < 2*AMREX_SPACEDIM; ori++)
            l_bc_extdir_vals_d[i][ori] = m_bc_extdir_vals[bccomp+i][ori];

    ParallelFor(bx, ncomp, [=] AMREX_GPU_DEVICE (int i, int j, int k, int n)
    {
        // Lo-x boundary
        if (i < dom_lo.x) {
            int iflip = dom_lo.x - 1 - i;
            if (bc_ptr[n].lo(0) == ERFBCType::ext_dir) {
                dest_arr(i,j,k) = l_bc_extdir_vals_d[n][0];
                if (l_use_terrain) {
                    dest_arr(i,j,k) = WFromOmega(i,j,k,dest_arr(i,j,k),
                                                 velx_arr,vely_arr,z_nd_arr,dxInv);
                }
            } else if (bc_ptr[n].lo(0) == ERFBCType::foextrap) {
                dest_arr(i,j,k) =  dest_arr(dom_lo.x,j,k);
            } else if (bc_ptr[n].lo(0) == ERFBCType::reflect_even) {
                dest_arr(i,j,k) =  dest_arr(iflip,j,k);
            } else if (bc_ptr[n].lo(0) == ERFBCType::reflect_odd) {
                dest_arr(i,j,k) = -dest_arr(iflip,j,k);
            }

        // Hi-x boundary
        } else if (i > dom_hi.x) {
            int iflip = 2*dom_hi.x + 1 - i;
            if (bc_ptr[n].hi(0) == ERFBCType::ext_dir) {
                dest_arr(i,j,k) = l_bc_extdir_vals_d[n][3];
                if (l_use_terrain) {
                    dest_arr(i,j,k) = WFromOmega(i,j,k,dest_arr(i,j,k),
                                                 velx_arr,vely_arr,z_nd_arr,dxInv);
                }
            } else if (bc_ptr[n].hi(0) == ERFBCType::foextrap) {
                dest_arr(i,j,k) =  dest_arr(dom_hi.x,j,k);
            } else if (bc_ptr[n].hi(0) == ERFBCType::reflect_even) {
                dest_arr(i,j,k) =  dest_arr(iflip,j,k);
            } else if (bc_ptr[n].hi(0) == ERFBCType::reflect_odd) {
                dest_arr(i,j,k) = -dest_arr(iflip,j,k);
            }
        }

        // Lo-y boundary
        if (j < dom_lo.y) {
            int jflip = dom_lo.y - 1 - j;
            if (bc_ptr[n].lo(1) == ERFBCType::ext_dir) {
                dest_arr(i,j,k) = l_bc_extdir_vals_d[n][1];
                if (l_use_terrain) {
                    dest_arr(i,j,k) = WFromOmega(i,j,k,dest_arr(i,j,k),
                                                 velx_arr,vely_arr,z_nd_arr,dxInv);
                }
            } else if (bc_ptr[n].lo(1) == ERFBCType::foextrap) {
                dest_arr(i,j,k) =  dest_arr(i,dom_lo.y,k);
            } else if (bc_ptr[n].lo(1) == ERFBCType::reflect_even) {
                dest_arr(i,j,k) =  dest_arr(i,jflip,k);
            } else if (bc_ptr[n].lo(1) == ERFBCType::reflect_odd) {
                dest_arr(i,j,k) = -dest_arr(i,jflip,k);
            }

        // Hi-y boundary
        } else if (j > dom_hi.y) {
            int jflip =  2*dom_hi.y + 1 - j;
            if (bc_ptr[n].hi(1) == ERFBCType::ext_dir) {
                dest_arr(i,j,k) = l_bc_extdir_vals_d[n][4];
                dest_arr(i,j,k) = WFromOmega(i,j,k,dest_arr(i,j,k),
                                                       velx_arr,vely_arr,z_nd_arr,dxInv);
            } else if (bc_ptr[n].hi(1) == ERFBCType::foextrap) {
                dest_arr(i,j,k) =  dest_arr(i,dom_hi.y,k);
            } else if (bc_ptr[n].hi(1) == ERFBCType::reflect_even) {
                dest_arr(i,j,k) =  dest_arr(i,jflip,k);
            } else if (bc_ptr[n].hi(1) == ERFBCType::reflect_odd) {
                dest_arr(i,j,k) = -dest_arr(i,jflip,k);
            }
        }
    });

    ParallelFor(bx, ncomp, [=] AMREX_GPU_DEVICE (int i, int j, int k, int n)
    {
        // Lo-z boundary
        if (k < dom_lo.z) {
            int kflip = dom_lo.z - k;
            if (bc_ptr[n].lo(2) == ERFBCType::ext_dir) {
                dest_arr(i,j,k) = l_bc_extdir_vals_d[n][2];
            } else if (bc_ptr[n].lo(2) == ERFBCType::foextrap) {
                dest_arr(i,j,k) =  dest_arr(i,j,dom_lo.z);
            } else if (bc_ptr[n].lo(2) == ERFBCType::reflect_even) {
                dest_arr(i,j,k) =  dest_arr(i,j,kflip);
            } else if (bc_ptr[n].lo(2) == ERFBCType::reflect_odd) {
                dest_arr(i,j,k) = -dest_arr(i,j,kflip);
            }

        // Hi-z boundary
        } else if (k > dom_hi.z+1) {
            int kflip =  2*(dom_hi.z + 1) - k;
            if (bc_ptr[n].hi(5) == ERFBCType::ext_dir) {
                dest_arr(i,j,k) = l_bc_extdir_vals_d[n][5];
            } else if (bc_ptr[n].hi(5) == ERFBCType::foextrap) {
                dest_arr(i,j,k) =  dest_arr(i,j,dom_hi.z+1);
            } else if (bc_ptr[n].hi(5) == ERFBCType::reflect_even) {
                dest_arr(i,j,k) =  dest_arr(i,j,kflip);
            } else if (bc_ptr[n].hi(5) == ERFBCType::reflect_odd) {
                dest_arr(i,j,k) = -dest_arr(i,j,kflip);
            }
        }

        // Populate face values on z-boundaries themselves only if EXT_DIR
        if (k == dom_lo.z && l_use_terrain && l_moving_terrain) {
            //************************************************************
            // time_mt = time - delta_t/2   delta_t = dt[lev] or dt_stage
            //************************************************************
            Real dhdt_val = dhdt(i,j,dx,time_mt,delta_t);
            dest_arr(i,j,k) = WFromOmega(i,j,k,dhdt_val,
                                         velx_arr,vely_arr,z_nd_arr,dxInv);

        } else if (k == dom_lo.z && bc_ptr[n].lo(2) == ERFBCType::ext_dir) {
            if (l_use_terrain)
                dest_arr(i,j,k) = WFromOmega(i,j,k,l_bc_extdir_vals_d[n][2],
                                             velx_arr,vely_arr,z_nd_arr,dxInv);
            else
               dest_arr(i,j,k) = l_bc_extdir_vals_d[n][2];

        } else if (k == dom_hi.z+1 && bc_ptr[n].hi(2) == ERFBCType::ext_dir) {
            if (l_use_terrain)
                dest_arr(i,j,k) = WFromOmega(i,j,k,l_bc_extdir_vals_d[n][5],
                                             velx_arr,vely_arr,z_nd_arr,dxInv);
            else
                dest_arr(i,j,k) = l_bc_extdir_vals_d[n][5];
        }
    });

    Gpu::streamSynchronize();
}
