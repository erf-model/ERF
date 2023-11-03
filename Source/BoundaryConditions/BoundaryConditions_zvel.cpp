#include "AMReX_PhysBCFunct.H"
#include <ERF_PhysBCFunct.H>
#include <prob_common.H>

using namespace amrex;

/*
 * Impose lateral boundary conditions on z-component of velocity
 *
 * @param[in] dest_arr  Array4 of the quantity to be filled
 * @param[in] bx        box associated with this data
 * @param[in] domain    computational domain
 * @param[in] z_phys_nd height coordinate at nodes
 * @param[in] bccomp    index into m_domain_bcs_type
 */
void ERFPhysBCFunct::impose_lateral_zvel_bcs (const Array4<Real>& dest_arr,
                                              const Array4<Real const>& xvel_arr,
                                              const Array4<Real const>& yvel_arr,
                                              const Box& bx, const Box& domain,
                                              const Array4<Real const>& z_phys_nd,
                                              const GpuArray<Real,AMREX_SPACEDIM> dxInv,
                                              int bccomp)
{
    BL_PROFILE_VAR("impose_lateral_zvel_bcs()",impose_lateral_zvel_bcs);
    const auto& dom_lo = amrex::lbound(domain);
    const auto& dom_hi = amrex::ubound(domain);

    // Based on BCRec for the domain, we need to make BCRec for this Box
    // bccomp is used as starting index for m_domain_bcs_type
    //      0 is used as starting index for bcrs
    int ncomp = 1;
    Vector<BCRec> bcrs_w(1);
    amrex::setBC(bx, domain, bccomp, 0, 1, m_domain_bcs_type, bcrs_w);

    // xlo: ori = 0
    // ylo: ori = 1
    // zlo: ori = 2
    // xhi: ori = 3
    // yhi: ori = 4
    // zhi: ori = 5

    amrex::Gpu::DeviceVector<BCRec> bcrs_w_d(ncomp);
#ifdef AMREX_USE_GPU
    Gpu::htod_memcpy_async(bcrs_w_d.data(), bcrs_w.data(), sizeof(BCRec));
#else
    std::memcpy(bcrs_w_d.data(), bcrs_w.data(), sizeof(BCRec));
#endif
    const amrex::BCRec* bc_ptr_w = bcrs_w_d.data();

    GpuArray<GpuArray<Real, AMREX_SPACEDIM*2>,AMREX_SPACEDIM+NVAR> l_bc_extdir_vals_d;

    bool l_use_terrain = (m_z_phys_nd != nullptr);

    for (int i = 0; i < ncomp; i++)
        for (int ori = 0; ori < 2*AMREX_SPACEDIM; ori++)
            l_bc_extdir_vals_d[i][ori] = m_bc_extdir_vals[bccomp+i][ori];

    GeometryData const& geomdata = m_geom.data();
    bool is_periodic_in_x = geomdata.isPeriodic(0);
    bool is_periodic_in_y = geomdata.isPeriodic(1);

    FArrayBox dhdtfab;

    // First do all ext_dir bcs
    if (!is_periodic_in_x)
    {
        Box bx_xlo(bx);  bx_xlo.setBig  (0,dom_lo.x-1);
        Box bx_xhi(bx);  bx_xhi.setSmall(0,dom_hi.x+1);
        ParallelFor(
            bx_xlo, ncomp, [=] AMREX_GPU_DEVICE (int i, int j, int k, int n) {
                int iflip = dom_lo.x - 1 - i;
                if (bc_ptr_w[n].lo(0) == ERFBCType::ext_dir) {
                    dest_arr(i,j,k) = l_bc_extdir_vals_d[n][0];
                    if (l_use_terrain) {
                        dest_arr(i,j,k) = WFromOmega(i,j,k,dest_arr(i,j,k),xvel_arr,yvel_arr,z_phys_nd,dxInv);
                    }
                } else if (bc_ptr_w[n].lo(0) == ERFBCType::foextrap) {
                    dest_arr(i,j,k) =  dest_arr(dom_lo.x,j,k);
                } else if (bc_ptr_w[n].lo(0) == ERFBCType::reflect_even) {
                    dest_arr(i,j,k) =  dest_arr(iflip,j,k);
                } else if (bc_ptr_w[n].lo(0) == ERFBCType::reflect_odd) {
                    dest_arr(i,j,k) = -dest_arr(iflip,j,k);
                }
            },
            bx_xhi, ncomp, [=] AMREX_GPU_DEVICE (int i, int j, int k, int n) {
                int iflip = 2*dom_hi.x + 1 - i;
                if (bc_ptr_w[n].hi(0) == ERFBCType::ext_dir) {
                    dest_arr(i,j,k) = l_bc_extdir_vals_d[n][3];
                    if (l_use_terrain) {
                        dest_arr(i,j,k) = WFromOmega(i,j,k,dest_arr(i,j,k),xvel_arr,yvel_arr,z_phys_nd,dxInv);
                    }
                } else if (bc_ptr_w[n].hi(0) == ERFBCType::foextrap) {
                    dest_arr(i,j,k) =  dest_arr(dom_hi.x,j,k);
                } else if (bc_ptr_w[n].hi(0) == ERFBCType::reflect_even) {
                    dest_arr(i,j,k) =  dest_arr(iflip,j,k);
                } else if (bc_ptr_w[n].hi(0) == ERFBCType::reflect_odd) {
                    dest_arr(i,j,k) = -dest_arr(iflip,j,k);
                }
            }
        );
    }

    // First do all ext_dir bcs
    if (!is_periodic_in_y)
    {
        Box bx_ylo(bx);  bx_ylo.setBig  (1,dom_lo.y-1);
        Box bx_yhi(bx);  bx_yhi.setSmall(1,dom_hi.y+1);
        ParallelFor(bx_ylo, ncomp, [=] AMREX_GPU_DEVICE (int i, int j, int k, int n) {
            int jflip = dom_lo.y - 1 - j;
            if (bc_ptr_w[n].lo(1) == ERFBCType::ext_dir) {
                dest_arr(i,j,k) = l_bc_extdir_vals_d[n][1];
                if (l_use_terrain) {
                    dest_arr(i,j,k) = WFromOmega(i,j,k,dest_arr(i,j,k),xvel_arr,yvel_arr,z_phys_nd,dxInv);
                }
            } else if (bc_ptr_w[n].lo(1) == ERFBCType::foextrap) {
                dest_arr(i,j,k) =  dest_arr(i,dom_lo.y,k);
            } else if (bc_ptr_w[n].lo(1) == ERFBCType::reflect_even) {
                dest_arr(i,j,k) =  dest_arr(i,jflip,k);
            } else if (bc_ptr_w[n].lo(1) == ERFBCType::reflect_odd) {
                dest_arr(i,j,k) = -dest_arr(i,jflip,k);
            }
        },
        bx_yhi, ncomp, [=] AMREX_GPU_DEVICE (int i, int j, int k, int n) {
            int jflip =  2*dom_hi.y + 1 - j;
            if (bc_ptr_w[n].hi(1) == ERFBCType::ext_dir) {
                dest_arr(i,j,k) = l_bc_extdir_vals_d[n][4];
                dest_arr(i,j,k) = WFromOmega(i,j,k,dest_arr(i,j,k),xvel_arr,yvel_arr,z_phys_nd,dxInv);
            } else if (bc_ptr_w[n].hi(1) == ERFBCType::foextrap) {
                dest_arr(i,j,k) =  dest_arr(i,dom_hi.y,k);
            } else if (bc_ptr_w[n].hi(1) == ERFBCType::reflect_even) {
                dest_arr(i,j,k) =  dest_arr(i,jflip,k);
            } else if (bc_ptr_w[n].hi(1) == ERFBCType::reflect_odd) {
                dest_arr(i,j,k) = -dest_arr(i,jflip,k);
            }
        });
    }

    Gpu::streamSynchronize();
}

/*
 * Impose vertical boundary conditions on z-component of velocity
 *
 * @param[in] dest_arr  Array4 of the quantity to be filled
 * @param[in] bx        box associated with this data
 * @param[in] domain    computational domain
 * @param[in] z_phys_nd height coordinate at nodes
 * @param[in] dxInv     inverse cell size array
 * @param[in] bccomp_u  index into m_domain_bcs_type corresponding to u
 * @param[in] bccomp_v  index into m_domain_bcs_type corresponding to v
 * @param[in] bccomp_w  index into m_domain_bcs_type corresponding to w
 * @param[in] terrain_type if 1 then the terrain is moving; otherwise fixed
 */

void ERFPhysBCFunct::impose_vertical_zvel_bcs (const Array4<Real>& dest_arr,
                                               const Array4<Real const>& xvel_arr,
                                               const Array4<Real const>& yvel_arr,
                                               const Box& bx, const Box& domain,
                                               const Array4<Real const>& z_phys_nd,
                                               const GpuArray<Real,AMREX_SPACEDIM> dxInv,
                                               int bccomp_u, int bccomp_v, int bccomp_w,
                                               TerrainType terrain_type)
{
    BL_PROFILE_VAR("impose_vertical_zvel_bcs()",impose_vertical_zvel_bcs);
    const auto& dom_lo = amrex::lbound(domain);
    const auto& dom_hi = amrex::ubound(domain);

    // xlo: ori = 0
    // ylo: ori = 1
    // zlo: ori = 2
    // xhi: ori = 3
    // yhi: ori = 4
    // zhi: ori = 5

    // Based on BCRec for the domain, we need to make BCRec for this Box
    // bccomp is used as starting index for m_domain_bcs_type
    //      0 is used as starting index for bcrs
    int ncomp = 1;
    Vector<BCRec> bcrs_u(1), bcrs_v(1), bcrs_w(1);
    amrex::setBC(bx, domain, bccomp_u, 0, 1, m_domain_bcs_type, bcrs_u);
    amrex::setBC(bx, domain, bccomp_v, 0, 1, m_domain_bcs_type, bcrs_v);
    amrex::setBC(bx, domain, bccomp_w, 0, 1, m_domain_bcs_type, bcrs_w);

    amrex::Gpu::DeviceVector<BCRec> bcrs_u_d(ncomp);
    amrex::Gpu::DeviceVector<BCRec> bcrs_v_d(ncomp);
    amrex::Gpu::DeviceVector<BCRec> bcrs_w_d(ncomp);
#ifdef AMREX_USE_GPU
    Gpu::htod_memcpy_async(bcrs_u_d.data(), bcrs_u.data(), sizeof(BCRec));
    Gpu::htod_memcpy_async(bcrs_v_d.data(), bcrs_v.data(), sizeof(BCRec));
    Gpu::htod_memcpy_async(bcrs_w_d.data(), bcrs_w.data(), sizeof(BCRec));
#else
    std::memcpy(bcrs_u_d.data(), bcrs_u.data(), sizeof(BCRec));
    std::memcpy(bcrs_v_d.data(), bcrs_v.data(), sizeof(BCRec));
    std::memcpy(bcrs_w_d.data(), bcrs_w.data(), sizeof(BCRec));
#endif
    const amrex::BCRec* bc_ptr_u = bcrs_u_d.data();
    const amrex::BCRec* bc_ptr_v = bcrs_v_d.data();
    const amrex::BCRec* bc_ptr_w = bcrs_w_d.data();

    bool l_use_terrain = (m_z_phys_nd != nullptr);
    bool l_moving_terrain = (terrain_type == TerrainType::Moving);

    GpuArray<GpuArray<Real, AMREX_SPACEDIM*2>,AMREX_SPACEDIM+NVAR> l_bc_extdir_vals_d;

    for (int i = 0; i < ncomp; i++)
        for (int ori = 0; ori < 2*AMREX_SPACEDIM; ori++)
            l_bc_extdir_vals_d[i][ori] = m_bc_extdir_vals[bccomp_w+i][ori];

    {
        Box bx_zlo(bx);  bx_zlo.setBig  (2,dom_lo.z-1);
        Box bx_zhi(bx);  bx_zhi.setSmall(2,dom_hi.z+2);

        // Populate face values on z-boundaries themselves only if EXT_DIR
        Box bx_zlo_face(bx); bx_zlo_face.setSmall(2,dom_lo.z  ); bx_zlo_face.setBig(2,dom_lo.z  );
        Box bx_zhi_face(bx); bx_zhi_face.setSmall(2,dom_hi.z+1); bx_zhi_face.setBig(2,dom_hi.z+1);

        ParallelFor(bx_zlo, [=] AMREX_GPU_DEVICE (int i, int j, int k) {
            int n = 0;
            int kflip = dom_lo.z - k;
            if (bc_ptr_w[n].lo(2) == ERFBCType::ext_dir) {
                dest_arr(i,j,k) = l_bc_extdir_vals_d[n][2];
            } else if (bc_ptr_w[n].lo(2) == ERFBCType::foextrap) {
                dest_arr(i,j,k) =  dest_arr(i,j,dom_lo.z);
            } else if (bc_ptr_w[n].lo(2) == ERFBCType::reflect_even) {
                dest_arr(i,j,k) =  dest_arr(i,j,kflip);
            } else if (bc_ptr_w[n].lo(2) == ERFBCType::reflect_odd) {
                dest_arr(i,j,k) = -dest_arr(i,j,kflip);
            }
        });

        if (l_use_terrain && l_moving_terrain)
        {
            //************************************************************
            // NOTE: z_t depends on the time interval in which it is
            //       evaluated so we can't arbitrarily define it at a
            //       given time, we must specify an interval
            //************************************************************
        } else if (l_use_terrain) {
            ParallelFor(bx_zlo_face, ncomp, [=] AMREX_GPU_DEVICE (int i, int j, int k, int n) {
                if (bc_ptr_w[n].lo(2) == ERFBCType::ext_dir) {
                    if (bc_ptr_u[n].lo(2) == ERFBCType::ext_dir &&
                        bc_ptr_v[n].lo(2) == ERFBCType::ext_dir) {
                        dest_arr(i,j,k) = WFromOmega(i,j,k,l_bc_extdir_vals_d[n][2],xvel_arr,yvel_arr,z_phys_nd,dxInv);

                    } else if (bc_ptr_u[n].lo(2) != ERFBCType::ext_dir &&
                               bc_ptr_v[n].lo(2) != ERFBCType::ext_dir) {
                        dest_arr(i,j,k) = WFromOmega(i,j,k,l_bc_extdir_vals_d[n][2],xvel_arr,yvel_arr,z_phys_nd,dxInv);
                    } else {
#ifndef AMREX_USE_GPU
                       amrex::Abort("Bad combination of boundary conditions");
#endif
                    }
                }
            });
        } else {
            ParallelFor(bx_zlo_face, ncomp, [=] AMREX_GPU_DEVICE (int i, int j, int k, int n) {
                if (bc_ptr_w[n].lo(2) == ERFBCType::ext_dir) {
                   dest_arr(i,j,k) = l_bc_extdir_vals_d[n][2];
                }
            });
        }

        ParallelFor(
          bx_zhi, ncomp, [=] AMREX_GPU_DEVICE (int i, int j, int k, int n) {
            int kflip =  2*(dom_hi.z + 1) - k;
            if (bc_ptr_w[n].hi(5) == ERFBCType::ext_dir) {
                dest_arr(i,j,k) = l_bc_extdir_vals_d[n][5];
            } else if (bc_ptr_w[n].hi(5) == ERFBCType::foextrap) {
                dest_arr(i,j,k) =  dest_arr(i,j,dom_hi.z+1);
            } else if (bc_ptr_w[n].hi(5) == ERFBCType::reflect_even) {
                dest_arr(i,j,k) =  dest_arr(i,j,kflip);
            } else if (bc_ptr_w[n].hi(5) == ERFBCType::reflect_odd) {
                dest_arr(i,j,k) = -dest_arr(i,j,kflip);
            }

          },
          bx_zhi_face, ncomp, [=] AMREX_GPU_DEVICE (int i, int j, int k, int n) {
            if (bc_ptr_w[n].hi(2) == ERFBCType::ext_dir) {
                if (l_use_terrain)
                    dest_arr(i,j,k) = WFromOmega(i,j,k,l_bc_extdir_vals_d[n][5],xvel_arr,yvel_arr,z_phys_nd,dxInv);
                else
                    dest_arr(i,j,k) = l_bc_extdir_vals_d[n][5];
            }
          }
        );
    }
    Gpu::streamSynchronize();
}
