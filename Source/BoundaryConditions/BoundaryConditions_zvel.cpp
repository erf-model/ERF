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
void ERFPhysBCFunct_w::impose_lateral_zvel_bcs (const Array4<Real      >& dest_arr,
                                                const Array4<Real const>& xvel_arr,
                                                const Array4<Real const>& yvel_arr,
                                                const Box& bx, const Box& domain,
                                                const Array4<Real const>& z_phys_nd,
                                                const GpuArray<Real,AMREX_SPACEDIM> dxInv,
                                                int bccomp)
{
    BL_PROFILE_VAR("impose_lateral_zvel_bcs()",impose_lateral_zvel_bcs);
    const auto& dom_lo = lbound(domain);
    const auto& dom_hi = ubound(domain);

    // Based on BCRec for the domain, we need to make BCRec for this Box
    // bccomp is used as starting index for m_domain_bcs_type
    //      0 is used as starting index for bcrs
    int ncomp = 1;
    Vector<BCRec> bcrs_w(1);
    setBC(enclosedCells(bx), domain, bccomp, 0, 1, m_domain_bcs_type, bcrs_w);

    // xlo: ori = 0
    // ylo: ori = 1
    // zlo: ori = 2
    // xhi: ori = 3
    // yhi: ori = 4
    // zhi: ori = 5

    Gpu::DeviceVector<BCRec> bcrs_w_d(ncomp);
    Gpu::copyAsync(Gpu::hostToDevice, bcrs_w.begin(), bcrs_w.end(), bcrs_w_d.begin());
    const BCRec* bc_ptr_w = bcrs_w_d.data();

    GpuArray<GpuArray<Real, AMREX_SPACEDIM*2>,AMREX_SPACEDIM+NVAR_max> l_bc_extdir_vals_d;

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
        Real* zvel_bc_ptr = m_w_bc_data;
        Box bx_xlo(bx);  bx_xlo.setBig  (0,dom_lo.x-1);
        Box bx_xhi(bx);  bx_xhi.setSmall(0,dom_hi.x+1);
        ParallelFor(
            bx_xlo, ncomp, [=] AMREX_GPU_DEVICE (int i, int j, int k, int n) {
                int iflip = dom_lo.x - 1 - i;
                if (bc_ptr_w[n].lo(0) == ERFBCType::ext_dir) {
                    dest_arr(i,j,k) = (zvel_bc_ptr) ? zvel_bc_ptr[k] : l_bc_extdir_vals_d[n][0];
                    if (l_use_terrain) {
                        dest_arr(i,j,k) = WFromOmega(i,j,k,dest_arr(i,j,k),xvel_arr,yvel_arr,z_phys_nd,dxInv);
                    }
                } else if (bc_ptr_w[n].lo(0) == ERFBCType::foextrap) {
                    dest_arr(i,j,k) =  dest_arr(dom_lo.x,j,k);
                } else if (bc_ptr_w[n].lo(0) == ERFBCType::open) {
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
                    dest_arr(i,j,k) = (zvel_bc_ptr) ? zvel_bc_ptr[k] : l_bc_extdir_vals_d[n][3];
                    if (l_use_terrain) {
                        dest_arr(i,j,k) = WFromOmega(i,j,k,dest_arr(i,j,k),xvel_arr,yvel_arr,z_phys_nd,dxInv);
                    }
                } else if (bc_ptr_w[n].hi(0) == ERFBCType::foextrap) {
                    dest_arr(i,j,k) =  dest_arr(dom_hi.x,j,k);
                } else if (bc_ptr_w[n].hi(0) == ERFBCType::open) {
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
        Real* zvel_bc_ptr = m_w_bc_data;
        Box bx_ylo(bx);  bx_ylo.setBig  (1,dom_lo.y-1);
        Box bx_yhi(bx);  bx_yhi.setSmall(1,dom_hi.y+1);
        ParallelFor(bx_ylo, ncomp, [=] AMREX_GPU_DEVICE (int i, int j, int k, int n) {
            int jflip = dom_lo.y - 1 - j;
            if (bc_ptr_w[n].lo(1) == ERFBCType::ext_dir) {
                dest_arr(i,j,k) = (zvel_bc_ptr) ? zvel_bc_ptr[k] : l_bc_extdir_vals_d[n][1];
                if (l_use_terrain) {
                    dest_arr(i,j,k) = WFromOmega(i,j,k,dest_arr(i,j,k),xvel_arr,yvel_arr,z_phys_nd,dxInv);
                }
            } else if (bc_ptr_w[n].lo(1) == ERFBCType::foextrap) {
                dest_arr(i,j,k) =  dest_arr(i,dom_lo.y,k);
            } else if (bc_ptr_w[n].lo(1) == ERFBCType::open) {
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
                dest_arr(i,j,k) = (zvel_bc_ptr) ? zvel_bc_ptr[k] : l_bc_extdir_vals_d[n][4];
                if (l_use_terrain) {
                    dest_arr(i,j,k) = WFromOmega(i,j,k,dest_arr(i,j,k),xvel_arr,yvel_arr,z_phys_nd,dxInv);
                }
            } else if (bc_ptr_w[n].hi(1) == ERFBCType::foextrap) {
                dest_arr(i,j,k) =  dest_arr(i,dom_hi.y,k);
            } else if (bc_ptr_w[n].hi(1) == ERFBCType::open) {
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

void ERFPhysBCFunct_w::impose_vertical_zvel_bcs (const Array4<Real>& dest_arr,
                                                 const Array4<Real const>& xvel_arr,
                                                 const Array4<Real const>& yvel_arr,
                                                 const Box& bx, const Box& domain,
                                                 const Array4<Real const>& z_phys_nd,
                                                 const GpuArray<Real,AMREX_SPACEDIM> dxInv,
                                                 int bccomp_u, int bccomp_v, int bccomp_w,
                                                 TerrainType terrain_type)
{
    BL_PROFILE_VAR("impose_vertical_zvel_bcs()",impose_vertical_zvel_bcs);
    const auto& dom_lo = lbound(domain);
    const auto& dom_hi = ubound(domain);

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
    setBC(enclosedCells(bx), domain, bccomp_u, 0, 1, m_domain_bcs_type, bcrs_u);
    setBC(enclosedCells(bx), domain, bccomp_v, 0, 1, m_domain_bcs_type, bcrs_v);
    setBC(enclosedCells(bx), domain, bccomp_w, 0, 1, m_domain_bcs_type, bcrs_w);

    // We use these for the asserts below
    const BCRec* bc_ptr_u_h = bcrs_u.data();
    const BCRec* bc_ptr_v_h = bcrs_v.data();
    const BCRec* bc_ptr_w_h = bcrs_w.data();

    bool l_use_terrain = (m_z_phys_nd != nullptr);
    bool l_moving_terrain = (terrain_type == TerrainType::Moving);

    GpuArray<GpuArray<Real, AMREX_SPACEDIM*2>,AMREX_SPACEDIM+NVAR_max> l_bc_extdir_vals_d;

    for (int i = 0; i < ncomp; i++) {
        for (int ori = 0; ori < 2*AMREX_SPACEDIM; ori++) {
            l_bc_extdir_vals_d[i][ori] = m_bc_extdir_vals[bccomp_w+i][ori];
        }
    }

    // *******************************************************
    // Bottom boundary
    // *******************************************************

    // *******************************************************
    // Moving terrain
    // *******************************************************
    if (l_use_terrain && l_moving_terrain)
    {
        //************************************************************
        // NOTE: z_t depends on the time interval in which it is
        //       evaluated so we can't arbitrarily define it at a
        //       given time, we must specify an interval
        //************************************************************

    // Static terrain
    } else if (l_use_terrain) {

        if (m_lev == 0) {
            AMREX_ALWAYS_ASSERT( (bc_ptr_u_h[0].lo(2) == ERFBCType::ext_dir && bc_ptr_v_h[0].lo(2) == ERFBCType::ext_dir) ||
                                 (bc_ptr_u_h[0].lo(2) != ERFBCType::ext_dir && bc_ptr_v_h[0].lo(2) != ERFBCType::ext_dir) );
        } else {
            // If we do not reach to the top or bottom boundary then the z-vel should be
            //    filled by interpolation from the coarser grid using ERF_FillPatcher.
        }
        if (bx.smallEnd(2) == dom_lo.z) {
            ParallelFor(makeSlab(bx,2,dom_lo.z), [=] AMREX_GPU_DEVICE (int i, int j, int k) {
                dest_arr(i,j,k) = WFromOmega(i,j,k,l_bc_extdir_vals_d[0][2],xvel_arr,yvel_arr,z_phys_nd,dxInv);
            });
        }

    // No terrain
    } else {
        if (bx.smallEnd(2) == dom_lo.z) {
            ParallelFor(makeSlab(bx,2,dom_lo.z), [=] AMREX_GPU_DEVICE (int i, int j, int k) {
               dest_arr(i,j,k) = l_bc_extdir_vals_d[0][2];
            });
        }
    }

    // *******************************************************
    // Top boundary
    // *******************************************************

    // NOTE: if we set SlipWall at top, that generates ERFBCType::ext_dir which sets w=0 here
    // NOTE: if we set  Outflow at top, that generates ERFBCType::foextrap which doesn't touch w here
    if (bx.bigEnd(2) == dom_hi.z+1) {
        if (bc_ptr_w_h[0].hi(2) == ERFBCType::ext_dir) {
            ParallelFor(makeSlab(bx,2,dom_hi.z+1), [=] AMREX_GPU_DEVICE (int i, int j, int k)
            {
                if (l_use_terrain) {
                    dest_arr(i,j,k) = WFromOmega(i,j,k,l_bc_extdir_vals_d[0][5],xvel_arr,yvel_arr,z_phys_nd,dxInv);
                } else {
                    dest_arr(i,j,k) = l_bc_extdir_vals_d[0][5];
                }
            });
        } else if (bc_ptr_w_h[0].hi(2) == ERFBCType::neumann_int) {
            ParallelFor(makeSlab(bx,2,dom_hi.z+1), [=] AMREX_GPU_DEVICE (int i, int j, int k)
            {
                dest_arr(i,j,k) = (4.0*dest_arr(i,j,dom_hi.z) - dest_arr(i,j,dom_hi.z-1))/3.0;
            });
        }
    }
    Gpu::streamSynchronize();
}

/*
 * Impose lateral boundary conditions on z-component of velocity
 * ASSUMING NO TERRAIN
 *
 * @param[in] dest_arr  Array4 of the quantity to be filled
 * @param[in] bx        box associated with this data
 * @param[in] domain    computational domain
 * @param[in] bccomp    index into m_domain_bcs_type
 */
void ERFPhysBCFunct_w_no_terrain::impose_lateral_zvel_bcs (const Array4<Real      >& dest_arr,
                                                           const Box& bx, const Box& domain,
                                                           int bccomp)
{
    BL_PROFILE_VAR("impose_lateral_zvel_bcs()",impose_lateral_zvel_bcs);
    const auto& dom_lo = lbound(domain);
    const auto& dom_hi = ubound(domain);

    // Based on BCRec for the domain, we need to make BCRec for this Box
    // bccomp is used as starting index for m_domain_bcs_type
    //      0 is used as starting index for bcrs
    int ncomp = 1;
    Vector<BCRec> bcrs_w(1);
    setBC(enclosedCells(bx), domain, bccomp, 0, 1, m_domain_bcs_type, bcrs_w);

    // xlo: ori = 0
    // ylo: ori = 1
    // zlo: ori = 2
    // xhi: ori = 3
    // yhi: ori = 4
    // zhi: ori = 5

    Gpu::DeviceVector<BCRec> bcrs_w_d(ncomp);
    Gpu::copyAsync(Gpu::hostToDevice, bcrs_w.begin(), bcrs_w.end(), bcrs_w_d.begin());
    const BCRec* bc_ptr_w = bcrs_w_d.data();

    GpuArray<GpuArray<Real, AMREX_SPACEDIM*2>,AMREX_SPACEDIM+NVAR_max> l_bc_extdir_vals_d;

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
                } else if (bc_ptr_w[n].lo(0) == ERFBCType::foextrap) {
                    dest_arr(i,j,k) =  dest_arr(dom_lo.x,j,k);
                } else if (bc_ptr_w[n].lo(0) == ERFBCType::open) {
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
                } else if (bc_ptr_w[n].hi(0) == ERFBCType::foextrap) {
                    dest_arr(i,j,k) =  dest_arr(dom_hi.x,j,k);
                } else if (bc_ptr_w[n].hi(0) == ERFBCType::open) {
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
            } else if (bc_ptr_w[n].lo(1) == ERFBCType::foextrap) {
                dest_arr(i,j,k) =  dest_arr(i,dom_lo.y,k);
            } else if (bc_ptr_w[n].lo(1) == ERFBCType::open) {
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
            } else if (bc_ptr_w[n].hi(1) == ERFBCType::foextrap) {
                dest_arr(i,j,k) =  dest_arr(i,dom_hi.y,k);
            } else if (bc_ptr_w[n].hi(1) == ERFBCType::open) {
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
 * @param[in] bccomp_w  index into m_domain_bcs_type corresponding to w
 */

void ERFPhysBCFunct_w_no_terrain::impose_vertical_zvel_bcs (const Array4<Real>& dest_arr,
                                                            const Box& bx, const Box& domain,
                                                            int bccomp_w)
{
    BL_PROFILE_VAR("impose_vertical_zvel_bcs()",impose_vertical_zvel_bcs);
    const auto& dom_lo = lbound(domain);
    const auto& dom_hi = ubound(domain);

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
    Vector<BCRec> bcrs_w(1);
    setBC(enclosedCells(bx), domain, bccomp_w, 0, 1, m_domain_bcs_type, bcrs_w);

    // We use these for the asserts below
    const BCRec* bc_ptr_w_h = bcrs_w.data();

    GpuArray<GpuArray<Real, AMREX_SPACEDIM*2>,AMREX_SPACEDIM+NVAR_max> l_bc_extdir_vals_d;

    for (int i = 0; i < ncomp; i++) {
        for (int ori = 0; ori < 2*AMREX_SPACEDIM; ori++) {
            l_bc_extdir_vals_d[i][ori] = m_bc_extdir_vals[bccomp_w+i][ori];
        }
    }

    // *******************************************************
    // Bottom boundary
    // *******************************************************

    if (bx.smallEnd(2) == dom_lo.z) {
        AMREX_ALWAYS_ASSERT(bc_ptr_w_h[0].lo(2) == ERFBCType::ext_dir);
        ParallelFor(makeSlab(bx,2,dom_lo.z), [=] AMREX_GPU_DEVICE (int i, int j, int k) {
           dest_arr(i,j,k) = l_bc_extdir_vals_d[0][2];
        });
    }

    // *******************************************************
    // Top boundary
    // *******************************************************

    // NOTE: if we set SlipWall at top, that generates ERFBCType::ext_dir which sets w=0 here
    // NOTE: if we set  Outflow at top, that generates ERFBCType::foextrap which doesn't touch w here
    if (bx.bigEnd(2) == dom_hi.z+1) {
        AMREX_ALWAYS_ASSERT(bc_ptr_w_h[0].hi(2) == ERFBCType::ext_dir ||
                            bc_ptr_w_h[0].hi(2) == ERFBCType::neumann_int);
        if (bc_ptr_w_h[0].hi(2) == ERFBCType::ext_dir) {
            ParallelFor(makeSlab(bx,2,dom_hi.z+1), [=] AMREX_GPU_DEVICE (int i, int j, int k)
            {
                dest_arr(i,j,k) = l_bc_extdir_vals_d[0][5];
            });
        } else if (bc_ptr_w_h[0].hi(2) == ERFBCType::neumann_int) {
            ParallelFor(makeSlab(bx,2,dom_hi.z+1), [=] AMREX_GPU_DEVICE (int i, int j, int k)
            {
                dest_arr(i,j,k) = (4.0*dest_arr(i,j,dom_hi.z) - dest_arr(i,j,dom_hi.z-1))/3.0;
            });
        }
    }
    Gpu::streamSynchronize();
}
