#include "AMReX_PhysBCFunct.H"
#include <ERF_PhysBCFunct.H>

using namespace amrex;

/*
 * Impose lateral boundary conditions on conserved scalars (at cell centers)
 *
 * @param[in,out] dest_arr cell-centered data to be filled
 * @param[in]     bx       box holding data to be filled
 * @param[in]     domain   simulation domain
 * @param[in]     icomp    index into the MultiFab -- this can be any value from 0 to NVAR-1
 * @param[in]     ncomp    the number of components -- this can be any value from 1 to NVAR
 *                         as long as icomp+ncomp <= NVAR-1.
 */

void ERFPhysBCFunct_cons::impose_lateral_cons_bcs (const Array4<Real>& dest_arr, const Box& bx, const Box& domain,
                                                   int icomp, int ncomp, int ngz)
{
    BL_PROFILE_VAR("impose_lateral_cons_bcs()",impose_lateral_cons_bcs);
    const auto& dom_lo = lbound(domain);
    const auto& dom_hi = ubound(domain);

    // xlo: ori = 0
    // ylo: ori = 1
    // zlo: ori = 2
    // xhi: ori = 3
    // yhi: ori = 4
    // zhi: ori = 5

    // Based on BCRec for the domain, we need to make BCRec for this Box
    //      0 is used as starting index for bcrs
    Vector<BCRec> bcrs(ncomp);

    GpuArray<GpuArray<Real, AMREX_SPACEDIM*2>,NBCVAR_max> l_bc_extdir_vals_d;

    const int* bxlo = bx.loVect();
    const int* bxhi = bx.hiVect();
    const int* dlo  = domain.loVect();
    const int* dhi  = domain.hiVect();

    for (int nc = 0; nc < ncomp; nc++)
    {
        int bc_comp = (icomp+nc >= RhoScalar_comp && icomp+nc < RhoScalar_comp+NSCALARS) ?
                       BCVars::RhoScalar_bc_comp : icomp+nc;
        if (bc_comp > BCVars::RhoScalar_bc_comp) bc_comp -= (NSCALARS-1);
        for (int dir = 0; dir < AMREX_SPACEDIM; dir++)
        {
            bcrs[nc].setLo(dir, ( bxlo[dir]<=dlo[dir]
                                 ? m_domain_bcs_type[bc_comp].lo(dir) : BCType::int_dir ));
            bcrs[nc].setHi(dir, ( bxhi[dir]>=dhi[dir]
                                 ? m_domain_bcs_type[bc_comp].hi(dir) : BCType::int_dir ));
        }

        for (int ori = 0; ori < 2*AMREX_SPACEDIM; ori++) {
            l_bc_extdir_vals_d[bc_comp][ori]  = m_bc_extdir_vals[bc_comp][ori];
        }
    }

    Gpu::DeviceVector<BCRec> bcrs_d(ncomp);
    Gpu::copyAsync(Gpu::hostToDevice, bcrs.begin(), bcrs.end(), bcrs_d.begin());
    const BCRec* bc_ptr = bcrs_d.data();

    GeometryData const& geomdata = m_geom.data();
    bool is_periodic_in_x = geomdata.isPeriodic(0);
    bool is_periodic_in_y = geomdata.isPeriodic(1);

    // First do all ext_dir bcs
    if (!is_periodic_in_x)
    {
        Box bx_xlo(bx);  bx_xlo.setBig  (0,dom_lo.x-1);
        Box bx_xhi(bx);  bx_xhi.setSmall(0,dom_hi.x+1);
        ParallelFor(
            bx_xlo, ncomp, [=] AMREX_GPU_DEVICE (int i, int j, int k, int n)
            {
                int dest_comp = icomp+n;
                int bc_comp   = (dest_comp >= RhoScalar_comp && dest_comp < RhoScalar_comp+NSCALARS) ?
                                 BCVars::RhoScalar_bc_comp : dest_comp;
                if (bc_comp > BCVars::RhoScalar_bc_comp) bc_comp -= (NSCALARS-1);
                int l_bc_type = bc_ptr[n].lo(0);

                if (l_bc_type == ERFBCType::ext_dir) {
                    dest_arr(i,j,k,dest_comp) = l_bc_extdir_vals_d[bc_comp][0];
                } else if (l_bc_type == ERFBCType::ext_dir_prim) {
                    Real rho = dest_arr(dom_lo.x,j,k,Rho_comp);
                    dest_arr(i,j,k,dest_comp) = rho * l_bc_extdir_vals_d[bc_comp][0];
                }
            },
            bx_xhi, ncomp, [=] AMREX_GPU_DEVICE (int i, int j, int k, int n)
            {
                int dest_comp = icomp+n;
                int bc_comp   = (dest_comp >= RhoScalar_comp && dest_comp < RhoScalar_comp+NSCALARS) ?
                                 BCVars::RhoScalar_bc_comp : dest_comp;
                if (bc_comp > BCVars::RhoScalar_bc_comp) bc_comp -= (NSCALARS-1);
                int h_bc_type = bc_ptr[n].hi(0);

                if (h_bc_type == ERFBCType::ext_dir) {
                    dest_arr(i,j,k,dest_comp) = l_bc_extdir_vals_d[bc_comp][3];
                } else if (h_bc_type == ERFBCType::ext_dir_prim) {
                    Real rho = dest_arr(dom_hi.x,j,k,Rho_comp);
                    dest_arr(i,j,k,dest_comp) = rho * l_bc_extdir_vals_d[bc_comp][3];
                }
            }
        );
    }

    if (!is_periodic_in_y)
    {
        Box bx_ylo(bx);  bx_ylo.setBig  (1,dom_lo.y-1);
        Box bx_yhi(bx);  bx_yhi.setSmall(1,dom_hi.y+1);
        ParallelFor(
            bx_ylo, ncomp, [=] AMREX_GPU_DEVICE (int i, int j, int k, int n)
            {
                int dest_comp = icomp+n;
                int bc_comp   = (dest_comp >= RhoScalar_comp && dest_comp < RhoScalar_comp+NSCALARS) ?
                                 BCVars::RhoScalar_bc_comp : dest_comp;
                if (bc_comp > BCVars::RhoScalar_bc_comp) bc_comp -= (NSCALARS-1);
                int l_bc_type = bc_ptr[n].lo(1);
                if (l_bc_type == ERFBCType::ext_dir) {
                    dest_arr(i,j,k,dest_comp) = l_bc_extdir_vals_d[bc_comp][1];
                } else if (l_bc_type == ERFBCType::ext_dir_prim) {
                    Real rho = dest_arr(i,dom_lo.y,k,Rho_comp);
                    dest_arr(i,j,k,dest_comp) = rho * l_bc_extdir_vals_d[bc_comp][1];
                }
            },
            bx_yhi, ncomp, [=] AMREX_GPU_DEVICE (int i, int j, int k, int n)
            {
                int dest_comp = icomp+n;
                int bc_comp   = (dest_comp >= RhoScalar_comp && dest_comp < RhoScalar_comp+NSCALARS) ?
                                 BCVars::RhoScalar_bc_comp : dest_comp;
                if (bc_comp > BCVars::RhoScalar_bc_comp) bc_comp -= (NSCALARS-1);
                int h_bc_type = bc_ptr[n].hi(1);
                if (h_bc_type == ERFBCType::ext_dir) {
                    dest_arr(i,j,k,dest_comp) = l_bc_extdir_vals_d[bc_comp][4];
                } else if (h_bc_type == ERFBCType::ext_dir_prim) {
                    Real rho = dest_arr(i,dom_hi.y,k,Rho_comp);
                    dest_arr(i,j,k,dest_comp) = rho * l_bc_extdir_vals_d[bc_comp][4];
                }
            }
        );
    }

    // Next do ghost cells in x-direction but not reaching out in y
    // The corners we miss here will be covered in the y-loop below or by periodicity
    if (!is_periodic_in_x)
    {
        // Populate ghost cells on lo-x and hi-x domain boundaries
        Box bx_xlo(bx);  bx_xlo.setBig  (0,dom_lo.x-1);
        Box bx_xhi(bx);  bx_xhi.setSmall(0,dom_hi.x+1);
        if (bx_xlo.smallEnd(2) != domain.smallEnd(2)) bx_xlo.growLo(2,ngz);
        if (bx_xlo.bigEnd(2)   != domain.bigEnd(2))   bx_xlo.growHi(2,ngz);
        if (bx_xhi.smallEnd(2) != domain.smallEnd(2)) bx_xhi.growLo(2,ngz);
        if (bx_xhi.bigEnd(2)   != domain.bigEnd(2))   bx_xhi.growHi(2,ngz);
        ParallelFor(
            bx_xlo, ncomp, [=] AMREX_GPU_DEVICE (int i, int j, int k, int n)
            {
                int dest_comp = icomp+n;
                int l_bc_type = bc_ptr[n].lo(0);
                int iflip = dom_lo.x - 1 - i;
                if (l_bc_type == ERFBCType::foextrap) {
                    dest_arr(i,j,k,dest_comp) =  dest_arr(dom_lo.x,j,k,dest_comp);
                } else if (l_bc_type == ERFBCType::open) {
                    dest_arr(i,j,k,dest_comp) =  dest_arr(dom_lo.x,j,k,dest_comp);
                } else if (l_bc_type == ERFBCType::reflect_even) {
                    dest_arr(i,j,k,dest_comp) =  dest_arr(iflip,j,k,dest_comp);
                } else if (l_bc_type == ERFBCType::reflect_odd) {
                    dest_arr(i,j,k,dest_comp) = -dest_arr(iflip,j,k,dest_comp);
                } else if (l_bc_type == ERFBCType::hoextrapcc) {
                    dest_arr(i,j,k,dest_comp) = 2.0*dest_arr(dom_lo.x,j,k,dest_comp) - dest_arr(dom_lo.x+1,j,k,dest_comp) ;
                }
            },
            bx_xhi, ncomp, [=] AMREX_GPU_DEVICE (int i, int j, int k, int n)
            {
                int dest_comp = icomp+n;
                int h_bc_type = bc_ptr[n].hi(0);
                int iflip =  2*dom_hi.x + 1 - i;
                if (h_bc_type == ERFBCType::foextrap) {
                    dest_arr(i,j,k,dest_comp) =  dest_arr(dom_hi.x,j,k,dest_comp);
                } else if (h_bc_type == ERFBCType::open) {
                    dest_arr(i,j,k,dest_comp) =  dest_arr(dom_hi.x,j,k,dest_comp);
                } else if (h_bc_type == ERFBCType::reflect_even) {
                    dest_arr(i,j,k,dest_comp) =  dest_arr(iflip,j,k,dest_comp);
                } else if (h_bc_type == ERFBCType::reflect_odd) {
                    dest_arr(i,j,k,dest_comp) = -dest_arr(iflip,j,k,dest_comp);
                } else if (h_bc_type == ERFBCType::hoextrapcc) {
                    dest_arr(i,j,k,dest_comp) = 2.0*dest_arr(dom_hi.x,j,k,dest_comp) - dest_arr(dom_hi.x-1,j,k,dest_comp) ;
                }
            }
        );
    }

    if (!is_periodic_in_y)
    {
        // Populate ghost cells on lo-y and hi-y domain boundaries
        Box bx_ylo(bx);  bx_ylo.setBig  (1,dom_lo.y-1);
        Box bx_yhi(bx);  bx_yhi.setSmall(1,dom_hi.y+1);
        if (bx_ylo.smallEnd(2) != domain.smallEnd(2)) bx_ylo.growLo(2,ngz);
        if (bx_ylo.bigEnd(2)   != domain.bigEnd(2))   bx_ylo.growHi(2,ngz);
        if (bx_yhi.smallEnd(2) != domain.smallEnd(2)) bx_yhi.growLo(2,ngz);
        if (bx_yhi.bigEnd(2)   != domain.bigEnd(2))   bx_yhi.growHi(2,ngz);
        ParallelFor(
            bx_ylo, ncomp, [=] AMREX_GPU_DEVICE (int i, int j, int k, int n)
            {
                int dest_comp = icomp+n;
                int l_bc_type = bc_ptr[n].lo(1);
                int jflip = dom_lo.y - 1 - j;
                if (l_bc_type == ERFBCType::foextrap) {
                    dest_arr(i,j,k,dest_comp) =  dest_arr(i,dom_lo.y,k,dest_comp);
                } else if (l_bc_type == ERFBCType::open) {
                    dest_arr(i,j,k,dest_comp) =  dest_arr(i,dom_lo.y,k,dest_comp);
                } else if (l_bc_type == ERFBCType::reflect_even) {
                    dest_arr(i,j,k,dest_comp) =  dest_arr(i,jflip,k,dest_comp);
                } else if (l_bc_type == ERFBCType::reflect_odd) {
                    dest_arr(i,j,k,dest_comp) = -dest_arr(i,jflip,k,dest_comp);
                } else if (l_bc_type == ERFBCType::hoextrapcc) {
                    dest_arr(i,j,k,dest_comp) = 2.0*dest_arr(i,dom_lo.y,k,dest_comp) - dest_arr(i,dom_lo.y+1,k,dest_comp) ;
                }

            },
            bx_yhi, ncomp, [=] AMREX_GPU_DEVICE (int i, int j, int k, int n)
            {
                int dest_comp = icomp+n;
                int h_bc_type = bc_ptr[n].hi(1);
                int jflip =  2*dom_hi.y + 1 - j;
                if (h_bc_type == ERFBCType::foextrap) {
                    dest_arr(i,j,k,dest_comp) =  dest_arr(i,dom_hi.y,k,dest_comp);
                } else if (h_bc_type == ERFBCType::open) {
                    dest_arr(i,j,k,dest_comp) =  dest_arr(i,dom_hi.y,k,dest_comp);
                } else if (h_bc_type == ERFBCType::reflect_even) {
                    dest_arr(i,j,k,dest_comp) =  dest_arr(i,jflip,k,dest_comp);
                } else if (h_bc_type == ERFBCType::reflect_odd) {
                    dest_arr(i,j,k,dest_comp) = -dest_arr(i,jflip,k,dest_comp);
                } else if (h_bc_type == ERFBCType::hoextrapcc) {
                    dest_arr(i,j,k,dest_comp) = 2.0*dest_arr(i,dom_hi.y,k,dest_comp) - dest_arr(i,dom_hi.y-1,k,dest_comp);
                }
            }
        );
    }
    Gpu::streamSynchronize();
}

/*
 * Impose vertical boundary conditions on conserved scalars (at cell centers)
 *
 * @param[in] dest_arr  the Array4 of the quantity to be filled
 * @param[in] bx        the box associated with this data
 * @param[in] domain    the computational domain
 * @param[in] z_phys_nd height coordinate at nodes
 * @param[in] dxInv     inverse cell size array
 * @param[in] icomp     the index of the first component to be filled
 * @param[in] ncomp     the number of components -- this can be any value from 1 to NVAR
 *                      as long as icomp+ncomp <= NVAR-1.
 */

void ERFPhysBCFunct_cons::impose_vertical_cons_bcs (const Array4<Real>& dest_arr, const Box& bx, const Box& domain,
                                                    const Array4<Real const>& z_phys_nd,
                                                    const GpuArray<Real,AMREX_SPACEDIM> dxInv,
                                                    int icomp, int ncomp)
{
    BL_PROFILE_VAR("impose_lateral_cons_bcs()",impose_lateral_cons_bcs);
    const auto& dom_lo = lbound(domain);
    const auto& dom_hi = ubound(domain);

    Box per_grown_domain(domain);
    int growx = (m_geom.isPeriodic(0)) ? 1 : 0;
    int growy = (m_geom.isPeriodic(1)) ? 1 : 0;
    per_grown_domain.grow(IntVect(growx,growy,0));
    const auto& perdom_lo = lbound(per_grown_domain);
    const auto& perdom_hi = ubound(per_grown_domain);

    GeometryData const& geomdata = m_geom.data();

    // xlo: ori = 0
    // ylo: ori = 1
    // zlo: ori = 2
    // xhi: ori = 3
    // yhi: ori = 4
    // zhi: ori = 5

    // Based on BCRec for the domain, we need to make BCRec for this Box
    //      0 is used as starting index for bcrs
    Vector<BCRec> bcrs(ncomp);
    GpuArray<GpuArray<Real, AMREX_SPACEDIM*2>,NBCVAR_max> l_bc_extdir_vals_d;
    GpuArray<GpuArray<Real, AMREX_SPACEDIM*2>,NBCVAR_max> l_bc_neumann_vals_d;

    const int* bxlo = bx.loVect();
    const int* bxhi = bx.hiVect();
    const int* dlo  = domain.loVect();
    const int* dhi  = domain.hiVect();

    for (int nc = 0; nc < ncomp; nc++)
    {
        int bc_comp = (icomp+nc >= RhoScalar_comp && icomp+nc < RhoScalar_comp+NSCALARS) ?
                       BCVars::RhoScalar_bc_comp : icomp+nc;
        if (bc_comp > BCVars::RhoScalar_bc_comp) bc_comp -= (NSCALARS-1);
        for (int dir = 0; dir < AMREX_SPACEDIM; dir++)
        {
            bcrs[nc].setLo(dir, ( bxlo[dir]<=dlo[dir]
                                 ? m_domain_bcs_type[bc_comp].lo(dir) : BCType::int_dir ));
            bcrs[nc].setHi(dir, ( bxhi[dir]>=dhi[dir]
                                 ? m_domain_bcs_type[bc_comp].hi(dir) : BCType::int_dir ));
        }

        for (int ori = 0; ori < 2*AMREX_SPACEDIM; ori++) {
            l_bc_extdir_vals_d[bc_comp][ori]  = m_bc_extdir_vals[bc_comp][ori];
            l_bc_neumann_vals_d[bc_comp][ori] = m_bc_neumann_vals[bc_comp][ori];
        }
    }

    Gpu::DeviceVector<BCRec> bcrs_d(icomp+ncomp);
    Gpu::copyAsync(Gpu::hostToDevice, bcrs.begin(), bcrs.end(), bcrs_d.begin());
    const BCRec* bc_ptr = bcrs_d.data();

    {
        Box bx_zlo(bx);  bx_zlo.setBig  (2,dom_lo.z-1);
        Box bx_zhi(bx);  bx_zhi.setSmall(2,dom_hi.z+1);
        ParallelFor(
            bx_zlo, ncomp, [=] AMREX_GPU_DEVICE (int i, int j, int k, int n)
            {
                int dest_comp = icomp+n;
                int bc_comp   = (dest_comp >= RhoScalar_comp && dest_comp < RhoScalar_comp+NSCALARS) ?
                                 BCVars::RhoScalar_bc_comp : dest_comp;
                if (bc_comp > BCVars::RhoScalar_bc_comp) bc_comp -= (NSCALARS-1);
                int l_bc_type = bc_ptr[n].lo(2);
                if (l_bc_type == ERFBCType::ext_dir) {
                    dest_arr(i,j,k,dest_comp) = l_bc_extdir_vals_d[bc_comp][2];
                } else if (l_bc_type == ERFBCType::ext_dir_prim) {
                    Real rho = dest_arr(i,j,dom_lo.z,Rho_comp);
                    dest_arr(i,j,k,dest_comp) = rho * l_bc_extdir_vals_d[bc_comp][2];
                }
            },
            bx_zhi, ncomp, [=] AMREX_GPU_DEVICE (int i, int j, int k, int n)
            {
                int dest_comp = icomp+n;
                int bc_comp   = (dest_comp >= RhoScalar_comp && dest_comp < RhoScalar_comp+NSCALARS) ?
                                 BCVars::RhoScalar_bc_comp : dest_comp;
                if (bc_comp > BCVars::RhoScalar_bc_comp) bc_comp -= (NSCALARS-1);
                int h_bc_type = bc_ptr[n].hi(2);
                if (h_bc_type == ERFBCType::ext_dir) {
                    dest_arr(i,j,k,dest_comp) = l_bc_extdir_vals_d[bc_comp][5];
                } else if (h_bc_type == ERFBCType::ext_dir_prim) {
                    Real rho = dest_arr(i,j,dom_hi.z,Rho_comp);
                    dest_arr(i,j,k,dest_comp) = rho * l_bc_extdir_vals_d[bc_comp][5];
                }

            }
        );
    }

    {
        Box bx_zlo(bx);  bx_zlo.setBig  (2,dom_lo.z-1);
        Box bx_zhi(bx);  bx_zhi.setSmall(2,dom_hi.z+1);
        // Populate ghost cells on lo-z and hi-z domain boundaries
        ParallelFor(
            bx_zlo, ncomp, [=] AMREX_GPU_DEVICE (int i, int j, int k, int n)
            {
                int dest_comp = icomp+n;
                int bc_comp   = (dest_comp >= RhoScalar_comp && dest_comp < RhoScalar_comp+NSCALARS) ?
                                 BCVars::RhoScalar_bc_comp : dest_comp;
                if (bc_comp > BCVars::RhoScalar_bc_comp) bc_comp -= (NSCALARS-1);
                int l_bc_type = bc_ptr[n].lo(2);
                int kflip = dom_lo.z - 1 - i;
                if (l_bc_type == ERFBCType::foextrap) {
                    dest_arr(i,j,k,dest_comp) =  dest_arr(i,j,dom_lo.z,dest_comp);
                } else if (l_bc_type == ERFBCType::open) {
                    dest_arr(i,j,k,dest_comp) =  dest_arr(i,j,dom_lo.z,dest_comp);
                } else if (l_bc_type == ERFBCType::reflect_even) {
                    dest_arr(i,j,k,dest_comp) =  dest_arr(i,j,kflip,dest_comp);
                } else if (l_bc_type == ERFBCType::reflect_odd) {
                    dest_arr(i,j,k,dest_comp) = -dest_arr(i,j,kflip,dest_comp);
                } else if (l_bc_type == ERFBCType::neumann) {
                    Real delta_z = (dom_lo.z - k) / dxInv[2];
                    dest_arr(i,j,k,dest_comp) = dest_arr(i,j,dom_lo.z,dest_comp) -
                        delta_z*l_bc_neumann_vals_d[bc_comp][2]*dest_arr(i,j,dom_lo.z,Rho_comp);
                } else if (l_bc_type == ERFBCType::hoextrapcc) {
                    dest_arr(i,j,k,dest_comp) = 2.0*dest_arr(i,j,dom_lo.z,dest_comp) - dest_arr(i,j,dom_lo.z+1,dest_comp);
                }
            },
            bx_zhi, ncomp, [=] AMREX_GPU_DEVICE (int i, int j, int k, int n)
            {
                int dest_comp = icomp+n;
                int bc_comp   = (dest_comp >= RhoScalar_comp && dest_comp < RhoScalar_comp+NSCALARS) ?
                                 BCVars::RhoScalar_bc_comp : dest_comp;
                if (bc_comp > BCVars::RhoScalar_bc_comp) bc_comp -= (NSCALARS-1);
                int h_bc_type = bc_ptr[n].hi(2);
                int kflip =  2*dom_hi.z + 1 - i;
                if (h_bc_type == ERFBCType::foextrap) {
                    dest_arr(i,j,k,dest_comp) =  dest_arr(i,j,dom_hi.z,dest_comp);
                } else if (h_bc_type == ERFBCType::open) {
                    dest_arr(i,j,k,dest_comp) =  dest_arr(i,j,dom_hi.z,dest_comp);
                } else if (h_bc_type == ERFBCType::reflect_even) {
                    dest_arr(i,j,k,dest_comp) =  dest_arr(i,j,kflip,dest_comp);
                } else if (h_bc_type == ERFBCType::reflect_odd) {
                    dest_arr(i,j,k,dest_comp) = -dest_arr(i,j,kflip,dest_comp);
                } else if (h_bc_type == ERFBCType::neumann) {
                    Real delta_z = (k - dom_hi.z) / dxInv[2];
                    if( (icomp+n) == Rho_comp ) {
                        dest_arr(i,j,k,dest_comp) = dest_arr(i,j,dom_hi.z,dest_comp) +
                            delta_z*l_bc_neumann_vals_d[bc_comp][5];
                    } else {
                        dest_arr(i,j,k,dest_comp) = dest_arr(i,j,dom_hi.z,dest_comp) +
                            delta_z*l_bc_neumann_vals_d[bc_comp][5]*dest_arr(i,j,dom_hi.z,Rho_comp);
                    }
                } else if (h_bc_type == ERFBCType::hoextrapcc){
                    dest_arr(i,j,k,dest_comp) = 2.0*dest_arr(i,j,dom_hi.z,dest_comp) - dest_arr(i,j,dom_hi.z-1,dest_comp);
                }
            }
        );
    }

    if (m_z_phys_nd) {
        const auto&  bx_lo = lbound(bx);
        const auto&  bx_hi = ubound(bx);

        // Neumann conditions (d<var>/dn = 0) must be aware of the surface normal with terrain.
        // An additional source term arises from d<var>/dx & d<var>/dy & met_h_xi/eta/zeta.
        //=====================================================================================
        // Only modify scalars, U, or V
        // Loop over each component
        for (int n = 0; n < ncomp; n++) {
            // Hit for Neumann condition at kmin
            int dest_comp = icomp+n;
            int l_bc_type = bc_ptr[n].lo(2);
            if(l_bc_type == ERFBCType::foextrap)
            {
                // Loop over ghost cells in bottom XY-plane (valid box)
                Box xybx = bx;

                int k0 = 0;
                if (xybx.smallEnd(2) < 0) {

                    xybx.setBig(2,dom_lo.z-1);
                    xybx.setSmall(2,bx.smallEnd()[2]);

                    // Get the dz cell size
                    Real dz = geomdata.CellSize(2);

                    // Fill all the Neumann srcs with terrain
                    ParallelFor(xybx, [=] AMREX_GPU_DEVICE (int i, int j, int k)
                    {
                        // Clip indices for ghost-cells
                        int ii = amrex::min(amrex::max(i,perdom_lo.x),perdom_hi.x);
                        int jj = amrex::min(amrex::max(j,perdom_lo.y),perdom_hi.y);

                        // Get metrics
                        Real met_h_xi   = Compute_h_xi_AtCellCenter  (ii,jj,k0,dxInv,z_phys_nd);
                        Real met_h_eta  = Compute_h_eta_AtCellCenter (ii,jj,k0,dxInv,z_phys_nd);
                        Real met_h_zeta = Compute_h_zeta_AtCellCenter(ii,jj,k0,dxInv,z_phys_nd);

                        // GradX at IJK location inside domain -- this relies on the assumption that we have
                        // used foextrap for cell-centered quantities outside the domain to define the gradient as zero
                        Real GradVarx, GradVary;
                        if (i < dom_lo.x-1 || i > dom_hi.x+1) {
                            GradVarx = 0.0;
                        } else if (i+1 > bx_hi.x) {
                            GradVarx =       dxInv[0] * (dest_arr(i  ,j,k0,dest_comp) - dest_arr(i-1,j,k0,dest_comp));
                        } else if (i-1 < bx_lo.x) {
                            GradVarx =       dxInv[0] * (dest_arr(i+1,j,k0,dest_comp) - dest_arr(i  ,j,k0,dest_comp));
                        } else {
                            GradVarx = 0.5 * dxInv[0] * (dest_arr(i+1,j,k0,dest_comp) - dest_arr(i-1,j,k0,dest_comp));
                        }

                        // GradY at IJK location inside domain -- this relies on the assumption that we have
                        // used foextrap for cell-centered quantities outside the domain to define the gradient as zero
                        if (j < dom_lo.y-1 || j > dom_hi.y+1) {
                            GradVary = 0.0;
                        } else if (j+1 > bx_hi.y) {
                            GradVary =       dxInv[1] * (dest_arr(i,j  ,k0,dest_comp) - dest_arr(i,j-1,k0,dest_comp));
                        } else if (j-1 < bx_lo.y) {
                            GradVary =       dxInv[1] * (dest_arr(i,j+1,k0,dest_comp) - dest_arr(i,j  ,k0,dest_comp));
                        } else {
                            GradVary = 0.5 * dxInv[1] * (dest_arr(i,j+1,k0,dest_comp) - dest_arr(i,j-1,k0,dest_comp));
                        }

                        // Prefactor
                        Real met_fac =  met_h_zeta / ( met_h_xi*met_h_xi + met_h_eta*met_h_eta + 1. );

                        // Accumulate in bottom ghost cell (EXTRAP already populated)
                        dest_arr(i,j,k,dest_comp) -= dz * met_fac * ( met_h_xi * GradVarx + met_h_eta * GradVary );
                    });
                } // box includes k0
            } // foextrap
        } // ncomp
    } // m_z_phys_nd
    Gpu::streamSynchronize();
}
