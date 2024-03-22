#include "AMReX_PhysBCFunct.H"
#include <ERF_PhysBCFunct.H>

using namespace amrex;

/*
 * Impose lateral boundary conditions on y-component of velocity
 *
 * @param[in] dest_arr Array4 of the quantity to be filled
 * @param[in] bx       box associated with this data
 * @param[in] domain   computational domain
 * @param[in] bccomp   index into m_domain_bcs_type
 */
void ERFPhysBCFunct_v::impose_lateral_yvel_bcs (const Array4<Real>& dest_arr,
                                                const Box& bx, const Box& domain,
                                                int bccomp)
{
    BL_PROFILE_VAR("impose_lateral_yvel_bcs()",impose_lateral_yvel_bcs);
    const auto& dom_lo = lbound(domain);
    const auto& dom_hi = ubound(domain);

    // Based on BCRec for the domain, we need to make BCRec for this Box
    // bccomp is used as starting index for m_domain_bcs_type
    //      0 is used as starting index for bcrs
    int ncomp = 1;
    Vector<BCRec> bcrs(ncomp);
    setBC(bx, domain, bccomp, 0, ncomp, m_domain_bcs_type, bcrs);

    // xlo: ori = 0
    // ylo: ori = 1
    // zlo: ori = 2
    // xhi: ori = 3
    // yhi: ori = 4
    // zhi: ori = 5

    Gpu::DeviceVector<BCRec> bcrs_d(ncomp);
#ifdef AMREX_USE_GPU
    Gpu::htod_memcpy_async(bcrs_d.data(), bcrs.data(), sizeof(BCRec)*ncomp);
#else
    std::memcpy(bcrs_d.data(), bcrs.data(), sizeof(BCRec)*ncomp);
#endif
    const BCRec* bc_ptr = bcrs_d.data();

    GpuArray<GpuArray<Real, AMREX_SPACEDIM*2>, AMREX_SPACEDIM+NVAR_max> l_bc_extdir_vals_d;

    for (int i = 0; i < ncomp; i++)
        for (int ori = 0; ori < 2*AMREX_SPACEDIM; ori++)
            l_bc_extdir_vals_d[i][ori] = m_bc_extdir_vals[bccomp+i][ori];

    GeometryData const& geomdata = m_geom.data();
    bool is_periodic_in_x = geomdata.isPeriodic(0);
    bool is_periodic_in_y = geomdata.isPeriodic(1);

    // First do all ext_dir bcs
    if (!is_periodic_in_x)
    {
        // Populate ghost cells on lo-x and hi-x domain boundaries
        Box bx_xlo(bx);  bx_xlo.setBig  (0,dom_lo.x-1);
        Box bx_xhi(bx);  bx_xhi.setSmall(0,dom_hi.x+1);
        ParallelFor(
            bx_xlo, ncomp, [=] AMREX_GPU_DEVICE (int i, int j, int k, int n) {
                int iflip = dom_lo.x - 1- i;
                if (bc_ptr[n].lo(0) == ERFBCType::ext_dir) {
                    dest_arr(i,j,k) = l_bc_extdir_vals_d[n][0];
                } else if (bc_ptr[n].lo(0) == ERFBCType::foextrap) {
                    dest_arr(i,j,k) =  dest_arr(dom_lo.x,j,k);
                } else if (bc_ptr[n].lo(0) == ERFBCType::open) {
                    dest_arr(i,j,k) =  dest_arr(dom_lo.x,j,k);
                } else if (bc_ptr[n].lo(0) == ERFBCType::reflect_even) {
                    dest_arr(i,j,k) =  dest_arr(iflip,j,k);
                } else if (bc_ptr[n].lo(0) == ERFBCType::reflect_odd) {
                    dest_arr(i,j,k) = -dest_arr(iflip,j,k);
                }
            },
            bx_xhi, ncomp, [=] AMREX_GPU_DEVICE (int i, int j, int k, int n) {
                int iflip =  2*dom_hi.x + 1 - i;
                if (bc_ptr[n].hi(0) == ERFBCType::ext_dir) {
                    dest_arr(i,j,k) = l_bc_extdir_vals_d[n][3];
                } else if (bc_ptr[n].hi(0) == ERFBCType::foextrap) {
                    dest_arr(i,j,k) =  dest_arr(dom_hi.x,j,k);
                } else if (bc_ptr[n].hi(0) == ERFBCType::open) {
                    dest_arr(i,j,k) =  dest_arr(dom_hi.x,j,k);
                } else if (bc_ptr[n].hi(0) == ERFBCType::reflect_even) {
                    dest_arr(i,j,k) =  dest_arr(iflip,j,k);
                } else if (bc_ptr[n].hi(0) == ERFBCType::reflect_odd) {
                    dest_arr(i,j,k) = -dest_arr(iflip,j,k);
                }
            }
        );
    }

    if (!is_periodic_in_y)
    {
        // Populate ghost cells on lo-y and hi-y domain boundaries
        Box bx_ylo(bx);  bx_ylo.setBig  (1,dom_lo.y-1);
        Box bx_yhi(bx);  bx_yhi.setSmall(1,dom_hi.y+2);
        Box bx_ylo_face(bx); bx_ylo_face.setSmall(1,dom_lo.y  ); bx_ylo_face.setBig(1,dom_lo.y  );
        Box bx_yhi_face(bx); bx_yhi_face.setSmall(1,dom_hi.y+1); bx_yhi_face.setBig(1,dom_hi.y+1);

        ParallelFor(
            bx_ylo, ncomp, [=] AMREX_GPU_DEVICE (int i, int j, int k, int n)
            {
                int jflip = dom_lo.y-j;
                if (bc_ptr[n].lo(1) == ERFBCType::ext_dir) {
                    dest_arr(i,j,k) = l_bc_extdir_vals_d[n][1];
                } else if (bc_ptr[n].lo(1) == ERFBCType::foextrap) {
                    dest_arr(i,j,k) =  dest_arr(i,dom_lo.y,k);
                } else if (bc_ptr[n].lo(1) == ERFBCType::open) {
                    dest_arr(i,j,k) =  dest_arr(i,dom_lo.y,k);
                } else if (bc_ptr[n].lo(1) == ERFBCType::reflect_even) {
                    dest_arr(i,j,k) =  dest_arr(i,jflip,k);
                } else if (bc_ptr[n].lo(1) == ERFBCType::reflect_odd) {
                    dest_arr(i,j,k) = -dest_arr(i,jflip,k);
                } else if (bc_ptr[n].lo(1) == ERFBCType::neumann_int) {
                    dest_arr(i,j,k) = (4.0*dest_arr(i,dom_lo.y+1,k) - dest_arr(i,dom_lo.y+2,k))/3.0;
                }
            },
            // We only set the values on the domain faces themselves if EXT_DIR
            bx_ylo_face, ncomp, [=] AMREX_GPU_DEVICE (int i, int j, int k, int n)
            {
                if (bc_ptr[n].lo(1) == ERFBCType::ext_dir) {
                    dest_arr(i,j,k) = l_bc_extdir_vals_d[n][1];
                } else if (bc_ptr[n].lo(1) == ERFBCType::neumann_int) {
                    dest_arr(i,j,k) = (4.0*dest_arr(i,dom_lo.y+1,k) - dest_arr(i,dom_lo.y+2,k))/3.0;
                }
            }
        );
        ParallelFor(
            bx_yhi, ncomp, [=] AMREX_GPU_DEVICE (int i, int j, int k, int n)
            {
                 int jflip =  2*(dom_hi.y + 1) - j;
                 if (bc_ptr[n].hi(1) == ERFBCType::ext_dir) {
                     dest_arr(i,j,k) = l_bc_extdir_vals_d[n][4];
                 } else if (bc_ptr[n].hi(1) == ERFBCType::foextrap) {
                     dest_arr(i,j,k) =  dest_arr(i,dom_hi.y+1,k);
                 } else if (bc_ptr[n].hi(1) == ERFBCType::open) {
                     dest_arr(i,j,k) =  dest_arr(i,dom_hi.y+1,k);
                 } else if (bc_ptr[n].hi(1) == ERFBCType::reflect_even) {
                     dest_arr(i,j,k) =  dest_arr(i,jflip,k);
                 } else if (bc_ptr[n].hi(1) == ERFBCType::reflect_odd) {
                     dest_arr(i,j,k) = -dest_arr(i,jflip,k);
                } else if (bc_ptr[n].hi(1) == ERFBCType::neumann_int) {
                    dest_arr(i,j,k) = (4.0*dest_arr(i,dom_hi.y,k) - dest_arr(i,dom_hi.y-1,k))/3.0;
                }
            },
            // We only set the values on the domain faces themselves if EXT_DIR
            bx_yhi_face, ncomp, [=] AMREX_GPU_DEVICE (int i, int j, int k, int n)
            {
                if (bc_ptr[n].hi(1) == ERFBCType::ext_dir) {
                    dest_arr(i,j,k) = l_bc_extdir_vals_d[n][4];
                } else if (bc_ptr[n].hi(1) == ERFBCType::neumann_int) {
                    dest_arr(i,j,k) = (4.0*dest_arr(i,dom_hi.y,k) - dest_arr(i,dom_hi.y-1,k))/3.0;
                }
            }
        );
    }
    Gpu::streamSynchronize();
}

/*
 * Impose vertical boundary conditions on y-component of velocity
 *
 * @param[in] dest_arr  the Array4 of the quantity to be filled
 * @param[in] bx        the box associated with this data
 * @param[in] domain    the computational domain
 * @param[in] z_phys_nd height coordinate at nodes
 * @param[in] dxInv     inverse cell size array
 * @param[in] bccomp    index into m_domain_bcs_type
 */

void ERFPhysBCFunct_v::impose_vertical_yvel_bcs (const Array4<Real>& dest_arr,
                                                 const Box& bx, const Box& domain,
                                                 const Array4<Real const>& z_phys_nd,
                                                 const GpuArray<Real,AMREX_SPACEDIM> dxInv,
                                                 int bccomp)
{
    BL_PROFILE_VAR("impose_vertical_yvel_bcs()",impose_vertical_yvel_bcs);
    const auto& dom_lo = lbound(domain);
    const auto& dom_hi = ubound(domain);

    // Based on BCRec for the domain, we need to make BCRec for this Box
    // bccomp is used as starting index for m_domain_bcs_type
    //      0 is used as starting index for bcrs
    int ncomp = 1;
    Vector<BCRec> bcrs(ncomp);
    setBC(bx, domain, bccomp, 0, ncomp, m_domain_bcs_type, bcrs);

    // xlo: ori = 0
    // ylo: ori = 1
    // zlo: ori = 2
    // xhi: ori = 3
    // yhi: ori = 4
    // zhi: ori = 5

    Gpu::DeviceVector<BCRec> bcrs_d(ncomp);
#ifdef AMREX_USE_GPU
    Gpu::htod_memcpy_async(bcrs_d.data(), bcrs.data(), sizeof(BCRec)*ncomp);
#else
    std::memcpy(bcrs_d.data(), bcrs.data(), sizeof(BCRec)*ncomp);
#endif
    const BCRec* bc_ptr = bcrs_d.data();

    GpuArray<GpuArray<Real, AMREX_SPACEDIM*2>, AMREX_SPACEDIM+NVAR_max> l_bc_extdir_vals_d;

    for (int i = 0; i < ncomp; i++)
        for (int ori = 0; ori < 2*AMREX_SPACEDIM; ori++)
            l_bc_extdir_vals_d[i][ori] = m_bc_extdir_vals[bccomp+i][ori];

    GeometryData const& geomdata = m_geom.data();

    {
        // Populate ghost cells on lo-z and hi-z domain boundaries
        Box bx_zlo(bx);  bx_zlo.setBig  (2,dom_lo.z-1);
        Box bx_zhi(bx);  bx_zhi.setSmall(2,dom_hi.z+1);
        ParallelFor(
            bx_zlo, ncomp, [=] AMREX_GPU_DEVICE (int i, int j, int k, int n) {
                int kflip = dom_lo.z - 1 - k;
                if (bc_ptr[n].lo(2) == ERFBCType::ext_dir) {
                    dest_arr(i,j,k) = l_bc_extdir_vals_d[n][2];
                } else if (bc_ptr[n].lo(2) == ERFBCType::foextrap) {
                    dest_arr(i,j,k) =  dest_arr(i,j,dom_lo.z);
                } else if (bc_ptr[n].lo(2) == ERFBCType::open) {
                    dest_arr(i,j,k) =  dest_arr(i,j,dom_lo.z);
                } else if (bc_ptr[n].lo(2) == ERFBCType::reflect_even) {
                    dest_arr(i,j,k) =  dest_arr(i,j,kflip);
                } else if (bc_ptr[n].lo(2) == ERFBCType::reflect_odd) {
                    dest_arr(i,j,k) = -dest_arr(i,j,kflip);
                }
            },
            bx_zhi, ncomp, [=] AMREX_GPU_DEVICE (int i, int j, int k, int n) {
                int kflip =  2*dom_hi.z + 1 - k;
                if (bc_ptr[n].hi(2) == ERFBCType::ext_dir) {
                    dest_arr(i,j,k) = l_bc_extdir_vals_d[n][5];
                } else if (bc_ptr[n].hi(2) == ERFBCType::foextrap) {
                    dest_arr(i,j,k) =  dest_arr(i,j,dom_hi.z);
                } else if (bc_ptr[n].hi(2) == ERFBCType::open) {
                    dest_arr(i,j,k) =  dest_arr(i,j,dom_hi.z);
                } else if (bc_ptr[n].hi(2) == ERFBCType::reflect_even) {
                    dest_arr(i,j,k) =  dest_arr(i,j,kflip);
                } else if (bc_ptr[n].hi(2) == ERFBCType::reflect_odd) {
                    dest_arr(i,j,k) = -dest_arr(i,j,kflip);
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
            if(bcrs[n].lo(2) == ERFBCType::foextrap) {
                // Loop over ghost cells in bottom XY-plane (valid box)
                Box xybx = bx;
                xybx.setBig(2,-1);
                xybx.setSmall(2,bx.smallEnd()[2]);
                int k0 = 0;

                // Get the dz cell size
                Real dz = geomdata.CellSize(2);

                // Fill all the Neumann srcs with terrain
                ParallelFor(xybx, [=] AMREX_GPU_DEVICE (int i, int j, int k)
                {
                    // Clip indices for ghost-cells
                    int ii = amrex::min(amrex::max(i,dom_lo.x),dom_hi.x);
                    int jj = amrex::min(amrex::max(j,dom_lo.y),dom_hi.y);

                    // Get metrics
                    Real met_h_xi   = Compute_h_xi_AtJface  (ii, jj, k0, dxInv, z_phys_nd);
                    Real met_h_eta  = Compute_h_eta_AtJface (ii, jj, k0, dxInv, z_phys_nd);
                    Real met_h_zeta = Compute_h_zeta_AtJface(ii, jj, k0, dxInv, z_phys_nd);

                    // GradX at IJK location inside domain -- this relies on the assumption that we have
                    // used foextrap for cell-centered quantities outside the domain to define the gradient as zero
                    Real GradVarx, GradVary;
                    if (i < dom_lo.x-1 || i > dom_hi.x+1)
                        GradVarx = 0.0;
                    else if (i+1 > bx_hi.x)
                        GradVarx =       dxInv[0] * (dest_arr(i  ,j,k0) - dest_arr(i-1,j,k0));
                    else if (i-1 < bx_lo.x)
                        GradVarx =       dxInv[0] * (dest_arr(i+1,j,k0) - dest_arr(i  ,j,k0));
                    else
                        GradVarx = 0.5 * dxInv[0] * (dest_arr(i+1,j,k0) - dest_arr(i-1,j,k0));

                    // GradY at IJK location inside domain -- this relies on the assumption that we have
                    // used foextrap for cell-centered quantities outside the domain to define the gradient as zero
                    if (j < dom_lo.y-1 || j > dom_hi.y+1)
                        GradVary = 0.0;
                    else if (j+1 > bx_hi.y)
                        GradVary =       dxInv[1] * (dest_arr(i,j  ,k0) - dest_arr(i,j-1,k0));
                    else if (j-1 < bx_lo.y)
                        GradVary =       dxInv[1] * (dest_arr(i,j+1,k0) - dest_arr(i,j  ,k0));
                    else
                        GradVary = 0.5 * dxInv[1] * (dest_arr(i,j+1,k0) - dest_arr(i,j-1,k0));

                    // Prefactor
                    Real met_fac =  met_h_zeta / ( met_h_xi*met_h_xi + met_h_eta*met_h_eta + 1. );

                    // Accumulate in bottom ghost cell (EXTRAP already populated)
                    dest_arr(i,j,k) -= dz * met_fac * ( met_h_xi * GradVarx + met_h_eta * GradVary );
                });
            } // foextrap
        } // ncomp
    } //m_z_phys_nd
    Gpu::streamSynchronize();
}
