#include "AMReX_PhysBCFunct.H"
#include <ERF_PhysBCFunct.H>

//
// dest_arr is the Array4 to be filled
// time is the time at which the data should be filled
// bccomp is the index into both domain_bcs_type_bcr and bc_extdir_vals
//     so this follows the BCVars enum
//
void ERFPhysBCFunct::impose_yvel_bcs (const Array4<Real>& dest_arr, const Box& bx, const Box& domain,
                                      const Array4<Real const>& z_nd,
                                      const GpuArray<Real,AMREX_SPACEDIM> dxInv,
                                      Real /*time*/, int bccomp)
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

    amrex::Gpu::DeviceVector<BCRec> bcrs_d(ncomp);
#ifdef AMREX_USE_GPU
    Gpu::htod_memcpy_async(bcrs_d.data(), bcrs.data(), sizeof(BCRec)*ncomp);
#else
    std::memcpy(bcrs_d.data(), bcrs.data(), sizeof(BCRec)*ncomp);
#endif
    const amrex::BCRec* bc_ptr = bcrs_d.data();

    amrex::GpuArray<amrex::GpuArray<amrex::Real, AMREX_SPACEDIM*2>,
                                                 AMREX_SPACEDIM+NVAR> l_bc_extdir_vals_d;

    for (int i = 0; i < ncomp; i++)
        for (int ori = 0; ori < 2*AMREX_SPACEDIM; ori++)
            l_bc_extdir_vals_d[i][ori] = m_bc_extdir_vals[bccomp+i][ori];

    ParallelFor(bx, ncomp, [=] AMREX_GPU_DEVICE (int i, int j, int k, int n)
    {
        // Populate ghost cells on lo-x and hi-x domain boundaries
        if (i < dom_lo.x) {
            int iflip = dom_lo.x - 1- i;
            if (bc_ptr[n].lo(0) == ERFBCType::ext_dir) {
                dest_arr(i,j,k) = l_bc_extdir_vals_d[n][0];
            } else if (bc_ptr[n].lo(0) == ERFBCType::foextrap) {
                dest_arr(i,j,k) =  dest_arr(dom_lo.x,j,k);
            } else if (bc_ptr[n].lo(0) == ERFBCType::reflect_even) {
                dest_arr(i,j,k) =  dest_arr(iflip,j,k);
            } else if (bc_ptr[n].lo(0) == ERFBCType::reflect_odd) {
                dest_arr(i,j,k) = -dest_arr(iflip,j,k);
            }

        } else if (i > dom_hi.x) {
            int iflip =  2*dom_hi.x + 1 - i;
            if (bc_ptr[n].hi(0) == ERFBCType::ext_dir) {
                dest_arr(i,j,k) = l_bc_extdir_vals_d[n][3];
            } else if (bc_ptr[n].hi(0) == ERFBCType::foextrap) {
                dest_arr(i,j,k) =  dest_arr(dom_hi.x,j,k);
            } else if (bc_ptr[n].hi(0) == ERFBCType::reflect_even) {
                dest_arr(i,j,k) =  dest_arr(iflip,j,k);
            } else if (bc_ptr[n].hi(0) == ERFBCType::reflect_odd) {
                dest_arr(i,j,k) = -dest_arr(iflip,j,k);
            }
        }

        // Populate ghost cells on lo-y and hi-y domain boundaries
        if (j < dom_lo.y) {
            int jflip = dom_lo.y-j;
            if (bc_ptr[n].lo(1) == ERFBCType::ext_dir) {
                dest_arr(i,j,k) = l_bc_extdir_vals_d[n][1];
            } else if (bc_ptr[n].lo(1) == ERFBCType::foextrap) {
                dest_arr(i,j,k) =  dest_arr(i,dom_lo.y,k);
            } else if (bc_ptr[n].lo(1) == ERFBCType::reflect_even) {
                dest_arr(i,j,k) =  dest_arr(i,jflip,k);
            } else if (bc_ptr[n].lo(1) == ERFBCType::reflect_odd) {
                dest_arr(i,j,k) = -dest_arr(i,jflip,k);
            }
        } else if (j > dom_hi.y+1) {
            int jflip =  2*(dom_hi.y + 1) - j;
            if (bc_ptr[n].hi(1) == ERFBCType::ext_dir) {
                dest_arr(i,j,k) = l_bc_extdir_vals_d[n][4];
            } else if (bc_ptr[n].hi(1) == ERFBCType::foextrap) {
                dest_arr(i,j,k) =  dest_arr(i,dom_hi.y+1,k);
            } else if (bc_ptr[n].hi(1) == ERFBCType::reflect_even) {
                dest_arr(i,j,k) =  dest_arr(i,jflip,k);
            } else if (bc_ptr[n].hi(1) == ERFBCType::reflect_odd) {
                dest_arr(i,j,k) = -dest_arr(i,jflip,k);
            }
        }
    });

    ParallelFor(bx, ncomp, [=] AMREX_GPU_DEVICE (int i, int j, int k, int n)
    {
        // Populate ghost cells on lo-z and hi-z domain boundaries
        if (k < dom_lo.z) {
            int kflip = dom_lo.z - 1 - k;
            if (bc_ptr[n].lo(2) == ERFBCType::ext_dir) {
                dest_arr(i,j,k) = l_bc_extdir_vals_d[n][2];
            } else if (bc_ptr[n].lo(2) == ERFBCType::foextrap) {
                dest_arr(i,j,k) =  dest_arr(i,j,dom_lo.z);
            } else if (bc_ptr[n].lo(2) == ERFBCType::reflect_even) {
                dest_arr(i,j,k) =  dest_arr(i,j,kflip);
            } else if (bc_ptr[n].lo(2) == ERFBCType::reflect_odd) {
                dest_arr(i,j,k) = -dest_arr(i,j,kflip);
            }

        } else if (k > dom_hi.z) {
            int kflip =  2*dom_hi.z + 1 - k;
            if (bc_ptr[n].hi(2) == ERFBCType::ext_dir) {
                dest_arr(i,j,k) = l_bc_extdir_vals_d[n][5];
            } else if (bc_ptr[n].hi(2) == ERFBCType::foextrap) {
                dest_arr(i,j,k) =  dest_arr(i,j,dom_hi.z);
            } else if (bc_ptr[n].hi(2) == ERFBCType::reflect_even) {
                dest_arr(i,j,k) =  dest_arr(i,j,kflip);
            } else if (bc_ptr[n].hi(2) == ERFBCType::reflect_odd) {
                dest_arr(i,j,k) = -dest_arr(i,j,kflip);
            }
        }

        // We only set the values on the domain faces themselves if EXT_DIR
        if (j == dom_lo.y   && bc_ptr[n].lo(1) == ERFBCType::ext_dir)
            dest_arr(i,j,k) = l_bc_extdir_vals_d[n][1];
        if (j == dom_hi.y+1 && bc_ptr[n].lo(4) == ERFBCType::ext_dir)
            dest_arr(i,j,k) = l_bc_extdir_vals_d[n][4];
    });

    if (m_z_phys_nd) {
        const auto&  bx_lo = amrex::lbound(bx);
        const auto&  bx_hi = amrex::ubound(bx);
        // Neumann conditions (d<var>/dn = 0) must be aware of the surface normal with terrain.
        // An additional source term arises from d<var>/dx & d<var>/dy & met_h_xi/eta/zeta.
        //=====================================================================================
        // Only modify scalars, U, or V
        // Loop over each component
        for (int n = 0; n < ncomp; n++) {
            // Hit for Neumann condition at kmin
            if( bcrs[n].lo(2) == ERFBCType::foextrap) {
                // Loop over ghost cells in bottom XY-plane (valid box)
                amrex::Box xybx = bx;
                xybx.setBig(2,-1);
                xybx.setSmall(2,bx.smallEnd()[2]);
                int k0 = 0;

                // Get the dz cell size
                amrex::GeometryData const& geomdata = m_geom.data();
                amrex::Real dz = geomdata.CellSize(2);

                // Fill all the Neumann srcs with terrain
                ParallelFor(xybx, [=] AMREX_GPU_DEVICE (int i, int j, int k)
                {
                    // Clip indices for ghost-cells
                    int ii = amrex::min(amrex::max(i,dom_lo.x),dom_hi.x);
                    int jj = amrex::min(amrex::max(j,dom_lo.y),dom_hi.y);

                    // Get metrics
                    amrex::Real met_h_xi,met_h_eta,met_h_zeta;

                    ComputeMetricAtJface(ii,jj,k0,met_h_xi,met_h_eta,met_h_zeta,dxInv,z_nd,TerrainMet::all);

                    // GradX at IJK location inside domain -- this relies on the assumption that we have
                    // used foextrap for cell-centered quantities outside the domain to define the gradient as zero
                    amrex::Real GradVarx, GradVary;
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
                    amrex::Real met_fac =  met_h_zeta / ( met_h_xi*met_h_xi + met_h_eta*met_h_eta + 1. );

                    // Accumulate in bottom ghost cell (EXTRAP already populated)
                    dest_arr(i,j,k) -= dz * met_fac * ( met_h_xi * GradVarx + met_h_eta * GradVary );
                });
            } // foextrap
        } // ncomp
    } //m_z_phys_nd

    Gpu::streamSynchronize();
}
