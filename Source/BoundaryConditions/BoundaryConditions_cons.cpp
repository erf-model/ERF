#include "AMReX_PhysBCFunct.H"
#include <ERF_PhysBCFunct.H>

using namespace amrex;

//
// mf is the multifab to be filled
// icomp is the index into the MultiFab -- if cell-centered this can be any value
//       from 0 to NVAR-1, if face-centered this must be 0
// ncomp is the number of components -- if cell-centered (var_idx = 0) this can be any value
//       from 1 to NVAR as long as icomp+ncomp <= NVAR-1.  If face-centered this
//       must be 1
// time is the time at which the data should be filled
// bccomp is the index into both domain_bcs_type_bcr and bc_extdir_vals for icomp = 0  --
//     so this follows the BCVars enum
//
void ERFPhysBCFunct::impose_cons_bcs (const Array4<Real>& dest_arr, const Box& bx, const Box& domain,
                                      const Array4<Real const>& z_nd,
                                      const GpuArray<Real,AMREX_SPACEDIM> dxInv,
                                      int icomp, int ncomp, Real /*time*/, int bccomp)
{
    const auto& dom_lo = amrex::lbound(domain);
    const auto& dom_hi = amrex::ubound(domain);

    // Based on BCRec for the domain, we need to make BCRec for this Box
    // bccomp is used as starting index for m_domain_bcs_type
    //      0 is used as starting index for bcrs
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

    GpuArray<GpuArray<Real, AMREX_SPACEDIM*2>,AMREX_SPACEDIM+NVAR> l_bc_extdir_vals_d;

    for (int i = 0; i < ncomp; i++)
        for (int ori = 0; ori < 2*AMREX_SPACEDIM; ori++)
            l_bc_extdir_vals_d[i][ori] = m_bc_extdir_vals[bccomp+i][ori];

    // First do all ext_dir bcs
    ParallelFor(bx, ncomp, [=] AMREX_GPU_DEVICE (int i, int j, int k, int n)
    {
        if (i < dom_lo.x) {
            if (bc_ptr[n].lo(0) == ERFBCType::ext_dir) {
                dest_arr(i,j,k,icomp+n) = l_bc_extdir_vals_d[n][0];
            }
        } else if (i > dom_hi.x) {
            if (bc_ptr[n].hi(0) == ERFBCType::ext_dir) {
                dest_arr(i,j,k,icomp+n) = l_bc_extdir_vals_d[n][3];
            }
        }
        if (j < dom_lo.y) {
            if (bc_ptr[n].lo(1) == ERFBCType::ext_dir) {
                dest_arr(i,j,k,icomp+n) = l_bc_extdir_vals_d[n][1];
            }
        } else if (j > dom_hi.y) {
            if (bc_ptr[n].hi(1) == ERFBCType::ext_dir) {
                dest_arr(i,j,k,icomp+n) = l_bc_extdir_vals_d[n][4];
            }
        }

        if (k < dom_lo.z) {
            if (bc_ptr[n].lo(2) == ERFBCType::ext_dir) {
                dest_arr(i,j,k,icomp+n) = l_bc_extdir_vals_d[n][2];
            }
        } else if (k > dom_hi.z) {
            if (bc_ptr[n].hi(2) == ERFBCType::ext_dir) {
                dest_arr(i,j,k,icomp+n) = l_bc_extdir_vals_d[n][5];
            }
        }
    });

    // Next do ghost cells in x-direction but not reaching out in y
    // The corners we miss here will be covered in the y-loop below or by periodicity
    ParallelFor(bx, ncomp, [=] AMREX_GPU_DEVICE (int i, int j, int k, int n)
    {
        // Populate ghost cells on lo-x and hi-x domain boundaries
        if (i < dom_lo.x && k >= dom_lo.z && k <= dom_hi.z) {
            int iflip = dom_lo.x - 1 - i;
            if (bc_ptr[n].lo(0) == ERFBCType::foextrap) {
                dest_arr(i,j,k,icomp+n) =  dest_arr(dom_lo.x,j,k,icomp+n);
            } else if (bc_ptr[n].lo(0) == ERFBCType::reflect_even) {
                dest_arr(i,j,k,icomp+n) =  dest_arr(iflip,j,k,icomp+n);
            } else if (bc_ptr[n].lo(0) == ERFBCType::reflect_odd) {
                dest_arr(i,j,k,icomp+n) = -dest_arr(iflip,j,k,icomp+n);
            }

        } else if (i > dom_hi.x && k >= dom_lo.z && k <= dom_hi.z) {
            int iflip =  2*dom_hi.x + 1 - i;
            if (bc_ptr[n].hi(0) == ERFBCType::foextrap) {
                dest_arr(i,j,k,icomp+n) =  dest_arr(dom_hi.x,j,k,icomp+n);
            } else if (bc_ptr[n].hi(0) == ERFBCType::reflect_even) {
                dest_arr(i,j,k,icomp+n) =  dest_arr(iflip,j,k,icomp+n);
            } else if (bc_ptr[n].hi(0) == ERFBCType::reflect_odd) {
                dest_arr(i,j,k,icomp+n) = -dest_arr(iflip,j,k,icomp+n);
            }
        }

        // Populate ghost cells on lo-y and hi-y domain boundaries
        if (j < dom_lo.y && k >= dom_lo.z && k <= dom_hi.z) {
            int jflip = dom_lo.y - 1 - j;
            if (bc_ptr[n].lo(1) == ERFBCType::foextrap) {
                dest_arr(i,j,k,icomp+n) =  dest_arr(i,dom_lo.y,k,icomp+n);
            } else if (bc_ptr[n].lo(1) == ERFBCType::reflect_even) {
                dest_arr(i,j,k,icomp+n) =  dest_arr(i,jflip,k,icomp+n);
            } else if (bc_ptr[n].lo(1) == ERFBCType::reflect_odd) {
                dest_arr(i,j,k,icomp+n) = -dest_arr(i,jflip,k,icomp+n);
            }
        } else if (j > dom_hi.y && k >= dom_lo.z && k <= dom_hi.z) {
            int jflip =  2*dom_hi.y + 1 - j;
            if (bc_ptr[n].hi(1) == ERFBCType::foextrap) {
                dest_arr(i,j,k,icomp+n) =  dest_arr(i,dom_hi.y,k,icomp+n);
            } else if (bc_ptr[n].hi(1) == ERFBCType::reflect_even) {
                dest_arr(i,j,k,icomp+n) =  dest_arr(i,jflip,k,icomp+n);
            } else if (bc_ptr[n].hi(1) == ERFBCType::reflect_odd) {
                dest_arr(i,j,k,icomp+n) = -dest_arr(i,jflip,k,icomp+n);
            }
        }
    });

    ParallelFor(bx, ncomp, [=] AMREX_GPU_DEVICE (int i, int j, int k, int n)
    {
        // Populate ghost cells on lo-z and hi-z domain boundaries
        if (k < dom_lo.z) {
            int kflip = dom_lo.z - 1 - i;
            if (bc_ptr[n].lo(2) == ERFBCType::foextrap) {
                dest_arr(i,j,k,icomp+n) =  dest_arr(i,j,dom_lo.z,icomp+n);
            } else if (bc_ptr[n].lo(2) == ERFBCType::reflect_even) {
                dest_arr(i,j,k,icomp+n) =  dest_arr(i,j,kflip,icomp+n);
            } else if (bc_ptr[n].lo(2) == ERFBCType::reflect_odd) {
                dest_arr(i,j,k,icomp+n) = -dest_arr(i,j,kflip,icomp+n);
            }

        } else if (k > dom_hi.z) {
            int kflip =  2*dom_hi.z + 1 - i;
            if (bc_ptr[n].hi(2) == ERFBCType::foextrap) {
                dest_arr(i,j,k,icomp+n) =  dest_arr(i,j,dom_hi.z,icomp+n);
            } else if (bc_ptr[n].hi(2) == ERFBCType::reflect_even) {
                dest_arr(i,j,k,icomp+n) =  dest_arr(i,j,kflip,icomp+n);
            } else if (bc_ptr[n].hi(2) == ERFBCType::reflect_odd) {
                dest_arr(i,j,k,icomp+n) = -dest_arr(i,j,kflip,icomp+n);
            }
        }
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
                GeometryData const& geomdata = m_geom.data();
                Real dz = geomdata.CellSize(2);

                // Fill all the Neumann srcs with terrain
                ParallelFor(xybx, [=] AMREX_GPU_DEVICE (int i, int j, int k)
                {
                    // Clip indices for ghost-cells
                    int ii = amrex::min(amrex::max(i,dom_lo.x),dom_hi.x);
                    int jj = amrex::min(amrex::max(j,dom_lo.y),dom_hi.y);

                    // Get metrics
                    Real met_h_xi,met_h_eta,met_h_zeta;

                    ComputeMetricAtCellCenter(ii,jj,k0,met_h_xi,met_h_eta,met_h_zeta,dxInv,z_nd,TerrainMet::all);

                    // GradX at IJK location inside domain -- this relies on the assumption that we have
                    // used foextrap for cell-centered quantities outside the domain to define the gradient as zero
                    Real GradVarx, GradVary;
                    if (i < dom_lo.x-1 || i > dom_hi.x+1)
                        GradVarx = 0.0;
                    else if (i+1 > bx_hi.x)
                        GradVarx =       dxInv[0] * (dest_arr(i  ,j,k0,n) - dest_arr(i-1,j,k0,n));
                    else if (i-1 < bx_lo.x)
                        GradVarx =       dxInv[0] * (dest_arr(i+1,j,k0,n) - dest_arr(i  ,j,k0,n));
                    else
                        GradVarx = 0.5 * dxInv[0] * (dest_arr(i+1,j,k0,n) - dest_arr(i-1,j,k0,n));

                    // GradY at IJK location inside domain -- this relies on the assumption that we have
                    // used foextrap for cell-centered quantities outside the domain to define the gradient as zero
                    if (j < dom_lo.y-1 || j > dom_hi.y+1)
                        GradVary = 0.0;
                    else if (j+1 > bx_hi.y)
                        GradVary =       dxInv[1] * (dest_arr(i,j  ,k0,n) - dest_arr(i,j-1,k0,n));
                    else if (j-1 < bx_lo.y)
                        GradVary =       dxInv[1] * (dest_arr(i,j+1,k0,n) - dest_arr(i,j  ,k0,n));
                    else
                        GradVary = 0.5 * dxInv[1] * (dest_arr(i,j+1,k0,n) - dest_arr(i,j-1,k0,n));

                    // Prefactor
                    Real met_fac =  met_h_zeta / ( met_h_xi*met_h_xi + met_h_eta*met_h_eta + 1. );

                    // Accumulate in bottom ghost cell (EXTRAP already populated)
                    dest_arr(i,j,k,n) -= dz * met_fac * ( met_h_xi * GradVarx + met_h_eta * GradVary );
                });
            } // foextrap
        } // ncomp
    } // m_z_phys_nd
    Gpu::streamSynchronize();
}
