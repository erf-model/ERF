#include "AMReX_PhysBCFunct.H"
#include "ERF_PhysBCFunct.H"
#include "ERF_Constants.H"
#include "ERF_TerrainMetrics.H"

using namespace amrex;

/*
 * Impose lateral boundary conditions on the base state
 *
 * @param[in,out] dest_arr cell-centered data to be filled
 * @param[in]     bx       box holding data to be filled
 * @param[in]     domain   simulation domain
 */

void ERFPhysBCFunct_base::impose_lateral_basestate_bcs (const Array4<Real>& dest_arr,
                                                        const Box& bx, const Box& domain,
                                                        int ncomp, const IntVect& nghost)
{
    BL_PROFILE_VAR("impose_lateral_base_bcs()",impose_lateral_base_bcs);
    //
    // Note that the "bx" that comes in here has already been grown in the lateral directions
    //     but not in the vertical
    //

    const int* bxlo = bx.loVect();
    const int* bxhi = bx.hiVect();

    const int* dlo  = domain.loVect();
    const int* dhi  = domain.hiVect();

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

    int bc_comp = BaseBCVars::rho0_bc_comp;

    for (int nc = 0; nc < ncomp; nc++)
    {
        for (int dir = 0; dir < AMREX_SPACEDIM; dir++)
        {
            bcrs[nc].setLo(dir, ( bxlo[dir]<=dlo[dir]
                                 ? m_domain_bcs_type[bc_comp].lo(dir) : BCType::int_dir ));
            bcrs[nc].setHi(dir, ( bxhi[dir]>=dhi[dir]
                                 ? m_domain_bcs_type[bc_comp].hi(dir) : BCType::int_dir ));
        }
    }

    Gpu::DeviceVector<BCRec> bcrs_d(ncomp);
    Gpu::copyAsync(Gpu::hostToDevice, bcrs.begin(), bcrs.end(), bcrs_d.begin());
    const BCRec* bc_ptr = bcrs_d.data();

    GeometryData const& geomdata = m_geom.data();
    bool is_periodic_in_x = geomdata.isPeriodic(0);
    bool is_periodic_in_y = geomdata.isPeriodic(1);

    // Do ghost cells in x-direction but not reaching out in y
    // The corners we miss here will be covered in the y-loop below or by periodicity
    if (!is_periodic_in_x)
    {
        // Populate ghost cells on lo-x and hi-x domain boundaries
        Box bx_xlo(bx);  bx_xlo.setBig  (0,dom_lo.x-nghost[0]);
        Box bx_xhi(bx);  bx_xhi.setSmall(0,dom_hi.x+nghost[0]);

        ParallelFor(
            bx_xlo, ncomp, [=] AMREX_GPU_DEVICE (int i, int j, int k, int n)
            {
                int dest_comp = n;
                int l_bc_type = bc_ptr[n].lo(0);
                int iflip = dom_lo.x - 1 - i;
                if (l_bc_type == ERFBCType::foextrap) {
                    dest_arr(i,j,k,dest_comp) =  dest_arr(dom_lo.x,j,k,dest_comp);
                } else if (l_bc_type == ERFBCType::open) {
                    dest_arr(i,j,k,dest_comp) =  dest_arr(dom_lo.x,j,k,dest_comp);
                } else if (l_bc_type == ERFBCType::reflect_even) {
                    dest_arr(i,j,k,dest_comp) =  dest_arr(iflip,j,k,dest_comp);
                }
            },
            bx_xhi, ncomp, [=] AMREX_GPU_DEVICE (int i, int j, int k, int n)
            {
                int dest_comp = n;
                int h_bc_type = bc_ptr[n].hi(0);
                int iflip =  2*dom_hi.x + 1 - i;
                if (h_bc_type == ERFBCType::foextrap) {
                    dest_arr(i,j,k,dest_comp) =  dest_arr(dom_hi.x,j,k,dest_comp);
                } else if (h_bc_type == ERFBCType::open) {
                    dest_arr(i,j,k,dest_comp) =  dest_arr(dom_hi.x,j,k,dest_comp);
                } else if (h_bc_type == ERFBCType::reflect_even) {
                    dest_arr(i,j,k,dest_comp) =  dest_arr(iflip,j,k,dest_comp);
                }
            }
        );
    }

    if (!is_periodic_in_y)
    {
        // Populate ghost cells on lo-y and hi-y domain boundaries
        Box bx_ylo(bx);  bx_ylo.setBig  (1,dom_lo.y-nghost[1]);
        Box bx_yhi(bx);  bx_yhi.setSmall(1,dom_hi.y+nghost[1]);
        if (bx_ylo.smallEnd(2) != domain.smallEnd(2)) bx_ylo.growLo(2,nghost[2]);
        if (bx_ylo.bigEnd(2)   != domain.bigEnd(2))   bx_ylo.growHi(2,nghost[2]);
        if (bx_yhi.smallEnd(2) != domain.smallEnd(2)) bx_yhi.growLo(2,nghost[2]);
        if (bx_yhi.bigEnd(2)   != domain.bigEnd(2))   bx_yhi.growHi(2,nghost[2]);
        ParallelFor(
            bx_ylo, ncomp, [=] AMREX_GPU_DEVICE (int i, int j, int k, int n)
            {
                int dest_comp = n;
                int l_bc_type = bc_ptr[n].lo(1);
                int jflip = dom_lo.y - 1 - j;
                if (l_bc_type == ERFBCType::foextrap) {
                    dest_arr(i,j,k,dest_comp) =  dest_arr(i,dom_lo.y,k,dest_comp);
                } else if (l_bc_type == ERFBCType::open) {
                    dest_arr(i,j,k,dest_comp) =  dest_arr(i,dom_lo.y,k,dest_comp);
                } else if (l_bc_type == ERFBCType::reflect_even) {
                    dest_arr(i,j,k,dest_comp) =  dest_arr(i,jflip,k,dest_comp);
                }

            },
            bx_yhi, ncomp, [=] AMREX_GPU_DEVICE (int i, int j, int k, int n)
            {
                int dest_comp = n;
                int h_bc_type = bc_ptr[n].hi(1);
                int jflip =  2*dom_hi.y + 1 - j;
                if (h_bc_type == ERFBCType::foextrap) {
                    dest_arr(i,j,k,dest_comp) =  dest_arr(i,dom_hi.y,k,dest_comp);
                } else if (h_bc_type == ERFBCType::open) {
                    dest_arr(i,j,k,dest_comp) =  dest_arr(i,dom_hi.y,k,dest_comp);
                } else if (h_bc_type == ERFBCType::reflect_even) {
                    dest_arr(i,j,k,dest_comp) =  dest_arr(i,jflip,k,dest_comp);
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
        if (bx_xlo.smallEnd(2) != domain.smallEnd(2)) bx_xlo.growLo(2,nghost[2]);
        if (bx_xlo.bigEnd(2)   != domain.bigEnd(2))   bx_xlo.growHi(2,nghost[2]);
        if (bx_xhi.smallEnd(2) != domain.smallEnd(2)) bx_xhi.growLo(2,nghost[2]);
        if (bx_xhi.bigEnd(2)   != domain.bigEnd(2))   bx_xhi.growHi(2,nghost[2]);
        ParallelFor(
            bx_xlo, ncomp, [=] AMREX_GPU_DEVICE (int i, int j, int k, int n)
            {
                int dest_comp = n;
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
                } else if (l_bc_type == ERFBCType::hoextrap) {
                    Real delta_i = (dom_lo.x - i);
                    dest_arr(i,j,k,dest_comp) = (1.0 + delta_i)*dest_arr(dom_lo.x,j,k,dest_comp) - delta_i*dest_arr(dom_lo.x+1,j,k,dest_comp) ;
                }
            },
            bx_xhi, ncomp, [=] AMREX_GPU_DEVICE (int i, int j, int k, int n)
            {
                int dest_comp = n;
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
                } else if (h_bc_type == ERFBCType::hoextrap) {
                    Real delta_i = (i - dom_hi.x);
                    dest_arr(i,j,k,dest_comp) = (1.0 + delta_i)*dest_arr(dom_hi.x,j,k,dest_comp) - delta_i*dest_arr(dom_hi.x-1,j,k,dest_comp) ;
                }
            }
        );
    }

    if (!is_periodic_in_y)
    {
        // Populate ghost cells on lo-y and hi-y domain boundaries
        Box bx_ylo(bx);  bx_ylo.setBig  (1,dom_lo.y-1);
        Box bx_yhi(bx);  bx_yhi.setSmall(1,dom_hi.y+1);
        if (bx_ylo.smallEnd(2) != domain.smallEnd(2)) bx_ylo.growLo(2,nghost[2]);
        if (bx_ylo.bigEnd(2)   != domain.bigEnd(2))   bx_ylo.growHi(2,nghost[2]);
        if (bx_yhi.smallEnd(2) != domain.smallEnd(2)) bx_yhi.growLo(2,nghost[2]);
        if (bx_yhi.bigEnd(2)   != domain.bigEnd(2))   bx_yhi.growHi(2,nghost[2]);
        ParallelFor(
            bx_ylo, ncomp, [=] AMREX_GPU_DEVICE (int i, int j, int k, int n)
            {
                int dest_comp = n;
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
                } else if (l_bc_type == ERFBCType::hoextrap) {
                    Real delta_j = (dom_lo.y - j);
                    dest_arr(i,j,k,dest_comp) = (1.0 + delta_j)*dest_arr(i,dom_lo.y,k,dest_comp) - delta_j*dest_arr(i,dom_lo.y+1,k,dest_comp) ;
                }

            },
            bx_yhi, ncomp, [=] AMREX_GPU_DEVICE (int i, int j, int k, int n)
            {
                int dest_comp = n;
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
                } else if (h_bc_type == ERFBCType::hoextrap) {
                    Real delta_j = (j - dom_hi.y);
                    dest_arr(i,j,k,dest_comp) = (1.0 + delta_j)*dest_arr(i,dom_hi.y,k,dest_comp) - delta_j*dest_arr(i,dom_hi.y-1,k,dest_comp);
                }
            }
        );
    }
    Gpu::streamSynchronize();
}

void ERFPhysBCFunct_base::impose_vertical_basestate_bcs (const Array4<Real>& dest_arr,
                                                         const Array4<Real const>& /*z_phys_nd*/,
                                                         const Box& bx,
                                                         const Box& domain,
                                                         int ncomp,
                                                         const IntVect& /*nghost*/)
{
    BL_PROFILE_VAR("impose_vertical_base_bcs()",impose_vertical_base_bcs);

    const auto& dom_lo = lbound(domain);
    const auto& dom_hi = ubound(domain);

    Box bx_zlo1(bx); bx_zlo1.setBig(2,dom_lo.z-1); if (bx_zlo1.ok()) bx_zlo1.setSmall(2,dom_lo.z-1);
    ParallelFor(
        bx_zlo1, [=] AMREX_GPU_DEVICE (int i, int j, int k)
        {
            dest_arr(i,j,k,BaseState::r0_comp)  = dest_arr(i,j,dom_lo.z,BaseState::r0_comp);
            dest_arr(i,j,k,BaseState::p0_comp)  = dest_arr(i,j,dom_lo.z,BaseState::p0_comp);
            dest_arr(i,j,k,BaseState::pi0_comp) = dest_arr(i,j,dom_lo.z,BaseState::pi0_comp);
            dest_arr(i,j,k,BaseState::th0_comp) = dest_arr(i,j,dom_lo.z,BaseState::th0_comp);
        }
    );

    Box bx_zlo(bx); bx_zlo.setBig(2,dom_lo.z-2);
    Box bx_zhi(bx); bx_zhi.setSmall(2,dom_hi.z+1);
    ParallelFor(
        bx_zlo, ncomp, [=] AMREX_GPU_DEVICE (int i, int j, int k, int n)
        {
            dest_arr(i,j,k,n) = dest_arr(i,j,dom_lo.z-1,n);
        },
        bx_zhi, ncomp, [=] AMREX_GPU_DEVICE (int i, int j, int k, int n)
        {
            dest_arr(i,j,k,n) = dest_arr(i,j,dom_hi.z,n);
        }
    );
}
