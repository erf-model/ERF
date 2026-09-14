#include "AMReX_PhysBCFunct.H"
#include "ERF_PhysBCFunct.H"
#include "ERF_NumericalConstants.H"
#include "ERF_TerrainMetrics.H"
#include "ERF_EOS.H"

using namespace amrex;

/**
 * Impose lateral boundary conditions on the base state
 *
 * @param[in,out] dest_arr cell-centered data to be filled
 * @param[in]     z_nd     nodal physical height
 * @param[in]     bx       box holding data to be filled
 * @param[in]     domain   simulation domain
 * @param[in]     ncomp    number of base-state components to fill
 * @param[in]     nghost   number of ghost cells in each coordinate direction
 */

void ERFPhysBCFunct_base::impose_lateral_basestate_bcs (const Array4<Real>& dest_arr,
                                                        const Array4<Real const>& z_nd,
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
        Box bx_xlo(bx);  bx_xlo.setBig  (0,dom_lo.x-1);
        Box bx_xhi(bx);  bx_xhi.setSmall(0,dom_hi.x+1);

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
                    Real delta_i = static_cast<Real>(dom_lo.x - i);
                    dest_arr(i,j,k,dest_comp) = (one + delta_i)*dest_arr(dom_lo.x,j,k,dest_comp) - delta_i*dest_arr(dom_lo.x+1,j,k,dest_comp) ;
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
                    Real delta_i = static_cast<Real>(i - dom_hi.x);
                    dest_arr(i,j,k,dest_comp) = (one + delta_i)*dest_arr(dom_hi.x,j,k,dest_comp) - delta_i*dest_arr(dom_hi.x-1,j,k,dest_comp) ;
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
                    Real delta_j = static_cast<Real>(dom_lo.y - j);
                    dest_arr(i,j,k,dest_comp) = (one + delta_j)*dest_arr(i,dom_lo.y,k,dest_comp) - delta_j*dest_arr(i,dom_lo.y+1,k,dest_comp) ;
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
                    Real delta_j = static_cast<Real>(j - dom_hi.y);
                    dest_arr(i,j,k,dest_comp) = (one + delta_j)*dest_arr(i,dom_hi.y,k,dest_comp) - delta_j*dest_arr(i,dom_hi.y-1,k,dest_comp);
                }
            }
        );
    }

    //
    // The copies above give a lateral ghost cell the base state of the cell it copies from,
    // which on a terrain-fitted mesh sits at a different height: the nodal mesh is extrapolated
    // past the domain, so z_cc of a ghost cell differs from z_cc of that cell whenever the
    // terrain has a slope at the boundary.  The base state left there is then neither the
    // hydrostatic profile at the ghost cell's own height nor in hydrostatic balance along the
    // ghost column, while detJ and the metric terms in that same column are height-consistent.
    //
    // Here we rebuild the base state of those cells at the height the mesh says the cell is at:
    // p_0 is carried across the height offset by the hydrostatic relation, and rho_0 and pi_0
    // follow from the equation of state so that it holds exactly in the ghost cell.
    //
    // This is done only where the boundary is foextrap or open: the reflecting and hoextrap
    // conditions above are not copies of a single cell and are left alone.  It reduces to the
    // identity, bit for bit, wherever the two heights agree -- every constant-dz mesh, and any
    // terrain that is flat at the boundary -- so it changes nothing there.
    //
    if (m_use_terrain && z_nd && ncomp == BaseState::num_comps)
    {
        auto is_ext_bc = [] (int bc_type) {
            return (bc_type == ERFBCType::foextrap || bc_type == ERFBCType::open);
        };

        const bool xlo_ext = !is_periodic_in_x && is_ext_bc(bcrs[0].lo(0));
        const bool xhi_ext = !is_periodic_in_x && is_ext_bc(bcrs[0].hi(0));
        const bool ylo_ext = !is_periodic_in_y && is_ext_bc(bcrs[0].lo(1));
        const bool yhi_ext = !is_periodic_in_y && is_ext_bc(bcrs[0].hi(1));

        if (xlo_ext || xhi_ext || ylo_ext || yhi_ext)
        {
            //
            // "bx" is grown laterally but not vertically, so we grow it in z the same way the
            // loops above do, then trim it to the domain: the ghost cells above and below the
            // domain are filled after us, by impose_vertical_basestate_bcs.
            //
            Box gbx(bx);
            if (gbx.smallEnd(2) != domain.smallEnd(2)) gbx.growLo(2,nghost[2]);
            if (gbx.bigEnd(2)   != domain.bigEnd(2))   gbx.growHi(2,nghost[2]);
            gbx.setSmall(2, amrex::max(gbx.smallEnd(2), dom_lo.z));
            gbx.setBig  (2, amrex::min(gbx.bigEnd(2)  , dom_hi.z));

            const Real l_rdOcp   = m_rdOcp;
            const Real l_gravity = m_gravity;

            ParallelFor(gbx, [=] AMREX_GPU_DEVICE (int i, int j, int k)
            {
                //
                // The cell this one was copied from: the nearest cell that is inside the domain
                // laterally, which is just what the x-pass followed by the y-pass above leaves
                // here, corners included.
                //
                const bool out_xlo = (i < dom_lo.x);
                const bool out_xhi = (i > dom_hi.x);
                const bool out_ylo = (j < dom_lo.y);
                const bool out_yhi = (j > dom_hi.y);

                if (!(out_xlo || out_xhi || out_ylo || out_yhi)) { return; }

                // Leave the cell alone unless every boundary it lies outside of is foextrap/open
                if (out_xlo && !xlo_ext) { return; }
                if (out_xhi && !xhi_ext) { return; }
                if (out_ylo && !ylo_ext) { return; }
                if (out_yhi && !yhi_ext) { return; }

                const int ir = amrex::min(amrex::max(i,dom_lo.x),dom_hi.x);
                const int jr = amrex::min(amrex::max(j,dom_lo.y),dom_hi.y);

                auto zcc = [=] (int ii, int jj, int kk) {
                    return Real(0.125) * ( z_nd(ii,jj  ,kk  ) + z_nd(ii+1,jj  ,kk  )
                                         + z_nd(ii,jj+1,kk  ) + z_nd(ii+1,jj+1,kk  )
                                         + z_nd(ii,jj  ,kk+1) + z_nd(ii+1,jj  ,kk+1)
                                         + z_nd(ii,jj+1,kk+1) + z_nd(ii+1,jj+1,kk+1) );
                };

                const Real dz_offset = zcc(i,j,k) - zcc(ir,jr,k);

                // Same height as the cell we copied from: the copy is already right
                if (dz_offset == Real(0.0)) { return; }

                const Real  r0_ref = dest_arr(ir,jr,k,BaseState::r0_comp );
                const Real  p0_ref = dest_arr(ir,jr,k,BaseState::p0_comp );
                const Real th0_ref = dest_arr(ir,jr,k,BaseState::th0_comp);
                const Real qv0_ref = dest_arr(ir,jr,k,BaseState::qv0_comp);

                //
                // A base state that has not been built yet -- this is called on a MultiFab that
                // was only just allocated and zeroed -- has nothing to transfer, and the scale
                // height below would divide by zero.  Leave the copy alone.
                //
                if (!(p0_ref > Real(0.0)) || !(r0_ref > Real(0.0)) || !(th0_ref > Real(0.0))) { return; }

                //
                // p_0 across the height offset, using the scale height of the cell we copied
                // from: p_0 * exp(-g dz / (R_d T_v)) with R_d T_v = p_0 / rho_0.
                //
                const Real p0_new = p0_ref * std::exp(-l_gravity * dz_offset * r0_ref / p0_ref);

                //
                // theta_0 and qv_0 carry over unchanged.  Reconstructing them at the height of
                // this cell would need the profile of the column above and below k, and the
                // cells we could read there are the ones this box happens to hold, so the
                // stencil -- and with it the answer -- would depend on where the grids are
                // split in z.  Holding them fixed leaves an error of dtheta_0/dz * dz, which
                // for the offsets a terrain-fitted mesh produces is a few thousandths of a
                // kelvin, against the several pascals in p_0 that the transfer above removes.
                //
                dest_arr(i,j,k,BaseState::p0_comp ) = p0_new;
                dest_arr(i,j,k,BaseState::th0_comp) = th0_ref;
                dest_arr(i,j,k,BaseState::qv0_comp) = qv0_ref;
                dest_arr(i,j,k,BaseState::pi0_comp) = getExnergivenP(p0_new, l_rdOcp);
                dest_arr(i,j,k,BaseState::r0_comp ) = getRhogivenThetaPress(th0_ref, p0_new,
                                                                           l_rdOcp, qv0_ref);
            });
        }
    }

    Gpu::streamSynchronize();
}

/**
 * Impose vertical boundary conditions on the base state
 *
 * @param[in,out] dest_arr  cell-centered base-state data to be filled
 * @param[in]     z_phys_nd height coordinate at nodes, unused for base-state fills
 * @param[in]     bx        box holding data to be filled
 * @param[in]     domain    simulation domain
 * @param[in]     ncomp     number of base-state components to fill
 * @param[in]     nghost    number of ghost cells, unused for base-state fills
 */
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

    Box bx_zlo(bx); bx_zlo.setBig(2,dom_lo.z-1);
    Box bx_zhi(bx); bx_zhi.setSmall(2,dom_hi.z+1);
    ParallelFor(
        bx_zlo, ncomp, [=] AMREX_GPU_DEVICE (int i, int j, int k, int n)
        {
            dest_arr(i,j,k,n) = dest_arr(i,j,dom_lo.z,n);
        },
        bx_zhi, ncomp, [=] AMREX_GPU_DEVICE (int i, int j, int k, int n)
        {
            dest_arr(i,j,k,n) = dest_arr(i,j,dom_hi.z,n);
        }
    );
}
