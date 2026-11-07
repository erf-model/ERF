#include "ERF.H"
#include "ERF_BndryPlaneIndex.H"

using namespace amrex;


/**
 * Impose boundary conditions using data read in as BndryRegisters from a previous ERF run
 *
 * The planes are read at level-0 resolution only, one cell wide just outside each lateral
 * face of the level-0 domain.  On a finer level they are interpolated linearly along the
 * face, in the two directions tangential to it, onto the fine ghost cells and faces; the
 * fine ghost cells take the same value at every distance from the face, as they do at
 * level 0.  The planes hold cell-centred values (the writer averages the face velocities to
 * cell centres), so at level 0 a cell-centred component copies the plane, and a velocity
 * component that is nodal along the face (y-velocity on an x face, x-velocity on a y face,
 * z-velocity on either) takes the average of the two plane cells either side of its node
 * (plane_index).
 *
 * @param[in]     lev       level the MultiFabs live on
 * @param[in,out] mfs       Vector of MultiFabs to be filled: cons, xvel, yvel, zvel
 * @param[in]     time      time at which the data should be filled
 * @param[in]     cons_only fill only mfs[Vars::cons]; the others are not touched
 */

void
ERF::fill_from_bndryregs (int lev, const Vector<MultiFab*>& mfs, const double time, bool cons_only)
{
    //
    // We now assume that if we read in on one face, we read in on all faces
    //
    AMREX_ALWAYS_ASSERT(m_r2d);

    // The plane data lives in level-0 index space
    const auto& pdom_lo = lbound(geom[0].Domain());
    const auto& pdom_hi = ubound(geom[0].Domain());

    // The cells to fill are those outside this level's domain
    const auto& dom_lo = lbound(geom[lev].Domain());
    const auto& dom_hi = ubound(geom[lev].Domain());

    IntVect rr(1);
    for (int l = 0; l < lev; ++l) {
        rr *= refRatio(l);
    }
    const int rx = rr[0];
    const int ry = rr[1];
    const int rz = rr[2];

    // Boundary-plane files are indexed by absolute simulation time.
    Vector<std::unique_ptr<PlaneVector>>& bndry_data = m_r2d->interp_in_time(time + start_time);

    const BCRec* bc_ptr = domain_bcs_type_d.data();

    // The ghost cells beside a face along a periodic direction take the plane at the
    // periodic image, as the interior cells across the seam do (bndry_plane_index)
    const bool xper = geom[lev].isPeriodic(0);
    const bool yper = geom[lev].isPeriodic(1);

    // xlo: ori = 0
    // ylo: ori = 1
    // zlo: ori = 2
    // xhi: ori = 3
    // yhi: ori = 4
    // zhi: ori = 5
    const auto& bdatxlo = (*bndry_data[0])[0].const_array();
    const auto& bdatylo = (*bndry_data[1])[0].const_array();
    const auto& bdatxhi = (*bndry_data[3])[0].const_array();
    const auto& bdatyhi = (*bndry_data[4])[0].const_array();

    int bccomp;

    for (int var_idx = 0; var_idx < Vars::NumTypes; ++var_idx)
    {
        if (cons_only && var_idx != Vars::cons) { continue; }

        MultiFab& mf = *mfs[var_idx];
        const int icomp = 0;
        const int ncomp = mf.nComp();
        // The ingestion check and the plane read use one slot: for a conserved
        // component its boundary-condition slot (every passive scalar shares one,
        // so an ingested scalar plane fills them all, and the planes are stored
        // by these slots, ReadBndryPlanes::plane_slot), for a velocity its own.
        const bool is_cons = (var_idx == Vars::cons);

        if (var_idx == Vars::xvel) {
           bccomp = BCVars::xvel_bc;
        } else if (var_idx == Vars::yvel) {
           bccomp = BCVars::yvel_bc;
        } else if (var_idx == Vars::zvel) {
           bccomp = BCVars::zvel_bc;
        } else  if (var_idx == Vars::cons) {
           bccomp = BCVars::cons_bc;
        }

        const bool nodal_x = (var_idx == Vars::xvel);
        const bool nodal_y = (var_idx == Vars::yvel);
        const bool nodal_z = (var_idx == Vars::zvel);

#ifdef AMREX_USE_OMP
#pragma omp parallel if (Gpu::notInLaunchRegion())
#endif
        for (MFIter mfi(mf); mfi.isValid(); ++mfi)
        {
            const Array4<Real>& dest_arr = mf.array(mfi);
            Box bx = mfi.growntilebox();

            // x-faces
            {
            Box bx_xlo(bx); bx_xlo.setBig(0,dom_lo.x-1);
            if (var_idx == Vars::xvel) bx_xlo.setBig(0,dom_lo.x);

            Box bx_xhi(bx); bx_xhi.setSmall(0,dom_hi.x+1);

            ParallelFor(
                bx_xlo, ncomp, [=] AMREX_GPU_DEVICE (int i, int j, int k, int n) {
                    const int bc_comp = is_cons ? cons_bc_comp(icomp+n) : bccomp+n;
                    if (bc_ptr[bc_comp].lo(0) == ERFBCType::ext_dir_ingested) {
                        const PlaneIndex pj = plane_index(j, ry, nodal_y, pdom_lo.y, pdom_hi.y, yper);
                        const PlaneIndex pk = plane_index(k, rz, nodal_z, pdom_lo.z, pdom_hi.z, false);
                        const int ib = pdom_lo.x-1;
                        const int c  = bc_comp;
                        dest_arr(i,j,k,icomp+n) =
                            (one-pk.w) * ((one-pj.w) * bdatxlo(ib,pj.lo,pk.lo,c) + pj.w * bdatxlo(ib,pj.hi,pk.lo,c))
                          +      pk.w  * ((one-pj.w) * bdatxlo(ib,pj.lo,pk.hi,c) + pj.w * bdatxlo(ib,pj.hi,pk.hi,c));
                    }
                },
                bx_xhi, ncomp, [=] AMREX_GPU_DEVICE (int i, int j, int k, int n) {
                    const int bc_comp = is_cons ? cons_bc_comp(icomp+n) : bccomp+n;
                    if (bc_ptr[bc_comp].hi(0) == ERFBCType::ext_dir_ingested) {
                        const PlaneIndex pj = plane_index(j, ry, nodal_y, pdom_lo.y, pdom_hi.y, yper);
                        const PlaneIndex pk = plane_index(k, rz, nodal_z, pdom_lo.z, pdom_hi.z, false);
                        const int ib = pdom_hi.x+1;
                        const int c  = bc_comp;
                        dest_arr(i,j,k,icomp+n) =
                            (one-pk.w) * ((one-pj.w) * bdatxhi(ib,pj.lo,pk.lo,c) + pj.w * bdatxhi(ib,pj.hi,pk.lo,c))
                          +      pk.w  * ((one-pj.w) * bdatxhi(ib,pj.lo,pk.hi,c) + pj.w * bdatxhi(ib,pj.hi,pk.hi,c));
                    }
                }
            );
            } // x-faces

            // y-faces
            {
            Box bx_ylo(bx); bx_ylo.setBig  (1,dom_lo.y-1);
            if (var_idx == Vars::yvel) bx_ylo.setBig(1,dom_lo.y);

            Box bx_yhi(bx); bx_yhi.setSmall(1,dom_hi.y+1);

            ParallelFor(
               bx_ylo, ncomp, [=] AMREX_GPU_DEVICE (int i, int j, int k, int n) {
                    const int bc_comp = is_cons ? cons_bc_comp(icomp+n) : bccomp+n;
                    if (bc_ptr[bc_comp].lo(1) == ERFBCType::ext_dir_ingested) {
                        const PlaneIndex pi = plane_index(i, rx, nodal_x, pdom_lo.x, pdom_hi.x, xper);
                        const PlaneIndex pk = plane_index(k, rz, nodal_z, pdom_lo.z, pdom_hi.z, false);
                        const int jb = pdom_lo.y-1;
                        const int c  = bc_comp;
                        dest_arr(i,j,k,icomp+n) =
                            (one-pk.w) * ((one-pi.w) * bdatylo(pi.lo,jb,pk.lo,c) + pi.w * bdatylo(pi.hi,jb,pk.lo,c))
                          +      pk.w  * ((one-pi.w) * bdatylo(pi.lo,jb,pk.hi,c) + pi.w * bdatylo(pi.hi,jb,pk.hi,c));
                    }
                },
                bx_yhi, ncomp, [=] AMREX_GPU_DEVICE (int i, int j, int k, int n) {
                    const int bc_comp = is_cons ? cons_bc_comp(icomp+n) : bccomp+n;
                    if (bc_ptr[bc_comp].hi(1) == ERFBCType::ext_dir_ingested) {
                        const PlaneIndex pi = plane_index(i, rx, nodal_x, pdom_lo.x, pdom_hi.x, xper);
                        const PlaneIndex pk = plane_index(k, rz, nodal_z, pdom_lo.z, pdom_hi.z, false);
                        const int jb = pdom_hi.y+1;
                        const int c  = bc_comp;
                        dest_arr(i,j,k,icomp+n) =
                            (one-pk.w) * ((one-pi.w) * bdatyhi(pi.lo,jb,pk.lo,c) + pi.w * bdatyhi(pi.hi,jb,pk.lo,c))
                          +      pk.w  * ((one-pi.w) * bdatyhi(pi.lo,jb,pk.hi,c) + pi.w * bdatyhi(pi.hi,jb,pk.hi,c));
                    }
                }
            );
            } // y-faces
        } // mf
    } // var_idx
}
