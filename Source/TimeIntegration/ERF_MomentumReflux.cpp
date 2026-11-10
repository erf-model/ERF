/**
 * \file ERF_MomentumReflux.cpp
 *
 * The reflux of the pressure impulse on the faces of a coarse level that lie on the coarse/fine
 * boundary and are normal to the momentum on them (erf.momentum_reflux, off by default;
 * erf-model/ERF#5017).
 *
 * The momenta live on faces, and the momentum flux of a face's control volume is at the cell
 * centres on either side of it.  Under two-way coupling the composite solution is the fine level
 * where there is one and the coarse level elsewhere, and the control volume of a coarse/fine face
 * then runs from the centre of the uncovered coarse cell beside it to the centre of the first fine
 * cell, a width h = (dx_c + dx_f) / 2.  The fine level does not evolve its coarse/fine faces (they
 * are set from the coarse level), and its first interior face takes its flux at the centre of the
 * first fine cell.  The coarse level, though, updates the coarse/fine face with its own flux at the
 * centre of the covered coarse cell, over its own width dx_c.  What the fine interior gives up at
 * its edge is not what the coarse/fine face receives, and nothing puts the difference back: the
 * composite momentum is not conserved (#5017).  In the cases there the difference is almost all
 * pressure.
 *
 * So each level adds up, over each coarse step of level 0, the pressure its x- and y-momenta saw in
 * each cell: the stage pressure of the slow right-hand side and the acoustic part of the substeps,
 * in the RK stage whose update the step ends with, each times the interval it is applied over.
 * At the end of the step (post_timestep), for each such coarse/fine face of the coarse level, the
 * impulse that the pressure difference across the composite control volume, from the uncovered
 * coarse cell to the first fine cell, over h, would have given it, less the one the coarse level
 * gave it, is
 *
 *     delta (rho u) = - (P_f - P_c(I-1)) / h + (P_c(I) - P_c(I-1)) / dx_c     (fine level on the high side)
 *
 * with P the time integral of the pressure and P_f that of the first fine cells, averaged over the
 * face.  The composite x-momentum (the control volume of each coarse/fine face being h wide) then
 * changes, up to the approximations below, only by what crosses the domain boundary.
 *
 * The impulse is not added to the momenta at once but as a momentum source on those faces of the
 * coarse level, at a constant rate over its steps in the next coarse step of level 0; the fine
 * level's coarse/fine faces follow the coarse level's as always.  Added at once, it acts as an
 * explicit acoustic update over a whole coarse step, at an acoustic Courant number far above one,
 * and where the coarse/fine faces of two directions meet (the vertical edges of a patch) it grows
 * from one step to the next.  As a source, the acoustic substeps take it in.  So the momentum at
 * the end of a step is short by the impulse of that step, which the next step gives; and that last
 * impulse is not in a checkpoint, so a restarted run does not have it.  It is dropped when the
 * coarse level's grids change, and it is not made in a step in which those of either level change.
 *
 * Only the pressure is refluxed: not the advective or diffusive momentum fluxes, nor the momentum
 * sources on these faces, nor the momentum on coarse/fine faces parallel to it, nor (rho w) on the
 * z-faces.  The pressure is taken as a cell value: the acoustic part is gamma R_d pi (rho theta)''
 * with the Exner function pi of the cell, where the substep uses its average over the face, and the
 * factor 1/(1+q_t) is that of the cell where the substep uses its face average.  Both differences
 * are of second order in the perturbation.  Map factors are taken to be one.  This is done for
 * fully compressible levels without terrain or EB.
 */
#include "ERF.H"
#include "ERF_Constants.H"
#include "ERF_EOS.H"

using namespace amrex;

/**
 * Whether the pressure impulse on the coarse/fine faces normal to the momentum is refluxed
 * (erf.momentum_reflux, false by default)
 */
bool
ERF::momentum_reflux_enabled ()
{
    static const bool enabled = [] {
        bool on = false;
        ParmParse("erf").query("momentum_reflux", on);
        return on;
    }();
    return enabled;
}

/**
 * The time integral of the pressure that the momenta of level lev have seen in each cell since the
 * start of the coarse step of level 0, to which the caller adds; nullptr unless the momentum reflux
 * is on.  Made again (zero) when the level's grids have changed, in which case its sum is
 * incomplete and the level is not refluxed in this step.
 *
 * @param lev the level
 */
MultiFab*
ERF::momentum_reflux_pressure (int lev)
{
    if (!momentum_reflux_enabled() || solverChoice.coupling_type != CouplingType::TwoWay ||
        finest_level == 0) {
        return nullptr;
    }

    static bool checked = false;
    if (!checked) {
        for (int l = 0; l <= max_level; ++l) {
            if (solverChoice.compressibility_type[l] != CompressibilityType::FullyCompressible) {
                Abort("erf.momentum_reflux is only implemented for fully compressible levels");
            }
        }
        if (SolverChoice::mesh_type == MeshType::VariableDz ||
            SolverChoice::terrain_type == TerrainType::EB ||
            SolverChoice::terrain_type == TerrainType::StaticFittedMesh ||
            SolverChoice::terrain_type == TerrainType::MovingFittedMesh) {
            Abort("erf.momentum_reflux is only implemented without terrain or EB");
        }
        checked = true;
    }

    if (static_cast<int>(m_mom_reflux_p.size()) <= lev) {
        m_mom_reflux_p.resize(lev+1);
        m_mom_reflux_used.resize(lev+1, 0);
        m_mom_reflux_stale.resize(lev+1, 0);
    }
    MultiFab& P = m_mom_reflux_p[lev];
    if (!P.isDefined() || P.boxArray() != grids[lev] || P.DistributionMap() != dmap[lev]) {
        if (m_mom_reflux_used[lev]) { m_mom_reflux_stale[lev] = 1; }
        P.define(grids[lev], dmap[lev], 1, 1);
        P.setVal(0);
    }
    m_mom_reflux_used[lev] = 1;
    return &P;
}

/**
 * Add the stage pressure of the slow right-hand side (the perturbation pressure, from the state
 * the stage's pressure gradient is made from) times the interval it is applied over
 *
 * @param lev   the level
 * @param P     the sum to add to (momentum_reflux_pressure)
 * @param cons  the state of the stage
 * @param p0    the base-state pressure
 * @param qt    the total moisture
 * @param moist whether qt is defined
 * @param dt_stage the interval
 */
void
ERF::momentum_reflux_add_slow (int /*lev*/, MultiFab& P, const MultiFab& cons, const MultiFab& p0,
                               const MultiFab& qt, bool moist, Real dt_stage)
{
#ifdef _OPENMP
#pragma omp parallel if (Gpu::notInLaunchRegion())
#endif
    for (MFIter mfi(P, TilingIfNotGPU()); mfi.isValid(); ++mfi) {
        const Box& bx = mfi.tilebox();
        const Array4<Real>&       p_arr  = P.array(mfi);
        const Array4<const Real>& c_arr  = cons.const_array(mfi);
        const Array4<const Real>& p0_arr = p0.const_array(mfi);
        const Array4<const Real>& qt_arr = (moist) ? qt.const_array(mfi) : Array4<const Real>{};
        ParallelFor(bx, [=] AMREX_GPU_DEVICE (int i, int j, int k) noexcept
        {
            Real qv = (moist) ? c_arr(i,j,k,RhoQ1_comp) / c_arr(i,j,k,Rho_comp) : zero;
            Real q  = (moist) ? qt_arr(i,j,k) : zero;
            Real pp = getPgivenRTh(c_arr(i,j,k,RhoTheta_comp), qv) - p0_arr(i,j,k);
            p_arr(i,j,k) += dt_stage * pp / (one + q);
        });
    }
}

/**
 * The impulse that the reflux of the pressure gives the coarse/fine faces of level crse_lev normal
 * to the x- and y-momenta, from the pressure sums of the step of level 0 just taken; stored for
 * momentum_reflux_add_sources, which gives it over the next one
 *
 * @param crse_lev the coarse level
 */
void
ERF::momentum_reflux (int crse_lev)
{
    BL_PROFILE("ERF::momentum_reflux()");

    const int c = crse_lev;
    const int f = crse_lev + 1;
    if (!momentum_reflux_enabled() || solverChoice.coupling_type != CouplingType::TwoWay ||
        f > finest_level || static_cast<int>(m_mom_reflux_p.size()) <= f) {
        return;
    }
    for (int l : {c, f}) {
        const MultiFab& P = m_mom_reflux_p[l];
        if (!P.isDefined() || P.boxArray() != grids[l] || P.DistributionMap() != dmap[l] ||
            m_mom_reflux_stale[l]) {
            if (verbose > 0) {
                Print() << "Momentum reflux skipped at level " << c << " (grids changed in the step)\n";
            }
            if (static_cast<int>(m_mom_reflux_src.size()) > c) {
                for (auto& D : m_mom_reflux_src[c]) { D.clear(); }
            }
            return;
        }
    }

    if (static_cast<int>(m_mom_reflux_src.size()) <= c) { m_mom_reflux_src.resize(c+1); }

    const IntVect rr = refRatio(c);
    MultiFab& Pc = m_mom_reflux_p[c];
    const MultiFab& Pf = m_mom_reflux_p[f];
    Pc.FillBoundary(geom[c].periodicity());

    MultiFab& rho_c = vars_new[c][Vars::cons];
    rho_c.FillBoundary(geom[c].periodicity());

    // Which coarse cells, ghost cells included, the fine level covers
    const iMultiFab fmask = makeFineMask(grids[c], dmap[c], IntVect(1), grids[f], rr,
                                         geom[c].periodicity(), 0, 1);

    const BoxArray cfba = coarsen(grids[f], rr);
    const Box& domain = geom[c].Domain();

    Real max_du = zero, sum_dmom = zero;

    for (int dir = 0; dir < 2; ++dir)
    {
        const int r  = rr[dir];
        const int t1 = (dir == 0) ? 1 : 0;
        const int t2 = 2;
        const int n1 = rr[t1];
        const int n2 = rr[t2];
        const Real inv_n = one / static_cast<Real>(n1 * n2);

        // The pressure sum of the first and last layer of fine cells in each covered coarse cell
        // (along dir), averaged over the fine cells across it
        MultiFab Lf(cfba, dmap[f], 2, 0);
#ifdef _OPENMP
#pragma omp parallel if (Gpu::notInLaunchRegion())
#endif
        for (MFIter mfi(Lf, TilingIfNotGPU()); mfi.isValid(); ++mfi) {
            const Box& bx = mfi.tilebox();
            const Array4<Real>&       l_arr = Lf.array(mfi);
            const Array4<const Real>& p_arr = Pf.const_array(mfi);
            ParallelFor(bx, [=] AMREX_GPU_DEVICE (int i, int j, int k) noexcept
            {
                IntVect lo(i*rr[0], j*rr[1], k*rr[2]);
                Real slo = zero, shi = zero;
                for (int b = 0; b < n1; ++b) {
                    for (int a = 0; a < n2; ++a) {
                        IntVect iv = lo;
                        iv[t1] += b;
                        iv[t2] += a;
                        slo += p_arr(iv);
                        iv[dir] += r - 1;
                        shi += p_arr(iv);
                    }
                }
                l_arr(i,j,k,0) = slo * inv_n;
                l_arr(i,j,k,1) = shi * inv_n;
            });
        }
        MultiFab Lc(grids[c], dmap[c], 2, 1);
        Lc.setVal(0);
        Lc.ParallelCopy(Lf, 0, 0, 2, IntVect(0), IntVect(1), geom[c].periodicity());

        const Real dxc   = geom[c].CellSize(dir);
        const Real dxf   = geom[f].CellSize(dir);
        const Real dxc_i = one / dxc;
        const Real h_i   = one / (myhalf * (dxc + dxf));
        const Real dV    = geom[c].CellSize(t1) * geom[c].CellSize(t2) / h_i;
        const bool per   = geom[c].isPeriodic(dir);
        const int  dlo   = domain.smallEnd(dir);
        const int  dhi   = domain.bigEnd(dir) + 1;

        // The impulse, (rho u) per unit volume, on the coarse/fine faces of level c
        MultiFab& D = m_mom_reflux_src[c][dir];
        D.define(convert(grids[c], IntVect::TheDimensionVector(dir)), dmap[c], 1, 0);
        ReduceOps<ReduceOpMax, ReduceOpSum> reduce_op;
        ReduceData<Real, Real> reduce_data(reduce_op);
        using ReduceTuple = typename decltype(reduce_data)::Type;
        for (MFIter mfi(D); mfi.isValid(); ++mfi) {
            const Box& bx = mfi.validbox();
            const Array4<Real>&       d_arr = D.array(mfi);
            const Array4<const Real>& r_arr = rho_c.const_array(mfi);
            const Array4<const Real>& p_arr = Pc.const_array(mfi);
            const Array4<const Real>& l_arr = Lc.const_array(mfi);
            const Array4<const int>&  m_msk = fmask.const_array(mfi);
            // A face shared by two grids takes the same impulse in both, and is counted in the
            // one it is the low face of
            const int bhi = bx.bigEnd(dir);
            reduce_op.eval(bx, reduce_data, [=] AMREX_GPU_DEVICE (int i, int j, int k) -> ReduceTuple
            {
                IntVect iv(i,j,k);
                IntVect im = iv; im[dir] -= 1;
                d_arr(iv) = zero;
                if (!per && (iv[dir] == dlo || iv[dir] == dhi)) { return {zero, zero}; }
                const int fhi = m_msk(iv);
                const int flo = m_msk(im);
                if (fhi == flo) { return {zero, zero}; }
                Real dmom_crse = -(p_arr(iv) - p_arr(im)) * dxc_i;
                Real dmom_comp = (fhi) ? -(l_arr(iv,0) - p_arr(im)) * h_i
                                       : -(p_arr(iv) - l_arr(im,1)) * h_i;
                Real dmom = dmom_comp - dmom_crse;
                d_arr(iv) = dmom;
                Real du = dmom / (myhalf * (r_arr(im,Rho_comp) + r_arr(iv,Rho_comp)));
                return {std::abs(du), (iv[dir] < bhi) ? dmom * dV : zero};
            });
        }
        ReduceTuple hv = reduce_data.value(reduce_op);
        max_du    = std::max(max_du, get<0>(hv));
        sum_dmom += get<1>(hv);
    }

    if (verbose > 0) {
        ParallelDescriptor::ReduceRealMax(max_du);
        ParallelDescriptor::ReduceRealSum(sum_dmom);
        Print() << "Momentum reflux at level " << c << ": max |du| = " << max_du
                << ", momentum to be added on the coarse/fine faces = " << sum_dmom << "\n";
    }
}

/**
 * Add the impulse of the last momentum reflux onto level lev (momentum_reflux) to its x- and
 * y-momentum sources, as a constant rate over the steps of level lev in this coarse step of
 * level 0
 *
 * @param lev      the level
 * @param xmom_src the x-momentum source
 * @param ymom_src the y-momentum source
 */
void
ERF::momentum_reflux_add_sources (int lev, MultiFab& xmom_src, MultiFab& ymom_src)
{
    if (static_cast<int>(m_mom_reflux_src.size()) <= lev) { return; }
    for (int dir = 0; dir < 2; ++dir) {
        const MultiFab& D = m_mom_reflux_src[lev][dir];
        MultiFab& src = (dir == 0) ? xmom_src : ymom_src;
        if (D.isDefined() && D.boxArray() == src.boxArray() &&
            D.DistributionMap() == src.DistributionMap()) {
            MultiFab::Saxpy(src, one / static_cast<Real>(dt[0]), D, 0, 0, 1, 0);
        }
    }
}

/**
 * Start the sums of the next coarse step (after all levels have been refluxed)
 */
void
ERF::momentum_reflux_reset ()
{
    for (int lev = 0; lev < static_cast<int>(m_mom_reflux_p.size()); ++lev) {
        if (m_mom_reflux_p[lev].isDefined()) { m_mom_reflux_p[lev].setVal(0); }
        m_mom_reflux_used[lev]  = 0;
        m_mom_reflux_stale[lev] = 0;
    }
}
