/**
 * \file ERF_InitForEnsemble.cpp
 */

#include <ERF.H>
#include <ERF_TileNoZ.H>

using namespace amrex;

void
ERF::create_random_perturbations(const int lev,
                                 MultiFab& cons_pert,
                                 MultiFab& xvel_pert,
                                 MultiFab& yvel_pert,
                                 MultiFab& zvel_pert)
{
    ignore_unused(cons_pert);
    ignore_unused(yvel_pert);
    ignore_unused(zvel_pert);

    auto& lev_new = vars_new[lev];
    for (MFIter mfi(lev_new[Vars::cons], TileNoZ()); mfi.isValid(); ++mfi) {
        const auto &xvel_pert_arr = xvel_pert.array(mfi);
        const Box &xbx = mfi.tilebox(IntVect(1,0,0));
        ParallelForRNG(xbx, [=] AMREX_GPU_DEVICE(int i, int j, int k, const RandomEngine& engine) noexcept
        {
            xvel_pert_arr(i, j, k) = Random(engine);
        });
    }
}

void
ERF::apply_gaussian_smoothing_to_perturbations(const int lev,
                                               MultiFab& cons_pert,
                                               MultiFab& xvel_pert,
                                               MultiFab& yvel_pert,
                                               MultiFab& zvel_pert)
{
    ignore_unused(cons_pert);
    ignore_unused(yvel_pert);
    ignore_unused(zvel_pert);

    const Geometry& gm = geom[lev];
    const Real dx = gm.CellSize(0);
    const Real dy = gm.CellSize(1);

    const Real dmesh = std::min(dx, dy);
    // ---- User choices ----
    const Real sigma = solverChoice.pert_correlated_radius; // e.g. 2 km correlation length
    const int  r     = static_cast<int>(3.0 * sigma / dmesh);  // stencil radius

    // ---- Precompute Gaussian weights on host ----
    const int wsize = 2*r + 1;
    Vector<Real> w_host(wsize * wsize);

    Real Z = 0.0;
    for (int m = -r; m <= r; ++m) {
        for (int n = -r; n <= r; ++n) {
            Real val = std::exp(-(m*m*dx*dx + n*n*dy*dy)/(2.0*sigma*sigma));
            w_host[(m+r)*wsize + (n+r)] = val;
            Z += val;
        }
    }
    for (auto& v : w_host) {
        v = v/Z;
    }

    Gpu::DeviceVector<Real> w_dev;
    w_dev.resize(w_host.size());
    Gpu::copy(Gpu::hostToDevice, w_host.begin(), w_host.end(), w_dev.begin());

    Real const* w = w_dev.data();

    // 1. Define ngrow_big using the actual dimension macro
    IntVect ngrow_big(AMREX_D_DECL(r, r, 0));

    // 2. Create the copy
    MultiFab xvel_pert_copy(xvel_pert.boxArray(),
                        xvel_pert.DistributionMap(),
                        1, ngrow_big);
    //MultiFab::Copy(xvel_pert_copy, xvel_pert, 0, 0, 1, 0);

    // 3. Use the built-in copy that includes ghost cell logic
    // Copy(dst, src, src_comp, dst_comp, num_comp, ngrow)
    // Setting ngrow to 0 ensures we only take valid data from the original
    xvel_pert_copy.ParallelCopy(xvel_pert, 0, 0, 1, IntVect(0), ngrow_big, gm.periodicity());

    for (MFIter mfi(xvel_pert, TileNoZ()); mfi.isValid(); ++mfi)
    {
        const Box& tbx = mfi.tilebox();

        auto const& in  = xvel_pert_copy.array(mfi);
        auto const& out = xvel_pert.array(mfi);

        ParallelFor(tbx,
        [=] AMREX_GPU_DEVICE (int i, int j, int k) noexcept
        {
            Real sum = 0.0;
            for (int m = -r; m <= r; ++m) {
                for (int n = -r; n <= r; ++n) {
                    Real wij = w[(m+r)*wsize + (n+r)];
                    sum += wij * in(i+m, j+n, k);
                }
            }
            out(i,j,k) = sum;
        });
    }
}

// Reads the plotfile data into cell cenetred multifab
// Does not fill ghost cells
void
read_plot_file(PlotFileData& pf,
               const std::string varnames,
               MultiFab& mf)
{
    // ------------------------------------------------------------
    // Open plotfile
    // ------------------------------------------------------------
    const Vector<std::string>& var_names_pf = pf.varNames();

    // ------------------------------------------------------------
    // Validate requested variables
    // ------------------------------------------------------------
    for (auto const& v : varnames) {
        bool found = false;
        for (auto const& vpf : var_names_pf) {
            if (v == vpf) {
                found = true;
                break;
            }
        }
        if (!found) {
            Abort("read_plot_file: invalid variable name: " + v);
        }
    }

    // ------------------------------------------------------------
    // Define destination MultiFab (single level only)
    // ------------------------------------------------------------
    const int level = 0;

    BoxArray ba = pf.boxArray(level);
    DistributionMapping dm(ba);

    int ncomp = varnames.size();
    int ngrow = 1;

    mf.define(ba, dm, ncomp, ngrow);

    // ------------------------------------------------------------
    // Copy plotfile data → mf
    // ------------------------------------------------------------
    for (int comp = 0; comp < ncomp; ++comp)
    {
        const MultiFab& src = pf.get(level, varnames[comp]);
        MultiFab::Copy(mf, src, 0, comp, 1, 0);
    }
}

void
create_nodal_mf_from_cc_mf (MultiFab& mf_nc,      // output nodal MF
                            MultiFab& mf_cc,      // input cell-centered MF (coarse)
                            const Geometry& geom)
{

    // -------------------------------------------------
    // 1. Build nodal MultiFab if not already defined
    // -------------------------------------------------
    if (!mf_nc.isDefined())
    {
        BoxArray ba_nd = amrex::convert(mf_cc.boxArray(),
                                        IntVect::TheNodeVector());

        mf_nc.define(ba_nd,
                     mf_cc.DistributionMap(),
                     mf_cc.nComp(),
                     0);   // nodal MF typically needs no ghosts
    }
    // -------------------------------------------------
    // 2. Fill interior + periodic ghost cells
    // -------------------------------------------------
    mf_cc.FillBoundary(geom.periodicity());
    // -------------------------------------------------
    // 3. Apply FOExtrap (Neumann) at domain boundaries
    // -------------------------------------------------
    const Box& domain = geom.Domain();

    for (MFIter mfi(mf_cc, TilingIfNotGPU()); mfi.isValid(); ++mfi)
    {
        const Box& gbx = mfi.growntilebox();   // includes ghost cells
        const Box& vbx = mfi.validbox();

        auto const& arr = mf_cc.array(mfi);
        int ncomp = mf_cc.nComp();

        ParallelFor(gbx, ncomp,
        [=] AMREX_GPU_DEVICE (int i, int j, int k, int n)
        {
            if (vbx.contains(i,j,k)) return;

            int ii = i;
            int jj = j;
            int kk = k;

            // Clamp to domain interior (FOExtrap)
            ii = amrex::max(domain.smallEnd(0),
                 amrex::min(i, domain.bigEnd(0)));

            jj = amrex::max(domain.smallEnd(1),
                 amrex::min(j, domain.bigEnd(1)));

            kk = amrex::max(domain.smallEnd(2),
                 amrex::min(k, domain.bigEnd(2)));

            arr(i,j,k,n) = arr(ii,jj,kk,n);
        });
    }

    AMREX_ALWAYS_ASSERT(mf_nc.nComp() == mf_cc.nComp());

    int ncomp = mf_cc.nComp();

    for (MFIter mfi(mf_nc, TilingIfNotGPU()); mfi.isValid(); ++mfi)
    {
        const Box& bx = mfi.tilebox();       // nodal valid box

        auto const& arr_nc = mf_nc.array(mfi);
        auto const& arr_cc = mf_cc.array(mfi);

        ParallelFor(bx, ncomp,
        [=] AMREX_GPU_DEVICE(int i, int j, int k, int n)
        {
            // 3D: eight surrounding cells
            arr_nc(i,j,k,n) =
                0.125 * ( arr_cc(i-1,j-1,k-1,n) + arr_cc(i  ,j-1,k-1,n)
                        + arr_cc(i-1,j  ,k-1,n) + arr_cc(i  ,j  ,k-1,n)
                        + arr_cc(i-1,j-1,k  ,n) + arr_cc(i  ,j-1,k  ,n)
                        + arr_cc(i-1,j  ,k  ,n) + arr_cc(i  ,j  ,k  ,n) );
        });
    }
}

IntVect
find_bound_idx(const Real& x, const Real& y, const Real& z,
               const BoxList& bl_coarse, const Geometry& geom_coarse,
               BoundType bound_type)
{
    const auto prob_lo_coarse  = geom_coarse.ProbLoArray();
    const auto dx_coarse       = geom_coarse.CellSizeArray();

    int i, j, k;

    if (bound_type == BoundType::Lo) {
        i = static_cast<int>(std::floor((x - prob_lo_coarse[0]) / dx_coarse[0]));
        j = static_cast<int>(std::floor((y - prob_lo_coarse[1]) / dx_coarse[1]));
        k = static_cast<int>(std::floor((z - prob_lo_coarse[2]) / dx_coarse[2]));
    } else { // BoundType::Hi
        i = static_cast<int>(std::ceil((x - prob_lo_coarse[0]) / dx_coarse[0]));
        j = static_cast<int>(std::ceil((y - prob_lo_coarse[1]) / dx_coarse[1]));
        k = static_cast<int>(std::ceil((z - prob_lo_coarse[2]) / dx_coarse[2]));
    }

    IntVect idx(i, j, k);

    for (const auto& b : bl_coarse) {
        if (b.contains(idx)) {
            return idx;
        }
    }

    amrex::Print() << x << " " << y << " " << z << " " << idx << std::endl;

    amrex::Print() << "Printing BoxList (coarse):\n";
for (const auto& b : bl_coarse) {
    amrex::Print() << b << "\n";
}

    amrex::Abort("Bound index not found in any box in BoxList!");
    return IntVect::TheZeroVector(); // unreachable if Abort
}

void
get_coarse_mf_on_fine_dmap(const Geometry& geom_coarse,
                           const Geometry& geom_fine,
                           const MultiFab& mf_nc_coarse,
                           const MultiFab& mf_cc_fine,
                           MultiFab& coarse_multifab_on_fine_dmap)
{
    BoxList bl_coarse = mf_nc_coarse.boxArray().boxList();
    BoxList bl_fine   = mf_cc_fine.boxArray().boxList();

    const auto prob_lo_fine  = geom_fine.ProbLoArray();
    const auto dx_fine       = geom_fine.CellSizeArray();

    for (auto& b : bl_fine) {
        // You look at the lo corner of b, and find out the lowest cell in
        // coarse mutlifab data you need for the interpolation. That gives
        // you the lo corner of the new b. Similarly, you can find out the
        // hi corner of the new b. For cells outside the coarse multifab data's
        // bounding data, it's up to you. You probably want to use a biased
        // interpolation stencil.

        // Get the cell indices of the bottom corner and top corner
        const IntVect& lo_fine = b.smallEnd();  // Lower corner (inclusive)
        const IntVect& hi_fine = b.bigEnd();    // Upper corner (inclusive)

        Real x = prob_lo_fine[0] + lo_fine[0] * dx_fine[0];
        Real y = prob_lo_fine[1] + lo_fine[1] * dx_fine[1];
        Real z = prob_lo_fine[2] + lo_fine[2] * dx_fine[2];

        auto idx_lo = find_bound_idx(x, y, z, bl_coarse, geom_coarse, BoundType::Lo);


        x = prob_lo_fine[0] + hi_fine[0] * dx_fine[0];
        y = prob_lo_fine[1] + hi_fine[1] * dx_fine[1];
        z = prob_lo_fine[2] + hi_fine[2] * dx_fine[2];

        auto idx_hi = find_bound_idx(x, y, z, bl_coarse, geom_coarse, BoundType::Hi);

        b.setSmall(idx_lo);
        b.setBig(idx_hi);

         /*Print() << "lo fine = " << lo_fine << std::endl;
         Print() << "hi fine = " << hi_fine << std::endl;
        Print() << " idx lo = " << idx_lo << std::endl;
        Print() << "idx_hi = " << idx_hi << std::endl;*/

    }

    BoxArray cba(std::move(bl_fine));
    cba.convert(IndexType::TheNodeType());  // <-- Make it nodal in all directions
    coarse_multifab_on_fine_dmap.define(cba, mf_cc_fine.DistributionMap(), mf_nc_coarse.nComp(), 0);
    coarse_multifab_on_fine_dmap.ParallelCopy(mf_nc_coarse);
}

void
populate_fine_cc_mf_from_coarse_nodal_mf(const Geometry& geom_coarse,
                                        const Geometry& geom_fine,
                                        const MultiFab& coarse_multifab_on_fine_dmap,
                                        const MultiFab& mf_cc_fine,
                                        MultiFab& mf_cc_from_coarse)
{
    AMREX_ALWAYS_ASSERT(coarse_multifab_on_fine_dmap.ixType().nodeCentered());

    if (!mf_cc_from_coarse.isDefined())
    {
        mf_cc_from_coarse.define(mf_cc_fine.boxArray(),
                                 mf_cc_fine.DistributionMap(),
                                 coarse_multifab_on_fine_dmap.nComp(),
                                 0);
    }

    AMREX_ALWAYS_ASSERT(mf_cc_from_coarse.boxArray() == mf_cc_fine.boxArray());
    AMREX_ALWAYS_ASSERT(mf_cc_from_coarse.DistributionMap() == mf_cc_fine.DistributionMap());
    AMREX_ALWAYS_ASSERT(mf_cc_from_coarse.nComp() == coarse_multifab_on_fine_dmap.nComp());

    const auto prob_lo_coarse = geom_coarse.ProbLoArray();
    const auto dx_coarse = geom_coarse.CellSizeArray();

    const auto prob_lo_fine = geom_fine.ProbLoArray();
    const auto dx_fine = geom_fine.CellSizeArray();

    Box nodal_domain = amrex::convert(geom_coarse.Domain(), IntVect::TheNodeVector());
    const auto nd_lo = nodal_domain.smallEnd();
    const auto nd_hi = nodal_domain.bigEnd();

    int ncomp = mf_cc_from_coarse.nComp();

    for (MFIter mfi(mf_cc_from_coarse, TilingIfNotGPU()); mfi.isValid(); ++mfi)
    {
        const Box& cbox = mfi.tilebox();
        auto const& arr_cc = mf_cc_from_coarse.array(mfi);
        auto const& arr_nc = coarse_multifab_on_fine_dmap.array(mfi);

        ParallelFor(cbox, ncomp,
        [=] AMREX_GPU_DEVICE(int i, int j, int k, int n)
        {
            Real x = prob_lo_fine[0] + (static_cast<Real>(i) + 0.5_rt) * dx_fine[0];
            Real y = prob_lo_fine[1] + (static_cast<Real>(j) + 0.5_rt) * dx_fine[1];
            Real z = prob_lo_fine[2] + (static_cast<Real>(k) + 0.5_rt) * dx_fine[2];

            Real fx = (x - prob_lo_coarse[0]) / dx_coarse[0];
            Real fy = (y - prob_lo_coarse[1]) / dx_coarse[1];
            Real fz = (z - prob_lo_coarse[2]) / dx_coarse[2];

            int i0 = static_cast<int>(amrex::Math::floor(fx));
            int j0 = static_cast<int>(amrex::Math::floor(fy));
            int k0 = static_cast<int>(amrex::Math::floor(fz));

            Real wx = fx - static_cast<Real>(i0);
            Real wy = fy - static_cast<Real>(j0);
            Real wz = fz - static_cast<Real>(k0);

            int i1 = i0 + 1;
            int j1 = j0 + 1;
            int k1 = k0 + 1;

            i0 = amrex::max(nd_lo[0], amrex::min(i0, nd_hi[0]-1));
            j0 = amrex::max(nd_lo[1], amrex::min(j0, nd_hi[1]-1));
            k0 = amrex::max(nd_lo[2], amrex::min(k0, nd_hi[2]-1));

            i1 = amrex::max(nd_lo[0], amrex::min(i1, nd_hi[0]));
            j1 = amrex::max(nd_lo[1], amrex::min(j1, nd_hi[1]));
            k1 = amrex::max(nd_lo[2], amrex::min(k1, nd_hi[2]));

            Real c000 = arr_nc(i0,j0,k0,n);
            Real c100 = arr_nc(i1,j0,k0,n);
            Real c010 = arr_nc(i0,j1,k0,n);
            Real c110 = arr_nc(i1,j1,k0,n);
            Real c001 = arr_nc(i0,j0,k1,n);
            Real c101 = arr_nc(i1,j0,k1,n);
            Real c011 = arr_nc(i0,j1,k1,n);
            Real c111 = arr_nc(i1,j1,k1,n);

            Real c00 = c000 * (1.0_rt - wx) + c100 * wx;
            Real c10 = c010 * (1.0_rt - wx) + c110 * wx;
            Real c01 = c001 * (1.0_rt - wx) + c101 * wx;
            Real c11 = c011 * (1.0_rt - wx) + c111 * wx;

            Real c0 = c00 * (1.0_rt - wy) + c10 * wy;
            Real c1 = c01 * (1.0_rt - wy) + c11 * wy;

            arr_cc(i,j,k,n) = c0 * (1.0_rt - wz) + c1 * wz;
        });
    }
}

void
ERF::create_background_state_for_ensemble ()
{
    bckgnd_state.resize(max_level+1);
    for (int lev = 0; lev < max_level+1; ++lev) {
        bckgnd_state[lev].resize(vars_new[lev].size()+1);
        for (int comp = 0; comp < vars_new[lev].size(); ++comp) {
            const MultiFab& src = vars_new[lev][comp];
            bckgnd_state[lev][comp].define(src.boxArray(), src.DistributionMap(),
                                           src.nComp(), src.nGrow());
         }
    }

    MultiFab mf_cc_coarse;
    PlotFileData pf_coarse(pltfile_bckgnd_coarse);
    Vector<std::string> varnames = {"density"};
    read_plot_file(pf_coarse, varnames, mf_cc_coarse);
    
    Geometry geom_coarse(pf_coarse.probDomain(0),
                         RealBox(pf_coarse.probLo(), pf_coarse.probHi()),
                         pf_coarse.coordSys(),
                         is_periodic);

    MultiFab mf_nc_coarse;
    create_nodal_mf_from_cc_mf(mf_nc_coarse, 
                               mf_cc_coarse, 
                               geom_coarse);

    MultiFab coarse_multifab_on_fine_dmap;
    get_coarse_mf_on_fine_dmap(geom_coarse, geom_fine,
                                mf_nc_coarse, mf_cc_fine,
                                coarse_multifab_on_fine_dmap);

    MultiFab mf_cc_fine_from_coarse;
    populate_fine_cc_mf_from_coarse_nodal_mf(geom_coarse, 
                                            geom_fine, 
                                            coarse_multifab_on_fine_dmap,
                                            mf_cc_fine, 
                                            mf_cc_fine_from_coarse);
}
