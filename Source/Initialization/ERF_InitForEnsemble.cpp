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
        ParallelForRNG(xbx, [=] AMREX_GPU_DEVICE(int i, int j, int k, const amrex::RandomEngine& engine) noexcept
        {
            xvel_pert_arr(i, j, k) = amrex::Random(engine);
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
