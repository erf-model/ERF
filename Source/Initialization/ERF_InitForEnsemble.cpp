/**
 * \file ERF_InitForEnsemble.cpp
 */

#include <ERF.H>
#include <ERF_TileNoZ.H>
#include <AMReX_PlotFileUtil.H>
#include <filesystem>

using namespace amrex;
namespace fs = std::filesystem;

void
ERF::create_random_perturbations(const int lev,
                                 MultiFab& mf_cc_pert)
{
    const MultiFab& src = vars_new[lev][Vars::cons];

    int ncomp = 5;
    mf_cc_pert.define(src.boxArray(), src.DistributionMap(),
                      ncomp, src.nGrow());

    // Loop over cell-centered boxes
    for (MFIter mfi(mf_cc_pert, TilingIfNotGPU()); mfi.isValid(); ++mfi)
    {
        const Box& bx = mfi.tilebox();

        auto const& pert_arr = mf_cc_pert.array(mfi);

        // Loop over all 5 components
        amrex::ParallelForRNG(bx, ncomp,
        [=] AMREX_GPU_DEVICE (int i, int j, int k, int n,
                             const amrex::RandomEngine& engine) noexcept
        {
            pert_arr(i,j,k,n) = amrex::Random(engine);
        });
    }
}

void NormalizeMultiFabRMS_PerComponent(MultiFab& mf_cc_pert)
{
    const int ncomp = mf_cc_pert.nComp();

    for (int n = 0; n < ncomp; ++n)
    {
        // 1. Set up AMReX reduction (sum of squares + count)
        ReduceOps<ReduceOpSum, ReduceOpSum> reduce_op;
        ReduceData<Real, Long> reduce_data(reduce_op);
        using ReduceTuple = typename decltype(reduce_data)::Type;

        // 2. Loop over tiles and accumulate
        for (MFIter mfi(mf_cc_pert, TilingIfNotGPU()); mfi.isValid(); ++mfi)
        {
            const Box& bx = mfi.tilebox();
            auto const& arr = mf_cc_pert.const_array(mfi);

            reduce_op.eval(bx, reduce_data,
                [=] AMREX_GPU_DEVICE (int i, int j, int k) noexcept -> ReduceTuple
                {
                    Real v = arr(i, j, k, n);
                    return { v * v, 1L };
                });
        }

        // 3. Retrieve results (includes implicit GPU sync + host copy)
        auto rv      = reduce_data.value(reduce_op);
        Real h_sumsq = amrex::get<0>(rv);
        Long h_count = amrex::get<1>(rv);

        // 4. Sum across MPI ranks
        ParallelDescriptor::ReduceRealSum(h_sumsq);
        ParallelDescriptor::ReduceLongSum(h_count);

        // 5. Compute RMS and normalize
        if (h_count > 0)
        {
            Real rms = std::sqrt(h_sumsq / static_cast<Real>(h_count));
            if (rms > 0.0) {
                mf_cc_pert.mult(1.0 / rms, n, 1);
            }
        }
    }
}

void
ERF::apply_gaussian_smoothing_to_perturbations(const int lev,
                                               MultiFab& mf_cc_pert)
{
    const Geometry& gm = geom[lev];
    const Real dx = gm.CellSize(0);
    const Real dy = gm.CellSize(1);

    const Real dmesh = std::min(dx, dy);

    // ---- User choice ----
    const Real sigma = solverChoice.ens_pert_correlated_radius;
    const int  r     = static_cast<int>(3.0 * sigma / dmesh);

    const int ncomp = mf_cc_pert.nComp();

    // ---- Precompute Gaussian weights ----
    const int wsize = 2*r + 1;
    Vector<Real> w_host(wsize * wsize);

    Real Z = 0.0;
    for (int m = -r; m <= r; ++m) {
        for (int n = -r; n <= r; ++n) {
            Real val = std::exp(-(m*m*dx*dx + n*n*dy*dy)
                                 /(2.0*sigma*sigma));
            w_host[(m+r)*wsize + (n+r)] = val;
            Z += val;
        }
    }

    for (auto& v : w_host) {
        v /= Z;
    }

    Gpu::DeviceVector<Real> w_dev(w_host.size());
    Gpu::copy(Gpu::hostToDevice, w_host.begin(), w_host.end(), w_dev.begin());

    Real const* w = w_dev.data();

    // ---- Create a grown copy (for stencil access) ----
    IntVect ngrow_big(AMREX_D_DECL(r, r, 0));

    MultiFab mf_copy(mf_cc_pert.boxArray(),
                     mf_cc_pert.DistributionMap(),
                     ncomp, ngrow_big);

    mf_copy.ParallelCopy(mf_cc_pert,
                         0, 0, ncomp,
                         IntVect(0), ngrow_big,
                         gm.periodicity());

    // ---- Apply smoothing ----
    for (MFIter mfi(mf_cc_pert, TilingIfNotGPU()); mfi.isValid(); ++mfi)
    {
        const Box& bx = mfi.tilebox();

        auto const& in  = mf_copy.const_array(mfi);
        auto const& out = mf_cc_pert.array(mfi);

        ParallelFor(bx, ncomp,
        [=] AMREX_GPU_DEVICE (int i, int j, int k, int n) noexcept
        {
            Real sum = 0.0;

            for (int m = -r; m <= r; ++m) {
                for (int nn = -r; nn <= r; ++nn) {
                    Real wij = w[(m+r)*wsize + (nn+r)];
                    sum += wij * in(i+m, j+nn, k, n);
                }
            }

            out(i,j,k,n) = sum;
        });
    }
    NormalizeMultiFabRMS_PerComponent(mf_cc_pert);
}

void ApplyNeumannBCs(const Geometry& geom,
                     MultiFab& mf_cc)
{

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
}


void ReadCustomDataFile(const std::string& filename_custom,
                        int& nx, int& ny, int& nz,
                        int& ng, int& ncomp,
                        std::array<Real,3>& problo_ext,
                        std::array<Real,3>& probhi_ext,
                        Vector<Real>& data_rho,
                        Vector<Real>& data_theta,
                        Vector<Real>& data_xvel,
                        Vector<Real>& data_yvel,
                        Vector<Real>& data_zvel)
{
    std::ifstream ifs(filename_custom, std::ios::binary);
    if (!ifs.is_open()) {
        Abort("Failed to open file " + filename_custom + " for reading");
    }

    // ----------------------------
    // Read header
    // ----------------------------
    ifs.read(reinterpret_cast<char*>(&nx), sizeof(int));
    ifs.read(reinterpret_cast<char*>(&ny), sizeof(int));
    ifs.read(reinterpret_cast<char*>(&nz), sizeof(int));

    ifs.read(reinterpret_cast<char*>(&ng), sizeof(int));
    ifs.read(reinterpret_cast<char*>(&ncomp), sizeof(int));

    ifs.read(reinterpret_cast<char*>(&problo_ext[0]), sizeof(Real));
    ifs.read(reinterpret_cast<char*>(&problo_ext[1]), sizeof(Real));
    ifs.read(reinterpret_cast<char*>(&problo_ext[2]), sizeof(Real));

    ifs.read(reinterpret_cast<char*>(&probhi_ext[0]), sizeof(Real));
    ifs.read(reinterpret_cast<char*>(&probhi_ext[1]), sizeof(Real));
    ifs.read(reinterpret_cast<char*>(&probhi_ext[2]), sizeof(Real));

    const std::size_t ncell = static_cast<std::size_t>(nx) * ny * nz;

    // ----------------------------
    // Allocate storage
    // ----------------------------
    data_rho.resize(ncell);
    data_theta.resize(ncell);
    data_xvel.resize(ncell);
    data_yvel.resize(ncell);
    data_zvel.resize(ncell);

    // ----------------------------
    // Read data
    // ----------------------------
    std::size_t idx = 0;

    for (int k = 0; k < nz; ++k)
    {
        for (int j = 0; j < ny; ++j)
        {
            for (int i = 0; i < nx; ++i)
            {
                // Skip coordinates
                Real x, y, z;
                ifs.read(reinterpret_cast<char*>(&x), sizeof(Real));
                ifs.read(reinterpret_cast<char*>(&y), sizeof(Real));
                ifs.read(reinterpret_cast<char*>(&z), sizeof(Real));

                // Read components (fixed order)
                ifs.read(reinterpret_cast<char*>(&data_rho[idx]),   sizeof(Real));
                ifs.read(reinterpret_cast<char*>(&data_theta[idx]), sizeof(Real));
                ifs.read(reinterpret_cast<char*>(&data_xvel[idx]),  sizeof(Real));
                ifs.read(reinterpret_cast<char*>(&data_yvel[idx]),  sizeof(Real));
                ifs.read(reinterpret_cast<char*>(&data_zvel[idx]),  sizeof(Real));

                ++idx;
            }
        }
    }

    ifs.close();
}

AMREX_GPU_HOST_DEVICE
AMREX_FORCE_INLINE
int idx(int i, int j, int k, int nx, int ny)
{
    return i + nx * (j + ny * k);
}

AMREX_GPU_HOST_DEVICE
AMREX_FORCE_INLINE
Real interp_trilinear(
    const Real* f,      // <-- raw pointer
    int i, int j, int k,
    Real tx, Real ty, Real tz,
    int nx, int ny, int nz)
{
    int i1 = amrex::min(i+1, nx-1);
    int j1 = amrex::min(j+1, ny-1);
    int k1 = amrex::min(k+1, nz-1);

    Real c000 = f[idx(i ,j ,k ,nx,ny)];
    Real c100 = f[idx(i1,j ,k ,nx,ny)];
    Real c010 = f[idx(i ,j1,k ,nx,ny)];
    Real c110 = f[idx(i1,j1,k ,nx,ny)];
    Real c001 = f[idx(i ,j ,k1,nx,ny)];
    Real c101 = f[idx(i1,j ,k1,nx,ny)];
    Real c011 = f[idx(i ,j1,k1,nx,ny)];
    Real c111 = f[idx(i1,j1,k1,nx,ny)];

    Real c00 = c000*(1-tx) + c100*tx;
    Real c10 = c010*(1-tx) + c110*tx;
    Real c01 = c001*(1-tx) + c101*tx;
    Real c11 = c011*(1-tx) + c111*tx;

    Real c0 = c00*(1-ty) + c10*ty;
    Real c1 = c01*(1-ty) + c11*ty;

    return c0*(1-tz) + c1*tz;
}

void
InterpolateToFineMF(
    const Vector<Real>& data_rho,
    const Vector<Real>& data_theta,
    const Vector<Real>& data_xvel,
    const Vector<Real>& data_yvel,
    const Vector<Real>& data_zvel,
    int nx, int ny, int nz,
    const std::array<Real,3>& problo,
    const std::array<Real,3>& probhi,
    MultiFab& mf_fine,
    const Geometry& geom_fine)
{
    // coarse spacing
    Real dx_c[3];
    dx_c[0] = (probhi[0] - problo[0]) / nx;
    dx_c[1] = (probhi[1] - problo[1]) / ny;
    dx_c[2] = (probhi[2] - problo[2]) / nz;

    const auto problo_f = geom_fine.ProbLoArray();
    const auto dx_f     = geom_fine.CellSizeArray();

    // Step 1: declare device vectors with correct size
    amrex::Gpu::DeviceVector<Real> d_rho(data_rho.size());
    amrex::Gpu::DeviceVector<Real> d_theta(data_theta.size());
    amrex::Gpu::DeviceVector<Real> d_xvel(data_xvel.size());
    amrex::Gpu::DeviceVector<Real> d_yvel(data_yvel.size());
    amrex::Gpu::DeviceVector<Real> d_zvel(data_zvel.size());

    // Step 2: copy data from host to device
    amrex::Gpu::copyAsync(amrex::Gpu::hostToDevice,
                      data_rho.begin(), data_rho.end(),
                      d_rho.begin());

    amrex::Gpu::copyAsync(amrex::Gpu::hostToDevice,
                      data_theta.begin(), data_theta.end(),
                      d_theta.begin());

    amrex::Gpu::copyAsync(amrex::Gpu::hostToDevice,
                      data_xvel.begin(), data_xvel.end(),
                      d_xvel.begin());

    amrex::Gpu::copyAsync(amrex::Gpu::hostToDevice,
                      data_yvel.begin(), data_yvel.end(),
                      d_yvel.begin());

    amrex::Gpu::copyAsync(amrex::Gpu::hostToDevice,
                      data_zvel.begin(), data_zvel.end(),
                      d_zvel.begin());

    const Real* rho_ptr   = d_rho.data();
    const Real* theta_ptr = d_theta.data();
    const Real* xvel_ptr  = d_xvel.data();
    const Real* yvel_ptr  = d_yvel.data();
    const Real* zvel_ptr  = d_zvel.data();
    // -------------------------------
    // GPU kernel over MultiFab
    // -------------------------------
    for (MFIter mfi(mf_fine); mfi.isValid(); ++mfi)
    {
        const Box& bx = mfi.validbox();
        auto arr = mf_fine.array(mfi);

        amrex::ParallelFor(bx,
        [=] AMREX_GPU_DEVICE (int i, int j, int k) noexcept
        {
            // physical location (fine cell center)
            Real x = problo_f[0] + (i + 0.5) * dx_f[0];
            Real y = problo_f[1] + (j + 0.5) * dx_f[1];
            Real z = problo_f[2] + (k + 0.5) * dx_f[2];

            // map to coarse index space
            Real rx = (x - problo[0]) / dx_c[0] - 0.5;
            Real ry = (y - problo[1]) / dx_c[1] - 0.5;
            Real rz = (z - problo[2]) / dx_c[2] - 0.5;

            int ic = static_cast<int>(floor(rx));
            int jc = static_cast<int>(floor(ry));
            int kc = static_cast<int>(floor(rz));

            Real tx = rx - ic;
            Real ty = ry - jc;
            Real tz = rz - kc;

            // clamp
            ic = amrex::max(0, amrex::min(ic, nx-1));
            jc = amrex::max(0, amrex::min(jc, ny-1));
            kc = amrex::max(0, amrex::min(kc, nz-1));

            //printf("The values are x, y, z, rx, ry, rz = %0.15g %0.15g %0.15g %0.15g %0.15g %0.15g\n", x, y, z, rx, ry, rz);

            // interpolate each component using device trilinear
            arr(i,j,k,0) = interp_trilinear(rho_ptr,   ic,jc,kc, tx,ty,tz, nx,ny,nz);
            arr(i,j,k,1) = interp_trilinear(theta_ptr, ic,jc,kc, tx,ty,tz, nx,ny,nz);
            arr(i,j,k,2) = interp_trilinear(xvel_ptr,  ic,jc,kc, tx,ty,tz, nx,ny,nz);
            arr(i,j,k,3) = interp_trilinear(yvel_ptr,  ic,jc,kc, tx,ty,tz, nx,ny,nz);
            arr(i,j,k,4) = interp_trilinear(zvel_ptr,  ic,jc,kc, tx,ty,tz, nx,ny,nz);
        });
    }
}

void
MakeFinalMultiFabs (const MultiFab& mf_cc_fine,
                    MultiFab& cons_pert,
                    MultiFab& xvel_pert,
                    MultiFab& yvel_pert,
                    MultiFab& zvel_pert)
{

    for (MFIter mfi(cons_pert, TilingIfNotGPU()); mfi.isValid(); ++mfi)
    {
        const Box& bx = mfi.tilebox();

        auto const& mf_cc_fine_arr  = mf_cc_fine.const_array(mfi);
        auto const& cons_pert_arr = cons_pert.array(mfi);

        amrex::ParallelFor(bx, [=] AMREX_GPU_DEVICE (int i, int j, int k)
        {
            Real tmp_rho = mf_cc_fine_arr(i,j,k,0);
            Real tmp_theta = mf_cc_fine_arr(i,j,k,1);
            cons_pert_arr(i,j,k,Rho_comp) = tmp_rho;
            cons_pert_arr(i,j,k,RhoTheta_comp) = tmp_rho*tmp_theta;
        });
    }

    // --- X-faces (component 2) ---
    for (MFIter mfi(xvel_pert, TilingIfNotGPU()); mfi.isValid(); ++mfi)
    {
        const Box& bx = mfi.tilebox();
        auto const& uface = xvel_pert.array(mfi);
        auto const& cc    = mf_cc_fine.const_array(mfi);

        amrex::ParallelFor(bx, [=] AMREX_GPU_DEVICE (int i, int j, int k)
        {
            uface(i,j,k) = 0.5 * (cc(i-1,j,k,2) + cc(i,j,k,2));
        });
    }

    // --- Y-faces (component 3) ---
    for (MFIter mfi(yvel_pert, TilingIfNotGPU()); mfi.isValid(); ++mfi)
    {
        const Box& bx = mfi.tilebox();
        auto const& vface = yvel_pert.array(mfi);
        auto const& cc    = mf_cc_fine.const_array(mfi);

        amrex::ParallelFor(bx, [=] AMREX_GPU_DEVICE (int i, int j, int k)
        {
            vface(i,j,k) = 0.5 * (cc(i,j-1,k,3) + cc(i,j,k,3));
        });
    }

    // --- Z-faces (component 4) ---
    for (MFIter mfi(zvel_pert, TilingIfNotGPU()); mfi.isValid(); ++mfi)
    {
        const Box& bx = mfi.tilebox();
        auto const& wface = zvel_pert.array(mfi);
        //auto const& cc    = mf_cc_fine.const_array(mfi);

        amrex::ParallelFor(bx, [=] AMREX_GPU_DEVICE (int i, int j, int k)
        {
            wface(i,j,k) = 0.0;//0.5 * (cc(i,j,k-1,4) + cc(i,j,k,4));
        });
    }
}

void
AddPertToBckgnd(MultiFab& mf_cc_fine,
                const MultiFab& mf_cc_pert)
{
    const int ncomp = mf_cc_fine.nComp();

    // Optional safety check (recommended)
    AMREX_ALWAYS_ASSERT(mf_cc_pert.nComp() == ncomp);

    for (MFIter mfi(mf_cc_fine, TilingIfNotGPU()); mfi.isValid(); ++mfi)
    {
        const Box& bx = mfi.tilebox();

        auto const& bg   = mf_cc_fine.array(mfi);
        auto const& pert = mf_cc_pert.const_array(mfi);

        amrex::ParallelFor(bx, ncomp,
        [=] AMREX_GPU_DEVICE (int i, int j, int k, int n) noexcept
        {
            Real ens_amp = 0.02*std::abs(bg(i,j,k,n));
            bg(i,j,k,n) += ens_amp*pert(i,j,k,n);
        });
    }
}

void
ERF::create_background_state_for_ensemble (int lev,
                                           MultiFab& mf_cc_pert,
                                           MultiFab& cons_pert,
                                           MultiFab& xvel_pert,
                                           MultiFab& yvel_pert,
                                           MultiFab& zvel_pert)
{

    ignore_unused(lev);
    int nx_crse, ny_crse, nz_crse, ng_crse, ncomp_crse;
    Vector<Vector<Real>> data_crse;
    std::array<Real,3> problo_ext, probhi_ext;

    Vector<Real> data_rho, data_theta, data_xvel, data_yvel, data_zvel;

    ReadCustomDataFile(solverChoice.coarse_bckgnd_data_file,
                       nx_crse, ny_crse, nz_crse, ng_crse, ncomp_crse,
                       problo_ext, probhi_ext,
                       data_rho, data_theta, data_xvel, data_yvel, data_zvel);

    Geometry& geom_fine = geom[0];
    // Create a cell-centered multifab on the fine mesh - ie. something with the same boxarray,
    // distributed mapping, nGrow, but with 5 components
    MultiFab mf_cc_fine;
    const MultiFab& src = vars_new[0][0];
    int ncomp = 5;
    mf_cc_fine.define(src.boxArray(), src.DistributionMap(),
                                           ncomp, src.nGrow());

    InterpolateToFineMF(data_rho, data_theta, data_xvel, data_yvel, data_zvel,
                        nx_crse, ny_crse, nz_crse,
                        problo_ext, probhi_ext,
                        mf_cc_fine,
                        geom_fine);

    ApplyNeumannBCs(geom_fine, mf_cc_fine);

    Vector<std::string> varnames = {"density","theta", "x_velocity","y_velocity","z_velocity"};

     // Add pertubrations stored in the "pert" variables in the function arguments
    // (multiplied by the corresponding amplitude)
    AddPertToBckgnd(mf_cc_fine, mf_cc_pert);
    ApplyNeumannBCs(geom_fine, mf_cc_fine);
    //WriteSingleLevelPlotfile("1_plt_final", mf_cc_fine, varnames, geom_fine, 0.0, 0);

    MakeFinalMultiFabs(mf_cc_fine, cons_pert, xvel_pert, yvel_pert, zvel_pert);
}
