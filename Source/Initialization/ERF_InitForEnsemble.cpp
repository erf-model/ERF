/**
 * \file ERF_InitForEnsemble.cpp
 */

#include <ERF.H>
#include <ERF_TileNoZ.H>
#include <AMReX_PlotFileUtil.H>

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
    const Real sigma = solverChoice.ens_pert_correlated_radius; // e.g. 2 km correlation length
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
        Abort("Failed to open file for reading");
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

void 
populate_mf_cc_fine_from_mf_cc_coarse (const Geometry& geom_coarse,
                                       const Geometry& geom_fine,
                                       const MultiFab& coarse_mf_cc_on_fine_dmap,
                                       const MultiFab& mf_cc_fine,
                                       MultiFab& mf_cc_from_coarse)
{
}

void
populate_mf_face_fine_from_mf_cc_coarse(const Geometry& geom_coarse,
                                           const Geometry& geom_fine,
                                           const MultiFab& mf_cc_coarse,
                                           const MultiFab& mf_face_fine,
                                           MultiFab& mf_face_from_coarse,
                                           int dir) // 0=x,1=y,2=z
{
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
MakeFaceCenteredVelocities (const MultiFab& mf_cc_fine,
                           MultiFab& mf_xvel,
                           MultiFab& mf_yvel,
                           MultiFab& mf_zvel)
{
    BL_PROFILE("MakeFaceCenteredVelocities");

    const BoxArray& ba = mf_cc_fine.boxArray();
    const DistributionMapping& dm = mf_cc_fine.DistributionMap();
    int ng = mf_xvel.nGrow();

    // --- Define face-centered MultiFabs ---
    BoxArray ba_x = amrex::convert(ba, IntVect(1,0,0));
    BoxArray ba_y = amrex::convert(ba, IntVect(0,1,0));
    BoxArray ba_z = amrex::convert(ba, IntVect(0,0,1));

    mf_xvel.define(ba_x, dm, 1, ng);
    mf_yvel.define(ba_y, dm, 1, ng);
    mf_zvel.define(ba_z, dm, 1, ng);

    // --- X-faces (component 2) ---
    for (MFIter mfi(mf_xvel, TilingIfNotGPU()); mfi.isValid(); ++mfi)
    {
        const Box& bx = mfi.tilebox();
        auto const& uface = mf_xvel.array(mfi);
        auto const& cc    = mf_cc_fine.const_array(mfi);

        amrex::ParallelFor(bx, [=] AMREX_GPU_DEVICE (int i, int j, int k)
        {
            uface(i,j,k) = 0.5 * (cc(i-1,j,k,2) + cc(i,j,k,2));
        });
    }

    // --- Y-faces (component 3) ---
    for (MFIter mfi(mf_yvel, TilingIfNotGPU()); mfi.isValid(); ++mfi)
    {
        const Box& bx = mfi.tilebox();
        auto const& vface = mf_yvel.array(mfi);
        auto const& cc    = mf_cc_fine.const_array(mfi);

        amrex::ParallelFor(bx, [=] AMREX_GPU_DEVICE (int i, int j, int k)
        {
            vface(i,j,k) = 0.5 * (cc(i,j-1,k,3) + cc(i,j,k,3));
        });
    }

    // --- Z-faces (component 4) ---
    for (MFIter mfi(mf_zvel, TilingIfNotGPU()); mfi.isValid(); ++mfi)
    {
        const Box& bx = mfi.tilebox();
        auto const& wface = mf_zvel.array(mfi);
        auto const& cc    = mf_cc_fine.const_array(mfi);

        amrex::ParallelFor(bx, [=] AMREX_GPU_DEVICE (int i, int j, int k)
        {
            wface(i,j,k) = 0.5 * (cc(i,j,k-1,4) + cc(i,j,k,4));
        });
    }
}

void
ERF::create_background_state_for_ensemble (int lev,
                                           MultiFab& cons_pert,
                                           MultiFab& xvel_pert,
                                           MultiFab& yvel_pert,
                                           MultiFab& zvel_pert)
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

    std::string filename_custom = "coarse_data.bin";
    int nx_crse, ny_crse, nz_crse, ng_crse, ncomp_crse;
    Vector<Vector<Real>> data_crse;
    std::array<Real,3> problo_ext, probhi_ext;

    Vector<Real> data_rho, data_theta, data_xvel, data_yvel, data_zvel;

    ReadCustomDataFile(filename_custom,
                       nx_crse, ny_crse, nz_crse, ng_crse, ncomp_crse,
                       problo_ext, probhi_ext,
                       data_rho, data_theta, data_xvel, data_yvel, data_zvel);

    Geometry& geom_fine = geom[0];
    // Create a cell-centered multifab on the fine mesh - ie. something with the same boxarray,
    // distributed mapping, nGrow, but with 5 components
    MultiFab mf_cc_fine;
    const MultiFab& src = vars_new[0][0];
    int ngrow = src.nGrow();
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
    WriteSingleLevelPlotfile("plt_final", mf_cc_fine, varnames, geom_fine, 0.0, 0);

    MakeFaceCenteredVelocities(mf_cc_fine, xvel_pert, yvel_pert, zvel_pert);
}
