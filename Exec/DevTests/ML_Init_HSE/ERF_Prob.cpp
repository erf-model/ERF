#include "ERF_Prob.H"
#include "ERF_TerrainMetrics.H"

using namespace amrex;

std::unique_ptr<ProblemBase>
amrex_probinit (const amrex_real* /*problo*/,
                const amrex_real* /*probhi*/)
{
    return std::make_unique<Problem>();
}

Problem::Problem ()
{
  // Parse params
  ParmParse pp("prob");
  pp.query("rho_0", parms.rho_0);
  pp.query("T_0", parms.T_0);
  pp.query("Th_0", parms.Th_0);
  pp.get("model_filename", parms.model_filename);
  pp.get("model_unscale_filename", parms.model_unscale_filename);

  init_base_parms(parms.rho_0, parms.T_0);
}

void
Problem::init_custom_pert (const Box& bx,
                           const Box& xbx,
                           const Box& ybx,
                           const Box& zbx,
                           Array4<Real const> const& /*state*/,
                           Array4<Real      > const& state_pert,
                           Array4<Real      > const& x_vel_pert,
                           Array4<Real      > const& y_vel_pert,
                           Array4<Real      > const& z_vel_pert,
                           Array4<Real      > const& /*r_hse*/,
                           Array4<Real      > const& /*p_hse*/,
                           Array4<Real const> const& z_nd,
                           Array4<Real const> const& /*z_cc*/,
                           GeometryData const& geomdata,
                           Array4<Real const> const& /*mf_m*/,
                           Array4<Real const> const& mf_u,
                           Array4<Real const> const& mf_v,
                           const SolverChoice& sc,
                           const int /*lev*/)
{
    const bool use_moisture = (sc.moisture_type != MoistureType::None);

    ParallelFor(bx, [=] AMREX_GPU_DEVICE(int i, int j, int k) noexcept
    {
        // Set scalar = 0 everywhere
        state_pert(i, j, k, RhoScalar_comp) = 0.0;

        if (use_moisture) {
            state_pert(i, j, k, RhoQ1_comp) = 0.0;
            state_pert(i, j, k, RhoQ2_comp) = 0.0;
        }
    });

    // Set the x-velocity
    ParallelFor(xbx, [=, parms_d=parms] AMREX_GPU_DEVICE(int i, int j, int k) noexcept
    {
        x_vel_pert(i, j, k) = 0.0;
    });

    // Set the y-velocity
    ParallelFor(ybx, [=] AMREX_GPU_DEVICE(int i, int j, int k) noexcept
    {
        y_vel_pert(i, j, k) = 0.0;
    });

    const auto dx = geomdata.CellSize();
    amrex::GpuArray<Real, AMREX_SPACEDIM> dxInv;
    dxInv[0] = 1. / dx[0];
    dxInv[1] = 1. / dx[1];
    dxInv[2] = 1. / dx[2];

    // Set the z-velocity from impenetrable condition
    if (sc.terrain_type == TerrainType::StaticFittedMesh) {
        ParallelFor(zbx, [=] AMREX_GPU_DEVICE(int i, int j, int k) noexcept
        {
            z_vel_pert(i, j, k) = WFromOmega(i, j, k, 0.0, x_vel_pert, y_vel_pert,
                                             mf_u, mf_v, z_nd, dxInv);
        });
    } else {
        ParallelFor(zbx, [=] AMREX_GPU_DEVICE(int i, int j, int k) noexcept
        {
            z_vel_pert(i, j, k) = 0.0;
        });
    }

    amrex::Gpu::streamSynchronize();
}

void
Problem::init_custom_terrain (const Geometry& geom,
                              FArrayBox& terrain_fab,
                              const Real& /*time*/)
{
    // Domain cell size and real bounds
    auto dx = geom.CellSizeArray();
    auto ProbLoArr = geom.ProbLoArray();
    auto ProbHiArr = geom.ProbHiArray();

    const amrex::Box& domain = geom.Domain();
    int domlo_x = domain.smallEnd(0); int domhi_x = domain.bigEnd(0) + 1;
    int domlo_z = domain.smallEnd(2);

    // User function parameters
    Real xcen = 0.5 * (ProbLoArr[0] + ProbHiArr[0]);

    Real asq    = 5000.0 * 5000.0;
    Real Hm     =  250.0;
    Real lambda = 4000.0;

    // Populate bottom plane
    int k0 = domlo_z;

    amrex::Box zbx = terrain_fab.box();
    if (zbx.smallEnd(2) <= k0)
    {
        amrex::Array4<Real> const& z_arr = terrain_fab.array();

        ParallelFor(zbx, [=] AMREX_GPU_DEVICE (int i, int j, int)
        {
            // Clip indices for ghost-cells
            int ii = amrex::min(amrex::max(i,domlo_x),domhi_x);

            // Location of nodes
            Real x = (ProbLoArr[0] + ii * dx[0] - xcen);

            Real cosx = cos(PI * x / lambda);

            z_arr(i,j,k0) = Hm * std::exp(-x*x/asq) * cosx * cosx;
        });
    }
}

void
Problem::erf_init_dens_hse (amrex::MultiFab& rho_hse,
                            std::unique_ptr<amrex::MultiFab>& /*z_phys_nd*/,
                            std::unique_ptr<amrex::MultiFab>& z_phys_cc,
                            amrex::Geometry const& geom)
{
    // Set up pytorch model
    //====================================
    torch::jit::script::Module module;
    try {
        // Deserialize the ScriptModule from a file using torch::jit::load().
        Print() << "Attempting to load model file: " << parms.model_filename << "\n";
        module = torch::jit::load(parms.model_filename);
    }
    catch (const c10::Error& e) {
        Abort("Error loading the model\n");
    }

    Print() << "Model loaded.\n";

    // set pytorch data type (default is float or torch::kFloat32)
    auto dtype0 = torch::kFloat64;

#ifdef AMREX_USE_CUDA
    torch::Device device0(torch::kCUDA);
    module.to(device0);
    Print() << "Copying model to GPU." << std::endl;

    // set tensor options
    auto tensoropt = torch::TensorOptions().dtype(dtype0).device(device0);
#else
    auto tensoropt = torch::TensorOptions().dtype(dtype0);
#endif

    // Custom minmax scaling (needed for quality model)
    int nscale = 128;
    Vector<Real> Theta(nscale), z(nscale);         // INPUTS
    Vector<Real> rho_min(nscale), rho_max(nscale); // OUTPUTS
    std::ifstream file(parms.model_unscale_filename, std::ios::binary);
    file.read(reinterpret_cast<char*>(Theta.data())  , Theta.size()   * sizeof(Real));
    file.read(reinterpret_cast<char*>(z.data())      , z.size()       * sizeof(Real));
    file.read(reinterpret_cast<char*>(rho_min.data()), rho_min.size() * sizeof(Real));
    file.read(reinterpret_cast<char*>(rho_max.data()), rho_max.size() * sizeof(Real));
    auto thmin = Theta[0]; auto thmax = Theta[nscale-1];
    auto zmin  = z[0];     auto zmax  = z[nscale-1];
    auto rmin  = interpolate_1d(Theta.data(), rho_min.data(), parms.Th_0, nscale);
    auto rmax  = interpolate_1d(Theta.data(), rho_max.data(), parms.Th_0, nscale);


    // Populate the HSE state
    //====================================
    auto Th_hse = parms.Th_0;
    auto dz  = geom.CellSize(2);
    auto zlo = geom.ProbLo(2);
    auto khi = geom.Domain().bigEnd(2);
    auto klo = geom.Domain().smallEnd(2);
#ifdef _OPENMP
#pragma omp parallel if (Gpu::notInLaunchRegion())
#endif
    for ( MFIter mfi(rho_hse); mfi.isValid(); ++mfi)
    {
        const Box& bx = mfi.growntilebox();

        // Dycore data
        Array4<Real> rho_hse_arr = rho_hse.array(mfi);
        Array4<Real> zcc_arr     = (z_phys_cc) ? z_phys_cc->array(mfi) : Array4<Real>{};

        // Auxiliary array for pytorch
        int nin   = 2; // theta, z
        int nout  = 2; // P, rho
        int ncell = bx.numPts();
        Gpu::ManagedVector<Real> ML_aux(ncell*nin);
        Real* AMREX_RESTRICT ML_auxPtr = ML_aux.dataPtr();

        // Get smallend and size of box
        const IntVect bx_lo = bx.smallEnd();
        const IntVect nbox  = bx.size();

        // Copy the ML inputs into auxiliary array
        ParallelFor(bx, nin, [=] AMREX_GPU_DEVICE(int i, int j, int k, int n) noexcept
        {
            // Flatten indexing
            int ii = i - bx_lo[0];
            int jj = j - bx_lo[1];
            int index = jj*nbox[0] + ii;
            int kk = k - bx_lo[2];
            index += kk*nbox[0]*nbox[1];

            // Ensure the z-value is FOEXTRAP
            int  k_hse  = std::min(std::max(k,klo),khi);
            Real z_hse  = (zcc_arr) ? zcc_arr(i,j,k_hse) : zlo + dz*(Real(k_hse) + 0.5);

            // array order is row-based [index][comp]
            if (n==0) {
                ML_auxPtr[index*nin + n] = (Th_hse - thmin) / (thmax - thmin);
            } else {
                ML_auxPtr[index*nin + n] = (z_hse - zmin) / (zmax - zmin);
            }
        });

        // Create torch tensor from array
        at::Tensor inputs_torch = torch::from_blob(ML_auxPtr, {ncell, nin}, tensoropt);

        // Evaluate torch model
        at::Tensor outputs_torch = module.forward({inputs_torch}).toTensor();
        outputs_torch = outputs_torch.to(dtype0);

        // get accessor to tensor (read-only)
#ifdef AMREX_USE_CUDA
        auto outputs_torch_acc = outputs_torch.packed_accessor64<Real,2>();
#else
        auto outputs_torch_acc = outputs_torch.accessor<Real,2>();
#endif

        // Copy ML outputs to the dycore vars
        ParallelFor(bx, [=] AMREX_GPU_DEVICE(int i, int j, int k) noexcept
        {
            // NOTE: We only assign the density here.
            //       The call pathway (from ERF_Init1D.cpp) then goes
            //       to "erf_enforce_hse" which does the pressure integration
            //       and sets theta according to rho and pressure.

            // Flatten indexing
            int ii = i - bx_lo[0];
            int jj = j - bx_lo[1];
            int index = jj*nbox[0] + ii;
            int kk = k - bx_lo[2];
            index += kk*nbox[0]*nbox[1];

            // ML output is {P, Rho}
            int n = 1;
            rho_hse_arr(i,j,k) = outputs_torch_acc[index][n]*(rmax - rmin) + rmin;
        });
    } // mfi

    rho_hse.FillBoundary(geom.periodicity());
}
