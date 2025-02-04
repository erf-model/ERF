#include "ERF_Prob.H"
#include <ERF_Constants.H>

using namespace amrex;

std::unique_ptr<ProblemBase>
amrex_probinit (
    const amrex_real* /*problo*/,
    const amrex_real* /*probhi*/)
{
    return std::make_unique<Problem>();
}

Problem::Problem ()
{
  // Parse params
  amrex::ParmParse pp("prob");
}

void
Problem::erf_init_dens_hse_moist (MultiFab& rho_hse,
                                  std::unique_ptr<MultiFab>& z_phys_nd,
                                  Geometry const& geom)
{
    const Real prob_lo_z = geom.ProbLo()[2];
    const Real dz        = geom.CellSize()[2];
    const int khi        = geom.Domain().bigEnd()[2];

    // use_terrain = 1
    if (z_phys_nd) {

        if (khi > 255) amrex::Abort("1D Arrays are hard-wired to only 256 high");

        for ( MFIter mfi(rho_hse, TileNoZ()); mfi.isValid(); ++mfi )
        {
            Array4<Real      > rho_arr  = rho_hse.array(mfi);
            //Array4<Real const> z_cc_arr = z_phys_cc->const_array(mfi);

            // Create a flat box with same horizontal extent but only one cell in vertical
            const Box& tbz = mfi.nodaltilebox(2);
            Box b2d = tbz; // Copy constructor
            b2d.grow(0,1); b2d.grow(1,1); // Grow by one in the lateral directions
            b2d.setRange(2,0);

            ParallelFor(b2d, [=] AMREX_GPU_DEVICE (int i, int j, int) {
              Array1D<Real,0,255> r;

              //init_isentropic_hse_terrain(i,j,rho_sfc,Thetabar,&(r(0)),&(p(0)),z_cc_arr,khi);

              for (int k = 0; k <= khi; k++) {
                 rho_arr(i,j,k) = r(k);
              }
              rho_arr(i,j,   -1) = rho_arr(i,j,0);
              rho_arr(i,j,khi+1) = rho_arr(i,j,khi);
            });
        } // mfi
    } else { // use_terrain = 0

        // These are at cell centers (unstaggered)
        Vector<Real> h_r(khi+2);
        Vector<Real> h_p(khi+2);
        Vector<Real> h_t(khi+2);
        Vector<Real> h_q_v(khi+2);

        amrex::Gpu::DeviceVector<Real> d_r(khi+2);

        amrex::Gpu::copyAsync(amrex::Gpu::hostToDevice, h_r.begin(), h_r.end(), d_r.begin());

        Real* r     = d_r.data();

#ifdef _OPENMP
#pragma omp parallel if (amrex::Gpu::notInLaunchRegion())
#endif
          for ( MFIter mfi(rho_hse,TilingIfNotGPU()); mfi.isValid(); ++mfi)
          {
              const Box& bx = mfi.growntilebox(1);
              const Array4<Real> rho_hse_arr   = rho_hse[mfi].array();
              ParallelFor(bx, [=] AMREX_GPU_DEVICE(int i, int j, int k) noexcept
              {
                  int kk = std::max(k,0);
                  rho_hse_arr(i,j,k) = r[kk];
              });
          } // mfi
    } // no terrain
}

void
Problem::init_custom_pert (
    const Box& bx,
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
    Array4<Real const> const& /*z_nd*/,
    Array4<Real const> const& /*z_cc*/,
    GeometryData const& geomdata,
    Array4<Real const> const& /*mf_m*/,
    Array4<Real const> const& /*mf_u*/,
    Array4<Real const> const& /*mf_v*/,
    const SolverChoice& sc)
{
    const bool use_moisture = (sc.moisture_type != MoistureType::None);

    const int khi = geomdata.Domain().bigEnd()[2];

    AMREX_ALWAYS_ASSERT(bx.length()[2] == khi+1);

  // This is what we do at k = 0 -- note we assume p = p_0 and T = T_0 at z=0
  const amrex::Real& dz        = geomdata.CellSize()[2];
  const amrex::Real& prob_lo_z = geomdata.ProbLo()[2];

  // Call the routine to calculate the 1d background condition

   Vector<Real> h_r(khi+2);
   Vector<Real> h_p(khi+2);
   Vector<Real> h_t(khi+2);
   Vector<Real> h_q_v(khi+2);

   amrex::Gpu::DeviceVector<Real> d_r(khi+2);
   amrex::Gpu::DeviceVector<Real> d_p(khi+2);
   amrex::Gpu::DeviceVector<Real> d_t(khi+2);
   amrex::Gpu::DeviceVector<Real> d_q_v(khi+2);

   amrex::Gpu::copyAsync(amrex::Gpu::hostToDevice, h_r.begin(), h_r.end(), d_r.begin());
   amrex::Gpu::copyAsync(amrex::Gpu::hostToDevice, h_p.begin(), h_p.end(), d_p.begin());
   amrex::Gpu::copyAsync(amrex::Gpu::hostToDevice, h_t.begin(), h_t.end(), d_t.begin());
   amrex::Gpu::copyAsync(amrex::Gpu::hostToDevice, h_q_v.begin(), h_q_v.end(), d_q_v.begin());

    Real* t   = d_t.data();
    Real* p   = d_p.data();

 // File to read
    const std::string filename = "ERF_IC.bin";

    // Open the binary file in input mode
    std::ifstream infile(filename, std::ios::binary);
    if (!infile) {
        std::cerr << "Error: Could not open file " << filename << std::endl;
    }


	int nx, ny, nz, ndata;
    float value;

    // Read the four integers
    infile.read(reinterpret_cast<char*>(&nx), sizeof(int));
    infile.read(reinterpret_cast<char*>(&ny), sizeof(int));
    infile.read(reinterpret_cast<char*>(&nz), sizeof(int));
    infile.read(reinterpret_cast<char*>(&ndata), sizeof(int));

	
	for(int i=0; i<nx; i++) {
		infile.read(reinterpret_cast<char*>(&value), sizeof(float));
	}
	
	for(int j=0; j<ny; j++) {
		infile.read(reinterpret_cast<char*>(&value), sizeof(float));
	}

	for(int k=0; k<nz; k++) {
		infile.read(reinterpret_cast<char*>(&value), sizeof(float));
	}

    // Vector to store the data
    std::vector<float> data;

    // Read the file
	for(int idx=0; idx<ndata; idata++){
		for(int k=0; k<nz; k++) {
			for(int j=0; j<ny; j++) {
				for(int i=0; i<nx; i++) {
					infile.read(reinterpret_cast<char*>(&value), sizeof(float));
				}
			}
		}
    }

    infile.close();

    // Output the data
    std::cout << "Read " << data.size() << " elements from the file." << std::endl;
    for (size_t i = 0; i < data.size(); ++i) {
        std::cout << "Element " << i << ": " << data[i] << std::endl;
    }	
	
  ParallelFor(bx, [=] AMREX_GPU_DEVICE(int i, int j, int k) noexcept
  {
    // Geometry (note we must include these here to get the data on device)
    const auto prob_lo  = geomdata.ProbLo();
    const auto dx       = geomdata.CellSize();
    const Real x        = prob_lo[0] + (i + 0.5) * dx[0];
    const Real y        = prob_lo[1] + (j + 0.5) * dx[1];
    const Real z        = prob_lo[2] + (k + 0.5) * dx[2];

    // This version perturbs rho but not p
    state_pert(i, j, k, RhoTheta_comp) = 0.0;//rho*theta_total - rho_back*t[k]*(1.0 + (R_v/R_d)*q_v_back);// rho*d_t[k]*(1.0 + R_v_by_R_d*q_v_hot);
    state_pert(i, j, k, Rho_comp)      = 0.0;//rho - rho_back*(1.0 + q_v_back);

    // Set scalar = 0 everywhere
    state_pert(i, j, k, RhoScalar_comp) = 0.0;//rho*scalar;

    // mean states
    if (use_moisture) {
        state_pert(i, j, k, RhoQ1_comp) = 0.0;//rho*q_v_hot;
        state_pert(i, j, k, RhoQ2_comp) = 0.0;
    }

  });

  // Set the x-velocity
  ParallelFor(xbx, [=] AMREX_GPU_DEVICE(int i, int j, int k) noexcept {
      const amrex::Real z = prob_lo_z + (k+0.5) * dz;
      x_vel_pert(i,j,k) = -12.0*std::max(0.0, (2.5e3 - z)/2.5e3);
  });

  // Set the y-velocity
  ParallelFor(ybx, [=] AMREX_GPU_DEVICE(int i, int j, int k) noexcept {
      y_vel_pert(i, j, k) = 0.0;
  });

  // Set the z-velocity
  ParallelFor(zbx, [=] AMREX_GPU_DEVICE(int i, int j, int k) noexcept {
      z_vel_pert(i, j, k) = 0.0;
  });

  Gpu::streamSynchronize();
}

void
Problem::erf_init_rayleigh(
    amrex::Vector<amrex::Vector<amrex::Real> >& rayleigh_ptrs,
    amrex::Geometry      const& geom,
    std::unique_ptr<MultiFab>& /*z_phys_nd*/,
    amrex::Real /*zdamp*/)
{
  const int khi = geom.Domain().bigEnd()[2];

  // We just use these values to test the Rayleigh damping
  for (int k = 0; k <= khi; k++)
  {
      rayleigh_ptrs[Rayleigh::ubar][k]     = 2.0;
      rayleigh_ptrs[Rayleigh::vbar][k]     = 1.0;
      rayleigh_ptrs[Rayleigh::wbar][k]     = 0.0;
      rayleigh_ptrs[Rayleigh::thetabar][k] = parms.T_0;
  }
}
