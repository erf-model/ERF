#include <TerrainMetrics.H>
#include <AMReX_ParmParse.H>
#include <math.h>
#define PI 3.141592653589793238462643383279502884197

using namespace amrex;

//*****************************************************************************************
// Compute the terrain grid from BTF or STF model
//
// NOTE: Multilevel is not yet working for either of these terrain-following coordinates,
//       but (we think) the issue is deep in ERF and this code will work once the deeper
//       problem is fixed. For now, make sure to run on a single level. -mmsanders
//*****************************************************************************************
void 
init_terrain_grid( int lev, Geometry& geom, MultiFab& z_phys_nd)
{
  auto dx = geom.CellSizeArray();
  auto ProbHiArr = geom.ProbHiArray();
  auto ProbLoArr = geom.ProbLoArray();

  // z_nd is nodal in all directions
  const Box& domain = geom.Domain();
  int imin = domain.smallEnd(0); int imax = domain.bigEnd(0) + 1;
  int jmin = domain.smallEnd(1); int jmax = domain.bigEnd(1) + 1;
  int kmin = domain.smallEnd(2); int kmax = domain.bigEnd(2) + 1;

  // User-selected method from inputs file (BTF default)
  ParmParse pp("erf");

  // "custom" refers to terrain analytically specified, such as WitchOfAgnesi 
  std::string terrain_type = "custom"; 

  int nz = domain.length(2)+1;
  amrex::Vector<Real> z_levels_h; 
  z_levels_h.resize(nz);

  // This is the default for z_levels
  for (int k = 0; k < nz; k++)
  {
      z_levels_h[k] = k * dx[2];
  }

  // But we can read them in from the inputs file as an alternative
  int n_zlevels = pp.queryarr("terrain_z_levels", z_levels_h, 0, nz);
  if (n_zlevels > 0 && n_zlevels != nz)
      amrex::Abort("You must specify a z_level for every value of k");

   amrex::Gpu::DeviceVector<Real> z_levels_d;
  z_levels_d.resize(nz);
#ifdef AMREX_USE_GPU
    Gpu::htod_memcpy(z_levels_d.data(), z_levels_h.data(), sizeof(Real)*nz);
#else
    std::memcpy(z_levels_d.data(), z_levels_h.data(), sizeof(Real)*nz);
#endif

  int terrain_smoothing = 0;
  pp.query("terrain_smoothing", terrain_smoothing);

  switch(terrain_smoothing) {
    case 0: // BTF Method
    {
      //********************************************************************************
      // Populate z_nd in valid box and account for domain boundary cells (zero grad).
      // Top and bottom z boundary populated at end (fillboundary doesn't touch these)
      //********************************************************************************
      int k0 = 0;
      Real ztop = ProbHiArr[2];

      // WoA terms to compute ghost cell values
      Real a = 0.5;
      Real num = 8 * a * a * a;
      Real xcen = 0.5 * (ProbLoArr[0] + ProbHiArr[0]);
      Real ycen = 0.5 * (ProbLoArr[1] + ProbHiArr[1]);

      int ngrow = z_phys_nd.nGrow();

      for ( amrex::MFIter mfi(z_phys_nd, amrex::TilingIfNotGPU()); mfi.isValid(); ++mfi )
      {
          //use valid box to get domain bounds in all cases
          const Box& vbx = mfi.validbox();

          Array4<Real> const& z_arr = z_phys_nd.array(mfi);
          auto const& z_lev = z_levels_d.data();

          // Does the box border the domain
          int limin = vbx.smallEnd(0); int limax = vbx.bigEnd(0);
          int ljmin = vbx.smallEnd(1); int ljmax = vbx.bigEnd(1);
          int lkmax = vbx.bigEnd(2);
          int bxflag  = (limin == imin || limax == imax) ? 1 : 0;
          int byflag  = (ljmin == jmin || ljmax == jmax) ? 1 : 0;

          Box xybx;
          if (lev > 0) {
              xybx = mfi.growntilebox(ngrow);
              xybx.setRange(2,0,lkmax+ngrow);
          }
          else
              xybx = vbx;

          ParallelFor(xybx, [=] AMREX_GPU_DEVICE (int i, int j, int k) {

              // Location of nodes
              Real x = (i * dx[0] - xcen);
              Real y = (j * dx[1] - ycen);
              Real z =  z_lev[k]; 

              // WoA Hill for ghost points
              Real height = num / (x*x + 4*a*a);

              // YZ plane ghost cells
              if(i<limin) z_arr(i,j,k) = height;
              if(i>limax) z_arr(i,j,k) = height;

              // XZ plane ghost cells
              if(j<ljmin) z_arr(i,j,k) = height;
              if(j>ljmax) z_arr(i,j,k) = height;

              // Corner cells along z
              if(i<limin && j<ljmin) z_arr(i,j,k) = height;
              if(i>limax && j>ljmax) z_arr(i,j,k) = height;
              if(i<limin && j>ljmax) z_arr(i,j,k) = height;
              if(i>limax && j<ljmin) z_arr(i,j,k) = height;

              // Fill non-ghost cells with BTF model from p2163 of Klemp2011
              z_arr(i,j,k) = z + (1. - (z/ztop)) * z_arr(i,j,k0);

              // Fill boundary cells if needed
              if(bxflag || byflag) terrain_fill_domain_bndry_XY(i, j, k, limin, ljmin, limax, ljmax, z_arr(i,j,k), z_arr);

          });
        }

        break;
    }

    case 1: // STF Method
    {
        // Get Multifab spanning domain with 1 level of ghost cells
        amrex::MultiFab h_mf(z_phys_nd.boxArray(), z_phys_nd.DistributionMap(), 1, 1);
        amrex::MultiFab h_mf_old(z_phys_nd.boxArray(), z_phys_nd.DistributionMap(), 1, 1);

        // Save max height for smoothing
        Real max_h;

        // Crete 2D MF without allocation
        amrex::MultiFab mf2d;
        {
            BoxList bl2d = h_mf.boxArray().boxList();
            for (auto& b : bl2d) {
                b.setRange(2,0);
            }
            BoxArray ba2d(std::move(bl2d));
            mf2d = amrex::MultiFab(ba2d, h_mf.DistributionMap(), 1, 0, amrex::MFInfo().SetAlloc(false));
        }

        // Get MultiArray4s from the multifabs
        amrex::MultiArray4<Real> const& ma_h_s     = h_mf.arrays();
        amrex::MultiArray4<Real> const& ma_h_s_old = h_mf_old.arrays();
        amrex::MultiArray4<Real> const& ma_z_phys  = z_phys_nd.arrays();

        // Get max value
        max_h = ParReduce(amrex::TypeList<amrex::ReduceOpMax>{}, amrex::TypeList<Real>{}, mf2d, amrex::IntVect(0),
                    [=] AMREX_GPU_DEVICE (int box_no, int i, int j, int) noexcept
                        -> amrex::GpuTuple<Real>
                {

                  // Bottom boundary
                  int k0 = 0;

                  // Get Array4s
                  auto& h     = ma_h_s[box_no];
                  auto& z_arr = ma_z_phys[box_no];

                  // Does the box border the domain?
                  int limin = lbound(h).x + 1; int limax = ubound(h).x - 1;
                  int ljmin = lbound(h).y + 1; int ljmax = ubound(h).y - 1;
                  int bxflag  = (limin == imin || limax == imax) ? 1 : 0;
                  int byflag  = (ljmin == jmin || ljmax == jmax) ? 1 : 0;


                  // Populate h with terrain
                  h(i,j,k0) = z_arr(i,j,k0);


                  // Fill boundary cells if needed
                  if(bxflag || byflag){
                    terrain_fill_domain_bndry_XY(i, j, k0, imin, jmin, imax, jmax, z_arr(i,j,k0), z_arr);
                    terrain_fill_domain_bndry_XY(i, j, k0, imin, jmin, imax, jmax, z_arr(i,j,k0), h    );
                  }

                  // Return height for max
                  return { z_arr(i,j,k0) };
                });

        // Fill ghost cells (neglects domain boundary)
        h_mf.FillBoundary(geom.periodicity());

        // Make h_mf copy for old values
        amrex::MultiFab::Copy(h_mf_old, h_mf,0,0,1,1);

        // Populate h_mf at k>0 with h_s, solving in ordered 2D slices
        for (int k = kmin+1; k <= kmax; k++) // skip terrain level
        {
            Real zz = k * dx[2];
            Real zz_minus = (k - 1) * dx[2];

            Real gamma_m = 0.5; // min allowed fractional grid spacing
            Real z_H = 2.44/(1-gamma_m);

            Real foo = cos((PI/2)*(zz/z_H));
            Real A;
            if(zz < z_H) { A = foo*foo*foo*foo*foo*foo; } // A controls rate of return to atm
            else         { A = 0; }
            Real foo_minus = cos((PI/2)*(zz_minus/z_H));
            Real A_minus;
            if(zz_minus < z_H) { A_minus = foo_minus*foo_minus*foo_minus*foo_minus*foo_minus*foo_minus; } // A controls rate of return to atm
            else               { A_minus = 0; }

            unsigned maxIter = 50; // M_k in paper
            unsigned iter = 0;
            Real threshold = gamma_m;
            Real diff = 1.e20;
            while (iter < maxIter && diff > threshold)
            {
                diff = ParReduce(amrex::TypeList<amrex::ReduceOpMin>{}, amrex::TypeList<Real>{}, mf2d, amrex::IntVect(0),
                    [=] AMREX_GPU_DEVICE (int box_no, int i, int j, int) noexcept
                        -> amrex::GpuTuple<Real>
                {
                    auto& h_s     = ma_h_s[box_no];
                    auto& h_s_old = ma_h_s_old[box_no];

                    Real h_m    = max_h; //high point of hill
                    Real beta_k = 0.2*std::min(zz/(2*h_m),1.0); //smoothing coefficient

                    if (iter == 0) {
                        h_s(i,j,k) = h_s_old(i,j,k-1) + beta_k*(h_s_old(i+1,j  ,k-1)
                                                              + h_s_old(i-1,j  ,k-1)
                                                              + h_s_old(i  ,j+1,k-1)
                                                              + h_s_old(i  ,j-1,k-1) - 4*h_s_old(i,j,k-1));
                    }
                    else {
                        h_s(i,j,k) = h_s_old(i,j,k) + beta_k*(h_s_old(i+1,j  ,k)
                                                            + h_s_old(i-1,j  ,k)
                                                            + h_s_old(i  ,j+1,k)
                                                            + h_s_old(i  ,j-1,k) - 4*h_s_old(i,j,k));
                    }

                    return { (zz + A * h_s(i,j,k) - (zz_minus + A_minus * h_s(i,j,k-1))) / (zz - zz_minus) };

                }); //ParReduce



                ParallelFor(mf2d, [=] AMREX_GPU_DEVICE (int box_no, int i, int j, int)
                {
                    ma_h_s_old[box_no](i,j,k) = ma_h_s[box_no](i,j,k);

                    // Fill boundary cells
                    terrain_fill_domain_bndry_XY(i, j, k, imin, jmin, imax, jmax, ma_h_s[box_no](i,j,k), ma_h_s_old[box_no]);

                });

                amrex::ParallelDescriptor::ReduceRealMin(diff);

                iter++;

                //fill ghost points
                h_mf_old.FillBoundary(geom.periodicity());

            } //while

            //Populate z_phys_nd by solving z_arr(i,j,k) = z + A*h_s(i,j,k)
            for ( amrex::MFIter mfi(z_phys_nd, amrex::TilingIfNotGPU()); mfi.isValid(); ++mfi )
            {
                // Only loop the XY-plane of valid box
                const Box& vbx = mfi.validbox();
                Box xybx = vbx;
                xybx.setRange(2,0);
                Array4<Real> const& h_s   = h_mf_old.array(mfi);
                Array4<Real> const& z_arr = z_phys_nd.array(mfi);

                // Does the box border the domain?
                int limin = vbx.smallEnd(0); int limax = vbx.bigEnd(0);
                int ljmin = vbx.smallEnd(1); int ljmax = vbx.bigEnd(1);
                int bxflag  = (limin == imin || limax == imax) ? 1 : 0;
                int byflag  = (ljmin == jmin || ljmax == jmax) ? 1 : 0;

                ParallelFor(xybx, [=] AMREX_GPU_DEVICE (int i, int j, int) {

                    // Location of nodes
                    Real z = k * dx[2];

                    // STF model from p2164 of Klemp2011
                    z_arr(i,j,k) = z + A*h_s(i,j,k);

                    // Fill boundary cells if needed
                    if(bxflag || byflag) terrain_fill_domain_bndry_XY(i, j, k, imin, jmin, imax, jmax, z_arr(i,j,k), z_arr);

                });
            }//mfi
        } //k

        amrex::Gpu::streamSynchronize();

        break;
      }
  } //switch

  // Fill ghost layers and corners
  z_phys_nd.FillBoundary(geom.periodicity());

  //********************************************************************************
  // Populate domain boundary cells in z-direction
  //********************************************************************************
  for ( MFIter mfi(z_phys_nd, TilingIfNotGPU()); mfi.isValid(); ++mfi )
  {
       // Only loop the XY-plane of grown box
       const Box& vbx = mfi.growntilebox(1);
       Box xybx = vbx;
       xybx.setRange(2,0);
       Array4<Real> const& z_arr = z_phys_nd.array(mfi);

       ParallelFor(xybx, [=] AMREX_GPU_DEVICE (int i, int j, int ) {

           // Extrapolate top/bottom layer
           z_arr(i,j,kmin-1) = 2.0*z_arr(i,j,kmin) - z_arr(i,j,kmin+1);
           z_arr(i,j,kmax+1) = 2.0*z_arr(i,j,kmax) - z_arr(i,j,kmax-1);

           });
  }

  /*
  // Debug
  for ( MFIter mfi(z_phys_nd, TilingIfNotGPU()); mfi.isValid(); ++mfi )
  {
       const Box& gbx = mfi.growntilebox(1);
       amrex::Print() << "Ann print: " << z_phys_nd[mfi].box() << std::endl;
       Array4<Real> z_arr = z_phys_nd.array(mfi);

       int rank = 0;
       amrex::Print(rank) << gbx << "\n";

       ParallelFor(gbx, [=] AMREX_GPU_DEVICE (int i, int j, int k) {
           amrex::Print(rank) << amrex::IntVect(i,j,k) << "\n";
           amrex::Print(rank) << z_arr(i,j,k) << "\n";
           amrex::Print(rank) << "\n";
       });
   }
  */
}
