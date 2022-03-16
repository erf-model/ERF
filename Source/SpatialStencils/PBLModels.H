#ifndef _PBLMODELS_H_
#define _PBLMODELS_H_

#include <TimeIntegration.H>
#include "DirectionSelector.H"

inline
void ComputeTurbulentViscosityPBL(const MultiFab& xvel, const MultiFab& yvel, const MultiFab& zvel,
                                  const MultiFab& cons_in, MultiFab& eddyViscosity,
                                  const Geometry& geom,
                                  const SolverChoice& solverChoice,
                                  const Gpu::DeviceVector<amrex::BCRec> domain_bcs_type_d)
{
  // MYNN Level 2.5 PBL Model
  if (solverChoice.pbl_type == PBLType::MYNN25) {

    const amrex::Real A1 = solverChoice.pbl_A1;
    const amrex::Real A2 = solverChoice.pbl_A2;
    const amrex::Real B1 = solverChoice.pbl_B1;
    const amrex::Real B2 = solverChoice.pbl_B2;
    const amrex::Real C1 = solverChoice.pbl_C1;
    const amrex::Real C2 = solverChoice.pbl_C2;
    const amrex::Real C3 = solverChoice.pbl_C3;
    const amrex::Real C4 = solverChoice.pbl_C4;
    const amrex::Real C5 = solverChoice.pbl_C5;

      // Quantites that are needed as a function of (i,j) constant over whole k column (and need whole k to calculate)
      //    L_obukhov (from MOST)
      //    l_T (integral over column)
      //    vertical surface heat flux \bar{w*theta}_0 (from MOST)
      // add to datastruct -> Require MOST with MYNN
    
    for ( MFIter mfi(eddyViscosity,TilingIfNotGPU()); mfi.isValid(); ++mfi) {
      
      const Box &bx = mfi.growntilebox(1);
      const Array4<Real const > &cell_data = cons_in.array(mfi);
      const Array4<Real> &K_turb = eddyViscosity.array(mfi);
      const Array4<Real const> &uvel = xvel.array(mfi);
      const Array4<Real const> &vvel = yvel.array(mfi);
      const Array4<Real const> &wvel = zvel.array(mfi);
      const amrex::BCRec* bc_ptr = domain_bcs_type_d.data();
      
      // Compute some quantities that are constant in each column
      // Sbox is shrunk to only include the interior of the domain in the vertical direction to compute integrals
      // NOTE: Here we requite that sbx covers the entire vertical domain
      const Box &dbx = geom.Domain();
      Box sbx(bx.smallEnd(), bx.bigEnd());
      sbx.grow(2,-1); 
      AMREX_ALWAYS_ASSERT(sbx.smallEnd(2) == dbx.smallEnd(2) && sbx.bigEnd(2) == dbx.bigEnd(2));
      
      const Box xybx = PerpendicularBox<ZDir>(bx, IntVect{0,0,0});
      FArrayBox qintegral(xybx,2);
      qintegral.setVal(0.0);
      FArrayBox qturb(bx,1);
      const Array4<amrex::Real> qint = qintegral.array();
      const Array4<amrex::Real> qvel= qturb.array();
      
      amrex::ParallelFor(Gpu::KernelInfo().setReduction(true),
			 bx, [=] AMREX_GPU_DEVICE (int i, int j, int k, Gpu::Handler const& handler) noexcept
      {
	amrex::Real Zval = geom.ProbLo(2) + (k + 0.5)*geom.CellSize(2);
	amrex::Real rho = cell_data(i,j,k,Rho_comp);
	qvel(i,j,k) = std::sqrt(cell_data(i,j,k,RhoQKE_comp) / cell_data(i,j,k,Rho_comp));
	if (sbx.contains(i,j,k)) {
	  Gpu::deviceReduceSum(&qint(i,j,0,0), Zval*qvel(i,j,k), handler);
	  Gpu::deviceReduceSum(&qint(i,j,0,1), qvel(i,j,k), handler);
	}
      });
      
      amrex::ParallelFor(bx, [=] AMREX_GPU_DEVICE (int i, int j, int k) noexcept
      {
	amrex::Real l_T = qint(i,j,0,0)/qint(i,j,0,1);
	K_turb(i,j,k,EddyDiff::Mom_v) = cell_data(i, j, k, Rho_comp) * l_T * qvel(i,j,k);
      });
    }
	/*
	// Get some relevant cell quantities
	amrex::Real rho  = cell_data(i, j, k, Rho_comp);
	amrex::Real rhoinv = 1.0/rho;
	amrex::Real q2   = cell_data(i, j, k, RhoQKE_comp) * rhoinv;
	amrex::Real qvel = std::sqrt(q2);	
	
	// Compute Length Scales
	amrex::Real L_obukhov = 0.0;
	amrex::Real l_s = 0.0;
	amrex::Real l_B = 0.0;
	amrex::Real l_T = 0.0;

	amrex::Real l_comb = 1.0 / (1.0/l_s + 1.0/l_T + 1.0/l_B);

	// Compute non-dimensional parameters
	amrex::Real GM = 0.0;
	amrex::Real GH = 0.0;
	
	amrex::Real E1 = 1.0 + 6.0*A1*A1*GM - 9.0*A1*A2*(1.0-C2)*GH;
	amrex::Real E2 = -3.0*A1*(4.0*A1 + 3.0*A2*(1.0-C5))*(1.0-C2)*GH;
	amrex::Real E3 = 6.0*A2*A1*GM;
	amrex::Real E4 = 1.0 - 12.0*A2*A1*(1.0-C2)*GH -3.0*A2*B2*(1.0-C3)*GH;
	amrex::Real R1 = A1*(1.0-3.0*C1);
	
	amrex::Real SM = (A2*E2 - R1*E4)/(E2*E3 - E1*E4);
	amrex::Real SH = (R1*E3 - A2*E1)/(E2*E3 - E1*E4);
	amrex::Real SQ = 3.0 * SM;

	// Finally, compute the eddy viscosity/diffusivities
	// TODO: Factor of 2 for momentum?
	K_turb(i,j,k,EddyDiff::Mom_v)   = rho * l_comb * qvel * SM;
	K_turb(i,j,k,EddyDiff::Theta_v) = rho * l_comb * qvel * SH;
	K_turb(i,j,k,EddyDiff::QKE_v)   = rho * l_comb * qvel * 3.0 * SQ;

	// TODO: How should this be done for other components (scalars, moisture)
	// Horizontal/LES component for QKE?  
      });
			     
    }*/
  }
}

#endif
