#ifndef _PBLMODELS_H_
#define _PBLMODELS_H_

#include <TimeIntegration.H>
#include "DirectionSelector.H"

inline
void ComputeTurbulentViscosityPBL(const MultiFab& xvel, const MultiFab& yvel, const MultiFab& zvel,
                                  const MultiFab& cons_in, MultiFab& eddyViscosity,
                                  const Geometry& geom,
                                  const SolverChoice& solverChoice/*,
								    const ABLMost& most*/)
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

      // add to datastruct -> Require MOST with MYNN
    
    for ( MFIter mfi(eddyViscosity,TilingIfNotGPU()); mfi.isValid(); ++mfi) {
      
      const Box &bx = mfi.growntilebox(1);
      const Array4<Real const > &cell_data = cons_in.array(mfi);
      const Array4<Real> &K_turb = eddyViscosity.array(mfi);
      const Array4<Real const> &uvel = xvel.array(mfi);
      const Array4<Real const> &vvel = yvel.array(mfi);
      const Array4<Real const> &wvel = zvel.array(mfi);
      
      // Compute some quantities that are constant in each column
      // Sbox is shrunk to only include the interior of the domain in the vertical direction to compute integrals
      // NOTE: Here we requite that sbx covers the entire vertical domain
      const Box &dbx = geom.Domain();
      Box sbx(bx.smallEnd(), bx.bigEnd());
      sbx.grow(2,-1); 
      AMREX_ALWAYS_ASSERT(sbx.smallEnd(2) == dbx.smallEnd(2) && sbx.bigEnd(2) == dbx.bigEnd(2));

      const GeometryData gdata = geom.data();
      
      const Box xybx = PerpendicularBox<ZDir>(bx, IntVect{0,0,0});
      FArrayBox qintegral(xybx,2);
      qintegral.setVal(0.0);
      FArrayBox qturb(bx,1);
      const Array4<amrex::Real> qint = qintegral.array();
      const Array4<amrex::Real> qvel= qturb.array();
      
      amrex::ParallelFor(Gpu::KernelInfo().setReduction(true),
			 bx, [=] AMREX_GPU_DEVICE (int i, int j, int k, Gpu::Handler const& handler) noexcept
      {
	const amrex::Real Zval = gdata.ProbLo(2) + (k + 0.5)*gdata.CellSize(2);
	const amrex::Real rho = cell_data(i,j,k,Rho_comp);
	qvel(i,j,k) = std::sqrt(cell_data(i,j,k,RhoQKE_comp) / cell_data(i,j,k,Rho_comp));
	if (sbx.contains(i,j,k)) {
	  Gpu::deviceReduceSum(&qint(i,j,0,0), Zval*qvel(i,j,k), handler);
	  Gpu::deviceReduceSum(&qint(i,j,0,1), qvel(i,j,k), handler);
	}
      });

      amrex::Real dz_inv = geom.InvCellSize(2);
      int izmin = geom.Domain().smallEnd(2);
      int izmax = geom.Domain().bigEnd(2);
      const amrex::Real l_obukhov = 4;//most.obukhov_len;
      const amrex::Real surface_heat_flux = 7;//most.surf_temp_flux;
      const amrex::Real theta0 = 300.0;//most.theta_mean (TODO: IS THIS ACTUALLY RHOTHETA)
      amrex::ParallelFor(bx, [=] AMREX_GPU_DEVICE (int i, int j, int k) noexcept
      {
	// Compute some partial derivatives that we will need (1st order at domain boundary)
	amrex::Real dthetadz, dudz, dvdz;
	if (k == izmax) {
	  dthetadz = (cell_data(i,j,k,RhoTheta_comp)/cell_data(i,j,k,Rho_comp) -
		      cell_data(i,j,k-1,RhoTheta_comp)/cell_data(i,j,k-1,Rho_comp))*dz_inv;
	  dudz = (uvel(i,j,k) - uvel(i,j,k-1))*dz_inv;
	  dvdz = (uvel(i,j,k) - uvel(i,j,k-1))*dz_inv;
	} else if (k == izmin){
	  dthetadz = (cell_data(i,j,k+1,RhoTheta_comp)/cell_data(i,j,k+1,Rho_comp) -
		      cell_data(i,j,k,RhoTheta_comp)/cell_data(i,j,k,Rho_comp))*dz_inv;
	  dudz = (uvel(i,j,k+1) - uvel(i,j,k))*dz_inv;
	  dvdz = (uvel(i,j,k+1) - uvel(i,j,k))*dz_inv;
	} else {
	  dthetadz = 0.5*(cell_data(i,j,k+1,RhoTheta_comp)/cell_data(i,j,k+1,Rho_comp) -
			  cell_data(i,j,k-1,RhoTheta_comp)/cell_data(i,j,k-1,Rho_comp))*dz_inv;
	  dudz = 0.5*(uvel(i,j,k+1) - uvel(i,j,k-1))*dz_inv;
	  dvdz = 0.5*(uvel(i,j,k+1) - uvel(i,j,k-1))*dz_inv;
	}
	  
	// First Length Scale
	AMREX_ASSERT(l_obukhov != 0);
	const amrex::Real zval = gdata.ProbLo(2) + (k + 0.5)*gdata.CellSize(2);
	const amrex::Real zeta = zval/l_obukhov;
	amrex::Real l_S;
	if (zeta >= 1.0) {
	  l_S = KAPPA*zval/3.7;
	} else if (zeta >= 0) {
	  l_S = KAPPA*zval/(1+2.7*zeta);
	} else {
	  l_S = KAPPA*zval*std::pow(1.0 - 100.0 * zeta, 0.2);
	}

	// Second Length Scale
	amrex::Real l_T;
	if (qint(i,j,0,1) > 0.0) {
	  l_T = 0.23*qint(i,j,0,0)/qint(i,j,0,1);
	} else {
	  l_T = std::numeric_limits<amrex::Real>::max();
	}
	
	// Third Length Scale
	amrex::Real l_B;
	if (dthetadz > 0) {
	  amrex::Real N_brunt_vaisala = CONST_GRAV/theta0 * std::sqrt(dthetadz);
	  if (zeta < 0) {
	    amrex::Real qc = CONST_GRAV/theta0 * surface_heat_flux * l_T;
	    qc = std::pow(qc,1.0/3.0);
	    l_B = (1.0 + 5.0*std::sqrt(qc/(N_brunt_vaisala * l_T))) * qvel(i,j,k)/N_brunt_vaisala;
	  } else {
	    l_B = qvel(i,j,k) / N_brunt_vaisala;
	  }
	} else {
	  l_B = std::numeric_limits<amrex::Real>::max();
	}

	// Overall Length Scale
	amrex::Real l_comb = 1.0 / (1.0/l_S + 1.0/l_T + 1.0/l_B);
	
	// Compute non-dimensional parameters
	amrex::Real l2_over_q2 = l_comb*l_comb/(qvel(i,j,k)*qvel(i,j,k))
	amrex::Real GM = l2_over_q2 * (dudz*dudz + dvdz*dvdz);
	amrex::Real GH = -l2_over_q2 / theta0 * dthetadz;
	
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
    }
  }
}

#endif
