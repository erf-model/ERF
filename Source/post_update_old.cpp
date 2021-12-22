#include "ET_Integration.H"

using namespace amrex;

void rescale_state (MultiFab& state_mf)
{
#ifdef _OPENMP
#pragma omp parallel
#endif
  for ( MFIter mfi(state_mf, TilingIfNotGPU()); mfi.isValid(); ++mfi )
  {
    const Box& bx = mfi.tilebox();
    const auto ncomp = state_mf.nComp();

    const auto& state_fab = state_mf.array(mfi);

    // For each grid, loop over all the valid points
    amrex::ParallelFor(bx,
    [=] AMREX_GPU_DEVICE (int i, int j, int k) noexcept
    {
	// for example:
	amrex::Real det = state_fab(i, j, k, Idx::gambar00)*state_fab(i, j, k, Idx::gambar11)*state_fab(i, j, k, Idx::gambar22) - state_fab(i, j, k, Idx::gambar00)*std::pow(state_fab(i, j, k, Idx::gambar12), 2) - std::pow(state_fab(i, j, k, Idx::gambar01), 2)*state_fab(i, j, k, Idx::gambar22) + 2*state_fab(i, j, k, Idx::gambar01)*state_fab(i, j, k, Idx::gambar02)*state_fab(i, j, k, Idx::gambar12) - std::pow(state_fab(i, j, k, Idx::gambar02), 2)*state_fab(i, j, k, Idx::gambar11);
	amrex::Real scale_factor = 1.0/std::pow(det,1.0/3.0);
	
	state_fab(i, j, k, Idx::gambar00) = state_fab(i, j, k, Idx::gambar00) * scale_factor;
	state_fab(i, j, k, Idx::gambar01) = state_fab(i, j, k, Idx::gambar01) * scale_factor;
	state_fab(i, j, k, Idx::gambar02) = state_fab(i, j, k, Idx::gambar02) * scale_factor;
	state_fab(i, j, k, Idx::gambar11) = state_fab(i, j, k, Idx::gambar11) * scale_factor;
	state_fab(i, j, k, Idx::gambar12) = state_fab(i, j, k, Idx::gambar12) * scale_factor;
	state_fab(i, j, k, Idx::gambar22) = state_fab(i, j, k, Idx::gambar22) * scale_factor;

	amrex::Real gambar00 = state_fab(i, j, k, Idx::gambar00);
        amrex::Real gambar01 = state_fab(i, j, k, Idx::gambar01);
        amrex::Real gambar02 = state_fab(i, j, k, Idx::gambar02);
        amrex::Real gambar10 = state_fab(i, j, k, Idx::gambar01);
        amrex::Real gambar11 = state_fab(i, j, k, Idx::gambar11);
        amrex::Real gambar12 = state_fab(i, j, k, Idx::gambar12);
        amrex::Real gambar20 = state_fab(i, j, k, Idx::gambar02);
        amrex::Real gambar21 = state_fab(i, j, k, Idx::gambar12);
        amrex::Real gambar22 = state_fab(i, j, k, Idx::gambar22);	

	amrex::Real Abar00 = state_fab(i, j, k, Idx::Abar00);
        amrex::Real Abar01 = state_fab(i, j, k, Idx::Abar01);
        amrex::Real Abar02 = state_fab(i, j, k, Idx::Abar02);
        amrex::Real Abar10 = state_fab(i, j, k, Idx::Abar01);
        amrex::Real Abar11 = state_fab(i, j, k, Idx::Abar11);
        amrex::Real Abar12 = state_fab(i, j, k, Idx::Abar12);
        amrex::Real Abar20 = state_fab(i, j, k, Idx::Abar02);
        amrex::Real Abar21 = state_fab(i, j, k, Idx::Abar12);
        amrex::Real Abar22 = state_fab(i, j, k, Idx::Abar22);

	amrex::Real gambarinv00 = (gambar11*gambar22 - gambar12*gambar21)/(gambar00*gambar11*gambar22 - gambar00*gambar12*gambar21 - gambar01*gambar10*gambar22 + gambar01*gambar12*gambar20 + gambar02*gambar10*gambar21 - gambar02*gambar11*gambar20);
        amrex::Real gambarinv01 = (-gambar01*gambar22 + gambar02*gambar21)/(gambar00*gambar11*gambar22 - gambar00*gambar12*gambar21 - gambar01*gambar10*gambar22 + gambar01*gambar12*gambar20 + gambar02*gambar10*gambar21 - gambar02*gambar11*gambar20);
        amrex::Real gambarinv02 = (gambar01*gambar12 - gambar02*gambar11)/(gambar00*gambar11*gambar22 - gambar00*gambar12*gambar21 - gambar01*gambar10*gambar22 + gambar01*gambar12*gambar20 + gambar02*gambar10*gambar21 - gambar02*gambar11*gambar20);
        amrex::Real gambarinv10 = (-gambar10*gambar22 + gambar12*gambar20)/(gambar00*gambar11*gambar22 - gambar00*gambar12*gambar21 - gambar01*gambar10*gambar22 + gambar01*gambar12*gambar20 + gambar02*gambar10*gambar21 - gambar02*gambar11*gambar20);
        amrex::Real gambarinv11 = gambar00*(gambar00*gambar22 - gambar02*gambar20)/((gambar00*gambar11 - gambar01*gambar10)*(gambar00*gambar22 - gambar02*gambar20) - (gambar00*gambar12 - gambar02*gambar10)*(gambar00*gambar21 - gambar01*gambar20));
        amrex::Real gambarinv12 = -gambar00*(gambar00*gambar12 - gambar02*gambar10)/((gambar00*gambar11 - gambar01*gambar10)*(gambar00*gambar22 - gambar02*gambar20) - (gambar00*gambar12 - gambar02*gambar10)*(gambar00*gambar21 - gambar01*gambar20));
        amrex::Real gambarinv20 = (gambar10*gambar21 - gambar11*gambar20)/(gambar00*gambar11*gambar22 - gambar00*gambar12*gambar21 - gambar01*gambar10*gambar22 + gambar01*gambar12*gambar20 + gambar02*gambar10*gambar21 - gambar02*gambar11*gambar20);
        amrex::Real gambarinv21 = -gambar00*(gambar00*gambar21 - gambar01*gambar20)/((gambar00*gambar11 - gambar01*gambar10)*(gambar00*gambar22 - gambar02*gambar20) - (gambar00*gambar12 - gambar02*gambar10)*(gambar00*gambar21 - gambar01*gambar20));
        amrex::Real gambarinv22 = gambar00*(gambar00*gambar11 - gambar01*gambar10)/((gambar00*gambar11 - gambar01*gambar10)*(gambar00*gambar22 - gambar02*gambar20) - (gambar00*gambar12 - gambar02*gambar10)*(gambar00*gambar21 - gambar01*gambar20));
    
	amrex::Real TrAbar = Abar00*gambarinv00 + Abar01*gambarinv01 + Abar02*gambarinv02 + Abar10*gambarinv10 + Abar11*gambarinv11 + Abar12*gambarinv12 + Abar20*gambarinv20 + Abar21*gambarinv21 + Abar22*gambarinv22;

	state_fab(i, j, k, Idx::Abar00) = Abar00 - 1.0/3.0*gambar00*TrAbar;
	state_fab(i, j, k, Idx::Abar01) = Abar01 - 1.0/3.0*gambar01*TrAbar;
	state_fab(i, j, k, Idx::Abar02) = Abar02 - 1.0/3.0*gambar02*TrAbar;
	state_fab(i, j, k, Idx::Abar11) = Abar11 - 1.0/3.0*gambar11*TrAbar;
	state_fab(i, j, k, Idx::Abar12) = Abar12 - 1.0/3.0*gambar12*TrAbar;
	state_fab(i, j, k, Idx::Abar22) = Abar22 - 1.0/3.0*gambar22*TrAbar;

	});
  }
}
