#include <AMReX_MultiFab.H>
#include <ERF_SrcHeaders.H>
#include <AMReX_ParmParse.H>

using namespace amrex;

void
ApplySurfaceTreatment_BulkCoeff_Mom (
  const Box& tbx,
  const Box& tby,
  const Array4<Real>& rho_u_rhs,
  const Array4<Real>& rho_v_rhs,
  const Array4<const Real>& rho_u,
  const Array4<const Real>& rho_v,
  const Array4<const Real>& cons_state,
  const Array4<const Real>& z_phys_nd,
  const Array4<const Real>& surface_state_arr)
{
    int ndrag = 1;
    Real Cd_sea = Real(0.001);
    Real Cd_land = Real(0.01);

    ParallelFor(tbx, [=] AMREX_GPU_DEVICE(int i, int j, int k)
    {
        if(k <= ndrag) {
            Real dz = z_phys_nd(i,j,1) - z_phys_nd(i,j,0);
            Real ls_mask = (surface_state_arr(i-1,j,0)+surface_state_arr(i,j,0))/two;
            Real weight = one - Real(k)/ndrag;  // linear decrease
            Real fac = k==0? one : zero;
            Real Cd = fac*Cd_sea*(one-ls_mask) + Cd_land*ls_mask;
            Real rho_for_u = (cons_state(i-1,j,k,0)+cons_state(i,j,k,0))/two;
            Real rho_for_v = (cons_state(i,j-1,k,0)+cons_state(i,j,k,0))/two;
            Real uvel = rho_u(i,j,k)/rho_for_u;
            Real vvel = rho_v(i,j,k)/rho_for_v;
            Real velmag = std::sqrt(uvel*uvel + vvel*vvel);
            rho_u_rhs(i, j, k) += -one*weight*Cd*velmag*rho_for_u*uvel/dz;
        }
    });


    ParallelFor(tby, [=] AMREX_GPU_DEVICE(int i, int j, int k)
    {
       if(k <= ndrag) {
            Real dz = z_phys_nd(i,j,1) - z_phys_nd(i,j,0);
            Real ls_mask = (surface_state_arr(i,j-1,0)+surface_state_arr(i,j,0))/two;
            Real fac = k==0? one : zero;
            Real Cd = fac*Cd_sea*(one-ls_mask) + Cd_land*ls_mask;
            Real weight = one - Real(k)/ndrag;  // linear decrease
            Real rho_for_u = (cons_state(i-1,j,k,0)+cons_state(i,j,k,0))/two;
            Real rho_for_v = (cons_state(i,j-1,k,0)+cons_state(i,j,k,0))/two;
            Real uvel = rho_u(i,j,k)/rho_for_u;
            Real vvel = rho_v(i,j,k)/rho_for_v;
            Real velmag = std::sqrt(uvel*uvel + vvel*vvel);
            rho_v_rhs(i, j, k) += -one*weight*Cd*velmag*rho_for_v*vvel/dz;
        }
    });
}

void
ApplySurfaceTreatment_BulkCoeff_CC (const Box& bx,
                         const Array4<Real>& cell_rhs,
                         const Array4<const Real>& cons_state,
                         const Array4<const Real>& z_phys_cc,
                         const Array4<const Real>& surface_state_arr)
{
     ParallelFor(bx, [=] AMREX_GPU_DEVICE(int i, int j, int k) noexcept {
        if(k == 0) {
            Real dz = z_phys_cc(i,j,1)-z_phys_cc(i,j,0);
            Real ls_mask = surface_state_arr(i,j,0);
            Real Ch = Real(0.0015)*(one-ls_mask);
            Real Ce = Real(0.0015)*(one-ls_mask);

            Real rho = cons_state(i,j,k, Rho_comp);
            Real rhotheta = cons_state(i,j,k, RhoTheta_comp);
            Real rhoqv = cons_state(i,j,k, RhoQ1_comp);
            Real theta = rhotheta/rho;
            Real qv = rhoqv/rho;
            Real temp = getTgivenRandRTh(rho, rhotheta, qv);
            Real uvel = cons_state(i,j,k,1)/cons_state(i,j,k,0);
            Real vvel = cons_state(i,j,k,2)/cons_state(i,j,k,0);
            Real velmag = std::sqrt(uvel*uvel + vvel*vvel);
            Real dT = max(zero, Real(301.0) - temp);
            Real dq = max(zero, Real(0.024) - qv);
            cell_rhs(i, j, k, RhoTheta_comp) += (theta/(Real(1005.0)*temp))*rho*Ch*velmag*dT/dz;
            cell_rhs(i, j, k, RhoQ1_comp) += rho*Ce*velmag*dq/dz;
        }
    });
}
