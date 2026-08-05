#ifndef ERF_RESOLVED_WALL_FLUX_H_
#define ERF_RESOLVED_WALL_FLUX_H_

#include <AMReX_Box.H>
#include <AMReX_GpuContainers.H>
#include <AMReX_MultiFab.H>

#include "ERF_Constants.H"
#include "ERF_EOS.H"
#include "ERF_MicrophysicsUtils.H"
#include "ERF_IndexDefines.H"
#include "ERF_DataStruct.H"
#include "BoundaryConditions/ERF_WallThermodynamics.H"

namespace erf_resolved_wall_flux {

template <int DIR, bool HIGH>
AMREX_FORCE_INLINE
void apply_face (const amrex::Box& bx,
                 const amrex::Box& domain,
                 const int quantity,
                 const int flux_comp,
                 const amrex::Array4<const amrex::Real>& state,
                 const amrex::Array4<const amrex::Real>& prim,
                 const amrex::Array4<const amrex::Real>& base_state,
                 const amrex::Array4<amrex::Real>& rhs,
                 const amrex::Array4<amrex::Real>& xflux,
                 const amrex::Array4<amrex::Real>& yflux,
                 const amrex::Array4<amrex::Real>& zflux,
                 const amrex::GpuArray<amrex::Real, AMREX_SPACEDIM>& dx_inv,
                 const erf_wall_thermodynamics::Boundary& walls,
                 const amrex::Real alpha_T,
                 const amrex::Real alpha_C,
                 const amrex::Real rdOcp)
{
    // A tile may belong to an interior BoxArray box.  Do not relocate that
    // tile to a global boundary: its Array4 views do not own data there.
    if constexpr (HIGH) {
        if (bx.bigEnd(DIR) != domain.bigEnd(DIR)) { return; }
    } else {
        if (bx.smallEnd(DIR) != domain.smallEnd(DIR)) { return; }
    }

    const int cell_face = HIGH ? domain.bigEnd(DIR) : domain.smallEnd(DIR);
    amrex::Box face_cells = bx;
    face_cells.setSmall(DIR, cell_face);
    face_cells.setBig(DIR, cell_face);
    face_cells &= domain;
    if (face_cells.isEmpty()) { return; }
    AMREX_ALWAYS_ASSERT(bx.contains(face_cells));

    const auto wall = walls[2*DIR + (HIGH ? 1 : 0)];
    const bool is_theta = quantity == RhoTheta_comp;
    const bool is_qv = quantity == RhoQ1_comp;
    const bool is_qc = quantity == RhoQ2_comp;
    if (!((is_theta && wall.thermal_mode == erf_wall_thermodynamics::ThermalMode::FixedPhysicalTemperature) ||
          (is_qv && wall.moisture_mode != erf_wall_thermodynamics::MoistureMode::LegacyNumerical) ||
          (is_qc && wall.moisture_mode != erf_wall_thermodynamics::MoistureMode::LegacyNumerical))) {
        return;
    }

    const auto& flux = (DIR == 0) ? xflux : ((DIR == 1) ? yflux : zflux);
    const amrex::Real alpha = is_theta ? alpha_T : alpha_C;
    const amrex::Real inv_half_cell = amrex::Real(2.0) * dx_inv[DIR];

    amrex::ParallelFor(face_cells, [=] AMREX_GPU_DEVICE (int i, int j, int k) noexcept
    {
        const int fi = HIGH && DIR == 0 ? i + 1 : i;
        const int fj = HIGH && DIR == 1 ? j + 1 : j;
        const int fk = HIGH && DIR == 2 ? k + 1 : k;
        const amrex::Real rho = state(i,j,k,Rho_comp);
        const amrex::Real p0 = base_state(i,j,k,BaseState::p0_comp);
        const amrex::Real wall_temperature = wall.physical_temperature_K;
        amrex::Real wall_value = amrex::Real(0.0);
        amrex::Real cell_value = prim(i,j,k,quantity-1);

        if (is_theta) {
            wall_value = wall_temperature * std::pow(p_0 / p0, rdOcp);
        } else if (is_qv && wall.moisture_mode == erf_wall_thermodynamics::MoistureMode::WetEquilibrium) {
            amrex::Real qsat = amrex::Real(0.0);
            erf_qsatw(wall_temperature, p0 * amrex::Real(0.01), qsat);
            wall_value = qsat;
        } else if (is_qc || (is_qv && wall.moisture_mode == erf_wall_thermodynamics::MoistureMode::DryImpermeable)) {
            wall_value = amrex::Real(0.0);
        }

        const amrex::Real old_flux = flux(fi,fj,fk,flux_comp);
        amrex::Real new_flux = amrex::Real(0.0);
        if (!(is_qc || (is_qv && wall.moisture_mode == erf_wall_thermodynamics::MoistureMode::DryImpermeable))) {
            amrex::ignore_unused(alpha, inv_half_cell);
            if constexpr (HIGH) {
                new_flux = -rho * alpha * (wall_value - cell_value) * inv_half_cell;
            } else {
                new_flux = -rho * alpha * (cell_value - wall_value) * inv_half_cell;
            }
        }
        flux(fi,fj,fk,flux_comp) = new_flux;

        // DiffusionSrcForState_N has already applied the old boundary flux to
        // the RHS.  Add only the boundary correction so the state update and
        // the retained flux array use exactly the same physical flux.
        const amrex::Real correction = (new_flux - old_flux) * dx_inv[DIR];
        amrex::ignore_unused(rhs);
        if constexpr (HIGH) {
            rhs(i,j,k,quantity) -= correction;
        } else {
            rhs(i,j,k,quantity) += correction;
        }
    });
}

inline void apply (const amrex::Box& bx,
                   const amrex::Box& domain,
                   const int quantity,
                   const int flux_comp,
                   const amrex::Array4<const amrex::Real>& state,
                   const amrex::Array4<const amrex::Real>& prim,
                   const amrex::Array4<const amrex::Real>& base_state,
                   const amrex::Array4<amrex::Real>& rhs,
                   const amrex::Array4<amrex::Real>& xflux,
                   const amrex::Array4<amrex::Real>& yflux,
                   const amrex::Array4<amrex::Real>& zflux,
                   const amrex::GpuArray<amrex::Real, AMREX_SPACEDIM>& dx_inv,
                   const erf_wall_thermodynamics::Boundary& walls,
                   const amrex::Real alpha_T,
                   const amrex::Real alpha_C,
                   const amrex::Real rdOcp)
{
    apply_face<0,false>(bx,domain,quantity,flux_comp,state,prim,base_state,rhs,xflux,yflux,zflux,dx_inv,walls,alpha_T,alpha_C,rdOcp);
    apply_face<0,true >(bx,domain,quantity,flux_comp,state,prim,base_state,rhs,xflux,yflux,zflux,dx_inv,walls,alpha_T,alpha_C,rdOcp);
    apply_face<1,false>(bx,domain,quantity,flux_comp,state,prim,base_state,rhs,xflux,yflux,zflux,dx_inv,walls,alpha_T,alpha_C,rdOcp);
    apply_face<1,true >(bx,domain,quantity,flux_comp,state,prim,base_state,rhs,xflux,yflux,zflux,dx_inv,walls,alpha_T,alpha_C,rdOcp);
    apply_face<2,false>(bx,domain,quantity,flux_comp,state,prim,base_state,rhs,xflux,yflux,zflux,dx_inv,walls,alpha_T,alpha_C,rdOcp);
    apply_face<2,true >(bx,domain,quantity,flux_comp,state,prim,base_state,rhs,xflux,yflux,zflux,dx_inv,walls,alpha_T,alpha_C,rdOcp);
}

} // namespace erf_resolved_wall_flux

#endif
