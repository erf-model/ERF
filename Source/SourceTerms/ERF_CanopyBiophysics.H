#ifndef ERF_CANOPY_BIOPHYSICS_H_
#define ERF_CANOPY_BIOPHYSICS_H_

#include <cmath>

#include <AMReX_MultiFab.H>

#include <ERF_Constants.H>
#include <ERF_DataStruct.H>
#include <ERF_IndexDefines.H>

/**
 * Add the portion of the ICLASS canopy-biophysics model exercised by the
 * fixed-leaf-potential-temperature configuration.  The source uses the leaf
 * area density produced by ForestDrag and adds sensible heat to rho-theta.
 * When moisture is active, it also adds the corresponding transpiration flux
 * using the reference model's default minimum stomatal conductance.
 * @param[in,out] cell_source staged source terms for the conserved state
 * @param[in] S_data conserved atmospheric state
 * @param[in] xvel x-face velocity
 * @param[in] yvel y-face velocity
 * @param[in] zvel z-face velocity
 * @param[in] frontal_area cell-centered leaf/frontal area density
 * @param[in] base_state hydrostatic pressure and Exner state
 * @param[in] solver_choice active physics and canopy configuration
 */
AMREX_INLINE void
AddCanopyBiophysicsHeatSources (amrex::MultiFab&       cell_source,
                                const amrex::MultiFab& S_data,
                                const amrex::MultiFab& xvel,
                                const amrex::MultiFab& yvel,
                                const amrex::MultiFab& zvel,
                                const amrex::MultiFab* frontal_area,
                                const amrex::MultiFab& base_state,
                                const SolverChoice&    solver_choice)
{
    BL_PROFILE("AddCanopyBiophysicsHeatSources");

    if (frontal_area == nullptr || !solver_choice.forest_biophysics_heat) {
        return;
    }

    const bool has_moisture = solver_choice.moisture_type != MoistureType::None;
    const amrex::Real leaf_theta = solver_choice.forest_leaf_theta_fixed;

    // Defaults used by the ICLASS canopy implementation for the fixed-theta
    // path.  Radiation and the bisection solve are not needed in this mode.
    constexpr amrex::Real leaf_width = amrex::Real(0.05); // m
    constexpr amrex::Real min_wind = amrex::Real(0.01);   // m/s
    constexpr amrex::Real boundary_coeff = amrex::Real(100.0);
    constexpr amrex::Real min_stomatal_conductance = amrex::Real(0.01); // mol/(m^2 s)
    constexpr amrex::Real gas_constant = amrex::Real(8.314);            // J/(mol K)

    // Make namespace-scope constants explicit lambda captures for CUDA.  The
    // host compiler can resolve these constexpr variables directly, whereas
    // NVCC does not expose them as device symbols from an implicit [=] capture.
    const amrex::Real zero_d   = zero;
    const amrex::Real one_d    = one;
    const amrex::Real two_d    = two;
    const amrex::Real myhalf_d = myhalf;
    const amrex::Real cp_d     = Cp_d;
    const amrex::Real rdo_rv_d = RdoRv;
    const amrex::Real lv_d     = L_v;

#ifdef _OPENMP
#pragma omp parallel if (amrex::Gpu::notInLaunchRegion())
#endif
    for (amrex::MFIter mfi(cell_source, amrex::TilingIfNotGPU()); mfi.isValid(); ++mfi) {
        const amrex::Box& box = mfi.tilebox();

        const auto& source = cell_source.array(mfi);
        const auto& state = S_data.const_array(mfi);
        const auto& u = xvel.const_array(mfi);
        const auto& v = yvel.const_array(mfi);
        const auto& w = zvel.const_array(mfi);
        const auto& leaf_area_density = frontal_area->const_array(mfi);
        const auto& base = base_state.const_array(mfi);

        amrex::ParallelFor(box, [=] AMREX_GPU_DEVICE (int i, int j, int k) noexcept
        {
            const amrex::Real lad = leaf_area_density(i, j, k);
            if (lad <= zero_d) {
                return;
            }

            const amrex::Real rho = state(i, j, k, Rho_comp);
            const amrex::Real theta_air = state(i, j, k, RhoTheta_comp) / rho;
            const amrex::Real exner = base(i, j, k, BaseState::pi0_comp);
            const amrex::Real pressure = base(i, j, k, BaseState::p0_comp);
            const amrex::Real temp_air = theta_air * exner;
            const amrex::Real temp_leaf = leaf_theta * exner;

            const amrex::Real u_cc = myhalf_d * (u(i, j, k) + u(i+1, j, k));
            const amrex::Real v_cc = myhalf_d * (v(i, j, k) + v(i, j+1, k));
            const amrex::Real w_cc = myhalf_d * (w(i, j, k) + w(i, j, k+1));
            const amrex::Real wind = std::sqrt(u_cc*u_cc + v_cc*v_cc + w_cc*w_cc);
            const amrex::Real resistance_boundary =
                boundary_coeff * std::sqrt(leaf_width / amrex::max(wind, min_wind));

            const amrex::Real sensible_heat =
                two_d * rho * cp_d * (temp_leaf - temp_air) / resistance_boundary;
            source(i, j, k, RhoTheta_comp) +=
                sensible_heat * lad / (cp_d * exner);

            if (has_moisture) {
                const amrex::Real qv = state(i, j, k, RhoQ1_comp) / rho;
                const amrex::Real temp_celsius = temp_leaf - amrex::Real(273.15);
                const amrex::Real saturation_pressure =
                    amrex::Real(611.2) * std::exp(amrex::Real(17.67) * temp_celsius /
                                                  (temp_celsius + amrex::Real(243.5)));
                const amrex::Real qv_saturation =
                    rdo_rv_d * saturation_pressure /
                    amrex::max(pressure - (one_d - rdo_rv_d) * saturation_pressure, one_d);
                const amrex::Real resistance_stomatal =
                    pressure /
                    amrex::max(min_stomatal_conductance * gas_constant * temp_air,
                               amrex::Real(1.0e-20));
                const amrex::Real latent_flux =
                    rho * lv_d * amrex::max(qv_saturation - qv, zero_d) /
                    (amrex::Real(1.075) * resistance_boundary + resistance_stomatal);

                source(i, j, k, RhoQ1_comp) += latent_flux * lad / lv_d;
            }
        });
    }
}

#endif
