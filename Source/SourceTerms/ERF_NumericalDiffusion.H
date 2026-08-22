#ifndef ERF_NUMERICALDIFFUSION_H_
#define ERF_NUMERICALDIFFUSION_H_

#include <AMReX.H>
#include <ERF_DataStruct.H>
#include <AMReX_MultiFab.H>

/**
 * Compute a fifth-order derivative approximation.
 *
 * @param[in] dnp2 Value at index n+2
 * @param[in] dnp1 Value at index n+1
 * @param[in] dn Value at index n
 * @param[in] dnm1 Value at index n-1
 * @param[in] dnm2 Value at index n-2
 * @param[in] dnm3 Value at index n-3
 * @return The calculated derivative
 */
AMREX_GPU_DEVICE
AMREX_FORCE_INLINE
amrex::Real
calc_fifth_order_deriv (const amrex::Real& dnp2,
                        const amrex::Real& dnp1,
                        const amrex::Real& dn,
                        const amrex::Real& dnm1,
                        const amrex::Real& dnm2,
                        const amrex::Real& dnm3)
{
    amrex::Real interp = amrex::Real(10.) * (dn   - dnm1)
                        - amrex::Real(5.) * (dnp1 - dnm2)
                             + (dnp2 - dnm3);
    return interp;
}

/**
 * Apply numerical diffusion to scalar variables.
 *
 * @param[in] bx Box of cells to process
 * @param[in] start_comp Starting component index
 * @param[in] num_comp Number of components to process
 * @param[in] dt Time step
 * @param[in] num_diff_coeff Numerical diffusion coefficient
 * @param[in] prim_data Primitive variable data
 * @param[in] cell_data Cell-centered data
 * @param[out] rhs Right-hand side array for updates
 * @param[in] mfx Map factors in x
 * @param[in] mfy Map factors in y
 */
void NumericalDiffusion_Scal (const amrex::Box& bx,
                              const int start_comp,
                              const int num_comp,
                              const double dt,
                              const amrex::Real num_diff_coeff,
                              const amrex::Array4<const amrex::Real>& prim_data,
                              const amrex::Array4<const amrex::Real>& cell_data,
                              const amrex::Array4<      amrex::Real>& rhs,
                              const amrex::Array4<const amrex::Real>& mfx,
                              const amrex::Array4<const amrex::Real>& mfy);

/**
 * Apply numerical diffusion to the x-momentum equation.
 *
 * @param[in] bx Box of cells to process
 * @param[in] dt Time step
 * @param[in] num_diff_coeff Numerical diffusion coefficient
 * @param[in] prim_data Primitive variable data
 * @param[in] cell_data Cell-centered data
 * @param[out] rhs Right-hand side array for updates
 * @param[in] mfx Map factors in x
 * @param[in] mfy Map factors in y
 */
void NumericalDiffusion_Xmom (const amrex::Box& bx,
                              const double dt,
                              const amrex::Real num_diff_coeff,
                              const amrex::Array4<const amrex::Real>& prim_data,
                              const amrex::Array4<const amrex::Real>& cell_data,
                              const amrex::Array4<      amrex::Real>& rhs,
                              const amrex::Array4<const amrex::Real>& mfx,
                              const amrex::Array4<const amrex::Real>& mfy);

/**
 * Apply numerical diffusion to the y-momentum equation.
 *
 * @param[in] bx Box of cells to process
 * @param[in] dt Time step
 * @param[in] num_diff_coeff Numerical diffusion coefficient
 * @param[in] prim_data Primitive variable data
 * @param[in] cell_data Cell-centered data
 * @param[out] rhs Right-hand side array for updates
 * @param[in] mfx Map factors in x
 * @param[in] mfy Map factors in y
 */
void NumericalDiffusion_Ymom (const amrex::Box& bx,
                              const double dt,
                              const amrex::Real num_diff_coeff,
                              const amrex::Array4<const amrex::Real>& prim_data,
                              const amrex::Array4<const amrex::Real>& cell_data,
                              const amrex::Array4<      amrex::Real>& rhs,
                              const amrex::Array4<const amrex::Real>& mfx,
                              const amrex::Array4<const amrex::Real>& mfy);
#endif
