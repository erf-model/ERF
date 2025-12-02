#ifndef ERF_SOLVE_TRIDIAG
#define ERF_SOLVE_TRIDIAG

AMREX_GPU_HOST_DEVICE
AMREX_FORCE_INLINE
void SolveTridiag (int i, int j, int klo, int khi,
                   const amrex::Array4<      amrex::Real>& soln_a,
                   const amrex::Array4<const amrex::Real>& coeffA_a, // lower diagonal
                   const amrex::Array4<      amrex::Real>& coeffB_a, // main diagonal
                   const amrex::Array4<      amrex::Real>& inv_coeffB_a,
                   const amrex::Array4<const amrex::Real>& coeffC_a, // upper diagonal
                   const amrex::Array4<const amrex::Real>& RHS_a)
{
    // Forward sweep

    amrex::Real bet = coeffB_a(i,j,klo);

    for (int k(klo+1); k<=khi; ++k) {
        amrex::Real gam = coeffC_a(i,j,k-1) / bet;
        bet = coeffB_a(i,j,k) - coeffA_a(i,j,k)*gam;
        AMREX_ASSERT(bet != 0.0);
        coeffB_a(i,j,k) = bet;
    }

    for (int k(klo); k<=khi; ++k) {
        inv_coeffB_a(i,j,k) = 1.0 / coeffB_a(i,j,k);
    }

    //
    // Tridiagonal solve
    //
    soln_a(i,j,klo) = RHS_a(i,j,klo) * inv_coeffB_a(i,j,klo);

    for (int k(klo+1); k<=khi; ++k) {
        soln_a(i,j,k) = (RHS_a(i,j,k)-coeffA_a(i,j,k)*soln_a(i,j,k-1)) * inv_coeffB_a(i,j,k);
    }

    for (int k(khi-1); k>=klo; --k) {
        soln_a(i,j,k) -= ( coeffC_a(i,j,k) * inv_coeffB_a(i,j,k) ) * soln_a(i,j,k+1);
    }
}
#endif
