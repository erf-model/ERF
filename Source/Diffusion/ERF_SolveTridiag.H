        // Forward sweep

        Real bet = coeffB_a(i,j,klo);

        for (int k(klo+1); k<=khi; ++k) {
            Real gam = coeffC_a(i,j,k-1) / bet;
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
