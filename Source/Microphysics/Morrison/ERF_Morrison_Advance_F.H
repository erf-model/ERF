          double dummy_reflectivity = Real(0);
          double* dummy_reflectivity_ptr = &dummy_reflectivity;

          // If using Fortran version, update the Fortran module variable as well
          set_morrison_ndcnst_c(m_ndcnst);

          const int ilom = grown_box.loVect()[0];
          const int ihim = grown_box.hiVect()[0];
          const int jlom = grown_box.loVect()[1];
          const int jhim = grown_box.hiVect()[1];
          const int klom = grown_box.loVect()[2];
          const int khim = grown_box.hiVect()[2];

#if defined(ERF_USE_MORR_FORT) && defined(AMREX_USE_GPU)
    Arena* Arena_Used = The_Pinned_Arena();
#else
    Arena* Arena_Used = The_Arena();
#endif

          // Create dummy arrays for cumulus tendencies (if needed)
          FArrayBox qrcuten_fab(grown_box, 1, Arena_Used);
          FArrayBox qscuten_fab(grown_box, 1, Arena_Used);
          FArrayBox qicuten_fab(grown_box, 1, Arena_Used);
          auto const& qrcuten_arr = qrcuten_fab.array();
          auto const& qscuten_arr = qscuten_fab.array();
          auto const& qicuten_arr = qicuten_fab.array();

          // Initialize tendencies to Real(0) (no cumulus parameterization in this example)
          ParallelFor(grown_box, [=] AMREX_GPU_DEVICE (int i, int j, int k) {
            qrcuten_arr(i,j,k) = Real(0);
            qscuten_arr(i,j,k) = Real(0);
            qicuten_arr(i,j,k) = Real(0);
          });

          // WRF-Chem related variables (optional)
          bool flag_qndrop = false;  // Flag to indicate droplet number prediction

          // Now create arrays for other optional variables
          FArrayBox rainprod_fab(grown_box, 1, Arena_Used);
          FArrayBox evapprod_fab(grown_box, 1, Arena_Used);
          FArrayBox qlsink_fab(grown_box, 1, Arena_Used);
          FArrayBox precr_fab(grown_box, 1, Arena_Used);
          FArrayBox preci_fab(grown_box, 1, Arena_Used);
          FArrayBox precs_fab(grown_box, 1, Arena_Used);
          FArrayBox precg_fab(grown_box, 1, Arena_Used);

          auto const& rainprod_arr = rainprod_fab.array();
          auto const& evapprod_arr = evapprod_fab.array();
          auto const& qlsink_arr = qlsink_fab.array();
          auto const& precr_arr = precr_fab.array();
          auto const& preci_arr = preci_fab.array();
          auto const& precs_arr = precs_fab.array();
          auto const& precg_arr = precg_fab.array();

          // Initialize WRF-Chem arrays to Real(0)
          ParallelFor(grown_box, [=] AMREX_GPU_DEVICE (int i, int j, int k) {
            rainprod_arr(i,j,k) = Real(0);
            evapprod_arr(i,j,k) = Real(0);
            qlsink_arr(i,j,k)   = Real(0);
            precr_arr(i,j,k)    = Real(0);
            preci_arr(i,j,k)    = Real(0);
            precs_arr(i,j,k)    = Real(0);
            precg_arr(i,j,k)    = Real(0);
          });

          mp_morr_two_moment_c
          (
              1,  // ITIMESTEP - Use 1 for simplicity

              // 3D arrays in Fortran expected order (assume column-major for Fortran)
              theta_arr.dataPtr(),      // TH
              qv_arr.dataPtr(),         // QV
              qcl_arr.dataPtr(),        // QC
              qpr_arr.dataPtr(),        // QR
              qci_arr.dataPtr(),        // QI
              qps_arr.dataPtr(),        // QS
              qpg_arr.dataPtr(),        // QG
              ni_arr.dataPtr(),         // NI
              ns_arr.dataPtr(),         // NS
              nr_arr.dataPtr(),         // NR
              ng_arr.dataPtr(),         // NG

              rho_arr.dataPtr(),        // RHO
              pii_arr.dataPtr(),        // PII (Exner function)
              pres_arr.dataPtr(),       // P (in hPa, convert if needed)
              dt,                       // DT_IN
              dz_arr.dataPtr(),         // DZ
              w_arr.dataPtr(),          // W (vertical velocity)

              // 2D arrays for precipitation accounting
              rain_accum_arr.dataPtr(), // RAINNC
              rainncv_arr.dataPtr(),    // RAINNCV
              sr_arr.dataPtr(),         // SR
              snow_accum_arr.dataPtr(), // SNOWNC
              snowncv_arr.dataPtr(),    // SNOWNCV
              graup_accum_arr.dataPtr(),// GRAUPELNC
              graupelncv_arr.dataPtr(), // GRAUPELNCV

              // Radar reflectivity
              dummy_reflectivity_ptr,  // refl_10cm
              true,                     // diagflag
              0,   // do_radar_ref

              // Cumulus tendencies
              qrcuten_arr.dataPtr(),    // qrcuten
              qscuten_arr.dataPtr(),    // qscuten
              qicuten_arr.dataPtr(),    // qicuten

              // WRF-Chem flags
              flag_qndrop,              // F_QNDROP
              nullptr,                  // qndrop (not used here)
              // ht_arr.dataPtr(),         // HT (terrain height - not used)

              // Domain dimensions
              ilo, ihi, jlo, jhi, klo, khi,  // IDS,IDE,JDS,JDE,KDS,KDE
              ilom, ihim, jlom, jhim, klom, khim,  // IMS,IME,JMS,JME,KMS,KME
              ilo, ihi, jlo, jhi, klo, khi,  // ITS,ITE,JTS,JTE,KTS,KTE

              // Optional WRF-Chem outputs
              false,                    // wetscav_on
              rainprod_arr.dataPtr(),   // rainprod
              evapprod_arr.dataPtr(),   // evapprod
              qlsink_arr.dataPtr(),     // QLSINK
              precr_arr.dataPtr(),      // PRECR
              preci_arr.dataPtr(),      // PRECI
              precs_arr.dataPtr(),      // PRECS
              precg_arr.dataPtr()       // PRECG
          );
