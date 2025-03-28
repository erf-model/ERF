#include <string>
#include <vector>
#include <memory>

#include <AMReX_FArrayBox.H>
#include <AMReX_Geometry.H>
#include <AMReX_TableData.H>
#include <AMReX_MultiFabUtil.H>

#include "ERF_Constants.H"
#include "ERF_MicrophysicsUtils.H"
#include "ERF_IndexDefines.H"
#include "ERF_DataStruct.H"
#include "ERF_NullMoist.H"
#include "ERF_Morrison.H"
#include <ERF_Morrison_Fortran_Interface.H>

    // wrapper to do all the updating
    void
    Morrison::Advance (const amrex::Real& dt_advance,
                       const SolverChoice& sc)
    {
        dt = dt_advance;
#if 0
        this->Cloud(sc);
        this->IceFall(sc);
        this->Precip(sc);
        this->PrecipFall(sc);
#else
        // Store timestep
        dt = dt_advance;

        // Loop through the grids
        for (amrex::MFIter mfi(*mic_fab_vars[MicVar_Morr::qcl],TileNoZ()); mfi.isValid(); ++mfi)
        {
          const amrex::Box& box = mfi.tilebox();

          // Get array data from class member variables
          auto const& theta_arr = mic_fab_vars[MicVar_Morr::theta]->array(mfi);
          auto const& qv_arr = mic_fab_vars[MicVar_Morr::qv]->array(mfi);
          auto const& qcl_arr = mic_fab_vars[MicVar_Morr::qcl]->array(mfi);
          auto const& qpr_arr = mic_fab_vars[MicVar_Morr::qpr]->array(mfi);
          auto const& qci_arr = mic_fab_vars[MicVar_Morr::qci]->array(mfi);
          auto const& qps_arr = mic_fab_vars[MicVar_Morr::qps]->array(mfi);
          auto const& qpg_arr = mic_fab_vars[MicVar_Morr::qpg]->array(mfi);
          auto const& ni_arr = mic_fab_vars[MicVar_Morr::ni]->array(mfi);
          auto const& ns_arr = mic_fab_vars[MicVar_Morr::ns]->array(mfi);
          auto const& nr_arr = mic_fab_vars[MicVar_Morr::nr]->array(mfi);
          auto const& ng_arr = mic_fab_vars[MicVar_Morr::ng]->array(mfi);
          auto const& rho_arr = mic_fab_vars[MicVar_Morr::rho]->array(mfi);
          auto const& pres_arr = mic_fab_vars[MicVar_Morr::pres]->array(mfi);
          auto const& tabs_arr = mic_fab_vars[MicVar_Morr::tabs]->array(mfi);
          auto const& rain_accum_arr = mic_fab_vars[MicVar_Morr::rain_accum]->array(mfi);
          auto const& snow_accum_arr = mic_fab_vars[MicVar_Morr::snow_accum]->array(mfi);
          auto const& graup_accum_arr = mic_fab_vars[MicVar_Morr::graup_accum]->array(mfi);
          auto const& w_arr = mic_fab_vars[MicVar_Morr::omega]->array(mfi);

          // Get radar reflectivity array if radar diagnostics enabled
          //        auto const& refl_arr = m_do_radar_ref ? m_radar->array(mfi) : nullptr;
          //        auto const& refl_arr = m_radar->array(mfi);

          // Extract box dimensions
          const int ilo = box.loVect()[0];
          const int ihi = box.hiVect()[0];
          const int jlo = box.loVect()[1];
          const int jhi = box.hiVect()[1];
          const int klo = box.loVect()[2];
          const int khi = box.hiVect()[2];

          amrex::Box grown_box(box); grown_box.grow(3);

          const int ilom = grown_box.loVect()[0];
          const int ihim = grown_box.hiVect()[0];
          const int jlom = grown_box.loVect()[1];
          const int jhim = grown_box.hiVect()[1];
          const int klom = grown_box.loVect()[2];
          const int khim = grown_box.hiVect()[2];

          // Calculate Exner function (PII) to convert potential temperature to temperature
          // PII = (P/P0)^(R/cp)
          amrex::FArrayBox pii_fab(grown_box, 1);
          auto const& pii_arr = pii_fab.array();

          const amrex::Real p0 = 100000.0; // Reference pressure (Pa)

          const amrex::Real rdcp = m_rdOcp; // R/cp ratio

          // Calculate Exner function
          amrex::ParallelFor(grown_box, [=] AMREX_GPU_DEVICE (int i, int j, int k) {
            // NOTE: the Morrison Fortran version uses Pa not hPa so we didn't divide p by 100
            //       so we don't need to multiply by 100 here
            pii_arr(i,j,k) = std::pow((pres_arr(i,j,k)) / p0, rdcp);
          });

          // Create arrays for height differences (dz)
          amrex::FArrayBox dz_fab(grown_box, 1);
          auto const& dz_arr = dz_fab.array();

          // Calculate height differences
          const amrex::Real dz_val = m_geom.CellSize(m_axis);
          amrex::ParallelFor(grown_box, [=] AMREX_GPU_DEVICE (int i, int j, int k) {
            dz_arr(i,j,k) = dz_val;
          });

          // Arrays to store precipitation rates
          amrex::FArrayBox    rainncv_fab(grown_box, 1);
          amrex::FArrayBox         sr_fab(grown_box, 1);     // Ratio of snow to total precipitation
          amrex::FArrayBox    snowncv_fab(grown_box, 1);
          amrex::FArrayBox graupelncv_fab(grown_box, 1);

          auto const& rainncv_arr = rainncv_fab.array();
          auto const& sr_arr = sr_fab.array();
          auto const& snowncv_arr = snowncv_fab.array();
          auto const& graupelncv_arr = graupelncv_fab.array();

          // Initialize precipitation rate arrays to zero
          rainncv_fab.setVal(0.0);
          sr_fab.setVal(0.0);
          snowncv_fab.setVal(0.0);
          graupelncv_fab.setVal(0.0);

          // Create terrain height array (not actually used by Morrison scheme)
          amrex::FArrayBox ht_fab(amrex::Box(amrex::IntVect(ilo, jlo, 0), amrex::IntVect(ihi, jhi, 0)), 1);
          auto const& ht_arr = ht_fab.array();
          ht_fab.setVal(0.0);  // Not used by Morrison scheme

          // Create dummy arrays for cumulus tendencies (if needed)
          amrex::FArrayBox qrcuten_fab(grown_box, 1);
          amrex::FArrayBox qscuten_fab(grown_box, 1);
          amrex::FArrayBox qicuten_fab(grown_box, 1);
          auto const& qrcuten_arr = qrcuten_fab.array();
          auto const& qscuten_arr = qscuten_fab.array();
          auto const& qicuten_arr = qicuten_fab.array();

          // Initialize tendencies to zero (no cumulus parameterization in this example)
          qrcuten_fab.setVal(0.0);
          qscuten_fab.setVal(0.0);
          qicuten_fab.setVal(0.0);

          // WRF-Chem related variables (optional)
          bool flag_qndrop = false;  // Flag to indicate droplet number prediction

          // Now create arrays for other optional variables
          amrex::FArrayBox rainprod_fab(grown_box, 1);
          amrex::FArrayBox evapprod_fab(grown_box, 1);
          amrex::FArrayBox qlsink_fab(grown_box, 1);
          amrex::FArrayBox precr_fab(grown_box, 1);
          amrex::FArrayBox preci_fab(grown_box, 1);
          amrex::FArrayBox precs_fab(grown_box, 1);
          amrex::FArrayBox precg_fab(grown_box, 1);

          auto const& rainprod_arr = rainprod_fab.array();
          auto const& evapprod_arr = evapprod_fab.array();
          auto const& qlsink_arr = qlsink_fab.array();
          auto const& precr_arr = precr_fab.array();
          auto const& preci_arr = preci_fab.array();
          auto const& precs_arr = precs_fab.array();
          auto const& precg_arr = precg_fab.array();

          // Initialize WRF-Chem arrays to zero
          rainprod_fab.setVal(0.0);
          evapprod_fab.setVal(0.0);
          qlsink_fab.setVal(0.0);
          precr_fab.setVal(0.0);
          preci_fab.setVal(0.0);
          precs_fab.setVal(0.0);
          precg_fab.setVal(0.0);

          // Prepare data pointers for Fortran call
          // These would be passed directly to the Fortran interface
          double dummy_reflectivity = 0.0;
          double* dummy_reflectivity_ptr = &dummy_reflectivity;
          // Example call (pseudo-code - actual interface would depend on your Fortran interop setup)

          //          amrex::Print()<<amrex::FArrayBox(qv_arr)<<std::endl;
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
              ht_arr.dataPtr(),         // HT (terrain height - not used)

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
          //          amrex::Print()<<amrex::FArrayBox(qv_arr)<<std::endl;
          // After the call, all fields are updated
          // We don't need to copy results back since we passed direct pointers
          // to our class member arrays
        }
#endif
    }
