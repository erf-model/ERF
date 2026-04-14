
    ParmParse pp_for_pert_vels("erf");
    InitType init_type;
    pp_for_pert_vels.query_enum_case_insensitive("init_type",init_type);

    ParmParse pp_for_pert_vels_prob("prob");
    Real rho_0 =   one; pp_for_pert_vels_prob.query("rho_0", rho_0);

    if (init_type == InitType::Uniform)
    {
        Real rot_time_period;
        pp_for_pert_vels.get("rotational_time_period", rot_time_period);
        Real coriolis_factor = amrex::Real(4.0) * PI / rot_time_period;

        Real Az;
        pp_for_pert_vels.get("dynamic_viscosity", Az); // dynamic viscosity [kg-m/s]
        Az = Az / rho_0; // kinematic viscosity [m^2/s]

        Vector<Real> abl_geo_wind(3);
        pp_for_pert_vels.queryarr("abl_geo_wind",abl_geo_wind);
        AMREX_ALWAYS_ASSERT_WITH_MESSAGE(
                                         (amrex::Math::abs(abl_geo_wind[1]) < amrex::Real(1.0e-15)) &&
                                         (amrex::Math::abs(abl_geo_wind[2]) < amrex::Real(1.0e-15)),
                                         "Ekman Spiral uses geostrophic forcing of the form (V_0, 0, 0)");
        const Real u_0 = abl_geo_wind[0];

        const Real DE = std::sqrt(two * Az / coriolis_factor);
        const Real a = one / DE;

        // Set the x-velocity
        ParallelFor(xbx, [=] AMREX_GPU_DEVICE(int i, int j, int k) noexcept {
            const auto dx = geomdata.CellSize();
            const Real z = (k + myhalf) * dx[2];

            // Set the x-velocity
            x_vel_pert(i, j, k) = u_0 * (one - std::exp(-a * z) * std::cos(-a * z));
        });

        // Set the y-velocity
        ParallelFor(ybx, [=] AMREX_GPU_DEVICE(int i, int j, int k) noexcept {
            const auto dx = geomdata.CellSize();
            const Real z = (k + myhalf) * dx[2];

            // Set the y-velocity
            y_vel_pert(i, j, k) = -u_0 * std::exp(-a * z) * std::sin(-a * z);
        });

        // Set the z-velocity
        ParallelFor(zbx, [=] AMREX_GPU_DEVICE(int i, int j, int k) noexcept {
            z_vel_pert(i, j, k) = zero;
        });
    } // end if uniform
