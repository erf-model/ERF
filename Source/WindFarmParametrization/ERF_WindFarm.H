#ifndef ERF_WINDFARM_H
#define ERF_WINDFARM_H

#include <ERF_DataStruct.H>
#include <AMReX_Geometry.H>
#include <AMReX_MultiFab.H>

#include "ERF_NullWindFarm.H"
#include "ERF_Fitch.H"
#include "ERF_EWP.H"
#include "ERF_SimpleAD.H"
#include "ERF_GeneralAD.H"

/**
 * @brief Manage a collection of wind turbine models.
 */
class WindFarm : public NullWindFarm {

public:

    WindFarm (){}

    virtual ~WindFarm () = default;

    WindFarm (int nlev,
              const WindFarmType& a_windfarm_type)
    {
        m_windfarm_model.resize(nlev);
        if (a_windfarm_type == WindFarmType::Fitch) {
            SetModel<Fitch>();
            amrex::Print() << "Fitch windfarm model!\n";
        }
        else if (a_windfarm_type == WindFarmType::EWP) {
            SetModel<EWP>();
            amrex::Print() << "EWP windfarm model!\n";
        }
        else if (a_windfarm_type == WindFarmType::SimpleAD) {
            SetModel<SimpleAD>();
            amrex::Print() << "Simplified actuator disk windfarm model!\n";
        }
        else if (a_windfarm_type == WindFarmType::GeneralAD) {
            SetModel<GeneralAD>();
            amrex::Print() << "Generalized actuator disk windfarm model!\n";
        }
         else {
            amrex::Abort("WindFarm: Dont know this windfarm_type!") ;
        }
    }

    /**
     * @brief Read turbine location and specification tables.
     * @param[in] windfarm_loc_table Path to the location table.
     * @param[in] windfarm_spec_table Path to the specification table.
     * @param[in] x_y Whether to use Cartesian coordinates.
     * @param[in] lat_lon Whether to use latitude and longitude coordinates.
     * @param[in] windfarm_x_shift X-coordinate shift.
     * @param[in] windfarm_y_shift Y-coordinate shift.
     */
    void read_tables (std::string windfarm_loc_table,
                      std::string windfarm_spec_table,
                      bool x_y, bool lat_lon,
                      const amrex::Real windfarm_x_shift = zero,
                      const amrex::Real windfarm_y_shift = zero);

    /**
     * @brief Initialize wind farm using latitude and longitude coordinates.
     * @param[in] windfarm_loc_table Path to the location table.
     * @param[in] windfarm_x_shift X-coordinate shift.
     * @param[in] windfarm_y_shift Y-coordinate shift.
     */
    void init_windfarm_lat_lon (const std::string windfarm_loc_table,
                                const amrex::Real windfarm_x_shift,
                                const amrex::Real windfarm_y_shift);

    /**
     * @brief Initialize wind farm using Cartesian coordinates.
     * @param[in] windfarm_loc_table Path to the location table.
     */
    void init_windfarm_x_y (const std::string windfarm_loc_table);

    /**
     * @brief Read the turbine location table.
     * @param[in] windfarm_loc_table Path to the location table.
     * @param[in] x_y Whether to use Cartesian coordinates.
     * @param[in] lat_lon Whether to use latitude and longitude coordinates.
     * @param[in] windfarm_x_shift X-coordinate shift.
     * @param[in] windfarm_y_shift Y-coordinate shift.
     */
    void read_windfarm_locations_table (const std::string windfarm_loc_table,
                                        bool x_y, bool lat_lon,
                                        const amrex::Real windfarm_x_shift = zero,
                                        const amrex::Real windfarm_y_shift = zero);

    /**
     * @brief Read the turbine specification table.
     * @param[in] windfarm_spec_table Path to the specification table.
     */
    void read_windfarm_spec_table (const std::string windfarm_spec_table);

    /**
     * @brief Read the turbine blade table.
     * @param[in] windfarm_blade_table Path to the blade table.
     */
    void read_windfarm_blade_table (const std::string windfarm_blade_table);

    /**
     * @brief Read the airfoil tables.
     * @param[in] windfarm_airfoil_tables Path to the airfoil tables.
     * @param[in] windfarm_blade_table Path to the blade table.
     */
    void read_windfarm_airfoil_tables (const std::string windfarm_airfoil_tables,
                                       const std::string windfarm_blade_table);

    /**
     * @brief Read extra turbine specifications.
     * @param[in] windfarm_spec_table_extra Path to the extra specification table.
     */
    void read_windfarm_spec_table_extra (const std::string windfarm_spec_table_extra);

    /**
     * @brief Fill the turbine count MultiFab.
     * @param[in] geom Geometry
     * @param[out] mf_Nturb MultiFab to fill with turbine counts.
     * @param[in,out] z_phys_nd Physical height MultiFab.
     */
    void fill_Nturb_multifab (const amrex::Geometry& geom,
                              amrex::MultiFab& mf_Nturb,
                              std::unique_ptr<amrex::MultiFab>& z_phys_nd);

    /**
     * @brief Fill the SMark MultiFab.
     * @param[in] geom Geometry
     * @param[out] mf_SMark MultiFab to fill with SMark values.
     * @param[in] sampling_distance_by_D Sampling distance normalized by rotor diameter.
     * @param[in] turb_disk_angle Turbine disk angle.
     * @param[in,out] z_phys_cc Physical height MultiFab.
     */
    void fill_SMark_multifab (const amrex::Geometry& geom,
                              amrex::MultiFab& mf_SMark,
                              const amrex::Real& sampling_distance_by_D,
                              const amrex::Real& turb_disk_angle,
                              std::unique_ptr<amrex::MultiFab>& z_phys_cc);

    /**
     * @brief Fill the SMark MultiFab for mesoscale models.
     * @param[in] geom Geometry
     * @param[out] mf_SMark MultiFab to fill with SMark values.
     * @param[in] mf_Nturb MultiFab containing turbine counts.
     * @param[in,out] z_phys_cc Physical height MultiFab.
     */
    void fill_SMark_multifab_mesoscale_models (const amrex::Geometry& geom,
                                          amrex::MultiFab& mf_SMark,
                                          const amrex::MultiFab& mf_Nturb,
                                          std::unique_ptr<amrex::MultiFab>& z_phys_cc);

    /**
     * @brief Write turbine locations to a VTK file.
     */
    void write_turbine_locations_vtk ();

    /**
     * @brief Write actuator disks to a VTK file.
     * @param[in] geom Geometry
     * @param[in] sampling_distance_by_D Sampling distance normalized by rotor diameter.
     */
    void write_actuator_disks_vtk (const amrex::Geometry& geom,
                                   const amrex::Real& sampling_distance_by_D);

    /**
     * @brief Advance the wind farm model in time.
     * @param[in] a_geom Geometry
     * @param[in] dt_advance Timestep.
     * @param[in,out] cons_in Conservative variables.
     * @param[in,out] mf_vars_windfarm Wind farm variables MultiFab.
     * @param[in,out] U_old Old x-velocity MultiFab.
     * @param[in,out] V_old Old y-velocity MultiFab.
     * @param[in,out] W_old Old z-velocity MultiFab.
     * @param[in] mf_Nturb Turbine count MultiFab.
     * @param[in] mf_SMark SMark MultiFab.
     * @param[in] time Current simulation time.
     */
    void advance (const amrex::Geometry& a_geom,
                  const double& dt_advance,
                  amrex::MultiFab& cons_in,
                  amrex::MultiFab& mf_vars_windfarm,
                  amrex::MultiFab& U_old,
                  amrex::MultiFab& V_old,
                  amrex::MultiFab& W_old,
                  const amrex::MultiFab& mf_Nturb,
                  const amrex::MultiFab& mf_SMark,
                  const double& time) override
    {
        m_windfarm_model[0]->advance(a_geom, dt_advance, cons_in, mf_vars_windfarm,
                                     U_old, V_old, W_old, mf_Nturb, mf_SMark, time);
    }

    /**
     * @brief Set general turbine specifications.
     * @param[in] a_rotor_rad Rotor radius.
     * @param[in] a_hub_height Hub height.
     * @param[in] a_thrust_coeff_standing Standing thrust coefficient.
     * @param[in] a_wind_speed Wind speed vector.
     * @param[in] a_thrust_coeff Thrust coefficient vector.
     * @param[in] a_power Power vector.
     */
    void set_turb_spec (const amrex::Real& a_rotor_rad, const amrex::Real& a_hub_height,
                        const amrex::Real& a_thrust_coeff_standing, const amrex::Vector<amrex::Real>& a_wind_speed,
                        const amrex::Vector<amrex::Real>& a_thrust_coeff,
                        const amrex::Vector<amrex::Real>& a_power) override
    {
        m_windfarm_model[0]->set_turb_spec(a_rotor_rad, a_hub_height, a_thrust_coeff_standing,
                  a_wind_speed, a_thrust_coeff, a_power);
    }

    /**
     * @brief Set turbine locations.
     * @param[in] a_xloc X locations.
     * @param[in] a_yloc Y locations.
     */
    void set_turb_loc (const amrex::Vector<amrex::Real>& a_xloc,
                       const amrex::Vector<amrex::Real>& a_yloc) override
    {
        m_windfarm_model[0]->set_turb_loc(a_xloc, a_yloc);
    }

    /**
     * @brief Set the turbine disk angle.
     * @param[in] a_turb_disk_angle Disk angle.
     */
    void set_turb_disk_angle (const amrex::Real& a_turb_disk_angle) override
    {
        m_windfarm_model[0]->set_turb_disk_angle(a_turb_disk_angle);
    }

    /**
     * @brief Set blade specifications.
     * @param[in] a_bld_rad_loc Blade radial locations.
     * @param[in] a_bld_twist Blade twist.
     * @param[in] a_bld_chord Blade chord.
     */
    void set_blade_spec (const amrex::Vector<amrex::Real>& a_bld_rad_loc,
                         const amrex::Vector<amrex::Real>& a_bld_twist,
                         const amrex::Vector<amrex::Real>& a_bld_chord) override
    {
        m_windfarm_model[0]->set_blade_spec(a_bld_rad_loc, a_bld_twist, a_bld_chord);
    }

    /**
     * @brief Set blade airfoil specifications.
     * @param[in] a_bld_airfoil_aoa Angle of attack.
     * @param[in] a_bld_airfoil_Cl Lift coefficient.
     * @param[in] a_bld_airfoil_Cd Drag coefficient.
     */
    void set_blade_airfoil_spec (const amrex::Vector<amrex::Vector<amrex::Real>>& a_bld_airfoil_aoa,
                                 const amrex::Vector<amrex::Vector<amrex::Real>>& a_bld_airfoil_Cl,
                                 const amrex::Vector<amrex::Vector<amrex::Real>>& a_bld_airfoil_Cd) override
    {
        m_windfarm_model[0]->set_blade_airfoil_spec(a_bld_airfoil_aoa, a_bld_airfoil_Cl, a_bld_airfoil_Cd);
    }

    /**
     * @brief Set extra turbine specifications.
     * @param[in] a_velocity Velocity.
     * @param[in] a_C_P Power coefficient.
     * @param[in] a_C_T Thrust coefficient.
     * @param[in] a_rotor_RPM Rotor RPM.
     * @param[in] a_blade_pitch Blade pitch.
     */
    void set_turb_spec_extra (const amrex::Vector<amrex::Real>& a_velocity,
                              const amrex::Vector<amrex::Real>& a_C_P,
                              const amrex::Vector<amrex::Real>& a_C_T,
                              const amrex::Vector<amrex::Real>& a_rotor_RPM,
                              const amrex::Vector<amrex::Real>& a_blade_pitch) override
    {
        m_windfarm_model[0]->set_turb_spec_extra(a_velocity, a_C_P, a_C_T, a_rotor_RPM, a_blade_pitch);
    }

protected:

    amrex::Vector<amrex::Real> xloc, yloc, zloc;
    amrex::Real my_turb_disk_angle;
    amrex::Real hub_height, rotor_rad, thrust_coeff_standing, nominal_power;
    amrex::Vector<amrex::Real> wind_speed, thrust_coeff, power;
    amrex::Vector<amrex::Real> bld_rad_loc, bld_twist, bld_chord;
    amrex::Vector<amrex::Vector<amrex::Real>> bld_airfoil_aoa, bld_airfoil_Cl, bld_airfoil_Cd;
    int n_bld_sections;
    amrex::Vector<amrex::Real> velocity, C_P, C_T, rotor_RPM, blade_pitch;

   /*! \brief Create and set the specified windfarm model */
    template<class NewWindFarmModel>
    void SetModel ()
    {
        for (int lev(0); lev<m_windfarm_model.size(); ++lev) {
            m_windfarm_model[lev] = std::make_unique<NewWindFarmModel>();
        }
    }

private:
    amrex::Vector<std::unique_ptr<NullWindFarm>> m_windfarm_model; /*!< windfarm model */
};

#endif

