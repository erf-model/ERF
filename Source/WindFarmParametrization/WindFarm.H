#ifndef WINDFARM_H
#define WINDFARM_H

#include "NullWindFarm.H"
#include "SimpleAD.H"
#include <DataStruct.H>
#include <AMReX_Geometry.H>
#include <AMReX_MultiFab.H>
#include <AMReX_Gpu.H>

class WindFarm : public NullWindFarm {

public:

    WindFarm(){}

	WindFarm (int nlev, 
			  const WindFarmType& a_windfarm_type) 
	{
		m_windfarm_model.resize(nlev);
		if (a_windfarm_type == WindFarmType::SimpleAD) {
            SetModel<SimpleAD>();
            amrex::Print() << "Simple actuator disk windfarm model!\n";
        } else {
            amrex::Abort("WindFarm: Dont know this windfarm_type!") ;
        }
	}

    void read_tables (std::string windfarm_loc_table,
                       std::string windfarm_spec_table,
                       bool x_y, bool lat_lon,
                       const amrex::Real latitude_lo=0.0,
                       const amrex::Real longitude_lo=0.0);

    void init_windfarm_lat_lon(const std::string windfarm_loc_table,
                               const amrex::Real latitude_lo,
                               const amrex::Real longitude_lo);

    void init_windfarm_x_y(const std::string windfarm_loc_table);

    void read_windfarm_locations_table(const std::string windfarm_loc_table,
										  bool x_y, bool lat_lon,
										  const amrex::Real latitude_lo=0.0,
										  const amrex::Real longitude_lo=0.0);

    void read_windfarm_spec_table(const std::string windfarm_spec_table);

	void fill_Nturb_multifab(const amrex::Geometry& geom, amrex::MultiFab& mf_Nturb);

    void write_turbine_locations_vtk();

    void write_actuator_disks_vtk();

	void advance (int lev,
                        const amrex::Geometry& geom,
                        const amrex::Real& dt_advance,
                        amrex::MultiFab& cons_in,
                        amrex::MultiFab& U_old, amrex::MultiFab& V_old, amrex::MultiFab& W_old,
                        amrex::MultiFab& mf_vars_windfarm, const amrex::MultiFab& mf_Nturb) override
	{
		m_windfarm_model[lev]->advance(lev, geom, dt_advance, cons_in, U_old, V_old, W_old,
                      				   mf_vars_windfarm, mf_Nturb);	
	}

	/*void source_terms_cellcentered (const amrex::Geometry& geom,
                                          const amrex::MultiFab& cons_in,
                                          const amrex::MultiFab& U_old, const amrex::MultiFab& V_old, const amrex::MultiFab& W_old,
                                          amrex::MultiFab& mf_vars_ewp, const amrex::MultiFab& mf_Nturb)
	{
		m_windfarm_model[lev]->source_terms_cellcentered(geom, cons_in, U_old, V_old, W_old, mf_vars_simpleAD, mf_Nturb);
	}

	void update (const amrex::Real& dt_advance,
                      amrex::MultiFab& cons_in,
                      amrex::MultiFab& U_old, amrex::MultiFab& V_old,
                      const amrex::MultiFab& mf_vars_ewp)
	{
		 m_windfarm_model[lev]->update(dt_advance, cons_in, U_old, V_old, mf_vars_simpleAD);
	}*/


protected:

   /*! \brief Create and set the specified windfarm model */
    template<class NewWindFarmModel>
    void SetModel ()
    {
		std::cout << "Setting the wind farm model with " << m_windfarm_model.size() << " levels " << "\n";
        for (int lev(0); lev<m_windfarm_model.size(); ++lev) {
            m_windfarm_model[lev] = std::make_unique<NewWindFarmModel>();
        }
    }

private:
    amrex::Vector<std::unique_ptr<NullWindFarm>> m_windfarm_model; /*!< windfarm model */

};

#endif

