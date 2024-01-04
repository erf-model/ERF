#include <string>
#include <ERF.H>
#include <ERFPC.H>

#ifdef ERF_USE_PARTICLES

using namespace amrex;

/*! Read tracer and hydro particles parameters */
void ERF::readTracersParams()
{
    ParmParse pp(pp_prefix);

    m_use_tracer_particles = 0;
    m_use_hydro_particles = 0;

    pp.query(std::string("use_"+ERFParticleNames::tracers).c_str(), m_use_tracer_particles);
    pp.query(std::string("use_"+ERFParticleNames::hydro).c_str(), m_use_hydro_particles);

    if (m_use_tracer_particles) {
        particleData.addName(ERFParticleNames::tracers);
    }

    if (m_use_hydro_particles) {
        particleData.addName(ERFParticleNames::hydro);
    }
    return;
}

/*! Initialize tracer and hydro particles */
void ERF::initializeTracers( ParGDBBase* a_gdb,
                             const Vector<std::unique_ptr<MultiFab>>& a_z_phys_nd )
{
    auto& namelist_unalloc( particleData.getNamesUnalloc() );

    for (auto it = namelist_unalloc.begin(); it != namelist_unalloc.end(); ++it) {

        std::string species_name( *it );

        if (species_name == ERFParticleNames::tracers) {

            AMREX_ASSERT(m_use_tracer_particles);
            ERFPC* pc = new ERFPC(a_gdb, ERFParticleNames::tracers);
            pc->InitializeParticles(a_z_phys_nd[0]);
            amrex::Print() << "Initialized " << pc->TotalNumberOfParticles() << " tracer particles.\n";
            particleData.pushBack(ERFParticleNames::tracers, pc);

        } else if (species_name == ERFParticleNames::hydro) {

            AMREX_ASSERT(m_use_hydro_particles);
            ERFPC* pc = new ERFPC(a_gdb, ERFParticleNames::hydro);
            pc->InitializeParticles(a_z_phys_nd[0]);
            amrex::Print() << "Initialized " << pc->TotalNumberOfParticles() << " hydro particles.\n";
            particleData.pushBack(ERFParticleNames::hydro, pc);

        }
    }

    if (m_use_tracer_particles) namelist_unalloc.remove( ERFParticleNames::tracers );
    if (m_use_hydro_particles)  namelist_unalloc.remove( ERFParticleNames::hydro );

    return;
}

/*! Evolve tracers and hydro particles for one time step*/
void ERF::evolveTracers( int                                        a_lev,
                         Real                                       a_dt_lev,
                         Vector<Vector<MultiFab>>&                  a_vars_new,
                         const Vector<std::unique_ptr<MultiFab>>&   a_z_phys_nd )
{
    if (m_use_tracer_particles) {
      particleData[ERFParticleNames::tracers]->EvolveParticles(  a_lev,
                                                                 a_dt_lev,
                                                                 a_vars_new,
                                                                 a_z_phys_nd );
    }
    if (m_use_hydro_particles) {
      particleData[ERFParticleNames::hydro]->EvolveParticles( a_lev,
                                                              a_dt_lev,
                                                              a_vars_new,
                                                              a_z_phys_nd );
    }
    return;
}

#endif
