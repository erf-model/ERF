#include <string>
#include <ERF.H>
#include <ERFPC.H>

#ifdef ERF_USE_PARTICLES

using namespace amrex;

/*! Read tracer particle parameters */
void ERF::readTracersParams ()
{
    ParmParse pp(pp_prefix);

    m_use_tracer_particles = 0;

    pp.query(std::string("use_"+ERFParticleNames::tracers).c_str(), m_use_tracer_particles);

    if (m_use_tracer_particles) {
        particleData.addName(ERFParticleNames::tracers);
    }
    return;
}

/*! Initialize tracer particles */
void ERF::initializeTracers ( ParGDBBase* a_gdb,
                              const Vector<std::unique_ptr<MultiFab>>& a_z_phys_nd,
                              const Real time)
{
    auto& namelist_unalloc( particleData.getNamesUnalloc() );

    for (auto it = namelist_unalloc.begin(); it != namelist_unalloc.end(); ++it) {

        std::string species_name( *it );

        if ( (species_name == ERFParticleNames::tracers) && (!particleData.HasSpecies(species_name)) ) {

            AMREX_ASSERT(m_use_tracer_particles);
            ERFPC* pc = new ERFPC(a_gdb, ERFParticleNames::tracers);
            pc->InitializeParticles(time,a_z_phys_nd[0]);
            if (pc->TotalNumberOfParticles() > 0) {
                amrex::Print() << "Initialized " << pc->TotalNumberOfParticles() << " tracer particles.\n";
                particleData.pushBack(ERFParticleNames::tracers, pc);
            }
        }
    }

    if (m_use_tracer_particles && particleData.HasSpecies(ERFParticleNames::tracers)) {
        namelist_unalloc.remove( ERFParticleNames::tracers );
    }

    return;
}

/*! Restart tracer particles */
void ERF::restartTracers ( ParGDBBase* a_gdb,
                           const std::string& a_fname )
{
    auto& namelist_unalloc( particleData.getNamesUnalloc() );

    for (auto it = namelist_unalloc.begin(); it != namelist_unalloc.end(); ++it) {

        std::string species_name( *it );

        if (species_name == ERFParticleNames::tracers)
        {
            AMREX_ASSERT(m_use_tracer_particles);
            std::string m_name = ERFParticleNames::tracers;
            std::string HeaderFileName(a_fname + "/" + m_name + "/Header");
            if (FileExists(HeaderFileName)) {
                ERFPC* pc = new ERFPC(a_gdb, m_name);
                pc->Restart(a_fname, m_name);
                amrex::Print() << "Restarted " << pc->TotalNumberOfParticles() << " tracer particles.\n";
                particleData.pushBack(m_name, pc);
            }
        }
    }

    if (m_use_tracer_particles && particleData.HasSpecies(ERFParticleNames::tracers)) {
        namelist_unalloc.remove( ERFParticleNames::tracers );
    }

    return;
}

/*! Evolve tracers particles for one time step*/
void ERF::evolveTracers ( int                                        a_lev,
                          Real                                       a_dt_lev,
                          Vector<Vector<MultiFab>>&                  a_vars_new,
                          const Vector<std::unique_ptr<MultiFab>>&   a_z_phys_nd )
{
    if ( m_use_tracer_particles && particleData.HasSpecies(ERFParticleNames::tracers) ) {
      particleData[ERFParticleNames::tracers]->EvolveParticles(  a_lev,
                                                                 a_dt_lev,
                                                                 a_vars_new,
                                                                 a_z_phys_nd );
    }
    return;
}

#endif
