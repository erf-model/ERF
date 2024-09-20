#ifndef ERF_PARTICLE_DATA_H_
#define ERF_PARTICLE_DATA_H_

#ifdef ERF_USE_PARTICLES

#include <map>
#include <vector>
#include <list>
#include <string>
#include <iostream>

#include <AMReX_ParmParse.H>
#include <AMReX_Print.H>
#include <AMReX_Vector.H>
#include <AMReX_Gpu.H>

#include <ERFPC.H>

/**
 * Container holding many of the particle-related data and options
 */

typedef std::map<std::string, ERFPC*> ParticleSpeciesMap;
typedef std::vector<std::string> ParticlesNamesVector;
typedef std::list<std::string> ParticlesNamesList;

class ParticleData
{
    public:

        /*! Constructor */
        ParticleData ()
        {
            BL_PROFILE("ParticleData::ParticleData()");

            amrex::ParmParse pp("particles");
            m_disable_particle_op = false;
            pp.query("disable_plt", m_disable_particle_op);

            m_particle_species.clear();
            m_namelist.clear();
            m_namelist_unalloc.clear();
        }

        /*! Destructor */
        ~ParticleData ()
        {
            BL_PROFILE("ParticleData::~ParticleData()");
            for (ParticlesNamesVector::size_type i = 0; i < m_namelist.size(); i++) {
                auto particles( m_particle_species[m_namelist[i]] );
                delete particles;
            }
            m_particle_species.clear();
            m_namelist.clear();
            m_namelist_unalloc.clear();
        }

        /*! Write particle info to plot files */
        void writePlotFile ( const std::string& a_fname ) const
        {
            BL_PROFILE("ParticleData::writePlotFile");
            if (!m_disable_particle_op) {
                for (ParticlesNamesVector::size_type i = 0; i < m_namelist.size(); i++) {
                    auto name( m_namelist[i] );
                    auto particles( m_particle_species.at(name) );
                    particles->Checkpoint( a_fname, name, true, particles->varNames() );
                }
            }
        }

        /*! Write checkpoint files */
        void Checkpoint ( const std::string& a_fname ) const
        {
            BL_PROFILE("ParticleData::Checkpoint()");
            for (ParticlesNamesVector::size_type i = 0; i < m_namelist.size(); i++) {
                auto name( m_namelist[i] );
                auto particles( m_particle_species.at(name) );
                particles->Checkpoint( a_fname, name, true, particles->varNames() );
            }
        }

        /*! Read from restart file */
        void Restart ( amrex::ParGDBBase* a_gdb, const std::string& a_fname )
        {
            BL_PROFILE("ParticleData::Restart()");
            AMREX_ASSERT(isEmpty());
            for (auto it = m_namelist_unalloc.begin(); it != m_namelist_unalloc.end(); ++it) {
                std::string species_name( *it );
                ERFPC* pc = new ERFPC( a_gdb, species_name );
                pc->Restart(a_fname, species_name );
                pushBack( species_name, pc );
            }
            m_namelist_unalloc.clear();
        }

        /*! Get mesh plot quantities from each particle container */
        void GetMeshPlotVarNames ( amrex::Vector<std::string>& a_names ) const
        {
            BL_PROFILE("ParticleData::GetMeshPlotVarNames()");
            a_names.clear();
            for (ParticlesNamesVector::size_type i = 0; i < m_namelist.size(); i++) {
                auto name( m_namelist[i] );
                auto particles( m_particle_species.at(name) );

                auto var_names = particles->meshPlotVarNames();
                for (int n = 0; n < var_names.size(); n++) {
                    a_names.push_back( std::string(name+"_"+var_names[n]) );
                }
            }
        }

        void GetMeshPlotVar (   const std::string& a_var_name,
                                amrex::MultiFab&   a_mf,
                                const int          a_lev )
        {
            BL_PROFILE("ParticleData::GetMeshPlotVar()");
            for (ParticlesNamesVector::size_type i = 0; i < m_namelist.size(); i++) {
                auto particle_name( m_namelist[i] );
                auto particles( m_particle_species.at(particle_name) );

                auto particle_var_names = particles->meshPlotVarNames();

                for (int n = 0; n < particle_var_names.size(); n++) {

                    std::string var_name = particle_name+"_"+particle_var_names[n];
                    if ( var_name == a_var_name ) {
                        particles->computeMeshVar(particle_var_names[n], a_mf, a_lev);
                        return;
                    }
                }
            }
            amrex::Abort("Requested var_name not found in ParticleData::GetMeshPlotVar");
        }

        /*! Redistribute/rebalance particles data */
        inline void Redistribute ()
        {
            BL_PROFILE("ParticleData::Redistribute()");
            for (ParticlesNamesVector::size_type i = 0; i < m_namelist.size(); i++) {
                m_particle_species[m_namelist[i]]->Redistribute();
            }
        }

        /*! Get species of a given name */
        inline ERFPC* GetSpecies ( const std::string& a_name )
        {
            BL_PROFILE("ParticleData::GetSpecies()");
            ParticleSpeciesMap::iterator it (m_particle_species.find(a_name));
            if (it == m_particle_species.end()) {
                amrex::Print()  << "ERROR: unable to find particle species with name \""
                                << a_name << "\"!";
                return nullptr;
            } else {
                return it->second;
            }
        }

        /*! accessor */
        inline ERFPC* operator[] ( const std::string& a_name )
        {
            BL_PROFILE("ParticleData::operator[]");
            ParticleSpeciesMap::iterator it (m_particle_species.find(a_name));
            if (it == m_particle_species.end()) {
                amrex::Print()  << "ERROR: unable to find particle species with name \""
                                << a_name << "\"!";
                return nullptr;
            } else {
                return it->second;
            }
        }

        /*! Get species of a given name (const version) */
        inline const ERFPC* GetSpecies ( const std::string& a_name ) const
        {
            BL_PROFILE("ParticleData::GetSpecies()");
            ParticleSpeciesMap::const_iterator it (m_particle_species.find(a_name));
            if (it == m_particle_species.end()) {
                amrex::Print()  << "ERROR: unable to find particle species with name \""
                                << a_name << "\"!";
                return nullptr;
            } else {
                return it->second;
            }
        }

        /*! accessor */
        inline const ERFPC* operator[] ( const std::string& a_name ) const
        {
            BL_PROFILE("ParticleData::operator[]");
            ParticleSpeciesMap::const_iterator it (m_particle_species.find(a_name));
            if (it == m_particle_species.end()) {
                amrex::Print()  << "ERROR: unable to find particle species with name \""
                                << a_name << "\"!";
                return nullptr;
            } else {
                return it->second;
            }
        }

        /*! Add a particle species to this container */
        inline void pushBack (const std::string&  a_name,
                             ERFPC* const        a_pc )
        {
            BL_PROFILE("ParticleData::pushBack()");
            AMREX_ASSERT(!contains(a_name));
            m_particle_species[a_name] = a_pc;
            m_namelist.push_back(a_name);
        }

        /*! Add a name; particle container will be initialized later */
        inline void addName (const std::string& a_name )
        {
            BL_PROFILE("ParticleData::addName()");
            m_namelist_unalloc.push_back(a_name);
        }

        /*! Returns list of names of particle species */
        inline const ParticlesNamesVector& getNames () const
        {
            BL_PROFILE("ParticleData::getNames()");
            return m_namelist;
        }

        /*! Returns list of names of particle species that are unallocated */
        inline ParticlesNamesList& getNamesUnalloc ()
        {
            BL_PROFILE("ParticleData::getNamesUnalloc()");
            return m_namelist_unalloc;
        }

        /*! queries if container has species of a certain name */
        inline bool contains ( const std::string& a_name ) const
        {
            BL_PROFILE("ParticleData::contains()");
            ParticleSpeciesMap::const_iterator it (m_particle_species.find(a_name));
            return (it != m_particle_species.end());
        }

        /*! query if container is empty */
        inline bool isEmpty () const
        {
            BL_PROFILE("ParticleData::isEmpty()");
            return (m_particle_species.size() == 0);
        }


    private:

        /*! Vector of all particle species */
        ParticleSpeciesMap m_particle_species;
        /*! Vector of particle species names */
        ParticlesNamesVector m_namelist;
        /*! List with names of unallocated species */
        ParticlesNamesList m_namelist_unalloc;

        /*! Disable particle output in plotfile? (because expensive) */
        bool m_disable_particle_op;
};

#endif
#endif
