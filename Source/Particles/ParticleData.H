#ifndef _PARTICLE_DATA_H_
#define _PARTICLE_DATA_H_

#include <string>
#include <iostream>

#include <IndexDefines.H>
#include <AMReX_Particles.H>
#include <AMReX_ParmParse.H>
#include <AMReX_Print.H>
#include <AMReX_Gpu.H>

#include <TracerPC.H>
#include <HydroPC.H>

/**
 * Container holding many of the particle-related data and options
 */

struct ParticleData {
  public:
    void init_particle_params()
    {
        amrex::ParmParse pp(pp_prefix);

        use_tracer_particles = 0;
        pp.query("use_tracer_particles", use_tracer_particles);

        use_hydro_particles = 0;
        pp.query("use_hydro_particles", use_hydro_particles);
    }

    void init_particles(amrex::ParGDBBase* a_gdb, const amrex::Vector<std::unique_ptr<amrex::MultiFab>>& z_phys_nd)

    {
        // Initialize tracer particles
        if (use_tracer_particles) {
            tracer_particles = std::make_unique<TracerPC>(a_gdb);
            if (z_phys_nd[0] != nullptr) {
                tracer_particles->InitParticles(*z_phys_nd[0]);
            } else {
                tracer_particles->InitParticles();
            }
            tracer_particles->Redistribute();
            amrex::Print() << "Initialized " << tracer_particles->TotalNumberOfParticles() << " tracer particles." << std::endl;
        }

        // Initialize hydro particles
        if (use_hydro_particles) {
            hydro_particles = std::make_unique<HydroPC>(a_gdb);
            if (z_phys_nd[0] != nullptr) {
                hydro_particles->InitParticles(*z_phys_nd[0]);
            } else {
                hydro_particles->InitParticles();
            }
            hydro_particles->Redistribute();
            amrex::Print() << "Initialized " << hydro_particles->TotalNumberOfParticles() << " hydro particles." << std::endl;
        }
    }

    void Checkpoint(const std::string& filename)
    {
        if (use_tracer_particles) {
            tracer_particles->Checkpoint(filename, "tracer_particles", true, tracer_particle_varnames);
        }
        if (use_hydro_particles) {
            hydro_particles->Checkpoint(filename, "hydro_particles", true, hydro_particle_varnames);
        }
    }

    void Restart(amrex::ParGDBBase* a_gdb, std::string& restart_file)
    {
        if (use_tracer_particles) {
            // tracer_particles = std::make_unique<TracerPC>(Geom(0), dmap[0], grids[0]);
            tracer_particles = std::make_unique<TracerPC>(a_gdb);
            std::string tracer_file("tracer_particles");
            tracer_particles->Restart(restart_file, tracer_file);
        }
        if (use_hydro_particles) {
            // hydro_particles = std::make_unique<HydroPC>(Geom(0), dmap[0], grids[0]);
            hydro_particles = std::make_unique<HydroPC>(a_gdb);
            std::string hydro_file("hydro_particles");
            hydro_particles->Restart(restart_file, hydro_file);
        }
    }

    void advance_particles(int lev, amrex::Real dt_lev,       amrex::Vector<amrex::Vector<amrex::MultiFab>>& vars_new,
                                                        const amrex::Vector<std::unique_ptr<amrex::MultiFab>>& z_phys_nd)
    {
        // Update tracer particles
        if (use_tracer_particles) {
            amrex::MultiFab* mac_vel = &vars_new[lev][Vars::xvel];
            bool use_terrain = (z_phys_nd[lev] != nullptr);
            if (use_terrain) {
                tracer_particles->AdvectWithUmac(mac_vel, lev, dt_lev, use_terrain, *z_phys_nd[lev]);
            } else {
                amrex::MultiFab dummy;
                tracer_particles->AdvectWithUmac(mac_vel, lev, dt_lev, use_terrain, dummy);
            }
        }

        // Update  hydro particles
        if (use_hydro_particles) {
            bool use_terrain = (z_phys_nd[lev] != nullptr);
            if (use_terrain) {
                hydro_particles->AdvectWithGravity(lev, dt_lev, use_terrain, *z_phys_nd[lev]);
            } else {
                amrex::MultiFab dummy;
                hydro_particles->AdvectWithGravity(lev, dt_lev, use_terrain, dummy);
            }
        }
    }

    void Redistribute()
    {
        if (use_tracer_particles) {
            tracer_particles->Redistribute();
        }
        if (use_hydro_particles) {
            hydro_particles->Redistribute();
        }
    }

    std::string pp_prefix {"erf"};

    int use_tracer_particles;
    std::unique_ptr<TracerPC> tracer_particles;
    amrex::Vector<std::string> tracer_particle_varnames = {AMREX_D_DECL("xvel", "yvel", "zvel")};

    int use_hydro_particles;
    std::unique_ptr<HydroPC> hydro_particles;
    amrex::Vector<std::string> hydro_particle_varnames = {AMREX_D_DECL("xvel", "yvel", "zvel"), "mass"};
};
#endif
