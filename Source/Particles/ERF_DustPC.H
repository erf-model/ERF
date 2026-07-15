/**
 * @file ERF_DustPC.H
 * @brief Dust super-particle container for source-receptor attribution.
 *
 * Subclasses amrex::ParticleContainer with SoA layout:
 *   mass [kg], v_settle [m/s], release_time [s], src_i_f, src_j_f
 *
 * Advection: nearest-cell Euler step from ERF face velocities.
 * Settling: v_settle from ERF_DustSettling.H compute_stokes_settling.
 * Deposition: particle removed when z < z_lo + 0.5*dz.
 *             mass added to dust_source_map at (src_i, src_j).
 *
 * Pattern: ERFPCEvolve.cpp AdvectWithFlow on ERF-Mining.
 *
 * References:
 *   AMReX particles: https://amrex-codes.github.io/amrex/docs_html/Particles.html
 *   ERFPCEvolve.cpp: Source/Particles/ERFPCEvolve.cpp
 *   ERF_DustSettling.H: Source/Dust/ERF_DustSettling.H
 *   Shao (2008). Physics and Modelling of Wind Erosion. Springer.
 *   Seinfeld & Pandis (2006), Ch. 9.
 *     https://doi.org/10.1002/9781118227572
 */
#if defined(ERF_USE_DUST) && defined(ERF_USE_PARTICLES)
#ifndef ERF_DUST_PC_H_
#define ERF_DUST_PC_H_
#include <AMReX_Particles.H>
#include <ERF_DustSettling.H>

struct DustParticleRealIdx {
    enum { mass=0, v_settle, release_time, src_i_f, src_j_f, ncomps };
};

class ERFDustPC : public amrex::ParticleContainer<
    0, 0, DustParticleRealIdx::ncomps, 0>
{
public:
    ERFDustPC(const amrex::Geometry& g,
              const amrex::DistributionMapping& dm,
              const amrex::BoxArray& ba)
        : amrex::ParticleContainer<0,0,DustParticleRealIdx::ncomps,0>(g,dm,ba)
    {}

    /**
     * @brief Release one super-particle per active emission cell.
     * mass = flux * dx_dust * dy_dust * dt.
     * v_settle from compute_stokes_settling (ERF_DustSettling.H).
     * Position: cell centre x,y; z = z_lo + 0.5*dz_atm.
     */
    void ReleaseParticles(const amrex::MultiFab& emission_flux,
                          const amrex::Geometry& geom_atm,
                          const amrex::Geometry& geom_dust,
                          amrex::Real dt, amrex::Real d_m,
                          amrex::Real rho_p);

    /**
     * @brief Advance particles one timestep.
     * Nearest-cell interpolation of xvel/yvel/zvel.
     * z velocity adjusted by -v_settle.
     * Euler step: pos += vel * dt.
     * Deposit when z < z_lo + 0.5*dz: add mass to source_map(src_i,src_j).
     * Pattern: AdvectWithFlow in ERFPCEvolve.cpp.
     */
    void AdvanceParticles(const amrex::MultiFab& xvel,
                          const amrex::MultiFab& yvel,
                          const amrex::MultiFab& zvel,
                          amrex::MultiFab& source_map,
                          const amrex::Geometry& geom_atm,
                          const amrex::Geometry& geom_dust,
                          amrex::Real dt);

    amrex::Vector<std::string> varNames() const {
        return {"mass","v_settle","release_time","src_i","src_j"};
    }
};
#endif
#endif
