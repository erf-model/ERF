/**
 * @file ERF_DustPC.cpp
 * @brief Implementation of ERFDustPC particle release and advection.
 *
 * Phase 19: Lagrangian super-particles for dust source-receptor attribution.
 * Particles are released at dust grid cells with emission flux > 0,
 * advected by nearest-cell interpolated ERF face velocities with Stokes settling,
 * and deposited when they reach z < z_lo + 0.5*dz_atm.
 *
 * References:
 *   AMReX particles: https://amrex-codes.github.io/amrex/docs_html/Particles.html
 *   ERFPCEvolve.cpp pattern for AdvectWithFlow.
 *   Shao (2008). Physics and Modelling of Wind Erosion. Springer.
 *   Seinfeld & Pandis (2006) Ch. 9. Atmospheric Chemistry and Physics.
 */

#if defined(ERF_USE_DUST) && defined(ERF_USE_PARTICLES)

#include <ERF_DustPC.H>
#include <AMReX_Random.H>
#include <cmath>

using namespace amrex;

void ERFDustPC::ReleaseParticles(const amrex::MultiFab& emission_flux,
                                  const amrex::Geometry& geom_atm,
                                  const amrex::Geometry& geom_dust,
                                  amrex::Real dt, amrex::Real d_m,
                                  amrex::Real rho_p)
{
    BL_PROFILE("ERFDustPC::ReleaseParticles");

    // Get dust grid spacing for mass calculation
    const auto& dx_dust = geom_dust.CellSize();
    const Real dx_dust_m = dx_dust[0];
    const Real dy_dust_m = dx_dust[1];
    
    // Get atmospheric geometry for z position
    const auto& plo_atm = geom_atm.ProbLoArray();
    const auto& dx_atm = geom_atm.CellSize();
    const Real dz_atm = dx_atm[2];
    const Real z_lo_atm = plo_atm[2];
    
    // Release height: z_lo + 0.5*dz_atm
    const Real z_release = z_lo_atm + 0.5 * dz_atm;

    // Settling velocity
    const Real v_settle = compute_stokes_settling(d_m, rho_p, 1.225, 
                                                   DustSettlingConst::MU_AIR_STD);

    // Loop over dust grid and release particles where emission_flux > 0
    for (MFIter mfi(emission_flux, true); mfi.isValid(); ++mfi) {
        const Box& bx = mfi.tilebox();
        auto flux_arr = emission_flux.const_array(mfi);
        
        // Get dust grid positions for this tile
        const auto& dust_domain = geom_dust.Domain();
        const auto& plo_dust = geom_dust.ProbLoArray();
        
        // Iterate over dust cells in this tile
        amrex::LoopOnCpu(bx, [&](int i, int j, int k) {
            // k should be 0 for 2D dust grid
            if (k == 0 && flux_arr(i, j, k) > 0.0) {
                // Position: cell center in x,y; z_release in z
                const Real x_pos = plo_dust[0] + (i + 0.5) * dx_dust_m;
                const Real y_pos = plo_dust[1] + (j + 0.5) * dy_dust_m;
                
                // Mass per particle: flux [kg/m^2/s] * area [m^2] * dt [s]
                const Real mass = flux_arr(i, j, k) * dx_dust_m * dy_dust_m * dt;
                
                // Create particle using DefineAndReturnParticleTile
                auto& particle_tile = DefineAndReturnParticleTile(
                    0, 0, mfi);
                
                ParticleType p;
                p.pos(0) = x_pos;
                p.pos(1) = y_pos;
                p.pos(2) = z_release;
                p.id() = amrex::ParticleBase::NextID();
                p.cpu() = ParallelDescriptor::MyProc();
                
                // Real attributes: mass, v_settle, release_time, src_i_f, src_j_f
                p.rdata(DustParticleRealIdx::mass) = mass;
                p.rdata(DustParticleRealIdx::v_settle) = v_settle;
                p.rdata(DustParticleRealIdx::release_time) = 0.0;  // Current time
                p.rdata(DustParticleRealIdx::src_i_f) = Real(i);
                p.rdata(DustParticleRealIdx::src_j_f) = Real(j);
                
                particle_tile.push_back(p);
            }
        });
    }

    // Redistribute particles to correct processors
    Redistribute();
}

void ERFDustPC::AdvanceParticles(const amrex::MultiFab& xvel,
                                 const amrex::MultiFab& yvel,
                                 const amrex::MultiFab& zvel,
                                 amrex::MultiFab& source_map,
                                 const amrex::Geometry& geom_atm,
                                 const amrex::Geometry& geom_dust,
                                 amrex::Real dt)
{
    BL_PROFILE("ERFDustPC::AdvanceParticles");

    const auto& plo_atm = geom_atm.ProbLoArray();
    const auto& dx_atm = geom_atm.CellSize();
    
    const Real z_lo = plo_atm[2];
    const Real dz_atm = dx_atm[2];
    const Real z_deposit = z_lo + 0.5 * dz_atm;
    
    // Get atmospheric domain for nearest-cell clamping
    const Box& dom_atm = geom_atm.Domain();
    const int i_lo = dom_atm.smallEnd(0), i_hi = dom_atm.bigEnd(0);
    const int j_lo = dom_atm.smallEnd(1), j_hi = dom_atm.bigEnd(1);
    const int k_lo = dom_atm.smallEnd(2), k_hi = dom_atm.bigEnd(2);

    // Iterate over particles
    for (ParIterType pti(*this, 0); pti.isValid(); ++pti) {
        auto& particles = pti.GetArrayOfStructs();
        int npart = particles.size();
        
        // Get velocity arrays for this tile
        auto xvel_const_arr = xvel.const_array(pti);
        auto yvel_const_arr = yvel.const_array(pti);
        auto zvel_const_arr = zvel.const_array(pti);
        auto source_map_arr = source_map.array(pti);

        for (int p = 0; p < npart; ++p) {
            ParticleType& particle = particles[p];
            
            // Skip dead particles
            if (particle.id() < 0) continue;
            
            // Get particle position
            Real x = particle.pos(0);
            Real y = particle.pos(1);
            Real z = particle.pos(2);
            
            // Get SoA attributes
            auto& soa = pti.GetStructOfArrays();
            auto mass_arr = soa.GetRealData(DustParticleRealIdx::mass);
            auto v_settle_arr = soa.GetRealData(DustParticleRealIdx::v_settle);
            auto src_i_arr = soa.GetRealData(DustParticleRealIdx::src_i_f);
            auto src_j_arr = soa.GetRealData(DustParticleRealIdx::src_j_f);
            
            const Real mass = mass_arr[p];
            const Real v_settle = v_settle_arr[p];
            const Real src_i_f = src_i_arr[p];
            const Real src_j_f = src_j_arr[p];
            
            // Find nearest atmospheric grid cell
            int i_cell = int(std::floor((x - plo_atm[0]) / dx_atm[0]));
            int j_cell = int(std::floor((y - plo_atm[1]) / dx_atm[1]));
            int k_cell = int(std::floor((z - plo_atm[2]) / dx_atm[2]));
            
            // Clamp to domain
            i_cell = amrex::max(i_lo, amrex::min(i_hi, i_cell));
            j_cell = amrex::max(j_lo, amrex::min(j_hi, j_cell));
            k_cell = amrex::max(k_lo, amrex::min(k_hi, k_cell));
            
            // Get velocities at this cell (averaged from faces)
            Real u = 0.0, v = 0.0, w = 0.0;
            
            // Check bounds before accessing
            if (i_cell >= i_lo && i_cell < i_hi && j_cell >= j_lo && j_cell < j_hi && 
                k_cell >= k_lo && k_cell < k_hi) {
                // Average x-velocities
                Real u1 = xvel_const_arr(i_cell, j_cell, k_cell);
                Real u2 = (i_cell + 1 < i_hi) ? xvel_const_arr(i_cell+1, j_cell, k_cell) : u1;
                u = 0.5 * (u1 + u2);
                
                // Average y-velocities
                Real v1 = yvel_const_arr(i_cell, j_cell, k_cell);
                Real v2 = (j_cell + 1 < j_hi) ? yvel_const_arr(i_cell, j_cell+1, k_cell) : v1;
                v = 0.5 * (v1 + v2);
                
                // Average z-velocities
                Real w1 = zvel_const_arr(i_cell, j_cell, k_cell);
                Real w2 = (k_cell + 1 < k_hi) ? zvel_const_arr(i_cell, j_cell, k_cell+1) : w1;
                w = 0.5 * (w1 + w2);
            }
            
            // Apply settling: reduce vertical velocity
            w -= v_settle;
            
            // Euler step: update position
            x += u * dt;
            y += v * dt;
            z += w * dt;
            
            // Check deposition condition
            if (z < z_deposit) {
                // Particle has deposited
                // Add mass to source_map at source cell (src_i, src_j)
                int src_i = int(src_i_f);
                int src_j = int(src_j_f);
                
                // Clamp to dust grid domain
                const Box& dom_dust = geom_dust.Domain();
                const int i_dust_lo = dom_dust.smallEnd(0);
                const int i_dust_hi = dom_dust.bigEnd(0);
                const int j_dust_lo = dom_dust.smallEnd(1);
                const int j_dust_hi = dom_dust.bigEnd(1);
                
                src_i = amrex::max(i_dust_lo, amrex::min(i_dust_hi, src_i));
                src_j = amrex::max(j_dust_lo, amrex::min(j_dust_hi, src_j));
                
                // Add mass to source_map (CPU: direct add, not GPU atomic for now)
                source_map_arr(src_i, src_j, 0) += mass;
                
                // Mark particle as dead
                particle.id() = -1;
            } else {
                // Update particle position
                particle.pos(0) = x;
                particle.pos(1) = y;
                particle.pos(2) = z;
            }
        }
    }

    // Redistribute to remove dead particles and redistribute across processors
    Redistribute();
}

#endif // ERF_USE_DUST && ERF_USE_PARTICLES
