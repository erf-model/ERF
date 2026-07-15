/**
 * @file ERF_DustPC.cpp
 * @brief Implementation of ERFDustPC particle release and advection.
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

    const auto& dx_dust  = geom_dust.CellSize();
    const Real dx_dust_m = dx_dust[0];
    const Real dy_dust_m = dx_dust[1];

    // Use ATM geometry for all positions so particles are inside geom_atm domain
    const auto& plo_atm  = geom_atm.ProbLoArray();
    const auto& dx_atm   = geom_atm.CellSize();
    const Real dz_atm    = dx_atm[2];
    const Real z_lo_atm  = plo_atm[2];
    const Real z_release = z_lo_atm + 1.5 * dz_atm;

    // Also use atm dx for x,y so positions map into atm domain
    const Real dx_xy = dx_atm[0];
    const Real dy_xy = dx_atm[1];

    const Real v_settle = compute_stokes_settling(d_m, rho_p, 1.225,
                                                   DustSettlingConst::MU_AIR_STD);

    // Copy emission_flux onto the particle container's own BoxArray/DM
    amrex::MultiFab flux_on_pc(ParticleBoxArray(0),
                               ParticleDistributionMap(0),
                               emission_flux.nComp(),
                               amrex::IntVect(0));
    flux_on_pc.setVal(0.0);
    flux_on_pc.ParallelCopy(emission_flux, 0, 0, emission_flux.nComp(),
                            emission_flux.nGrowVect(),
                            amrex::IntVect(0),
                            geom_atm.periodicity());

    for (MFIter mfi(flux_on_pc, true); mfi.isValid(); ++mfi) {
        const Box& bx   = mfi.tilebox();
        auto flux_arr   = flux_on_pc.const_array(mfi);

        auto& particle_tile = DefineAndReturnParticleTile(0, mfi);
        auto& aos = particle_tile.GetArrayOfStructs();
        auto& soa = particle_tile.GetStructOfArrays();

        amrex::LoopOnCpu(bx, [&](int i, int j, int k) {
            if (k != 0 || flux_arr(i, j, k) <= 0.0) return;

            // Compute x,y from ATM geometry so position is inside geom_atm
            const Real x_pos = plo_atm[0] + (i + 0.5) * dx_xy;
            const Real y_pos = plo_atm[1] + (j + 0.5) * dy_xy;
            const Real mass  = flux_arr(i, j, k) * dx_dust_m * dy_dust_m * dt;

            ParticleType p;
            p.pos(0) = x_pos;
            p.pos(1) = y_pos;
            p.pos(2) = z_release;
            p.id()   = ParticleType::NextID();
            p.cpu()  = ParallelDescriptor::MyProc();
            aos().push_back(p);

            soa.GetRealData(DustParticleRealIdx::mass        ).push_back(mass);
            soa.GetRealData(DustParticleRealIdx::v_settle    ).push_back(v_settle);
            soa.GetRealData(DustParticleRealIdx::release_time).push_back(Real(0.0));
            soa.GetRealData(DustParticleRealIdx::src_i_f     ).push_back(Real(i));
            soa.GetRealData(DustParticleRealIdx::src_j_f     ).push_back(Real(j));
        });
    }

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

    const auto& plo_atm  = geom_atm.ProbLoArray();
    const auto& dx_atm   = geom_atm.CellSize();
    const Real z_lo      = plo_atm[2];
    const Real dz_atm    = dx_atm[2];
    const Real z_deposit = z_lo + 0.5 * dz_atm;

    const Box& dom_atm  = geom_atm.Domain();
    const int i_lo = dom_atm.smallEnd(0), i_hi = dom_atm.bigEnd(0);
    const int j_lo = dom_atm.smallEnd(1), j_hi = dom_atm.bigEnd(1);
    const int k_lo = dom_atm.smallEnd(2), k_hi = dom_atm.bigEnd(2);

    const Box& dom_dust = geom_dust.Domain();
    const int i_dust_lo = dom_dust.smallEnd(0), i_dust_hi = dom_dust.bigEnd(0);
    const int j_dust_lo = dom_dust.smallEnd(1), j_dust_hi = dom_dust.bigEnd(1);

    for (ParIterType pti(*this, 0); pti.isValid(); ++pti) {
        auto& particles     = pti.GetArrayOfStructs();
        const int npart     = particles.size();

        auto xvel_arr       = xvel.const_array(pti);
        auto yvel_arr       = yvel.const_array(pti);
        auto zvel_arr       = zvel.const_array(pti);
        auto source_map_arr = source_map.array(pti);

        auto& soa         = pti.GetStructOfArrays();
        auto& mass_vec    = soa.GetRealData(DustParticleRealIdx::mass);
        auto& vsett_vec   = soa.GetRealData(DustParticleRealIdx::v_settle);
        auto& srci_vec    = soa.GetRealData(DustParticleRealIdx::src_i_f);
        auto& srcj_vec    = soa.GetRealData(DustParticleRealIdx::src_j_f);

        for (int p = 0; p < npart; ++p) {
            ParticleType& particle = particles[p];
            if (particle.id() < 0) continue;

            Real x = particle.pos(0);
            Real y = particle.pos(1);
            Real z = particle.pos(2);

            const Real mass     = mass_vec[p];
            const Real v_settle = vsett_vec[p];
            const Real src_i_f  = srci_vec[p];
            const Real src_j_f  = srcj_vec[p];

            int i_cell = amrex::max(i_lo, amrex::min(i_hi,
                int(std::floor((x - plo_atm[0]) / dx_atm[0]))));
            int j_cell = amrex::max(j_lo, amrex::min(j_hi,
                int(std::floor((y - plo_atm[1]) / dx_atm[1]))));
            int k_cell = amrex::max(k_lo, amrex::min(k_hi,
                int(std::floor((z - plo_atm[2]) / dx_atm[2]))));

            Real u = 0.0, v_p = 0.0, w = 0.0;

            if (i_cell >= i_lo && i_cell < i_hi &&
                j_cell >= j_lo && j_cell < j_hi &&
                k_cell >= k_lo && k_cell < k_hi) {

                Real u1 = xvel_arr(i_cell,   j_cell, k_cell);
                Real u2 = (i_cell+1 <= i_hi) ? xvel_arr(i_cell+1, j_cell, k_cell) : u1;
                u   = 0.5*(u1+u2);

                Real v1 = yvel_arr(i_cell, j_cell,   k_cell);
                Real v2 = (j_cell+1 <= j_hi) ? yvel_arr(i_cell, j_cell+1, k_cell) : v1;
                v_p = 0.5*(v1+v2);

                Real w1 = zvel_arr(i_cell, j_cell, k_cell);
                Real w2 = (k_cell+1 <= k_hi) ? zvel_arr(i_cell, j_cell, k_cell+1) : w1;
                w   = 0.5*(w1+w2);
            }

            w -= v_settle;

            x += u   * dt;
            y += v_p * dt;
            z += w   * dt;

            if (z < z_deposit) {
                int src_i = amrex::max(i_dust_lo, amrex::min(i_dust_hi, int(src_i_f)));
                int src_j = amrex::max(j_dust_lo, amrex::min(j_dust_hi, int(src_j_f)));
                source_map_arr(src_i, src_j, 0) += mass;
                particle.id() = -1;
            } else {
                particle.pos(0) = x;
                particle.pos(1) = y;
                particle.pos(2) = z;
            }
        }
    }

    Redistribute();
}

#endif // ERF_USE_DUST && ERF_USE_PARTICLES