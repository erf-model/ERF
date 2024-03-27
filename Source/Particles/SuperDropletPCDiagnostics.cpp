#include <limits>
#include "SuperDropletPC.H"

#ifdef ERF_USE_PARTICLES

using namespace amrex;

/*! Compute diagnostics (max, min, avg radius, mass, etc) */
void SuperDropletPC::Diagnostics()
{
    BL_PROFILE("SuperDropletPC::Diagnostics()");

    Long num_total_particles = 0.0;
    Real avg_par_mass = 0.0;
    Real avg_par_radius = 0.0;
    Real max_par_mass = 0.0;
    Real max_par_radius = 0.0;
    Real min_par_mass = DBL_MAX;
    Real min_par_radius = DBL_MAX;

    for (ParIterType pti(*this, m_lev); pti.isValid(); ++pti) {

        auto& ptile = ParticlesAt(m_lev, pti);
        auto& aos = ptile.GetArrayOfStructs();
        auto& soa = ptile.GetStructOfArrays();
        const size_t np = aos.numParticles();
        auto* p_pbox = aos().data();

        /* SoA attributes */
        auto* mass_ptr = soa.GetRealData(SuperDropletsRealIdxSoA::mass).data();

        /* Runtime-added SoA attributes */
        int rt_offset = SuperDropletsRealIdxSoA::ncomps;
        auto* radius_ptr = soa.GetRealData(rt_offset+SuperDropletsRealIdxSoA_RT::radius).data();

        Gpu::Buffer<Real> avg_particle_mass({0}), avg_particle_radius({0});
        auto* avg_particle_mass_ptr = avg_particle_mass.data();
        auto* avg_particle_radius_ptr = avg_particle_radius.data();

        Gpu::Buffer<Real> min_particle_mass({DBL_MAX}), min_particle_radius({DBL_MAX});
        auto* min_particle_mass_ptr = min_particle_mass.data();
        auto* min_particle_radius_ptr = min_particle_radius.data();

        Gpu::Buffer<Real> max_particle_mass({0}), max_particle_radius({0});
        auto* max_particle_mass_ptr = max_particle_mass.data();
        auto* max_particle_radius_ptr = max_particle_radius.data();

        Gpu::Buffer<Long> num_particles({0});
        Long* num_particles_ptr = num_particles.data();

        ParallelFor(np, [=] AMREX_GPU_DEVICE (int i)
        {
            ParticleType& p = p_pbox[i];
            if (p.id() <= 0) { return; }

            auto radius = radius_ptr[i];
            auto mass = mass_ptr[i];

            Gpu::Atomic::Add( avg_particle_mass_ptr, mass );
            Gpu::Atomic::Add( avg_particle_radius_ptr, radius );

            Gpu::Atomic::Min( min_particle_mass_ptr, mass );
            Gpu::Atomic::Min( min_particle_radius_ptr, radius );

            Gpu::Atomic::Max( max_particle_mass_ptr, mass );
            Gpu::Atomic::Max( max_particle_radius_ptr, radius );

            Gpu::Atomic::Add( num_particles_ptr, Long(1) );

        });

        Gpu::streamSynchronize();

        max_par_mass = std::max( max_par_mass, *(max_particle_mass.copyToHost()) );
        max_par_radius = std::max( max_par_radius, *(max_particle_radius.copyToHost()) );

        min_par_mass = std::min( min_par_mass, *(min_particle_mass.copyToHost()) );
        min_par_radius = std::min( min_par_radius, *(min_particle_radius.copyToHost()) );

        avg_par_mass += *(avg_particle_mass.copyToHost());
        avg_par_radius += *(avg_particle_radius.copyToHost());

        num_total_particles += *(num_particles.copyToHost());
    }

    ParallelDescriptor::ReduceLongSum(&num_total_particles,1,ParallelDescriptor::IOProcessorNumber());

    ParallelDescriptor::ReduceRealMin(&min_par_mass,1,ParallelDescriptor::IOProcessorNumber());
    ParallelDescriptor::ReduceRealMin(&min_par_radius,1,ParallelDescriptor::IOProcessorNumber());

    ParallelDescriptor::ReduceRealMax(&max_par_mass,1,ParallelDescriptor::IOProcessorNumber());
    ParallelDescriptor::ReduceRealMax(&max_par_radius,1,ParallelDescriptor::IOProcessorNumber());

    ParallelDescriptor::ReduceRealSum(&avg_par_mass,1,ParallelDescriptor::IOProcessorNumber());
    ParallelDescriptor::ReduceRealSum(&avg_par_radius,1,ParallelDescriptor::IOProcessorNumber());
    avg_par_mass /= static_cast<Real>(num_total_particles);
    avg_par_radius /= static_cast<Real>(num_total_particles);

    Print() << "SuperDropletPC(" << m_name << ") diagnostics:\n"
            << "    particle mass (min,max,avg): "
            << min_par_mass << ", " << max_par_mass << ", " << avg_par_mass << "\n"
            << "    particle radius (min,max,avg): "
            << min_par_radius << ", " << max_par_radius << ", " << avg_par_radius << "\n";

}

#endif

