#include <stdio.h>
#include <limits>
#include <fstream>
#include "ERF_Constants.H"
#include "SuperDropletPC.H"

#ifdef ERF_USE_PARTICLES

using namespace amrex;

/*! Compute diagnostics (max, min, avg radius, mass, etc) */
void SuperDropletPC::Diagnostics( const int& a_iter,
                                  const bool a_flag )
{
    BL_PROFILE("SuperDropletPC::Diagnostics()");
    using PTDType = typename SuperDropletPC::ParticleTileType::ConstParticleTileDataType;

    auto num_total_particles = TotalNumberOfParticles();
    auto num_superdroplets = static_cast<Real>(NumSuperDroplets());

    auto min_par_radius = ReduceMin( *this,
                                     [=] AMREX_GPU_HOST_DEVICE (const PTDType& ptd, const int i) -> Real
                                     { return ptd.m_runtime_rdata[SuperDropletsRealIdxSoA_RT::radius][i]; } );

    auto max_par_radius = ReduceMax( *this,
                                     [=] AMREX_GPU_HOST_DEVICE (const PTDType& ptd, const int i) -> Real
                                     { return ptd.m_runtime_rdata[SuperDropletsRealIdxSoA_RT::radius][i]; } );

    auto avg_par_radius = ReduceSum( *this,
                                     [=] AMREX_GPU_HOST_DEVICE (const PTDType& ptd, const int i) -> Real
                                     {
                                         auto n = ptd.m_runtime_rdata[SuperDropletsRealIdxSoA_RT::multiplicity][i];
                                         auto r = ptd.m_runtime_rdata[SuperDropletsRealIdxSoA_RT::radius][i];
                                         return n*r;
                                     } );

    auto min_multiplic  = ReduceMin( *this,
                                     [=] AMREX_GPU_HOST_DEVICE (const PTDType& ptd, const int i) -> Real
                                     { return ptd.m_runtime_rdata[SuperDropletsRealIdxSoA_RT::multiplicity][i]; } );

    auto max_multiplic  = ReduceMax( *this,
                                     [=] AMREX_GPU_HOST_DEVICE (const PTDType& ptd, const int i) -> Real
                                     { return ptd.m_runtime_rdata[SuperDropletsRealIdxSoA_RT::multiplicity][i]; } );

    auto avg_multiplic  = ReduceSum( *this,
                                     [=] AMREX_GPU_HOST_DEVICE (const PTDType& ptd, const int i) -> Real
                                     { return ptd.m_runtime_rdata[SuperDropletsRealIdxSoA_RT::multiplicity][i]; } );

    auto min_par_mass   = ReduceMin( *this,
                                     [=] AMREX_GPU_HOST_DEVICE (const PTDType& ptd, const int i) -> Real
                                     { return ptd.m_rdata[SuperDropletsRealIdxSoA::mass][i]; } );

    auto max_par_mass   = ReduceMax( *this,
                                     [=] AMREX_GPU_HOST_DEVICE (const PTDType& ptd, const int i) -> Real
                                     { return ptd.m_rdata[SuperDropletsRealIdxSoA::mass][i]; } );

    auto avg_par_mass   = ReduceSum( *this,
                                     [=] AMREX_GPU_HOST_DEVICE (const PTDType& ptd, const int i) -> Real
                                     {
                                         auto n = ptd.m_runtime_rdata[SuperDropletsRealIdxSoA_RT::multiplicity][i];
                                         auto m = ptd.m_rdata[SuperDropletsRealIdxSoA::mass][i];
                                         return n*m;
                                     } );

    auto min_par_vx = ReduceMin( *this,
                                 [=] AMREX_GPU_HOST_DEVICE (const PTDType& ptd, const int i) -> Real
                                 { return ptd.m_rdata[SuperDropletsRealIdxSoA::vx][i]; } );

    auto max_par_vx = ReduceMax( *this,
                                 [=] AMREX_GPU_HOST_DEVICE (const PTDType& ptd, const int i) -> Real
                                 { return ptd.m_rdata[SuperDropletsRealIdxSoA::vx][i]; } );

    auto avg_par_vx = ReduceSum( *this,
                                 [=] AMREX_GPU_HOST_DEVICE (const PTDType& ptd, const int i) -> Real
                                 {
                                     auto n = ptd.m_runtime_rdata[SuperDropletsRealIdxSoA_RT::multiplicity][i];
                                     auto m = ptd.m_rdata[SuperDropletsRealIdxSoA::vx][i];
                                     return n*m;
                                 } );

    auto min_par_vy = ReduceMin( *this,
                                 [=] AMREX_GPU_HOST_DEVICE (const PTDType& ptd, const int i) -> Real
                                 { return ptd.m_rdata[SuperDropletsRealIdxSoA::vy][i]; } );

    auto max_par_vy = ReduceMax( *this,
                                 [=] AMREX_GPU_HOST_DEVICE (const PTDType& ptd, const int i) -> Real
                                 { return ptd.m_rdata[SuperDropletsRealIdxSoA::vy][i]; } );

    auto avg_par_vy = ReduceSum( *this,
                                 [=] AMREX_GPU_HOST_DEVICE (const PTDType& ptd, const int i) -> Real
                                 {
                                     auto n = ptd.m_runtime_rdata[SuperDropletsRealIdxSoA_RT::multiplicity][i];
                                     auto m = ptd.m_rdata[SuperDropletsRealIdxSoA::vy][i];
                                     return n*m;
                                 } );

    auto min_par_vz = ReduceMin( *this,
                                 [=] AMREX_GPU_HOST_DEVICE (const PTDType& ptd, const int i) -> Real
                                 { return ptd.m_rdata[SuperDropletsRealIdxSoA::vz][i]; } );

    auto max_par_vz = ReduceMax( *this,
                                 [=] AMREX_GPU_HOST_DEVICE (const PTDType& ptd, const int i) -> Real
                                 { return ptd.m_rdata[SuperDropletsRealIdxSoA::vz][i]; } );

    auto avg_par_vz = ReduceSum( *this,
                                 [=] AMREX_GPU_HOST_DEVICE (const PTDType& ptd, const int i) -> Real
                                 {
                                     auto n = ptd.m_runtime_rdata[SuperDropletsRealIdxSoA_RT::multiplicity][i];
                                     auto m = ptd.m_rdata[SuperDropletsRealIdxSoA::vz][i];
                                     return n*m;
                                 } );

    ParallelDescriptor::ReduceRealMin(&min_par_mass,1,ParallelDescriptor::IOProcessorNumber());
    ParallelDescriptor::ReduceRealMin(&min_par_radius,1,ParallelDescriptor::IOProcessorNumber());
    ParallelDescriptor::ReduceRealMin(&min_multiplic,1,ParallelDescriptor::IOProcessorNumber());
    ParallelDescriptor::ReduceRealMin(&min_par_vx,1,ParallelDescriptor::IOProcessorNumber());
    ParallelDescriptor::ReduceRealMin(&min_par_vy,1,ParallelDescriptor::IOProcessorNumber());
    ParallelDescriptor::ReduceRealMin(&min_par_vz,1,ParallelDescriptor::IOProcessorNumber());

    ParallelDescriptor::ReduceRealMax(&max_par_mass,1,ParallelDescriptor::IOProcessorNumber());
    ParallelDescriptor::ReduceRealMax(&max_par_radius,1,ParallelDescriptor::IOProcessorNumber());
    ParallelDescriptor::ReduceRealMax(&max_multiplic,1,ParallelDescriptor::IOProcessorNumber());
    ParallelDescriptor::ReduceRealMax(&max_par_vx,1,ParallelDescriptor::IOProcessorNumber());
    ParallelDescriptor::ReduceRealMax(&max_par_vy,1,ParallelDescriptor::IOProcessorNumber());
    ParallelDescriptor::ReduceRealMax(&max_par_vz,1,ParallelDescriptor::IOProcessorNumber());

    ParallelDescriptor::ReduceRealSum(&avg_par_mass,1,ParallelDescriptor::IOProcessorNumber());
    ParallelDescriptor::ReduceRealSum(&avg_par_radius,1,ParallelDescriptor::IOProcessorNumber());
    ParallelDescriptor::ReduceRealSum(&avg_multiplic,1,ParallelDescriptor::IOProcessorNumber());
    ParallelDescriptor::ReduceRealSum(&avg_par_vx,1,ParallelDescriptor::IOProcessorNumber());
    ParallelDescriptor::ReduceRealSum(&avg_par_vy,1,ParallelDescriptor::IOProcessorNumber());
    ParallelDescriptor::ReduceRealSum(&avg_par_vz,1,ParallelDescriptor::IOProcessorNumber());
    avg_par_mass /= num_total_particles;
    avg_par_radius /= num_total_particles;
    avg_multiplic /= num_superdroplets;
    avg_par_vx /= num_total_particles;
    avg_par_vy /= num_total_particles;
    avg_par_vz /= num_total_particles;

    Print() << "SuperDropletPC(" << m_name << ") attributes (min, max, avg):\n"
            << "    mass: "
            << min_par_mass << ", " << max_par_mass << ", " << avg_par_mass << "\n"
            << "    radius: "
            << min_par_radius << ", " << max_par_radius << ", " << avg_par_radius << "\n"
            << "    multiplicity: "
            << min_multiplic << ", " << max_multiplic << ", " << avg_multiplic << "\n"
            << "    velocity components:\n"
            << "        x: "
            << min_par_vx << ", " << max_par_vx << ", " << avg_par_vx << "\n"
            << "        y: "
            << min_par_vy << ", " << max_par_vy << ", " << avg_par_vy << "\n"
            << "        z: "
            << min_par_vz << ", " << max_par_vz << ", " << avg_par_vz << "\n";

    Long num_unconverged_particles = m_num_unconverged_particles;
    m_num_unconverged_particles = 0;
    ParallelDescriptor::ReduceLongSum(  &num_unconverged_particles,
                                        1,
                                        ParallelDescriptor::IOProcessorNumber() );
    int max_substeps_actual = m_max_substeps_actual;
    m_max_substeps_actual = 1;
    ParallelDescriptor::ReduceIntMax( &max_substeps_actual,
                                      1,
                                      ParallelDescriptor::IOProcessorNumber() );

    if (max_substeps_actual > 1) {
        Print() << "SuperDropletPC::MassChange(): max substeps = "
                << max_substeps_actual << "\n";
    }
    if (num_unconverged_particles > 0) {
        Print() << "SuperDropletPC::MassChange(): Warning - "
                << num_unconverged_particles
                << " particles did not converge during Newton solve.\n";
    }

    if (a_flag) {
        MassDensityDistribution( a_iter, min_par_radius, max_par_radius );
    }
}

/*! Compute and write the mass density distribution (as a function of the log of
    the droplet radius. The file written is a text file with two columns:
    R and g(ln R). */
void SuperDropletPC::MassDensityDistribution( const int& a_iter,
                                              const ParticleReal& a_r_min,
                                              const ParticleReal& a_r_max )
{
    int Nr = m_distribution_grid_size;

    const Geometry& geom = m_gdb->Geom(m_lev);
    const auto dxi = geom.InvCellSizeArray();

    const ParticleReal inv_cell_volume = dxi[0]*dxi[1]*dxi[2];
    const ParticleReal inv_bin_size
        = 1.0 / (  static_cast<ParticleReal>(m_coalescence_bin_size[0])
                 * static_cast<ParticleReal>(m_coalescence_bin_size[1])
                 * static_cast<ParticleReal>(m_coalescence_bin_size[2]) );
    const ParticleReal inv_bin_volume = inv_cell_volume*inv_bin_size;

    Vector<Real> ln_R, g_ln_R;
    ln_R.resize(Nr);
    g_ln_R.resize(Nr);

    // Set ln R grid
    for (int n = 0; n < Nr; n++) {
        ln_R[n] = std::log(a_r_min) + n*(std::log(a_r_max)-std::log(a_r_min))/(Nr-1);
    }

    const auto np = NumSuperDroplets();
    const ParticleReal sigma_0 = 0.62;
    const ParticleReal sigma = sigma_0 * std::exp(-0.2*std::log(static_cast<ParticleReal>(np)));
    const ParticleReal lambda = 1.0 / (2.0*sigma*sigma);
    const ParticleReal gamma = 1.0/(std::sqrt(2.0*PI)*sigma) * inv_bin_volume;

    // compute g(ln R)
    using PTDType = typename SuperDropletPC::ParticleTileType::ConstParticleTileDataType;
    for (int n = 0; n < Nr; n++) {
        const auto lnR = ln_R[n];
        g_ln_R[n] = ReduceSum(  *this,
                                [=] AMREX_GPU_HOST_DEVICE (const PTDType& ptd, const int i) -> Real
                                {
                                    auto ri = ptd.m_runtime_rdata[SuperDropletsRealIdxSoA_RT::radius][i];
                                    auto ni = ptd.m_runtime_rdata[SuperDropletsRealIdxSoA_RT::multiplicity][i];
                                    auto mi  = ptd.m_rdata[SuperDropletsRealIdxSoA::mass][i];

                                    auto lnRi = std::log(ri);
                                    return gamma*ni*mi*std::exp(-lambda*(lnR-lnRi)*(lnR-lnRi));
                                } );
    }

    // Sum g(ln R) over MPI subdomains
    ParallelDescriptor::ReduceRealSum(g_ln_R.dataPtr(),Nr);

    // Write to file
    char iter_str[12]; sprintf(iter_str, "%05d", a_iter+1);
    std::string output_filename =   m_name
                                    + "_mass_density_distribution_"
                                    + std::string(iter_str) + ".txt";
    Print() << "Writing " << output_filename << "\n";
    if (ParallelDescriptor::IOProcessor()) {
        std::ofstream outfile;
        outfile.open(output_filename.c_str(), std::ios::out|std::ios::trunc);
        if (!outfile.good()) { amrex::FileOpenFailed(output_filename); }

        for (int n = 0; n < Nr; n++) {
            outfile << std::exp(ln_R[n]) << " " << g_ln_R[n] << "\n";
        }

        outfile.flush();
        outfile.close();
        if (!outfile.good()) { amrex::Abort("problem writing output file"); }
    }
}

#endif

