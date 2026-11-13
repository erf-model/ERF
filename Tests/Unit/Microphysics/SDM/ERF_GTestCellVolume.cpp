#include <gtest/gtest.h>
#include <cmath>
#include <limits>
#include <map>
#include <memory>

#include <ERF_SuperDropletPC.H>

using namespace amrex;

// A column of 2 x 2 x 4 cells with computational spacing 1 m: two 0.5 m layers
// under two 1.5 m layers. Particle counts and densities must follow the physical
// cell volumes, not the computational volume of 1 m^3.
// Keep device lambdas in free/public functions, not GoogleTest's private TestBody (NVCC).

struct SDMCellTotal
{
    int count = 0;
    ParticleReal mult = 0;
    ParticleReal condensate = 0;
};

const Box sdm_cv_domain (IntVect(0), IntVect(AMREX_D_DECL(1,1,3)));

Geometry sdm_cv_geom ()
{
    const RealBox rb(Real(0), Real(0), Real(0), Real(2), Real(2), Real(4));
    return Geometry(sdm_cv_domain, rb, CoordSys::cartesian, {1,1,0});
}

Real sdm_cv_dz (int k) { return (k < 2) ? Real(0.5) : Real(1.5); }

// Node heights: two 0.5 m layers under two layers of a_top_dz
Vector<std::unique_ptr<MultiFab>> sdm_cv_heights (const BoxArray& ba, const DistributionMapping& dm,
                                                  Real a_top_dz = Real(1.5))
{
    BoxArray nodes(ba);
    nodes.surroundingNodes();
    Vector<std::unique_ptr<MultiFab>> z_nd(1);
    z_nd[0] = std::make_unique<MultiFab>(nodes, dm, 1, 2);
    for (MFIter mfi(*z_nd[0]); mfi.isValid(); ++mfi) {
        const auto h = z_nd[0]->array(mfi);
        ParallelFor(mfi.fabbox(), [=] AMREX_GPU_DEVICE (int i, int j, int k) noexcept {
            amrex::ignore_unused(i, j);
            h(i,j,k) = (k <= 2) ? Real(0.5)*Real(k) : Real(1) + a_top_dz*Real(k-2);
        });
    }
    return z_nd;
}

// Number, multiplicity and condensate (mass less the aerosol mass a_dry) per cell
std::map<int, SDMCellTotal> sdm_cv_totals (SuperDropletPC& pc, ParticleReal a_dry = 0)
{
    std::map<int, SDMCellTotal> totals;
    for (SuperDropletPC::ParIterType pti(pc, 0); pti.isValid(); ++pti) {
        const auto& tile = pti.GetParticleTile();
        const int np = tile.numParticles();
        Vector<SuperDropletPC::ParticleType> particles(np);
        Gpu::copy(Gpu::deviceToHost, tile.GetArrayOfStructs().dataPtr(),
                  tile.GetArrayOfStructs().dataPtr()+np, particles.begin());
        const auto& soa = tile.GetStructOfArrays();
        const auto& xi = soa.GetRealData(static_cast<int>(SDPCDefn::SuperDropletsRealIdx::ncomps)
                                        + SDPCDefn::SuperDropletsRealIdxSoA_RT::multiplicity);
        const auto& m = soa.GetRealData(SDPCDefn::SuperDropletsRealIdx::mass);
        Vector<ParticleReal> mult(np), mass(np);
        Gpu::copy(Gpu::deviceToHost, xi.begin(), xi.end(), mult.begin());
        Gpu::copy(Gpu::deviceToHost, m.begin(), m.end(), mass.begin());
        for (int n = 0; n < np; ++n) {
            if (particles[n].id() <= 0) { continue; }
            const int i = static_cast<int>(std::floor(particles[n].pos(0)));
            const int j = static_cast<int>(std::floor(particles[n].pos(1)));
            const int k = static_cast<int>(std::floor(particles[n].pos(2)));
            auto& t = totals[i + 2*(j + 2*k)];
            t.count += 1;
            t.mult += mult[n];
            t.condensate += mult[n] * (mass[n] - a_dry);
        }
    }
    return totals;
}

void sdm_cv_set (SuperDropletPC& pc, int a_active, ParticleReal a_mult)
{
    for (SuperDropletPC::ParIterType pti(pc, 0); pti.isValid(); ++pti) {
        const auto data = pti.GetParticleTile().getParticleTileData();
        ParallelFor(pti.numParticles(), [=] AMREX_GPU_DEVICE (int i) noexcept {
            if (a_active >= 0) {
                data.m_runtime_idata[SDPCDefn::SuperDropletsIntIdxSoA_RT::active][i] = a_active;
            }
            if (a_mult != ParticleReal(0)) {
                data.m_runtime_rdata[SDPCDefn::SuperDropletsRealIdxSoA_RT::multiplicity][i] = a_mult;
            }
        });
    }
}

void sdm_cv_initialization_test ()
{
    const Geometry geom = sdm_cv_geom();
    const BoxArray ba(sdm_cv_domain);
    const DistributionMapping dm(ba);
    auto z_nd = sdm_cv_heights(ba, dm);

    const std::string name = "sdm_cell_volume_init";
    ParmParse pp(name);
    pp.add("num_initializations", 1);
    pp.add("num_injections", 0);
    pp.add("initial_number_density", Real(100));
    pp.add("initial_particles_per_cell", 4);
    pp.add("multiplicity_type", "constant");
    pp.add("initial_species_mean_mass_H2O", Real(1.e-12));
    pp.add("initial_aerosol_mean_mass_NaCl", Real(1.e-19));
    pp.add("place_randomly_in_cells", false);
    SuperDropletPC pc(geom, dm, ba, {Species::Name::H2O}, {Species::Name::NaCl}, 1.0, name);
    pc.setVariableDz(true);
    pc.InitializeParticles(0, 0, z_nd[0]);

    // without ice the ice attributes are defined, and zero
    for (SuperDropletPC::ParIterType pti(pc, 0); pti.isValid(); ++pti) {
        const auto& soa = pti.GetStructOfArrays();
        for (int c : {idx_ice_Tfz(1,1), idx_ice_a(1,1), idx_ice_c(1,1), idx_ice_mrime(1,1), idx_ice_nmono(1,1)}) {
            const auto& attr = soa.GetRealData(c);
            Vector<ParticleReal> values(attr.size());
            Gpu::copy(Gpu::deviceToHost, attr.begin(), attr.end(), values.begin());
            for (auto v : values) { EXPECT_EQ(v, ParticleReal(0)); }
        }
    }

    // ceil(100 m^-3 x cell volume) particles per cell, on 4 super-droplets
    auto totals = sdm_cv_totals(pc);
    ASSERT_EQ(static_cast<int>(totals.size()), 16);
    for (const auto& [cell, t] : totals) {
        const Real vol = sdm_cv_dz(cell/4);
        EXPECT_EQ(t.count, 4);
        EXPECT_EQ(t.mult, ParticleReal(100*vol));
    }

    // the condensate of each cell is its density times its volume
    MultiFab condensate(ba, dm, 1, 0);
    condensate.setVal(Real(1.e-3));
    pc.SetAttributes(condensate, *z_nd[0], 0);
    totals = sdm_cv_totals(pc, ParticleReal(1.e-19));
    for (const auto& [cell, t] : totals) {
        const Real expected = Real(1.e-3) * sdm_cv_dz(cell/4);
        EXPECT_NEAR(t.condensate, expected,
                    Real(64)*Real(std::numeric_limits<ParticleReal>::epsilon())*expected);
    }

    // super-droplets that entered through an open boundary take the initial multiplicity of their cell
    sdm_cv_set(pc, -1, ParticleReal(-1));
    pc.setInflowMultiplicity(z_nd);
    totals = sdm_cv_totals(pc);
    for (const auto& [cell, t] : totals) {
        EXPECT_EQ(t.mult, t.count * ParticleReal(100*sdm_cv_dz(cell/4)/4));
    }
}

void sdm_cv_zero_density_test ()
{
    const Geometry geom = sdm_cv_geom();
    const BoxArray ba(sdm_cv_domain);
    const DistributionMapping dm(ba);
    auto z_nd = sdm_cv_heights(ba, dm);

    const std::string name = "sdm_cell_volume_zero";
    ParmParse pp(name);
    pp.add("num_initializations", 1);
    pp.add("num_injections", 0);
    pp.add("initial_number_density", Real(0));
    pp.add("initial_particles_per_cell", 2);
    pp.add("multiplicity_type", "constant");
    pp.add("initial_species_mean_mass_H2O", Real(1.e-12));
    pp.add("initial_aerosol_mean_mass_NaCl", Real(1.e-19));
    SuperDropletPC pc(geom, dm, ba, {Species::Name::H2O}, {Species::Name::NaCl}, 1.0, name);
    pc.setVariableDz(true);
    pc.InitializeParticles(0, 0, z_nd[0]);
    EXPECT_EQ(pc.TotalNumberOfParticles(), Real(0));
}

void sdm_cv_created_counts_test ()
{
    const Geometry geom = sdm_cv_geom();
    const BoxArray ba(sdm_cv_domain);
    const DistributionMapping dm(ba);
    auto z_nd = sdm_cv_heights(ba, dm);

    // a second initialization without particles creates no super-droplets, so the
    // condensate of a cell is split over those of the first one only
    {
        const std::string name = "sdm_cell_volume_two_inits";
        ParmParse pp(name);
        pp.add("num_initializations", 2);
        pp.add("num_injections", 0);
        pp.add("initial_number_density", Real(100));
        pp.add("initial_particles_per_cell", 4);
        pp.add("multiplicity_type", "constant");
        pp.add("initial_species_mean_mass_H2O", Real(1.e-12));
        pp.add("initial_aerosol_mean_mass_NaCl", Real(1.e-19));
        pp.add("place_randomly_in_cells", false);
        ParmParse pp1(name + ".1");
        pp1.add("initial_number_density", Real(0));
        pp1.add("initial_particles_per_cell", 2);
        SuperDropletPC pc(geom, dm, ba, {Species::Name::H2O}, {Species::Name::NaCl}, 1.0, name);
        pc.setVariableDz(true);
        pc.InitializeParticles(0, 0, z_nd[0]);

        MultiFab condensate(ba, dm, 1, 0);
        condensate.setVal(Real(1.e-3));
        pc.SetAttributes(condensate, *z_nd[0], 0);
        const auto totals = sdm_cv_totals(pc, ParticleReal(1.e-19));
        ASSERT_EQ(static_cast<int>(totals.size()), 16);
        for (const auto& [cell, t] : totals) {
            const Real expected = Real(1.e-3) * sdm_cv_dz(cell/4);
            EXPECT_EQ(t.count, 4);
            EXPECT_NEAR(t.condensate, expected,
                        Real(64)*Real(std::numeric_limits<ParticleReal>::epsilon())*expected);
        }
    }

    // fewer particles than super-droplets per cell: one super-droplet per particle
    {
        const std::string name = "sdm_cell_volume_sparse";
        ParmParse pp(name);
        pp.add("num_initializations", 1);
        pp.add("num_injections", 0);
        pp.add("initial_number_density", Real(4));
        pp.add("initial_particles_per_cell", 4);
        pp.add("multiplicity_type", "constant");
        pp.add("initial_species_mean_mass_H2O", Real(1.e-12));
        pp.add("initial_aerosol_mean_mass_NaCl", Real(1.e-19));
        pp.add("place_randomly_in_cells", false);
        SuperDropletPC pc(geom, dm, ba, {Species::Name::H2O}, {Species::Name::NaCl}, 1.0, name);
        pc.setVariableDz(true);
        pc.InitializeParticles(0, 0, z_nd[0]);
        const auto totals = sdm_cv_totals(pc);
        ASSERT_EQ(static_cast<int>(totals.size()), 16);
        for (const auto& [cell, t] : totals) {
            const int num_par = static_cast<int>(std::ceil(Real(4)*sdm_cv_dz(cell/4)));  // 2 or 6
            EXPECT_EQ(t.count, std::min(num_par, 4));
            EXPECT_EQ(t.mult, ParticleReal(num_par));
        }
    }
}

void sdm_cv_recycle_test ()
{
    // columns of 3 m (mean cell volume 0.75 m^3) under a computational spacing of 1 m
    const Geometry geom = sdm_cv_geom();
    const BoxArray ba(sdm_cv_domain);
    const DistributionMapping dm(ba);
    auto z_nd = sdm_cv_heights(ba, dm, Real(1));

    const std::string name = "sdm_cell_volume_recycle";
    ParmParse pp(name);
    pp.add("num_initializations", 1);
    pp.add("num_injections", 0);
    pp.add("initial_number_density", Real(100));
    pp.add("initial_particles_per_cell", 4);
    pp.add("multiplicity_type", "constant");
    pp.add("initial_species_mean_mass_H2O", Real(1.e-12));
    pp.add("initial_aerosol_mean_mass_NaCl", Real(1.e-19));
    pp.add("recycle_zmax", Real(3));
    SuperDropletPC pc(geom, dm, ba, {Species::Name::H2O}, {Species::Name::NaCl}, 1.0, name);
    pc.setVariableDz(true);
    pc.InitializeParticles(0, 0, z_nd[0]);

    // recycled super-droplets land uniformly in physical space and all take the multiplicity
    // of the mean cell volume, 100 m^-3 x 0.75 m^3 / 4, so the number density stays uniform
    sdm_cv_set(pc, 0, 0);
    pc.Recycle(0, z_nd, 0, 1.0, true);
    int recycled = 0;
    for (SuperDropletPC::ParIterType pti(pc, 0); pti.isValid(); ++pti) {
        const auto& soa = pti.GetStructOfArrays();
        const auto& xi = soa.GetRealData(static_cast<int>(SDPCDefn::SuperDropletsRealIdx::ncomps)
                                        + SDPCDefn::SuperDropletsRealIdxSoA_RT::multiplicity);
        const auto& act = soa.GetIntData(static_cast<int>(SDPCDefn::SuperDropletsIntIdx::ncomps)
                                        + SDPCDefn::SuperDropletsIntIdxSoA_RT::active);
        Vector<ParticleReal> mult(xi.size());
        Vector<int> active(act.size());
        Gpu::copy(Gpu::deviceToHost, xi.begin(), xi.end(), mult.begin());
        Gpu::copy(Gpu::deviceToHost, act.begin(), act.end(), active.begin());
        for (int n = 0; n < static_cast<int>(mult.size()); ++n) {
            EXPECT_EQ(active[n], 1);
            EXPECT_EQ(mult[n], ParticleReal(18.75));
            recycled++;
        }
    }
    EXPECT_EQ(recycled, 64);
}

void sdm_cv_injection_test (bool a_perbox)
{
    const Geometry geom = sdm_cv_geom();
    const BoxArray ba(sdm_cv_domain);
    const DistributionMapping dm(ba);
    auto z_nd = sdm_cv_heights(ba, dm);

    const std::string name = a_perbox ? "sdm_cell_volume_perbox" : "sdm_cell_volume_inject";
    ParmParse pp(name);
    pp.add("num_initializations", 0);
    pp.add("initial_number_density", Real(0));
    pp.add("initial_particles_per_cell", 0);
    pp.add("num_injections", 1);
    pp.add("injection.distribution_type", "uniform");
    pp.addarr("injection.particle_box_lo", std::vector<Real>{0, 0, 0});
    pp.addarr("injection.particle_box_hi", std::vector<Real>{2, 2, 4});
    pp.add("injection.multiplicity_type", "constant");
    pp.add("injection.aerosol_mean_mass_NaCl", Real(1.e-19));
    if (a_perbox) {
        // 2500 m^-3 s^-1 super-droplets at multiplicity 10, scattered over the box
        pp.add("injection.rate", Real(25000));
        pp.add("injection.sd_rate", Real(2500));
    } else {
        // 1 m^-3 s^-1 over steps of 1, 0.5, 1.5 and 1 s: 0.25 to 2.25 particles per cell and
        // step, on a number of super-droplets per cell that changes with the step
        pp.add("injection.rate", Real(1));
        pp.add("injection.sd_rate", Real(2));
        pp.add("injection.particles_per_cell", 2);
    }
    SuperDropletPC pc(geom, dm, ba, {Species::Name::H2O}, {Species::Name::NaCl}, 1.0, name);
    pc.setVariableDz(true);

    if (a_perbox) {
        pc.InjectParticles(0.0, z_nd[0], 1.0);
        auto totals = sdm_cv_totals(pc);
        Vector<Real> layer(4, 0);
        for (const auto& [cell, t] : totals) {
            EXPECT_EQ(t.mult, ParticleReal(10 * t.count));
            layer[cell/4] += Real(t.count);
        }
        // 40000 super-droplets, in proportion to the layer volumes 2, 2, 6, 6 m^3
        for (int k = 0; k < 4; ++k) {
            const Real expected = Real(40000) * Real(4)*sdm_cv_dz(k) / Real(16);
            EXPECT_NEAR(layer[k], expected, Real(0.05)*expected);
        }
    } else {
        Real time = 0;
        for (Real dt : {Real(1), Real(0.5), Real(1.5), Real(1)}) {
            pc.InjectParticles(time, z_nd[0], dt);
            time += dt;
        }
        // 4 m^-3 after 4 s in every cell, on whole super-droplets of unit multiplicity
        auto totals = sdm_cv_totals(pc);
        ASSERT_EQ(static_cast<int>(totals.size()), 16);
        for (const auto& [cell, t] : totals) {
            EXPECT_EQ(t.mult, ParticleReal(4*sdm_cv_dz(cell/4)));
            EXPECT_EQ(t.mult, ParticleReal(t.count));
        }
    }
}

// Motivation: On a stretched mesh the initial super-droplets, their condensate, and the
// multiplicity of super-droplets entering through an open boundary must follow each
// cell's physical volume; with the computational volume the thin layers hold twice and
// the thick layers two thirds of the specified number density.
TEST(SDMCellVolume, InitializationFollowsPhysicalCellVolume)
{
    sdm_cv_initialization_test();
}

// Motivation: A zero number density must add no super-droplets on a stretched mesh,
// as on a uniform one, whatever the number of super-droplets per cell.
TEST(SDMCellVolume, ZeroDensityAddsNothing)
{
    sdm_cv_zero_density_test();
}

// Motivation: A cell's condensate must be split over the super-droplets that exist in it,
// not over those of an initialization without particles, and no super-droplet may get a
// zero multiplicity where a cell has fewer particles than super-droplets.
TEST(SDMCellVolume, CountsFollowCreatedSuperDroplets)
{
    sdm_cv_created_counts_test();
}

// Motivation: Recycled super-droplets are placed uniformly in physical space, so they must
// all take one multiplicity, that of the mean physical cell volume; per-cell multiplicities
// would make the number density follow the cell volume, and the computational volume
// overstates it where the columns are shorter than the computational height.
TEST(SDMCellVolume, RecyclingUsesMeanPhysicalCellVolume)
{
    sdm_cv_recycle_test();
}

// Motivation: Per-cell injection must deliver the specified number per unit physical
// volume in thin and thick cells alike, including steps that add less than one particle
// per super-droplet and a time step that changes the number of super-droplets per cell.
TEST(SDMCellVolume, PerCellInjectionFollowsPhysicalCellVolume)
{
    sdm_cv_injection_test(false);
}

// Motivation: Per-box injection must bin its scattered super-droplets by physical height,
// so each layer receives a share in proportion to its physical volume.
TEST(SDMCellVolume, PerBoxInjectionFollowsPhysicalCellVolume)
{
    sdm_cv_injection_test(true);
}
