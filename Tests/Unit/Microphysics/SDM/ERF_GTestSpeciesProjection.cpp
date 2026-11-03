#include <gtest/gtest.h>
#include <limits>
#include <map>
#include <AMReX_BoxIterator.H>

#include <ERFPCParticleToMesh.H>
#include <AMReX_MultiFabUtil.H>
#include <AMReX_ParGDB.H>
#include <ERF_SuperDropletsMoist.H>

using namespace amrex;

// Keep device lambdas in free/public functions, not GoogleTest's private TestBody (NVCC).
void composite_cic_mass_test ()
{
    const Box coarse_domain(IntVect(0), IntVect(3));
    const RealBox rb(Real(0), Real(0), Real(0), Real(4), Real(4), Real(4));
    const Vector<Geometry> geom{
        Geometry(coarse_domain, rb, CoordSys::cartesian, {1, 0, 0}),
        Geometry(refine(coarse_domain, 2), rb, CoordSys::cartesian, {1, 0, 0})};
    Box fine_box = geom[1].Domain();
    fine_box.setBig(2, 3);
    const Vector<BoxArray> ba{BoxArray(coarse_domain), BoxArray(fine_box)};
    const Vector<DistributionMapping> dm{DistributionMapping(ba[0]), DistributionMapping(ba[1])};
    const Vector<IntVect> rr{IntVect(2)};
    ParGDB gdb(geom, dm, ba, rr);

    // Fine wall, fine-side interface, coarse-side interface; x also crosses
    // the periodic seam. No advection, microphysics, or plotting can repair this deposit.
    for (Real zpos : {Real(0.1), Real(1.9), Real(2.1)}) {
        ERFPC pc(&gdb, "species_projection_test", ERFPCOptions{false, false});
        pc.InitRandom(1, 42, {});
        for (int lev = 0; lev <= pc.finestLevel(); ++lev) {
            for (ERFPC::ParIterType pti(pc, lev); pti.isValid(); ++pti) {
                auto* particles = pti.GetArrayOfStructs().dataPtr();
                ParallelFor(pti.numParticles(), [=] AMREX_GPU_DEVICE (int i) noexcept {
                    particles[i].pos(0) = ParticleReal(0.1);
                    particles[i].pos(1) = ParticleReal(1.5);
                    particles[i].pos(2) = ParticleReal(zpos);
                });
            }
        }
        pc.Redistribute();
        Real mass = Real(0);
        for (int lev = 0; lev < 2; ++lev) {
            MultiFab deposit(ba[lev], dm[lev], 1, 1);
            BoxArray nodes(ba[lev]);
            nodes.surroundingNodes();
            MultiFab height(nodes, dm[lev], 1, 2);
            const Real dz = geom[lev].CellSize(2);
            for (MFIter mfi(height); mfi.isValid(); ++mfi) {
                const auto h = height.array(mfi);
                ParallelFor(mfi.fabbox(), [=] AMREX_GPU_DEVICE (int i, int j, int k) noexcept {
                    h(i,j,k) = Real(k) * dz;
                });
            }
            pc.ERFPCParticleToMesh(deposit, height, lev, 0,
                [=] AMREX_GPU_DEVICE (const ERFPC::ParticleTileType::ConstParticleTileDataType&, int) {
                    return ParticleReal(1);
                }, true);
            if (lev == 0) {
                const iMultiFab mask = makeFineMask(deposit, ba[1], rr[0], 1, 0);
                for (MFIter mfi(deposit); mfi.isValid(); ++mfi) {
                    const auto d = deposit.array(mfi);
                    const auto m = mask.const_array(mfi);
                    ParallelFor(mfi.validbox(), [=] AMREX_GPU_DEVICE (int i, int j, int k) noexcept {
                        if (!m(i,j,k)) { d(i,j,k) = Real(0); }
                    });
                }
            }
            const auto dx = geom[lev].CellSizeArray();
            mass += deposit.sum(0) * dx[0] * dx[1] * dx[2];
        }
        EXPECT_NEAR(mass, Real(1), Real(64) * std::numeric_limits<Real>::epsilon())
            << "particle z=" << zpos;
    }
}

// Motivation: CIC weight outside physical walls or on covered coarse cells must
// not lose or duplicate represented particle mass. A unit-mass particle checks
// composite deposition on both sides of a level edge and across a periodic seam.
TEST(SpeciesProjection, CompositeCICConservesParticleMassAtWallsAndLevelEdges)
{
    composite_cic_mass_test();
}

void composite_initialization_test ()
{
    const Box domain(IntVect(0), IntVect(3));
    const RealBox rb(Real(0), Real(0), Real(0), Real(4), Real(4), Real(4));
    Vector<Geometry> geom{Geometry(domain, rb, CoordSys::cartesian, {0,0,0}),
                          Geometry(refine(domain, 2), rb, CoordSys::cartesian, {0,0,0})};
    Box fine = geom[1].Domain();
    fine.setBig(2, 3);
    Vector<BoxArray> ba{BoxArray(domain), BoxArray(fine)};
    Vector<DistributionMapping> dm{DistributionMapping(ba[0]), DistributionMapping(ba[1])};
    ParGDB gdb(geom, dm, ba, {IntVect(2)});
    const std::string name = "composite_initialization_test";
    ParmParse pp(name);
    pp.add("num_initializations", 1);
    pp.add("num_injections", 0);
    pp.add("initial_number_density", Real(16));
    pp.add("initial_super_droplet_density", Real(16));
    pp.add("multiplicity_type", "constant");
    pp.add("initial_species_mean_mass_H2O", Real(1.e-12));
    pp.add("initial_species_mean_mass_water", Real(1.e-12));
    pp.add("initial_aerosol_mean_mass_NaCl", Real(1.e-19));
    pp.add("density_scaling", true);
    pp.add("place_randomly_in_cells", false);
    SuperDropletPC pc(&gdb, {Species::Name::H2O, Species::Name::water},
                     {Species::Name::NaCl}, 0.02, name);
    ParticleReal total = 0;
    for (int lev = 0; lev < 2; ++lev) {
        BoxArray nodes(ba[lev]);
        nodes.surroundingNodes();
        auto height = std::make_unique<MultiFab>(nodes, dm[lev], 1, 2);
        const Real dz = geom[lev].CellSize(2);
        for (MFIter mfi(*height); mfi.isValid(); ++mfi) {
            const auto h = height->array(mfi);
            ParallelFor(mfi.fabbox(), [=] AMREX_GPU_DEVICE (int i, int j, int k) noexcept {
                h(i,j,k) = Real(k) * dz;
            });
        }
        pc.InitializeParticles(lev, 0, height);
        MultiFab rho(ba[lev], dm[lev], 1, 1);
        rho.setVal(Real(lev + 2));
        pc.DensityScaling(rho, lev);
        EXPECT_EQ(pc.NumberOfParticlesAtLevel(lev), 512);
        for (SuperDropletPC::ParIterType pti(pc, lev); pti.isValid(); ++pti) {
            const auto& xi = pti.GetStructOfArrays().GetRealData(
                static_cast<int>(SDPCDefn::SuperDropletsRealIdx::ncomps) +
                SDPCDefn::SuperDropletsRealIdxSoA_RT::multiplicity);
            Vector<ParticleReal> multiplicity(xi.size());
            Gpu::copy(Gpu::deviceToHost, xi.begin(), xi.end(), multiplicity.begin());
            for (auto value : multiplicity) {
                EXPECT_EQ(value, ParticleReal(lev + 2));
                total += value;
            }
        }
        MultiFab condensate(ba[lev], dm[lev], 1, 0);
        condensate.setVal(Real(1.e-4));
        pc.SetAttributes(condensate, lev);
        for (SuperDropletPC::ParIterType pti(pc, lev); pti.isValid(); ++pti) {
            const auto& soa = pti.GetStructOfArrays();
            const auto& xi = soa.GetRealData(static_cast<int>(SDPCDefn::SuperDropletsRealIdx::ncomps) +
                                            SDPCDefn::SuperDropletsRealIdxSoA_RT::multiplicity);
            const auto& mass = soa.GetRealData(SDPCDefn::SuperDropletsRealIdx::mass);
            Vector<ParticleReal> multiplicity(xi.size()), per_drop(mass.size());
            Gpu::copy(Gpu::deviceToHost, xi.begin(), xi.end(), multiplicity.begin());
            Gpu::copy(Gpu::deviceToHost, mass.begin(), mass.end(), per_drop.begin());
            const Real expected = Real(1.e-4) / Real(16);
            for (int i = 0; i < static_cast<int>(multiplicity.size()); ++i) {
                EXPECT_NEAR(multiplicity[i] * (per_drop[i] - ParticleReal(1.e-12) - ParticleReal(1.e-19)),
                            expected, expected * Real(32) * std::numeric_limits<ParticleReal>::epsilon());
            }
        }
    }
    EXPECT_EQ(total, ParticleReal(2560));
}

// Motivation: Covered coarse cells must not seed a second initial population,
// and density scaling must use the flow on each owning level before particles
// move. Distinct level densities give an exact composite multiplicity oracle;
// volume-based sampling also checks the level-local condensate attribute divisor.
TEST(SpeciesProjection, CompositeInitializationScalesEachLevelOnce)
{
    composite_initialization_test();
}

void set_test_particle (SuperDropletPC& pc, ParticleReal multiplicity,
                        Real x, Real y, Real z, int tag)
{
    pc.InitRandom(1, 42, {});
    for (int lev = 0; lev <= pc.finestLevel(); ++lev) {
        for (SuperDropletPC::ParIterType pti(pc, lev); pti.isValid(); ++pti) {
            const auto data = pti.GetParticleTile().getParticleTileData();
            ParallelFor(pti.numParticles(), [=] AMREX_GPU_DEVICE (int i) noexcept {
                data.m_aos[i].pos(0) = ParticleReal(x);
                data.m_aos[i].pos(1) = ParticleReal(y);
                data.m_aos[i].pos(2) = ParticleReal(z);
                for (int c = 0; c < data.m_num_runtime_real; ++c) {
                    data.m_runtime_rdata[c][i] = ParticleReal(c+1) * ParticleReal(1.e-12);
                }
                data.m_runtime_rdata[SDPCDefn::SuperDropletsRealIdxSoA_RT::multiplicity][i] = multiplicity;
                data.m_runtime_rdata[SDPCDefn::SuperDropletsRealIdxSoA_RT::radius][i] = ParticleReal(1.e-6);
                data.m_runtime_idata[SDPCDefn::SuperDropletsIntIdxSoA_RT::active][i] = tag;
            });
        }
    }
    pc.Redistribute();
}

void limited_donor_test ()
{
    const Box domain(IntVect(0), IntVect(1));
    const RealBox rb(Real(0), Real(0), Real(0), Real(4), Real(4), Real(4));
    Vector<Geometry> geom{Geometry(domain, rb, CoordSys::cartesian, {0,0,0}),
                          Geometry(refine(domain, 2), rb, CoordSys::cartesian, {0,0,0})};
    Vector<BoxArray> ba{BoxArray(domain), BoxArray(geom[1].Domain())};
    Vector<DistributionMapping> dm{DistributionMapping(ba[0]), DistributionMapping(ba[1])};
    ParGDB gdb(geom, dm, ba, {IntVect(2)});
    const std::string name = "limited_donor_test";
    ParmParse pp(name);
    pp.add("num_initializations", 0);
    pp.add("num_injections", 0);
    pp.add("split_merge_amr", true);
    for (int mult : {1, 2, 3, 8}) {
        std::map<int, ParticleReal> reference;
        for (int repeat = 0; repeat < 4; ++repeat) {
            SuperDropletPC pc(&gdb, {Species::Name::H2O, Species::Name::water, Species::Name::agua},
                             {Species::Name::NaCl}, 0.02, name);
            set_test_particle(pc, ParticleReal(mult), Real(0.5), Real(0.5), Real(0.5), 2);
            pc.SplitMergeAtLevelBoundary();
            std::map<int, ParticleReal> cells;
            ParticleReal total = 0;
            ParticleReal species_mass = 0;
            int stored = 0;
            for (SuperDropletPC::ParIterType pti(pc, 1); pti.isValid(); ++pti) {
                const auto& tile = pti.GetParticleTile();
                const int np = tile.numParticles();
                stored += np;
                Vector<SuperDropletPC::ParticleType> particles(np);
                Gpu::copy(Gpu::deviceToHost, tile.GetArrayOfStructs().dataPtr(),
                          tile.GetArrayOfStructs().dataPtr()+np, particles.begin());
                Vector<ParticleReal> multiplicity(np), mass(np);
                const auto& soa = tile.GetStructOfArrays();
                const auto& xi = soa.GetRealData(static_cast<int>(SDPCDefn::SuperDropletsRealIdx::ncomps)
                                                + SDPCDefn::SuperDropletsRealIdxSoA_RT::multiplicity);
                Gpu::copy(Gpu::deviceToHost, xi.begin(), xi.end(), multiplicity.begin());
                for (int c = 0; c < soa.NumRealComps(); ++c) {
                    if (c < SDPCDefn::SuperDropletsRealIdx::ncomps) { continue; }
                    const int runtime = c - SDPCDefn::SuperDropletsRealIdx::ncomps;
                    if (runtime == SDPCDefn::SuperDropletsRealIdxSoA_RT::multiplicity ||
                        runtime == SDPCDefn::SuperDropletsRealIdxSoA_RT::radius) { continue; }
                    const auto& values = soa.GetRealData(c);
                    Gpu::copy(Gpu::deviceToHost, values.begin(), values.end(), mass.begin());
                    for (int i = 0; i < np; ++i) {
                        EXPECT_EQ(mass[i], ParticleReal(runtime+1) * ParticleReal(1.e-12));
                    }
                }
                const auto& sp = soa.GetRealData(idx_s(1,1,3));
                Gpu::copy(Gpu::deviceToHost, sp.begin(), sp.end(), mass.begin());
                for (int i = 0; i < np; ++i) {
                    EXPECT_GT(particles[i].id(), 0);
                    EXPECT_GE(multiplicity[i], ParticleReal(1));
                    const int b = int(particles[i].pos(0)) + 4 *
                                  (int(particles[i].pos(1)) + 4 * int(particles[i].pos(2)));
                    cells[b] = multiplicity[i];
                    total += multiplicity[i];
                    species_mass += multiplicity[i] * mass[i];
                }
            }
            EXPECT_EQ(total, ParticleReal(mult));
            const ParticleReal per_drop = ParticleReal(ridx_s(1,1,3)+1) * ParticleReal(1.e-12);
            EXPECT_NEAR(species_mass, ParticleReal(mult) * per_drop,
                        Real(16) * std::numeric_limits<ParticleReal>::epsilon() * mult * per_drop);
            EXPECT_EQ(stored, std::min(mult, 4));
            ASSERT_EQ(cells.size(), std::min(mult, 4));
            EXPECT_TRUE(cells.count(0));
            const int pieces = std::min(mult, 4);
            EXPECT_EQ(cells.at(0), ParticleReal(mult - (pieces-1) * (mult/pieces)));
            const int expected_bins[] = {0,1,4,16};
            for (int i = 1; i < pieces; ++i) {
                ASSERT_TRUE(cells.count(expected_bins[i]));
                EXPECT_EQ(cells.at(expected_bins[i]), ParticleReal(mult/pieces));
            }
            if (repeat == 0) { reference = cells; }
            else { EXPECT_EQ(cells, reference); }
        }
    }
}

// Motivation: Several empty fine cells can exhaust one donor's multiplicity.
// Repopulation must conserve its species mass, retain the donor, discard invalid
// daughters and choose the same surviving cells on repeated identical inputs,
// without changing per-particle species, aerosol or runtime metadata.
TEST(SpeciesProjection, LimitedDonorIsConservativeCompactAndDeterministic)
{
    limited_donor_test();
}

class ProjectionModel : public SuperDropletsMoist
{
public:
    void setup (SuperDropletPC& pc, const MultiFab& cons, SDMSimulationDim dim)
    {
        m_super_droplets = &pc;
        m_geom = pc.Geom(0);
        m_dimensionality = dim;
        m_kinematic_mode = false;
        m_istart_sp = 1;
        m_num_species = 3;
        m_num_aerosols = 1;
        m_qstate_nonmoist_size = 4;
        EnsureMicFabVars(0, cons);
        for (auto& f : m_mic_fab_vars[0]) { f->setVal(Real(-777)); }
        field(MicVar_SD::rho).setVal(Real(2));
    }
    MultiFab& field (int n) { return *m_mic_fab_vars[0][n]; }
};

void diagnostic_projection_test ()
{
    const Box domain(IntVect(0), IntVect(3));
    const RealBox rb(Real(0), Real(0), Real(0), Real(4), Real(4), Real(4));
    Geometry geom(domain, rb, CoordSys::cartesian, {1,1,0});
    for (auto dim : {SDMSimulationDim::one_d_z, SDMSimulationDim::three_d,
                     SDMSimulationDim::two_d_xz, SDMSimulationDim::two_d_yz}) {
        const bool one_d = dim == SDMSimulationDim::one_d_z;
        BoxArray ba(domain);
        if (!one_d) {
            Box left = domain, right = domain;
            left.setBig(0,1);
            right.setSmall(0,2);
            ba = BoxArray(Vector<Box>{left, right}.data(), 2);
        }
        DistributionMapping dm(ba);
        ParmParse pp("diagnostic_projection_test");
        pp.add("num_initializations",0);
        pp.add("num_injections",0);
        SuperDropletPC pc(geom, dm, ba,
                         {Species::Name::H2O, Species::Name::water, Species::Name::agua},
                         {Species::Name::NaCl}, 0.02, "diagnostic_projection_test");
        set_test_particle(pc, ParticleReal(8), one_d ? Real(0.5) : Real(1.9),
                          Real(0.5), Real(1.5), 1);
        MultiFab cons(ba, dm, RhoQ1_comp+10, 1);
        cons.setVal(Real(-999));
        cons.setVal(Real(2), Rho_comp, 1, 0);
        cons.setVal(Real(0.004), RhoQ7_comp, 1, 0);
        cons.setVal(Real(0.008), RhoQ9_comp, 1, 0);
        BoxArray nodes(ba);
        nodes.surroundingNodes();
        MultiFab height(nodes, dm, 1, 2);
        for (MFIter mfi(height); mfi.isValid(); ++mfi) {
            const auto h = height.array(mfi);
            ParallelFor(mfi.fabbox(), [=] AMREX_GPU_DEVICE (int i, int j, int k) noexcept {
                h(i,j,k) = Real(k);
            });
        }
        ProjectionModel model;
        model.setup(pc, cons, dim);
        model.computeQcSpecies(height);
        MultiFab reference(ba, dm, 1, 1);
        MultiFab::Copy(reference, model.field(s_qc_idx(1,1)), 0, 0, 1, 1);
        model.field(s_qc_idx(1,1)).setVal(Real(-888));
        model.field(MicVar_SD::rho).setVal(Real(-888));
        model.Update_Diagnostic_State(cons, height);
        MultiFab::Subtract(reference, model.field(s_qc_idx(1,1)), 0, 0, 1, 1);
        EXPECT_EQ(reference.norm0(0, 1, IntVect(1)), Real(0));
        const Real represented = Real(8) * Real(ridx_s(1,1,3)+1) * Real(1.e-12);
        EXPECT_NEAR(cons.sum(RhoQ8_comp), represented * (one_d ? Real(16) : Real(1)),
                    represented * Real(256) * std::numeric_limits<ParticleReal>::epsilon());
        for (MFIter mfi(cons); mfi.isValid(); ++mfi) {
            const Box box = mfi.fabbox();
            FArrayBox c(box, 1, The_Pinned_Arena()), v(box, 1, The_Pinned_Arena());
            FArrayBox t(box, 1, The_Pinned_Arena()), r(box, 1, The_Pinned_Arena());
            const auto& cf = model.field(s_qc_idx(1,1))[mfi];
            const auto& vf = model.field(s_qv_idx(1,1))[mfi];
            const auto& tf = model.field(s_qt_idx(1,1))[mfi];
            const auto& rf = model.field(MicVar_SD::rho)[mfi];
            Gpu::dtoh_memcpy(c.dataPtr(), cf.dataPtr(), cf.nBytes());
            Gpu::dtoh_memcpy(v.dataPtr(), vf.dataPtr(), vf.nBytes());
            Gpu::dtoh_memcpy(t.dataPtr(), tf.dataPtr(), tf.nBytes());
            Gpu::dtoh_memcpy(r.dataPtr(), rf.dataPtr(), rf.nBytes());
            for (BoxIterator it(box); it.ok(); ++it) {
                const auto cell = it();
                EXPECT_EQ(r.const_array()(cell), Real(2));
                EXPECT_EQ(v.const_array()(cell), Real(0.002));
                EXPECT_GE(c.const_array()(cell), Real(0));
                EXPECT_TRUE(std::isfinite(t.const_array()(cell)));
                EXPECT_EQ(t.const_array()(cell), v.const_array()(cell)+c.const_array()(cell));
                if (one_d && mfi.validbox().contains(cell)) {
                    EXPECT_EQ(c.const_array()(cell), c.const_array()(0,0,cell[2]));
                }
            }
        }
        EXPECT_EQ(cons.min(RhoQ1_comp, 0), Real(-999)); // water is untouched
        EXPECT_EQ(cons.min(RhoQ7_comp, 0), Real(0.004));
    }
}

// Motivation: Redistribution-only refresh must preserve normal SDM projection
// semantics, including one_d_z horizontal replication, without mixing density
// and mixing-ratio representations in ghosts. Sentinel-filled, split/periodic
// boxes expose stale data while checking that water and conserved vapor stay
// unchanged in every supported dimensional mode.
TEST(SpeciesProjection, DiagnosticRefreshMatchesNormalDimensionsAndGhostRepresentation)
{
    diagnostic_projection_test();
}
