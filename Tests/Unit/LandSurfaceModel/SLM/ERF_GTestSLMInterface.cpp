#include <array>
#include <cmath>
#include <filesystem>
#include <fstream>
#include <memory>
#include <string>
#include <utility>
#include <vector>

#include <AMReX_Box.H>
#include <AMReX_BoxArray.H>
#include <AMReX_DistributionMapping.H>
#include <AMReX_Geometry.H>
#include <AMReX_MultiFab.H>
#include <AMReX_ParmParse.H>
#include <AMReX_RealBox.H>
#include <AMReX_Utility.H>

#include <gtest/gtest.h>

#include "ERF_LandSurface.H"
#include "ERF_SLM.H"

namespace {

amrex::MultiFab* get_slm_data (LandSurface& land_surface, const int lev,
                               const char* name)
{
    return land_surface.Get_Data_Ptr(lev, land_surface.Get_DataIdx(lev, name));
}

/**
 * Sets a horizontal surface-temperature gradient for a test state.
 */
void set_surface_temperature_gradient (LandSurface& land_surface,
                                       const int lev,
                                       const amrex::Real base_temperature)
{
    auto* tsurf = get_slm_data(land_surface, lev, "tsurf");
    for (amrex::MFIter mfi(*tsurf, amrex::TilingIfNotGPU()); mfi.isValid(); ++mfi) {
        const amrex::Box box = mfi.validbox();
        auto tsurf_arr = tsurf->array(mfi);
        amrex::ParallelFor(box, [=] AMREX_GPU_DEVICE (int i, int j, int k) noexcept {
            tsurf_arr(i, j, k) = base_temperature + i + 2.0 * j;
        });
    }
}

struct SLMTestState {
    amrex::Box domain;
    amrex::BoxArray ba;
    amrex::DistributionMapping dm;
    amrex::Geometry geom;
    amrex::MultiFab cons;
    amrex::MultiFab uvel;
    amrex::MultiFab vvel;

    explicit SLMTestState (const amrex::Box& test_domain = amrex::Box(
        amrex::IntVect(AMREX_D_DECL(0, 0, 0)),
        amrex::IntVect(AMREX_D_DECL(1, 1, 7))))
        : domain(test_domain),
          ba(domain),
          dm(ba),
          geom(make_geometry(domain)),
          cons(ba, dm, 8, amrex::IntVect(0)),
          uvel(amrex::BoxArray(amrex::convert(domain, amrex::IntVect(1, 0, 0))),
               dm,
               1, amrex::IntVect(0)),
          vvel(amrex::BoxArray(amrex::convert(domain, amrex::IntVect(0, 1, 0))),
               dm,
               1, amrex::IntVect(0))
    {
        cons.setVal(0.0);
        uvel.setVal(0.0);
        vvel.setVal(0.0);
    }

    static amrex::Geometry make_geometry (const amrex::Box& box)
    {
        amrex::RealBox real_box({AMREX_D_DECL(0.0, 0.0, 0.0)},
                                 {AMREX_D_DECL(2.0, 2.0, 8.0)});
        std::array<int, AMREX_SPACEDIM> is_periodic{};
        return amrex::Geometry(box, &real_box, 0, is_periodic.data());
    }
};

void configure_slm_parameters ()
{
    amrex::ParmParse pp("slm");
    pp.add("nsoil", 7);
    pp.addarr("soil_dz", std::vector<amrex::Real>{
        amrex::Real(0.1), amrex::Real(0.2), amrex::Real(0.3),
        amrex::Real(0.4), amrex::Real(0.5), amrex::Real(0.6),
        amrex::Real(0.7)});
    pp.add("landtype0", 10);
    pp.add("LAI0", amrex::Real(2.0));
    pp.add("clay0", amrex::Real(13.0));
    pp.add("sand0", amrex::Real(17.0));
    pp.addarr("sw0", std::vector<amrex::Real>{
        amrex::Real(0.60), amrex::Real(0.61), amrex::Real(0.62),
        amrex::Real(0.63), amrex::Real(0.64), amrex::Real(0.65),
        amrex::Real(0.66)});
    pp.addarr("st0", std::vector<amrex::Real>{
        amrex::Real(300.15), amrex::Real(300.14), amrex::Real(300.13),
        amrex::Real(300.12), amrex::Real(300.11), amrex::Real(300.10),
        amrex::Real(300.09)});
    pp.addarr("relax_hgt", std::vector<amrex::Real>{
        amrex::Real(0.0), amrex::Real(0.0), amrex::Real(0.0),
        amrex::Real(0.0), amrex::Real(0.0), amrex::Real(0.0),
        amrex::Real(1.0)});
}

amrex::BoxArray make_surface_box_array (const amrex::BoxArray& soil_boxes)
{
    amrex::BoxList surface_boxes = soil_boxes.boxList();
    for (auto& box : surface_boxes) {
        box.setRange(2, 0);
    }
    return amrex::BoxArray(std::move(surface_boxes));
}

struct SyntheticRadiationResult {
    amrex::Real soil_temperature;
    amrex::Real surface_temperature;
    std::array<amrex::Real, 6> radiation_fields;
};

SyntheticRadiationResult run_synthetic_radiation_case (
    SolverChoice& solver_choice,
    const std::array<amrex::Real, 6>& radiation)
{
    SLMTestState state;
    LandSurface land_surface;
    land_surface.ReSize(1);
    land_surface.SetModel<SLM>();
    land_surface.Define(0, solver_choice);

    std::unique_ptr<amrex::MultiFab> z_phys_nd;
    amrex::Vector<amrex::BCRec> domain_bcs_type;
    amrex::IntVect ref_ratio(1);
    amrex::Vector<amrex::Vector<std::string>> nc_init_file;
    land_surface.Init(0, state.cons, state.uvel, state.vvel,
                      state.geom, state.geom, domain_bcs_type, ref_ratio,
                      amrex::Real(1.0), z_phys_nd, nc_init_file);
    state.cons.setVal(1.0);
    state.cons.setVal(amrex::Real(1.0), Rho_comp, 1);
    state.cons.setVal(amrex::Real(300.0), RhoTheta_comp, 1);
    state.cons.setVal(amrex::Real(0.01), RhoQ1_comp, 1);
    state.uvel.setVal(0.1);
    state.vvel.setVal(0.1);

    const auto lsm_box = get_slm_data(land_surface, 0, "tsurf")->boxArray();
    const auto lsm_dm = get_slm_data(land_surface, 0, "tsurf")->DistributionMap();
    const auto surface_box = make_surface_box_array(lsm_box);

    amrex::Vector<amrex::Vector<std::unique_ptr<amrex::MultiFab>>> sst(1);
    amrex::Vector<amrex::Vector<std::unique_ptr<amrex::iMultiFab>>> landmask(1);
    sst[0].resize(1);
    landmask[0].resize(1);
    sst[0][0] = std::make_unique<amrex::MultiFab>(surface_box, lsm_dm, 1, amrex::IntVect(0));
    landmask[0][0] = std::make_unique<amrex::iMultiFab>(surface_box, lsm_dm, 1, amrex::IntVect(0));
    sst[0][0]->setVal(301.0);
    landmask[0][0]->setVal(1);

    amrex::MultiFab precip(surface_box, lsm_dm, 1, amrex::IntVect(0));
    precip.setVal(0.0);

    get_slm_data(land_surface, 0, "SW_dw_dir_vis")->setVal(radiation[0]);
    get_slm_data(land_surface, 0, "SW_dw_dir_nir")->setVal(radiation[1]);
    get_slm_data(land_surface, 0, "SW_dw_dif_vis")->setVal(radiation[2]);
    get_slm_data(land_surface, 0, "SW_dw_dif_nir")->setVal(radiation[3]);
    get_slm_data(land_surface, 0, "LW_dw")->setVal(radiation[4]);
    get_slm_data(land_surface, 0, "cos_zenith")->setVal(radiation[5]);

    land_surface.Update_Lsm_Vars_Lev(0, state.cons, state.uvel, state.vvel);
    land_surface.set_LSM_terrain_inputs(0, sst, landmask);
    land_surface.set_LSM_precip_input(0, &precip);
    land_surface.Advance(0, amrex::Real(0.3), amrex::Real(0.3), amrex::Real(0.0));

    return {
        get_slm_data(land_surface, 0, "tsoil")->max(0),
        get_slm_data(land_surface, 0, "tsurf")->max(0),
        {
            get_slm_data(land_surface, 0, "SW_dw_dir_vis")->max(0),
            get_slm_data(land_surface, 0, "SW_dw_dir_nir")->max(0),
            get_slm_data(land_surface, 0, "SW_dw_dif_vis")->max(0),
            get_slm_data(land_surface, 0, "SW_dw_dif_nir")->max(0),
            get_slm_data(land_surface, 0, "LW_dw")->max(0),
            get_slm_data(land_surface, 0, "cos_zenith")->max(0)}};
}

void set_synthetic_reference_atmosphere (
    LandSurface& land_surface,
    const amrex::Real pressure_mb,
    const amrex::Real temperature_k,
    const amrex::Real humidity_g_per_kg,
    const amrex::Real u_velocity,
    const amrex::Real v_velocity)
{
    // Values and units follow the first row of sndref:
    // pressure [mb], temperature [K], q [g/kg], and wind [m/s].
    get_slm_data(land_surface, 0, "ref_t")->setVal(temperature_k);
    get_slm_data(land_surface, 0, "ref_q")->setVal(humidity_g_per_kg / 1000.0);
    get_slm_data(land_surface, 0, "ref_p")->setVal(pressure_mb);
    get_slm_data(land_surface, 0, "ref_d")->setVal(
        pressure_mb * 100.0 / (amrex::Real(287.0) * temperature_k));
    get_slm_data(land_surface, 0, "ref_u")->setVal(u_velocity);
    get_slm_data(land_surface, 0, "ref_v")->setVal(v_velocity);
}

class SLMInterfaceTest : public ::testing::Test {
protected:
    static amrex::Vector<amrex::BCRec> make_domain_bcs_type ()
    {
        const amrex::BCRec interior_bc(
            AMREX_D_DECL(amrex::BCType::int_dir, amrex::BCType::int_dir,
                         amrex::BCType::int_dir),
            AMREX_D_DECL(amrex::BCType::int_dir, amrex::BCType::int_dir,
                         amrex::BCType::int_dir));
        return amrex::Vector<amrex::BCRec>(AMREX_SPACEDIM + NBCVAR_max,
                                           interior_bc);
    }

    static void SetUpTestSuite ()
    {
        configure_slm_parameters();
    }

    void SetUp () override
    {
        state = std::make_unique<SLMTestState>();
        solver_choice = SolverChoice{};
        solver_choice.terrain_type = TerrainType::None;
        solver_choice.init_type = InitType::Input_Sounding;

        land_surface.ReSize(1);
        land_surface.SetModel<SLM>();
        land_surface.Define(0, solver_choice);

        std::unique_ptr<amrex::MultiFab> z_phys_nd;
        amrex::Vector<amrex::BCRec> domain_bcs_type = make_domain_bcs_type();
        amrex::IntVect ref_ratio(1);
        amrex::Vector<amrex::Vector<std::string>> nc_init_file;
        land_surface.Init(0, state->cons, state->uvel, state->vvel,
                          state->geom, state->geom, domain_bcs_type, ref_ratio,
                          amrex::Real(1.0), z_phys_nd, nc_init_file);

        // Use a physically valid synthetic atmospheric state.  Setting every
        // conservative component to one produces an approximately 1 K
        // atmosphere and violates the saturation-pressure routine's range.
        state->cons.setVal(1.0);
        state->cons.setVal(amrex::Real(1.0), Rho_comp, 1);
        state->cons.setVal(amrex::Real(300.0), RhoTheta_comp, 1);
        state->cons.setVal(amrex::Real(0.01), RhoQ1_comp, 1);
        state->uvel.setVal(0.1);
        state->vvel.setVal(0.1);
    }

    void initialize_multilevel (const int nlevels,
                                const InitType init_type = InitType::Input_Sounding)
    {
        solver_choice.init_type = init_type;
        multilevel_states.clear();
        multilevel_states.reserve(nlevels);

        for (int lev = 0; lev < nlevels; ++lev) {
            const int xy_hi = (1 << (lev + 1)) - 1;
            const amrex::Box level_box(
                amrex::IntVect(AMREX_D_DECL(0, 0, 0)),
                amrex::IntVect(AMREX_D_DECL(xy_hi, xy_hi, 7)));
            multilevel_states.push_back(std::make_unique<SLMTestState>(level_box));
        }

        land_surface.ReSize(nlevels);
        land_surface.SetModel<SLM>();
        for (int lev = 0; lev < nlevels; ++lev) {
            land_surface.Define(lev, solver_choice);
        }

        for (int lev = 0; lev < nlevels; ++lev) {
            std::unique_ptr<amrex::MultiFab> z_phys_nd;
            amrex::Vector<amrex::BCRec> domain_bcs_type = make_domain_bcs_type();
            amrex::IntVect ref_ratio = (lev == 0) ? amrex::IntVect(1) : amrex::IntVect(2);
            amrex::Vector<amrex::Vector<std::string>> nc_init_file;
            land_surface.Init(lev, multilevel_states[lev]->cons,
                              multilevel_states[lev]->uvel,
                              multilevel_states[lev]->vvel,
                              multilevel_states[lev]->geom,
                              multilevel_states[0]->geom,
                              domain_bcs_type, ref_ratio, amrex::Real(1.0),
                              z_phys_nd, nc_init_file);
        }
    }

    void seed_transfer_state (const int lev,
                              const amrex::Real tsurf,
                              const amrex::Real soil_temperature,
                              const amrex::Real soil_moisture,
                              const amrex::Real flux)
    {
        get_slm_data(land_surface, lev, "tsurf")->setVal(tsurf);
        get_slm_data(land_surface, lev, "tsoil")->setVal(soil_temperature);
        get_slm_data(land_surface, lev, "wsoil")->setVal(soil_moisture);
        get_slm_data(land_surface, lev, "vegtype")->setVal(10.0);
        get_slm_data(land_surface, lev, "soiltype")->setVal(1.0);
        get_slm_data(land_surface, lev, "lai")->setVal(2.0);
        get_slm_data(land_surface, lev, "veg_frac")->setVal(1.0);
        get_slm_data(land_surface, lev, "veg_frac_min")->setVal(0.0);
        get_slm_data(land_surface, lev, "veg_frac_max")->setVal(0.0);
        for (int var = 0; var < land_surface.Get_Flux_Size(); ++var) {
            land_surface.Get_Flux_Ptr(lev, var)->setVal(flux);
        }
    }

    amrex::Real surface_temperature_at (const int lev,
                                        const amrex::IntVect& iv)
    {
        auto* tsurf = get_slm_data(land_surface, lev, "tsurf");
        for (amrex::MFIter mfi(*tsurf, amrex::TilingIfNotGPU()); mfi.isValid(); ++mfi) {
            if (mfi.validbox().contains(iv)) {
                return tsurf->const_array(mfi)(iv[0], iv[1], iv[2]);
            }
        }
        return amrex::Real(-1.0);
    }

    std::unique_ptr<SLMTestState> state;
    std::vector<std::unique_ptr<SLMTestState>> multilevel_states;
    SolverChoice solver_choice;
    LandSurface land_surface;
};

TEST_F(SLMInterfaceTest, ExposesCurrentStateAndFluxContract)
{
    const std::vector<std::string> expected_data_names{
        "tsurf", "ustar", "tstar", "qstar", "tveg", "mv", "tsoil",
        "wsoil", "sand", "clay", "soil_thickness", "surface_u", "surface_v",
        "surface_vapor", "surface_heat", "precip_soil", "ref_precip",
        "SW_dw_dir_vis", "SW_dw_dir_nir", "SW_dw_dif_vis", "SW_dw_dif_nir",
        "LW_dw", "cos_zenith", "ref_t", "ref_u", "ref_v", "ref_d", "ref_q",
        "ref_p", "node_z", "soilt_nudge", "soilw_nudge", "lai", "vegtype",
        "soiltype", "veg_frac", "veg_frac_min", "veg_frac_max", "emis_sfc",
        "alb_nir_sfc", "alb_vis_sfc", "alb_nir_sfc_diff", "alb_vis_sfc_diff",
        "soil_transp_frac"};
    const std::vector<std::string> expected_flux_names{
        "t_flux", "q_flux", "tau13", "tau23", "olen"};

    ASSERT_EQ(land_surface.Get_Data_Size(), static_cast<int>(expected_data_names.size()));
    ASSERT_EQ(land_surface.Get_Flux_Size(), static_cast<int>(expected_flux_names.size()));

    for (int var = 0; var < land_surface.Get_Data_Size(); ++var) {
        EXPECT_EQ(land_surface.Get_DataName(var), expected_data_names[var]);
        std::string name = expected_data_names[var];
        EXPECT_EQ(land_surface.Get_DataIdx(0, name), var);
        EXPECT_NE(land_surface.Get_Data_Ptr(0, var), nullptr);
    }

    for (int var = 0; var < land_surface.Get_Flux_Size(); ++var) {
        EXPECT_EQ(land_surface.Get_FluxName(var), expected_flux_names[var]);
        std::string name = expected_flux_names[var];
        EXPECT_EQ(land_surface.Get_FluxIdx(0, name), var);
        EXPECT_NE(land_surface.Get_Flux_Ptr(0, var), nullptr);
    }
}

TEST_F(SLMInterfaceTest, PreservesSoilGeometryAndFieldLocations)
{
    const amrex::Box expected_soil_box(
        amrex::IntVect(AMREX_D_DECL(0, 0, -7)),
        amrex::IntVect(AMREX_D_DECL(1, 1, -1)));
    const amrex::Box expected_flux_box(
        amrex::IntVect(AMREX_D_DECL(0, 0, 0)),
        amrex::IntVect(AMREX_D_DECL(1, 1, 0)));

    const auto* soil_temperature = get_slm_data(land_surface, 0, "tsoil");
    const auto* surface_temperature = get_slm_data(land_surface, 0, "tsurf");
    const auto* soil_flux = land_surface.Get_Flux_Ptr(0, LsmFlux_SLM::t_flux);
    const auto* obukhov_length = land_surface.Get_Flux_Ptr(0, LsmFlux_SLM::olen);

    ASSERT_NE(soil_temperature, nullptr);
    ASSERT_NE(surface_temperature, nullptr);
    ASSERT_NE(soil_flux, nullptr);
    ASSERT_NE(obukhov_length, nullptr);

    EXPECT_EQ(soil_temperature->boxArray()[0], expected_soil_box);
    EXPECT_EQ(surface_temperature->boxArray()[0], expected_soil_box);
    EXPECT_EQ(soil_flux->boxArray()[0], expected_flux_box);
    EXPECT_EQ(obukhov_length->boxArray()[0], expected_flux_box);
    EXPECT_EQ(land_surface.Get_Lsm_Geom(0).Domain(), expected_soil_box);
}

TEST_F(SLMInterfaceTest, InitializesConfiguredSoilAndSurfaceFields)
{
    const auto* soil_temperature = get_slm_data(land_surface, 0, "tsoil");
    const auto* soil_moisture = get_slm_data(land_surface, 0, "wsoil");
    const auto* soil_depth = get_slm_data(land_surface, 0, "soil_thickness");
    const auto* node_z = get_slm_data(land_surface, 0, "node_z");
    const auto* sand = get_slm_data(land_surface, 0, "sand");
    const auto* clay = get_slm_data(land_surface, 0, "clay");
    const auto* tsurf = get_slm_data(land_surface, 0, "tsurf");
    const auto* lai = get_slm_data(land_surface, 0, "lai");
    const auto* vegtype = get_slm_data(land_surface, 0, "vegtype");
    const auto* veg_frac = get_slm_data(land_surface, 0, "veg_frac");
    const auto* veg_frac_min = get_slm_data(land_surface, 0, "veg_frac_min");
    const auto* veg_frac_max = get_slm_data(land_surface, 0, "veg_frac_max");
    const auto* emissivity = get_slm_data(land_surface, 0, "emis_sfc");
    const auto* albedo_nir = get_slm_data(land_surface, 0, "alb_nir_sfc");
    const auto* albedo_vis = get_slm_data(land_surface, 0, "alb_vis_sfc");
    const auto* albedo_nir_diff = get_slm_data(land_surface, 0, "alb_nir_sfc_diff");
    const auto* albedo_vis_diff = get_slm_data(land_surface, 0, "alb_vis_sfc_diff");
    const auto* olen = land_surface.Get_Flux_Ptr(0, LsmFlux_SLM::olen);

    ASSERT_NE(soil_temperature, nullptr);
    ASSERT_NE(soil_moisture, nullptr);
    ASSERT_NE(soil_depth, nullptr);
    ASSERT_NE(node_z, nullptr);
    ASSERT_NE(sand, nullptr);
    ASSERT_NE(clay, nullptr);
    ASSERT_NE(tsurf, nullptr);
    ASSERT_NE(lai, nullptr);
    ASSERT_NE(vegtype, nullptr);
    ASSERT_NE(veg_frac, nullptr);
    ASSERT_NE(veg_frac_min, nullptr);
    ASSERT_NE(veg_frac_max, nullptr);
    ASSERT_NE(emissivity, nullptr);
    ASSERT_NE(albedo_nir, nullptr);
    ASSERT_NE(albedo_vis, nullptr);
    ASSERT_NE(albedo_nir_diff, nullptr);
    ASSERT_NE(albedo_vis_diff, nullptr);
    ASSERT_NE(olen, nullptr);

    {
        amrex::MFIter mfi(*soil_temperature, amrex::TilingIfNotGPU());
        ASSERT_TRUE(mfi.isValid());
        const auto& fab = (*soil_temperature)[mfi];
        const auto soil_temperature_arr = soil_temperature->const_array(mfi);
        const auto soil_moisture_arr = soil_moisture->const_array(mfi);
        const auto soil_depth_arr = soil_depth->const_array(mfi);
        const auto node_z_arr = node_z->const_array(mfi);
        const auto sand_arr = sand->const_array(mfi);
        const auto clay_arr = clay->const_array(mfi);
        const auto tsurf_arr = tsurf->const_array(mfi);
        const auto lai_arr = lai->const_array(mfi);
        const auto vegtype_arr = vegtype->const_array(mfi);
        const auto veg_frac_arr = veg_frac->const_array(mfi);
        const auto veg_frac_min_arr = veg_frac_min->const_array(mfi);
        const auto veg_frac_max_arr = veg_frac_max->const_array(mfi);
        const auto emissivity_arr = emissivity->const_array(mfi);
        const auto albedo_nir_arr = albedo_nir->const_array(mfi);
        const auto albedo_vis_arr = albedo_vis->const_array(mfi);
        const auto albedo_nir_diff_arr = albedo_nir_diff->const_array(mfi);
        const auto albedo_vis_diff_arr = albedo_vis_diff->const_array(mfi);
        const auto olen_arr = olen->const_array(mfi);

    ASSERT_TRUE(fab.box().contains(amrex::IntVect(AMREX_D_DECL(0, 0, -1))));
    for (int layer = 0; layer < 7; ++layer) {
        const int k = -1 - layer;
        EXPECT_NEAR(soil_temperature_arr(0, 0, k),
                    amrex::Real(300.15 - 0.01 * layer), amrex::Real(1.0e-12));
        EXPECT_NEAR(soil_moisture_arr(0, 0, k),
                    amrex::Real(0.60 + 0.01 * layer), amrex::Real(1.0e-12));
        EXPECT_NEAR(soil_depth_arr(0, 0, k),
                    amrex::Real(0.1 + 0.1 * layer), amrex::Real(1.0e-12));
        amrex::Real expected_node_z = amrex::Real(0.5) * (amrex::Real(0.1) + amrex::Real(0.1 * layer));
        for (int previous = 0; previous < layer; ++previous) {
            expected_node_z += amrex::Real(0.1 + 0.1 * previous);
        }
        EXPECT_NEAR(node_z_arr(0, 0, k), expected_node_z, amrex::Real(1.0e-12));
        EXPECT_NEAR(sand_arr(0, 0, k), amrex::Real(17.0), amrex::Real(1.0e-12));
        EXPECT_NEAR(clay_arr(0, 0, k), amrex::Real(13.0), amrex::Real(1.0e-12));
    }
    EXPECT_NEAR(tsurf_arr(0, 0, -1), amrex::Real(300.15), amrex::Real(1.0e-12));
    EXPECT_NEAR(lai_arr(0, 0, -1), amrex::Real(2.0), amrex::Real(1.0e-12));
    EXPECT_EQ(vegtype_arr(0, 0, -1), 10);
    EXPECT_NEAR(veg_frac_arr(0, 0, -1), amrex::Real(1.0), amrex::Real(1.0e-12));
    EXPECT_NEAR(veg_frac_min_arr(0, 0, -1), amrex::Real(0.0), amrex::Real(1.0e-12));
    EXPECT_NEAR(veg_frac_max_arr(0, 0, -1), amrex::Real(0.0), amrex::Real(1.0e-12));
    EXPECT_TRUE(std::isfinite(emissivity_arr(0, 0, -1)));
    EXPECT_TRUE(std::isfinite(albedo_nir_arr(0, 0, -1)));
    EXPECT_TRUE(std::isfinite(albedo_vis_arr(0, 0, -1)));
    EXPECT_TRUE(std::isfinite(albedo_nir_diff_arr(0, 0, -1)));
    EXPECT_TRUE(std::isfinite(albedo_vis_diff_arr(0, 0, -1)));
        EXPECT_EQ(olen_arr(0, 0, 0), amrex::Real(1.0e34));
    }

    for (int var = 0; var < land_surface.Get_Data_Size(); ++var) {
        EXPECT_FALSE(land_surface.Get_Data_Ptr(0, var)->contains_nan());
    }
    for (int var = 0; var < land_surface.Get_Flux_Size(); ++var) {
        EXPECT_FALSE(land_surface.Get_Flux_Ptr(0, var)->contains_nan());
    }
}

TEST_F(SLMInterfaceTest, RadiationExportUsesCurrentSLMFields)
{
    SLM* slm = land_surface.get_model_lev<SLM>(0);
    ASSERT_NE(slm, nullptr);

    const auto& radiation_fields = slm->export_to_RRTMGP();
    ASSERT_EQ(static_cast<int>(radiation_fields.size()), 6);
    EXPECT_EQ(radiation_fields[0], get_slm_data(land_surface, 0, "tsurf"));
    EXPECT_EQ(radiation_fields[1], get_slm_data(land_surface, 0, "emis_sfc"));
    EXPECT_EQ(radiation_fields[2], get_slm_data(land_surface, 0, "alb_vis_sfc"));
    EXPECT_EQ(radiation_fields[3], get_slm_data(land_surface, 0, "alb_nir_sfc"));
    EXPECT_EQ(radiation_fields[4], get_slm_data(land_surface, 0, "alb_vis_sfc_diff"));
    EXPECT_EQ(radiation_fields[5], get_slm_data(land_surface, 0, "alb_nir_sfc_diff"));

    const auto& output_map = slm->get_rad_output_map();
    EXPECT_EQ(output_map.at("cos_zenith_angle"), "cos_zenith");
    EXPECT_EQ(output_map.at("sw_flux_dn_dir_vis"), "SW_dw_dir_vis");
    EXPECT_EQ(output_map.at("sw_flux_dn_dir_nir"), "SW_dw_dir_nir");
    EXPECT_EQ(output_map.at("sw_flux_dn_dif_vis"), "SW_dw_dif_vis");
    EXPECT_EQ(output_map.at("sw_flux_dn_dif_nir"), "SW_dw_dif_nir");
    EXPECT_EQ(output_map.at("lw_flux_dn"), "LW_dw");
}

TEST_F(SLMInterfaceTest, ExposesWRFInputFieldMapping)
{
    const auto& mapping = land_surface.Get_WRFInputNames();
    const std::vector<std::pair<std::string, std::string>> expected{
        {"DZS", "soil_thickness"}, {"ZS", "node_z"},
        {"TSLB", "tsoil"}, {"SMOIS", "wsoil"},
        {"LAI", "lai"}, {"IVGTYP", "vegtype"},
        {"ISLTYP", "soiltype"}, {"TSK", "tsurf"},
        {"VEGFRA", "veg_frac"}, {"SHDMIN", "veg_frac_min"},
        {"SHDMAX", "veg_frac_max"}};

    ASSERT_EQ(mapping.size(), expected.size());
    for (const auto& [wrf_name, slm_name] : expected) {
        ASSERT_NE(mapping.find(wrf_name), mapping.end());
        EXPECT_EQ(mapping.at(wrf_name), slm_name);
    }
}

TEST_F(SLMInterfaceTest, ConsumesSyntheticAtmosphericRadiationFields)
{
    // Representative values based on the CASS radiation fixture:
    // direct SW visible/NIR, diffuse SW visible/NIR, LW, and cos(zenith).
    const std::array<amrex::Real, 6> radiation{
        amrex::Real(250.0), amrex::Real(230.0), amrex::Real(40.0),
        amrex::Real(10.0), amrex::Real(400.0), amrex::Real(0.6)};

    const auto result = run_synthetic_radiation_case(solver_choice, radiation);

    EXPECT_TRUE(std::isfinite(result.soil_temperature));
    EXPECT_TRUE(std::isfinite(result.surface_temperature));
    EXPECT_NEAR(result.radiation_fields[0],
                radiation[0], amrex::Real(1.0e-12));
    EXPECT_NEAR(result.radiation_fields[1],
                radiation[1], amrex::Real(1.0e-12));
    EXPECT_NEAR(result.radiation_fields[2],
                radiation[2], amrex::Real(1.0e-12));
    EXPECT_NEAR(result.radiation_fields[3],
                radiation[3], amrex::Real(1.0e-12));
    EXPECT_NEAR(result.radiation_fields[4],
                radiation[4], amrex::Real(1.0e-12));
    EXPECT_NEAR(result.radiation_fields[5],
                radiation[5], amrex::Real(1.0e-12));
}

TEST_F(SLMInterfaceTest, LongwaveRadiationChangesSLMThermalResponse)
{
    const std::array<amrex::Real, 6> no_longwave{
        amrex::Real(0.0), amrex::Real(0.0), amrex::Real(0.0),
        amrex::Real(0.0), amrex::Real(0.0), amrex::Real(0.0)};
    const std::array<amrex::Real, 6> cass_longwave{
        amrex::Real(0.0), amrex::Real(0.0), amrex::Real(0.0),
        amrex::Real(0.0), amrex::Real(400.0), amrex::Real(0.0)};

    const auto without_longwave =
        run_synthetic_radiation_case(solver_choice, no_longwave);
    const auto with_longwave =
        run_synthetic_radiation_case(solver_choice, cass_longwave);

    EXPECT_TRUE(std::isfinite(without_longwave.soil_temperature));
    EXPECT_TRUE(std::isfinite(with_longwave.soil_temperature));
    EXPECT_GT(std::abs(with_longwave.soil_temperature -
                       without_longwave.soil_temperature),
              amrex::Real(1.0e-10));
}

TEST_F(SLMInterfaceTest, PositiveZenithAngleEnablesShortwaveResponse)
{
    const std::array<amrex::Real, 6> zero_cosine{
        amrex::Real(250.0), amrex::Real(230.0), amrex::Real(40.0),
        amrex::Real(10.0), amrex::Real(400.0), amrex::Real(0.0)};
    const std::array<amrex::Real, 6> positive_cosine{
        amrex::Real(250.0), amrex::Real(230.0), amrex::Real(40.0),
        amrex::Real(10.0), amrex::Real(400.0), amrex::Real(0.6)};

    const auto without_shortwave =
        run_synthetic_radiation_case(solver_choice, zero_cosine);
    const auto with_shortwave =
        run_synthetic_radiation_case(solver_choice, positive_cosine);

    EXPECT_TRUE(std::isfinite(without_shortwave.soil_temperature));
    EXPECT_TRUE(std::isfinite(with_shortwave.soil_temperature));
    EXPECT_GT(std::abs(with_shortwave.soil_temperature -
                       without_shortwave.soil_temperature),
              amrex::Real(1.0e-10));
}

TEST_F(SLMInterfaceTest, AcceptsSyntheticReferenceAtmosphereRow)
{
    // First row of swlw: zero shortwave, 367 W/m2 downwelling longwave,
    // with the reference-input convention cos(zenith) = 1.
    const std::array<amrex::Real, 6> reference_radiation{
        amrex::Real(0.0), amrex::Real(0.0), amrex::Real(0.0),
        amrex::Real(0.0), amrex::Real(367.0), amrex::Real(1.0)};
    get_slm_data(land_surface, 0, "SW_dw_dir_vis")->setVal(reference_radiation[0]);
    get_slm_data(land_surface, 0, "SW_dw_dir_nir")->setVal(reference_radiation[1]);
    get_slm_data(land_surface, 0, "SW_dw_dif_vis")->setVal(reference_radiation[2]);
    get_slm_data(land_surface, 0, "SW_dw_dif_nir")->setVal(reference_radiation[3]);
    get_slm_data(land_surface, 0, "LW_dw")->setVal(reference_radiation[4]);
    get_slm_data(land_surface, 0, "cos_zenith")->setVal(reference_radiation[5]);

    const auto lsm_box = get_slm_data(land_surface, 0, "tsurf")->boxArray();
    const auto lsm_dm = get_slm_data(land_surface, 0, "tsurf")->DistributionMap();
    const auto surface_box = make_surface_box_array(lsm_box);
    amrex::Vector<amrex::Vector<std::unique_ptr<amrex::MultiFab>>> sst(1);
    amrex::Vector<amrex::Vector<std::unique_ptr<amrex::iMultiFab>>> landmask(1);
    sst[0].resize(1);
    landmask[0].resize(1);
    sst[0][0] = std::make_unique<amrex::MultiFab>(surface_box, lsm_dm, 1, amrex::IntVect(0));
    landmask[0][0] = std::make_unique<amrex::iMultiFab>(surface_box, lsm_dm, 1, amrex::IntVect(0));
    sst[0][0]->setVal(292.39);
    landmask[0][0]->setVal(1);
    amrex::MultiFab precip(surface_box, lsm_dm, 1, amrex::IntVect(0));
    precip.setVal(0.0);

    state->cons.setVal(1.0);
    state->cons.setVal(amrex::Real(1.0), Rho_comp, 1);
    state->cons.setVal(amrex::Real(300.0), RhoTheta_comp, 1);
    state->cons.setVal(amrex::Real(0.01), RhoQ1_comp, 1);
    state->uvel.setVal(0.1);
    state->vvel.setVal(0.1);
    land_surface.Update_Lsm_Vars_Lev(0, state->cons, state->uvel, state->vvel);
    // First row of Exec/DevTests/LandSurfaceModel_SLM/sndref. Set these
    // after the ERF-state update because that update populates the reference
    // atmospheric fields.
    set_synthetic_reference_atmosphere(
        land_surface, amrex::Real(857.0), amrex::Real(292.75),
        amrex::Real(15.541), amrex::Real(2.287272), amrex::Real(0.878001));
    land_surface.set_LSM_terrain_inputs(0, sst, landmask);
    land_surface.set_LSM_precip_input(0, &precip);
    land_surface.Advance(0, amrex::Real(10.0), amrex::Real(10.0), amrex::Real(0.0));

    EXPECT_NEAR(get_slm_data(land_surface, 0, "ref_t")->max(0),
                amrex::Real(292.75), amrex::Real(1.0e-12));
    EXPECT_NEAR(get_slm_data(land_surface, 0, "ref_q")->max(0),
                amrex::Real(0.015541), amrex::Real(1.0e-12));
    EXPECT_NEAR(get_slm_data(land_surface, 0, "ref_p")->max(0),
                amrex::Real(857.0), amrex::Real(1.0e-12));
    EXPECT_FALSE(get_slm_data(land_surface, 0, "tsurf")->contains_nan());
    EXPECT_FALSE(land_surface.Get_Flux_Ptr(0, LsmFlux_SLM::t_flux)->contains_nan());
}

TEST_F(SLMInterfaceTest, ConsumesReferenceStyleSplitRadiationFields)
{
    // Daytime swlw row at t = 209.333328: swdn = 401 W/m2, lwdn = 352 W/m2.
    // SLM's reference-input convention is 70% direct, 30% diffuse, with
    // each component split equally between visible and NIR.
    const amrex::Real swdn = 401.0;
    const std::array<amrex::Real, 6> expected_radiation{
        amrex::Real(0.5) * (swdn * amrex::Real(0.7)),
        amrex::Real(0.5) * (swdn * amrex::Real(0.7)),
        amrex::Real(0.5) * (swdn * amrex::Real(0.3)),
        amrex::Real(0.5) * (swdn * amrex::Real(0.3)),
        amrex::Real(352.0), amrex::Real(1.0)};
    const auto result = run_synthetic_radiation_case(solver_choice, expected_radiation);

    EXPECT_TRUE(std::isfinite(result.soil_temperature));
    EXPECT_TRUE(std::isfinite(result.surface_temperature));
    for (int comp = 0; comp < 6; ++comp) {
        EXPECT_NEAR(result.radiation_fields[comp], expected_radiation[comp],
                    amrex::Real(1.0e-12));
    }
}

TEST_F(SLMInterfaceTest, AcceptsSyntheticSSTAndPrecipitationRow)
{
    // First sfc row: SST = 292.390 K and precipitation = 0 mm/s.
    // A later file value is used to cover positive precipitation as well.
    const auto lsm_box = get_slm_data(land_surface, 0, "tsurf")->boxArray();
    const auto lsm_dm = get_slm_data(land_surface, 0, "tsurf")->DistributionMap();
    const auto surface_box = make_surface_box_array(lsm_box);
    amrex::Vector<amrex::Vector<std::unique_ptr<amrex::MultiFab>>> sst(1);
    amrex::Vector<amrex::Vector<std::unique_ptr<amrex::iMultiFab>>> landmask(1);
    sst[0].resize(1);
    landmask[0].resize(1);
    sst[0][0] = std::make_unique<amrex::MultiFab>(surface_box, lsm_dm, 1, amrex::IntVect(0));
    landmask[0][0] = std::make_unique<amrex::iMultiFab>(surface_box, lsm_dm, 1, amrex::IntVect(0));
    sst[0][0]->setVal(292.39);
    landmask[0][0]->setVal(1);
    amrex::MultiFab precip(surface_box, lsm_dm, 1, amrex::IntVect(0));
    precip.setVal(0.001157);

    land_surface.set_LSM_terrain_inputs(0, sst, landmask);
    land_surface.set_LSM_precip_input(0, &precip);

    // Terrain SST is stored on the SLM atmosphere/interface staging plane.
    auto* tsurf_mf = get_slm_data(land_surface, 0, "tsurf");
    auto* precip_mf = get_slm_data(land_surface, 0, "ref_precip");
    amrex::MFIter mfi(*tsurf_mf, amrex::TilingIfNotGPU());
    ASSERT_TRUE(mfi.isValid());
    const auto tsurf_arr = tsurf_mf->const_array(mfi);
    const auto precip_arr = precip_mf->const_array(mfi);
    EXPECT_NEAR(tsurf_arr(0, 0, 0), amrex::Real(292.39), amrex::Real(1.0e-12));
    EXPECT_NEAR(precip_arr(0, 0, -1), amrex::Real(0.001157), amrex::Real(1.0e-12));
}

TEST_F(SLMInterfaceTest, FirstAdvancePreservesLandSurfaceInterface)
{
    const auto lsm_box = get_slm_data(land_surface, 0, "tsurf")->boxArray();
    const auto lsm_dm = get_slm_data(land_surface, 0, "tsurf")->DistributionMap();
    const auto surface_box = make_surface_box_array(lsm_box);

    amrex::Vector<amrex::Vector<std::unique_ptr<amrex::MultiFab>>> sst(1);
    amrex::Vector<amrex::Vector<std::unique_ptr<amrex::iMultiFab>>> landmask(1);
    sst[0].resize(1);
    landmask[0].resize(1);
    sst[0][0] = std::make_unique<amrex::MultiFab>(surface_box, lsm_dm, 1, amrex::IntVect(0));
    landmask[0][0] = std::make_unique<amrex::iMultiFab>(surface_box, lsm_dm, 1, amrex::IntVect(0));
    sst[0][0]->setVal(301.0);
    landmask[0][0]->setVal(1);

    amrex::MultiFab precip(surface_box, lsm_dm, 1, amrex::IntVect(0));
    precip.setVal(0.0);

    land_surface.Update_Lsm_Vars_Lev(0, state->cons, state->uvel, state->vvel);
    land_surface.set_LSM_terrain_inputs(0, sst, landmask);
    land_surface.set_LSM_precip_input(0, &precip);
    land_surface.Advance(0, amrex::Real(0.3), amrex::Real(0.3), amrex::Real(0.0));
    land_surface.Update_State_Vars_Lev(0, state->cons);

    EXPECT_FALSE(state->cons.contains_nan(0, state->cons.nComp(), 0));
    EXPECT_NE(get_slm_data(land_surface, 0, "tsurf"), nullptr);
    EXPECT_NE(land_surface.Get_Flux_Ptr(0, LsmFlux_SLM::t_flux), nullptr);
}

TEST_F(SLMInterfaceTest, AdvancesMixedLandAndNonLandCells)
{
    const auto lsm_box = get_slm_data(land_surface, 0, "tsurf")->boxArray();
    const auto lsm_dm = get_slm_data(land_surface, 0, "tsurf")->DistributionMap();
    const auto surface_box = make_surface_box_array(lsm_box);

    amrex::Vector<amrex::Vector<std::unique_ptr<amrex::MultiFab>>> sst(1);
    amrex::Vector<amrex::Vector<std::unique_ptr<amrex::iMultiFab>>> landmask(1);
    sst[0].resize(1);
    landmask[0].resize(1);
    sst[0][0] = std::make_unique<amrex::MultiFab>(surface_box, lsm_dm, 1, amrex::IntVect(0));
    landmask[0][0] = std::make_unique<amrex::iMultiFab>(surface_box, lsm_dm, 1, amrex::IntVect(0));
    sst[0][0]->setVal(301.0);

    for (amrex::MFIter mfi(*landmask[0][0], amrex::TilingIfNotGPU()); mfi.isValid(); ++mfi) {
        auto mask = landmask[0][0]->array(mfi);
        const auto& bx = mfi.validbox();
        amrex::ParallelFor(bx, [=] AMREX_GPU_DEVICE (int i, int j, int k) noexcept {
            mask(i, j, k) = ((i + j) % 2 == 0) ? 1 : 0;
        });
    }

    amrex::MultiFab precip(surface_box, lsm_dm, 1, amrex::IntVect(0));
    precip.setVal(0.0);
    land_surface.Update_Lsm_Vars_Lev(0, state->cons, state->uvel, state->vvel);
    land_surface.set_LSM_terrain_inputs(0, sst, landmask);
    land_surface.set_LSM_precip_input(0, &precip);
    land_surface.Advance(0, amrex::Real(0.3), amrex::Real(0.3), amrex::Real(0.0));

    const auto* t_flux = land_surface.Get_Flux_Ptr(0, LsmFlux_SLM::t_flux);
    const auto* q_flux = land_surface.Get_Flux_Ptr(0, LsmFlux_SLM::q_flux);
    const auto* olen = land_surface.Get_Flux_Ptr(0, LsmFlux_SLM::olen);
    ASSERT_FALSE(t_flux->contains_nan());
    ASSERT_FALSE(q_flux->contains_nan());
    ASSERT_FALSE(olen->contains_nan());

    for (amrex::MFIter mfi(*t_flux, amrex::TilingIfNotGPU()); mfi.isValid(); ++mfi) {
        const auto flux = t_flux->const_array(mfi);
        const auto mask = landmask[0][0]->const_array(mfi);
        const auto& bx = mfi.validbox();
        amrex::Loop(bx, [=] (int i, int j, int k) noexcept {
            if (mask(i, j, k) == 0) {
                EXPECT_EQ(flux(i, j, k), amrex::Real(0.0));
            }
        });
    }
}

TEST_F(SLMInterfaceTest, AcceptsNonIntegralAndRestartLikeTimes)
{
    const amrex::Real dt_values[] = {amrex::Real(0.3), amrex::Real(1.25)};
    const amrex::Real time_values[] = {amrex::Real(0.3), amrex::Real(17.75)};
    const amrex::Real start_values[] = {amrex::Real(0.0), amrex::Real(12.5)};

    for (int i = 0; i < 2; ++i) {
        land_surface.Update_Lsm_Vars_Lev(0, state->cons, state->uvel, state->vvel);
        land_surface.Advance(0, dt_values[i], time_values[i], start_values[i]);
        land_surface.Update_State_Vars_Lev(0, state->cons);
        EXPECT_FALSE(state->cons.contains_nan(0, state->cons.nComp(), 0));
    }
}

TEST_F(SLMInterfaceTest, CopiesSyntheticTerrainAndPrecipitationInputs)
{
    const auto lsm_box = get_slm_data(land_surface, 0, "tsurf")->boxArray();
    const auto lsm_dm = get_slm_data(land_surface, 0, "tsurf")->DistributionMap();
    const auto surface_box = make_surface_box_array(lsm_box);

    amrex::Vector<amrex::Vector<std::unique_ptr<amrex::MultiFab>>> sst(1);
    amrex::Vector<amrex::Vector<std::unique_ptr<amrex::iMultiFab>>> landmask(1);
    sst[0].resize(1);
    landmask[0].resize(1);
    sst[0][0] = std::make_unique<amrex::MultiFab>(surface_box, lsm_dm, 1, amrex::IntVect(0));
    landmask[0][0] = std::make_unique<amrex::iMultiFab>(surface_box, lsm_dm, 1, amrex::IntVect(0));
    sst[0][0]->setVal(301.25);
    landmask[0][0]->setVal(0);

    amrex::MultiFab precip(surface_box, lsm_dm, 1, amrex::IntVect(0));
    precip.setVal(4.5);

    land_surface.set_LSM_terrain_inputs(0, sst, landmask);
    land_surface.set_LSM_precip_input(0, &precip);

    auto* tsurf_mf = get_slm_data(land_surface, 0, "tsurf");
    auto* precip_mf = get_slm_data(land_surface, 0, "ref_precip");
    amrex::MFIter mfi(*tsurf_mf, amrex::TilingIfNotGPU());
    ASSERT_TRUE(mfi.isValid());
    const auto tsurf_arr = tsurf_mf->const_array(mfi);
    const auto precip_arr = precip_mf->const_array(mfi);
    EXPECT_NEAR(tsurf_arr(0, 0, 0), amrex::Real(301.25), amrex::Real(1.0e-12));
    EXPECT_NEAR(precip_arr(0, 0, -1), amrex::Real(4.5), amrex::Real(1.0e-12));
}

TEST_F(SLMInterfaceTest, InitializesIndependentMultilevelSLMState)
{
    const amrex::Box coarse_box(
        amrex::IntVect(AMREX_D_DECL(0, 0, 0)),
        amrex::IntVect(AMREX_D_DECL(1, 1, 7)));
    const amrex::Box fine_box(
        amrex::IntVect(AMREX_D_DECL(0, 0, 0)),
        amrex::IntVect(AMREX_D_DECL(3, 3, 7)));

    SLMTestState coarse(coarse_box);
    SLMTestState fine(fine_box);
    LandSurface multilevel;
    SolverChoice& sc = solver_choice;
    multilevel.ReSize(2);
    multilevel.SetModel<SLM>();
    multilevel.Define(0, sc);
    multilevel.Define(1, sc);

    std::unique_ptr<amrex::MultiFab> coarse_z;
    std::unique_ptr<amrex::MultiFab> fine_z;
    multilevel.Init(0, coarse.cons, coarse.uvel, coarse.vvel,
                    coarse.geom, amrex::Real(0.3), coarse_z);
    multilevel.Init(1, fine.cons, fine.uvel, fine.vvel,
                    fine.geom, amrex::Real(0.3), fine_z);

    ASSERT_NE(get_slm_data(multilevel, 0, "tsurf"), nullptr);
    ASSERT_NE(get_slm_data(multilevel, 1, "tsurf"), nullptr);
    ASSERT_NE(multilevel.Get_Flux_Ptr(0, LsmFlux_SLM::t_flux), nullptr);
    ASSERT_NE(multilevel.Get_Flux_Ptr(1, LsmFlux_SLM::t_flux), nullptr);
    EXPECT_EQ(multilevel.Get_Lsm_Geom(0).Domain().length(0), 2);
    EXPECT_EQ(multilevel.Get_Lsm_Geom(1).Domain().length(0), 4);
    EXPECT_EQ(multilevel.Get_DataName(0), "tsurf");
    EXPECT_EQ(multilevel.Get_FluxName(0), "t_flux");

    for (int lev = 0; lev < 2; ++lev) {
        for (int var = 0; var < multilevel.Get_Data_Size(); ++var) {
            ASSERT_NE(multilevel.Get_Data_Ptr(lev, var), nullptr);
            EXPECT_EQ(multilevel.Get_DataIdx(lev, multilevel.Get_DataName(var)), var);
        }
        for (int var = 0; var < multilevel.Get_Flux_Size(); ++var) {
            ASSERT_NE(multilevel.Get_Flux_Ptr(lev, var), nullptr);
            EXPECT_EQ(multilevel.Get_FluxIdx(lev, multilevel.Get_FluxName(var)), var);
        }
        EXPECT_EQ(get_slm_data(multilevel, lev, "tsurf")->nGrowVect(),
                  amrex::IntVect(1));
        EXPECT_EQ(multilevel.Get_Flux_Ptr(lev, LsmFlux_SLM::t_flux)->nGrowVect(),
                  amrex::IntVect(1, 1, 0));
    }

    get_slm_data(multilevel, 0, "tsurf")->setVal(271.0);
    EXPECT_NE(get_slm_data(multilevel, 1, "tsurf")->max(0),
              amrex::Real(271.0));
}

TEST_F(SLMInterfaceTest, TransfersRawInputStateFromLevelZero)
{
    initialize_multilevel(2, InitType::WRFInput);
    seed_transfer_state(0, amrex::Real(285.0), amrex::Real(290.0),
                        amrex::Real(0.4), amrex::Real(3.0));

    land_surface.Initialize_From_Source(1, 0, amrex::IntVect(2),
                                         LSMTransferMode::RawInput);

    EXPECT_NEAR(get_slm_data(land_surface, 1, "tsurf")->max(0),
                amrex::Real(285.0), amrex::Real(1.0e-12));
    EXPECT_NEAR(get_slm_data(land_surface, 1, "tsoil")->max(0),
                amrex::Real(290.0), amrex::Real(1.0e-12));
    EXPECT_NEAR(get_slm_data(land_surface, 1, "wsoil")->max(0),
                amrex::Real(0.4), amrex::Real(1.0e-12));
    EXPECT_NEAR(land_surface.Get_Flux_Ptr(1, LsmFlux_SLM::t_flux)->max(0),
                amrex::Real(3.0), amrex::Real(1.0e-12));
}

TEST_F(SLMInterfaceTest, TransfersProcessedStateFromLevelZero)
{
    initialize_multilevel(2);
    seed_transfer_state(0, amrex::Real(286.0), amrex::Real(291.0),
                        amrex::Real(0.45), amrex::Real(4.0));

    land_surface.Initialize_From_Source(1, 0, amrex::IntVect(2),
                                         LSMTransferMode::ProcessedState);

    EXPECT_NEAR(get_slm_data(land_surface, 1, "tsurf")->max(0),
                amrex::Real(286.0), amrex::Real(1.0e-12));
    EXPECT_NEAR(get_slm_data(land_surface, 1, "tsoil")->max(0),
                amrex::Real(291.0), amrex::Real(1.0e-12));
    EXPECT_NEAR(get_slm_data(land_surface, 1, "wsoil")->max(0),
                amrex::Real(0.45), amrex::Real(1.0e-12));
    EXPECT_NEAR(land_surface.Get_Flux_Ptr(1, LsmFlux_SLM::t_flux)->max(0),
                amrex::Real(4.0), amrex::Real(1.0e-12));
    EXPECT_FALSE(get_slm_data(land_surface, 1, "tsurf")->contains_nan());
    EXPECT_FALSE(land_surface.Get_Flux_Ptr(1, LsmFlux_SLM::olen)->contains_nan());
}

TEST_F(SLMInterfaceTest, UsesCumulativeAndParentSourceRatios)
{
    initialize_multilevel(3);
    seed_transfer_state(0, amrex::Real(287.0), amrex::Real(292.0),
                        amrex::Real(0.5), amrex::Real(5.0));
    seed_transfer_state(1, amrex::Real(288.0), amrex::Real(293.0),
                        amrex::Real(0.55), amrex::Real(6.0));
    set_surface_temperature_gradient(land_surface, 0, amrex::Real(287.0));
    set_surface_temperature_gradient(land_surface, 1, amrex::Real(288.0));

    land_surface.Initialize_From_Source(2, 0, amrex::IntVect(4),
                                         LSMTransferMode::ProcessedState);
    EXPECT_NEAR(surface_temperature_at(
                    2, amrex::IntVect(AMREX_D_DECL(4, 4, -1))),
                amrex::Real(290.0), amrex::Real(1.0e-12));
    EXPECT_NEAR(land_surface.Get_Flux_Ptr(2, LsmFlux_SLM::t_flux)->max(0),
                amrex::Real(5.0), amrex::Real(1.0e-12));

    land_surface.Initialize_From_Source(2, 1, amrex::IntVect(2),
                                         LSMTransferMode::ProcessedState);
    EXPECT_NEAR(surface_temperature_at(
                    2, amrex::IntVect(AMREX_D_DECL(6, 6, -1))),
                amrex::Real(297.0), amrex::Real(1.0e-12));
    EXPECT_NEAR(land_surface.Get_Flux_Ptr(2, LsmFlux_SLM::t_flux)->max(0),
                amrex::Real(6.0), amrex::Real(1.0e-12));
}

TEST_F(SLMInterfaceTest, RemakePreservesProcessedState)
{
    seed_transfer_state(0, amrex::Real(289.0), amrex::Real(294.0),
                        amrex::Real(0.6), amrex::Real(7.0));

    std::unique_ptr<amrex::MultiFab> z_phys_nd;
    amrex::Vector<amrex::BCRec> domain_bcs_type = make_domain_bcs_type();
    amrex::IntVect ref_ratio(1);
    amrex::Vector<amrex::Vector<std::string>> nc_init_file;
    land_surface.Remake_Level(0, state->cons, state->uvel, state->vvel,
                              state->geom, state->geom, domain_bcs_type,
                              ref_ratio, amrex::Real(1.0), z_phys_nd,
                              nc_init_file);

    EXPECT_NEAR(get_slm_data(land_surface, 0, "tsurf")->max(0),
                amrex::Real(289.0), amrex::Real(1.0e-12));
    EXPECT_NEAR(get_slm_data(land_surface, 0, "tsoil")->max(0),
                amrex::Real(294.0), amrex::Real(1.0e-12));
    EXPECT_NEAR(get_slm_data(land_surface, 0, "wsoil")->max(0),
                amrex::Real(0.6), amrex::Real(1.0e-12));
    EXPECT_NEAR(land_surface.Get_Flux_Ptr(0, LsmFlux_SLM::t_flux)->max(0),
                amrex::Real(7.0), amrex::Real(1.0e-12));
}

TEST_F(SLMInterfaceTest, CheckpointWritesAndReadsSyntheticState)
{
    const std::filesystem::path checkpoint =
        std::filesystem::temp_directory_path() / "erf_slm_unit_checkpoint";
    std::filesystem::remove_all(checkpoint);
    std::filesystem::create_directories(checkpoint / "Level_0");

    land_surface.Update_Lsm_Vars_Lev(0, state->cons, state->uvel, state->vvel);
    land_surface.Advance(0, amrex::Real(0.3), amrex::Real(0.3), amrex::Real(0.0));
    land_surface.WriteCheckpoint(0, checkpoint.string());

    const std::filesystem::path header = checkpoint / "Level_0" / "SLM_Header";
    ASSERT_TRUE(std::filesystem::exists(header));
    {
        std::ifstream header_stream(header);
        std::string header_text((std::istreambuf_iterator<char>(header_stream)),
                                std::istreambuf_iterator<char>());
        EXPECT_NE(header_text.find("Checkpoint file for SLM"), std::string::npos);
    }

    land_surface.ReadCheckpoint(0, checkpoint.string());
    EXPECT_FALSE(get_slm_data(land_surface, 0, "tsurf")->contains_nan());
    EXPECT_FALSE(land_surface.Get_Flux_Ptr(0, LsmFlux_SLM::t_flux)->contains_nan());

    std::filesystem::remove_all(checkpoint);
}

TEST_F(SLMInterfaceTest, ReadsDirectSLMCheckpointIntoFreshSLM)
{
    const std::filesystem::path checkpoint =
        std::filesystem::temp_directory_path() / "erf_slm_unit_checkpoint_fresh";
    std::filesystem::remove_all(checkpoint);
    std::filesystem::create_directories(checkpoint / "Level_0");

    land_surface.Update_Lsm_Vars_Lev(0, state->cons, state->uvel, state->vvel);
    land_surface.Advance(0, amrex::Real(0.3), amrex::Real(0.3), amrex::Real(0.0));
    land_surface.WriteCheckpoint(0, checkpoint.string());

    const std::filesystem::path header = checkpoint / "Level_0" / "SLM_Header";
    ASSERT_TRUE(std::filesystem::exists(header));

    LandSurface restored;
    restored.ReSize(1);
    restored.SetModel<SLM>();
    restored.Define(0, solver_choice);
    std::unique_ptr<amrex::MultiFab> z_phys_nd;
    restored.Init(0, state->cons, state->uvel, state->vvel,
                  state->geom, amrex::Real(1.0), z_phys_nd);
    restored.ReadCheckpoint(0, checkpoint.string());

    // ERF writes and restores the public LsmData*/LsmFlux* files around this
    // direct SLM checkpoint call.  This test intentionally checks only the
    // direct SLM API: it must read its own checkpoint into a fresh object and
    // leave the exposed fields finite and geometrically valid.
    EXPECT_EQ(restored.Get_Lsm_Geom(0).Domain(),
              land_surface.Get_Lsm_Geom(0).Domain());
    for (int var = 0; var < restored.Get_Data_Size(); ++var) {
        EXPECT_FALSE(restored.Get_Data_Ptr(0, var)->contains_nan());
    }
    for (int var = 0; var < restored.Get_Flux_Size(); ++var) {
        EXPECT_FALSE(restored.Get_Flux_Ptr(0, var)->contains_nan());
    }

    std::filesystem::remove_all(checkpoint);
}

} // namespace
