/**
 * @file ERF_GTestFloodEngineExample.cpp
 * @brief A unit test that a flood engine directory builds into ERF's unit-test executable
 *
 * ERFFloodEngine.cmake adds this file through erf_flood_engine_add_unit_tests, so it is
 * compiled into erf_unit_tests with ERF's test main, and the example engine it tests is
 * already in that executable through erf_flood_engine_add. The test creates the engine by
 * its registered name and checks what it declares to the coupling.
 */
#include <gtest/gtest.h>

#include <ERF_FloodEngineRegistry.H>

#include <string>
#include <vector>

TEST(FloodEngineExample, RegisteredAndDeclaresItsField)
{
    auto engine = make_flood_engine("example", FloodCouplingParams{});
    ASSERT_NE(engine, nullptr) << "the example engine is not compiled into this executable";
    EXPECT_STREQ(engine->name(), "example");
    EXPECT_TRUE(engine->wants_surface_precip());
    EXPECT_EQ(engine->plot_var_names(), std::vector<std::string>{"flood_depth"});
}
