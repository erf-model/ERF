// Unit tests for reading erf.pbl_type (TurbChoice::init_params):
//  * YSU selects the YSU scheme (formerly YSUNew, which replaced the original YSU);
//  * the deprecated name YSUNew, in any case, is read as YSU, so that no other code
//    sees PBLType::YSUNew.

#include <string>

#include <AMReX_ParmParse.H>

#include <ERF_DataStruct.H>

#include <gtest/gtest.h>

namespace {

PBLType
read_pbl_type (const std::string& prefix, const std::string& value)
{
    // YSU requires the surface_layer BC at zlo, which init_params reads
    std::string zlo_type = "surface_layer";
    amrex::ParmParse("zlo").queryAdd("type", zlo_type);
    amrex::ParmParse pp(prefix);
    pp.add("pbl_type", value);
    TurbChoice tc;
    tc.init_params(0, 0, prefix);
    return tc.pbl_type;
}

} // namespace

TEST(PBLTypeInputs, YSUSelectsYSU)
{
    EXPECT_EQ(read_pbl_type("pbl_type_gtest_ysu", "YSU"), PBLType::YSU);
}

TEST(PBLTypeInputs, YSUNewIsReadAsYSU)
{
    EXPECT_EQ(read_pbl_type("pbl_type_gtest_ysunew", "YSUNew"), PBLType::YSU);
    EXPECT_EQ(read_pbl_type("pbl_type_gtest_ysunew_lower", "ysunew"), PBLType::YSU);
}
