#ifdef ERF_USE_NETCDF

#include <filesystem>
#include <string>
#include <type_traits>
#include <utility>

#include <gtest/gtest.h>

#include "ERF_NCInterface.H"
#include "../ERF_GTestTempDir.H"

namespace {

// The file owns its NetCDF ID and closes it on destruction, so a copy would close it twice
static_assert(!std::is_copy_constructible_v<ncutils::NCFile>);
static_assert(!std::is_copy_assignable_v<ncutils::NCFile>);
static_assert(std::is_move_constructible_v<ncutils::NCFile>);

// Writes a NetCDF-4 file with two dimensions in the root group and one in group "child".
// NetCDF-4 dimension IDs are unique across the file, so the child's dimension is not ID 0
std::string write_grouped_file ()
{
    const std::string fname = erf_gtest_temp_path("erf_unit_nc_groups").string() + ".nc";
    auto ncf = ncutils::NCFile::create(fname, NC_CLOBBER | NC_NETCDF4);
    ncf.def_dim("root_a", 2);
    ncf.def_dim("root_b", 3);
    int child_id = -1;
    EXPECT_EQ(nc_def_grp(ncf.ncid, "child", &child_id), NC_NOERR);
    auto child = ncf.group("child");
    child.def_dim("child_c", 5);
    child.def_var("child_var", NC_DOUBLE, {"child_c"});
    ncf.close();
    return fname;
}

} // namespace

TEST(NCInterface, AllDimsInGroupUsesTheGroupDimensionIDs)
{
    if (!amrex::ParallelDescriptor::IOProcessor()) { GTEST_SKIP(); }
    const std::string fname = write_grouped_file();
    {
        auto ncf = ncutils::NCFile::open(fname, NC_NOWRITE);

        const auto root_dims = ncf.all_dims();
        ASSERT_EQ(root_dims.size(), 2u);
        EXPECT_EQ(root_dims[0].name(), "root_a");
        EXPECT_EQ(root_dims[1].name(), "root_b");

        const auto child_dims = ncf.group("child").all_dims();
        ASSERT_EQ(child_dims.size(), 1u);
        EXPECT_EQ(child_dims[0].name(), "child_c");
        EXPECT_EQ(child_dims[0].len(), 5u);

        const auto child_vars = ncf.group("child").all_vars();
        ASSERT_EQ(child_vars.size(), 1u);
        EXPECT_EQ(child_vars[0].name(), "child_var");
    }
    std::filesystem::remove(fname);
}

TEST(NCInterface, MovedFromFileDoesNotCloseTheHandle)
{
    if (!amrex::ParallelDescriptor::IOProcessor()) { GTEST_SKIP(); }
    const std::string fname = write_grouped_file();
    {
        auto first = ncutils::NCFile::open(fname, NC_NOWRITE);
        ncutils::NCFile second(std::move(first));
        EXPECT_EQ(second.dim("root_b").len(), 3u);
        second.close();
        // Leaving the scope destroys both; closing the ID again would abort
    }
    std::filesystem::remove(fname);
}

#endif // ERF_USE_NETCDF
