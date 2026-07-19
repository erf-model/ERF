#include <gtest/gtest.h>

#include "ERF_PlotfileStream.H"

namespace
{

enum class TestFormat { Amrex, Netcdf };

TEST (PlotfileStream, FinalStreamFormatsRemainIndependent)
{
    EXPECT_EQ(erf_plotfile::stream_format(1, TestFormat::Amrex, TestFormat::Netcdf),
              TestFormat::Amrex);
    EXPECT_EQ(erf_plotfile::stream_format(2, TestFormat::Amrex, TestFormat::Netcdf),
              TestFormat::Netcdf);
}

} // namespace
