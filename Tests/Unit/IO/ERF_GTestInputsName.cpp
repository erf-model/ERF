// Contract of CheckForDuplicateInputs: a deck may be split into files included with
// FILE = <name> (the AMReX ParmParse include, which is not a parameter), but every key is set in
// one place across the deck and everything it includes. A key set twice in one file, in the deck
// and an included file, in two included files, or through a nested include is a duplicate. The
// scan is tested through erf_inputs_detail::scan_inputs, which reports duplicates instead of
// aborting.

#include <filesystem>
#include <fstream>
#include <string>
#include <unordered_map>
#include <vector>

#include <gtest/gtest.h>

#include <ERF_InputsName.H>

#include "../ERF_GTestTempDir.H"

namespace {

struct Deck {
    std::filesystem::path dir = erf_gtest_temp_path("erf_gtest_inputs_name");
    Deck () { std::filesystem::create_directories(dir); }
    ~Deck () { std::filesystem::remove_all(dir); }
    // a file of the deck, by its name in the scratch directory; FILE lines name included files the same way
    std::string write (const std::string& name, const std::string& text) const
    {
        std::ofstream(dir / name) << text;
        return (dir / name).string();
    }
    std::string path (const std::string& name) const { return (dir / name).string(); }
};

bool has_duplicates (const std::string& file)
{
    std::unordered_map<std::string, erf_inputs_detail::EntryInfo> seen;
    std::vector<std::string> open_files;
    return erf_inputs_detail::scan_inputs(file, seen, open_files);
}

} // namespace

TEST(InputsName, SeveralFileIncludesAreNotDuplicates)
{
    Deck d;
    d.write("flow.inputs", "erf.fixed_dt = 0.5\nerf.use_gravity = true\n");
    d.write("network.inputs", "erf.conductors.spans = L1  # generated\n");
    const std::string deck = d.write("inputs", "# shared settings and a generated block\nFILE = " + d.path("flow.inputs") +
                                     "\nmax_step = 10\nFILE = \"" + d.path("network.inputs") + "\"   # quoted\nerf.v = 1\n");
    EXPECT_FALSE(has_duplicates(deck));
    // and the public check returns
    EXPECT_NO_FATAL_FAILURE(CheckForDuplicateInputs(deck));
}

TEST(InputsName, AKeySetTwiceAnywhereInTheDeckIsADuplicate)
{
    Deck d;
    d.write("a.inputs", "erf.fixed_dt = 0.5\n");
    d.write("b.inputs", "erf.fixed_dt = 0.25\n");
    d.write("inner.inputs", "max_step = 20\n");
    d.write("outer.inputs", "FILE = " + d.path("inner.inputs") + "\n");
    // in one file, as before
    EXPECT_TRUE(has_duplicates(d.write("one", "max_step = 10\nerf.v = 1\nmax_step = 20\n")));
    // in the deck and a file it includes
    EXPECT_TRUE(has_duplicates(d.write("deck_and_include", "erf.fixed_dt = 1.0\nFILE = " + d.path("a.inputs") + "\n")));
    // in two included files
    EXPECT_TRUE(has_duplicates(d.write("two_includes", "FILE = " + d.path("a.inputs") + "\nFILE = " + d.path("b.inputs") + "\n")));
    // through a nested include
    EXPECT_TRUE(has_duplicates(d.write("nested", "max_step = 10\nFILE = " + d.path("outer.inputs") + "\n")));
    // the same file included twice sets its keys twice
    EXPECT_TRUE(has_duplicates(d.write("twice", "FILE = " + d.path("a.inputs") + "\nFILE = " + d.path("a.inputs") + "\n")));
    // a comment that mentions a key is not a setting
    EXPECT_FALSE(has_duplicates(d.write("comment", "max_step = 10   # was max_step = 5\n# max_step = 7\n")));
}
