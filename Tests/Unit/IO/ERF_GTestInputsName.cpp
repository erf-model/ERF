// Contract of CheckForDuplicateInputs: a deck may be split into files included with
// FILE = <name> (the AMReX ParmParse include, which is not a parameter), and a file may re-set a
// key that a file it included earlier set (base + override; ParmParse keeps the last value). A key
// set twice in one file, in two files included side by side, or in a file and then again in a file
// it includes afterwards is a duplicate. The scan reads files as ParmParse does: names relative to
// the working directory with $AMREX_INPUTS_FILE_PREFIX in front, FILE with several values is a
// parameter, UNSET drops keys, [name] tables prefix their keys, and #if regions that are off are
// skipped. It is tested through erf_inputs_detail::scan_inputs, which reports duplicates and
// errors instead of aborting.

#include <cstdlib>
#include <filesystem>
#include <fstream>
#include <iostream>
#include <sstream>
#include <string>

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
        std::filesystem::create_directories((dir / name).parent_path());
        std::ofstream(dir / name) << text;
        return (dir / name).string();
    }
    std::string path (const std::string& name) const { return (dir / name).string(); }
};

// runs a test body in another working directory, as ERF runs a deck from its own directory
struct InDirectory {
    std::filesystem::path old = std::filesystem::current_path();
    explicit InDirectory (const std::filesystem::path& d) { std::filesystem::current_path(d); }
    ~InDirectory () { std::filesystem::current_path(old); }
};

void set_prefix (const char* value)
{
#ifdef _WIN32
    _putenv_s("AMREX_INPUTS_FILE_PREFIX", value == nullptr ? "" : value);
#else
    if (value == nullptr) { unsetenv("AMREX_INPUTS_FILE_PREFIX"); }
    else                  { setenv("AMREX_INPUTS_FILE_PREFIX", value, 1); }
#endif
}

erf_inputs_detail::ScanState scan (const std::string& file)
{
    erf_inputs_detail::ScanState state;
    erf_inputs_detail::scan_inputs(file, state);
    return state;
}

bool has_duplicates (const std::string& file)
{
    auto state = scan(file);
    EXPECT_EQ(state.error, "") << file;
    return state.found_duplicates;
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
    // and the public check returns (an abort there ends the test binary, a failure ctest reports)
    CheckForDuplicateInputs(deck);
}

TEST(InputsName, AKeySetTwiceIsADuplicate)
{
    Deck d;
    d.write("a.inputs", "erf.fixed_dt = 0.5\n");
    d.write("b.inputs", "erf.fixed_dt = 0.25\n");
    d.write("inner.inputs", "max_step = 20\n");
    d.write("outer.inputs", "FILE = " + d.path("inner.inputs") + "\n");
    // in one file, as before
    EXPECT_TRUE(has_duplicates(d.write("one", "max_step = 10\nerf.v = 1\nmax_step = 20\n")));
    // in a file and then in a file it includes afterwards: the include would win silently
    EXPECT_TRUE(has_duplicates(d.write("deck_then_include", "erf.fixed_dt = 1.0\nFILE = " + d.path("a.inputs") + "\n")));
    // in two files included side by side
    EXPECT_TRUE(has_duplicates(d.write("two_includes", "FILE = " + d.path("a.inputs") + "\nFILE = " + d.path("b.inputs") + "\n")));
    // in a file and then through a nested include
    EXPECT_TRUE(has_duplicates(d.write("nested", "max_step = 10\nFILE = " + d.path("outer.inputs") + "\n")));
    // the same file included twice sets its keys twice
    EXPECT_TRUE(has_duplicates(d.write("twice", "FILE = " + d.path("a.inputs") + "\nFILE = " + d.path("a.inputs") + "\n")));
    // a comment that mentions a key is not a setting
    EXPECT_FALSE(has_duplicates(d.write("comment", "max_step = 10   # was max_step = 5\n# max_step = 7\n")));
}

TEST(InputsName, AFileMayOverrideWhatItIncludedEarlier)
{
    Deck d;
    d.write("base.inputs", "erf.moisture_model = SAM\nmax_step = 100\n");
    d.write("mid.inputs", "FILE = " + d.path("base.inputs") + "\nmax_step = 50\n");
    // the base + override idiom of the Bubble WDM6 and SEB decks
    EXPECT_FALSE(has_duplicates(d.write("override", "FILE = " + d.path("base.inputs") + "\nerf.moisture_model = WDM6\n")));
    // through two levels: the deck overrides what the middle file included and overrode
    EXPECT_FALSE(has_duplicates(d.write("deep", "FILE = " + d.path("mid.inputs") + "\nmax_step = 10\nerf.moisture_model = WDM6\n")));
    // overriding the same key twice in one file is a duplicate in that file
    EXPECT_TRUE(has_duplicates(d.write("twice", "FILE = " + d.path("base.inputs") + "\nmax_step = 5\nmax_step = 6\n")));
    // a second include that sets the key again is not an override
    EXPECT_TRUE(has_duplicates(d.write("sibling", "FILE = " + d.path("mid.inputs") + "\nFILE = " + d.path("base.inputs") + "\n")));
    // nor is one whose first setting lies deeper, on another branch: the deck includes A (which
    // includes B, setting the key) and then C, a sibling of A, sets it again
    d.write("B.inputs", "erf.moisture_model = SAM\n");
    d.write("A.inputs", "FILE = " + d.path("B.inputs") + "\n");
    d.write("C.inputs", "erf.moisture_model = WDM6\n");
    EXPECT_TRUE(has_duplicates(d.write("branches", "FILE = " + d.path("A.inputs") + "\nFILE = " + d.path("C.inputs") + "\n")));
}

TEST(InputsName, UnsetIsADirectiveThatDropsKeys)
{
    Deck d;
    d.write("base.inputs", "erf.fixed_dt = 0.5\nerf.les_type = Smagorinsky\n");
    // UNSET lines are not keys, however many there are
    EXPECT_FALSE(has_duplicates(d.write("unsets", "FILE = " + d.path("base.inputs") + "\nUNSET = erf.fixed_dt\nUNSET = erf.les_type\n")));
    // a key dropped with UNSET may be set again, in the same file too; several keys on one line
    EXPECT_FALSE(has_duplicates(d.write("reset", "a = 1\nb = 2\nUNSET = a b\na = 3\nb = 4\n")));
    // a key not dropped is still a duplicate
    EXPECT_TRUE(has_duplicates(d.write("partial", "a = 1\nb = 2\nUNSET = a\na = 3\nb = 4\n")));
}

TEST(InputsName, FileWithSeveralValuesIsAParameter)
{
    Deck d;
    // ParmParse includes only for one value; FILE = a b is a parameter named FILE
    auto state = scan(d.write("param", "FILE = a b\n"));
    EXPECT_EQ(state.error, "");
    EXPECT_EQ(state.seen.count("FILE"), 1u);
    EXPECT_TRUE(has_duplicates(d.write("param_twice", "FILE = a b\nFILE = c d\n")));
    // one quoted value with a space is one include
    d.write("with space.inputs", "x = 1\n");
    EXPECT_FALSE(has_duplicates(d.write("quoted", "FILE = \"" + d.path("with space.inputs") + "\"\nx2 = 1\n")));
}

TEST(InputsName, RelativeIncludesAreOpenedFromTheWorkingDirectory)
{
    Deck d;
    d.write("flow.inputs", "erf.fixed_dt = 0.5\n");
    d.write("sub/inner.inputs", "max_step = 1\n");
    d.write("inputs", "FILE = flow.inputs\nFILE = sub/inner.inputs\nerf.fixed_dt = 0.25\n");
    InDirectory here(d.dir);
    EXPECT_FALSE(has_duplicates("inputs"));
    EXPECT_EQ(scan("inputs").seen.at("erf.fixed_dt").value, "0.25");
    // a name that is not there is an error, not a crash
    d.write("missing", "FILE = nowhere.inputs\n");
    EXPECT_NE(scan("missing").error.find("Could not open inputs file: nowhere.inputs"), std::string::npos);
}

TEST(InputsName, AnIncludeLoopIsFoundWhateverTheSpelling)
{
    Deck d;
    d.write("self", "a = 1\nFILE = ./self\n");
    d.write("a", "FILE = sub/b\n");
    // names are opened from the working directory, not from the including file's directory
    d.write("sub/b", "FILE = sub/../a\n");
    d.write("direct", "FILE = direct\n");
    InDirectory here(d.dir);
    EXPECT_NE(scan("self").error.find("includes itself"), std::string::npos);
    EXPECT_NE(scan("a").error.find("includes itself"), std::string::npos);
    EXPECT_NE(scan("direct").error.find("includes itself"), std::string::npos);
    // an error deep in the includes leaves no include open: the state can scan again
    auto state = scan("a");
    EXPECT_TRUE(state.includes.empty());
    EXPECT_TRUE(state.open_files.empty());
}

TEST(InputsName, TheInputsFilePrefixIsUsedForEveryFile)
{
    Deck d;
    d.write("flow.inputs", "erf.fixed_dt = 0.5\n");
    d.write("inputs", "FILE = flow.inputs\nmax_step = 10\n");
    // the deck directory without a trailing '/', which ParmParse adds
    set_prefix(d.dir.string().c_str());
    auto state = scan("inputs");
    set_prefix(nullptr);
    EXPECT_EQ(state.error, "");
    EXPECT_EQ(state.seen.count("erf.fixed_dt"), 1u);
    // compared as paths: ParmParse joins the prefix and the name with '/', and the scratch
    // directory uses the native separator, a backslash on Windows
    const std::string& file = state.seen.at("erf.fixed_dt").file;
    EXPECT_TRUE(std::filesystem::path(file) == d.dir / "flow.inputs") << file;
}

TEST(InputsName, TablesAndPreprocessorRegionsAreReadAsParmParseDoes)
{
    Deck d;
    // [name] tables: the same key in two tables is two keys, and name.key at the root is the same key
    EXPECT_FALSE(has_duplicates(d.write("tables", "[erf]\nfixed_dt = 1\n[other]\nfixed_dt = 2\n")));
    EXPECT_TRUE(has_duplicates(d.write("table_root", "erf.fixed_dt = 1\n[erf]\nfixed_dt = 2\n")));
    // an included file starts at the root, and the including file's table resumes after it
    d.write("root.inputs", "fixed_dt = 1\n");
    EXPECT_FALSE(has_duplicates(d.write("table_include", "[erf]\nfixed_dt = 2\nFILE = " + d.path("root.inputs") + "\nmax_step = 3\n")));
    EXPECT_EQ(scan(d.path("table_include")).seen.count("erf.max_step"), 1u);
    // only one branch of an #if region is read
    EXPECT_FALSE(has_duplicates(d.write("dim", "#if AMREX_SPACEDIM == 3\na = 1\n#elif AMREX_SPACEDIM == 2\na = 2\n#else\na = 3\n#endif\n")));
    EXPECT_FALSE(has_duplicates(d.write("gpu", "#ifdef AMREX_USE_GPU\nb = 1\n#else\nb = 2\n#endif\n")));
    // and a duplicate inside the branch that is read is still found
    EXPECT_TRUE(has_duplicates(d.write("dim_dup", "#if AMREX_SPACEDIM == 3\na = 1\na = 2\n#endif\n")));
    EXPECT_EQ(scan(d.path("dim")).seen.at("a").value, AMREX_SPACEDIM == 3 ? "1" : (AMREX_SPACEDIM == 2 ? "2" : "3"));
}

TEST(InputsName, TheMessageNamesBothPlacesOnSeparateLines)
{
    Deck d;
    d.write("repeat.inputs", "erf.fixed_dt = 0.01\n");
    const std::string deck = d.write("inputs", "erf.fixed_dt = 0.0005\nFILE = " + d.path("repeat.inputs") + "\n");
    std::ostringstream out;
    auto* old = std::cout.rdbuf(out.rdbuf());
    auto state = scan(deck);
    std::cout.rdbuf(old);
    EXPECT_TRUE(state.found_duplicates);
    if (amrex::ParallelDescriptor::IOProcessor()) {
        const std::string text = out.str();
        EXPECT_NE(text.find("    First line  : 1 of " + deck + "\n    Value       : 0.0005\n"), std::string::npos) << text;
        EXPECT_NE(text.find("    Second line : 1 of " + d.path("repeat.inputs") + "\n    Value       : 0.01\n"), std::string::npos) << text;
    }
}
