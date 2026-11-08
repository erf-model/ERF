// The frame model against SubDyn: the reader on SubDyn's own input files, its refusals, and the
// stiffness of a lattice tower condensed to its peak joint compared with SubDyn's KBBt
// (Tests/test_files/FrameSubDynTower, from OpenFAST 5.0.0's subdyn_driver in double precision).
//
// FrameSubDyn.ReadsEveryTableOfASubDynFile: case B's counts, beam theory, divisions, member
//   types and spin, spring, interface joint and section values.
// FrameSubDyn.RefusesWhatTheFrameCannotModelNamingIt: cables, rigid links, tapered members,
//   FEMMod 2, NDiv 0, non-cantilever joints, unknown joints, a missing SSI file, a negative
//   modulus and a short table are each refused with a message naming the item.
// FrameSubDyn.RefusesSpringsMassesAndSupportsThatDoNotReadWhole: an SSI entry with an unknown name, an
//   SSI file without entries or with a non-finite mass, a support that restrains nothing, and a
//   concentrated-mass row of neither 5 nor 11 values or with a non-finite inertia are refused.
// FrameSubDyn.TheStiffnessAtThePeakIsSubDyns: for cases A and B (Euler-Bernoulli with arbitrary
//   sections; Timoshenko with circular, rectangular and spun arbitrary sections, two elements per
//   member and a coupled spring base) the 6 x 6 stiffness at the peak, and for case T (a 30 m
//   lattice tower with a 12 m cross-arm) at the cross-arm's centre, equals SubDyn's to the 7 digits it
//   prints.
// The suite also runs in the parallel executable: the reader is collective (the I/O rank reads
//   the file and broadcasts it), so every rank must get the same frame and the same refusals.

#include <gtest/gtest.h>

#include <AMReX_ParallelDescriptor.H>

#include <array>
#include <cmath>
#include <cstdio>
#include <filesystem>
#include <fstream>
#include <sstream>
#include <string>
#include <vector>

#include "ERF_Frame.H"
#include "../ERF_GTestTempDir.H"

using namespace erf_towers;

namespace {

std::string case_file (const std::string& c) { return std::string(ERF_FRAME_TEST_FILES) + "/case" + c + "/tower" + c + ".dat"; }

std::string slurp (const std::string& path)
{
    std::ifstream f(path);
    std::stringstream s;
    s << f.rdbuf();
    return s.str();
}

/**
 * Write case A with pieces of text replaced, and beside it the files given (name, contents), in a
 * directory of its own; read it, and return the reader's or the checks' message.
 */
std::string read_edits (const std::vector<std::pair<std::string, std::string>>& edits, const std::string& name,
                        const std::vector<std::pair<std::string, std::string>>& files = {})
{
    std::string text = slurp(case_file("A"));
    for (const auto& [from, to] : edits) {
        const auto p = text.find(from);
        EXPECT_NE(p, std::string::npos) << "the edit '" << from << "' does not apply to case A";
        if (p == std::string::npos) { return "edit not applied"; }
        text.replace(p, from.size(), to);
    }
    // one rank writes the files, every rank reads them (read_subdyn is collective: the I/O rank reads and broadcasts)
    const std::filesystem::path dir = erf_gtest_temp_path("frame_subdyn_" + name);
    const std::string path = (dir / "tower.dat").string();
    if (amrex::ParallelDescriptor::IOProcessor()) {
        std::filesystem::create_directories(dir);
        { std::ofstream f(path); f << text; }
        for (const auto& [fname, contents] : files) { std::ofstream f(dir / fname); f << contents; }
    }
    amrex::ParallelDescriptor::Barrier();
    FrameInputs in;
    std::string err = read_subdyn(path, in);
    if (err.empty()) { err = in.validate(); }
    amrex::ParallelDescriptor::Barrier();
    if (amrex::ParallelDescriptor::IOProcessor()) { std::filesystem::remove_all(dir); }
    return err;
}

/** read_edits() with one piece of text replaced. */
std::string read_edited (const std::string& from, const std::string& to, const std::string& name)
{
    return read_edits({{from, to}}, name);
}

/** SubDyn's KBBt as tower<case>_kbbt.txt holds it: 6 rows of 6, '#' lines skipped. */
std::array<double,36> subdyn_kbbt (const std::string& c)
{
    std::ifstream f(std::string(ERF_FRAME_TEST_FILES) + "/case" + c + "/tower" + c + "_kbbt.txt");
    std::array<double,36> k{};
    std::string line;
    std::size_t n = 0;
    while (std::getline(f, line) && n < 36) {
        if (line.empty() || line[0] == '#') { continue; }
        std::istringstream is(line);
        double v = 0.0;
        while (n < 36 && (is >> v)) { k[n++] = v; }
    }
    EXPECT_EQ(n, 36u);
    return k;
}

/** The frame's stiffness condensed to joint id: the inverse of the 6 x 6 flexibility from six unit loads. */
std::array<double,36> condensed_stiffness (const Frame& f, int id)
{
    const std::size_t node = f.node_of_joint(id);
    std::vector<double> flex(36, 0.0);
    for (std::size_t k = 0; k < 6; ++k) {
        std::vector<double> load(f.num_dofs(), 0.0);
        load[6 * node + k] = 1.0;
        const FrameSolution s = f.solve(load, 0.0);
        for (std::size_t i = 0; i < 6; ++i) { flex[6 * i + k] = s.displacement[6 * node + i]; }
    }
    DenseCholesky c;
    EXPECT_LT(c.factor(flex, 6), 0) << "the flexibility at joint " << id << " is not positive definite";
    std::array<double,36> k{};
    for (std::size_t j = 0; j < 6; ++j) {
        std::vector<double> e(6, 0.0);
        e[j] = 1.0;
        c.solve(e);
        for (std::size_t i = 0; i < 6; ++i) { k[6 * i + j] = e[i]; }
    }
    return k;
}

} // namespace

TEST(FrameSubDyn, ReadsEveryTableOfASubDynFile)
{
    FrameInputs in;
    const std::string err = read_subdyn(case_file("B"), in);
    ASSERT_TRUE(err.empty()) << err;
    EXPECT_TRUE(in.validate().empty()) << in.validate();
    EXPECT_EQ(in.theory, BeamTheory::Timoshenko);
    EXPECT_EQ(in.divisions, 2);
    EXPECT_EQ(in.joints.size(), 19u);
    EXPECT_EQ(in.members.size(), 56u);
    EXPECT_EQ(in.supports.size(), 4u);
    ASSERT_EQ(in.interface_joints.size(), 1u);
    EXPECT_EQ(in.interface_joints[0], 17);
    EXPECT_DOUBLE_EQ(in.joints[16].x[2], 25.0);
    // joint 1 stands free on its spring, the others are fixed
    const FrameSupport& s1 = in.supports[0];
    EXPECT_EQ(s1.joint, 1);
    for (const bool b : s1.fixed) { EXPECT_FALSE(b); }
    EXPECT_DOUBLE_EQ(s1.stiffness[0], 2.0e8);    // Kxx
    EXPECT_DOUBLE_EQ(s1.stiffness[5], 5.0e8);    // Kzz
    EXPECT_DOUBLE_EQ(s1.stiffness[10], 1.0e7);   // Kxty
    EXPECT_DOUBLE_EQ(s1.stiffness[20], 1.0e8);   // Ktztz
    for (const bool b : in.supports[1].fixed) { EXPECT_TRUE(b); }
    // legs circular, struts rectangular, diagonals arbitrary and spun 30 deg
    EXPECT_EQ(in.members[0].shape, SectionShape::Circular);
    EXPECT_EQ(in.members[12].shape, SectionShape::Rectangular);
    EXPECT_EQ(in.members[24].shape, SectionShape::Arbitrary);
    EXPECT_NEAR(in.members[24].spin, 30.0 * 3.14159265358979323846 / 180.0, 1e-15);
    const FrameSection* leg = in.section(in.members[0].section, SectionShape::Circular);
    ASSERT_NE(leg, nullptr);
    EXPECT_DOUBLE_EQ(leg->D, 0.3);
    EXPECT_DOUBLE_EQ(leg->t, 0.012);
    std::string ferr;
    auto f = Frame::create(in, ferr);
    ASSERT_TRUE(f) << ferr;
    EXPECT_EQ(f->num_nodes(), 19u + 56u);        // a mid-node per member
    EXPECT_EQ(f->elements().size(), 112u);
    EXPECT_EQ(f->num_free_dofs(), 6u * 75u - 18u);
}

TEST(FrameSubDyn, RefusesWhatTheFrameCannotModelNamingIt)
{
    struct Edit { std::string from, to, expect; };
    const std::string row = "    1          1           5            1             1        4         0";
    const std::vector<Edit> edits = {
        {"             1   FEMMod", "             2   FEMMod", "FEMMod = 2"},
        {"             1   NDiv", "             0   NDiv", "NDiv must be >= 1"},
        {row, "    1          1           5            1             1        2         0", "cables"},
        {row, "    1          1           5            1             1        3         0", "rigid links"},
        {row, "    1          1           5            1             2        4         0", "different property sets"},
        {row, "    1          1          99            1             1        4         0", "joint 99"},
        {row, "    1          1           5            7             7        4         0", "property set 7"},
        {"   1          1           1           1           1           1           1",
         "   1          0           0           0           0           0           0        \"no_such_ssi.dat\"", "cannot read the SSI file"},
        {"             2   NXPropSets", "             3   NXPropSets", "ARBITRARY BEAM CROSS-SECTION PROPERTIES"},
        {"0.000000000000000E+00   1   0.0   0.0   0.0   0.0", "0.000000000000000E+00   2   0.0   0.0   0.0   0.0", "JointType 2"},
    };
    for (std::size_t i = 0; i < edits.size(); ++i) {
        const std::string err = read_edited(edits[i].from, edits[i].to, "edit" + std::to_string(i));
        EXPECT_NE(err.find(edits[i].expect), std::string::npos) << "edit " << i << ": " << err;
    }
    // a negative Young's modulus in the first arbitrary property set
    {
        FrameInputs in;
        ASSERT_TRUE(read_subdyn(case_file("A"), in).empty());
        in.sections[0].E = -1.0;
        EXPECT_NE(in.validate().find("YoungE > 0"), std::string::npos) << in.validate();
    }
    FrameInputs in;
    EXPECT_NE(read_subdyn("no_such_subdyn_file.dat", in).find("cannot read"), std::string::npos);
}

TEST(FrameSubDyn, RefusesSpringsMassesAndSupportsThatDoNotReadWhole)
{
    // joint 1 free vertically on a spring file beside the frame's
    const std::string fixed = "   1          1           1           1           1           1           1";
    const std::string on_spring = "   1          1           1           0           1           1           1        \"ssi.dat\"";
    auto ssi = [&] (const std::string& contents, const std::string& name) {
        return read_edits({{fixed, on_spring}}, name, {{"ssi.dat", contents}});
    };
    EXPECT_TRUE(ssi("! a vertical spring\n   2.0E+08   Kzz\n   5.0E+02   Mzz\n", "ssi_ok").empty());
    // a misspelt name once left its entry at 0 without a word
    EXPECT_NE(ssi("   2.0E+08   Kzzz\n", "ssi_name").find("'Kzzz' is not one of"), std::string::npos);
    EXPECT_NE(ssi("! nothing but comments\n", "ssi_empty").find("holds no stiffness or mass entry"), std::string::npos);
    EXPECT_NE(ssi("   2.0E+08   Kzz\n   nan   Mzz\n", "ssi_nan").find("non-finite SSI mass"), std::string::npos);
    // a support free in every DOF without a spring holds nothing
    EXPECT_NE(read_edited(fixed, "   1          0           0           0           0           0           0", "free")
                  .find("restrains nothing"), std::string::npos);
    // a concentrated-mass row of 8 values once dropped its products of inertia
    const std::pair<std::string, std::string> one{"             0   NCmass", "             1   NCmass"};
    const std::string next = "\n---------------------------- OUTPUT: SUMMARY";
    EXPECT_NE(read_edits({one, {next, "\n   5   100.0   1.0   1.0   1.0   0.2   0.0   0.0" + next}}, "cmass8")
                  .find("has 8 values"), std::string::npos);
    EXPECT_TRUE(read_edits({one, {next, "\n   5   100.0   1.0   1.0   1.0" + next}}, "cmass5").empty());
    EXPECT_NE(read_edits({one, {next, "\n   5   100.0   nan   1.0   1.0   0.0   0.0   0.0   0.0   0.0   0.0" + next}}, "cmass_nan")
                  .find("non-finite inertia"), std::string::npos);
}

TEST(FrameSubDyn, TheStiffnessAtThePeakIsSubDyns)
{
    for (const std::string c : {"A", "B", "T"}) {
        FrameInputs in;
        const std::string err = read_subdyn(case_file(c), in);
        ASSERT_TRUE(err.empty()) << err;
        std::string ferr;
        auto f = Frame::create(in, ferr);
        ASSERT_TRUE(f) << ferr;
        // the interface joint: the peak of cases A and B, the cross-arm's centre of case T
        const auto ours = condensed_stiffness(*f, c == "T" ? 21 : 17);
        const auto theirs = subdyn_kbbt(c);
        // each entry to 2e-6 of the geometric mean of its row's and column's diagonal: SubDyn prints 7 digits
        for (std::size_t i = 0; i < 6; ++i) {
            for (std::size_t j = 0; j < 6; ++j) {
                const double scale = std::sqrt(theirs[7 * i] * theirs[7 * j]);
                EXPECT_NEAR(ours[6 * i + j], theirs[6 * i + j], 2.0e-6 * scale) << "case " << c << ", KBBt(" << i + 1 << "," << j + 1 << ")";
            }
        }
    }
}
