/**
 * @file ERF_GTestStructureModelExample.cpp
 * @brief A unit test that a structure model directory builds into ERF's unit-test executable
 *
 * ERFStructureModel.cmake adds this file through erf_structure_model_add_unit_tests, so it
 * is compiled into erf_unit_tests with ERF's test main, and the example model it tests is
 * already in that executable through erf_structure_model_add. The test creates the model
 * by its registered name and checks the box support() returns: the point with the reach
 * of its kernel, 3 epsilon, on each side.
 */
#include <gtest/gtest.h>

#include <ERF_StructureModelRegistry.H>
#include <AMReX_ParmParse.H>

#include <string>
#include <vector>

using namespace amrex;

namespace {

/// Sets the structure.example inputs for one test and removes them afterwards
class ScopedExampleInputs
{
public:
    ScopedExampleInputs (const std::vector<Real>& pos, const std::vector<Real>& force, Real eps)
        : m_pp("structure.example")
    {
        m_pp.addarr("position", pos);
        m_pp.addarr("force", force);
        m_pp.add("epsilon", eps);
    }
    ~ScopedExampleInputs ()
    {
        m_pp.remove("position");
        m_pp.remove("force");
        m_pp.remove("epsilon");
    }
    ScopedExampleInputs (const ScopedExampleInputs&) = delete;
    ScopedExampleInputs& operator= (const ScopedExampleInputs&) = delete;
private:
    ParmParse m_pp;
};

} // namespace

TEST(StructureModelExample, SupportIsThePointAndItsKernelReach)
{
    const std::vector<Real> pos{Real(200.0), Real(150.0), Real(100.0)};
    const Real eps = Real(10.0);
    ScopedExampleInputs inputs(pos, {Real(50.0), Real(0.0), Real(0.0)}, eps);

    auto model = make_structure_model("example", 0);
    ASSERT_NE(model, nullptr) << "the example model is not compiled into this executable";
    EXPECT_STREQ(model->name(), "example");
    EXPECT_TRUE(model->supports_restart());

    // Exact in single and double precision: every value is an integer below 2^24
    const std::vector<RealBox> boxes = model->support();
    ASSERT_EQ(boxes.size(), 1u);
    for (int d = 0; d < 3; ++d) {
        EXPECT_EQ(boxes[0].lo(d), pos[d] - Real(3.0) * eps) << "direction " << d;
        EXPECT_EQ(boxes[0].hi(d), pos[d] + Real(3.0) * eps) << "direction " << d;
    }
}
