#include <AMReX_ParmParse.H>

#include <ERF_WallScalarBC.H>

#include <gtest/gtest.h>

#include <string>
#include <vector>

using amrex::ParmParse;
using amrex::Real;

namespace {

class ScopedWallParams
{
public:
    explicit ScopedWallParams (const char* prefix)
        : m_pp(prefix) {}

    ~ScopedWallParams ()
    {
        for (const auto& name : m_names) {
            m_pp.remove(name);
        }
    }

    void add (const char* name, Real value)
    {
        m_pp.add(name, value);
        m_names.emplace_back(name);
    }

    ParmParse& pp () { return m_pp; }

private:
    ParmParse m_pp;
    std::vector<std::string> m_names;
};

using erf_wall_scalar_bc::SolidWallKind;
using erf_wall_scalar_bc::WallScalarBCIntent;

// Motivation: Each wall face must parse its own density policy, so wall order
// cannot change whether theta is conserved or primitive.
TEST(WallBCMetadata, ParmParseSameFaceDensityControlsRepresentation)
{
    ScopedWallParams xlo("wall_parser_order_a_xlo");
    xlo.add("density", Real(1.0));
    xlo.add("theta", Real(300.0));
    ScopedWallParams xhi("wall_parser_order_a_xhi");
    xhi.add("theta", Real(301.0));

    const auto first = erf_wall_scalar_bc::parse_wall_face_scalars(
        xlo.pp(), "xlo", SolidWallKind::NoSlip, false);
    const auto second = erf_wall_scalar_bc::parse_wall_face_scalars(
        xhi.pp(), "xhi", SolidWallKind::NoSlip, false);
    ASSERT_TRUE(first.ok()) << first.error;
    ASSERT_TRUE(second.ok()) << second.error;
    EXPECT_EQ(first.scalars.theta.intent, WallScalarBCIntent::DirichletConserved);
    EXPECT_EQ(first.scalars.theta.stored_value, Real(300.0));
    EXPECT_EQ(second.scalars.theta.intent, WallScalarBCIntent::DirichletPrimitive);
    EXPECT_EQ(second.scalars.theta.stored_value, Real(301.0));

    ScopedWallParams xlo_reverse("wall_parser_order_b_xlo");
    xlo_reverse.add("density", Real(0.9));
    xlo_reverse.add("theta", Real(302.0));
    ScopedWallParams xhi_reverse("wall_parser_order_b_xhi");
    xhi_reverse.add("theta", Real(303.0));

    const auto reverse_second = erf_wall_scalar_bc::parse_wall_face_scalars(
        xhi_reverse.pp(), "xhi", SolidWallKind::NoSlip, false);
    const auto reverse_first = erf_wall_scalar_bc::parse_wall_face_scalars(
        xlo_reverse.pp(), "xlo", SolidWallKind::NoSlip, false);
    ASSERT_TRUE(reverse_first.ok()) << reverse_first.error;
    ASSERT_TRUE(reverse_second.ok()) << reverse_second.error;
    EXPECT_EQ(reverse_first.scalars.theta.intent, WallScalarBCIntent::DirichletConserved);
    EXPECT_EQ(reverse_first.scalars.theta.stored_value, Real(271.8));
    EXPECT_EQ(reverse_second.scalars.theta.intent, WallScalarBCIntent::DirichletPrimitive);
    EXPECT_EQ(reverse_second.scalars.theta.stored_value, Real(303.0));
}

// Motivation: Explicit zero values are valid inputs; ParmParse presence, not
// numeric sentinels, must preserve density, qv, and theta gradients.
TEST(WallBCMetadata, ParmParseExplicitZerosArePresent)
{
    ScopedWallParams params("wall_parser_zero_values");
    params.add("density", Real(1.0));
    params.add("qv", Real(0.0));
    params.add("theta_grad", Real(0.0));

    const auto parsed = erf_wall_scalar_bc::parse_wall_face_scalars(
        params.pp(), "xlo", SolidWallKind::NoSlip, false);
    ASSERT_TRUE(parsed.ok()) << parsed.error;
    EXPECT_EQ(parsed.scalars.density.intent, WallScalarBCIntent::DirichletConserved);
    EXPECT_EQ(parsed.scalars.density.stored_value, Real(1.0));
    EXPECT_EQ(parsed.scalars.qv.intent, WallScalarBCIntent::DirichletConserved);
    EXPECT_EQ(parsed.scalars.qv.stored_value, Real(0.0));
    EXPECT_EQ(parsed.scalars.theta.intent, WallScalarBCIntent::Neumann);
    EXPECT_EQ(parsed.scalars.theta.stored_value, Real(0.0));
}

// Motivation: An explicitly dry no-slip wall without same-face density must
// use primitive storage.
TEST(WallBCMetadata, ParmParseDryNoSlipWallWithoutDensityIsPrimitive)
{
    ScopedWallParams params("wall_parser_dry_primitive");
    params.add("qv", Real(0.0));

    const auto parsed = erf_wall_scalar_bc::parse_wall_face_scalars(
        params.pp(), "xlo", SolidWallKind::NoSlip, false);
    ASSERT_TRUE(parsed.ok()) << parsed.error;
    EXPECT_EQ(parsed.scalars.qv.intent, WallScalarBCIntent::DirichletPrimitive);
    EXPECT_EQ(parsed.scalars.qv.stored_value, Real(0.0));
    EXPECT_EQ(parsed.scalars.density.intent, WallScalarBCIntent::Unspecified);
}

// Motivation: Anelastic solid walls must reject density and density_grad on
// either wall kind with the face-specific diagnostic.
TEST(WallBCMetadata, ParmParseAnelasticDensityIsRejected)
{
    for (const auto kind : {SolidWallKind::NoSlip, SolidWallKind::Slip}) {
        ScopedWallParams density("wall_parser_anelastic_density");
        density.add("density", Real(1.0));
        const auto density_result = erf_wall_scalar_bc::parse_wall_face_scalars(
            density.pp(), "xlo", kind, true);
        EXPECT_FALSE(density_result.ok());
        EXPECT_NE(density_result.error.find("xlo.density"), std::string::npos);

        ScopedWallParams gradient("wall_parser_anelastic_gradient");
        gradient.add("density_grad", Real(0.0));
        const auto gradient_result = erf_wall_scalar_bc::parse_wall_face_scalars(
            gradient.pp(), "zhi", kind, true);
        EXPECT_FALSE(gradient_result.ok());
        EXPECT_NE(gradient_result.error.find("zhi.density_grad"), std::string::npos);
    }
}

// Motivation: NoSlipWall does not implement density_grad; fixed density
// remains supported for compressible walls.
TEST(WallBCMetadata, ParmParseNoSlipDensityGradientIsUnsupported)
{
    for (const Real gradient_value : {Real(0.0), Real(0.25)}) {
        ScopedWallParams params("wall_parser_noslip_gradient");
        params.add("density_grad", gradient_value);
        const auto parsed = erf_wall_scalar_bc::parse_wall_face_scalars(
            params.pp(), "xlo", SolidWallKind::NoSlip, false);
        EXPECT_FALSE(parsed.ok());
        EXPECT_NE(parsed.error.find("xlo.density_grad"), std::string::npos);
        EXPECT_NE(parsed.error.find("NoSlipWall"), std::string::npos);
    }

    ScopedWallParams fixed_density("wall_parser_noslip_density");
    fixed_density.add("density", Real(1.0));
    const auto parsed = erf_wall_scalar_bc::parse_wall_face_scalars(
        fixed_density.pp(), "xhi", SolidWallKind::NoSlip, false);
    ASSERT_TRUE(parsed.ok()) << parsed.error;
    EXPECT_EQ(parsed.scalars.density.intent, WallScalarBCIntent::DirichletConserved);
}

// Motivation: SlipWall density gradients remain supported, including zero,
// while value/gradient conflicts stay rejected and qv remains unsupported.
TEST(WallBCMetadata, ParmParseSlipDensityGradientRemainsSupported)
{
    ScopedWallParams zero_gradient("wall_parser_slip_gradient");
    zero_gradient.add("density_grad", Real(0.0));
    const auto parsed = erf_wall_scalar_bc::parse_wall_face_scalars(
        zero_gradient.pp(), "zhi", SolidWallKind::Slip, false);
    ASSERT_TRUE(parsed.ok()) << parsed.error;
    EXPECT_EQ(parsed.scalars.density.intent, WallScalarBCIntent::Neumann);
    EXPECT_EQ(parsed.scalars.density.stored_value, Real(0.0));

    ScopedWallParams conflict("wall_parser_slip_conflict");
    conflict.add("density", Real(1.0));
    conflict.add("density_grad", Real(0.1));
    const auto conflict_result = erf_wall_scalar_bc::parse_wall_face_scalars(
        conflict.pp(), "zhi", SolidWallKind::Slip, false);
    EXPECT_FALSE(conflict_result.ok());
    EXPECT_NE(conflict_result.error.find("zhi.density"), std::string::npos);
}

// Motivation: theta and theta_grad are mutually exclusive on both wall kinds,
// including explicit zero gradients.
TEST(WallBCMetadata, ParmParseThetaValueGradientConflict)
{
    for (const auto kind : {SolidWallKind::NoSlip, SolidWallKind::Slip}) {
        ScopedWallParams params("wall_parser_theta_conflict");
        params.add("theta", Real(300.0));
        params.add("theta_grad", Real(0.0));
        const auto parsed = erf_wall_scalar_bc::parse_wall_face_scalars(
            params.pp(), "yhi", kind, false);
        EXPECT_FALSE(parsed.ok());
        EXPECT_NE(parsed.error.find("yhi.theta"), std::string::npos);
        EXPECT_NE(parsed.error.find("yhi.theta_grad"), std::string::npos);
    }
}

// Motivation: Omitted wall keywords must remain Unspecified rather than
// inheriting a temporary or sentinel value.
TEST(WallBCMetadata, ParmParseAbsenceRemainsUnspecified)
{
    ScopedWallParams params("wall_parser_absence");
    const auto parsed = erf_wall_scalar_bc::parse_wall_face_scalars(
        params.pp(), "xhi", SolidWallKind::NoSlip, false);
    ASSERT_TRUE(parsed.ok()) << parsed.error;
    EXPECT_EQ(parsed.scalars.density.intent, WallScalarBCIntent::Unspecified);
    EXPECT_EQ(parsed.scalars.theta.intent, WallScalarBCIntent::Unspecified);
    EXPECT_EQ(parsed.scalars.qv.intent, WallScalarBCIntent::Unspecified);
}

// Motivation: qv is supported only on NoSlipWall and must not alter SlipWall
// metadata.
TEST(WallBCMetadata, ParmParseSlipWallDoesNotApplyQv)
{
    ScopedWallParams params("wall_parser_slip_qv");
    params.add("qv", Real(0.0));
    const auto parsed = erf_wall_scalar_bc::parse_wall_face_scalars(
        params.pp(), "xlo", SolidWallKind::Slip, false);
    ASSERT_TRUE(parsed.ok()) << parsed.error;
    EXPECT_EQ(parsed.scalars.qv.intent, WallScalarBCIntent::Unspecified);
}

// Motivation: The pure resolver must preserve explicit zero values for
// primitive and gradient inputs.
TEST(WallBCMetadata, PureResolverPreservesExplicitZero)
{
    const auto primitive = erf_wall_scalar_bc::resolve_wall_scalar_value(
        true, Real(0.0), false, Real(0.0));
    const auto gradient = erf_wall_scalar_bc::resolve_wall_scalar_gradient(
        true, Real(0.0));
    EXPECT_EQ(primitive.intent, WallScalarBCIntent::DirichletPrimitive);
    EXPECT_EQ(primitive.stored_value, Real(0.0));
    EXPECT_EQ(gradient.intent, WallScalarBCIntent::Neumann);
    EXPECT_EQ(gradient.stored_value, Real(0.0));
}

} // namespace
