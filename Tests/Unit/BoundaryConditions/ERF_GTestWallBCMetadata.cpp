#include <ERF_WallScalarBC.H>

#include <gtest/gtest.h>

using amrex::Real;

namespace {

// Regression motivation:
// Before Stage 0, solid-wall representation was selected by one last-parsed
// density flag and scalar presence was inferred from sentinel values. These
// tests exercise the production resolver's independent oracle: explicit
// ParmParse presence plus same-face conserved/primitive conversion. A failure
// means a wall can again inherit another face's representation or lose an
// explicitly supplied zero.
TEST(WallBCMetadata, SameFaceDensityControlsRepresentation)
{
    const auto conserved = erf_wall_scalar_bc::resolve_wall_scalar_value(true, Real(300.0), true, Real(1.0));
    EXPECT_EQ(conserved.intent, erf_wall_scalar_bc::WallScalarBCIntent::DirichletConserved);
    EXPECT_EQ(conserved.stored_value, Real(300.0));

    const auto primitive = erf_wall_scalar_bc::resolve_wall_scalar_value(true, Real(301.0), false, Real(0.0));
    EXPECT_EQ(primitive.intent, erf_wall_scalar_bc::WallScalarBCIntent::DirichletPrimitive);
    EXPECT_EQ(primitive.stored_value, Real(301.0));
}

TEST(WallBCMetadata, ParseOrderCannotChangeMixedFaces)
{
    const auto early_density = erf_wall_scalar_bc::resolve_wall_scalar_value(true, Real(300.0), true, Real(1.0));
    const auto later_primitive = erf_wall_scalar_bc::resolve_wall_scalar_value(true, Real(301.0), false, Real(0.0));
    const auto later_density = erf_wall_scalar_bc::resolve_wall_scalar_value(true, Real(302.0), true, Real(0.9));

    EXPECT_EQ(early_density.intent, erf_wall_scalar_bc::WallScalarBCIntent::DirichletConserved);
    EXPECT_EQ(later_primitive.intent, erf_wall_scalar_bc::WallScalarBCIntent::DirichletPrimitive);
    EXPECT_EQ(later_primitive.stored_value, Real(301.0));
    EXPECT_EQ(later_density.intent, erf_wall_scalar_bc::WallScalarBCIntent::DirichletConserved);
    EXPECT_EQ(later_density.stored_value, Real(271.8));
}

TEST(WallBCMetadata, HistoricalSentinelsAndExplicitZerosArePresent)
{
    const auto unit_density = erf_wall_scalar_bc::resolve_wall_scalar_value(true, Real(300.0), true, Real(1.0));
    const auto dry_primitive = erf_wall_scalar_bc::resolve_wall_scalar_value(true, Real(0.0), false, Real(0.0));
    const auto dry_conserved = erf_wall_scalar_bc::resolve_wall_scalar_value(true, Real(0.0), true, Real(1.2));
    const auto zero_gradient = erf_wall_scalar_bc::resolve_wall_scalar_gradient(true, Real(0.0));

    EXPECT_EQ(unit_density.intent, erf_wall_scalar_bc::WallScalarBCIntent::DirichletConserved);
    EXPECT_EQ(dry_primitive.intent, erf_wall_scalar_bc::WallScalarBCIntent::DirichletPrimitive);
    EXPECT_EQ(dry_primitive.stored_value, Real(0.0));
    EXPECT_EQ(dry_conserved.intent, erf_wall_scalar_bc::WallScalarBCIntent::DirichletConserved);
    EXPECT_EQ(dry_conserved.stored_value, Real(0.0));
    EXPECT_EQ(zero_gradient.intent, erf_wall_scalar_bc::WallScalarBCIntent::Neumann);
    EXPECT_EQ(zero_gradient.stored_value, Real(0.0));
}

TEST(WallBCMetadata, ValueAndGradientConflictNamesBothKeywords)
{
    EXPECT_TRUE(erf_wall_scalar_bc::wall_scalar_value_gradient_conflict(true, true));
    const std::string message = erf_wall_scalar_bc::wall_scalar_value_gradient_error("zhi", "theta");
    EXPECT_NE(message.find("zhi.theta"), std::string::npos);
    EXPECT_NE(message.find("zhi.theta_grad"), std::string::npos);
    EXPECT_FALSE(erf_wall_scalar_bc::wall_scalar_value_gradient_conflict(true, false));
}

TEST(WallBCMetadata, AnelasticDensityIsRejectedWithFaceAndKeyword)
{
    const std::string value_message = erf_wall_scalar_bc::anelastic_wall_density_error("xlo", true, false);
    const std::string gradient_message = erf_wall_scalar_bc::anelastic_wall_density_error("zhi", false, true);
    EXPECT_NE(value_message.find("xlo.density"), std::string::npos);
    EXPECT_NE(gradient_message.find("zhi.density_grad"), std::string::npos);
    EXPECT_TRUE(erf_wall_scalar_bc::anelastic_wall_density_error("xhi", false, false).empty());
}

TEST(WallBCMetadata, UnspecifiedValuesRemainUnspecified)
{
    const auto value = erf_wall_scalar_bc::resolve_wall_scalar_value(false, Real(999.0), false, Real(999.0));
    const auto gradient = erf_wall_scalar_bc::resolve_wall_scalar_gradient(false, Real(999.0));
    EXPECT_EQ(value.intent, erf_wall_scalar_bc::WallScalarBCIntent::Unspecified);
    EXPECT_EQ(gradient.intent, erf_wall_scalar_bc::WallScalarBCIntent::Unspecified);
}

} // namespace
