#include "ERF_SBMPhysicalBoundary.H"
#include "ERF_SBMCanonicalIdentity.H"
#include "ERF_SBMConstraintGroups.H"
#include "ERF_SBMFixtureValidation.H"

#include <cmath>
#include <locale>
#include <sstream>
#include <utility>

namespace erf_sbm {

namespace {

constexpr std::array<int, physical_face_count> identity_order{0, 3, 1, 4, 2, 5};

} // namespace

bool
SBMPhysicalBoundaryConfig::all_periodic () const noexcept
{
    for (const auto face_kind : kind) {
        if (face_kind != AdvectionFaceKind::Periodic) return false;
    }
    return true;
}

bool
SBMPhysicalBoundaryConfig::periodic_in (const int dir) const noexcept
{
    return kind[static_cast<std::size_t>(dir)] == AdvectionFaceKind::Periodic &&
           kind[static_cast<std::size_t>(dir + AMREX_SPACEDIM)] ==
               AdvectionFaceKind::Periodic;
}

std::string
SBMPhysicalBoundaryConfig::identity () const
{
    if (all_periodic()) return "boundary=triply-periodic-v1";
    std::ostringstream line;
    line << "boundary=physical-faces-v1";
    for (const int face : identity_order) {
        const auto face_kind = kind[static_cast<std::size_t>(face)];
        line << ';' << face_name(face) << '=' << face_kind_name(face_kind);
        if (face_kind == AdvectionFaceKind::PrescribedInflow) {
            line << '[';
            const auto& values = inflow_z[static_cast<std::size_t>(face)];
            for (std::size_t n = 0; n < values.size(); ++n) {
                line << (n == 0 ? "" : ",") << canonical_real(values[n]);
            }
            line << ']';
        }
    }
    return line.str();
}

const char*
face_name (const int face) noexcept
{
    static constexpr const char* names[physical_face_count] =
        {"xlo", "ylo", "zlo", "xhi", "yhi", "zhi"};
    return face >= 0 && face < physical_face_count ? names[face] : "unknown";
}

const char*
face_kind_name (const AdvectionFaceKind kind) noexcept
{
    switch (kind) {
    case AdvectionFaceKind::Periodic: return "periodic";
    case AdvectionFaceKind::ImpermeableWall: return "wall";
    case AdvectionFaceKind::AdvectiveOutflow: return "outflow";
    case AdvectionFaceKind::PrescribedInflow: return "inflow";
    }
    return "unknown";
}

bool
map_erf_boundary (const ERF_BC bc, AdvectionFaceKind& kind, std::string& reason)
{
    reason.clear();
    switch (bc) {
    case ERF_BC::periodic:
        kind = AdvectionFaceKind::Periodic;
        return true;
    case ERF_BC::slip_wall:
    case ERF_BC::no_slip_wall:
        kind = AdvectionFaceKind::ImpermeableWall;
        return true;
    case ERF_BC::surface_layer:
        reason = "a surface_layer boundary acts through the diffusion operator, "
                 "and SBM does not yet support diffusion";
        return false;
    case ERF_BC::symmetry:
        reason = "ERF does not hold the normal dry-air carrier of a symmetry face "
                 "at zero, so it is not an impermeable spectral wall";
        return false;
    case ERF_BC::outflow:
    case ERF_BC::ho_outflow:
        kind = AdvectionFaceKind::AdvectiveOutflow;
        return true;
    case ERF_BC::inflow:
    case ERF_BC::inflow_outflow:
        kind = AdvectionFaceKind::PrescribedInflow;
        return true;
    case ERF_BC::open:
        reason = "open boundaries advect with a characteristic scheme that the "
                 "spectral transport does not share";
        return false;
    case ERF_BC::undefined:
        break;
    }
    reason = "the boundary type is undefined";
    return false;
}

bool
validate_inflow_state (const SBMLayout& layout,
                       const std::vector<amrex::Real>& z,
                       std::string& reason)
{
    reason.clear();
    if (z.size() != static_cast<std::size_t>(layout.ncomp())) {
        reason = "it has " + std::to_string(z.size()) + " values but the SBM layout has " +
                 std::to_string(layout.ncomp()) + " components";
        return false;
    }
    for (std::size_t n = 0; n < z.size(); ++n) {
        if (!std::isfinite(z[n]) || z[n] < amrex::Real(0.0)) {
            reason = "component " + std::to_string(n) + " is not finite and nonnegative";
            return false;
        }
    }
    const auto canonicality = validate_fixture_initial_state_canonicality(layout, z);
    if (!canonicality.canonical) {
        reason = "population " + std::to_string(canonicality.population_id) + " bin " +
                 std::to_string(canonicality.bin) +
                 " is not canonical for its spectral representation";
        return false;
    }
    for (const auto& group : make_constraint_groups(layout)) {
        std::string failed_constraint;
        if (!group.admissible(z, nullptr, &failed_constraint)) {
            reason = "population " + std::to_string(group.population_id) + " bin " +
                     std::to_string(group.bin) + " (" + group.semantic_id +
                     ") fails " + failed_constraint;
            return false;
        }
    }
    return true;
}

void
SBMBoundaryLedger::reset (const int components)
{
    ncomp = components;
    outward.assign(static_cast<std::size_t>(physical_face_count * components), 0.0);
    gross.assign(static_cast<std::size_t>(physical_face_count * components), 0.0);
    face_stage_terms = 0;
}

double
SBMBoundaryLedger::total_outward (const int component) const
{
    double total = 0.0;
    for (int face = 0; face < physical_face_count; ++face) {
        total += outward[static_cast<std::size_t>(face * ncomp + component)];
    }
    return total;
}

double
SBMBoundaryLedger::total_gross (const int component) const
{
    double total = 0.0;
    for (int face = 0; face < physical_face_count; ++face) {
        total += gross[static_cast<std::size_t>(face * ncomp + component)];
    }
    return total;
}

std::string
SBMBoundaryLedger::serialize () const
{
    std::ostringstream text;
    text.imbue(std::locale::classic());
    text << "ERF-SBM-BOUNDARY-LEDGER-v1\n"
         << "ncomp=" << ncomp << '\n'
         << "face_stage_terms=" << face_stage_terms << '\n';
    for (int face = 0; face < physical_face_count; ++face) {
        for (int component = 0; component < ncomp; ++component) {
            const auto index = static_cast<std::size_t>(face * ncomp + component);
            text << face_name(face) << ' ' << component << ' '
                 << canonical_decimal(outward[index]) << ' '
                 << canonical_decimal(gross[index]) << '\n';
        }
    }
    return text.str();
}

bool
SBMBoundaryLedger::parse (const std::string& text, const int components,
                          SBMBoundaryLedger& ledger, std::string& reason)
{
    reason.clear();
    std::istringstream input(text);
    input.imbue(std::locale::classic());
    std::string header;
    std::string ncomp_line;
    std::string terms_line;
    if (!std::getline(input, header) || header != "ERF-SBM-BOUNDARY-LEDGER-v1" ||
        !std::getline(input, ncomp_line) ||
        ncomp_line != "ncomp=" + std::to_string(components) ||
        !std::getline(input, terms_line) ||
        terms_line.rfind("face_stage_terms=", 0) != 0) {
        reason = "the boundary ledger header does not match the active layout";
        return false;
    }
    SBMBoundaryLedger parsed;
    parsed.reset(components);
    std::istringstream terms(terms_line.substr(std::string("face_stage_terms=").size()));
    terms.imbue(std::locale::classic());
    if (!(terms >> parsed.face_stage_terms)) {
        reason = "the boundary ledger operation count is unreadable";
        return false;
    }
    for (int face = 0; face < physical_face_count; ++face) {
        for (int component = 0; component < components; ++component) {
            std::string name;
            int index = -1;
            double out_value = 0.0;
            double gross_value = 0.0;
            if (!(input >> name >> index >> out_value >> gross_value) ||
                name != face_name(face) || index != component ||
                !std::isfinite(out_value) || !std::isfinite(gross_value) ||
                gross_value < 0.0) {
                reason = "the boundary ledger entry for " + std::string(face_name(face)) +
                         " component " + std::to_string(component) + " is invalid";
                return false;
            }
            const auto slot = static_cast<std::size_t>(face * components + component);
            parsed.outward[slot] = out_value;
            parsed.gross[slot] = gross_value;
        }
    }
    std::string trailing;
    if (input >> trailing) {
        reason = "the boundary ledger has trailing data";
        return false;
    }
    ledger = std::move(parsed);
    return true;
}

} // namespace erf_sbm
