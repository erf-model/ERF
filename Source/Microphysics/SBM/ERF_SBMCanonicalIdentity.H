#ifndef ERF_SBM_CANONICAL_IDENTITY_H_
#define ERF_SBM_CANONICAL_IDENTITY_H_

#include <AMReX_REAL.H>

#include <cmath>
#include <iomanip>
#include <limits>
#include <locale>
#include <sstream>
#include <string>

namespace erf_sbm {

// Scientific metadata is part of the restart/schema identity.  Keep this
// policy in one place so grids and layouts cannot silently serialize the same
// value with different precision or locale rules.
inline std::string canonical_real(const amrex::Real value,
                                  const bool allow_unbounded_nan = false)
{
    if (allow_unbounded_nan && std::isnan(value)) return "unbounded";
    std::ostringstream out;
    out.imbue(std::locale::classic());
    out.setf(std::ios_base::scientific, std::ios_base::floatfield);
    out << std::setprecision(std::numeric_limits<amrex::Real>::max_digits10) << value;
    return out.str();
}

} // namespace erf_sbm

#endif
