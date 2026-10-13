#ifndef ERF_SURFACE_TEMPERATURE_H_
#define ERF_SURFACE_TEMPERATURE_H_

#include <iomanip>
#include <sstream>

#include <AMReX_Math.H>
#include <AMReX_GpuAtomic.H>
#include <AMReX_ParallelDescriptor.H>
#include <AMReX_Print.H>

#include <ERF_Constants.H>
#include <ERF_EOS.H>

namespace erf_surface_temperature {

/**
 * Pressure at a boundary face from the adjacent atmospheric cell state.
 *
 * The state pressure is diagnosed at the cell centre. Adding rho*g*delta_z
 * gives the pressure at the boundary face, where delta_z is the signed
 * distance from the face to the cell centre.
 */
AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE
amrex::Real
pressure_at_boundary_from_cell (const amrex::Real rho,
                                const amrex::Real rho_theta,
                                const amrex::Real qv,
                                const amrex::Real delta_z)
{
    return getPgivenRTh(rho_theta, qv) + rho * CONST_GRAV * delta_z;
}

AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE
amrex::Real
pressure_at_surface (const amrex::Real rho,
                     const amrex::Real rho_theta,
                     const amrex::Real qv,
                     const amrex::Real delta_z)
{
    return pressure_at_boundary_from_cell(rho, rho_theta, qv, delta_z);
}

AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE
bool
valid_pressure (const amrex::Real pressure)
{
    return amrex::Math::isfinite(pressure) && pressure > amrex::Real(0.0);
}

AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE
bool
valid_temperature (const amrex::Real temperature)
{
    return amrex::Math::isfinite(temperature) && temperature > amrex::Real(0.0);
}

AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE
bool
theta_to_temperature (const amrex::Real theta,
                      const amrex::Real pressure,
                      const amrex::Real rdOcp,
                      amrex::Real& temperature)
{
    if (!valid_pressure(pressure) || !valid_temperature(theta) ||
        !amrex::Math::isfinite(rdOcp) || rdOcp <= amrex::Real(0.0)) {
        return false;
    }
    const amrex::Real result = getTgivenPandTh(pressure, theta, rdOcp);
    if (!valid_temperature(result)) {
        return false;
    }
    temperature = result;
    return true;
}

AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE
bool
temperature_to_theta (const amrex::Real temperature,
                      const amrex::Real pressure,
                      const amrex::Real rdOcp,
                      amrex::Real& theta)
{
    if (!valid_pressure(pressure) || !valid_temperature(temperature) ||
        !amrex::Math::isfinite(rdOcp) || rdOcp <= amrex::Real(0.0)) {
        return false;
    }
    const amrex::Real result = getThgivenTandP(temperature, pressure, rdOcp);
    if (!valid_temperature(result)) {
        return false;
    }
    theta = result;
    return true;
}

/**
 * Why a surface conversion could not be made, in the order the kernels test for it.
 *
 * Every kernel that converts between the atmospheric state and a surface temperature or
 * humidity guards the same sequence of ways that state can be unusable.  Which guard
 * tripped is the difference between "the surface layer gave up" and knowing that the
 * state handed to it was already bad -- so the guards are named rather than merged into
 * a single flag.
 */
enum struct ConversionError : int {
    none = 0,
    bad_density,
    no_moisture_comp,
    bad_rho_qv,
    bad_qv,
    bad_theta_to_temperature,
    bad_temperature_to_theta,
    bad_pressure,
    bad_skin_temperature,
    bad_q_air,
    bad_q_surf
};

/**
 * The first failing cell on this rank, the state that produced it, and how many cells
 * failed.
 *
 * The count is the part that is easy to leave out and expensive to do without: one bad
 * column and a whole plane of NaN state are different problems with different causes,
 * and nothing else in the report tells them apart.  Cells the fill is not authoritative
 * for are counted separately -- they do not abort the run, but a halo full of failures
 * beside one bad owned cell says where to look.
 */
struct ConversionFailure
{
    int code       = static_cast<int>(ConversionError::none);
    int count      = 0;
    int halo_count = 0;
    int i = 0, j = 0, k = 0;
    amrex::Real rho       = amrex::Real(0.0);
    amrex::Real rho_theta = amrex::Real(0.0);
    amrex::Real rho_qv    = amrex::Real(0.0);
    amrex::Real t_surf    = amrex::Real(0.0);
    amrex::Real pressure  = amrex::Real(0.0);
    amrex::Real delta_z   = amrex::Real(0.0);
    amrex::Real q_air     = amrex::Real(0.0);
};

/**
 * Record one cell's failure.
 *
 * Every failing cell adds to a count; the first authoritative one to arrive also claims
 * the description, so what gets reported is a failure that matters rather than whichever
 * thread happened to finish last.  Pass the values the caller has reached by this point:
 * the report prints the ones that mean something for the code it is given.
 *
 * @param[in,out] fail  This rank's record
 * @param[in] fatal     Whether this is a cell the fill is authoritative for
 * @param[in] code      Which guard tripped
 * @param[in] i,j,k     The cell whose state was read
 */
AMREX_GPU_DEVICE AMREX_FORCE_INLINE
void
record_conversion_failure (ConversionFailure* fail,
                           const bool fatal,
                           const ConversionError code,
                           const int i, const int j, const int k,
                           const amrex::Real rho       = amrex::Real(0.0),
                           const amrex::Real rho_theta = amrex::Real(0.0),
                           const amrex::Real rho_qv    = amrex::Real(0.0),
                           const amrex::Real t_surf    = amrex::Real(0.0),
                           const amrex::Real pressure  = amrex::Real(0.0),
                           const amrex::Real delta_z   = amrex::Real(0.0),
                           const amrex::Real q_air     = amrex::Real(0.0))
{
    if (!fatal) {
        amrex::Gpu::Atomic::Add(&fail->halo_count, 1);
        return;
    }

    amrex::Gpu::Atomic::Add(&fail->count, 1);

    // The first failure to arrive owns the description; the rest only add to the count.
    constexpr int unclaimed = static_cast<int>(ConversionError::none);
    if (amrex::Gpu::Atomic::CAS(&fail->code, unclaimed, static_cast<int>(code)) != unclaimed) {
        return;
    }

    fail->i = i; fail->j = j; fail->k = k;
    fail->rho       = rho;
    fail->rho_theta = rho_theta;
    fail->rho_qv    = rho_qv;
    fail->t_surf    = t_surf;
    fail->pressure  = pressure;
    fail->delta_z   = delta_z;
    fail->q_air     = q_air;
}

AMREX_FORCE_INLINE
const char*
conversion_error_text (const ConversionError code)
{
    switch (code) {
    case ConversionError::bad_density:
        return "the density or (rho theta) of the state cell is not usable";
    case ConversionError::no_moisture_comp:
        return "the state carries no water vapor component to convert";
    case ConversionError::bad_rho_qv:
        return "the (rho qv) of the state cell is not usable";
    case ConversionError::bad_qv:
        return "the water vapor mixing ratio derived from the state cell is not usable";
    case ConversionError::bad_theta_to_temperature:
        return "the surface potential temperature could not be converted to absolute temperature";
    case ConversionError::bad_temperature_to_theta:
        return "the surface temperature could not be converted to potential temperature";
    case ConversionError::bad_pressure:
        return "the surface pressure derived from the state cell is not usable";
    case ConversionError::bad_skin_temperature:
        return "the skin temperature is not usable";
    case ConversionError::bad_q_air:
        return "the reference-height water vapor mixing ratio is not usable";
    case ConversionError::bad_q_surf:
        return "the surface mixing ratio came out non-finite";
    default:
        return "the conversion failed for a reason that was not recorded";
    }
}

/**
 * Describe what this rank recorded, then abort if any rank recorded anything.
 *
 * Collective: every rank must call this.  Each rank that failed describes its own first
 * failure before the abort, so a failure confined to one rank still names the cell and
 * the numbers rather than only the check that tripped.
 *
 * @param[in] who  Name of the fill, used to open the message
 * @param[in] lev  Level being filled
 * @param[in] fail This rank's record
 */
inline void
abort_on_conversion_failure (const char* who,
                             const int lev,
                             ConversionFailure const& fail)
{
    int any_failed = (fail.code != static_cast<int>(ConversionError::none)) ? 1 : 0;
    amrex::ParallelDescriptor::ReduceIntMax(any_failed);
    if (any_failed == 0) { return; }

    if (fail.code != static_cast<int>(ConversionError::none)) {
        const auto code = static_cast<ConversionError>(fail.code);

        const bool is_conversion = (code == ConversionError::bad_theta_to_temperature ||
                                    code == ConversionError::bad_temperature_to_theta);
        const bool has_state = (code != ConversionError::no_moisture_comp);
        const bool has_vapor = (is_conversion ||
                                code == ConversionError::bad_rho_qv ||
                                code == ConversionError::bad_qv     ||
                                code == ConversionError::bad_q_air  ||
                                code == ConversionError::bad_q_surf);
        const bool has_surface = (is_conversion ||
                                  code == ConversionError::bad_pressure ||
                                  code == ConversionError::bad_skin_temperature ||
                                  code == ConversionError::bad_q_surf);
        const bool has_q_air = (code == ConversionError::bad_q_air ||
                                code == ConversionError::bad_q_surf);

        // Assemble the whole report before writing it, so that ranks that failed do not
        // interleave their lines with one another
        std::ostringstream msg;
        msg << who << ": surface conversion failed at level " << lev
            << " [rank " << amrex::ParallelDescriptor::MyProc() << "]\n"
            << "  " << fail.count << " of the cells this fill owns failed";
        if (fail.halo_count > 0) {
            msg << " (and " << fail.halo_count << " halo cells, which are not fatal)";
        }
        msg << "\n  first failure at (" << fail.i << "," << fail.j << "," << fail.k << "): "
            << conversion_error_text(code) << "\n";

        if (has_state) {
            msg << std::setprecision(10) << "    rho = " << fail.rho
                << "   (rho theta) = " << fail.rho_theta;
            if (has_vapor)   { msg << "   (rho qv) = " << fail.rho_qv; }
            if (has_q_air)   { msg << "   q_air = " << fail.q_air; }
            if (has_surface) {
                msg << "\n    surface temperature = " << fail.t_surf
                    << "   pressure = " << fail.pressure
                    << "   delta_z = " << fail.delta_z;
            }
            msg << "\n";
        }

        amrex::AllPrint() << msg.str() << std::flush;
    }

    // Let every rank finish writing before one of them tears the job down
    amrex::ParallelDescriptor::Barrier();
    amrex::Abort(std::string(who) + ": could not convert the surface state (see above)");
}

} // namespace erf_surface_temperature

#endif
