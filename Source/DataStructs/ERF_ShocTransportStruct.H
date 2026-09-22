#ifndef ERF_SHOC_TRANSPORT_STRUCT_H_
#define ERF_SHOC_TRANSPORT_STRUCT_H_

#include <AMReX_ParmParse.H>

#include <algorithm>
#include <cctype>
#include <string>

enum class ShocTransportMode
{
    StateUpdate
};

enum class ShocMomentumTransport
{
    None,
    StateUpdate
};

inline const char*
shoc_transport_mode_name (ShocTransportMode mode)
{
    switch (mode) {
    case ShocTransportMode::StateUpdate: return "state_update";
    }
    return "unknown";
}

inline const char*
shoc_momentum_transport_name (ShocMomentumTransport mode)
{
    switch (mode) {
    case ShocMomentumTransport::None:        return "none";
    case ShocMomentumTransport::StateUpdate: return "state_update";
    }
    return "unknown";
}

inline bool
shoc_uses_momentum_state_update (ShocMomentumTransport mode)
{
    return (mode == ShocMomentumTransport::StateUpdate);
}

inline bool
shoc_disables_momentum_transport (ShocMomentumTransport mode)
{
    return (mode == ShocMomentumTransport::None);
}

inline std::string
shoc_to_lower (std::string value)
{
    std::transform(value.begin(), value.end(), value.begin(),
                   [] (unsigned char c) { return static_cast<char>(std::tolower(c)); });
    return value;
}

inline bool
parse_shoc_transport_mode_string (std::string value,
                                  ShocTransportMode& mode,
                                  std::string& error_message)
{
    value = shoc_to_lower(value);
    if (value == "state_update") {
        mode = ShocTransportMode::StateUpdate;
        return true;
    }
    if (value == "host_diffusion") {
        error_message =
            "erf.shoc.transport_mode = host_diffusion has been removed for native SHOC. Use erf.shoc.transport_mode = state_update.";
        return false;
    }
    if (value == "tendencies") {
        error_message =
            "erf.shoc.transport_mode = tendencies has been removed for native SHOC. Use erf.shoc.transport_mode = state_update. Native SHOC now applies its coupled column increment before the dycore to avoid splitting theta and momentum from moisture.";
        return false;
    }

    error_message = "erf.shoc.transport_mode must be 'state_update'";
    return false;
}

inline bool
parse_shoc_momentum_transport_string (std::string value,
                                      ShocMomentumTransport& mode,
                                      std::string& error_message)
{
    value = shoc_to_lower(value);
    if (value == "none") {
        mode = ShocMomentumTransport::None;
        return true;
    }
    if (value == "state_update") {
        mode = ShocMomentumTransport::StateUpdate;
        return true;
    }
    if (value == "host_diffusion") {
        error_message =
            "erf.shoc.momentum_transport = host_diffusion has been removed for native SHOC. Use 'state_update' to retain SHOC momentum transport, or 'none' to disable SHOC momentum transport.";
        return false;
    }

    error_message = "erf.shoc.momentum_transport must be 'none' or 'state_update'";
    return false;
}

/**
 * @brief Read the native SHOC transport modes from the erf.shoc input namespace.
 *
 * This is the single parse point shared by SolverChoice and the SHOC runtime
 * options, so the dycore and the SHOC driver can never disagree about which
 * component owns transport. Values that are absent from the inputs leave the
 * incoming modes untouched.
 */
inline void
read_shoc_transport_modes (ShocTransportMode& transport_mode,
                           ShocMomentumTransport& momentum_transport)
{
    amrex::ParmParse pp("erf.shoc");

    std::string transport_mode_string;
    if (pp.query("transport_mode", transport_mode_string)) {
        std::string error_message;
        if (!parse_shoc_transport_mode_string(transport_mode_string,
                                              transport_mode,
                                              error_message)) {
            amrex::Abort(error_message.c_str());
        }
    }

    std::string momentum_transport_string;
    if (pp.query("momentum_transport", momentum_transport_string)) {
        std::string error_message;
        if (!parse_shoc_momentum_transport_string(momentum_transport_string,
                                                  momentum_transport,
                                                  error_message)) {
            amrex::Abort(error_message.c_str());
        }
    }
}
#endif
