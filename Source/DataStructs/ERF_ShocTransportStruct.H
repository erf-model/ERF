#ifndef ERF_SHOC_TRANSPORT_STRUCT_H_
#define ERF_SHOC_TRANSPORT_STRUCT_H_

#include <algorithm>
#include <cctype>
#include <string>

enum class ShocTransportMode
{
    StateUpdate,
    HostDiffusion
};

enum class ShocMomentumTransport
{
    None,
    StateUpdate,
    HostDiffusion
};

inline const char*
shoc_transport_mode_name (ShocTransportMode mode)
{
    switch (mode) {
    case ShocTransportMode::StateUpdate:   return "state_update";
    case ShocTransportMode::HostDiffusion: return "host_diffusion";
    }
    return "unknown";
}

inline const char*
shoc_momentum_transport_name (ShocMomentumTransport mode)
{
    switch (mode) {
    case ShocMomentumTransport::None:          return "none";
    case ShocMomentumTransport::StateUpdate:   return "state_update";
    case ShocMomentumTransport::HostDiffusion: return "host_diffusion";
    }
    return "unknown";
}

inline bool
shoc_uses_state_update (ShocTransportMode mode)
{
    return (mode == ShocTransportMode::StateUpdate);
}

inline bool
shoc_uses_host_diffusion (ShocTransportMode mode)
{
    return (mode == ShocTransportMode::HostDiffusion);
}

inline bool
shoc_uses_momentum_state_update (ShocMomentumTransport mode)
{
    return (mode == ShocMomentumTransport::StateUpdate);
}

inline bool
shoc_uses_momentum_host_diffusion (ShocMomentumTransport mode)
{
    return (mode == ShocMomentumTransport::HostDiffusion);
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
        mode = ShocTransportMode::HostDiffusion;
        return true;
    }
    if (value == "tendencies") {
        error_message =
            "erf.shoc.transport_mode = tendencies has been removed for native SHOC. Use erf.shoc.transport_mode = state_update. Native SHOC now applies its coupled column increment before the dycore to avoid splitting theta and momentum from moisture.";
        return false;
    }

    error_message = "erf.shoc.transport_mode must be 'state_update' or 'host_diffusion'";
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
        mode = ShocMomentumTransport::HostDiffusion;
        return true;
    }

    error_message = "erf.shoc.momentum_transport must be 'none', 'state_update', or 'host_diffusion'";
    return false;
}
#endif
