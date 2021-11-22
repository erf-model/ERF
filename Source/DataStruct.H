#ifndef _DATA_STRUCT_H_
#define _DATA_STRUCT_H_

#include <ERF_Constants.H>
#include <AMReX_ParmParse.H>

enum class ABLDriverType {
    None, PressureGradient, GeostrophicWind
};

enum class Coord {
    x, y, z
};

enum class NextOrPrev {
    prev = -1,
    local = 0,
    next = 1
};

enum class AdvectedQuantity {
    unity, u, v, w, theta, scalar
};

enum class AdvectingQuantity {
    rho_u, rho_v, rho_w
};

enum class AdvectionDir {
    x, y, z
};

enum class DiffusionDir {
    x, y, z
};

enum class MomentumEqn {
    x, y, z
};

enum class LESType {
    None, Smagorinsky, Deardorff
};

enum class TurbulenceModel {
    DNS, Smagorinsky, Deardorff
};

#endif
