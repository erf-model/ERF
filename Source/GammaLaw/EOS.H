#ifndef _EOS_H_
#define _EOS_H_

#include <AMReX.H>
#include <AMReX_REAL.H>
#include <AMReX_Utility.H>
#include <AMReX_Arena.H>
#include <AMReX_ParmParse.H>

namespace EOS {

constexpr amrex::Real gamma = 1.4;
constexpr amrex::Real RU = 8.31446261815324e7;
constexpr amrex::Real RUC  = 1.98721558317399615845;
constexpr amrex::Real PATM = 1.01325e+06;
constexpr amrex::Real AIRMW = 28.97;

void init();

void close();


AMREX_GPU_HOST_DEVICE
AMREX_FORCE_INLINE
void
T2Ei(amrex::Real T, amrex::Real Ei)
{
  amrex::Real wbar = AIRMW;
  const amrex::Real Cv = RU / (wbar * (gamma - 1.0));
  Ei = Cv * T;
}

AMREX_GPU_HOST_DEVICE
AMREX_FORCE_INLINE
void
T2Cv(amrex::Real T, amrex::Real& Cv)
{
  amrex::Real wbar = AIRMW;
  Cv = RU / (wbar * (gamma - 1.0));
}

AMREX_GPU_HOST_DEVICE
AMREX_FORCE_INLINE
void
T2Cp(amrex::Real T, amrex::Real& Cp)
{
  amrex::Real cv;
  EOS::T2Cv(T, cv);
  Cp = gamma * cv;
}

AMREX_GPU_HOST_DEVICE
AMREX_FORCE_INLINE
void
T2Cpi(amrex::Real T, amrex::Real Cp)
{
  Cp = gamma * RU / (AIRMW * (gamma - 1.0));
}

AMREX_GPU_HOST_DEVICE
AMREX_FORCE_INLINE
void
RET2P(
  amrex::Real  R,
  amrex::Real& E,
  amrex::Real& T,
  amrex::Real& P)
{
  P = (gamma - 1.0) * R * E;
}

AMREX_GPU_HOST_DEVICE
AMREX_FORCE_INLINE
void
RP2Cs(amrex::Real R, amrex::Real P,  amrex::Real& Cs)
{
  Cs = std::sqrt(gamma * P / R);
}

AMREX_GPU_HOST_DEVICE
AMREX_FORCE_INLINE
void
RT2Cs(amrex::Real R, amrex::Real T, amrex::Real& Cs)
{
  amrex::Real P;
  amrex::Real wbar = AIRMW;
  amrex::Real Cv = RU / (wbar * (gamma - 1.0));
  amrex::Real E = Cv * T;
  EOS::RET2P(R, E, T, P);
  Cs = std::sqrt(gamma * P / R);
}

AMREX_GPU_HOST_DEVICE
AMREX_FORCE_INLINE
void
E2T(amrex::Real E, amrex::Real& T)
{
  amrex::Real poverrho;
  amrex::Real wbar = AIRMW;
  poverrho = (gamma - 1.0) * E;
  T = poverrho * wbar / RU;
}

AMREX_GPU_HOST_DEVICE
AMREX_FORCE_INLINE
void
H2T(amrex::Real H, amrex::Real& T)
{
  amrex::Real wbar = AIRMW;
  const amrex::Real Cv = RU / (wbar * (gamma - 1.0));
  T = H / (Cv * gamma);
}

AMREX_GPU_HOST_DEVICE
AMREX_FORCE_INLINE
void
RT2P(amrex::Real R, amrex::Real T, amrex::Real& P)
{
  amrex::Real Cv, E;
  EOS::T2Cv(T, Cv);
  E = Cv * T;
  P = (gamma - 1.0) * R * E;
}

AMREX_GPU_HOST_DEVICE
AMREX_FORCE_INLINE
void
RP2T(amrex::Real R, amrex::Real P, amrex::Real& T)
{
  amrex::Real wbar = AIRMW;
  T = P * wbar / (R * RU);
}

AMREX_GPU_HOST_DEVICE
AMREX_FORCE_INLINE
void
PT2RE(
  amrex::Real P, amrex::Real T, amrex::Real& R, amrex::Real& E)
{
  amrex::Real wbar = AIRMW;
  R = P * wbar / (RU * T);
  E = P / (R * (gamma - 1.0));
}

AMREX_GPU_HOST_DEVICE
AMREX_FORCE_INLINE
void
PT2R(amrex::Real P, amrex::Real T, amrex::Real &R) 
{
  R = P * (RU * T); 
}

AMREX_GPU_HOST_DEVICE
AMREX_FORCE_INLINE
void
RP2E(amrex::Real Rho, amrex::Real P, amrex::Real& E)
{
  E = P / (Rho * (gamma - 1.0));
}

AMREX_GPU_HOST_DEVICE
AMREX_FORCE_INLINE
void
T2Hi(amrex::Real T, amrex::Real Hi)
{
  amrex::Real wbar = AIRMW;
  const amrex::Real Cv = RU / (wbar * (gamma - 1.0));
  Hi = Cv * T * gamma;
}

AMREX_GPU_HOST_DEVICE
AMREX_FORCE_INLINE
void
S(amrex::Real& s)
{
  s = 1.0;
}

AMREX_GPU_HOST_DEVICE
AMREX_FORCE_INLINE
void
T2G(amrex::Real T, amrex::Real& G)
{
  G = gamma;
}

AMREX_GPU_HOST_DEVICE
AMREX_FORCE_INLINE
void
T2H(amrex::Real T, amrex::Real& H)
{
  amrex::Real wbar = AIRMW;
  const amrex::Real Cv = RU / (wbar * (gamma - 1.0)); 
  H = Cv * T * gamma;
}

AMREX_GPU_HOST_DEVICE
AMREX_FORCE_INLINE
void
RPE2dpdr_e(amrex::Real R, amrex::Real P, amrex::Real E, amrex::Real& dpdr_e)
{
  dpdr_e = (gamma - 1.0) * E;
}

AMREX_GPU_HOST_DEVICE
AMREX_FORCE_INLINE
void
RG2dpde(amrex::Real R, amrex::Real G, amrex::Real& dpde)
{
  dpde = (gamma - 1.0) * R;
}

}; // namespace EOS

#endif
