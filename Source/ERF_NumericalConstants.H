/*! @file ERF_NumericalConstants.H
 *  \brief Dimensionless numeric literals and pure mathematical constants.
 *
 *  Nothing in this file is physics.  These are spelling aids that let the code
 *  write `myhalf` instead of `amrex::Real(0.5)` and so stay precision-agnostic
 *  between double and single (AMREX_USE_FLOAT) builds, plus the mathematical
 *  constant pi and the two combinations of it that recur often enough to be
 *  worth naming.
 *
 *  Kept separate from ERF_Constants.H so that a file needing only numeric
 *  literals does not have to pull in the thermodynamic constants, and so that
 *  ERF_Constants.H reads as a list of physical constants.  ERF_Constants.H
 *  includes this header, so every existing consumer keeps compiling unchanged.
 */

#ifndef ERF_NUMERICAL_CONSTANTS_H_
#define ERF_NUMERICAL_CONSTANTS_H_

#include <AMReX_REAL.H>

constexpr amrex::Real zero         = amrex::Real(0.0);
constexpr amrex::Real one          = amrex::Real(1.0);
constexpr amrex::Real two          = amrex::Real(2.0);
constexpr amrex::Real three        = amrex::Real(3.0);
constexpr amrex::Real four         = amrex::Real(4.0);
constexpr amrex::Real myhalf       = amrex::Real(0.5);
constexpr amrex::Real fourth       = amrex::Real(0.25);
constexpr amrex::Real third        = amrex::Real(1.0/3.0);
constexpr amrex::Real eightthirds  = amrex::Real(8.0/3.0);

constexpr amrex::Real PI           = amrex::Real(3.14159265358979323846264338327950288);
constexpr amrex::Real PIoTwo       = PI/two;
constexpr amrex::Real four_thirds_pi = four/three*PI;

#endif
