/**
 * @file ERF_Dust.H
 * @brief Main include file for the ERF-Dust module.
 *
 * Provides all public interfaces and inline functions for dust emission modeling.
 * Includes all dust headers in proper order.
 */

#ifndef ERF_DUST_H
#define ERF_DUST_H

#ifdef ERF_USE_DUST

#include <ERF_DustParams.H>
#include <ERF_DustGrid.H>
#include <ERF_DustLayer.H>
#include <ERF_DustPrerequisites.H>
#include <ERF_DustSurfaceReader.H>
#include <ERF_PhreeqcReader.H>
#include <ERF_DustThreshold.H>
#include <ERF_DustEmission.H>
#include <ERF_DustBlastSchedule.H>
#include <ERF_DustSuppression.H>
#include <ERF_DustWindExtract.H>
#include <ERF_DustAtmCoupling.H>
#include <ERF_DustSettling.H>
#include <ERF_DustDeposition.H>
#include <ERF_DustAtmReturn.H>
#include <ERF_DustPM.H>
#include <ERF_DustNAAQSOutput.H>

#endif // ERF_USE_DUST

#endif // ERF_DUST_H
