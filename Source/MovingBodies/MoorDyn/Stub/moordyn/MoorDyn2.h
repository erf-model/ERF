#ifndef ERF_MOORDYN_STUB_MOORDYN2_H
#define ERF_MOORDYN_STUB_MOORDYN2_H

/*
 * Stand-in for MoorDyn-C's <moordyn/MoorDyn2.h>: the subset of the v2 C API that ERF's line
 * coupling calls, with the same names, signatures, error codes and log levels, implemented by
 * ERF_MoorDynStub.cpp. It lets the coupling be built and tested where MoorDyn is not installed.
 * MoorDyn line and point ids are 1-based as in MoorDyn; node indices are 0-based. Bodies, rods,
 * waves other than external ones, line failures and the per-component force getters (weight,
 * bending, seabed) are not declared. ERF_MoorDynStub.cpp lists where the stub's results differ
 * from MoorDyn-C 2.7.1.
 */

#include <stddef.h>
#include <stdint.h>

#define DECLDIR

#define MOORDYN_DBG_LEVEL 0
#define MOORDYN_MSG_LEVEL 1
#define MOORDYN_WRN_LEVEL 2
#define MOORDYN_ERR_LEVEL 3
#define MOORDYN_NO_OUTPUT 4096

#define MOORDYN_SUCCESS 0
#define MOORDYN_INVALID_INPUT_FILE -1
#define MOORDYN_INVALID_OUTPUT_FILE -2
#define MOORDYN_INVALID_INPUT -3
#define MOORDYN_NAN_ERROR -4
#define MOORDYN_MEM_ERROR -5
#define MOORDYN_INVALID_VALUE -6
#define MOORDYN_NON_IMPLEMENTED -7
#define MOORDYN_UNHANDLED_ERROR -255

#ifdef __cplusplus
extern "C" {
#endif

typedef struct __MoorDyn* MoorDyn;
typedef struct __MoorDynLine* MoorDynLine;
typedef struct __MoorDynPoint* MoorDynPoint;

MoorDyn DECLDIR MoorDyn_Create (const char* infilename);
int DECLDIR MoorDyn_NCoupledDOF (MoorDyn system, unsigned int* n);
int DECLDIR MoorDyn_SetVerbosity (MoorDyn system, int verbosity);
int DECLDIR MoorDyn_SetLogFile (MoorDyn system, const char* log_path);
int DECLDIR MoorDyn_SetLogLevel (MoorDyn system, int verbosity);
int DECLDIR MoorDyn_Init (MoorDyn system, const double* x, const double* xd);
int DECLDIR MoorDyn_Init_NoIC (MoorDyn system, const double* x, const double* xd);
int DECLDIR MoorDyn_Step (MoorDyn system, const double* x, const double* xd, double* f, double* t, double* dt);
int DECLDIR MoorDyn_Close (MoorDyn system);

int DECLDIR MoorDyn_ExternalWaveKinInit (MoorDyn system, unsigned int* n);
int DECLDIR MoorDyn_ExternalWaveKinGetN (MoorDyn system, unsigned int* n);
int DECLDIR MoorDyn_ExternalWaveKinGetCoordinates (MoorDyn system, double* r);
int DECLDIR MoorDyn_ExternalWaveKinSet (MoorDyn system, const double* U, const double* Ud, double t);

int DECLDIR MoorDyn_GetNumberPoints (MoorDyn system, unsigned int* n);
MoorDynPoint DECLDIR MoorDyn_GetPoint (MoorDyn system, unsigned int c);
int DECLDIR MoorDyn_GetNumberLines (MoorDyn system, unsigned int* n);
MoorDynLine DECLDIR MoorDyn_GetLine (MoorDyn system, unsigned int l);
int DECLDIR MoorDyn_GetDt (MoorDyn system, double* dt);
int DECLDIR MoorDyn_SetDt (MoorDyn system, double dt);
int DECLDIR MoorDyn_Save (MoorDyn system, const char* filepath);
int DECLDIR MoorDyn_Load (MoorDyn system, const char* filepath);
int DECLDIR MoorDyn_Serialize (MoorDyn system, size_t* size, uint64_t* data);
int DECLDIR MoorDyn_Deserialize (MoorDyn system, const uint64_t* data);

int DECLDIR MoorDyn_GetLineN (MoorDynLine l, unsigned int* n);
int DECLDIR MoorDyn_GetLineNumberNodes (MoorDynLine l, unsigned int* n);
int DECLDIR MoorDyn_GetLineUnstretchedLength (MoorDynLine l, double* ul);
int DECLDIR MoorDyn_GetLineNodePos (MoorDynLine l, unsigned int i, double pos[3]);
int DECLDIR MoorDyn_GetLineNodeVel (MoorDynLine l, unsigned int i, double vel[3]);
int DECLDIR MoorDyn_GetLineNodeTen (MoorDynLine l, unsigned int i, double t[3]);
int DECLDIR MoorDyn_GetLineNodeDrag (MoorDynLine l, unsigned int i, double f[3]);
int DECLDIR MoorDyn_GetLineNodeForce (MoorDynLine l, unsigned int i, double f[3]);
int DECLDIR MoorDyn_GetLineFairTen (MoorDynLine l, double* t);
int DECLDIR MoorDyn_GetLineMaxTen (MoorDynLine l, double* t);

int DECLDIR MoorDyn_GetPointType (MoorDynPoint point, int* t);
int DECLDIR MoorDyn_GetPointPos (MoorDynPoint point, double pos[3]);
int DECLDIR MoorDyn_GetPointForce (MoorDynPoint point, double f[3]);

#ifdef __cplusplus
}
#endif

#endif
