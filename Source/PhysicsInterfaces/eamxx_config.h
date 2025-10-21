#ifndef SCREAM_CONFIG_H
#define SCREAM_CONFIG_H

// If defined, Real is double; if not, Real is float.
#define SCREAM_DOUBLE_PRECISION

// If defined, enable floating point exceptions.
/* #undef SCREAM_FPE */

// The number of scalars in a scream::pack::Pack and Mask.
#define SCREAM_PACK_SIZE 1

// The number of scalars in a scream::pack::SmallPack and SmallMask.
#define SCREAM_SMALL_PACK_SIZE 1

// The number of scalars in a possibly-no-pack. Use this packsize when a routine does better with pksize=1 on some architectures (SKX).
#define SCREAM_POSSIBLY_NO_PACK_SIZE 1

// How many levels to use for the vertical grid
/* #define SCREAM_NUM_VERTICAL_LEV 72 */

// Whether this is a CUDA/HIP build
/* #undef EAMXX_ENABLE_GPU */

// Whether scream uses leap years or not
#define SCREAM_HAS_LEAP_YEAR

// What level of testing we are doing. 0=autotesting, 1=nightly, 2=experimental
#define SCREAM_TEST_LEVEL 0

// Whether getrusage can be used to get memory usage
/* #undef SCREAM_ENABLE_GETRUSAGE */
// Whether /proc/self/statm can be used to get memory usage
/* #undef SCREAM_ENABLE_STATM */

#if defined(SCREAM_ENABLE_STATM) || defined(SCREAM_ENABLE_GETRUSAGE)
#define SCREAM_HAS_MEMORY_USAGE
#endif

#define SCREAM_MPI_ON_DEVICE 1

// Data directory for the scream project
/*#define SCREAM_DATA_DIR "/home/alattanz/git/E3SM/components/eamxx/build/scream-input/atm/scream"*/

// Whether or not to run RRTMGP debug checks
/* #undef SCREAM_RRTMGP_DEBUG */

// Whether or not small kernels are used in ALL targets that support them
/* #undef SCREAM_SMALL_KERNELS */
// Whether or not small kernels are used in P3
/* #undef SCREAM_P3_SMALL_KERNELS */
// Whether or not small kernels are used in SHOC
/* #undef SCREAM_SHOC_SMALL_KERNELS */

#endif
