#include "ERF_NoahMPGitHash.H"

#include <AMReX_buildInfo.H>

// The Noah-MP hash arrives by a different route in each build system. AMReX's
// CMake generate_buildinfo() has only two hash slots (ERF and AMReX), so CMake
// runs its own git describe and hands the result to this one file as a compile
// definition; GNU Make instead appends Noah-MP to the makebuildinfo_C.py --GIT
// list, where it lands in slot 3. Both yield "" when Noah-MP is not in the
// build. Keeping this in its own translation unit means a changed hash under
// CMake recompiles this file alone rather than all of ERF.
const char*
NoahMPGitHash ()
{
#if defined(ERF_NOAHMP_GIT_HASH)
    return ERF_NOAHMP_GIT_HASH;
#else
    return amrex::buildInfoGetGitHash(3);
#endif
}
