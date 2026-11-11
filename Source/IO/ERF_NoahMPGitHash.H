/**
 * \file ERF_NoahMPGitHash.H
 */
#ifndef ERF_NOAHMP_GIT_HASH_H_
#define ERF_NOAHMP_GIT_HASH_H_

/**
 * The git hash of the Noah-MP submodule this executable was built against.
 *
 * Returns "" when Noah-MP is not part of the build, so callers guard on
 * strlen() > 0 exactly as they do for the other build-info hashes.
 */
const char* NoahMPGitHash ();

#endif
