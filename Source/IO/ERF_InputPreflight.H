#ifndef ERF_INPUT_PREFLIGHT_H_
#define ERF_INPUT_PREFLIGHT_H_

namespace erf_preflight {
// MPI, if enabled, must already be initialized. Emits one report on rank zero
// and returns the same status on every rank. Never constructs an ERF simulation.
int run (int argc, char* argv[]);
}

#endif
