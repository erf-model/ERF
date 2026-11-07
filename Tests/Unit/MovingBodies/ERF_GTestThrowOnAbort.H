#ifndef ERF_GTEST_THROW_ON_ABORT_H
#define ERF_GTEST_THROW_ON_ABORT_H

/**
 * Test helpers that turn amrex::Abort and the AMREX_ALWAYS_ASSERT family into a thrown
 * amrex::RuntimeError (a std::runtime_error), so a unit test can prove that an abort fires and
 * read its message. Only for code that aborts before it changes shared state: the throw unwinds
 * through the caller instead of ending the run.
 */

#include <functional>
#include <stdexcept>
#include <string>

#include <AMReX.H>

namespace erf_gtest {

/** While an object of this type lives, amrex::Abort and failed assertions throw instead of ending the run. */
class ThrowOnAbort
{
public:
    ThrowOnAbort () : m_old(amrex::system::throw_exception) { amrex::system::throw_exception = true; }
    ~ThrowOnAbort () { amrex::system::throw_exception = m_old; }
    ThrowOnAbort (const ThrowOnAbort&) = delete;
    ThrowOnAbort& operator= (const ThrowOnAbort&) = delete;
private:
    bool m_old;
};

/**
 * Run f with aborts turned into exceptions. Returns the abort's message, or an empty string when
 * f returned without aborting.
 */
inline std::string abort_message (const std::function<void()>& f)
{
    ThrowOnAbort guard;
    try {
        f();
    } catch (const std::runtime_error& e) {
        return e.what();
    }
    return {};
}

} // namespace erf_gtest

#endif
