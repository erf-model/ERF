#ifndef ERF_CHECKPOINT_LSM_DATA_H_
#define ERF_CHECKPOINT_LSM_DATA_H_

#include <cstddef>
#include <string>
#include <vector>

/**
 * File names of the LSM data in a checkpoint.
 *
 * lsm_data[i] is written as "LsmData<i>". The Noah-MP snow and snowh fields
 * were added in the middle of that numbering, which would shift every later
 * index (the soil fields) by two and make an older checkpoint restart with the
 * wrong fields. They are therefore stored under their own names and left out of
 * the numbering, so every other field keeps the index it always had. A restart
 * from a checkpoint that predates them keeps their initial values.
 */
namespace erf_checkpoint_lsm_data {

inline bool
is_separately_named (const std::vector<std::string>& names, const std::size_t ivar)
{
    return ivar < names.size() && (names[ivar] == "snow" || names[ivar] == "snowh");
}

// The part of the file name after the level directory, for lsm_data[ivar].
inline std::string
file_stem (const std::vector<std::string>& names, const std::size_t ivar)
{
    if (is_separately_named(names, ivar)) {
        return "LsmData_" + names[ivar];
    }
    std::size_t index = ivar;
    for (std::size_t k = 0; k < ivar; ++k) {
        if (is_separately_named(names, k)) { --index; }
    }
    return "LsmData" + std::to_string(index);
}

// True when a checkpoint written before the field existed is allowed to lack it.
inline bool
may_be_absent (const std::vector<std::string>& names, const std::size_t ivar)
{
    return is_separately_named(names, ivar);
}

} // namespace erf_checkpoint_lsm_data

#endif
