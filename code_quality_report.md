# Code Quality Report: `simple_radiation` and `lsf_and_nudging` merges

Date: 2026-03-17  
Branch context: `dg/sdm_w_cold_processes`  
Merge commits reviewed:
- `aaa30c0f0` — Merge branch `simple_radiation`
- `c1bfe88cb` — Merge branch `lsf_and_nudging`

This report focuses on code hygiene (indentation/whitespace, naming, dead code) and flags correctness risks that stand out during review.

---

## 1) `simple_radiation` (`aaa30c0f0`)

### 1.1 Whitespace & indentation
- Tabs and trailing whitespace were introduced:
  - `Source/PhysicsInterfaces/Radiation/Simple/ERF_RadiationSimple.cpp:159`
  - `Source/PhysicsInterfaces/Radiation/Simple/ERF_RadiationSimple.cpp:160`
  - `Source/PhysicsInterfaces/Radiation/Simple/ERF_RadiationSimple.cpp:161` (trailing whitespace)
  - `Source/PhysicsInterfaces/Radiation/Simple/ERF_RadiationSimple.cpp:162`

Recommendation:
- Convert tabs to spaces and remove trailing whitespace in the affected block.

### 1.2 API usage & const-correctness
- Unused parameter:
  - `Source/PhysicsInterfaces/Radiation/Simple/ERF_RadiationSimple.H:17` (`const int& lev` is unused)
- Overly permissive parameter type:
  - `Source/PhysicsInterfaces/Radiation/Simple/ERF_RadiationSimple.H:18` takes `SolverChoice& sc` but appears to only read it; prefer `const SolverChoice&`.

Recommendation:
- Remove the unused `lev` parameter (or use it meaningfully), and tighten `sc` to `const SolverChoice&` if mutation is not required.

### 1.3 Style consistency
- Very long conditional on a single line:
  - `Source/PhysicsInterfaces/Radiation/Simple/ERF_RadiationSimple.H:23`
- Constructor ends with an extra semicolon:
  - `Source/PhysicsInterfaces/Radiation/Simple/ERF_RadiationSimple.H:24`

Recommendation:
- Wrap long conditionals to match project style and remove the extra semicolon.

---

## 2) `lsf_and_nudging` (`c1bfe88cb`)

### 2.1 Correctness risks (high priority)

#### 2.1.1 Time-indexing can go out-of-bounds
- In `LargeScaleForcingData::get_forcing_time_coeffs`, `next` is set to `num_times` when `curr == num_times - 1`:
  - `Source/DataStructs/ERF_LargeScaleForcingData.H:339`
  - `Source/DataStructs/ERF_LargeScaleForcingData.H:342`
- Callers appear to index `*_int_lsf_d[itime_next]` directly:
  - `Source/TimeIntegration/ERF_AdvanceDycore.cpp:245`
  - `Source/TimeIntegration/ERF_AdvanceDycore.cpp:420`

Risk:
- If `itime_next == num_times`, this is an out-of-range access for vectors sized `num_times`.

Recommendation:
- Clamp `next` to `num_times - 1` (and set `coeff_next = 0`, `coeff_curr = 1`) for “use last set” behavior.

#### 2.1.2 Tautological assertion / likely bug from shadowing
- Local variable shadows the member `num_times` and asserts `num_times == num_times` (always true):
  - `Source/DataStructs/ERF_LargeScaleForcingData.H:218`
  - `Source/DataStructs/ERF_LargeScaleForcingData.H:219`

Recommendation:
- Remove the tautology; assert against the member if needed (or rename the local to avoid shadowing).

#### 2.1.3 Potential invalid access for `k=0`
- When interpolation does not “set” at `k=0`, code uses previous element `k-1`:
  - `Source/DataStructs/ERF_LargeScaleForcingData.H:289`
  - `Source/DataStructs/ERF_LargeScaleForcingData.H:290`

Risk:
- For `k=0` this becomes `[-1]`.

Recommendation:
- Handle `k=0` explicitly (e.g., initialize from first valid forcing level or a defined default).

### 2.2 Missing / implicit includes
- File uses `std::istringstream` and `std::ifstream`, but does not include `<sstream>` / `<fstream>`:
  - `Source/DataStructs/ERF_LargeScaleForcingData.H:32` (istringstream usage)
  - `Source/DataStructs/ERF_LargeScaleForcingData.H:48` (ifstream usage)

Risk:
- Build relies on transitive includes; can break depending on compiler/headers.

Recommendation:
- Add the required standard headers explicitly.

### 2.3 Debug/dead code left in place
- Unconditional debug printing blocks (even if they `break` after first iteration):
  - `Source/DataStructs/ERF_LargeScaleForcingData.H:191`
  - `Source/DataStructs/ERF_LargeScaleForcingData.H:304`
- Large commented-out nudging implementation:
  - `Source/SourceTerms/ERF_MakeMomSources.cpp:546`
- Disabled block via `#if 0`:
  - `Source/SourceTerms/ERF_MakeMomSources.cpp:712`
- Stray placeholder comment:
  - `Source/TimeIntegration/ERF_AdvanceDycore.cpp:229` (`/**/`)

Recommendation:
- Remove dead/commented blocks or convert to guarded debug (`#ifdef ERF_DEBUG`), and remove stray placeholders.

### 2.4 Unused variables / unused parameters (warnings)
- Unused local `inv_scale`:
  - `Source/TimeIntegration/ERF_AdvanceDycore.cpp:237`
- Unused parameter `lsf_data` in `make_sources`:
  - `Source/SourceTerms/ERF_MakeSources.cpp:56` (parameter is present but not referenced elsewhere in the file)

Recommendation:
- Remove unused locals; for unused parameters either wire up usage, drop the parameter, or mark as intentionally unused (consistent with project conventions).

### 2.5 Whitespace, redundancy, and minor inconsistencies
- Trailing whitespace:
  - `Source/IO/ERF_Write1DProfiles_stag.cpp:49`
  - `Source/IO/ERF_Write1DProfiles_stag.cpp:777`
- Redundant conditional:
  - `Source/IO/ERF_Write1DProfiles_stag.cpp:34` and `Source/IO/ERF_Write1DProfiles_stag.cpp:51` (`if (NumDataLogs() > 1)` nested under itself)
- Profiling label mismatch:
  - `Source/IO/ERF_Write1DProfiles_stag.cpp:27` uses `"ERF::write_1D_profiles()"` inside `ERF::write_1D_profiles_stag`

Recommendation:
- Remove trailing whitespace, simplify redundant checks, and align profile labels with function names for trace clarity.

---

## Suggested remediation order
1. Fix `LargeScaleForcingData` indexing/edge cases (`get_forcing_time_coeffs`, `k=0` handling) to eliminate out-of-bounds risk.
2. Add missing standard includes to `ERF_LargeScaleForcingData.H`.
3. Remove or gate debug prints and dead/commented blocks.
4. Clean whitespace/tabs/trailing spaces in the newly added/modified files.
5. Resolve unused locals/params to eliminate compiler warnings.

