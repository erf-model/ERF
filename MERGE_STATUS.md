# Merge Status: `origin/development` → `ERF-Hazard`

**Date:** 2026-07-21  
**Performed by:** GitHub Copilot Coding Agent  
**Merge commit (prior sync):** `ecf1ecaa` — "Merge branch 'erf-model:development' into ERF-Hazard" (2026-07-20)

---

## Executive Summary

`origin/development` has been **fully merged** into `ERF-Hazard`.  
The merge was completed on 2026-07-20 via commit `ecf1ecaa6990564f2f2717d8ae948faf951c4eb9`.

Running `git merge origin/development --no-commit --no-ff` on `ERF-Hazard` now returns:

```
Already up to date.
```

`ERF-Hazard` is **418 commits ahead** of `origin/development` and **0 commits behind**.  
There are no outstanding merge conflicts.

---

## What Was Merged (Group A — Auto-Resolved)

The previous merge commit (`ecf1ecaa`) successfully incorporated all development changes into `ERF-Hazard`, including:

### Build System
| File | Resolution |
|---|---|
| `CMake/BuildERFExe.cmake` | ✅ Development's `ERF_NearSurfaceDiagnostics` target added; ERF-Hazard's Fire/Dust/FireDust targets preserved |
| `CMakeLists.txt` | ✅ Both sides kept: ERF-Hazard's `ERF_ENABLE_FIRE`/`ERF_ENABLE_DUST` options + development additions |

### Documentation
| File | Resolution |
|---|---|
| `Docs/BuildDocs.sh` | ✅ Both sides merged |
| `Docs/sphinx_doc/Plotfiles.rst` | ✅ Development's land-surface diagnostics + ERF-Hazard's fire/dust plotfile catalog |
| `Docs/sphinx_doc/Inputs.rst` | ✅ Development parameters + ERF-Hazard fire model parameter section |
| `Docs/sphinx_doc/LandSurfaceDiagnostics.rst` | ✅ New land-surface diagnostics page from development |
| `Docs/sphinx_doc/CouplingToNoahMP.rst` | ✅ Both sides kept |
| `Docs/sphinx_doc/index.rst` | ✅ Both sides kept (Fire, DustModule theory entries added) |
| `Docs/sphinx_doc/scripts/check_plotfile2d_catalog.py` | ✅ Catalog checker from development included |
| `Docs/sphinx_doc/scripts/test_check_plotfile2d_catalog.py` | ✅ Tests for catalog checker included |

### Source — New Files Only (ERF-Hazard additions, not in development)
All **99 new source files** under `Source/Fire/`, `Source/Dust/`, `Source/FireDust/`, and `Source/Particles/ERF_DustPC.*` are ERF-Hazard additions with no counterpart in development — no conflicts possible.

### Test Cases (Exec/)
**214 new test input files** under `Exec/CanonicalTests/Dust/` and `Exec/CanonicalTests/Fire/` — all ERF-Hazard additions, no conflicts.

---

## Group B — Files That Required Careful Merge (Source/ modifications)

The following existing ERF source files were modified by **both** ERF-Hazard and development. These were resolved in the previous merge commit (`ecf1ecaa`). Manual review is recommended to verify correctness:

### Core ERF files
| File | ERF-Hazard changes | Development changes | Review status |
|---|---|---|---|
| `Source/ERF.H` | Added `m_fire_layer`, `m_DustLayer`, `m_fire_dust_coupling` members | Plotfile2D, land-surface diagnostic additions | ⚠️ Review recommended |
| `Source/ERF.cpp` | Fire/dust initialization calls | Land-surface diagnostics, near-surface fields | ⚠️ Review recommended |
| `Source/ERF_Constructors.cpp` | Fire/dust constructor setup | Minor cleanup | ⚠️ Review recommended |
| `Source/ERF_IndexDefines.H` | `RhoSmoke_comp`, `RhoDust_comp` index additions | Other index additions | ⚠️ Review recommended |

### I/O files
| File | ERF-Hazard changes | Development changes | Review status |
|---|---|---|---|
| `Source/IO/ERF_Plotfile2D.cpp` | Fire/dust plotfile fields | Land-surface 2D fields | ⚠️ Review recommended |
| `Source/IO/ERF_Plotfile2DCatalog.H` | Fire/dust catalog entries | Land-surface catalog entries | ⚠️ Review recommended |
| `Source/IO/ERF_Plotfile2DCatalog.cpp` | Fire/dust catalog entries | Land-surface catalog entries | ⚠️ Review recommended |
| `Source/IO/ERF_Plotfile2DFill.H` | Fire/dust fill routines | Land-surface fill routines | ⚠️ Review recommended |
| `Source/IO/ERF_Plotfile2DFill.cpp` | Fire/dust fill routines | Land-surface fill routines | ⚠️ Review recommended |
| `Source/IO/ERF_Plotfile2DMetadata.cpp` | Fire/dust metadata | Land-surface metadata | ⚠️ Review recommended |
| `Source/IO/ERF_Plotfile2DSampledLevel.cpp` | Fire/dust sampled level | Land-surface sampled level | ⚠️ Review recommended |
| `Source/IO/ERF_Plotfile2DWaterPath.H` | New file from development | — | ✅ No ERF-Hazard conflict |
| `Source/IO/ERF_Plotfile2DWaterPath.cpp` | New file from development | — | ✅ No ERF-Hazard conflict |
| `Source/IO/ERF_WriteScalarProfiles.cpp` | Fire/dust scalar profiles | Minor additions | ⚠️ Review recommended |
| `Source/IO/ERF_Checkpoint.cpp` | Fire/dust checkpoint state | Minor additions | ⚠️ Review recommended |

### Time Integration
| File | ERF-Hazard changes | Development changes | Review status |
|---|---|---|---|
| `Source/TimeIntegration/ERF_Advance.cpp` | Fire/dust advance calls | CPM flexibility, coupling improvements | ⚠️ Review recommended |
| `Source/TimeIntegration/ERF_AdvanceDycore.cpp` | Fire-related dycore changes | Minor changes | ⚠️ Review recommended |
| `Source/TimeIntegration/ERF_SlowRhsPost.cpp` | Fire/dust slow RHS terms | Minor additions | ⚠️ Review recommended |
| `Source/TimeIntegration/ERF_TI_slow_rhs_post.H` | Fire/dust slow RHS post header | Minor additions | ⚠️ Review recommended |
| `Source/TimeIntegration/ERF_TI_slow_rhs_pre.H` | Fire/dust slow RHS pre header | Minor additions | ⚠️ Review recommended |

### Diffusion / PBL
| File | ERF-Hazard changes | Development changes | Review status |
|---|---|---|---|
| `Source/Diffusion/ERF_ComputeTurbulentViscosity.cpp` | Fire-modified turbulence | Development turbulence improvements | ⚠️ Review recommended |
| `Source/Diffusion/ERF_EddyViscosity.H` | Fire-modified eddy viscosity | Development changes | ⚠️ Review recommended |
| `Source/PBL/ERF_ComputeDiffusivityMRF.cpp` | Dust-coupled MRF diffusivity | Development MRF improvements | ⚠️ Review recommended |
| `Source/PBL/ERF_PBLModels.H` | Dust PBL header | Development PBL header | ⚠️ Review recommended |

### Land Surface / Noah-MP
| File | ERF-Hazard changes | Development changes | Review status |
|---|---|---|---|
| `Source/LandSurfaceModel/ERF_LandSurface.H` | Fire/dust land surface additions | Development land-surface additions | ⚠️ Review recommended |
| `Source/LandSurfaceModel/Noah-MP/ERF_NOAHMP_Advance.cpp` | Fire/dust Noah-MP coupling | Development Noah-MP changes | ⚠️ Review recommended |
| `Source/LandSurfaceModel/Noah-MP/ERF_NOAHMP_Fields.H` | Fire/dust Noah-MP fields | Development fields | ⚠️ Review recommended |
| `Source/LandSurfaceModel/Noah-MP/ERF_NOAHMP_ResultPolicy.H` | Fire/dust result policy | Development result policy | ⚠️ Review recommended |
| `Source/LandSurfaceModel/Noah-MP/Make.package` | Fire/dust Noah-MP package | Development Make.package | ⚠️ Review recommended |
| `Source/DataStructs/ERF_TurbStruct.H` | Fire/dust turbulence struct | Development turbulence struct | ⚠️ Review recommended |

### Diagnostics
| File | ERF-Hazard changes | Development changes | Review status |
|---|---|---|---|
| `Source/Diagnostics/ERF_NearSurfaceDiagnostics.H` | New file added by development merge | — | ✅ From development |
| `Source/Diagnostics/ERF_NearSurfaceDiagnostics.cpp` | New file added by development merge | — | ✅ From development |
| `Source/Diagnostics/Make.package` | Added NearSurfaceDiagnostics entry | Development Make.package | ⚠️ Review recommended |

### Tests
| File | ERF-Hazard changes | Development changes | Review status |
|---|---|---|---|
| `Tests/Unit/CMakeLists.txt` | Fire unit test targets | Development unit test infrastructure | ⚠️ Review recommended |
| `Tests/Unit/Fire/CMakeLists.txt` | Fire unit tests | — | ✅ ERF-Hazard only |
| `Tests/Unit/IO/ERF_GTestPlotfile2D.cpp` | Fire/dust plotfile2D tests | Development plotfile2D tests | ⚠️ Review recommended |
| `Tests/Unit/LandSurfaceModel/Noah-MP/ERF_GTestNoahMPResultPolicy.cpp` | Fire/dust Noah-MP tests | Development Noah-MP tests | ⚠️ Review recommended |

---

## Action Items for Manual Review

The following Source/ files were merged automatically but should be manually reviewed by a developer familiar with both the ERF core and the ERF-Hazard fire/dust additions:

1. **`Source/ERF.H`** — Check that `m_fire_layer`, `m_DustLayer`, `m_fire_dust_coupling` members do not clash with new development members
2. **`Source/IO/ERF_Plotfile2DCatalog.cpp`** — Verify fire/dust catalog entries are compatible with development's expanded 2D land-surface catalog
3. **`Source/TimeIntegration/ERF_Advance.cpp`** — Verify fire/dust advance calls are correctly placed relative to development's CPM and coupling improvements
4. **`Source/PBL/ERF_ComputeDiffusivityMRF.cpp`** — Verify dust-coupled MRF diffusivity is correct with development's MRF improvements

---

## Branch Comparison Summary

| Metric | Value |
|---|---|
| Base branch | `origin/development` (`e52b32df`) |
| Target branch | `ERF-Hazard` (`fad88b2c`) |
| ERF-Hazard commits ahead | 418 |
| Development commits not yet in ERF-Hazard | 0 |
| Files changed (ERF-Hazard vs development) | 378 |
| New ERF-Hazard source files (Fire/Dust/FireDust) | 99 |
| New ERF-Hazard test cases | 214 |
| New ERF-Hazard documentation files | 20 |
| Merge conflicts outstanding | **0** |

---

## Conclusion

The `origin/development` branch has been fully incorporated into `ERF-Hazard`. No further merge action is required at this time. The ERF-Hazard branch is ready for a PR into development once the Source/ file manual reviews (listed above) have been completed by the development team.
