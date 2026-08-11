---
name: SDM+AMR missing features from ice branch
description: AMR-related features present in dg/ice_sdm_w_AMR but missing from dg/sdm_w_AMR, prioritized by correctness impact
type: project
---

Analysis date: 2026-04-14. Compared dg/sdm_w_AMR (HEAD 25b973bff) against dg/ice_sdm_w_AMR (HEAD 3ddd563e4).

## Critical (affect correctness)

1. **Particle splitting/merging on regrid**
   - `SplitParticlesForRefinement(lev, old_ba, ref_ratio)` — split particles when region gets refined so fine-level cells have adequate SD count
   - `MergeParticlesAfterDerefining(lev, old_ba, ref_ratio)` — merge particles when region de-refined
   - `SplitMergeAtLevelBoundary()` — handle particles crossing level boundaries during timestep
   - Declared as virtual no-ops in ERFPC base class, overridden in SuperDropletPC
   - Without this, fine-level cells are sparsely populated
   - Key commit: 245558fe3

2. **Fine-level coverage mask (`m_fine_level_mask`)**
   - Built in Advance() using `makeFineMask` when `current_lev < finest_level`
   - Skips Euler-field updates (phaseChange, Copy_Micro_to_State) in L0 cells covered by L1
   - Without this, L0 phase change double-counts cells that also have fine-level particles
   - Used extensively in PhaseChange and Copy_Micro_to_State

3. **Pre-regrid AverageDown + RestoreMoistVarsInState**
   - `AverageDownMicroVars(finest_lev)` — averages micro FABs (qc, qr) fine→coarse; fills L0 halo cells from covered neighbors
   - `RestoreMoistVarsInState_Lev(lev, cons)` — copies qc/qr from micro FABs back into conserved state
   - Called before `regrid()` in ERF_TimeStep.cpp so tagging sees particle-derived qc on L0
   - Called in plotfile writer so plotfiles show correct qc on all levels
   - Virtual methods added to base Microphysics class (ERF_Microphysics.H)
   - Key commit: 245558fe3, d89cdb952

## Important (affect solution quality)

4. **FillPatch before microphysics advance**
   - `FillPatchCrseLevel` (lev=0) or `FillPatchFineLevel` (lev>0) before `advance_microphysics`
   - Fills ghost cells in conserved state that are stale after dycore advance
   - In ERF_Advance.cpp

5. **Level-scaled coalescence bins**
   - `coal_bin_size *= refRatio(k)` for each level below current
   - Keeps physical bin volume constant across AMR levels
   - Without this, coalescence suppressed on fine levels (smaller cells → fewer particles per bin)
   - In ERF_SuperDropletPCCoalescence.cpp

## Robustness / QoL

6. **Redistribute before initial plotfile write** (5116e5067)
   - After L0-only particle init, Redistribute() needed before initial plotfile so container has valid structures at fine levels

7. **CountParticlesPerLevelAndHalo** — diagnostics after split/merge

8. **Terrain-aware FixKIndexAMR** (67b805ff5) — use z_phys_nd terrain array to correct k-index, not just flat-grid formula

9. **Refinement box index alignment** (ebae4ef9d) — auto-snap box indices to ref_ratio alignment instead of aborting
