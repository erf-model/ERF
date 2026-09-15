# ERF SBM P2 qualification report

This report is the evidence ledger for the P2 implementation branch.  It
distinguishes exercised evidence from capabilities that remain fail-closed.
The attached implementation prompt is the specification; this report does
not authorize P3 physics.

Overall qualification status: **P2 PARTIAL**.  The qualified single-level
transport/diffusion path and the auxiliary AMR/restart service contracts pass;
G4/G5 remain `STRONGLY_SUPPORTED` because this repository has no dedicated
full ERF two-level continuous-vs-restart fixture.  The production capability
gate remains fail-closed for that missing qualification and for unsupported
geometry/boundary combinations.

## Git

| Field | Value |
|---|---|
| Branch | `sbm-p2-fct-amr-lifecycle` |
| Base SHA | `f7dd387dfa148bb2cd2d5a58f15f7aab88ba7ab5` (`sbm-p1-qualification`) |
| Final implementation SHA | `da35e21aa` (completed code, test, and evidence commit) |
| Remote | `https://github.com/pressel/ERF.git` |
| Policy | Push this branch only; do not merge into `development` |

## Gate status

| Gate | Status | Evidence |
|---|---|---|
| G0 source archaeology | PASS | `P2_SOURCE_TRACE.md`; ERF WENO-Z3 and pinned YAFluxRegister semantics traced from source |
| G1 constraint groups | PASS | `SBMP2.ConstraintGroupsCoverOneAndTwoMomentLayouts`, endpoint and synthetic-property tests |
| G2 grouped FCT | PASS | host reference, production `GroupedFCT_WENOZ3`, 15-test MPI suite |
| G3 diffusion/boundaries | PASS | density-weighted diffusion and boundary budget tests; production explicit timestep bound |
| G4 AMR lifecycle/reflux | STRONGLY_SUPPORTED | per-level manager lifecycle, provider register wiring, reflux validation, host AMR tests, and the two-rank two-level 2M production fixture pass; full subcycling/restart equivalence evidence remains limited |
| G5 regrid/restart | STRONGLY_SUPPORTED | remake/destroy lifecycle and strict schema/projection services pass; full continuous-vs-restart AMR fixture is not present |
| G6 qualification/docs | PASS | focused MPI suite, P1 regression suite, CTest/build/doc checks recorded below |

## Supported capability matrix

Actually qualified in this branch:

* runtime-sized 1M and 2M layout/constraint algebra and production grouped
  transport;
* `DonorCell` reference transport and production `GroupedFCT_WENOZ3` transport
  on the qualified periodic Cartesian path;
* complete one-lambda population/bin FCT groups, including attached subset
  constraints;
* compressible RK3 and anelastic Heun stage contracts;
* explicit orthogonal density-weighted two-point diffusion with an explicit
  timestep bound;
* periodic and generic prescribed-inflow/advective-outflow/impermeable-wall
  transfer-budget descriptors;
* static Cartesian auxiliary level creation, conservative coarse/fine service
  views, provider-owned accepted-transfer register wiring, post-reflux
  validation, remake, and destruction;
* exact P2 schema comparison and compact-projection-before-overwrite restart
  ordering.

The generic boundary descriptors and the AMR/restart services are qualified at
the service/unit level.  The production SBM gate remains periodic Cartesian,
static, double precision, and rejects unsupported host boundary/geometry
combinations rather than guessing a boundary spectrum.

## Claim ledger

| Claim ID | Requirement | Implementation path | Test/reproducer | Result | Evidence | Remaining limitation |
|---|---|---|---|---|---|---|
| P2-FCT-01 | Conservative grouped FCT | `ERF_SBMFCT.cpp`, `ERF_SBMTransportPrototype.cpp` | `GroupedFCTUsesOneFaceLimiterAndConservesEveryComponent`; production WENO test | accepted transfer is conservative and common-limited | VERIFIED | Production limiter uses local FAB face budgets; host reference carries global deterministic budgets |
| P2-INV-02 | 1M positivity | constraint groups and post-state checks | `ConstraintGroupsCoverOneAndTwoMomentLayouts`; P1 manufactured runs | pass | VERIFIED | Physical warm-cloud source terms are P3+ |
| P2-INV-03 | 2M realizability | physical `(M,C)` plus bounded `(L,H)` scratch | endpoint, two-moment FCT, production build | pass | VERIFIED | Multi-population production transport remains intentionally outside the current single-population host adapter |
| P2-PROP-04 | Property constraints | `make_constraint_groups`, POD production descriptors | subset and synthetic population tests | `0 <= subset <= carrier` | VERIFIED | No physical ice process |
| P2-TIME-05 | Heun temporal weighting | `StageWeightContract`, ledger | `FCTStageWeightsAndAcceptedCorrectionAreExact` | pass | VERIFIED | End-to-end anelastic grouped-FCT fixture not yet in CTest |
| P2-TIME-06 | Compressible old baseline | `AuxiliaryStateManager`, production stage recurrence | P1 stage/transport regressions and production WENO test | pass | STRONGLY_SUPPORTED | Existing real regression is donor-focused |
| P2-DIFF-07 | Diffuse `X/rho` | `ERF_SBMDiffusion`, production face flux | `DensityWeightedDiffusionUsesIntensiveRatioAndPhysicalGeometry` | pass | VERIFIED | Coefficient is explicit SBM input, not automatic PBL diffusivity |
| P2-PROJ-08 | Accepted bulk projection | `SBMBulkProjection`, accepted ledger | P1 transfer-closure tests; production build | pass | VERIFIED | Compact projection is liquid-population only by design |
| P2-OWN-09 | Unique face ownership | face-centered ledger and host duplicate rejection | duplicate ownership negative control; 2-rank production test | pass | STRONGLY_SUPPORTED | Full AMR shared-interface multi-rank fixture remains limited |
| P2-CHUNK-10 | Chunk-safe semantics | runtime layout, atomic groups, chunk policy | chunk-size invariance assertion | pass to roundoff | STRONGLY_SUPPORTED | Production FAB path currently keeps descriptor groups resident rather than implementing a multi-pass chunk kernel |
| P2-AMR-11 | Restriction/prolongation | manager and `ERF_SBMAMR` | manager lifecycle plus volume/register tests and `SBM_P2_AMR_2M` | pass | STRONGLY_SUPPORTED | Full subcycling and hierarchy scientific diagnostics remain limited |
| P2-AMR-12 | Reflux fail closed | `validate_post_reflux`, `ERF::post_timestep` | inadmissible correction negative control | collective diagnostic failure/no clip | VERIFIED | Full hierarchy negative control is service-level |
| P2-RESTART-13 | Strict restart schema | `ERF_SBMRestart`, checkpoint hooks | schema mismatch/projection negative controls | pass | STRONGLY_SUPPORTED | Full continuous-vs-restart hierarchy fixture remains to be added |
| P2-NONSBM-14 | Ordinary ERF unchanged | provider allocation/gates conditional on SBM | existing full CTest and non-SBM regressions | pass on existing suite | STRONGLY_SUPPORTED | No binary-diff claim is made |

## Numerical evidence

The focused unit executable asserts admissibility and conservation; its
reported tolerances and the directly measured values are:

* 1M grouped-FCT accepted mass total: `4.0` before and after correction,
  error below `1e-14` in the two-cell host reference;
* production grouped-WENO 1M and 2M cases at 4/16/64 bins: mass and (2M)
  number totals agree to the test tolerance `2e-12*max(1,|initial|)`;
  every final 2M endpoint margin is at least the checked roundoff allowance
  `-2e-14`;
* 2M and subset states: all endpoint/support/subset constraints admissible;
* density-weighted diffusion: integrated transfer `-0.45` and updated states
  `2.1125`, `7.55` for the non-unit geometry case;
* AMR transfer adapter: `I=12`, `A=3`, `dt=2` maps to per-area flux `2`;
* post-reflux negative control: rejected with level `2`, cell `0`, and a
  populated group/constraint diagnostic;
* restart projection comparison: `1` versus `1+1e-14` accepted, `1` versus
  `1+1e-4` rejected under the scale-aware comparison.

The six real ERF P1 manufactured cases additionally enforce finite scaled
mass, projection, face-projection, and transfer-closure tolerances through
`Tests/SBMQualificationCheck.cpp`; exact per-case diagnostic values are kept
in the generated `BuildTests/Tests/test_files/SBM_P1_*` directories.

## Convergence evidence

No universal third-order end-to-end ERF claim is made.  The production P2
regression demonstrates WENO-Z3 reuse and materially different smooth
transport from donor transport while the existing P1 cases remain the donor
reference.  A formal grid-refinement order table is not generated by the
current routine fixture and is therefore not marked verified.

## MPI/AMR evidence

* MPI ranks: 1 and 2 for the focused P2 unit suite; the production P2 WENO/FCT
  1M and 2M tests pass at both rank counts.
* AMR refinement: manager/reference tests use refinement ratio `(2,2,2)` and
  unequal fine/coarse cell volumes in the host transfer tests.
* Changed decomposition: manager remake uses a changed fine BoxArray split.
* Production coarse/fine crossing: `SBM_P2_AMR_2M` runs two levels with ratio
  `(2,2,2)`, `TwoWay` coupling, grouped WENO/FCT, 2M physical storage, and two
  MPI ranks through coarse step 2; the CTest checker verifies layout identity,
  2M component count, and successful completion after reflux.
* True level-subcycling and full continuous-vs-restart equivalence are not yet
  covered by a dedicated production fixture, so G4 remains
  `STRONGLY_SUPPORTED` and G5 remains `STRONGLY_SUPPORTED`.

## Memory policy

The diagnostic reports logical grown-FAB payloads, not allocator metadata.
For a representative single 4-cell tile, double precision, auxiliary cell
state with four states (old/evaluation/output/scratch), one accepted and one
stage face ledger, and the compact two-component baseline snapshot, the
logical payloads are approximately:

| bins | 1M authoritative state | 2M authoritative state | transport scratch | accepted face ledger | AMR register | chunk |
|---:|---:|---:|---:|---:|---:|---:|
| 4 | `4 * 4 * 8 * 4` bytes | `4 * 8 * 8 * 4` bytes | one bounded state allocation | two face-transfer objects | provider-owned YAFluxRegister | 256 |
| 16 | `4 * 16 * 8 * 4` bytes | `4 * 32 * 8 * 4` bytes | one bounded state allocation | two face-transfer objects | provider-owned YAFluxRegister | 256 |
| 64 | `4 * 64 * 8 * 4` bytes | `4 * 128 * 8 * 4` bytes | one bounded state allocation | two face-transfer objects | provider-owned YAFluxRegister | 256 |

The table is a component-count scaling ledger, not a claim of allocator-level
performance.  For the existing 4-cell qualification tile, the corresponding
authoritative cell-state payloads are 512/2048/8192 bytes for 1M and
1024/4096/16384 bytes for 2M at 4/16/64 bins.  The production diagnostic also
records grown-box and face-ledger bytes for the six inherited P1 cases.  There
is no `MAX_BINS` production constant and no permanent duplicate full 2M
endpoint state.

## Tests/builds

Commands run with the established Spack wrappers:

```text
/Users/pres026/Spack/var/spack/environments/erf-fresh/.spack-env/view/bin/mpiexec -np 1 BuildTests/Tests/Unit/erf_unit_tests --gtest_filter='SBMP2.*' --gtest_color=no
/Users/pres026/Spack/var/spack/environments/erf-fresh/.spack-env/view/bin/mpiexec -np 2 BuildTests/Tests/Unit/erf_unit_tests --gtest_filter='SBMP2.*' --gtest_color=no
/Users/pres026/Spack/var/spack/environments/erf-fresh/.spack-env/view/bin/mpiexec -np 2 BuildTests/Tests/Unit/erf_unit_tests --gtest_filter='SBMP0.*:SBMP1.*' --gtest_color=no
ctest --test-dir BuildTests -L sbm --output-on-failure -j8
ctest --test-dir BuildTests --output-on-failure -j8
cmake --build BuildTests -- -j8
cmake --build Build -- -j8
make -j8 CC=/Users/pres026/Spack/var/spack/environments/erf-fresh/.spack-env/view/bin/mpicc CXX=/Users/pres026/Spack/var/spack/environments/erf-fresh/.spack-env/view/bin/mpicxx
python3 Docs/sphinx_doc/scripts/test_check_plotfile2d_catalog.py
```

Results at the current evidence checkpoint:

| Check | Result |
|---|---|
| Focused P2 unit tests | 15/15 at 1 MPI rank and 2 MPI ranks |
| P1 focused SBM unit subset | 13/13 on 2 MPI ranks |
| Real single-level P1 qualification | 6/6 |
| AMR qualification | 5/5: four targeted service/lifecycle tests plus `SBM_P2_AMR_2M` two-rank two-level production fixture |
| Restart/regrid qualification | 2/2 targeted schema/remake service tests; full continuous/restart fixture absent |
| SBM CTest subset | 8/8 |
| Full CTest | 779/779 |
| CMake production build | PASS with Spack `mpicxx`, `-j8` |
| GNUmake build | PASS with Spack `mpicxx`, `-j8` |
| Docs/catalog checks | 15/15 |
| Full Sphinx/Doxygen build | BLOCKED: `sphinx-build` and `doxygen` are not installed |

## Remaining limitations / P3+

The following remain explicitly unsupported or unqualified:

* P3 thermodynamic adapter;
* condensation/evaporation;
* aerosol activation/regeneration;
* collision/coalescence;
* sedimentation;
* physical ice processes;
* implicit moisture diffusion;
* SHOC/macrophysics condensate coupling;
* moving terrain;
* embedded boundaries;
* non-orthogonal diffusion;
* dynamic spectral-grid conversion;
* GPU performance qualification;
* cloud-chamber/LES physical validation.

## Review link

https://github.com/pressel/ERF/tree/sbm-p2-fct-amr-lifecycle
