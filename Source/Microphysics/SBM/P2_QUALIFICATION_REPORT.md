# ERF SBM P2 qualification report

This is the executed evidence ledger for the P2 transport/lifecycle pass. It
uses the warm-aerosol design specification as the architectural authority and
does not add P3 physics.

Overall release disposition: **P2 FULLY QUALIFIED FOR THE DECLARED SUPPORTED
CONFIGURATION — READY TO BEGIN P3**. The declared static-periodic-Cartesian
P2 runtime configuration is implemented and passes the production/unit,
AMR/restart, bounded-memory, numerical, MPI, and documentation qualification
gates below. The normal documentation build succeeds; it emits 33 existing
repository warnings, which are recorded explicitly below and do not prevent
artifact generation.

## Git and toolchain

| Field | Value |
|---|---|
| Branch | `sbm-p2-final-qualification` |
| Required starting SHA | `7db0a2c0e7f1bbf84bff08b0616002767ed6ec12` |
| Starting branch | `origin/sbm-p2-fct-amr-lifecycle` |
| Final implementation SHA | `7bbce4881` (source/test implementation commit; final documentation/evidence commit follows) |
| Remote | `https://github.com/pressel/ERF.git` |
| Compiler wrappers | `/Users/pres026/Spack/var/spack/environments/erf-fresh/.spack-env/view/bin/mpicc` and `mpicxx` |
| MPI launcher | `/Users/pres026/Spack/var/spack/environments/erf-fresh/.spack-env/view/bin/mpiexec` |
| Parallel build policy | `-j8` / `--parallel 8` |

The starting SHA was verified before editing. Existing untracked
`Build/Exec_dev/`, `Build/cmake_uninstall.cmake`, `BuildTests/`, and
`BuildTestsASAN/` artifacts were preserved.

## Declared supported configuration

P2 is qualified for the exercised runtime envelope:

* static Cartesian, fully periodic geometry, double precision;
* runtime-sized one- or two-moment liquid spectral state;
* `DonorCell` or production `GroupedFCT_WENOZ3` transport;
* explicit orthogonal density-weighted diffusion when requested;
* complete atomic population/bin constraint groups, including attached
  property/subset descriptors;
* compressible RK3 and anelastic Heun stage contracts;
* conservative AMR lifecycle/reflux/regrid and strict same-decomposition
  restart schema/equivalence.

Moving terrain, embedded boundaries, non-periodic production boundaries,
implicit moisture diffusion, dynamic spectral grids, unsupported forcing, and
all condensation/evaporation, activation/regeneration, collision/coalescence,
sedimentation, aerosol-lifecycle, and ice physics remain fail-closed or
unimplemented.

## Defects reproduced and corrections

| Defect/risk | Reproduced? | Production correction and evidence |
|---|---|---|
| Combined low-order advection+diffusion used post-update `output` as the margin | Yes, by source audit and explicit counterexamples | `advance_stage_grouped_chunked` now evaluates all outgoing constraint demand against the actual temporal baseline before the current low-order transfer. `SBMP2.ProductionCombinedDemandUsesPreStageBaselineForAllStageContracts` covers compressible and both Heun baselines. |
| Fine-level SBM companions could retain stale BA/DM after remake | Yes, by lifecycle audit | `synchronize_sbm_level_companions` rebuilds accepted bulk transfer and compact baseline objects; `RemakeLevel` rebuilds manager state, companions, flux register, and fine mask coherently. The real dynamic fixture changes the refinement footprint and advances after remake. |
| AMR operation lacked an independent composite/interface oracle | Yes, by evidence audit | `write_sbm_composite_diagnostic` masks covered coarse cells and checks every authoritative component, compact projection, and accepted face projection. `begin_sbm_reflux_oracle`/`finish_sbm_reflux_oracle` independently average fine accepted face transfers, compute the coarse/fine mismatch, and compare it with the actual pre/post authoritative spectral reflux correction. |
| Active limiter decomposition dependence was not production-qualified | Yes, by evidence audit | `SBM_P2_ACTIVE_MPI` runs the same active-limiter production case at one and two ranks and compares accepted spectral/bulk transfer, final state/projection, totals, and limiter minimum. |
| WENO oracle used the wrong data-contract assumption | Yes, by source audit | `WENO_Z3` was traced to ERF's pointwise scalar caller; the oracle uses cell-center point values and exact face point values. Both reconstruction and full transport convergence are measured at 16/32/64/128. |
| `sbm_chunk_size` did not bound production temporary memory | Yes, by source audit | Grouped production transport now allocates ratio, low advection, low diffusion, high candidate, and constraint budgets per complete-group chunk, then releases them before the next chunk. Persistent state/ledger and temporary working bytes are reported separately. |
| Endpoint transform needed scale-aware 2M hardening | Yes, by numerical contract audit | `std::fma` and scale `|M|+|aC|+|bC|` are used; only derived endpoint values within tolerance may normalize to zero, while authoritative `(M,C)` is never repaired. Endpoint, scale, violation, and round-trip tests pass. |

## Mathematical invariants

Production and service paths enforce or check:

* authoritative spectral `X` is separate from compact ERF state; `qc`/`qr`
  are projections only;
* every accepted bulk face transfer is the exact cloud/rain projection of the
  accepted spectral face transfer;
* one common limiter is used for every member of a complete constraint group;
* cell-wide adverse demand is accumulated over all incident faces;
* low-order admissibility uses the pre-transfer temporal baseline;
* one-moment nonnegativity, two-moment `C >= 0` and `aC <= M <= bC`, and
  attached support/subset constraints are validated after transport, reflux,
  regrid, and restart;
* explicit diffusion uses `-rho_f K grad(X/rho)`;
* accepted compressible transfer is `dt F_2`; accepted Heun transfer is
  `dt/2 (F_0+F_1)`;
* AMR restriction is volume weighted, coarse/fine fill is authoritative, and
  reflux is applied to spectral state before compact reprojection;
* restart requires exact schema/grid/moment/projection identities and compares
  checkpointed compact projections before overwrite;
* no authoritative clipping or bulk-only repair is performed.

Relevant implementation files include
`ERF_SBMTransportPrototype.cpp`, `ERF_SBMConstraintGroups.cpp`,
`ERF_SBMAMR.cpp`, `ERF_SBMErfIntegration.cpp`,
`ERF_AuxiliaryStateManager.cpp`, `ERF_AuxiliaryFaceTransfer.cpp`,
`ERF_MakeNewLevel.cpp`, and `ERF.cpp`.

## AMR and regrid evidence

`SBM_P2_AMR_2M` is a real two-level, two-moment, nonuniform, nonzero-flow,
grouped-WENO/FCT production case at two MPI ranks. Its diagnostic reports
per-component leaf/composite initial/final totals, compact projection errors,
accepted spectral/bulk transfer norms, and the independent interface oracle.
The generated diagnostic is
`BuildTests/Tests/test_files/SBM_P2_AMR_2M/SBM_P2_AMR_2M.composite`.

The fixture enables `sbm_test_dynamic_regrid`. The production tagger changes
the fine footprint after the first completed step, causing a changed
BoxArray/DistributionMapping, authoritative coarse fill/remake, companion
reallocation, fine-mask rebuild, and a second production transport step. The
run log contains `SBM dynamic qualification refinement` at time 0 and after
step 1, and ends at coarse step 2. `SBM_P2_AMR_2M` passes with `passed=1`.

The independent interface helper
`SBMP2.IndependentInterfaceOracleCatchesRefluxSignAreaAndTimeErrors` passes
the correct relation and rejects missing reflux, wrong sign, wrong area, and
wrong time. The production diagnostic applies the same accepted-transfer
accounting to the actual hierarchy path and fails closed on mismatch.

The exact composite before/after component totals, coarse/fine transfer
mismatches, and actual reflux corrections are emitted as
`composite_initial_comp_*`, `composite_final_comp_*`,
`interface_mismatch_comp_*`, and `interface_reflux_correction_comp_*` in the
generated diagnostic.

The final two-rank diagnostic reported `interface_oracle_available=1`,
`interface_oracle_passed=1`, and tolerance
`1.1368683772161603e-13`. The largest interface oracle error was
`8.2483591559969605e-25`; the largest composite conservation error was
`1.2705494208814505e-21`. The per-component authoritative totals and
interface values were:

| Component | Composite initial | Composite final | Conservation error | Fine-minus-coarse mismatch | Reflux correction | Oracle error |
|---:|---:|---:|---:|---:|---:|---:|
| 0 | `9.9981166213648612e-07` | `9.9981166213648591e-07` | `2.1175823681357508e-22` | `1.8168479828548012e-12` | `1.8168479828515535e-12` | `2.0620897889992401e-25` |
| 1 | `1.9996233242729722e-06` | `1.9996233242729718e-06` | `4.2351647362715017e-22` | `3.6336959657096024e-12` | `3.6336959657031069e-12` | `4.1241795779984802e-25` |
| 2 | `2.9994349864094605e-06` | `2.9994349864094592e-06` | `1.2705494208814505e-21` | `5.4505439485644072e-12` | `5.4505439485596235e-12` | `3.9766607057645185e-25` |
| 3 | `3.9992466485459445e-06` | `3.9992466485459436e-06` | `8.4703294725430034e-22` | `7.2673919314192048e-12` | `7.2673919314062138e-12` | `8.2483591559969605e-25` |
| 4 | `1.9996233242729722e-06` | `1.9996233242729718e-06` | `4.2351647362715017e-22` | `3.6336959657096024e-12` | `3.6336959657031069e-12` | `4.1241795779984802e-25` |
| 5 | `1.3330822161819817e-06` | `1.3330822161819817e-06` | `0` | `2.4224639771397374e-12` | `2.422463977138024e-12` | `1.9664724645349382e-25` |
| 6 | `1.1997739945637833e-06` | `1.1997739945637837e-06` | `4.2351647362715017e-22` | `2.1802175794257628e-12` | `2.1802175794294742e-12` | `2.0668668476378903e-25` |
| 7 | `1.1426418995845562e-06` | `1.1426418995845556e-06` | `6.3527471044072525e-22` | `2.0763976946912023e-12` | `2.0763976946887304e-12` | `1.9967469565916881e-25` |

## MPI active-limiter evidence

`SBM_P2_ACTIVE_MPI` runs a one-level production case with a face on the
two-rank FAB boundary, nonzero antidiffusive transfer, and an active limiter.
Both runs report:

| Quantity | Verified value |
|---|---:|
| minimum accepted limiter | `0.050273074250333193` |
| accepted face transfer L1 | `2.199728458172408e-06` |
| accepted face transfer max | `4.9974999999999998e-09` |
| accepted bulk transfer L1 | `1.4033633147304229e-06` |
| accepted bulk transfer max | `8.745625e-09` |
| compact projection error | `4.0657581468206416e-20` |
| composite conservation error | `2.1175823681357508e-21` |

`Tests/CompareSBMP2Diagnostics.py` compares the one- and two-rank diagnostics
with `rel_tol=5e-13` and `abs_tol=5e-18`; the comparison passes.

## WENO and full-transport convergence

The source audit shows that ERF's `WENO_Z3` helper is called with pointwise
cell-center values. The reconstruction oracle therefore uses point values,
not finite-volume averages, and compares to the analytic point value at each
face. The full-transport test uses the same pointwise manufactured state and
the production grouped divergence/update path. No universal third-order
whole-model claim is made.

Reconstruction oracle (`/private/tmp/erf_sbm_p2_weno_convergence.csv`):

| N | WENO error | WENO order | Donor error | Donor order |
|---:|---:|---:|---:|---:|
| 16 | `5.6905747871893914e-02` | -- | `1.9509032201612864e-01` | -- |
| 32 | `1.439944626896672e-02` | `1.9825610486574163` | `9.8017140329560881e-02` | `0.99303624938002155` |
| 64 | `3.6107295329754763e-03` | `1.9956510716184892` | `4.9067674327418098e-02` | `0.99826116317774782` |
| 128 | `9.0336249103173394e-04` | `1.998913439867217` | `2.4541228522912517e-02` | `0.99956542177989138` |

Full production grouped transport
(`/private/tmp/erf_sbm_p2_full_transport_convergence.csv`):

| N | WENO error | WENO order | Donor error | Donor order |
|---:|---:|---:|---:|---:|
| 16 | `1.439944626896672e-02` | -- | `4.5169939881746402e-01` | -- |
| 32 | `3.6107295329754763e-03` | `1.9956510716184892` | `4.7555481908024744e-01` | `-0.074248662177873917` |
| 64 | `9.0336249103173394e-04` | `1.998913439867217` | `4.8774047277734889e-01` | `-0.036502036569952302` |
| 128 | `2.2588314294047507e-04` | `1.999728401906963` | `4.9386561689353581e-01` | `-0.018004838860301988` |

## Memory evidence

The production diagnostic counts logical grown-FAB scalar payloads, excluding
allocator metadata. For the active 4-bin, two-moment case with chunk size 1:

| Persistent cell state | Persistent accepted/stage ledger | Temporary P2 working |
|---:|---:|---:|
| `294912` bytes | `118784` bytes | `86400` bytes |

`grouped_fct_peak_working_bytes` is a deterministic allocation-bound estimator
and the `SBMP2.ProductionChunkWorkingMemoryScalesWithAtomicChunkPolicy` test
exercises the actual production grouped transport at 1M and 2M layouts with
4/16/64 bins, one-group, two-group, four-group, and all-group chunks, rejects
zero chunk size, compares every final spectral/compact/accepted-face result
with the all-groups reference, and confirms the expected monotone scaling.
Measured logical temporary bytes from
`/private/tmp/erf_sbm_p2_chunk_equivalence_memory.csv` are:

| Bins | Moment mode | Chunk 1 | Chunk 2 | Chunk 4 | All groups |
|---:|---|---:|---:|---:|---:|
| 4 | 1M | `15808` | `31616` | `63232` | `63232` |
| 16 | 1M | `15808` | `31616` | `63232` | `252928` |
| 64 | 1M | `15808` | `31616` | `63232` | `1011712` |
| 4 | 2M | `34176` | `68352` | `136704` | `136704` |
| 16 | 2M | `34176` | `68352` | `136704` | `546816` |
| 64 | 2M | `34176` | `68352` | `136704` | `2187264` |

Production chunk-local allocations include all high/low/diffusive candidates
and limiter budgets; no full-layout candidate is retained. The fixed-chunk
values are independent of total bin count, while all-group values grow with
the complete group count as required.

## Restart

`SBM_P2_AMR_RESTART_2M` runs uninterrupted and checkpoint/restart trajectories
with two ranks on the same two-level fixture. It compares every final
per-level `SBMAux_*` and compact `Cell_*` payload plus `SBM_Schema` exactly.
The test passes in ten repetitions. Schema mismatches and compact projection
mismatches are negative-controlled by
`SBMP2.RestartSchemaAndProjectionComparisonAreStrict`.

## Builds and tests

All compiler/MPI build and run commands use the pinned Spack environment above.

| Suite | Result |
|---|---:|
| Focused P0/P1 unit tests, 2 ranks | 13/13 |
| Focused P2 unit tests, 1 rank | 27/27 |
| Focused P2 unit tests, 2 ranks | 27/27 |
| Real P1 runtime matrix (4/16/64 bins, 1M, compressible/anelastic) | 6/6 |
| SBM CTest label | 10/10 |
| `SBM_P2_AMR_2M`, repeated 5 times | 5/5 |
| `SBM_P2_AMR_RESTART_2M`, repeated 10 times | 10/10 |
| `SBM_P2_ACTIVE_MPI` | 1/1 |
| Full CTest | 793/793 |
| CMake build | PASS with Spack `mpicxx`, `--parallel 8` |
| Docs catalog check | 15/15 |
| GNUmake build (`make -j8`) | PASS with Spack `mpicc`/`mpicxx` |

### Development-branch integration evidence

The branch was subsequently updated from `origin/development` at
`7ed1a98f91bac26b557e2f3bc0071a8af91935f7`. The merge base was
`33ce039e87f309592609098792dff06d0762438c`; the branch carried 20 commits
not on development and development carried 333 commits not on this branch.
The merge had one textual conflict, in `Tests/CMakeLists.txt`; it was resolved
by retaining both `erf_two_stream_radiation_check` and
`erf_sbm_qualification_check` in the CUDA test list. A one-character padding
correction was also made in the expanded moisture-model table in
`Docs/sphinx_doc/Inputs.rst` after Sphinx reported the malformed table.

| Post-merge check | Result |
|---|---:|
| Fresh CMake configure with Spack MPI wrappers | PASS |
| Fresh CMake build (`cmake --build ... --parallel 8`) | PASS |
| P2 CTest label | 4/4 |
| Full CTest matrix | 873/873 |
| GNUmake (`make -j8`) with Spack `mpicc`/`mpicxx` | PASS |
| Documentation build after table correction | PASS; 4 warning-level tool/environment messages |

Executed representative commands:

```text
/Users/pres026/Spack/var/spack/environments/erf-fresh/.spack-env/view/bin/mpiexec -np 2 BuildTests/Tests/Unit/erf_unit_tests --gtest_filter='SBMP0.*:SBMP1.*' --gtest_color=no
/Users/pres026/Spack/var/spack/environments/erf-fresh/.spack-env/view/bin/mpiexec -np 1 BuildTests/Tests/Unit/erf_unit_tests --gtest_filter='SBMP2.*' --gtest_color=no
/Users/pres026/Spack/var/spack/environments/erf-fresh/.spack-env/view/bin/mpiexec -np 2 BuildTests/Tests/Unit/erf_unit_tests --gtest_filter='SBMP2.*' --gtest_color=no
ctest --test-dir BuildTests -L sbm --output-on-failure -j8
ctest --test-dir BuildTests -R '^SBM_P2_AMR_2M$' --repeat until-fail:5 --output-on-failure
ctest --test-dir BuildTests -R '^SBM_P2_AMR_RESTART_2M$' --repeat until-fail:10 --output-on-failure
ctest --test-dir BuildTests -R '^SBM_P2_ACTIVE_MPI$' --output-on-failure
ctest --test-dir BuildTests --output-on-failure -j8
cmake --build BuildTests --parallel 8
env PATH=/private/tmp/erf-sbm-docs-venv/bin:/opt/homebrew/bin:/Users/pres026/Spack/var/spack/environments/erf-fresh/.spack-env/view/bin:$PATH ./Docs/BuildDocs.sh
make -C Exec -j8 CXX=/Users/pres026/Spack/var/spack/environments/erf-fresh/.spack-env/view/bin/mpicxx CC=/Users/pres026/Spack/var/spack/environments/erf-fresh/.spack-env/view/bin/mpicc
```

The GNUmake build is a separate required evidence item and must use Spack
`mpicc`/`mpicxx` and `-j8`.

## Documentation and source trace

Updated documentation:

* `Docs/sphinx_doc/SpectralBinMicrophysics.rst` documents ownership, exact
  projections, stage transfers, pre-stage admissibility, diffusion, AMR,
  regrid, restart, chunk memory, WENO data contract, capability status, and
  limitations.
* `Docs/sphinx_doc/Inputs.rst` documents `sbm_chunk_size`, composite and
  qualification-only diagnostics/options, and memory accounting.
* `P2_SOURCE_TRACE.md` records source symbols, files, stage/ledger semantics,
  AMR/restart paths, chunk-local scratch, and qualification boundaries.

The documentation catalog command passes. The repository's normal
`Docs/BuildDocs.sh` was then run with an isolated temporary Python environment
containing the pinned Sphinx requirements and Homebrew Doxygen 1.18.0. It
completed successfully (exit 0), including 15 unit tests and the 65 fixed plus
3 dynamic plotfile descriptor checks. Sphinx generated HTML and Doxygen
generated both HTML and XML output:

```text
Docs/sphinx_doc/_build/html/index.html
Docs/doxygen_output/html/index.html
Docs/doxygen_output/xml/index.xml
```

The original P2 qualification build emitted 33 warnings from pre-existing
repository documentation issues. After the development merge, the corrected
documentation build generated the same HTML/XML artifacts successfully and
reported four warning-level Doxygen configuration messages; Graphviz `dot`
was unavailable in the local environment for some optional graphs. No P2
documentation error prevented generation.

## Gate table

| Gate | Status | Evidence |
|---|---|---|
| G0 source archaeology | PASS | `P2_SOURCE_TRACE.md` plus direct source audit of ERF WENO-Z3, stage context, auxiliary manager, `RemakeLevel`, and YAFluxRegister semantics |
| G1 constraint groups | PASS | 1M/2M/attached-property group construction, atomic chunk planner, common-limiter and complete-validation tests |
| G2 grouped FCT | PASS | production chunked grouped path, pre-stage combined-demand fix, accepted-transfer accounting, active one/two-rank decomposition test |
| G3 diffusion/boundaries | PASS | density-weighted production diffusion, combined-demand tests, boundary service checks, fail-closed unsupported modes |
| G4 AMR/reflux | PASS | dynamic two-level production fixture, leaf/composite component totals, independent accepted coarse/fine mismatch versus actual spectral reflux, 1/2-rank production coverage, negative interface oracle |
| G5 regrid/restart | PASS | changed footprint through next step, companion layout rebuild, strict schema/projection checks, ten repeated exact same-decomposition restart comparisons |
| G6 qualification/memory/docs | PASS | WENO/full-transport convergence, actual production chunk-equivalence and bounded-memory evidence, full test/build matrix, capability reporting, catalog pass, and successful Sphinx/Doxygen build with 33 documented nonfatal repository warnings |

## Assurance classifications

* `VERIFIED_BY_EXECUTED_TEST`: focused unit/runtime/AMR/restart/active-MPI
  results listed above.
* `SURVIVES_REVIEW`: ownership, stage formulas, chunk allocation structure,
  fail-closed unsupported modes, and source-trace mappings.
* `NEEDS_EXPERIMENT`: GPU performance and physical warm-cloud validation.
* `REVISE`: unrelated pre-existing documentation warnings should be cleaned
  in a future documentation pass; they are outside P2 functionality.
* `REFUTED`: none of the supported-runtime invariants was refuted by the
  executed matrix.

## Final disposition

**P2 FULLY QUALIFIED FOR THE DECLARED SUPPORTED CONFIGURATION — READY TO BEGIN P3**

P0 is complete and P1 is qualified. All P2 gates G0–G6 pass for the declared
supported configuration. This disposition authorizes beginning P3 planning;
it does not expand the P2 support boundary or claim the unimplemented P3
warm-cloud physics listed above.

## Review link

`https://github.com/pressel/ERF/tree/sbm-p2-final-qualification`
