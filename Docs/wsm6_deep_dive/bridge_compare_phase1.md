# WSM6 Bridge Compare (Phase 1)

## Purpose
This file is the short patch map between ERF's current WSM6 bridge and the canonical `physics_mmm/mp_wsm6.F90` interfaces, so follow-on work can be done with explainable incremental diffs.

## Entrypoint Mapping
- Canonical init: `mp_wsm6_init(den0,denr,dens,cl,cpv,hail_opt,errmsg,errflg)`
- Canonical run: `mp_wsm6_run(t,q,qc,qi,qr,qs,qg,den,p,delz,delt,g,cpd,cpv,rd,rv,t0c,ep1,ep2,qmin,xls,xlv0,xlf0,den0,denr,cliq,cice,psat,rain,rainncv,sr,snow,snowncv,graupel,graupelncv,rainprod2d,evapprod2d,its,ite,kts,kte,errmsg,errflg)`
- ERF C bridge entrypoints: `mp_wsm6_init_c`, `mp_wsm6_run_c`
- ERF Fortran shim entrypoints: `mp_wsm6_init`, `mp_wsm6_run`

## Status Table
| Area | Canonical (`physics_mmm`) | ERF phase-1 state | Status |
|---|---|---|---|
| Symbol names | `mp_wsm6_init`, `mp_wsm6_run` | Same names in ERF-owned module | Matched |
| C bindings | Not provided in canonical module | `mp_wsm6_init_c`, `mp_wsm6_run_c` in ERF isohelper | Matched (ERF bridge layer) |
| Init args | `den0,denr,dens,cl,cpv,hail_opt,+errmsg/errflg` | Same arg set in shim | Matched |
| Run core thermo/hydro args | `t,q,qc,qi,qr,qs,qg,den,p,delz,...` | Same core set in shim (`qv` naming on ERF side maps to canonical `q`) | Matched |
| 2-D precip diagnostics | `rain,rainncv,sr,snow,snowncv,graupel,graupelncv` | Same set in shim and C bridge | Matched |
| Optional 2-D diagnostics | `rainprod2d,evapprod2d` optional | Not in phase-1 shim/C bridge yet | Deferred |
| Run index contract | Canonical uses `(its,ite,kts,kte)` | ERF C bridge passes full `(ids..kte)` and shim currently keeps all | Stubbed for transition |
| Physics internals | Full WSM6 microphysics | ERF shim currently does validation + nonnegative clamps only | Stubbed |
| Error contract | `errmsg,errflg` out args | Preserved in shim; isohelper converts nonzero to `stop 1` | Matched |

## Files In ERF For This Bridge
- `Source/Microphysics/WSM6/ERF_module_mp_wsm6.F90`
- `Source/Microphysics/WSM6/ERF_module_mp_wsm6_isohelper.F90`
- `Source/Microphysics/WSM6/ERF_WSM6_Fortran_Interface.H`
- `Source/Microphysics/WSM6/ERF_AdvanceWSM6.cpp`

## Replace Order (Minimal-Risk)
1. Keep `mp_wsm6_*_c` signatures fixed; do not change C++ call sites yet.
2. Replace shim `mp_wsm6_init` body with canonical init internals.
3. Replace shim `mp_wsm6_run` body with canonical run internals, still keeping ERF C bridge ABI unchanged.
4. Add `rainprod2d/evapprod2d` as optional bridge args only after run body replacement is stable.
5. Reduce index argument surface to canonical `(its,ite,kts,kte)` once full run parity is validated.

## Quick Diff Commands
```bash
# Canonical vs ERF shim (module-level implementation)
diff -u /home/jmsexton/codes/physics_mmm/mp_wsm6.F90 \
        /home/jmsexton/codes/ERF/Source/Microphysics/WSM6/ERF_module_mp_wsm6.F90

# ERF C bridge surface used by C++
diff -u /home/jmsexton/codes/ERF/Source/Microphysics/WSM6/ERF_WSM6_Fortran_Interface.H \
        /home/jmsexton/codes/ERF/Source/Microphysics/WSM6/ERF_module_mp_wsm6_isohelper.F90
```
