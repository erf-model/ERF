# SubDyn reference: a small lattice tower

Reference stiffness for ERF's frame model (`Source/MovingBodies/Towers/ERF_Frame.H`), computed by
OpenFAST's SubDyn. The unit tests `FrameSubDyn.*` read these SubDyn input files with ERF's reader and
compare ERF's condensed stiffness at the peak joint with SubDyn's.

## The tower

- 4 legs from the base corners (±3, ±3, 0) m, tapering to (±1, ±1, 21) m.
- Panel levels at z = 0, 7, 14 and 21 m, with horizontal struts at the three upper levels.
- X-bracing in every face of every panel.
- A peak joint at (0, 0, 25) m, joined to the 4 top corners.
- A cross-arm reaching (±7, 0, 21) m.
- Joint and member numbering is in the docstring of `make_cases.py`.

| Case | Beam theory | Elements per member (NDiv) | Sections | Base |
|---|---|---|---|---|
| `caseA` | Euler-Bernoulli (FEMMod 1) | 1 | arbitrary (MType 4), spin 0 | 4 legs fixed |
| `caseB` | Timoshenko (FEMMod 3) | 2 | circular legs (1c), rectangular struts and peak members (1r), arbitrary diagonals spun 30 deg, arbitrary cross-arm | 3 legs fixed; joint 1 free on a 6 x 6 spring (`towerB_SSI_joint1.dat`, with a Kxty coupling) |
| `caseC` | as `caseB` | 2 | as `caseB`, plus a 400 kg concentrated mass at each arm tip (centre 1.5 m below, offset across, products of inertia) | as `caseB`, the spring with a 6 x 6 mass (with an Mxty coupling) |
| `caseT` | Euler-Bernoulli (FEMMod 1) | 1 | arbitrary (MType 4) as `caseA`; a 30 m tower, 6 m wide at the base and 1.5 m at the top, 4 panels, a 12 m cross-arm with a centre joint | 4 legs fixed |
| `caseG` | Euler-Bernoulli (FEMMod 1) | 1 | arbitrary (MType 4): equal-leg angles, 150 x 12 mm legs and cross-arm chords, 90 x 7 mm elsewhere; the tower ERF generates (`ERF_LatticeFrame.H`) for 6 panels, crossed bracing, a 12 m cross-arm 1.5 m deep and a 3 m peak (95 joints, 286 members) | 4 legs fixed |

Case T is a 30 m lattice tower with a 12 m cross-arm, the frame a tower type gives as its
`frame_file`, in tower-local axes: origin at the base centre, x along the line, y along the
cross-arm, z up. Its interface joint is the cross-arm's centre (joint 21), the point lines hang
from. `caseT/towerT_members.dat`, written by `make_cases.py`, is its member design file for ERF's
member checks (`ERF_MemberChecks.H`): angles whose areas are within 2.5 % of the arbitrary
sections'. No unit test reads it; it is the `member_file` a run's tower type names for this frame. Case T has no joints where its diagonals cross, so its 9 m X
diagonals are over the bracing slenderness limit and govern on the elastic buckling branch, at a
utilisation above 1: the file tests the checks, not a sound design.

Cases A, B and C have one interface joint, the peak (joint 17), locked to the transition piece, and
the driver's reference point is at that joint. Cases T and G have theirs at the cross-arm's centre,
(0, 0, 30) m.

Case G's `towerG.dat` is not written by `make_cases.py` but by ERF itself: `lattice_frame()` with the
dimensions of `case_g()` in `Tests/Unit/MovingBodies/ERF_GTestMemberChecks.cpp`, then `write_subdyn()`.
SubDyn reading it to the same stiffness and frequencies checks the generator and the writer together.

## The oracles

`tower<case>_kbbt.txt` is the 6 x 6 matrix `KBBt` that SubDyn writes to `tower<case>.SD.sum.yaml`.
It is the stiffness of the whole tower condensed to the peak joint, with the base supports applied.

- DOF order: x, y, z translations, then rotations about x, y, z.
- Units: N/m, N/rad and N m/rad.
- Precision: SubDyn prints 7 significant digits (format ES15.6E2).
- Gravity does not enter `KBBt`: SubDyn's beams have no geometric stiffness.

`tower<case>_frequencies.txt` lists SubDyn's `Full_frequencies` from the same summary, in Hz. These
are all the natural frequencies of the tower with its base reactions applied and its interface joint
free (90 for case A, 432 for cases B and C, 114 for case T, 546 for case G), ascending, again to 7
significant digits.

## Regenerating

Run with OpenFAST 5.0.0 built in double precision:

```bash
python3 make_cases.py
cd caseA && subdyn_driver towerA.dvr && cd ..
cd caseB && subdyn_driver towerB.dvr && cd ..
cd caseC && subdyn_driver towerC.dvr && cd ..
cd caseT && subdyn_driver towerT.dvr && cd ..
cd caseG && subdyn_driver towerG.dvr && cd ..
```

`make_cases.py` rewrites the `.dat`, `.dvr` and SSI files byte for byte (case G: its `.dvr` only;
regenerate `towerG.dat` with ERF's `write_subdyn()` of `lattice_frame(case_g())`). Copy `KBBt` from each
`.SD.sum.yaml` into `tower<case>_kbbt.txt` (cases A, B, T and G), and `Full_frequencies` into
`tower<case>_frequencies.txt`; the `#` header lines of those files are comments.

## Data and acknowledgements

The input files are written by `make_cases.py` (case G by ERF itself); the oracle files
(`tower<case>_kbbt.txt`, `tower<case>_frequencies.txt`) are SubDyn's summary output for them, from
OpenFAST 5.0.0 (NREL, https://github.com/OpenFAST/openfast, Apache License 2.0), built in double
precision. Thanks to the OpenFAST and SubDyn developers at NREL for the module these numbers come
from. Everything here is redistributed in the repository; nothing has to be downloaded to run the
tests.
