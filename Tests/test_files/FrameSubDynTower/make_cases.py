#!/usr/bin/env python3
"""Generate SubDyn reference inputs (cases A, B, C and T, and case G's driver file) for small lattice towers.

Joint numbering (1-based):
  level k = 0..3 at z = 7k, half-width h = 3 - 2 z / 21
  corner c = 0..3 : (+h,+h), (-h,+h), (-h,-h), (+h,-h)
  joint id = 4 k + c + 1   (1..16)
  17 = peak (0,0,25), 18 = arm tip (+7,0,21), 19 = arm tip (-7,0,21)
Member numbering:
  1-12  legs        (c major, k minor: id = 3 c + k + 1)
  13-24 struts      (levels 1..3, corner c -> c+1)
  25-48 diagonals   (panel k, face c->c+1: d1 = (k,c)->(k+1,c+1), d2 = (k,c+1)->(k+1,c))
  49-52 peak members (top corner c -> 17)
  53-56 cross-arm   (13->18, 16->18, 14->19, 15->19)
"""
import os
import sys

HERE = os.path.dirname(os.path.abspath(__file__))
SIGNS = [(1, 1), (-1, 1), (-1, -1), (1, -1)]


def jid(k, c):
    return 4 * k + (c % 4) + 1


def joints():
    out = []
    for k in range(4):
        z = 7.0 * k
        h = 3.0 - 2.0 * z / 21.0
        for c in range(4):
            sx, sy = SIGNS[c]
            out.append((jid(k, c), sx * h, sy * h, z))
    out.append((17, 0.0, 0.0, 25.0))
    out.append((18, 7.0, 0.0, 21.0))
    out.append((19, -7.0, 0.0, 21.0))
    return out


def members():
    """Return list of (group, j1, j2) in member-id order."""
    m = []
    for c in range(4):
        for k in range(3):
            m.append(("leg", jid(k, c), jid(k + 1, c)))
    for k in range(1, 4):
        for c in range(4):
            m.append(("strut", jid(k, c), jid(k, c + 1)))
    for k in range(3):
        for c in range(4):
            m.append(("diag", jid(k, c), jid(k + 1, c + 1)))
            m.append(("diag", jid(k, c + 1), jid(k + 1, c)))
    for c in range(4):
        m.append(("peak", jid(3, c), 17))
    m.append(("arm", 13, 18))
    m.append(("arm", 16, 18))
    m.append(("arm", 14, 19))
    m.append(("arm", 15, 19))
    return m


def joints_t():
    """Case T: a 30 m lattice tower with a 12 m cross-arm, in tower-local axes (x along the
    line, y along the cross-arm, z up, origin at the base centre): 5 levels at z = 7.5 k, half-width
    3 - 2.25 z / 30; 21 = the cross-arm's centre (0, 0, 30) the line hangs from; 22, 23 = the arm tips
    (0, +-6, 30)."""
    out = []
    for k in range(5):
        z = 7.5 * k
        h = 3.0 - 2.25 * z / 30.0
        for c in range(4):
            sx, sy = SIGNS[c]
            out.append((jid(k, c), sx * h, sy * h, z))
    out.append((21, 0.0, 0.0, 30.0))
    out.append((22, 0.0, 6.0, 30.0))
    out.append((23, 0.0, -6.0, 30.0))
    return out


def members_t():
    m = []
    for c in range(4):
        for k in range(4):
            m.append(("leg", jid(k, c), jid(k + 1, c)))
    for k in range(1, 5):
        for c in range(4):
            m.append(("strut", jid(k, c), jid(k, c + 1)))
    for k in range(4):
        for c in range(4):
            m.append(("diag", jid(k, c), jid(k + 1, c + 1)))
            m.append(("diag", jid(k, c + 1), jid(k + 1, c)))
    for c in range(4):
        m.append(("centre", jid(4, c), 21))
    # each arm tip to the two top corners and the two corners one level down on its side
    for tip, corners in ((22, (0, 1)), (23, (2, 3))):
        for k in (4, 3):
            for c in corners:
                m.append(("arm", jid(k, c), tip))
    return m


STEEL = "2.0E+11  7.7E+10  7850.0"


def write_members_t(path):
    """Case T's member design file (ERF_MemberChecks.H): equal-leg angles whose areas are within 2.5 % of
    the arbitrary sections' (legs 4.0e-3 m^2: 150 x 14 mm; the rest 1.0e-3 m^2: 65 x 8 mm), 345 MPa
    steel; the legs bolted in both faces, the rest by one leg with framing eccentricity at both ends,
    the struts redundant."""
    L = ["# Member design data of case T (towerT.dat, the 30 m lattice tower), for ERF's member checks",
         "# MemberID  Role  Fy(Pa)  b(m)  t(m)  NetArea(-)  Bolted  Ends  Restraint"]
    for n, (g, _, _) in enumerate(members_t(), start=1):
        if g == "leg":
            L.append("%d leg 3.45e8 0.15 0.014 1.0 both concentric none" % n)
        else:
            role = "redundant" if g == "strut" else "bracing"
            L.append("%d %s 3.45e8 0.065 0.008 0.85 one both none" % (n, role))
    with open(path, "w") as f:
        f.write("\n".join(L) + "\n")


def write_dat(case, path, ssi_name=None):
    A = case in ("A", "T")
    femmod = 1 if A else 3
    ndiv = 1 if A else 2
    J = joints_t() if case == "T" else joints()
    M = members_t() if case == "T" else members()
    interface = 21 if case == "T" else 17
    L = []
    w = L.append
    w("----------- SubDyn MultiMember Support Structure Input File ---------------------------")
    w(("Lattice tower case T for ERF's frame model, tower-local axes (4 tapered legs, struts, X-bracing, cross-arm)."
       if case == "T" else
       "Lattice tower oracle case %s for the ERF frame solver (4 tapered legs, struts, X-bracing, peak, cross-arm)." % case))
    w("-------------------------- SIMULATION CONTROL -----------------------------------------")
    w("%-16s Echo        - Echo input data to \"<rootname>.SD.ech\" (flag)" % ("False" if A else "True"))
    w("\"DEFAULT\"        SDdeltaT    - Local Integration Step. If \"default\", the glue-code integration step will be used.")
    w("             3   IntMethod   - Integration Method [1/2/3/4 = RK4/AB4/ABM4/AM2].")
    w("False            SttcSolve   - Solve dynamics about static equilibrium point")
    w("-------------------- FEA and CRAIG-BAMPTON PARAMETERS ---------------------------------")
    w("             %d   FEMMod      - FEM switch: element model in the FEM. [1= Euler-Bernoulli(E-B);  2=Tapered E-B (unavailable);  3= 2-node Timoshenko;  4= 2-node tapered Timoshenko (unavailable)]" % femmod)
    w("             %d   NDiv        - Number of sub-elements per member" % ndiv)
    w("             0   Nmodes      - Number of internal modes to retain. If Nmodes=0 --> Guyan Reduction. If Nmodes<0 --> retain all modes.")
    w("             0   JDampings   - Damping Ratios for each retained mode (% of critical) If Nmodes>0, list Nmodes structural damping ratios for each retained mode (% of critical), or a single damping ratio to be applied to all retained modes. (last entered value will be used for all remaining modes).")
    w("             0   GuyanDampMod - Guyan damping {0=none, 1=Rayleigh Damping, 2=user specified 6x6 matrix}")
    w("  0.000, 0.000   RayleighDamp - Mass and stiffness proportional damping coefficients (Rayleigh Damping) [only if GuyanDampMod=1]")
    w("             6   GuyanDampSize - Guyan damping matrix (6nTPx6nTP if fixed bottom or 6(nTP-1)-by-6(nTP-1) if floating) [only if GuyanDampMod=2]")
    for _ in range(6):
        w("   0.0000e+00   0.0000e+00   0.0000e+00   0.0000e+00   0.0000e+00   0.0000e+00")
    w("------- INITIAL RIGID-BODY POSITION [used only for floating structure with more than one transition pieces] -------")
    w("RBSurge    RBSway     RBHeave    RBRoll     RBPitch    RBYaw")
    w("  (m)        (m)        (m)      (deg)      (deg)      (deg)")
    w("  0.0        0.0        0.0       0.0        0.0        0.0")
    w("---- STRUCTURE JOINTS: joints connect structure members (~Hydrodyn Input File) --------")
    w("            %2d   NJoints     - Number of joints (-)" % len(J))
    w("JointID          JointXss               JointYss               JointZss     JointType JointDirX  JointDirY JointDirZ JointStiff    ![Coordinates of Member joints in SS-Coordinate System][JointType={1:cantilever, 2:universal joint, 3:revolute joint, 4:spherical joint}]")
    w("  (-)               (m)                    (m)                    (m)         (-)        (-)        (-)       (-)     (Nm/rad) ")
    for (i, x, y, z) in J:
        w("  %3d   %22.15E %22.15E %22.15E   1   0.0   0.0   0.0   0.0" % (i, x, y, z))
    w("------------------- BASE REACTION JOINTS: 1/0 for Locked/Free DOF @ each Reaction Node ---------------------")
    w("             4   NReact      - Number of Joints with reaction forces; be sure to remove all rigid motion DOFs of the structure  (else det([K])=[0])")
    w("RJointID   RctTDXss    RctTDYss    RctTDZss    RctRDXss    RctRDYss    RctRDZss     SSIfile ![Global Coordinate System]")
    w("  (-)       (flag)      (flag)      (flag)      (flag)      (flag)      (flag)      (string)      ")
    for j in (1, 2, 3, 4):
        if (not A) and j == 1:
            w("   %d          0           0           0           0           0           0        \"%s\"" % (j, ssi_name))
        else:
            w("   %d          1           1           1           1           1           1" % j)
    w("------- INTERFACE JOINTS: 1/0 for Locked (to the TP)/Free DOF @each Interface Joint (only Locked-to-TP implemented thus far (=rigid TP)) ---------")
    w("             1   NInterf     - Number of interface joints locked to the Transition Piece (TP):  be sure to remove all rigid motion dofs")
    w("IJointID   TPID   ItfTDXss    ItfTDYss    ItfTDZss    ItfRDXss    ItfRDYss    ItfRDZss     ![Global Coordinate System]")
    w("  (-)      (-)     (flag)      (flag)      (flag)      (flag)      (flag)      (flag)")
    w("  %2d        1        1           1           1           1           1           1" % interface)
    w("----------------------------------- MEMBERS -------------------------------------------")
    w("            %2d   NMembers    - Number of members (-)" % len(M))
    w("MemberID   MJointID1   MJointID2   MPropSetID1   MPropSetID2  MType  COSMID/MSpin   ![MType={1:beam circ., 2:cable, 3:rigid, 4:beam arb., 5:spring}. COMSID={-1:none}]")
    w("  (-)         (-)         (-)          (-)           (-)        (-)    (-)/(deg)")
    for n, (g, j1, j2) in enumerate(M, start=1):
        if A:
            pid, mt, spin = (1, "4", 0) if g == "leg" else (2, "4", 0)
        else:
            if g == "leg":
                pid, mt, spin = 1, "1c", 0
            elif g == "diag":
                pid, mt, spin = 3, "4", 30
            elif g in ("strut", "peak"):
                pid, mt, spin = 2, "1r", 0
            else:  # arm
                pid, mt, spin = 4, "4", 0
        w("  %3d        %3d         %3d          %3d           %3d        %-3s    %4g          # %s" % (n, j1, j2, pid, pid, mt, spin, g))
    w("------------------ CIRCULAR BEAM CROSS-SECTION PROPERTIES -----------------------------")
    w("             %d   NPropSetsCyl - Number of structurally unique circular cross-sections (if 0 the following table is ignored)" % (0 if A else 1))
    w("PropSetID     YoungE          ShearG          MatDens          XsecD           XsecT")
    w("  (-)         (N/m2)          (N/m2)          (kg/m3)           (m)             (m)  ")
    if not A:
        w("   1        %s    0.3     0.012      # legs" % STEEL)
    w("----------------- RECTANGULAR BEAM CROSS-SECTION PROPERTIES ---------------------------")
    w("             %d   NPropSetsRec - Number of structurally unique rectangular cross-sections (if 0 the following table is ignored)" % (0 if A else 1))
    w("PropSetID     YoungE          ShearG          MatDens          XsecSa         XsecSb          XsecT")
    w("  (-)         (N/m2)          (N/m2)          (kg/m3)           (m)            (m)             (m)")
    if not A:
        w("   2        %s    0.15    0.10    0.008      # struts and peak members" % STEEL)
    w("----------------- ARBITRARY BEAM CROSS-SECTION PROPERTIES -----------------------------")
    w("             %d   NXPropSets   - Number of structurally unique arbitrary cross-sections (if 0 the following table is ignored)" % 2)
    w("PropSetID     YoungE          ShearG          MatDens          XsecA          XsecAsx       XsecAsy       XsecJxx       XsecJyy        XsecJ0    XsecJt")
    w("  (-)         (N/m2)          (N/m2)          (kg/m3)          (m2)            (m2)          (m2)          (m4)          (m4)          (m4)       (m4)")
    legs = "4.0E-03  2.0E-03  2.0E-03  1.2E-05  0.6E-05  1.8E-05  4.0E-07"
    if A:
        w("   1        %s    %s   # legs" % (STEEL, legs))
        w("   2        %s    1.0E-03  5.0E-04  5.0E-04  8.0E-07  5.0E-07  1.3E-06  3.0E-09   # struts, diagonals, peak, arm" % STEEL)
    else:
        w("   3        %s    1.0E-03  3.0E-04  4.0E-04  8.0E-07  5.0E-07  1.3E-06  3.0E-09   # diagonals (MSpin 30 deg)" % STEEL)
        w("   4        %s    %s   # cross-arm (legs' set)" % (STEEL, legs))
    w("-------------------------- CABLE PROPERTIES -------------------------------------------")
    w("             0   NCablePropSets   - Number of cable cable properties")
    w("PropSetID     EA          MatDens        T0         CtrlChannel")
    w("  (-)         (N)         (kg/m)        (N)             (-)")
    w("----------------------- RIGID LINK PROPERTIES -----------------------------------------")
    w("             0   NRigidPropSets - Number of rigid link properties")
    w("PropSetID   MatDens   ")
    w("  (-)       (kg/m)")
    w("----------------------- SPRING ELEMENT PROPERTIES -------------------------------------")
    w("             0   NSpringPropSets - Number of spring properties")
    w("PropSetID   k11     k12     k13     k14     k15     k16     k22     k23     k24     k25     k26     k33     k34     k35     k36     k44      k45      k46      k55      k56      k66    ")
    w("  (-)      (N/m)   (N/m)   (N/m)  (N/rad) (N/rad) (N/rad)  (N/m)   (N/m)  (N/rad) (N/rad) (N/rad)  (N/m)  (N/rad) (N/rad) (N/rad) (Nm/rad) (Nm/rad) (Nm/rad) (Nm/rad) (Nm/rad) (Nm/rad)          ")
    w("---------------------- MEMBER COSINE MATRICES COSM(i,j) -------------------------------")
    w("             0   NCOSMs      - Number of unique cosine matrices (i.e., of unique member alignments including principal axis rotations); ignored if NXPropSets=0   or 9999 in any element below")
    w("COSMID    COSM11    COSM12    COSM13    COSM21    COSM22    COSM23    COSM31    COSM32    COSM33")
    w(" (-)       (-)       (-)       (-)       (-)       (-)       (-)       (-)       (-)       (-)     ")
    w("------------------------ JOINT ADDITIONAL CONCENTRATED MASSES--------------------------")
    cmass = CMASSES if case == "C" else []
    w("%14d   NCmass      - Number of joints with concentrated masses; Global Coordinate System" % len(cmass))
    w("CMJointID       JMass            JMXX             JMYY             JMZZ          JMXY        JMXZ         JMYZ        MCGX      MCGY        MCGZ")
    w("  (-)            (kg)          (kg*m^2)         (kg*m^2)         (kg*m^2)      (kg*m^2)    (kg*m^2)     (kg*m^2)       (m)      (m)          (m)")
    for row in cmass:
        w("   %3d   " % row[0] + "  ".join("%14.6E" % v for v in row[1:]))
    w("---------------------------- OUTPUT: SUMMARY & OUTFILE --------------------------------")
    w("True             SumPrint    - Output a Summary File (flag)")
    w("0                OutCBModes  - Output Guyan and Craig-Bampton modes {0: No output, 1: JSON output}, (flag)")
    w("0                OutFEMModes - Output first 30 FEM modes {0: No output, 1: JSON output} (flag)")
    w("False            OutCOSM     - Output cosine matrices with the selected output member forces (flag)")
    w("False            OutAll      - [T/F] Output all members' end forces")
    w("             1   OutSwtch    - [1/2/3] Output requested channels to: 1=<rootname>.SD.out;  2=<rootname>.out (generated by FAST);  3=both files.")
    w("True             TabDelim    - Generate a tab-delimited output in the <rootname>.SD.out file")
    w("             1   OutDec      - Decimation of output in the <rootname>.SD.out file")
    w("\"ES11.4e2\"       OutFmt      - Output format for numerical results in the <rootname>.SD.out file")
    w("\"A11\"            OutSFmt     - Output format for header strings in the <rootname>.SD.out file")
    w("------------------------- MEMBER OUTPUT LIST ------------------------------------------")
    w("             0   NMOutputs   - Number of members whose forces/displacements/velocities/accelerations will be output (-) [Must be <= 99].")
    w("MemberID   NOutCnt    NodeCnt ![NOutCnt=how many nodes to get output for [< 10]; NodeCnt are local ordinal numbers from the start of the member, and must be >=1 and <= NDiv+1] If NMOutputs=0 leave blank as well.")
    w("  (-)        (-)        (-)")
    w("------------------------- SSOutList: The next line(s) contains a list of output parameters that will be output in <rootname>.SD.out or <rootname>.out. ------")
    w("END of output channels and end of file. (the word \"END\" must appear in the first 3 columns of this line)")
    with open(path, "w") as f:
        f.write("\n".join(x.rstrip() for x in L) + "\n")


def write_dvr(case, path, datname, root):
    L = []
    w = L.append
    zref = 30 if case in ("T", "G") else 25
    w("SubDyn Driver file for stand-alone applications")
    w(("Lattice tower case %s: Guyan KBBt at the cross-arm centre (TP ref point 0,0,30)." % case if case in ("T", "G") else
       "Lattice tower oracle case %s: Guyan KBBt at the peak joint (TP ref point 0,0,25)." % case))
    w("False               Echo           - Echo the input file data (flag)")
    w("---------------------- ENVIRONMENTAL CONDITIONS -------------------------------------------------")
    w("9.80665             Gravity        - Gravity (m/s^2).")
    w("0                   WtrDpth        - Water Depth (m) positive value.")
    w("---------------------- SubDyn -------------------------------------------------------------------")
    w("\"%s\" SDInputFile    - Absolute or relative path." % datname)
    w("\"%s\"     OutRootName    - Basename for output files." % root)
    w("2                                   NSteps         - Number of time steps in the simulations (-)")
    w("0.01                                TimeInterval   - TimeInterval for the simulation (sec)")
    w("1                                   NTPs           - Number of transition pieces")
    w("0                                   TP_RefPoint_X  - X location of the TP reference points in global coordinates (m) {require NTPs entries}")
    w("0                                   TP_RefPoint_Y  - Y location of the TP reference points in global coordinates (m) {require NTPs entries}")
    w("%-36dTP_RefPoint_Z  - Z location of the TP reference points in global coordinates (m) {require NTPs entries}" % zref)
    w("0                                   SubRotateZ     - Rotation angle of the structure geometry in [deg] about the global Z axis.")
    w("---------------------- INPUTS -------------------------------------------------------------------")
    w("0                   InputsMod      - Inputs model {0: all inputs are zero for every timestep, 1: steady state inputs, 2: read inputs from a file (InputsFile)} (switch)")
    w("\"unused\"            InputsFile     - Name of the inputs file if InputsMod = 2.")
    w("---------------------- STEADY INPUTS (for InputsMod = 1) ----------------------------------------")
    w("0.0  0.0  0.0  0.0  0.0  0.0   uTPInSteady       - input displacements and rotations ( m, rads )")
    w("0.0  0.0  0.0  0.0  0.0  0.0   uDotTPInSteady    - input translational and rotational velocities ( m/s, rads/s)")
    w("0.0  0.0  0.0  0.0  0.0  0.0   uDotDotTPInSteady - input translational and rotational accelerations( m/s^2, rads/s^2)")
    w("---------------------- LOADS --------------------------------------------------------------------")
    w("0    nAppliedLoads  - Number of applied loads at given nodes")
    w("ALJointID    Fx     Fy    Fz     Mx     My     Mz   UnsteadyFile")
    w("   (-)       (N)    (N)   (N)   (Nm)   (Nm)   (Nm)     (-)")
    w("END of driver input file")
    with open(path, "w") as f:
        f.write("\n".join(x.rstrip() for x in L) + "\n")


# case C: case B plus the masses an insulator string and a line clamp put on each arm tip, hanging
# 1.5 m below it (mass, Jxx Jyy Jzz Jxy Jxz Jyz about its centre, centre offset x y z), and a
# foundation mass on joint 1's spring
CMASSES = [(18, 400.0, 20.0, 30.0, 25.0, 2.0, -1.0, 1.5, 0.1, -0.2, -1.5),
           (19, 400.0, 20.0, 30.0, 25.0, 2.0, -1.0, 1.5, -0.1, 0.2, -1.5)]
SSI_MASS = {"Mxx": 2.0e3, "Myy": 2.0e3, "Mzz": 2.0e3, "Mtxtx": 5.0e2, "Mtyty": 5.0e2, "Mtztz": 4.0e2, "Mxty": 1.0e2}


def write_ssi(path, with_mass=False):
    K = {"Kxx": 2.0e8, "Kyy": 2.0e8, "Kzz": 5.0e8,
         "Ktxtx": 3.0e8, "Ktyty": 3.0e8, "Ktztz": 1.0e8, "Kxty": 1.0e7}
    names = ['Kxx', 'Kxy', 'Kyy', 'Kxz', 'Kyz', 'Kzz', 'Kxtx', 'Kytx', 'Kztx', 'Ktxtx',
             'Kxty', 'Kyty', 'Kzty', 'Ktxty', 'Ktyty', 'Kxtz', 'Kytz', 'Kztz', 'Ktxtz', 'Ktytz', 'Ktztz']
    L = ["!---------------- SSI spring at joint 1 (+3,+3,0): K entries, M %s -------------------!" % ("entries" if with_mass else "all zero"),
         "!Upper-triangular names, value first then name; Kxty couples x-translation with y-rotation (K(1,5))"]
    for n in names:
        L.append("   %.6E        %s" % (K.get(n, 0.0), n))
    for n in names:
        m = "M" + n[1:]
        L.append("   %.6E        %s" % (SSI_MASS.get(m, 0.0) if with_mass else 0.0, m))
    with open(path, "w") as f:
        f.write("\n".join(x.rstrip() for x in L) + "\n")


if __name__ == "__main__":
    for case in ("A", "B", "C", "T"):
        d = os.path.join(HERE, "case" + case)
        os.makedirs(d, exist_ok=True)
        root = "tower" + case
        ssi = "tower%s_SSI_joint1.dat" % case if case in ("B", "C") else None
        write_dat(case, os.path.join(d, root + ".dat"), ssi)
        write_dvr(case, os.path.join(d, root + ".dvr"), root + ".dat", root)
        if ssi:
            write_ssi(os.path.join(d, ssi), with_mass=(case == "C"))
        if case == "T":
            write_members_t(os.path.join(d, "towerT_members.dat"))
    # case G's .dat is written by ERF itself (write_subdyn of the generated tower); only its driver file is here
    os.makedirs(os.path.join(HERE, "caseG"), exist_ok=True)
    write_dvr("G", os.path.join(HERE, "caseG", "towerG.dvr"), "towerG.dat", "towerG")
    print("joints", len(joints()), "members", len(members()))
