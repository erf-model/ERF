
namespace amrex {

const char* buildInfoGetBuildDate() {

  static const char BUILD_DATE[] = "2022-04-06 10:55:49.632412";
  return BUILD_DATE;
}

const char* buildInfoGetBuildDir() {

  static const char BUILD_DIR[] = "/home/alattanz/git/ERF/Exec/CouetteFlow";
  return BUILD_DIR;
}

const char* buildInfoGetBuildMachine() {

  static const char BUILD_MACHINE[] = "Linux LAPTOP-PO22H9R0 4.4.0-22000-Microsoft #1-Microsoft Fri Jun 04 16:28:00 PST 2021 x86_64 x86_64 x86_64 GNU/Linux";
  return BUILD_MACHINE;
}

const char* buildInfoGetAMReXDir() {

  static const char AMREX_DIR[] = "/home/alattanz/git/amrex/";
  return AMREX_DIR;
}

const char* buildInfoGetComp() {

  static const char COMP[] = "gnu";
  return COMP;
}

const char* buildInfoGetCompVersion() {

  static const char COMP_VERSION[] = "9.4.0";
  return COMP_VERSION;
}

// deprecated
const char* buildInfoGetFcomp() {

  static const char FCOMP[] = "";
  return FCOMP;
}

// deprecated
const char* buildInfoGetFcompVersion() {

  static const char FCOMP_VERSION[] = "";
  return FCOMP_VERSION;
}

const char* buildInfoGetCXXName() {

  static const char CXX_comp_name[] = "mpicxx";
  return CXX_comp_name;
}

const char* buildInfoGetFName() {

  static const char F_comp_name[] = "";
  return F_comp_name;
}

const char* buildInfoGetCXXFlags() {

  static const char CXX_flags[] = " -Werror=return-type -g -O3  -pthread   -Werror=return-type -g -O3  -pthread   -DBL_USE_MPI -DAMREX_USE_MPI -DBL_SPACEDIM=3 -DAMREX_SPACEDIM=3 -DBL_FORT_USE_UNDERSCORE -DAMREX_FORT_USE_UNDERSCORE -DBL_Linux -DAMREX_Linux -DNDEBUG -DOMPI_SKIP_MPICXX -Itmp_build_dir/s/3d.gnu.MPI.EXE -I. -I../../Source -I../../Source/SpatialStencils -I../../Source/TimeIntegration -I../../Source/IO -I../../Exec/CouetteFlow -I/home/alattanz/git/amrex//Src/Base -I/home/alattanz/git/amrex//Src/Base/Parser -I/home/alattanz/git/amrex//Src/Base -I/home/alattanz/git/amrex//Src/Base/Parser -I/home/alattanz/git/amrex//Src/Boundary -I/home/alattanz/git/amrex//Src/AmrCore -I/home/alattanz/git/amrex//Tools/C_scripts ";
  return CXX_flags;
}

const char* buildInfoGetFFlags() {

  static const char F_flags[] = "";
  return F_flags;
}

const char* buildInfoGetLinkFlags() {

  static const char link_flags[] = " -Werror=return-type -g -O3  -pthread   -DBL_USE_MPI -DAMREX_USE_MPI -DBL_SPACEDIM=3 -DAMREX_SPACEDIM=3 -DBL_FORT_USE_UNDERSCORE -DAMREX_FORT_USE_UNDERSCORE -DBL_Linux -DAMREX_Linux -DNDEBUG -DOMPI_SKIP_MPICXX -Itmp_build_dir/s/3d.gnu.MPI.EXE -I. -I../../Source -I../../Source/SpatialStencils -I../../Source/TimeIntegration -I../../Source/IO -I../../Exec/CouetteFlow -I/home/alattanz/git/amrex//Src/Base -I/home/alattanz/git/amrex//Src/Base/Parser -I/home/alattanz/git/amrex//Src/Base -I/home/alattanz/git/amrex//Src/Base/Parser -I/home/alattanz/git/amrex//Src/Boundary -I/home/alattanz/git/amrex//Src/AmrCore -I/home/alattanz/git/amrex//Tools/C_scripts  -L. -L/usr/lib/gcc/x86_64-linux-gnu/9/";
  return link_flags;
}

const char* buildInfoGetLibraries() {

  static const char libraries[] = " -pthread -I/usr/lib/x86_64-linux-gnu/openmpi/lib -L/usr/lib/x86_64-linux-gnu/openmpi/lib -lmpi_usempif08 -lmpi_usempi_ignore_tkr -lmpi_mpifh -lmpi -lmpi_cxx -lmpi -lgfortran -lquadmath";
  return libraries;
}

const char* buildInfoGetAux(int i) {

  //static const char AUX1[] = "${AUX[1]}";

  static const char EMPT[] = "";

  switch(i)
  {

    default: return EMPT;
  }
}

int buildInfoGetNumModules() {
  // int num_modules = X;
  int num_modules = 0;

  return num_modules;
}

const char* buildInfoGetModuleName(int i) {

  //static const char MNAME1[] = "${MNAME[1]}";

  static const char EMPT[] = "";

  switch(i)
  {

    default: return EMPT;
  }
}

const char* buildInfoGetModuleVal(int i) {

  //static const char MVAL1[] = "${MVAL[1]}";

  static const char EMPT[] = "";

  switch(i)
  {

    default: return EMPT;
  }
}

const char* buildInfoGetGitHash(int i) {

  //static const char HASH1[] = "${GIT[1]}";
  static const char HASH1[] = "22.03-15-g04621415c";

  static const char EMPT[] = "";

  switch(i)
  {
    case 1: return HASH1;

    default: return EMPT;
  }
}

const char* buildInfoGetBuildGitHash() {

  //static const char HASH[] = "${GIT}";
  static const char HASH[] = "";


  return HASH;
}

const char* buildInfoGetBuildGitName() {

  //static const char NAME[] = "";
  static const char NAME[] = "";


  return NAME;
}

#ifdef AMREX_USE_CUDA
const char* buildInfoGetCUDAVersion() {

  static const char CUDA_VERSION[] = "";
  return CUDA_VERSION;
}
#endif

}
