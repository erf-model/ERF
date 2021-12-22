#include <ERF.H>
#include <AMReX_buildInfo.H>

extern std::string inputs_name;

void
ERF::writeJobInfo(const std::string& dir) const
{
  // job_info file with details about the run
  std::ofstream jobInfoFile;
  std::string FullPathJobInfoFile = dir;
  FullPathJobInfoFile += "/job_info";
  jobInfoFile.open(FullPathJobInfoFile.c_str(), std::ios::out);

  std::string PrettyLine = "==================================================="
                           "============================\n";
  std::string OtherLine = "----------------------------------------------------"
                          "----------------------------\n";
  std::string SkipSpace = "        ";

  // job information
  jobInfoFile << PrettyLine;
  jobInfoFile << " ERF Job Information\n";
  jobInfoFile << PrettyLine;

  jobInfoFile << "inputs file: " << inputs_name << "\n\n";

  jobInfoFile << "number of MPI processes: "
              << amrex::ParallelDescriptor::NProcs() << "\n";
#ifdef _OPENMP
  jobInfoFile << "number of threads:       " << omp_get_max_threads() << "\n";
#endif

  jobInfoFile << "\n";
  jobInfoFile << "CPU time used since start of simulation (CPU-hours): "
              << getCPUTime() / 3600.0;

  jobInfoFile << "\n\n";

  // plotfile information
  jobInfoFile << PrettyLine;
  jobInfoFile << " Plotfile Information\n";
  jobInfoFile << PrettyLine;

  time_t now = time(0);

  // Convert now to tm struct for local timezone
  tm* localtm = localtime(&now);
  jobInfoFile << "output data / time: " << asctime(localtm);

  std::string currentDir = amrex::FileSystem::CurrentPath();
  jobInfoFile << "output dir:         " << currentDir << "\n";

  jobInfoFile << "\n\n";

  // build information
  jobInfoFile << PrettyLine;
  jobInfoFile << " Build Information\n";
  jobInfoFile << PrettyLine;

  jobInfoFile << "build date:    " << amrex::buildInfoGetBuildDate() << "\n";
  jobInfoFile << "build machine: " << amrex::buildInfoGetBuildMachine() << "\n";
  jobInfoFile << "build dir:     " << amrex::buildInfoGetBuildDir() << "\n";
  jobInfoFile << "AMReX dir:     " << amrex::buildInfoGetAMReXDir() << "\n";

  jobInfoFile << "\n";

  jobInfoFile << "COMP:          " << amrex::buildInfoGetComp() << "\n";
  jobInfoFile << "COMP version:  " << amrex::buildInfoGetCompVersion() << "\n";

  jobInfoFile << "\n";

  for (int n = 1; n <= amrex::buildInfoGetNumModules(); n++) {
    jobInfoFile << amrex::buildInfoGetModuleName(n) << ": "
                << amrex::buildInfoGetModuleVal(n) << "\n";
  }

  jobInfoFile << "\n";

  const char* githash1 = amrex::buildInfoGetGitHash(1);
  const char* githash2 = amrex::buildInfoGetGitHash(2);
  if (strlen(githash1) > 0) {
    jobInfoFile << "ERF       git hash: " << githash1 << "\n";
  }
  if (strlen(githash2) > 0) {
    jobInfoFile << "AMReX       git hash: " << githash2 << "\n";
  }

  const char* buildgithash = amrex::buildInfoGetBuildGitHash();
  const char* buildgitname = amrex::buildInfoGetBuildGitName();
  if (strlen(buildgithash) > 0) {
    jobInfoFile << buildgitname << " git hash: " << buildgithash << "\n";
  }

  jobInfoFile << "\n\n";

  // grid information
  jobInfoFile << PrettyLine;
  jobInfoFile << " Grid Information\n";
  jobInfoFile << PrettyLine;

  int f_lev = finest_level;

  for (int i = 0; i <= f_lev; i++) {
    jobInfoFile << " level: " << i << "\n";
    jobInfoFile << "   number of boxes = " << grids[i].size() << "\n";
    jobInfoFile << "   maximum zones   = ";
    for (int n = 0; n < AMREX_SPACEDIM; n++) {
      jobInfoFile << geom[i].Domain().length(n) << " ";
    }
    jobInfoFile << "\n\n";
  }

  jobInfoFile << " Boundary conditions\n";
  amrex::Vector<std::string> lo_bc_out(AMREX_SPACEDIM);
  amrex::Vector<std::string> hi_bc_out(AMREX_SPACEDIM);
  amrex::ParmParse pp("erf");
  //
  // Check for integer BC type specification in inputs file (older style)
  //
  if ( pp.contains("lo_bc") )
  {
      pp.getarr("lo_bc", lo_bc_out, 0, AMREX_SPACEDIM);
      pp.getarr("hi_bc", hi_bc_out, 0, AMREX_SPACEDIM);
  }

  if (geom[0].isPeriodic(0))
  {
      jobInfoFile << "   -x: " << "Interior" << "\n";
      jobInfoFile << "   +x: " << "Interior" << "\n";
  } else {
      jobInfoFile << "   -x: " << lo_bc_out[0] << "\n";
      jobInfoFile << "   +x: " << hi_bc_out[0] << "\n";
  }
  if (geom[0].isPeriodic(1))
  {
      jobInfoFile << "   -y: " << "Interior" << "\n";
      jobInfoFile << "   +y: " << "Interior" << "\n";
  } else {
      jobInfoFile << "   -y: " << lo_bc_out[1] << "\n";
      jobInfoFile << "   +y: " << hi_bc_out[1] << "\n";
  }
  if (geom[0].isPeriodic(2))
  {
      jobInfoFile << "   -z: " << "Interior" << "\n";
      jobInfoFile << "   +z: " << "Interior" << "\n";
  } else {
      jobInfoFile << "   -z: " << lo_bc_out[2] << "\n";
      jobInfoFile << "   +z: " << hi_bc_out[2] << "\n";
  }

  jobInfoFile << "\n\n";

  // runtime parameters
  jobInfoFile << PrettyLine;
  jobInfoFile << " Inputs File Parameters\n";
  jobInfoFile << PrettyLine;

  amrex::ParmParse::dumpTable(jobInfoFile, true);
  jobInfoFile.close();
}

void
ERF::writeBuildInfo(std::ostream& os)
{
  std::string PrettyLine = std::string(78, '=') + "\n";
  std::string OtherLine = std::string(78, '-') + "\n";
  std::string SkipSpace = std::string(8, ' ');

  // build information
  os << PrettyLine;
  os << " ERF Build Information\n";
  os << PrettyLine;

  os << "build date:    " << amrex::buildInfoGetBuildDate() << "\n";
  os << "build machine: " << amrex::buildInfoGetBuildMachine() << "\n";
  os << "build dir:     " << amrex::buildInfoGetBuildDir() << "\n";
  os << "AMReX dir:     " << amrex::buildInfoGetAMReXDir() << "\n";

  os << "\n";

  os << "COMP:          " << amrex::buildInfoGetComp() << "\n";
  os << "COMP version:  " << amrex::buildInfoGetCompVersion() << "\n";

  os << "C++ compiler:  " << amrex::buildInfoGetCXXName() << "\n";
  os << "C++ flags:     " << amrex::buildInfoGetCXXFlags() << "\n";

  os << "\n";

  os << "Link flags:    " << amrex::buildInfoGetLinkFlags() << "\n";
  os << "Libraries:     " << amrex::buildInfoGetLibraries() << "\n";

  os << "\n";

  for (int n = 1; n <= amrex::buildInfoGetNumModules(); n++) {
    os << amrex::buildInfoGetModuleName(n) << ": "
       << amrex::buildInfoGetModuleVal(n) << "\n";
  }

  os << "\n";
  const char* githash1 = amrex::buildInfoGetGitHash(1);
  const char* githash2 = amrex::buildInfoGetGitHash(2);
  if (strlen(githash1) > 0) {
    os << "ERF       git hash: " << githash1 << "\n";
  }
  if (strlen(githash2) > 0) {
    os << "AMReX       git hash: " << githash2 << "\n";
  }

  const char* buildgithash = amrex::buildInfoGetBuildGitHash();
  const char* buildgitname = amrex::buildInfoGetBuildGitName();
  if (strlen(buildgithash) > 0) {
    os << buildgitname << " git hash: " << buildgithash << "\n";
  }

  os << "\n";
  os << " ERF Compile time variables: \n";

  os << "\n";
  os << " ERF Defines: \n";
#ifdef _OPENMP
  os << std::setw(35) << std::left << "_OPENMP " << std::setw(6) << "ON"
     << std::endl;
#else
  os << std::setw(35) << std::left << "_OPENMP " << std::setw(6) << "OFF"
     << std::endl;
#endif

#ifdef MPI_VERSION
  os << std::setw(35) << std::left << "MPI_VERSION " << std::setw(6)
     << MPI_VERSION << std::endl;
#else
  os << std::setw(35) << std::left << "MPI_VERSION " << std::setw(6)
     << "UNDEFINED" << std::endl;
#endif

#ifdef MPI_SUBVERSION
  os << std::setw(35) << std::left << "MPI_SUBVERSION " << std::setw(6)
     << MPI_SUBVERSION << std::endl;
#else
  os << std::setw(35) << std::left << "MPI_SUBVERSION " << std::setw(6)
     << "UNDEFINED" << std::endl;
#endif

  os << "\n\n";
}
