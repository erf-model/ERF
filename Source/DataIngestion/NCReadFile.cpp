#include <iomanip>
#include <iostream>
#include <string>
#include <ctime>

#ifdef _OPENMP
#include <omp.h>
#endif

#include <AMReX_Utility.H>
#include <AMReX_buildInfo.H>
#include <AMReX_ParmParse.H>

#include "ERF.H"
#include "NCInterface.H"
#include "NCReadFile.H"
#include "IndexDefines.H"

void
ERF::ncrestart(istream& is)
{
    // trying to read from checkpoint; if nonexisting, set it to 0.
    if (input_version == -1) {
      if (amrex::ParallelDescriptor::IOProcessor()) {
        std::ifstream ERFHeaderFile;
        std::string FullPathERFHeaderFile = papa.theRestartFile();
        FullPathERFHeaderFile += "/ERFHeader";
        ERFHeaderFile.open(FullPathERFHeaderFile.c_str(), std::ios::in);
        if (ERFHeaderFile.good()) {
          char foo[256];
          ERFHeaderFile.getline(foo, 256, ':');
          ERFHeaderFile >> input_version;
          ERFHeaderFile.close();
        } else {
          input_version = 0;
        }
      }
      amrex::ParallelDescriptor::Bcast(
        &input_version, 1, amrex::ParallelDescriptor::IOProcessorNumber());
    }

    AMREX_ASSERT(input_version >= 0);

    erf.parent = &papa;

    // This is urgly in NetCDF input/output, we should store everything in NetCDF file, however,
    // the directory layout is the way need to know the level before it can read the file.
    // TODO: we need find a way to reorganize the output directory structure later.
    is >> erf.level;

    std::string LevelDir, FullPath;
    LevelDir = amrex::Concatenate("Level_", erf.level, 1);
    erf.LevelDirectoryNames(papa.theRestartFile(), LevelDir, FullPath);
    if( ! erf.levelDirectoryCreated) {
      amrex::Print() << "ERF::ncrestart: NetCDF restart file doesn't exist!" << '\n';
    }

    auto ncf = ncutils::NCFile::open(FullPath+"/SD_New_MF.nc");

    std::vector<int> data;
    ncf.get_attr("Level", data);
    erf.level = data[0];

    ncf.get_attr("finest_level",data);
    papa.SetFinestLevel(data[0]);
    int nflev = data[0];

    ncf.get_attr("CoordSys_ID", data);
    int coord = data[0];
    erf.geom.SetCoord(static_cast<amrex::CoordSys::CoordType>(data[0]));

    std::vector<amrex::Real> offset(AMREX_SPACEDIM);
    std::vector<amrex::Real> cellsize(AMREX_SPACEDIM);
    ncf.var("Coord.Offset").get(offset.data(), {0}, {AMREX_SPACEDIM});
    ncf.var("Coord.CellSize").get(cellsize.data(), {0}, {AMREX_SPACEDIM});

    erf.geom.SetOffset(offset.data());
    // TODO (this is important)
    // we need to setup the cellsize to coordinate here!!!
//    geom.SetCellsize(cellsize.data());

    amrex::IntVect lo(AMREX_SPACEDIM);
    amrex::IntVect hi(AMREX_SPACEDIM);
    amrex::IntVect type(AMREX_SPACEDIM);
    ncf.var("Domain.SmallEnd").get(lo.begin(), {0}, {AMREX_SPACEDIM});
    ncf.var("Domain.BigEnd").get(hi.begin(), {0}, {AMREX_SPACEDIM});
    ncf.var("Domain.BType").get(type.begin(), {0}, {AMREX_SPACEDIM});

    amrex::Box bx = amrex::Box(lo,hi,type);

    amrex::Real plo[AMREX_SPACEDIM];
    amrex::Real phi[AMREX_SPACEDIM];
    ncf.var("PDomain.Lo").get(plo, {0}, {AMREX_SPACEDIM});
    ncf.var("PDomain.Hi").get(phi, {0}, {AMREX_SPACEDIM});

    amrex::RealBox rb = amrex::RealBox(plo, phi);

    erf.geom.Domain(bx);
    erf.geom.ProbDomain(rb);

    amrex::IntArray periodic;
    ncf.var("Geom.Periodic").get(periodic.begin(), {0}, {AMREX_SPACEDIM});

    erf.geom.setPeriodicity({{AMREX_D_DECL(periodic[0],periodic[1],periodic[2])}});

    // setup geom
    erf.geom.define(bx, rb, coord, periodic);

    erf.fine_ratio = amrex::IntVect::TheUnitVector(); erf.fine_ratio.scale(-1);
    erf.crse_ratio = amrex::IntVect::TheUnitVector(); erf.crse_ratio.scale(-1);

    if (erf.level > 0) {
       erf.crse_ratio = papa.refRatio(erf.level-1);
    }

    if (erf.level < papa.maxLevel()) {
       erf.fine_ratio = papa.refRatio(erf.level);
    }

    auto grids_dim = ncf.dim("num_grids");
    int num_grids = static_cast<int>(grids_dim.len());
    erf.grids.resize(num_grids);

    for (auto ib=0; ib<num_grids; ++ib) {
       int blo[AMREX_SPACEDIM];
       int bhi[AMREX_SPACEDIM];
       int btype[AMREX_SPACEDIM];

       ncf.var("Grids.SmallEnd").get(blo, {ib, 0}, {1, AMREX_SPACEDIM});
       ncf.var("Grids.BigEnd").get(bhi, {ib, 0}, {1, AMREX_SPACEDIM});
       ncf.var("Grids.BType").get(btype, {ib, 0}, {1, AMREX_SPACEDIM});

       amrex::Box b = amrex::Box(amrex::IntVect(blo), amrex::IntVect(bhi), amrex::IntVect(btype));
       erf.grids.set(ib, b);
    }

    int ndesc = erf.desc_lst.size();
    int nstate = ndesc;

    amrex::Vector<int> state_in_checkpoint(ndesc, 1);
    if (ndesc > nstate) {
        erf.set_state_in_checkpoint(state_in_checkpoint);
    } else {
        BL_ASSERT(nstate == ndesc);
    }

    erf.dmap.define(erf.grids);
    papa.SetBoxArray(erf.level, erf.grids);
    papa.SetDistributionMap(erf.level, erf.dmap);
    erf.m_factory.reset(new amrex::FArrayBoxFactory());

    erf.state.resize(ndesc);

    amrex::Vector<int> npts;
    amrex::Vector<int> nsize;
    amrex::Vector<int> nblocks;
    amrex::Vector<std::string> lo_names;
    amrex::Vector<std::string> hi_names;
    amrex::Vector<std::string> typ_names;
    amrex::Vector<std::string> var_names;
    std::map<std::pair<int, int>, int> var_map;
    int inc = 0;
    for (int typ = 0; typ < erf.desc_lst.size(); typ++) {
        for (int comp = 0; comp < erf.desc_lst[typ].nComp(); comp++) {
           std::string name;
           std::pair<int, int> p {std::pair<int, int>(typ, comp)};
           var_map[p] = inc;
           if (erf.desc_lst[typ].name(comp) != "") {
              name = erf.desc_lst[typ].name(comp);
           } else {
              name = EnumToString(static_cast<StateType>(typ));
           }
             var_names.push_back(name.c_str());
             ++inc;
       }
    }

    for (int k = 0; k < var_names.size(); ++k) {
       int num_pts = static_cast<int>(ncf.dim("num_points_"+std::to_string(k)).len());
       int num_size = static_cast<int>(ncf.dim("nsize_"+std::to_string(k)).len());
       npts.push_back(num_pts);
       nsize.push_back(num_size);
    }

    for (int k = 0; k < erf.desc_lst.size(); ++k) {
      std::string typ_name = EnumToString(static_cast<StateType>(k));
      lo_names.push_back("lo_"+typ_name);
      hi_names.push_back("hi_"+typ_name);
      typ_names.push_back("typ_"+typ_name);
      int num_box = static_cast<int>(ncf.dim("num_boxs_"+std::to_string(k)).len());
      nblocks.push_back(num_box);
    }

    for (int i = 0; i < ndesc; ++i) {
        if (state_in_checkpoint[i]) {

           amrex::Real time, dt_new, dt_old;
           ncf.get_attr("cumtime",data);
           time = static_cast<amrex::Real>(data[0]);
           ncf.get_attr("DtNew",data);
           dt_new = static_cast<amrex::Real>(data[0]);
           ncf.get_attr("DtOld",data);
           dt_old = static_cast<amrex::Real>(data[0]);

           auto state_domain = erf.geom.Domain();
           auto state_grids  = erf.grids;

           amrex::IndexType typ(erf.desc_lst[i].getType());
           if (!typ.cellCentered()) {
             state_domain.convert(typ);
             state_grids.convert(typ);
           }

           amrex::Arena* arena = nullptr;

           // initialize state
           erf.state[i].define(state_domain, state_grids, erf.dmap, erf.desc_lst[i], time, dt_new, *(erf.m_factory));

           amrex::MultiFab& new_data = erf.state[i].newData();

           // initialize multifab
           new_data.define(state_grids,erf.dmap,erf.desc_lst[i].nComp(),erf.desc_lst[i].nExtra(),
                            amrex::MFInfo().SetTag("StateData").SetArena(arena),
                            *(erf.m_factory));

           new_data.setVal(0.);

           erf.state[i].setTimeLevel(time, dt_old, dt_new);

          for (amrex::MFIter new_mfi(new_data); new_mfi.isValid(); ++new_mfi) {
              auto ncomp   = new_data.nComp();

              amrex::IntVect smallend(AMREX_SPACEDIM);
              amrex::IntVect bigend(AMREX_SPACEDIM);
              amrex::IntVect itype(AMREX_SPACEDIM);

              ncf.var(lo_names[i] ).get(smallend.begin(), {new_mfi.index(), 0}, {1, AMREX_SPACEDIM});
              ncf.var(hi_names[i] ).get(bigend.begin()  , {new_mfi.index(), 0}, {1, AMREX_SPACEDIM});
              ncf.var(typ_names[i]).get(itype.begin()   , {new_mfi.index(), 0}, {1, AMREX_SPACEDIM});

              amrex::Box fab_box = amrex::Box(smallend, bigend, itype);

              amrex::FArrayBox &fab = new_data[new_mfi.index()];

              for(int k(0); k < fab.nComp(); ++k) {
                 auto dataPtr = fab.dataPtr(k);
                 std::pair<int, int> p {std::pair<int, int>(i, k)};
                 int index = var_map[p];
                 ncf.var(var_names[index]).get(dataPtr, {new_mfi.index(), 0}, {1, npts[index]});
              }
          }
       }
    }

    if (papa.useFixedCoarseGrids()) {
      erf.constructAreaNotToTag();
    }

    erf.post_step_regrid = 0;

    erf.finishConstructor();

  /*
    Deal here with new state descriptor types added, with corresponding
    input_version > 0, if applicable
   */
  for (int i = 0; i < erf.desc_lst.size(); ++i) {
    if (state_in_checkpoint[i] == 0) {
      const amrex::Real ctime = erf.state[i-1].curTime();
      erf.state[i].define(
        erf.geom.Domain(), erf.grids, erf.dmap, erf.desc_lst[i], ctime, papa.dtLevel(erf.level),
        *(erf.m_factory));
      erf.state[i] = erf.state[i - 1];
    }
  }
  erf.buildMetrics();

  // get the elapsed CPU time to now;
  if (erf.level == 0 && amrex::ParallelDescriptor::IOProcessor()) {
    // get elapsed CPU time
    std::ifstream CPUFile;
    std::string FullPathCPUFile = papa.theRestartFile();
    FullPathCPUFile += "/CPUtime";
    CPUFile.open(FullPathCPUFile.c_str(), std::ios::in);
    CPUFile >> erf.previousCPUTimeUsed;
    CPUFile.close();
    amrex::Print() << "read CPU time: " << erf.previousCPUTimeUsed << "\n";
  }

  /* Deprecated: erf.track_grid_losses is not a supported option
  if (erf.track_grid_losses && erf.level == 0) {

    // get the current value of the diagnostic quantities
    std::ifstream DiagFile;
    std::string FullPathDiagFile = papa.theRestartFile();
    FullPathDiagFile += "/Diagnostics";
    DiagFile.open(FullPathDiagFile.c_str(), std::ios::in);

    DiagFile.close();
  }
  */

  /*Not implemented for CUDA
      if (level == 0)
      {
    // get problem-specific stuff -- note all processors do this,
    // eliminating the need for a broadcast
    std::string dir = parent->theRestartFile();

    char * dir_for_pass = new char[dir.size() + 1];
    std::copy(dir.begin(), dir.end(), dir_for_pass);
    dir_for_pass[dir.size()] = '\0';

    int len = dir.size();

    Vector<int> int_dir_name(len);
    for (int j = 0; j < len; j++)
        int_dir_name[j] = (int) dir_for_pass[j];

    AMREX_FORT_PROC_CALL(PROBLEM_RESTART,problem_restart)(int_dir_name.dataPtr(),
    &len);

    delete [] dir_for_pass;

      }*/
}

