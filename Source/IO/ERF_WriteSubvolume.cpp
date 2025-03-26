#include <ERF_EOS.H>
#include <ERF.H>
#include <AMReX_ParmParse.H>

using namespace amrex;

// Write plotfile to disk
void
ERF::WriteSubvolume ()
{
    ParmParse pp("erf.subvol");

    Vector<Real> origin;
    Vector< int> ncell;
    Vector<Real> delta;

    // **************************************************************
    // Read in the origin, number of cells in each dir, and resolution
    // **************************************************************
    pp.getarr("origin",origin,0,AMREX_SPACEDIM);
    pp.getarr("nxnynz", ncell,0,AMREX_SPACEDIM);
    pp.getarr("dxdydz", delta,0,AMREX_SPACEDIM);

    int lev_for_sub = 0;

    bool found = false;
    for (int i = 0; i <= finest_level; i++) {
        if (!found) {
            if (almostEqual(delta[0],geom[i].CellSize(0)) &&
                almostEqual(delta[1],geom[i].CellSize(1)) &&
                almostEqual(delta[2],geom[i].CellSize(2)) ) {

                // amrex::Print() << "XDIR " << delta[0] << " " << geom[i].CellSize(0) << std::endl;
                // amrex::Print() << "YDIR " << delta[1] << " " << geom[i].CellSize(1) << std::endl;
                // amrex::Print() << "ZDIR " << delta[2] << " " << geom[i].CellSize(2) << std::endl;
                amrex::Print() << "Resolution specified matches that of level " << i << std::endl;
                found = true;
                lev_for_sub = i;
            }
        }
    }

    if (!found) {
        amrex::Abort("Resolution specified for subvol does not match the resolution of any of the levels.");
    }

    // **************************************************************
    // Now that we know which level we're at, we can figure out which (i,j,k) the origin corresponds to
    // Note we use 1.0001 as a fudge factor since the division of two reals --> integer will do a floor
    // **************************************************************
    int i0 = static_cast<int>((origin[0] - geom[lev_for_sub].ProbLo(0)) * 1.0001 / delta[0]);
    int j0 = static_cast<int>((origin[1] - geom[lev_for_sub].ProbLo(1)) * 1.0001 / delta[1]);
    int k0 = static_cast<int>((origin[2] - geom[lev_for_sub].ProbLo(2)) * 1.0001 / delta[2]);

    found = false;
    if (almostEqual(geom[lev_for_sub].ProbLo(0)+i0*delta[0],origin[0]) &&
        almostEqual(geom[lev_for_sub].ProbLo(1)+j0*delta[1],origin[1]) &&
        almostEqual(geom[lev_for_sub].ProbLo(2)+k0*delta[2],origin[2]) )
    {
        amrex::Print() << "Specified origin is the lower left corner of cell " << IntVect(i0,j0,k0) << std::endl;
        found = true;
    }

    if (!found) {
        amrex::Abort("Origin specified does not correspond to a node at this level.");
    }

    Box domain(geom[lev_for_sub].Domain());

    Box bx(IntVect(i0,j0,k0),IntVect(i0+ncell[0]-1,j0+ncell[1]-1,k0+ncell[2]-1));
    amrex::Print() << "Box requested is " << bx << std::endl;

    if (!domain.contains(bx))
    {
        amrex::Abort("Box requested is larger than the existing domain");
    }

    Vector<int> cs(3);
    int count = pp.countval("chunk_size");
    if (count > 0) {
        pp.queryarr("chunk_size",cs,0,AMREX_SPACEDIM);
    } else {
        cs[0] = max_grid_size[0][0];
        cs[1] = max_grid_size[0][1];
        cs[2] = max_grid_size[0][2];
    }
    IntVect chunk_size(cs[0],cs[1],cs[2]);

    BoxArray ba(bx);
    ba.maxSize(chunk_size);

    amrex::Print() << "BoxArray is " << ba << std::endl;

    int ncomp_mf = AMREX_SPACEDIM;

    DistributionMapping dm(ba);

    MultiFab mf(ba, dm, ncomp_mf, 0);

    MultiFab mf_cc_vel(grids[lev_for_sub], dmap[lev_for_sub], ncomp_mf, 0);
    average_face_to_cellcenter(mf_cc_vel,0,
                               Array<const MultiFab*,3>{&vars_new[lev_for_sub][Vars::xvel],
                                                        &vars_new[lev_for_sub][Vars::yvel],
                                                        &vars_new[lev_for_sub][Vars::zvel]});

    mf.ParallelCopy(mf_cc_vel,0,0,AMREX_SPACEDIM,0,0);

    std::string subvol_filename = Concatenate(subvol_file, istep[0], 5);

    Vector<std::string> varnames;
    varnames.push_back("x_velocity");
    varnames.push_back("y_velocity");
    varnames.push_back("z_velocity");

    Real time = t_new[lev_for_sub];

    amrex::Print() <<"Writing subvolume into " << subvol_filename << std::endl;
    WriteSingleLevelPlotfile(subvol_filename,mf,varnames,geom[lev_for_sub],time,istep[0]);
}
