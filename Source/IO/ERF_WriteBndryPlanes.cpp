/**
 * \file ERF_WriteBndryPlanes.cpp
 */
#include "AMReX_Gpu.H"
#include "AMReX_ParmParse.H"
#include "AMReX_PlotFileUtil.H"
#include "AMReX_MultiFabUtil.H"
#include <algorithm>
#include <cmath>
#include <filesystem>
#include <fstream>
#include <limits>
#include <sstream>
#include <string>
#include <system_error>
#include "ERF_WriteBndryPlanes.H"
#include "ERF_IndexDefines.H"
#include "ERF_Derive.H"

using namespace amrex;

/**
 * Copies the contents of one boundary register to another in either x or y dimensions.
 *
 * @param oit Orientation iterator used for determining which dimension we are to copy
 * @param b1 Boundary register containing the source data
 * @param b2 Boundary register containing the destination data
 */
void br_shift (OrientationIter oit, const BndryRegister& b1, BndryRegister& b2)
{
    auto ori = oit();
    int ncomp = b1[ori].nComp();
    if (ori.coordDir() < 2) {
#ifdef AMREX_USE_OMP
#pragma omp parallel if (Gpu::notInLaunchRegion())
#endif
        for (FabSetIter bfsi(b1[ori]); bfsi.isValid(); ++bfsi) {
            int idx = bfsi.index();
            const Box& bx1 = b1[ori].boxArray()[idx];
            const Box& bx2 = b2[ori].boxArray()[idx];

            // Copy onto a boundary register based on a box starting at (0,0,0)
            Array4<Real>       dest_arr = b2[ori][idx].array();
            Array4<Real const>  src_arr = b1[ori][idx].const_array();
            int ioff = bx1.smallEnd(0) - bx2.smallEnd(0);
            int joff = bx1.smallEnd(1) - bx2.smallEnd(1);
            ParallelFor(bx2, ncomp,
            [=] AMREX_GPU_DEVICE(int i, int j, int k, int n) noexcept
            {
                dest_arr(i,j,k,n) = src_arr(i+ioff,j+joff,k,n);
            });
        }
    }
}

// Default to level 0
int WriteBndryPlanes::bndry_lev = 0;

/**
 * Constructor for the WriteBndryPlanes class for writing the contents of boundary planes to an output file
 *
 * @param grids Vector of BoxArrays containing the grids at each level in the AMR
 * @param geom Vector of Geometry containing the geometry at each level in the AMR
 */
WriteBndryPlanes::WriteBndryPlanes (Vector<BoxArray>& grids,
                                    Vector<Geometry>& geom) : m_geom(geom)
{
    ParmParse pp("erf");

    // Get the radius inside the domain
    pp.queryAdd("in_rad",m_in_rad);

    // User-specified region is given in physical coordinates, not index space
    std::vector<Real> box_lo(3), box_hi(3);
    pp.getarr("bndry_output_box_lo",box_lo,0,2);
    pp.getarr("bndry_output_box_hi",box_hi,0,2);

    // If the target area is contained at a finer level, use the finest data possible.
    // target_box is expressed in the index space of bndry_lev, so it must only be set
    // together with bndry_lev; each level's candidate box is kept in a local variable.
    for (int ilev = 0; ilev < grids.size(); ilev++) {

        const Real* xLo   = m_geom[ilev].ProbLo();
        auto const dxi    = geom[ilev].InvCellSizeArray();
        const Box& domain = m_geom[ilev].Domain();

        // We create the smallest box that contains all of the cell centers
        // in the physical region specified
        int ilo = static_cast<int>(Math::floor((box_lo[0] - xLo[0]) * dxi[0] + Real(.5)));
        int jlo = static_cast<int>(Math::floor((box_lo[1] - xLo[1]) * dxi[1] + Real(.5)));
        int ihi = static_cast<int>(Math::floor((box_hi[0] - xLo[0]) * dxi[0] + Real(.5)))-1;
        int jhi = static_cast<int>(Math::floor((box_hi[1] - xLo[1]) * dxi[1] + Real(.5)))-1;

        // Map this to index space -- for now we do no interpolation
        Box box_this_lev(IntVect(ilo,jlo,0),
                         IntVect(ihi,jhi,domain.bigEnd(2)));

        // Test if the target box at this level fits in the grids at this level
        Box gbx = box_this_lev; gbx.grow(IntVect(1,1,0));

        // Ensure that the box is no larger than can fit in the (periodically grown) domain
        // at level 0
        if (ilev == 0) {
            Box per_grown_domain(domain);
            int growx = (geom[0].isPeriodic(0)) ? 1 : 0;
            int growy = (geom[0].isPeriodic(1)) ? 1 : 0;
            per_grown_domain.grow(IntVect(growx,growy,0));
            /*
            if (!per_grown_domain.contains(gbx))
                Error("WriteBndryPlanes: Requested box is too large to fill");
            */

            // Default to level 0
            target_box = box_this_lev;
        }

        // Use this level's data only if its grids fully contain the requested region
        if (grids[ilev].contains(gbx)) {
            bndry_lev  = ilev;
            target_box = box_this_lev;
        }
    }

    // The folder "m_filename" will contain the time series of data and the time.dat file
    pp.get("bndry_output_planes_file", m_filename);

    m_time_file = m_filename + "/time.dat";

    if (pp.contains("bndry_output_var_names"))
    {
        int num_vars = pp.countval("bndry_output_var_names");
        m_var_names.resize(num_vars);
        pp.queryarr("bndry_output_var_names",m_var_names,0,num_vars);
    }
}

/**
 * Number of leading rows of time.dat that a run starting at step start_step keeps.
 *
 * A restart replays the steps after its checkpoint, so the rows for those steps -- written by a
 * run that went on past the checkpoint -- are dropped and written again. Rows also stop being
 * kept at the first one whose step or time does not increase, or that is later than the restart,
 * or that is at the restart step but not at the restart time: those were not written by the run
 * being continued. A fresh start begins a new series.
 *
 * @param steps Step column of time.dat, in file order
 * @param times Time column of time.dat, in file order
 * @param start_step Step the run starts at
 * @param start_time Time the run starts at, as written in time.dat
 * @param restarting Whether the run starts from a checkpoint
 */
int WriteBndryPlanes::rows_to_keep (const Vector<int>& steps, const Vector<double>& times,
                                    const int start_step, const double start_time,
                                    const bool restarting)
{
    if (!restarting) { return 0; }
    // The times of the run being continued are the restart's own to the last bit (time.dat and
    // the checkpoint both hold 17 digits); a few units in the last place only allow for a
    // time.dat written by hand, and stay far below any plane spacing also at a calendar clock
    const double tol = 16.0 * std::numeric_limits<double>::epsilon() * std::max(1.0, std::abs(start_time));
    int nkeep = 0;
    for (int i = 0; i < steps.size(); ++i) {
        if (steps[i] > start_step || times[i] > start_time + tol) { break; }
        if (i > 0 && (steps[i] <= steps[i-1] || times[i] <= times[i-1])) { break; }
        if (steps[i] == start_step && std::abs(times[i] - start_time) > tol) { break; }
        ++nkeep;
    }
    return nkeep;
}

/**
 * Decide which rows of time.dat a run starting at step start_step continues from, and say whether
 * the series already holds the plane at that step. Only the I/O rank reads the file; every rank
 * gets the answer. The file itself is rewritten by write_planes, just before the run's first row,
 * so a run that stops before its first plane leaves the series as it found it.
 *
 * @param start_step Step the run starts at (the restart step, or 0)
 * @param start_time Time the run starts at, as written in time.dat
 * @param restarting Whether the run starts from a checkpoint
 */
bool WriteBndryPlanes::start_series (const int start_step, const double start_time,
                                     const bool restarting)
{
    int has_start_plane = 0;
    if (ParallelDescriptor::IOProcessor()) {
        std::ifstream in(m_time_file);
        if (in.good()) {
            // Every line counts, so that a line that does not parse -- a partial last row left
            // by a run that stopped while writing it, say -- is dropped with the rows after it
            Vector<std::string> rows;
            Vector<int> steps;
            Vector<double> times;
            int nlines = 0;
            std::string line;
            while (std::getline(in, line)) {
                ++nlines;
                std::istringstream row(line);
                int step;
                double time;
                if (static_cast<int>(rows.size()) + 1 < nlines || !(row >> step >> time)) { continue; }
                rows.push_back(line);
                steps.push_back(step);
                times.push_back(time);
            }
            in.close();

            const int nkeep = rows_to_keep(steps, times, start_step, start_time, restarting);
            has_start_plane = (nkeep > 0 && steps[nkeep-1] == start_step) ? 1 : 0;

            // A row whose step or time does not increase, anywhere in the file, or a row at or
            // before the restart that does not continue the kept rows, cannot come from one run
            bool broken = false;
            for (int i = 1; i < rows.size(); ++i) {
                if (steps[i] <= steps[i-1] || times[i] <= times[i-1]) { broken = true; }
            }
            if (restarting && nkeep < rows.size() && steps[nkeep] <= start_step) { broken = true; }

            if (restarting && broken) {
                Warning("WriteBndryPlanes: " + m_time_file + " does not hold one run's series: a step "
                        "or time does not increase, or a row at or before the restart step " +
                        std::to_string(start_step) + " is not at a time of this run. Only its first " +
                        std::to_string(nkeep) + " rows are kept. A series that an earlier ERF version "
                        "wrote across a restart has such a row, and its step-0 plane holds the "
                        "restart state; write that series again.");
            } else if (nkeep < nlines) {
                if (restarting) {
                    Print() << "WriteBndryPlanes: keeping " << nkeep << " of the " << nlines
                            << " rows of " << m_time_file << " up to the restart step "
                            << start_step << ", from this run's first plane on" << std::endl;
                } else {
                    Print() << "WriteBndryPlanes: a fresh start begins a new series; the "
                            << nlines << " rows of " << m_time_file
                            << " are dropped with this run's first plane" << std::endl;
                }
            }

            // Rewrite even when every row is kept: a last row without its newline would
            // otherwise run into the next row written
            m_kept_rows.clear();
            for (int i = 0; i < nkeep; ++i) {
                m_kept_rows.push_back(rows[i]);
            }
            m_rewrite_time_file = true;
        }
    }
    ParallelDescriptor::Bcast(&has_start_plane, 1, ParallelDescriptor::IOProcessorNumber(),
                              ParallelDescriptor::Communicator());
    return has_start_plane == 1;
}

/**
 * Function to write the specified grid data to an output file
 *
 * @param t_step Timestep number
 * @param time Current time
 * @param vars_new Grid data for all variables across the AMR hierarchy
 * @param is_moist Whether moist state components are included in the output
 */
void WriteBndryPlanes::write_planes (const int t_step, const double time_d,
                                     Vector<Vector<MultiFab>>& vars_new,
                                     bool is_moist)
{
    BL_PROFILE("ERF::WriteBndryPlanes::write_planes");

    MultiFab& S    = vars_new[bndry_lev][Vars::cons];
    MultiFab& xvel = vars_new[bndry_lev][Vars::xvel];
    MultiFab& yvel = vars_new[bndry_lev][Vars::yvel];
    MultiFab& zvel = vars_new[bndry_lev][Vars::zvel];

    const std::string chkname =
        m_filename + Concatenate("/bndry_output", t_step);

    //Print() << "Writing boundary planes at time " << time_d << std::endl;

    const std::string level_prefix = "Level_";
    PreBuildDirectorHierarchy(chkname, level_prefix, 1, true);

    // note: by using the entire domain box we end up using 1 processor
    // to hold all boundaries
    BoxArray ba(target_box);
    DistributionMapping dm{ba};

    IntVect new_hi = target_box.bigEnd() - target_box.smallEnd();
    Box target_box_shifted(IntVect(0,0,0),new_hi);
    BoxArray ba_shifted(target_box_shifted);

    for (int i = 0; i < m_var_names.size(); i++)
    {
        std::string var_name = m_var_names[i];
        std::string filename = MultiFabFileFullPrefix(bndry_lev, chkname, level_prefix, var_name);

        int ncomp;
        if (var_name == "velocity") {
            ncomp = AMREX_SPACEDIM;
        } else {
            ncomp = 1;
        }

        BndryRegister bndry        (ba        , dm, m_in_rad, m_out_rad, m_extent_rad, ncomp);
        BndryRegister bndry_shifted(ba_shifted, dm, m_in_rad, m_out_rad, m_extent_rad, ncomp);

        int nghost = 0;
        if (var_name == "density")
        {
            bndry.copyFrom(S, nghost, Rho_comp, 0, ncomp, m_geom[bndry_lev].periodicity());

        } else if (var_name == "temperature") {
            MultiFab Temp(S.boxArray(),S.DistributionMap(),ncomp,0);
            for (MFIter mfi(Temp, TilingIfNotGPU()); mfi.isValid(); ++mfi)
            {
                const Box& bx = mfi.tilebox();
                Real time = static_cast<Real>(time_d);
                if (is_moist) {
                    // NOTE: we send in S[mfi] where we should send in (*z_phys_cc[lev])[mfi] because we know this routine doesn't use it
                    derived::erf_dermoisttemp(bx, Temp[mfi], 0, 1, S[mfi], S[mfi], m_geom[bndry_lev], time, nullptr, bndry_lev);
                } else {
                    // NOTE: we send in S[mfi] where we should send in (*z_phys_cc[lev])[mfi] because we know this routine doesn't use it
                    derived::erf_dertemp(bx, Temp[mfi], 0, 1, S[mfi], S[mfi], m_geom[bndry_lev], time, nullptr, bndry_lev);
                }
            }
            bndry.copyFrom(Temp, nghost, 0, 0, ncomp, m_geom[bndry_lev].periodicity());
        } else if (var_name == "theta") {
            MultiFab Temp(S.boxArray(),S.DistributionMap(),ncomp,0);
            for (MFIter mfi(Temp, TilingIfNotGPU()); mfi.isValid(); ++mfi)
            {
                const Box& bx = mfi.tilebox();
                derived::erf_derrhodivide(bx, Temp[mfi], S[mfi], RhoTheta_comp);
            }
            bndry.copyFrom(Temp, nghost, 0, 0, ncomp, m_geom[bndry_lev].periodicity());
        } else if (var_name == "scalar") {
            MultiFab Temp(S.boxArray(),S.DistributionMap(),ncomp,0);
            for (MFIter mfi(Temp, TilingIfNotGPU()); mfi.isValid(); ++mfi)
            {
                const Box& bx = mfi.tilebox();
                derived::erf_derrhodivide(bx, Temp[mfi], S[mfi], RhoScalar_comp);
            }
            bndry.copyFrom(Temp, nghost, 0, 0, ncomp, m_geom[bndry_lev].periodicity());
        } else if (var_name == "ke") {
            MultiFab Temp(S.boxArray(),S.DistributionMap(),ncomp,0);
            for (MFIter mfi(Temp, TilingIfNotGPU()); mfi.isValid(); ++mfi)
            {
                const Box& bx = mfi.tilebox();
                derived::erf_derrhodivide(bx, Temp[mfi], S[mfi], RhoKE_comp);
            }
            bndry.copyFrom(Temp, nghost, 0, 0, ncomp, m_geom[bndry_lev].periodicity());
        } else if (var_name == "qv") {
            MultiFab Temp(S.boxArray(),S.DistributionMap(),ncomp,0);
            if (S.nComp() > RhoQ2_comp) {
                for (MFIter mfi(Temp, TilingIfNotGPU()); mfi.isValid(); ++mfi)
                {
                    const Box& bx = mfi.tilebox();
                    derived::erf_derrhodivide(bx, Temp[mfi], S[mfi], RhoQ1_comp);
                }
            } else {
                Temp.setVal(0.);
            }
            bndry.copyFrom(Temp, nghost, 0, 0, ncomp, m_geom[bndry_lev].periodicity());
        } else if (var_name == "qc") {
            MultiFab Temp(S.boxArray(),S.DistributionMap(),ncomp,0);
            if (S.nComp() > RhoQ2_comp) {
                for (MFIter mfi(Temp, TilingIfNotGPU()); mfi.isValid(); ++mfi)
                {
                    const Box& bx = mfi.tilebox();
                    derived::erf_derrhodivide(bx, Temp[mfi], S[mfi], RhoQ2_comp);
                }
            } else {
                Temp.setVal(0.);
            }
            bndry.copyFrom(Temp, nghost, 0, 0, ncomp, m_geom[bndry_lev].periodicity());
        } else if (var_name == "velocity") {
            MultiFab Vel(S.boxArray(), S.DistributionMap(), 3, m_out_rad);
            average_face_to_cellcenter(Vel,0,Array<const MultiFab*,3>{&xvel,&yvel,&zvel});
            bndry.copyFrom(Vel, nghost, 0, 0, ncomp, m_geom[bndry_lev].periodicity());
        } else {
            //Print() << "Trying to write planar output for " << var_name << std::endl;
            Error("Don't know how to output this variable");
        }

        for (OrientationIter oit; oit != nullptr; ++oit) {
            auto ori = oit();
            if (ori.coordDir() < 2) {
                std::string facename = Concatenate(filename + '_', ori, 1);
                br_shift(oit, bndry, bndry_shifted);
                bndry_shifted[ori].write(facename);
            }
        }

    } // loop over num_vars

    // Writing time.dat
    if (ParallelDescriptor::IOProcessor()) {
        // The rows this run continues from replace the file before its first row, through a
        // temporary file so that a crash leaves either the old file or the new one
        if (m_rewrite_time_file) {
            const std::string tmp_file = m_time_file + ".tmp";
            std::ofstream out(tmp_file, std::ios::out | std::ios::trunc);
            for (const auto& row : m_kept_rows) {
                out << row << '\n';
            }
            out.close();
            if (out.fail()) {
                Abort("WriteBndryPlanes: cannot write " + tmp_file);
            }
            std::error_code ec;
            std::filesystem::rename(tmp_file, m_time_file, ec);
            if (ec) {
                Abort("WriteBndryPlanes: cannot replace " + m_time_file + " with " + tmp_file + ": " + ec.message());
            }
            m_rewrite_time_file = false;
            m_kept_rows.clear();
        }
        std::ofstream oftime(m_time_file, std::ios::out | std::ios::app);
        oftime << std::setprecision(17) << t_step << ' ' << time_d << '\n';
        oftime.close();
    }
}
