#include "AMReX_Gpu.H"
#include "AMReX_ParmParse.H"
#include "AMReX_PlotFileUtil.H"
#include "AMReX_MultiFabUtil.H"
#include "ERF_WriteBndryPlanes.H"
#include "IndexDefines.H"
#include "Derive.H"

// Default to level 0
int WriteBndryPlanes::bndry_lev = 0;

WriteBndryPlanes::WriteBndryPlanes(amrex::Vector<amrex::BoxArray>& grids,
                                   amrex::Vector<amrex::Geometry>& geom): m_geom(geom)
{
    amrex::ParmParse pp("erf");

    // User-specified region is given in physical coordinates, not index space
    std::vector<amrex::Real> box_lo(3), box_hi(3);
    pp.getarr("bndry_output_box_lo",box_lo,0,2);
    pp.getarr("bndry_output_box_hi",box_hi,0,2);

    // If the target area is contained at a finer level, use the finest data possible
    for (int ilev = 0; ilev < grids.size(); ilev++) {

        auto const dxi = geom[ilev].InvCellSizeArray();
        const amrex::Box& domain = m_geom[ilev].Domain();

        // We create the smallest box that contains all of the cell centers
        // in the physical region specified
        int ilo = static_cast<int>(amrex::Math::floor(box_lo[0] * dxi[0])+.5);
        int jlo = static_cast<int>(amrex::Math::floor(box_lo[1] * dxi[1])+.5);
        int ihi = static_cast<int>(amrex::Math::floor(box_hi[0] * dxi[0])+.5)-1;
        int jhi = static_cast<int>(amrex::Math::floor(box_hi[1] * dxi[1])+.5)-1;

        // Map this to index space -- for now we do no interpolation
        target_box.setSmall(amrex::IntVect(ilo,jlo,0));
        target_box.setBig(amrex::IntVect(ihi,jhi,domain.bigEnd(2)));

        // Test if the target box at this level fits in the grids at this level
        amrex::Box gbx = target_box; gbx.grow(amrex::IntVect(1,1,0));
        if (grids[ilev].contains(gbx)) bndry_lev = ilev;

        if (ilev == 0) {
            AMREX_ALWAYS_ASSERT(grids[ilev].contains(gbx));
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

void WriteBndryPlanes::write_planes(const int t_step, const amrex::Real time,
                                    amrex::Vector<amrex::Vector<amrex::MultiFab>>& vars_new)
{
    BL_PROFILE("ERF::WriteBndryPlanes::write_planes");

    auto& lev_new = vars_new[bndry_lev];
    amrex::MultiFab& S    = vars_new[bndry_lev][Vars::cons];
    amrex::MultiFab& xvel = vars_new[bndry_lev][Vars::xvel];
    amrex::MultiFab& yvel = vars_new[bndry_lev][Vars::yvel];
    amrex::MultiFab& zvel = vars_new[bndry_lev][Vars::zvel];
    //FillPatch(bndry_lev, time, S   , 0, 1, Vars::cons);
    //FillPatch(bndry_lev, time, xvel, 0, 1, Vars::xvel);
    //FillPatch(bndry_lev, time, yvel, 0, 1, Vars::yvel);
    //FillPatch(bndry_lev, time, zvel, 0, 1, Vars::zvel);

    const std::string chkname =
        m_filename + amrex::Concatenate("/bndry_output", t_step);

    //amrex::Print() << "Writing boundary planes at time " << time << std::endl;

    const std::string level_prefix = "Level_";
    amrex::PreBuildDirectorHierarchy(chkname, level_prefix, 1, true);

    // note: by using the entire domain box we end up using 1 processor
    // to hold all boundaries
    amrex::BoxArray ba(target_box);
    amrex::DistributionMapping dm{ba};

    for (int i = 0; i < m_var_names.size(); i++)
    {
        std::string var_name = m_var_names[i];

        if (var_name == "density")
        {
            int ncomp = 1;
            amrex::BndryRegister bndry(ba, dm, m_in_rad, m_out_rad, m_extent_rad, ncomp);
            bndry.copyFrom(S, Cons::Rho, 0, 0, ncomp, m_geom[bndry_lev].periodicity());
            std::string filename = amrex::MultiFabFileFullPrefix(bndry_lev, chkname, level_prefix, m_var_names[i]);
            for (amrex::OrientationIter oit; oit != nullptr; ++oit) {
                auto ori = oit();
                if (ori.coordDir() < 2) {
                    std::string facename = amrex::Concatenate(filename + '_', ori, 1);
                    bndry[ori].write(facename);
                }
            }

        } else if (var_name == "temperature") {

            int ncomp = 1;
            amrex::MultiFab Temp(S.boxArray(),S.DistributionMap(),ncomp,0);
            for (amrex::MFIter mfi(Temp, amrex::TilingIfNotGPU()); mfi.isValid(); ++mfi)
            {
                const amrex::Box& bx = mfi.tilebox();
                auto& dfab    = Temp[mfi];
                auto& sfab    = S[mfi];
                derived::erf_dertemp(bx, dfab, 0, 1, sfab, m_geom[bndry_lev], time, nullptr, bndry_lev);
            }
            amrex::BndryRegister bndry(ba, dm, m_in_rad, m_out_rad, m_extent_rad, ncomp);
            bndry.copyFrom(Temp, 0, 0, 0, ncomp, m_geom[bndry_lev].periodicity());
            std::string filename = amrex::MultiFabFileFullPrefix(bndry_lev, chkname, level_prefix, m_var_names[i]);
            for (amrex::OrientationIter oit; oit != nullptr; ++oit) {
                auto ori = oit();
                if (ori.coordDir() < 2) {
                    std::string facename = amrex::Concatenate(filename + '_', ori, 1);
                    bndry[ori].write(facename);
                }
            }

        } else if (var_name == "velocity") {

            amrex::MultiFab Vel(S.boxArray(), S.DistributionMap(), 3, m_out_rad);
            int ncomp_V = Vel.nComp();
            average_face_to_cellcenter(Vel,0,amrex::Array<const amrex::MultiFab*,3>{&xvel,&yvel,&zvel});

            amrex::BndryRegister bndry_V(ba, dm, m_in_rad, m_out_rad, m_extent_rad, ncomp_V);
            bndry_V.copyFrom(Vel, 0, 0, 0, ncomp_V, m_geom[bndry_lev].periodicity());
            std::string filename = amrex::MultiFabFileFullPrefix(bndry_lev, chkname, level_prefix, m_var_names[i]);
            for (amrex::OrientationIter oit; oit != nullptr; ++oit) {
                auto ori = oit();
                if (ori.coordDir() < 2) {
                    std::string facename = amrex::Concatenate(filename + '_', ori, 1);
                    bndry_V[ori].write(facename);
                }
            }
        } else {
            //amrex::Print() << "Trying to write planar output for " << var_name << std::endl;
            amrex::Error("Don't know how to output this variable");
        }
    } // loop over num_vars

    // Writing time.dat
    if (amrex::ParallelDescriptor::IOProcessor()) {
        std::ofstream oftime(m_time_file, std::ios::out | std::ios::app);
        oftime << t_step << ' ' << time << '\n';
        oftime.close();
    }
}

