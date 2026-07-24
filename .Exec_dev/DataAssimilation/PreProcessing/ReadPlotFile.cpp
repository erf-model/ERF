#include <fstream>
#include <sstream>

#include <AMReX.H>
#include <AMReX_Vector.H>
#include <AMReX_MultiFab.H>
#include <AMReX_PlotFileUtil.H>
#include <AMReX_DistributionMapping.H>
#include <AMReX_BoxArray.H>

#include "ReadPlotFile.H"

#include <vector>
#include <algorithm>  // for std::remove_if, std::isspace

using namespace amrex;
// ------------------------------------------------------------
// Read variable names from a file
// ------------------------------------------------------------
Vector<std::string>
ReadVarNames(const std::string& filename)
{
    Vector<std::string> varnames;
    std::ifstream infile(filename);
    if (!infile.is_open()) {
        amrex::Abort("Cannot open variable file: " + filename);
    }

    std::string line;
    while (std::getline(infile, line)) {
        // trim whitespace
        line.erase(line.begin(),
                   std::find_if(line.begin(), line.end(),
                                [](unsigned char ch){ return !std::isspace(ch); }));
        line.erase(std::find_if(line.rbegin(), line.rend(),
                                [](unsigned char ch){ return !std::isspace(ch); }).base(),
                   line.end());

        if (!line.empty()) varnames.push_back(line);
    }

    if (varnames.empty()) {
        amrex::Abort("Variable file is empty: " + filename);
    }

    return varnames;
}

// Reads the plotfile data into cell cenetred multifab
// Does not fill ghost cells
void
ReadPlotFile(const std::string& var_filename,
             PlotFileData& pf,
             amrex::MultiFab& mf)
{
    // ------------------------------------------------------------
    // Read variable list from file
    // ------------------------------------------------------------
    Vector<std::string> varnames;
    {
        std::ifstream infile(var_filename);
        if (!infile.is_open()) {
            Abort("ReadPlotFile: Cannot open variable file: " + var_filename);
        }

        std::string line;
        while (std::getline(infile, line)) {
            if (!line.empty()) {
                // trim whitespace
                line.erase(0, line.find_first_not_of(" \t\r\n"));
                line.erase(line.find_last_not_of(" \t\r\n") + 1);

                if (!line.empty()) varnames.push_back(line);
            }
        }
    }

    if (varnames.empty()) {
        Abort("ReadPlotFile: No variables found in " + var_filename);
    }

    // ------------------------------------------------------------
    // Open plotfile
    // ------------------------------------------------------------
    const Vector<std::string>& var_names_pf = pf.varNames();

    // ------------------------------------------------------------
    // Validate requested variables
    // ------------------------------------------------------------
    for (auto const& v : varnames) {
        bool found = false;
        for (auto const& vpf : var_names_pf) {
            if (v == vpf) {
                found = true;
                break;
            }
        }
        if (!found) {
            Abort("ReadPlotFile: invalid variable name: " + v);
        }
    }

    // ------------------------------------------------------------
    // Define destination MultiFab (single level only)
    // ------------------------------------------------------------
    const int level = 0;

    BoxArray ba = pf.boxArray(level);
    DistributionMapping dm(ba);

    int ncomp = varnames.size();
    int ngrow = 1;

    mf.define(ba, dm, ncomp, ngrow);

    // ------------------------------------------------------------
    // Copy plotfile data → mf
    // ------------------------------------------------------------
    for (int comp = 0; comp < ncomp; ++comp)
    {
        const MultiFab& src = pf.get(level, varnames[comp]);
        MultiFab::Copy(mf, src, 0, comp, 1, 0);
    }
}

void AverageFineToCoarse(
    const MultiFab& mf_fine,
    MultiFab&       mf_coarse,
    const Geometry& geom_fine,
    const Geometry& geom_coarse)
{
    AMREX_ALWAYS_ASSERT(mf_fine.nComp() == mf_coarse.nComp());

    const int ncomp = mf_fine.nComp();

    // Fine grid spacing
    const Real* dx_f = geom_fine.CellSize();
    const Real* plo_f = geom_fine.ProbLo();

    // Coarse grid spacing
    const Real* dx_c = geom_coarse.CellSize();
    const Real* plo_c = geom_coarse.ProbLo();


    // Store all fine cell centers and values
    struct FineCell
    {
        Real x, y, z;
        Vector<Real> val;
    };

    Vector<FineCell> fine_cells;


    // Collect fine cell centers
    for (MFIter mfi(mf_fine); mfi.isValid(); ++mfi)
    {
        const Box& bx = mfi.validbox();
        auto const& fine_arr = mf_fine.const_array(mfi);

        for (int k = bx.smallEnd(2); k <= bx.bigEnd(2); ++k)
        {
            for (int j = bx.smallEnd(1); j <= bx.bigEnd(1); ++j)
            {
                for (int i = bx.smallEnd(0); i <= bx.bigEnd(0); ++i)
                {
                    FineCell cell;

                    cell.x = plo_f[0] + (i + 0.5)*dx_f[0];
                    cell.y = plo_f[1] + (j + 0.5)*dx_f[1];
                    cell.z = plo_f[2] + (k + 0.5)*dx_f[2];

                    cell.val.resize(ncomp);

                    for (int n = 0; n < ncomp; ++n)
                    {
                        cell.val[n] = fine_arr(i,j,k,n);
                    }

                    fine_cells.push_back(std::move(cell));
                }
            }
        }
    }


    Print() << "Number of fine cells = "
            << fine_cells.size() << "\n";


    // Loop over coarse cells
    for (MFIter mfi(mf_coarse); mfi.isValid(); ++mfi)
    {
        const Box& bx = mfi.validbox();
        auto const& coarse_arr = mf_coarse.array(mfi);


        for (int k = bx.smallEnd(2); k <= bx.bigEnd(2); ++k)
        {
            for (int j = bx.smallEnd(1); j <= bx.bigEnd(1); ++j)
            {
                for (int i = bx.smallEnd(0); i <= bx.bigEnd(0); ++i)
                {

                    Real xlo = plo_c[0] + i*dx_c[0];
                    Real xhi = xlo + dx_c[0];

                    Real ylo = plo_c[1] + j*dx_c[1];
                    Real yhi = ylo + dx_c[1];

                    Real zlo = plo_c[2] + k*dx_c[2];
                    Real zhi = zlo + dx_c[2];


                    Vector<Real> sum(ncomp,0.0);
                    int count = 0;


                    // Search fine cells inside coarse cell
                    for (const auto& fc : fine_cells)
                    {
                        if (fc.x >= xlo && fc.x < xhi &&
                            fc.y >= ylo && fc.y < yhi &&
                            fc.z >= zlo && fc.z < zhi)
                        {
                            for (int n = 0; n < ncomp; ++n)
                            {
                                sum[n] += fc.val[n];
                            }

                            count++;
                        }
                    }


                    // Average
                    if (count > 0)
                    {
                        for (int n = 0; n < ncomp; ++n)
                        {
                            coarse_arr(i,j,k,n)
                                = sum[n]/Real(count);
                        }
                    }
                    else
                    {
                        // no fine cell found
                        for (int n = 0; n < ncomp; ++n)
                        {
                            coarse_arr(i,j,k,n) = 0.0;
                        }
                    }
                }
            }
        }
    }
}


void ApplyNeumannBCs(const Geometry& geom,
                     MultiFab& mf_cc)
{

     // -------------------------------------------------
    // 2. Fill interior + periodic ghost cells
    // -------------------------------------------------
    mf_cc.FillBoundary(geom.periodicity());
    // -------------------------------------------------
    // 3. Apply FOExtrap (Neumann) at domain boundaries
    // -------------------------------------------------
    const Box& domain = geom.Domain();

    for (MFIter mfi(mf_cc, TilingIfNotGPU()); mfi.isValid(); ++mfi)
    {
        const Box& gbx = mfi.growntilebox();   // includes ghost cells
        const Box& vbx = mfi.validbox();

        auto const& arr = mf_cc.array(mfi);
        int ncomp = mf_cc.nComp();

        ParallelFor(gbx, ncomp,
        [=] AMREX_GPU_DEVICE (int i, int j, int k, int n)
        {
            if (vbx.contains(i,j,k)) return;

            int ii = i;
            int jj = j;
            int kk = k;

            // Clamp to domain interior (FOExtrap)
            ii = amrex::max(domain.smallEnd(0),
                 amrex::min(i, domain.bigEnd(0)));

            jj = amrex::max(domain.smallEnd(1),
                 amrex::min(j, domain.bigEnd(1)));

            kk = amrex::max(domain.smallEnd(2),
                 amrex::min(k, domain.bigEnd(2)));

            arr(i,j,k,n) = arr(ii,jj,kk,n);
        });
    }
}

void WriteCustomDataFile(const Geometry& geom,
                         const MultiFab& mf_cc,
                         const std::string& filename_custom)
{
    std::cout << "coarse box number is " << mf_cc.local_size() << std::endl;

    AMREX_ALWAYS_ASSERT(mf_cc.nComp() >= 1);
    AMREX_ALWAYS_ASSERT(mf_cc.local_size() == 1); // single box assumed

    const int ncomp = mf_cc.nComp();
    const int ng    = mf_cc.nGrow();

    const Box grown_domain = amrex::grow(geom.Domain(), ng);
    const auto lo = lbound(grown_domain);
    const auto hi = ubound(grown_domain);

    const auto problo = geom.ProbLoArray();
    const auto probhi = geom.ProbHiArray();
    const auto dx     = geom.CellSizeArray();

    // extend physical domain to include ghosts
    amrex::GpuArray<double,3> problo_ext, probhi_ext;
    for (int d=0; d<3; ++d) {
        problo_ext[d] = problo[d] - ng*dx[d];
        probhi_ext[d] = probhi[d] + ng*dx[d];
    }

    std::ofstream ofs(filename_custom, std::ios::binary);
    if (!ofs.is_open()) Abort("Failed to open file for writing");

    // header
    const int nx = hi.x - lo.x + 1;
    const int ny = hi.y - lo.y + 1;
    const int nz = hi.z - lo.z + 1;

    ofs.write(reinterpret_cast<const char*>(&nx), sizeof(int));
    ofs.write(reinterpret_cast<const char*>(&ny), sizeof(int));
    ofs.write(reinterpret_cast<const char*>(&nz), sizeof(int));
    ofs.write(reinterpret_cast<const char*>(&ng), sizeof(int));
    ofs.write(reinterpret_cast<const char*>(&ncomp), sizeof(int));

    for (int d=0; d<3; ++d) ofs.write(reinterpret_cast<const char*>(&problo_ext[d]), sizeof(double));
    for (int d=0; d<3; ++d) ofs.write(reinterpret_cast<const char*>(&probhi_ext[d]), sizeof(double));

    // data
    const auto& arr = mf_cc.const_array(0);
    for (int k=lo.z; k<=hi.z; ++k)
    for (int j=lo.y; j<=hi.y; ++j)
    for (int i=lo.x; i<=hi.x; ++i)
    {
        double x = problo[0] + (i + 0.5) * dx[0];
        double y = problo[1] + (j + 0.5) * dx[1];
        double z = problo[2] + (k + 0.5) * dx[2];

        ofs.write(reinterpret_cast<const char*>(&x), sizeof(double));
        ofs.write(reinterpret_cast<const char*>(&y), sizeof(double));
        ofs.write(reinterpret_cast<const char*>(&z), sizeof(double));

        for (int n=0; n<ncomp; ++n)
        {
            double val = arr(i,j,k,n);
            ofs.write(reinterpret_cast<const char*>(&val), sizeof(double));
        }
    }

    ofs.close();
}
