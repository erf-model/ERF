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


void
CreateNodalMultiFabFromCellCenteredMultiFab (MultiFab& mf_nc,      // output nodal MF
                                             MultiFab& mf_cc,      // input cell-centered MF (coarse)
                                             const Geometry& geom)
{

    // -------------------------------------------------
    // 1. Build nodal MultiFab if not already defined
    // -------------------------------------------------
    if (!mf_nc.isDefined())
    {
        BoxArray ba_nd = amrex::convert(mf_cc.boxArray(),
                                        IntVect::TheNodeVector());

        mf_nc.define(ba_nd,
                     mf_cc.DistributionMap(),
                     mf_cc.nComp(),
                     0);   // nodal MF typically needs no ghosts
    }
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

    AMREX_ALWAYS_ASSERT(mf_nc.nComp() == mf_cc.nComp());

    int ncomp = mf_cc.nComp();

    for (MFIter mfi(mf_nc, TilingIfNotGPU()); mfi.isValid(); ++mfi)
    {
        const Box& bx = mfi.tilebox();       // nodal valid box

        auto const& arr_nc = mf_nc.array(mfi);
        auto const& arr_cc = mf_cc.array(mfi);

        ParallelFor(bx, ncomp,
        [=] AMREX_GPU_DEVICE(int i, int j, int k, int n)
        {
            // 3D: eight surrounding cells
            arr_nc(i,j,k,n) =
                0.125 * ( arr_cc(i-1,j-1,k-1,n) + arr_cc(i  ,j-1,k-1,n)
                        + arr_cc(i-1,j  ,k-1,n) + arr_cc(i  ,j  ,k-1,n)
                        + arr_cc(i-1,j-1,k  ,n) + arr_cc(i  ,j-1,k  ,n)
                        + arr_cc(i-1,j  ,k  ,n) + arr_cc(i  ,j  ,k  ,n) );
        });
    }
}




void
CreateCellCenteredMultiFabFromNodalMultiFab(MultiFab& mf_cc_tmp,   // cell-centered output
                                            MultiFab& mf_nc)   // node-centered input
{

    // -------------------------------------------------
    // 1. Define cell-centered MF with same structure
    // -------------------------------------------------
    BoxArray ba_cc = amrex::convert(mf_nc.boxArray(),
                                    IntVect::TheCellVector());

    mf_cc_tmp.define(ba_cc,
                     mf_nc.DistributionMap(),
                     mf_nc.nComp(),
                     0);   // no ghosts needed for simple averaging

    int ncomp = mf_cc_tmp.nComp();

    for (MFIter mfi(mf_cc_tmp, TilingIfNotGPU()); mfi.isValid(); ++mfi)
    {
        const Box& cbox = mfi.tilebox();   // cell-centered valid box
        auto const& arr_cc = mf_cc_tmp.array(mfi);
        auto const& arr_nc = mf_nc.array(mfi);

        ParallelFor(cbox, ncomp,
        [=] AMREX_GPU_DEVICE(int i, int j, int k, int n)
        {
            // Average 8 surrounding nodes
            arr_cc(i,j,k,n) = 0.125 * ( arr_nc(i,j,k,n)       // corner 0
                              + arr_nc(i+1,j,k,n)     // corner 1
                              + arr_nc(i,j+1,k,n)     // corner 2
                              + arr_nc(i+1,j+1,k,n)   // corner 3
                              + arr_nc(i,j,k+1,n)     // corner 4
                              + arr_nc(i+1,j,k+1,n)   // corner 5
                              + arr_nc(i,j+1,k+1,n)   // corner 6
                              + arr_nc(i+1,j+1,k+1,n) // corner 7
                              );
        });
    }
}


enum class BoundType { Lo, Hi };
enum class MultiFabType { CC, NC };

IntVect
find_bound_idx(const Real& x, const Real& y, const Real& z,
               const BoxList& bl_coarse, const Geometry& geom_coarse,
               BoundType bound_type)
{
    const auto prob_lo_coarse  = geom_coarse.ProbLoArray();
    const auto dx_coarse       = geom_coarse.CellSizeArray();

    int i, j, k;

    if (bound_type == BoundType::Lo) {
        i = static_cast<int>(std::floor((x - prob_lo_coarse[0]) / dx_coarse[0]));
        j = static_cast<int>(std::floor((y - prob_lo_coarse[1]) / dx_coarse[1]));
        k = static_cast<int>(std::floor((z - prob_lo_coarse[2]) / dx_coarse[2]));
    } else { // BoundType::Hi
        i = static_cast<int>(std::ceil((x - prob_lo_coarse[0]) / dx_coarse[0]));
        j = static_cast<int>(std::ceil((y - prob_lo_coarse[1]) / dx_coarse[1]));
        k = static_cast<int>(std::ceil((z - prob_lo_coarse[2]) / dx_coarse[2]));
    }

    IntVect idx(i, j, k);

    for (const auto& b : bl_coarse) {
        if (b.contains(idx)) {
            return idx;
        }
    }

    amrex::Print() << x << " " << y << " " << z << " " << idx << std::endl;

    amrex::Print() << "Printing BoxList (coarse):\n";
for (const auto& b : bl_coarse) {
    amrex::Print() << b << "\n";
}

    amrex::Abort("Bound index not found in any box in BoxList!");
    return IntVect::TheZeroVector(); // unreachable if Abort
}


void
GetCoarseMultiFabOnFineDMap(const Geometry& geom_coarse,
                            const Geometry& geom_fine,
                            const MultiFab& mf_nc_coarse,
                            const MultiFab& mf_cc_fine,
                            MultiFab& coarse_multifab_on_fine_dmap)
{
    BoxList bl_coarse = mf_nc_coarse.boxArray().boxList();
    BoxList bl_fine   = mf_cc_fine.boxArray().boxList();

    const auto prob_lo_fine  = geom_fine.ProbLoArray();
    const auto dx_fine       = geom_fine.CellSizeArray();

    for (auto& b : bl_fine) {
        // You look at the lo corner of b, and find out the lowest cell in
        // coarse mutlifab data you need for the interpolation. That gives
        // you the lo corner of the new b. Similarly, you can find out the
        // hi corner of the new b. For cells outside the coarse multifab data's
        // bounding data, it's up to you. You probably want to use a biased
        // interpolation stencil.

        // Get the cell indices of the bottom corner and top corner
        const IntVect& lo_fine = b.smallEnd();  // Lower corner (inclusive)
        const IntVect& hi_fine = b.bigEnd();    // Upper corner (inclusive)

        Real x = prob_lo_fine[0] + lo_fine[0] * dx_fine[0];
        Real y = prob_lo_fine[1] + lo_fine[1] * dx_fine[1];
        Real z = prob_lo_fine[2] + lo_fine[2] * dx_fine[2];

        auto idx_lo = find_bound_idx(x, y, z, bl_coarse, geom_coarse, BoundType::Lo);


        x = prob_lo_fine[0] + hi_fine[0] * dx_fine[0];
        y = prob_lo_fine[1] + hi_fine[1] * dx_fine[1];
        z = prob_lo_fine[2] + hi_fine[2] * dx_fine[2];

        auto idx_hi = find_bound_idx(x, y, z, bl_coarse, geom_coarse, BoundType::Hi);

        b.setSmall(idx_lo);
        b.setBig(idx_hi);

         /*Print() << "lo fine = " << lo_fine << std::endl;
         Print() << "hi fine = " << hi_fine << std::endl;
        Print() << " idx lo = " << idx_lo << std::endl;
        Print() << "idx_hi = " << idx_hi << std::endl;*/

    }

    BoxArray cba(std::move(bl_fine));
    cba.convert(IndexType::TheNodeType());  // <-- Make it nodal in all directions
    coarse_multifab_on_fine_dmap.define(cba, mf_cc_fine.DistributionMap(), mf_nc_coarse.nComp(), 0);
    coarse_multifab_on_fine_dmap.ParallelCopy(mf_nc_coarse);
}


void
PopulateFineCellCenteredFromCoarseNodal(const Geometry& geom_coarse,
                                        const Geometry& geom_fine,
                                        const MultiFab& coarse_multifab_on_fine_dmap,
                                        const MultiFab& mf_cc_fine,
                                        MultiFab& mf_cc_from_coarse)
{
    AMREX_ALWAYS_ASSERT(coarse_multifab_on_fine_dmap.ixType().nodeCentered());

    if (!mf_cc_from_coarse.isDefined())
    {
        mf_cc_from_coarse.define(mf_cc_fine.boxArray(),
                                 mf_cc_fine.DistributionMap(),
                                 coarse_multifab_on_fine_dmap.nComp(),
                                 0);
    }

    AMREX_ALWAYS_ASSERT(mf_cc_from_coarse.boxArray() == mf_cc_fine.boxArray());
    AMREX_ALWAYS_ASSERT(mf_cc_from_coarse.DistributionMap() == mf_cc_fine.DistributionMap());
    AMREX_ALWAYS_ASSERT(mf_cc_from_coarse.nComp() == coarse_multifab_on_fine_dmap.nComp());

    const auto prob_lo_coarse = geom_coarse.ProbLoArray();
    const auto dx_coarse = geom_coarse.CellSizeArray();

    const auto prob_lo_fine = geom_fine.ProbLoArray();
    const auto dx_fine = geom_fine.CellSizeArray();

    Box nodal_domain = amrex::convert(geom_coarse.Domain(), IntVect::TheNodeVector());
    const auto nd_lo = nodal_domain.smallEnd();
    const auto nd_hi = nodal_domain.bigEnd();

    int ncomp = mf_cc_from_coarse.nComp();

    for (MFIter mfi(mf_cc_from_coarse, TilingIfNotGPU()); mfi.isValid(); ++mfi)
    {
        const Box& cbox = mfi.tilebox();
        auto const& arr_cc = mf_cc_from_coarse.array(mfi);
        auto const& arr_nc = coarse_multifab_on_fine_dmap.array(mfi);

        ParallelFor(cbox, ncomp,
        [=] AMREX_GPU_DEVICE(int i, int j, int k, int n)
        {
            Real x = prob_lo_fine[0] + (static_cast<Real>(i) + 0.5_rt) * dx_fine[0];
            Real y = prob_lo_fine[1] + (static_cast<Real>(j) + 0.5_rt) * dx_fine[1];
            Real z = prob_lo_fine[2] + (static_cast<Real>(k) + 0.5_rt) * dx_fine[2];

            Real fx = (x - prob_lo_coarse[0]) / dx_coarse[0];
            Real fy = (y - prob_lo_coarse[1]) / dx_coarse[1];
            Real fz = (z - prob_lo_coarse[2]) / dx_coarse[2];

            int i0 = static_cast<int>(amrex::Math::floor(fx));
            int j0 = static_cast<int>(amrex::Math::floor(fy));
            int k0 = static_cast<int>(amrex::Math::floor(fz));

            Real wx = fx - static_cast<Real>(i0);
            Real wy = fy - static_cast<Real>(j0);
            Real wz = fz - static_cast<Real>(k0);

            int i1 = i0 + 1;
            int j1 = j0 + 1;
            int k1 = k0 + 1;

            i0 = amrex::max(nd_lo[0], amrex::min(i0, nd_hi[0]-1));
            j0 = amrex::max(nd_lo[1], amrex::min(j0, nd_hi[1]-1));
            k0 = amrex::max(nd_lo[2], amrex::min(k0, nd_hi[2]-1));

            i1 = amrex::max(nd_lo[0], amrex::min(i1, nd_hi[0]));
            j1 = amrex::max(nd_lo[1], amrex::min(j1, nd_hi[1]));
            k1 = amrex::max(nd_lo[2], amrex::min(k1, nd_hi[2]));

            Real c000 = arr_nc(i0,j0,k0,n);
            Real c100 = arr_nc(i1,j0,k0,n);
            Real c010 = arr_nc(i0,j1,k0,n);
            Real c110 = arr_nc(i1,j1,k0,n);
            Real c001 = arr_nc(i0,j0,k1,n);
            Real c101 = arr_nc(i1,j0,k1,n);
            Real c011 = arr_nc(i0,j1,k1,n);
            Real c111 = arr_nc(i1,j1,k1,n);

            Real c00 = c000 * (1.0_rt - wx) + c100 * wx;
            Real c10 = c010 * (1.0_rt - wx) + c110 * wx;
            Real c01 = c001 * (1.0_rt - wx) + c101 * wx;
            Real c11 = c011 * (1.0_rt - wx) + c111 * wx;

            Real c0 = c00 * (1.0_rt - wy) + c10 * wy;
            Real c1 = c01 * (1.0_rt - wy) + c11 * wy;

            arr_cc(i,j,k,n) = c0 * (1.0_rt - wz) + c1 * wz;
        });
    }
}

void WriteCustomDataFile(const Geometry& geom,
                         const MultiFab& mf_cc,
                         const std::string& filename_custom)
{
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
