#include <AMReX_BoxArray.H>
#include <AMReX_DistributionMapping.H>
#include <AMReX_Geometry.H>
#include <AMReX_Interpolater.H>
#include <AMReX_iMultiFab.H>

#include <gtest/gtest.h>

#include "ERF_FillPatcher.H"

// Motivation (erf-model/ERF#4303): the relax band of ERFFillPatcher (cf_width > 0) was built by
// growing the complement of the fine boxes in x and y only, so a level that stops below the top of
// the domain got no band along its top coarse/fine face. The band must follow every coarse/fine
// face, the top one included, and must not appear along the physical top or bottom of the domain.

using namespace amrex;

namespace {

constexpr int ncx = 8, ncy = 2, ncz = 8;    // coarse cells; the fine level has ratio 2
constexpr int fi_lo = 4, fi_hi = 11;        // fine patch in x (interior coarse/fine faces)
constexpr int cf_width = 2;                 // width of the relax band, in fine cells

// The fine patch spans the domain in y (periodic) and runs from the bottom to fine cell fk_hi
iMultiFab relax_mask (int fk_hi, IntVect const& ixtype, int& set_val, int& relax_val)
{
    RealBox rb({0., 0., 0.}, {Real(ncx), Real(ncy), Real(ncz)});
    Array<int,3> per{0, 1, 0};
    Box cdom(IntVect(0), IntVect(ncx-1, ncy-1, ncz-1));
    Geometry cgeom(cdom, rb, 0, per);
    Geometry fgeom(amrex::refine(cdom, 2), rb, 0, per);

    BoxArray cba = convert(BoxArray(cdom), ixtype);
    BoxArray fba = convert(BoxArray(Box(IntVect(fi_lo, 0, 0), IntVect(fi_hi, 2*ncy-1, fk_hi))), ixtype);
    DistributionMapping cdm(cba), fdm(fba);

    InterpBase* interp = (ixtype == IntVect(0)) ? static_cast<InterpBase*>(&cell_cons_interp)
                                                : static_cast<InterpBase*>(&face_cons_linear_interp);
    ERFFillPatcher fp(fba, fdm, fgeom, cba, cdm, cgeom, -cf_width, 0, 1, interp);
    set_val   = fp.GetSetMaskVal();
    relax_val = fp.GetRelaxMaskVal();

    iMultiFab host(fba, fdm, 1, 0, MFInfo().SetArena(The_Pinned_Arena()));
    iMultiFab::Copy(host, *fp.GetMask(), 0, 0, 1, 0);
    Gpu::streamSynchronize();
    return host;
}

int at (iMultiFab const& host, int i, int j, int k)
{
    for (MFIter mfi(host); mfi.isValid(); ++mfi) {
        if (mfi.validbox().contains(IntVect(i, j, k))) { return host.const_array(mfi)(i, j, k); }
    }
    return -1;
}

bool in_x_band (int i) { return i < fi_lo + cf_width || i > fi_hi - cf_width; }

} // namespace

TEST(FillPatcherRelaxBand, CoversTheTopFaceOfALevelBelowTheDomainTop)
{
    constexpr int fk_hi = 9;   // the fine level stops 6 fine cells below the domain top
    int set_val = 0, relax_val = 0;
    const iMultiFab m = relax_mask(fk_hi, IntVect(0), set_val, relax_val);
    for (int k = 0; k <= fk_hi; ++k) {
        for (int i = fi_lo; i <= fi_hi; ++i) {
            // the band along the lateral faces and the top face; none along the domain bottom
            const bool band = in_x_band(i) || k > fk_hi - cf_width;
            EXPECT_EQ(at(m, i, 1, k), band ? relax_val : 0) << "i=" << i << " k=" << k;
        }
    }
}

TEST(FillPatcherRelaxBand, NoBandAlongThePhysicalTop)
{
    constexpr int fk_hi = 2*ncz - 1;   // the fine level reaches the domain top
    int set_val = 0, relax_val = 0;
    const iMultiFab m = relax_mask(fk_hi, IntVect(0), set_val, relax_val);
    for (int k = 0; k <= fk_hi; ++k) {
        for (int i = fi_lo; i <= fi_hi; ++i) {
            EXPECT_EQ(at(m, i, 1, k), in_x_band(i) ? relax_val : 0) << "i=" << i << " k=" << k;
        }
    }
}

TEST(FillPatcherRelaxBand, TopCoarseFineFaceOfZFacesKeepsTheSetRegion)
{
    constexpr int fk_hi = 9;   // z faces 0..fk_hi+1; the top one is the coarse/fine face
    int set_val = 0, relax_val = 0;
    const iMultiFab m = relax_mask(fk_hi, IntVect(0,0,1), set_val, relax_val);
    const int i = (fi_lo + fi_hi) / 2;   // away from the lateral band
    // w is set on the coarse/fine face itself and relaxed on the other faces of the cf_width
    // cells below it; the domain bottom has no band
    EXPECT_EQ(at(m, i, 1, fk_hi+1), set_val);
    for (int k = 0; k <= fk_hi; ++k) {
        const bool band = k >= fk_hi + 1 - cf_width;
        EXPECT_EQ(at(m, i, 1, k), band ? relax_val : 0) << "k=" << k;
    }
}
