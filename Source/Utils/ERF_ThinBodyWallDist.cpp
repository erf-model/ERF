#include <AMReX.H>
#include <AMReX_MultiFab.H>

using namespace amrex;

/*
 * Brute force approach to computing wall distance for thin immersed bodies.
 * This should be called after poisson_wall_dist, which will calculate the
 * wall distance relative to terrain.
 */
void
thinbody_wall_dist (std::unique_ptr<MultiFab>& wdist,
                    Vector<IntVect>& xfaces,
                    Vector<IntVect>& yfaces,
                    Vector<IntVect>& zfaces,
                    const Geometry& geomdata,
                    std::unique_ptr<MultiFab>& z_phys_cc)
{
    BL_PROFILE("thinbody_wall_dist()");

    const Real* prob_lo = geomdata.ProbLo();
    const Real* dx = geomdata.CellSize();

    const bool use_terrain = (z_phys_cc != nullptr);
    if (use_terrain) {
        Error("Thinbody wall dist calc not implemented for terrain yet");
    }

    for (MFIter mfi(*wdist); mfi.isValid(); ++mfi) {
        const Box& bx = mfi.validbox();

        const auto& z_cc = (use_terrain) ? z_phys_cc->const_array(mfi) : Array4<const Real>{};
        auto wd_arr      = wdist->array(mfi);

        if (!use_terrain) {
            ParallelFor(bx, [=] AMREX_GPU_DEVICE(int i, int j, int k) {
                Real xr = prob_lo[0] + (i + 0.5) * dx[0];
                Real yr = prob_lo[1] + (j + 0.5) * dx[1];
                Real zr = prob_lo[2] + (k + 0.5) * dx[2];

                for (int iface=0; iface < xfaces.size(); ++iface) {
                    int ii = xfaces[iface][0];
                    int jj = xfaces[iface][1];
                    int kk = xfaces[iface][2];
                    Real xfc = prob_lo[0] +  ii      * dx[0];
                    Real yfc = prob_lo[1] + (jj+0.5) * dx[1];
                    Real zfc = prob_lo[2] + (kk+0.5) * dx[2];
                    Real y0  = prob_lo[1] +  jj   * dx[1];
                    Real y1  = prob_lo[1] + (jj+1)* dx[1];
                    Real z0  = prob_lo[2] +  kk   * dx[2];
                    Real z1  = prob_lo[2] + (kk+1)* dx[2];
                    Real wd2 = wd_arr(i, j, k) * wd_arr(i, j, k);
                    wd2 = min(wd2, (xfc-xr)*(xfc-xr) + (yfc-yr)*(yfc-yr) + (zfc-zr)*(zfc-zr));
                    wd2 = min(wd2, (xfc-xr)*(xfc-xr) + ( y0-yr)*( y0-yr) + ( z0-zr)*( z0-zr));
                    wd2 = min(wd2, (xfc-xr)*(xfc-xr) + ( y0-yr)*( y0-yr) + ( z1-zr)*( z1-zr));
                    wd2 = min(wd2, (xfc-xr)*(xfc-xr) + ( y1-yr)*( y1-yr) + ( z0-zr)*( z0-zr));
                    wd2 = min(wd2, (xfc-xr)*(xfc-xr) + ( y1-yr)*( y1-yr) + ( z1-zr)*( z1-zr));
                    wd_arr(i, j, k) = std::sqrt(wd2);
                }

                for (int iface=0; iface < yfaces.size(); ++iface) {
                    int ii = yfaces[iface][0];
                    int jj = yfaces[iface][1];
                    int kk = yfaces[iface][2];
                    Real xfc = prob_lo[0] + (ii+0.5) * dx[0];
                    Real yfc = prob_lo[1] +  jj      * dx[1];
                    Real zfc = prob_lo[2] + (kk+0.5) * dx[2];
                    Real x0  = prob_lo[0] +  ii   * dx[0];
                    Real x1  = prob_lo[0] + (ii+1)* dx[0];
                    Real z0  = prob_lo[2] +  kk   * dx[2];
                    Real z1  = prob_lo[2] + (kk+1)* dx[2];
                    Real wd2 = wd_arr(i, j, k) * wd_arr(i, j, k);
                    wd2 = min(wd2, (xfc-xr)*(xfc-xr) + (yfc-yr)*(yfc-yr) + (zfc-zr)*(zfc-zr));
                    wd2 = min(wd2, ( x0-xr)*( x0-xr) + (yfc-yr)*(yfc-yr) + ( z0-zr)*( z0-zr));
                    wd2 = min(wd2, ( x0-xr)*( x0-xr) + (yfc-yr)*(yfc-yr) + ( z1-zr)*( z1-zr));
                    wd2 = min(wd2, ( x1-xr)*( x1-xr) + (yfc-yr)*(yfc-yr) + ( z0-zr)*( z0-zr));
                    wd2 = min(wd2, ( x1-xr)*( x1-xr) + (yfc-yr)*(yfc-yr) + ( z1-zr)*( z1-zr));
                    wd_arr(i, j, k) = std::sqrt(wd2);
                }

                for (int iface=0; iface < zfaces.size(); ++iface) {
                    int ii = zfaces[iface][0];
                    int jj = zfaces[iface][1];
                    int kk = zfaces[iface][2];
                    Real xfc = prob_lo[0] + (ii+0.5) * dx[0];
                    Real yfc = prob_lo[1] + (jj+0.5) * dx[1];
                    Real zfc = prob_lo[2] +  kk      * dx[2];
                    Real x0  = prob_lo[0] +  ii   * dx[0];
                    Real x1  = prob_lo[0] + (ii+1)* dx[0];
                    Real y0  = prob_lo[1] +  jj   * dx[1];
                    Real y1  = prob_lo[1] + (jj+1)* dx[1];
                    Real wd2 = wd_arr(i, j, k) * wd_arr(i, j, k);
                    wd2 = min(wd2, (xfc-xr)*(xfc-xr) + (yfc-yr)*(yfc-yr) + (zfc-zr)*(zfc-zr));
                    wd2 = min(wd2, ( x0-xr)*( x0-xr) + ( y0-yr)*( y0-yr) + (zfc-zr)*(zfc-zr));
                    wd2 = min(wd2, ( x0-xr)*( x0-xr) + ( y0-yr)*( y0-yr) + (zfc-zr)*(zfc-zr));
                    wd2 = min(wd2, ( x1-xr)*( x1-xr) + ( y1-yr)*( y1-yr) + (zfc-zr)*(zfc-zr));
                    wd2 = min(wd2, ( x1-xr)*( x1-xr) + ( y1-yr)*( y1-yr) + (zfc-zr)*(zfc-zr));
                    wd_arr(i, j, k) = std::sqrt(wd2);
                }
            });
        }
    }
}
