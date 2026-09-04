/**
 * \file ERF_IBFaceSet.cpp
 * \brief Detection, storage, output and reporting of the wall faces of
 *        resolved buildings (phase 1 of the immersed-boundary surface energy
 *        balance). The conventions are documented in ERF_IBFaceSet.H.
 */
#include "ERF_IBFaceSet.H"
#include "ERF_IBSEBSolar.H"
#include <ERF_EOS.H>
#include <ERF_IndexDefines.H>
#include <AMReX_ParallelDescriptor.H>
#include <AMReX_Gpu.H>
#include <fstream>
#include <iomanip>
#include <sstream>
#include <algorithm>

using namespace amrex;

namespace {

/**
 * Host-readable view of one fab of a MultiFab.
 *
 * The face detection is a host loop at initialisation. On CPU builds the fab
 * is read in place; on GPU builds a pinned copy is made so the loop does not
 * touch device memory. The copy keeps the fab's own box, so the returned
 * Array4 is indexed with the same global (i, j, k) as the device fab.
 */
struct HostFab
{
#ifdef AMREX_USE_GPU
    FArrayBox fab;
    explicit HostFab (const FArrayBox& d) : fab(d.box(), d.nComp(), The_Pinned_Arena())
    {
        Gpu::dtoh_memcpy_async(fab.dataPtr(), d.dataPtr(), d.size() * sizeof(Real));
        Gpu::streamSynchronize();
    }
    Array4<const Real> array () const { return fab.const_array(); }
#else
    const FArrayBox* fab;
    explicit HostFab (const FArrayBox& d) : fab(&d) {}
    Array4<const Real> array () const { return fab->const_array(); }
#endif
};

/** Resize a device vector to a host vector and copy it over. */
template <class T>
void upload (Gpu::DeviceVector<T>& d, const std::vector<T>& h)
{
    d.resize(h.size());
    Gpu::copy(Gpu::hostToDevice, h.begin(), h.end(), d.begin());
}

/** Resize a device vector to n entries all equal to value. */
template <class T>
void fill (Gpu::DeviceVector<T>& d, size_t n, T value)
{
    std::vector<T> h(n, value);
    upload(d, h);
}

} // namespace

/**
 * Detect the faces owned by this rank and allocate the per-face arrays.
 *
 * The scan visits every valid cell of every local fab in MFIter order. A
 * fluid cell (blanking < 0.5) contributes one face for each of its six
 * neighbours that is solid (blanking >= 0.5); the neighbour may sit in a
 * ghost cell, which is why the blanking must have its ghost cells filled.
 * Neighbours outside the domain in a non-periodic direction are skipped, so
 * a building against a non-periodic boundary has no face there. Solid cells
 * contribute nothing but the column mask used for the building ids.
 *
 * Faces are appended in scan order, which makes them contiguous per fab;
 * m_fab_start records where each fab's faces begin. The per-direction,
 * per-building and total counts are reduced over all ranks here, once, so
 * report() does not have to reduce static numbers every time.
 */
void
IBFaceSet::build (const MultiFab& blanking, const Geometry& geom)
{
    const Box& domain = geom.Domain();
    const int nx  = domain.length(0);
    const int ny  = domain.length(1);
    const int ilo = domain.smallEnd(0);
    const int jlo = domain.smallEnd(1);
    const auto dx = geom.CellSizeArray();
    const auto plo = geom.ProbLoArray();
    const Real face_area[3] = { dx[1] * dx[2], dx[0] * dx[2], dx[0] * dx[1] };
    m_nx = nx; m_ny = ny;
    m_x_lo = plo[0]; m_y_lo = plo[1];
    for (int d = 0; d < 3; ++d) { m_dx[d] = dx[d]; }
    m_per_x = geom.isPeriodic(0); m_per_y = geom.isPeriodic(1);
    m_max_path = 4.0 * (geom.ProbLength(0) + geom.ProbLength(1));
    m_z_ground = plo[2];

    std::vector<int>  h_i, h_j, h_k, h_dir, h_side, h_nbi, h_nbj;
    std::vector<Real> h_area, h_xf, h_yf, h_zf;
    // Largest blanking in each column of this rank; reduced below so every
    // rank labels the same solid columns. The highest solid cell of each
    // column gives the column top for the ray cast.
    std::vector<Real> colmax(static_cast<size_t>(nx) * ny, 0.0);
    std::vector<Real> colk(static_cast<size_t>(nx) * ny, -1.0);

    m_fab_start.clear();
    for (MFIter mfi(blanking); mfi.isValid(); ++mfi) {
        m_fab_start.push_back(static_cast<int>(h_i.size()));
        const Box& bx = mfi.validbox();
        const HostFab hf(blanking[mfi]);
        auto const& b = hf.array();
        const Box& fbox = Box(b);
        const auto lo = lbound(bx);
        const auto hi = ubound(bx);
        for (int k = lo.z; k <= hi.z; ++k) {
        for (int j = lo.y; j <= hi.y; ++j) {
        for (int i = lo.x; i <= hi.x; ++i) {
            if (b(i, j, k) >= 0.5) {
                const size_t c = static_cast<size_t>(i - ilo) * ny + (j - jlo);
                colmax[c] = std::max(colmax[c], b(i, j, k));
                colk[c]   = std::max(colk[c], static_cast<Real>(k - domain.smallEnd(2)));
                continue;
            }
            for (int d = 0; d < 3; ++d) {
                for (int s = -1; s <= 1; s += 2) {
                    IntVect nb(i, j, k);
                    nb[d] += s;
                    if (!fbox.contains(nb)) { continue; }
                    if (!geom.isPeriodic(d) && !domain.contains(nb)) { continue; }
                    if (b(nb) < 0.5) { continue; }
                    h_i.push_back(i); h_j.push_back(j); h_k.push_back(k);
                    h_dir.push_back(d); h_side.push_back(s);
                    // Face centre: the fluid cell centre moved half a cell toward the solid.
                    Real c[3] = { plo[0] + (i + 0.5) * dx[0], plo[1] + (j + 0.5) * dx[1], plo[2] + (k + 0.5) * dx[2] };
                    c[d] += 0.5 * s * dx[d];
                    h_xf.push_back(c[0]); h_yf.push_back(c[1]); h_zf.push_back(c[2]);
                    // Column of the solid neighbour, wrapped for periodic directions.
                    int ci = nb[0] - ilo, cj = nb[1] - jlo;
                    ci = (ci % nx + nx) % nx;
                    cj = (cj % ny + ny) % ny;
                    h_nbi.push_back(ci); h_nbj.push_back(cj);
                    h_area.push_back(face_area[d]);
                }
            }
        }}}
    }
    m_fab_start.push_back(static_cast<int>(h_i.size()));
    m_nface = static_cast<int>(h_i.size());

    // Buildings: 4-connected solid columns, numbered in scan order. The column
    // mask is reduced so every rank labels the same columns with the same ids;
    // the labelling itself is a plain depth-first flood fill on the host.
    ParallelDescriptor::ReduceRealMax(colmax.data(), nx * ny);
    ParallelDescriptor::ReduceRealMax(colk.data(), nx * ny);
    // Column tops for the ray cast: one cell above the highest solid cell,
    // the ground where the column is fluid.
    std::vector<Real> h_col_top(static_cast<size_t>(nx) * ny, plo[2]);
    m_col_top_max = plo[2];
    for (size_t c = 0; c < h_col_top.size(); ++c) {
        if (colk[c] >= 0.0) { h_col_top[c] = plo[2] + (colk[c] + 1.0) * dx[2]; }
        m_col_top_max = std::max(m_col_top_max, h_col_top[c]);
    }
    upload(d_col_top, h_col_top);
    std::vector<int> label(static_cast<size_t>(nx) * ny, 0);
    std::vector<int> stack;
    m_nbld = 0;
    for (int ci = 0; ci < nx; ++ci) {
        for (int cj = 0; cj < ny; ++cj) {
            const size_t p = static_cast<size_t>(ci) * ny + cj;
            if (colmax[p] < 0.5 || label[p] > 0) { continue; }
            ++m_nbld;
            label[p] = m_nbld;
            stack.clear();
            stack.push_back(static_cast<int>(p));
            while (!stack.empty()) {
                const int q = stack.back(); stack.pop_back();
                const int qi = q / ny, qj = q % ny;
                const int ni[4] = {qi - 1, qi + 1, qi, qi};
                const int nj[4] = {qj, qj, qj - 1, qj + 1};
                for (int n = 0; n < 4; ++n) {
                    if (ni[n] < 0 || ni[n] >= nx || nj[n] < 0 || nj[n] >= ny) { continue; }
                    const size_t r = static_cast<size_t>(ni[n]) * ny + nj[n];
                    if (colmax[r] >= 0.5 && label[r] == 0) {
                        label[r] = m_nbld;
                        stack.push_back(static_cast<int>(r));
                    }
                }
            }
        }
    }
    std::vector<int> h_bid(m_nface);
    for (int n = 0; n < m_nface; ++n) {
        h_bid[n] = label[static_cast<size_t>(h_nbi[n]) * ny + h_nbj[n]];
    }
    // Footprints, for the debug summary: columns and bounding box per building.
    m_bld_ncol.assign(m_nbld + 1, 0);
    m_bld_ilo.assign(m_nbld + 1, nx); m_bld_ihi.assign(m_nbld + 1, -1);
    m_bld_jlo.assign(m_nbld + 1, ny); m_bld_jhi.assign(m_nbld + 1, -1);
    for (int ci = 0; ci < nx; ++ci) {
        for (int cj = 0; cj < ny; ++cj) {
            const int b = label[static_cast<size_t>(ci) * ny + cj];
            if (b == 0) { continue; }
            m_bld_ncol[b] += 1;
            m_bld_ilo[b] = std::min(m_bld_ilo[b], ci + ilo); m_bld_ihi[b] = std::max(m_bld_ihi[b], ci + ilo);
            m_bld_jlo[b] = std::min(m_bld_jlo[b], cj + jlo); m_bld_jhi[b] = std::max(m_bld_jhi[b], cj + jlo);
        }
    }

    // Static totals: per direction, per building, and the summed area.
    std::vector<long> nd(3, 0);
    std::vector<long> bn(m_nbld + 1, 0);
    std::vector<Real> ba(m_nbld + 1, 0.0);
    Real area = 0.0;
    for (int n = 0; n < m_nface; ++n) {
        nd[h_dir[n]] += 1;
        bn[h_bid[n]] += 1;
        ba[h_bid[n]] += h_area[n];
        area += h_area[n];
    }
    ParallelDescriptor::ReduceLongSum(nd.data(), 3);
    ParallelDescriptor::ReduceLongSum(bn.data(), m_nbld + 1);
    ParallelDescriptor::ReduceRealSum(ba.data(), m_nbld + 1);
    ParallelDescriptor::ReduceRealSum(area);
    for (int d = 0; d < 3; ++d) { m_nface_dir[d] = nd[d]; }
    m_bld_nface = bn;
    m_bld_area  = ba;
    m_area_total = area;

    // Device arrays. The state starts uniform at T_skin_init (the slab too:
    // phase 5 sets the interior boundary); the view fractions and fluxes start
    // at zero and are filled by the later phases.
    const size_t nf = static_cast<size_t>(m_nface);
    upload(d_i, h_i); upload(d_j, h_j); upload(d_k, h_k);
    upload(d_dir, h_dir); upload(d_side, h_side); upload(d_bid, h_bid);
    upload(d_area, h_area);
    upload(d_xf, h_xf); upload(d_yf, h_yf); upload(d_zf, h_zf);
    fill(d_mat,      nf, 0);
    fill(d_T_skin,   nf, m_params.T_skin_init);
    fill(d_T_slab,   nf * static_cast<size_t>(n_layers()), m_params.T_skin_init);
    // View fractions: placeholders (a roof sees the whole sky, a wall half sky
    // and half ground) that compute_view_fractions() replaces at initialisation.
    {
        std::vector<Real> fs(nf), fg(nf), fb(nf, 0.0);
        for (size_t n = 0; n < nf; ++n) {
            const bool roof = (h_dir[n] == 2);
            fs[n] = roof ? 1.0 : 0.5;
            fg[n] = roof ? 0.0 : 0.5;
        }
        upload(d_f_sky, fs); upload(d_f_ground, fg); upload(d_f_bldg, fb);
    }
    fill(d_shadow,        nf, Real(0.0));
    fill(d_SW_direct_in,  nf, Real(0.0));
    fill(d_SW_diffuse_in, nf, Real(0.0));
    fill(d_LW_down_in,    nf, Real(0.0));
    fill(d_T_air,         nf, m_params.T_skin_init);
    fill(d_SW_abs,   nf, Real(0.0));
    fill(d_LW_net,   nf, Real(0.0));
    fill(d_H,        nf, Real(0.0));
    fill(d_LE,       nf, Real(0.0));
    fill(d_G,        nf, Real(0.0));
    fill(d_Q_ext,    nf, Real(0.0));
    Gpu::streamSynchronize();

    if (m_params.debug) { print_debug_summary(); }
}

/**
 * Hemisphere sampling for the view fractions. One kernel over the faces,
 * each looping over its rays; the count of rays ending on the sky, the
 * ground and a building over the total gives the three fractions. Roofs
 * point up and never see the ground; a wall on flat open ground sees half
 * sky and half ground, which the regtest checks.
 */
void
IBFaceSet::compute_view_fractions ()
{
    const int n_az = m_params.view_n_az, n_el = m_params.view_n_el;
    const Real* col_top = d_col_top.data();
    const int   nx = m_nx, ny = m_ny;
    const Real  x_lo = m_x_lo, y_lo = m_y_lo, dx = m_dx[0], dy = m_dx[1];
    const bool  per_x = m_per_x, per_y = m_per_y;
    const Real  z_ground = m_z_ground, z_max = m_col_top_max, max_path = m_max_path;
    const int*  pd = d_dir.data();  const int* ps = d_side.data();
    const Real* pxf = d_xf.data();  const Real* pyf = d_yf.data();  const Real* pzf = d_zf.data();
    Real* pfs = d_f_sky.data(); Real* pfg = d_f_ground.data(); Real* pfb = d_f_bldg.data();
    ParallelFor(m_nface, [=] AMREX_GPU_DEVICE (int f) noexcept {
        int n_sky = 0, n_gnd = 0, n_bld = 0;
        for (int ie = 0; ie < n_el; ++ie) {
            for (int ia = 0; ia < n_az; ++ia) {
                Real sx, sy, sz;
                ibseb::hemisphere_direction(pd[f], -ps[f], ia, ie, n_az, n_el, sx, sy, sz);
                const int hit = ibseb::ray_hit(pxf[f], pyf[f], pzf[f], sx, sy, sz, col_top, nx, ny,
                                               x_lo, y_lo, dx, dy, per_x, per_y, z_ground, z_max, max_path);
                if (hit == ibseb::RAY_SKY) ++n_sky; else if (hit == ibseb::RAY_GROUND) ++n_gnd; else ++n_bld;
            }
        }
        const Real inv = 1.0 / static_cast<Real>(n_az * n_el);
        pfs[f] = n_sky * inv; pfg[f] = n_gnd * inv; pfb[f] = n_bld * inv;
    });
    Gpu::streamSynchronize();

    if (m_params.debug) {
        std::vector<Real> fs(m_nface), fg(m_nface), fb(m_nface);
        Gpu::copy(Gpu::deviceToHost, d_f_sky.begin(), d_f_sky.end(), fs.begin());
        Gpu::copy(Gpu::deviceToHost, d_f_ground.begin(), d_f_ground.end(), fg.begin());
        Gpu::copy(Gpu::deviceToHost, d_f_bldg.begin(), d_f_bldg.end(), fb.begin());
        Gpu::streamSynchronize();
        Real s[3] = {0.0, 0.0, 0.0};
        for (int n = 0; n < m_nface; ++n) { s[0] += fs[n]; s[1] += fg[n]; s[2] += fb[n]; }
        ParallelDescriptor::ReduceRealSum(s, 3);
        const Real nf_all = static_cast<Real>(m_nface_dir[0] + m_nface_dir[1] + m_nface_dir[2]);
        Print() << "[IBSEB DEBUG] lev=" << m_lev << " view fractions: " << n_az * n_el
                << " rays per face, mean f_sky=" << s[0] / nf_all << " f_ground=" << s[1] / nf_all
                << " f_bldg=" << s[2] / nf_all << "\n";
    }
}

/**
 * Longwave of the current step. One kernel per fab over that fab's faces,
 * reading the air temperature of the fluid cell from the conserved state
 * (potential temperature and density through the equation of state).
 */
void
IBFaceSet::compute_longwave (const MultiFab& cons)
{
    constexpr Real sigma = 5.670374419e-8;
    const bool  gray   = (m_params.lw_mode == "gray");
    const Real  lw_fix = m_params.lw_down, eps_sky = m_params.sky_emissivity;
    const Real  eps    = m_params.emissivity, eps_g = m_params.emissivity_ground;
    const Real  Tg4    = std::pow(m_params.T_ground, 4);
    const int*  pi = d_i.data();  const int* pj = d_j.data();  const int* pk = d_k.data();
    const Real* pfs = d_f_sky.data(); const Real* pfg = d_f_ground.data(); const Real* pfb = d_f_bldg.data();
    const Real* pT = d_T_skin.data();
    Real* pTa = d_T_air.data(); Real* pin = d_LW_down_in.data(); Real* pnet = d_LW_net.data();
    for (MFIter mfi(cons); mfi.isValid(); ++mfi) {
        const int f0 = m_fab_start[mfi.LocalIndex()];
        const int f1 = m_fab_start[mfi.LocalIndex() + 1];
        auto const& c = cons.const_array(mfi);
        ParallelFor(f1 - f0, [=] AMREX_GPU_DEVICE (int m) noexcept {
            const int f = f0 + m;
            const Real rho = c(pi[f], pj[f], pk[f], Rho_comp);
            const Real Ta  = getTgivenRandRTh(rho, c(pi[f], pj[f], pk[f], RhoTheta_comp));
            const Real lw_sky = gray ? eps_sky * sigma * Ta * Ta * Ta * Ta : lw_fix;
            const Real Ts4 = pT[f] * pT[f] * pT[f] * pT[f];
            const Real lw_in = pfs[f] * lw_sky + pfg[f] * eps_g * sigma * Tg4 + pfb[f] * sigma * Ts4;
            pTa[f]  = Ta;
            pin[f]  = lw_in;
            pnet[f] = eps * (lw_in - sigma * Ts4);
        });
    }
    Gpu::streamSynchronize();
}

/**
 * Shortwave of the current step. The sun and the irradiances are scalars
 * computed on the host from the provider inputs; the per-face work is one
 * kernel: incidence cosine, ray cast against the column tops, direct and
 * diffuse incident, absorbed.
 *
 * With sun_mode = fixed the direct-normal irradiance and the horizontal
 * diffuse are the inputs as given. With sun_mode = solar the sun follows the
 * site and time and the clear-sky formulas give both irradiances; the sun
 * below the horizon gives zero everywhere.
 *
 * Diffuse on a face: ``f_sky * diffuse_h + f_ground * albedo_ground *
 * (direct_h + diffuse_h)``, the second term being the ground-reflected part,
 * with the direct on a horizontal surface ``direct_h = dni * cos z``. The
 * view fractions are the placeholders of build() until phase 3.
 */
void
IBFaceSet::compute_shortwave (Real time)
{
    // ---- Sun and irradiances (host scalars) ----
    SunState s;
    if (m_params.sun_mode == "fixed") {
        s.zenith  = m_params.sun_zenith_deg  * M_PI / 180.0;
        s.azimuth = m_params.sun_azimuth_deg * M_PI / 180.0;
        const Real cz = std::cos(s.zenith);
        s.dni       = (cz > 0.0) ? m_params.sw_direct_normal : 0.0;
        s.diffuse_h = (cz > 0.0) ? m_params.sw_diffuse : 0.0;
    } else {
        const Real t_utc = m_params.time_zero_utc_s + time;
        const Real decl  = ibseb::solar_declination(m_params.day_of_year);
        const Real ha    = ibseb::solar_hour_angle(t_utc, m_params.longitude_deg,
                                                   m_params.utc_offset_hours, m_params.day_of_year);
        s.zenith  = ibseb::solar_zenith(m_params.latitude_deg, decl, ha);
        s.azimuth = ibseb::solar_azimuth(m_params.latitude_deg, decl, ha, s.zenith);
        const Real cz = std::cos(s.zenith);
        const Real e0 = ibseb::earth_sun_distance_factor(m_params.day_of_year);
        s.dni       = ibseb::clear_sky_dni(cz, m_params.solar_constant, m_params.sw_transmission, e0);
        s.diffuse_h = ibseb::clear_sky_diffuse_h(cz, m_params.solar_constant, m_params.sw_transmission,
                                                 e0, m_params.sw_diffuse_coeff);
    }
    ibseb::sun_vector(s.zenith, s.azimuth, s.sx, s.sy, s.sz);
    m_sun = s;

    // ---- Per face ----
    const Real sx = s.sx, sy = s.sy, sz = s.sz;
    const Real dni = s.dni, dif_h = s.diffuse_h;
    const Real dir_h = dni * amrex::max(Real(0.0), sz);
    const Real albedo = m_params.albedo, alb_g = m_params.albedo_ground;
    const Real* col_top = d_col_top.data();
    const int   nx = m_nx, ny = m_ny;
    const Real  x_lo = m_x_lo, y_lo = m_y_lo, dx = m_dx[0], dy = m_dx[1];
    const bool  per_x = m_per_x, per_y = m_per_y;
    const Real  z_max = m_col_top_max, max_path = m_max_path;
    const int*  pd  = d_dir.data();  const int* ps = d_side.data();
    const Real* pxf = d_xf.data();   const Real* pyf = d_yf.data();  const Real* pzf = d_zf.data();
    const Real* pfs = d_f_sky.data(); const Real* pfg = d_f_ground.data();
    Real* psh = d_shadow.data();  Real* pdir = d_SW_direct_in.data();
    Real* pdif = d_SW_diffuse_in.data();  Real* pabs = d_SW_abs.data();
    ParallelFor(m_nface, [=] AMREX_GPU_DEVICE (int f) noexcept {
        // Outward normal: opposite to the side the solid is on.
        Real n[3] = {0.0, 0.0, 0.0};
        n[pd[f]] = -static_cast<Real>(ps[f]);
        const Real cosi = n[0] * sx + n[1] * sy + n[2] * sz;
        Real shadow = 0.0, direct = 0.0;
        if (sz > 0.0 && cosi > 0.0 && dni > 0.0) {
            shadow = ibseb::ray_blocked(pxf[f], pyf[f], pzf[f], sx, sy, sz, col_top, nx, ny,
                                        x_lo, y_lo, dx, dy, per_x, per_y, z_max, max_path) ? 1.0 : 0.0;
            direct = dni * cosi * (1.0 - shadow);
        }
        const Real diffuse = pfs[f] * dif_h + pfg[f] * alb_g * (dir_h + dif_h);
        psh[f]  = shadow;
        pdir[f] = direct;
        pdif[f] = diffuse;
        pabs[f] = (1.0 - albedo) * (direct + diffuse);
    });
    Gpu::streamSynchronize();

    if (m_params.debug) {
        Print() << "[IBSEB DEBUG] lev=" << m_lev << " sun: zenith=" << s.zenith * 180.0 / M_PI
                << " deg azimuth=" << s.azimuth * 180.0 / M_PI << " deg s=(" << s.sx << "," << s.sy << "," << s.sz
                << ") DNI=" << s.dni << " W/m2 diffuse_h=" << s.diffuse_h << " W/m2\n";
    }
}

/**
 * Scatter one per-face array as a per-cell mean, the same accumulation as
 * scatter_diagnostics() with the face count as the divisor.
 */
void
IBFaceSet::scatter_field (const Gpu::DeviceVector<Real>& v, MultiFab& out) const
{
    MultiFab cnt(out.boxArray(), out.DistributionMap(), 1, 0);
    cnt.setVal(0.0);
    out.setVal(0.0);
    const int* pi = d_i.data();  const int* pj = d_j.data();  const int* pk = d_k.data();
    const Real* pv = v.data();
    for (MFIter mfi(out); mfi.isValid(); ++mfi) {
        const int f0 = m_fab_start[mfi.LocalIndex()];
        const int f1 = m_fab_start[mfi.LocalIndex() + 1];
        auto const& n = cnt.array(mfi);
        auto const& o = out.array(mfi);
        ParallelFor(f1 - f0, [=] AMREX_GPU_DEVICE (int m) noexcept {
            const int f = f0 + m;
            Gpu::Atomic::AddNoRet(&n(pi[f], pj[f], pk[f]), Real(1.0));
            Gpu::Atomic::AddNoRet(&o(pi[f], pj[f], pk[f]), pv[f]);
        });
        ParallelFor(mfi.validbox(), [=] AMREX_GPU_DEVICE (int i, int j, int k) noexcept {
            if (n(i, j, k) > 0.0) { o(i, j, k) /= n(i, j, k); }
        });
    }
    Gpu::streamSynchronize();
}

/**
 * Per-rank face dump. Every rank writes its own file so no gather is
 * needed; readers concatenate ``<prefix>.rank*.csv``.
 */
void
IBFaceSet::dump_faces (const std::string& prefix) const
{
    const int nf = m_nface;
    auto host = [&](const Gpu::DeviceVector<Real>& d) {
        std::vector<Real> h(nf);
        Gpu::copy(Gpu::deviceToHost, d.begin(), d.end(), h.begin());
        return h;
    };
    auto hosti = [&](const Gpu::DeviceVector<int>& d) {
        std::vector<int> h(nf);
        Gpu::copy(Gpu::deviceToHost, d.begin(), d.end(), h.begin());
        return h;
    };
    const auto i = hosti(d_i), j = hosti(d_j), k = hosti(d_k), dir = hosti(d_dir), side = hosti(d_side), bid = hosti(d_bid);
    const auto xf = host(d_xf), yf = host(d_yf), zf = host(d_zf), area = host(d_area);
    const auto fs = host(d_f_sky), fg = host(d_f_ground), fb = host(d_f_bldg);
    const auto sh = host(d_shadow), sdir = host(d_SW_direct_in), sdif = host(d_SW_diffuse_in), sabs = host(d_SW_abs);
    const auto T = host(d_T_skin), Ta = host(d_T_air), lwin = host(d_LW_down_in), lwn = host(d_LW_net);
    Gpu::streamSynchronize();
    std::ofstream f(prefix + ".rank" + std::to_string(ParallelDescriptor::MyProc()) + ".csv");
    f << "i,j,k,dir,side,bid,x_m,y_m,z_m,area_m2,f_sky,f_ground,f_bldg,shadow,SW_direct_in,SW_diffuse_in,SW_abs,T_skin,T_air,LW_in,LW_net\n";
    f << std::setprecision(10);
    for (int n = 0; n < nf; ++n) {
        f << i[n] << "," << j[n] << "," << k[n] << "," << dir[n] << "," << side[n] << "," << bid[n] << ","
          << xf[n] << "," << yf[n] << "," << zf[n] << "," << area[n] << ","
          << fs[n] << "," << fg[n] << "," << fb[n] << "," << sh[n] << ","
          << sdir[n] << "," << sdif[n] << "," << sabs[n] << "," << T[n] << ","
          << Ta[n] << "," << lwin[n] << "," << lwn[n] << "\n";
    }
}

/**
 * Debug description of the set, in the style of the fire module's
 * ``[FIRE DEBUG]`` lines: what was read, what every rank holds, what the
 * buildings look like, and what the arrays cost. Per-rank lines use
 * AllPrint() so every rank reports; the rest is printed by the I/O rank.
 */
void
IBFaceSet::print_debug_summary () const
{
    const int nl = n_layers();
    Print() << "[IBSEB DEBUG] lev=" << m_lev << " inputs: n_slab_layers=" << nl
            << " T_skin_init=" << m_params.T_skin_init << " K"
            << " T_interior=" << m_params.T_interior << " K"
            << " csv_file=" << m_params.csv_file << " csv_int=" << m_params.csv_int << "\n";
    Print() << "[IBSEB DEBUG] lev=" << m_lev << " global faces=" << (m_nface_dir[0] + m_nface_dir[1] + m_nface_dir[2])
            << " (x=" << m_nface_dir[0] << " y=" << m_nface_dir[1] << " z=" << m_nface_dir[2] << ")"
            << " area=" << m_area_total << " m2 buildings=" << m_nbld << "\n";
    for (int b = 1; b <= m_nbld; ++b) {
        Print() << "[IBSEB DEBUG]   building " << b << ": columns=" << m_bld_ncol[b]
                << " i=[" << m_bld_ilo[b] << "," << m_bld_ihi[b] << "]"
                << " j=[" << m_bld_jlo[b] << "," << m_bld_jhi[b] << "]"
                << " faces=" << m_bld_nface[b] << " area=" << m_bld_area[b] << " m2\n";
    }
    // Per-rank ownership: one line per rank with its fab ranges.
    const int nfab = static_cast<int>(m_fab_start.size()) - 1;
    const size_t bytes = static_cast<size_t>(m_nface) *
        (7 * sizeof(int) + (11 + static_cast<size_t>(nl)) * sizeof(Real));
    std::ostringstream os;
    os << "[IBSEB DEBUG] lev=" << m_lev << " rank " << ParallelDescriptor::MyProc()
       << ": faces=" << m_nface << " over " << nfab << " fab(s), device memory "
       << bytes / 1024.0 << " kB; fab ranges:";
    for (int f = 0; f < nfab; ++f) { os << " [" << m_fab_start[f] << "," << m_fab_start[f + 1] << ")"; }
    AllPrint() << os.str() << "\n";
}

/**
 * Scatter the face count and the mean skin temperature into cell-centred
 * fields. One kernel per fab over that fab's faces (contiguous, see
 * fab_start()), with atomic adds since a corner cell collects several faces,
 * then a kernel over the box divides the temperature sum by the count.
 */
void
IBFaceSet::scatter_diagnostics (MultiFab& nfaces, MultiFab& tskin) const
{
    nfaces.setVal(0.0);
    tskin.setVal(0.0);
    const int* pi = d_i.data();  const int* pj = d_j.data();  const int* pk = d_k.data();
    const Real* pT = d_T_skin.data();
    for (MFIter mfi(nfaces); mfi.isValid(); ++mfi) {
        const int f0 = m_fab_start[mfi.LocalIndex()];
        const int f1 = m_fab_start[mfi.LocalIndex() + 1];
        auto const& n = nfaces.array(mfi);
        auto const& t = tskin.array(mfi);
        ParallelFor(f1 - f0, [=] AMREX_GPU_DEVICE (int m) noexcept {
            const int f = f0 + m;
            Gpu::Atomic::AddNoRet(&n(pi[f], pj[f], pk[f]), Real(1.0));
            Gpu::Atomic::AddNoRet(&t(pi[f], pj[f], pk[f]), pT[f]);
        });
        const Box& bx = mfi.validbox();
        ParallelFor(bx, [=] AMREX_GPU_DEVICE (int i, int j, int k) noexcept {
            if (n(i, j, k) > 0.0) { t(i, j, k) /= n(i, j, k); }
        });
    }
    Gpu::streamSynchronize();
}

/**
 * Write the skin and slab temperatures into the six-slot state field. The
 * slot of a face is ``(dir*2 + (side>0)) * (1 + n_layers)``; no two faces of a
 * cell share a slot, so the writes are plain stores.
 */
void
IBFaceSet::save_state (MultiFab& state) const
{
    AMREX_ALWAYS_ASSERT(state.nComp() == state_ncomp());
    state.setVal(0.0);
    const int nl = n_layers();
    const int* pi = d_i.data();  const int* pj = d_j.data();  const int* pk = d_k.data();
    const int* pd = d_dir.data(); const int* ps = d_side.data();
    const Real* pT = d_T_skin.data(); const Real* pS = d_T_slab.data();
    for (MFIter mfi(state); mfi.isValid(); ++mfi) {
        const int f0 = m_fab_start[mfi.LocalIndex()];
        const int f1 = m_fab_start[mfi.LocalIndex() + 1];
        auto const& st = state.array(mfi);
        ParallelFor(f1 - f0, [=] AMREX_GPU_DEVICE (int m) noexcept {
            const int f  = f0 + m;
            const int c0 = (pd[f] * 2 + (ps[f] > 0 ? 1 : 0)) * (1 + nl);
            st(pi[f], pj[f], pk[f], c0) = pT[f];
            for (int l = 0; l < nl; ++l) { st(pi[f], pj[f], pk[f], c0 + 1 + l) = pS[f * nl + l]; }
        });
    }
    Gpu::streamSynchronize();
    if (m_params.debug) {
        Print() << "[IBSEB DEBUG] lev=" << m_lev << " face state saved: " << state_ncomp()
                << " components (6 slots x " << (1 + n_layers()) << ")\n";
    }
}

/**
 * Read the skin and slab temperatures back from the six-slot state field,
 * the inverse of save_state(), into a list that build() has just recreated.
 */
void
IBFaceSet::load_state (const MultiFab& state)
{
    AMREX_ALWAYS_ASSERT(state.nComp() == state_ncomp());
    const int nl = n_layers();
    const int* pi = d_i.data();  const int* pj = d_j.data();  const int* pk = d_k.data();
    const int* pd = d_dir.data(); const int* ps = d_side.data();
    Real* pT = d_T_skin.data(); Real* pS = d_T_slab.data();
    for (MFIter mfi(state); mfi.isValid(); ++mfi) {
        const int f0 = m_fab_start[mfi.LocalIndex()];
        const int f1 = m_fab_start[mfi.LocalIndex() + 1];
        auto const& st = state.const_array(mfi);
        ParallelFor(f1 - f0, [=] AMREX_GPU_DEVICE (int m) noexcept {
            const int f  = f0 + m;
            const int c0 = (pd[f] * 2 + (ps[f] > 0 ? 1 : 0)) * (1 + nl);
            pT[f] = st(pi[f], pj[f], pk[f], c0);
            for (int l = 0; l < nl; ++l) { pS[f * nl + l] = st(pi[f], pj[f], pk[f], c0 + 1 + l); }
        });
    }
    Gpu::streamSynchronize();
    if (m_params.debug) {
        std::vector<Real> h_T(m_nface);
        Gpu::copy(Gpu::deviceToHost, d_T_skin.begin(), d_T_skin.end(), h_T.begin());
        Gpu::streamSynchronize();
        Real tmin = 1.0e30, tmax = -1.0e30;
        for (int n = 0; n < m_nface; ++n) { tmin = std::min(tmin, h_T[n]); tmax = std::max(tmax, h_T[n]); }
        ParallelDescriptor::ReduceRealMin(tmin);
        ParallelDescriptor::ReduceRealMax(tmax);
        Print() << "[IBSEB DEBUG] lev=" << m_lev << " face state loaded: T_skin in ["
                << tmin << ", " << tmax << "] K\n";
    }
}

/**
 * Summary line and per-building CSV rows. The dynamic quantities (skin
 * temperature range and area-weighted mean per building) are gathered on the
 * host from the device arrays and reduced; the static counts and areas were
 * reduced in build().
 */
void
IBFaceSet::report (Real time, int step, bool write_csv) const
{
    // Skin temperature: min, max and per-building mean over all ranks.
    std::vector<Real> h_T(m_nface), h_A(m_nface);
    std::vector<int>  h_b(m_nface);
    Gpu::copy(Gpu::deviceToHost, d_T_skin.begin(), d_T_skin.end(), h_T.begin());
    Gpu::copy(Gpu::deviceToHost, d_area.begin(),   d_area.end(),   h_A.begin());
    Gpu::copy(Gpu::deviceToHost, d_bid.begin(),    d_bid.end(),    h_b.begin());
    Gpu::streamSynchronize();

    std::vector<Real> h_S(m_nface), h_sh(m_nface), h_L(m_nface);
    Gpu::copy(Gpu::deviceToHost, d_SW_abs.begin(), d_SW_abs.end(), h_S.begin());
    Gpu::copy(Gpu::deviceToHost, d_shadow.begin(), d_shadow.end(), h_sh.begin());
    Gpu::copy(Gpu::deviceToHost, d_LW_net.begin(), d_LW_net.end(), h_L.begin());
    Gpu::streamSynchronize();

    std::vector<Real> bsum(m_nbld + 1, 0.0), bsw(m_nbld + 1, 0.0), bsh(m_nbld + 1, 0.0), blw(m_nbld + 1, 0.0);
    Real tmin = 1.0e30, tmax = -1.0e30, sw_sum = 0.0, sw_max = 0.0, sh_area = 0.0, lw_sum = 0.0;
    for (int n = 0; n < m_nface; ++n) {
        bsum[h_b[n]] += h_T[n] * h_A[n];
        bsw[h_b[n]]  += h_S[n] * h_A[n];
        bsh[h_b[n]]  += h_sh[n] * h_A[n];
        blw[h_b[n]]  += h_L[n] * h_A[n];
        tmin = std::min(tmin, h_T[n]);
        tmax = std::max(tmax, h_T[n]);
        sw_sum  += h_S[n] * h_A[n];
        sw_max   = std::max(sw_max, h_S[n]);
        sh_area += h_sh[n] * h_A[n];
        lw_sum  += h_L[n] * h_A[n];
    }
    ParallelDescriptor::ReduceRealSum(bsum.data(), m_nbld + 1);
    ParallelDescriptor::ReduceRealSum(bsw.data(),  m_nbld + 1);
    ParallelDescriptor::ReduceRealSum(bsh.data(),  m_nbld + 1);
    ParallelDescriptor::ReduceRealSum(blw.data(),  m_nbld + 1);
    ParallelDescriptor::ReduceRealMin(tmin);
    ParallelDescriptor::ReduceRealMax(tmax);
    ParallelDescriptor::ReduceRealSum(sw_sum);
    ParallelDescriptor::ReduceRealMax(sw_max);
    ParallelDescriptor::ReduceRealSum(sh_area);
    ParallelDescriptor::ReduceRealSum(lw_sum);
    const Real sw_mean  = (m_area_total > 0.0) ? sw_sum / m_area_total : 0.0;
    const Real sh_frac  = (m_area_total > 0.0) ? sh_area / m_area_total : 0.0;
    const Real lw_mean  = (m_area_total > 0.0) ? lw_sum / m_area_total : 0.0;

    const long nf_all = m_nface_dir[0] + m_nface_dir[1] + m_nface_dir[2];
    Print() << "[IBSEB] lev=" << m_lev << " step=" << step << " t=" << time
            << " faces=" << nf_all
            << " (x=" << m_nface_dir[0] << " y=" << m_nface_dir[1] << " z=" << m_nface_dir[2] << ")"
            << " buildings=" << m_nbld
            << " area=" << std::setprecision(10) << m_area_total << " m2"
            << " T_skin_min=" << ((nf_all > 0) ? tmin : Real(0.0))
            << " T_skin_max=" << ((nf_all > 0) ? tmax : Real(0.0))
            << " SW_abs_mean=" << sw_mean << " SW_abs_max=" << sw_max
            << " shadow_frac=" << sh_frac << " LW_net_mean=" << lw_mean << "\n";

    if (m_params.debug) {
        for (int b = 1; b <= m_nbld; ++b) {
            const Real tmean = (m_bld_area[b] > 0.0) ? bsum[b] / m_bld_area[b] : Real(0.0);
            const Real swm = (m_bld_area[b] > 0.0) ? bsw[b] / m_bld_area[b] : Real(0.0);
            const Real shf = (m_bld_area[b] > 0.0) ? bsh[b] / m_bld_area[b] : Real(0.0);
            const Real lwm = (m_bld_area[b] > 0.0) ? blw[b] / m_bld_area[b] : Real(0.0);
            Print() << "[IBSEB DEBUG]   building " << b << ": faces=" << m_bld_nface[b]
                    << " area=" << m_bld_area[b] << " m2 T_skin_mean=" << tmean << " K"
                    << " SW_abs_mean=" << swm << " W/m2 shadow_frac=" << shf
                    << " LW_net_mean=" << lwm << " W/m2\n";
        }
    }

    if (write_csv && !m_params.dump_faces_file.empty()) { dump_faces(m_params.dump_faces_file); }

    if (!write_csv || !ParallelDescriptor::IOProcessor()) { return; }
    bool need_header = true;
    {
        std::ifstream probe(m_params.csv_file, std::ios::ate);
        if (probe.good() && probe.tellg() > 0) { need_header = false; }
    }
    std::ofstream csv(m_params.csv_file, std::ios::app);
    if (need_header) {
        csv << "time_s,step,level,building,n_faces,area_m2,T_skin_mean_K,SW_abs_mean_Wm2,shadow_frac,LW_net_mean_Wm2\n";
    }
    for (int b = 1; b <= m_nbld; ++b) {
        const Real tmean = (m_bld_area[b] > 0.0) ? bsum[b] / m_bld_area[b] : Real(0.0);
        const Real swm   = (m_bld_area[b] > 0.0) ? bsw[b]  / m_bld_area[b] : Real(0.0);
        const Real shf   = (m_bld_area[b] > 0.0) ? bsh[b]  / m_bld_area[b] : Real(0.0);
        const Real lwm   = (m_bld_area[b] > 0.0) ? blw[b]  / m_bld_area[b] : Real(0.0);
        csv << std::setprecision(10) << time << "," << step << "," << m_lev << "," << b << ","
            << m_bld_nface[b] << "," << m_bld_area[b] << "," << tmean << "," << swm << "," << shf << "," << lwm << "\n";
    }
}
