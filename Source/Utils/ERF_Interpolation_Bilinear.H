#ifndef ERF_INTERPOLATE_BILINEAR_H_
#define ERF_INTERPOLATE_BILINEAR_H_

#include "AMReX_GpuPrint.H"
#include "ERF_DataStruct.H"

/**
 * Bilinear interpolation routines. There are different versions for
 * different scenarios that take in a different argument list
 */

AMREX_FORCE_INLINE
AMREX_GPU_HOST_DEVICE
int get_single_index(int i, int j, int k,
             int nx, int ny)
{

    int si = k*nx*ny + j*nx + i;
    return si;
}

AMREX_FORCE_INLINE
AMREX_GPU_HOST_DEVICE
int clamp_interp_index (int idx, int npts)
{
    if (npts <= 1) {
        return 0;
    }

    return amrex::max(0, amrex::min(idx, npts - 2));
}

AMREX_FORCE_INLINE
AMREX_GPU_HOST_DEVICE
amrex::Real get_interp_weight (amrex::Real x, amrex::Real xlo, amrex::Real xhi)
{
    return (xhi > xlo) ? (x - xlo) / (xhi - xlo) : amrex::Real(0);
}

AMREX_FORCE_INLINE
AMREX_GPU_HOST_DEVICE
void bilinear_interpolation(const amrex::Real* xvec, const amrex::Real* yvec, const amrex::Real* zvec,
                const amrex::Real dxvec, const amrex::Real dyvec,
                const int nx, const int ny, const int nz,
                const amrex::Real x, const amrex::Real y, const amrex::Real z,
                const amrex::Real* varvec,
                amrex::Real& tmp_var)
{

    if(z < amrex::Real(0)){
        tmp_var = amrex::Real(0);
        return;
    }

    amrex::Real x_interp = x;
    amrex::Real y_interp = y;
    amrex::Real z_interp = z;

    if (x_interp < xvec[0]) x_interp = xvec[0];
    if (x_interp > xvec[nx-1]) x_interp = xvec[nx-1];

    if (y_interp < yvec[0]) y_interp = yvec[0];
    if (y_interp > yvec[ny-1]) y_interp = yvec[ny-1];

    if (z_interp < zvec[0]) z_interp = zvec[0];
    if (z_interp > zvec[nz-1]) z_interp = zvec[nz-1];

    int iloc = clamp_interp_index(static_cast<int>(std::floor((x_interp - xvec[0]) / dxvec)), nx);
    int jloc = clamp_interp_index(static_cast<int>(std::floor((y_interp - yvec[0]) / dyvec)), ny);

    int kloc = 0;
    if (nz > 1) {
        kloc = nz - 2;
        for (int k = 1; k < nz; ++k) {
            if (zvec[k] >= z_interp) {
                kloc = k - 1;
                break;
            }
        }
    }

    int ihi = (nx > 1) ? iloc + 1 : iloc;
    int jhi = (ny > 1) ? jloc + 1 : jloc;
    int khi = (nz > 1) ? kloc + 1 : kloc;

    amrex::Real xlo = xvec[iloc];
    amrex::Real ylo = yvec[jloc];
    amrex::Real zlo = zvec[kloc];

    amrex::Real w_x = get_interp_weight(x_interp, xlo, xvec[ihi]);
    amrex::Real w_y = get_interp_weight(y_interp, ylo, yvec[jhi]);
    amrex::Real w_z = get_interp_weight(z_interp, zlo, zvec[khi]);

    int ind0 = get_single_index(iloc,jloc,kloc,nx,ny);
    int ind1 = get_single_index(ihi ,jloc,kloc,nx,ny);
    int ind2 = get_single_index(ihi ,jhi ,kloc,nx,ny);
    int ind3 = get_single_index(iloc,jhi ,kloc,nx,ny);
    int ind4 = get_single_index(iloc,jloc,khi ,nx,ny);
    int ind5 = get_single_index(ihi ,jloc,khi ,nx,ny);
    int ind6 = get_single_index(ihi ,jhi ,khi ,nx,ny);
    int ind7 = get_single_index(iloc,jhi ,khi ,nx,ny);

    tmp_var = (1-w_x)*(1-w_y)*(1-w_z)*varvec[ind0] + w_x*(1-w_y)*(1-w_z)*varvec[ind1] +
          w_x*w_y*(1-w_z)*varvec[ind2] + (1-w_x)*w_y*(1-w_z)*varvec[ind3] +
          (1-w_x)*(1-w_y)*w_z*varvec[ind4] + w_x*(1-w_y)*w_z*varvec[ind5] +
          w_x*w_y*w_z*varvec[ind6] + (1-w_x)*w_y*w_z*varvec[ind7];

    //std::cout << "Variable value is " << tmp_var << "\n";
}

AMREX_GPU_DEVICE
inline
amrex::Real interpolate_from_coarse(const amrex::Array4<const amrex::Real>& crse,
                       int n,
                       amrex::Real x, amrex::Real y, amrex::Real z,
                       const amrex::Real* prob_lo_crse,
                       const amrex::Real* dx_crse)
{
    // Convert to coarse grid index space
    amrex::Real ic_real = (x - prob_lo_crse[0]) / dx_crse[0];
    amrex::Real jc_real = (y - prob_lo_crse[1]) / dx_crse[1];
    amrex::Real kc_real = (z - prob_lo_crse[2]) / dx_crse[2];

    int ic = static_cast<int>(amrex::Math::floor(ic_real));
    int jc = static_cast<int>(amrex::Math::floor(jc_real));
    int kc = static_cast<int>(amrex::Math::floor(kc_real));

    amrex::Real dx = ic_real - ic;
    amrex::Real dy = jc_real - jc;
    amrex::Real dz = kc_real - kc;

    // Trilinear interpolation
    amrex::Real v000 = crse(ic,     jc,     kc,     n);
    amrex::Real v100 = crse(ic + 1, jc,     kc,     n);
    amrex::Real v010 = crse(ic,     jc + 1, kc,     n);
    amrex::Real v110 = crse(ic + 1, jc + 1, kc,     n);
    amrex::Real v001 = crse(ic,     jc,     kc + 1, n);
    amrex::Real v101 = crse(ic + 1, jc,     kc + 1, n);
    amrex::Real v011 = crse(ic,     jc + 1, kc + 1, n);
    amrex::Real v111 = crse(ic + 1, jc + 1, kc + 1, n);

    amrex::Real v = v000 * (1 - dx) * (1 - dy) * (1 - dz) +
            v100 * dx       * (1 - dy) * (1 - dz) +
            v010 * (1 - dx) * dy       * (1 - dz) +
            v110 * dx       * dy       * (1 - dz) +
            v001 * (1 - dx) * (1 - dy) * dz       +
            v101 * dx       * (1 - dy) * dz       +
            v011 * (1 - dx) * dy       * dz       +
            v111 * dx       * dy       * dz;

    return v;
}

AMREX_FORCE_INLINE
AMREX_GPU_HOST_DEVICE
void bilinear_interpolation_2d(const amrex::Real* xvec, const amrex::Real* yvec,
                            const amrex::Real dxvec, const amrex::Real dyvec,
                            const int nx, const int ny,
                            amrex::Real x, amrex::Real y,
                            const amrex::Real* varvec,
                            amrex::Real& tmp_var)
{
        if(x<xvec[0]) {
            x = xvec[0];
        }
        if(x>xvec[nx-1]) {
            x = xvec[nx-1];
        }

        if(y<yvec[0]) {
            y = yvec[0];
        }
        if(y>yvec[ny-1]) {
            y = yvec[ny-1];
        }

        int iloc = clamp_interp_index(static_cast<int>(std::floor((x-xvec[0])/dxvec)), nx);
        int jloc = clamp_interp_index(static_cast<int>(std::floor((y-yvec[0])/dyvec)), ny);
        int ihi = (nx > 1) ? iloc + 1 : iloc;
        int jhi = (ny > 1) ? jloc + 1 : jloc;

        amrex::Real xlo = xvec[iloc];
        amrex::Real ylo = yvec[jloc];

        amrex::Real w_x = get_interp_weight(x, xlo, xvec[ihi]);
        amrex::Real w_y = get_interp_weight(y, ylo, yvec[jhi]);

        int ind0 = get_single_index(iloc,jloc,0,nx,ny);
        int ind1 = get_single_index(ihi ,jloc,0,nx,ny);
        int ind2 = get_single_index(ihi ,jhi ,0,nx,ny);
        int ind3 = get_single_index(iloc,jhi ,0,nx,ny);

        tmp_var = (1-w_x)*(1-w_y)*varvec[ind0] + w_x*(1-w_y)*varvec[ind1] +
                  w_x*w_y*varvec[ind2] + (1-w_x)*w_y*varvec[ind3];

        //std::cout << "Variable value is " << tmp_var << "\n";
}


#endif
