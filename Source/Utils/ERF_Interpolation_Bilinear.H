#ifndef ERF_INTERPOLATE_BILINEAR_H_
#define ERF_INTERPOLATE_BILINEAR_H_

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
void bilinear_interpolation(const amrex::Real* xvec, const amrex::Real* yvec, const amrex::Real* zvec,
                const amrex::Real dxvec, const amrex::Real dyvec,
                const int nx, const int ny, const int nz,
                const amrex::Real x, const amrex::Real y, const amrex::Real z,
                const amrex::Real* varvec,
                amrex::Real& tmp_var)
{

    if(z < 0.0){
        tmp_var = 0.0;
        return;
    }

    int iloc=-1, jloc=-1, kloc=-1;
    for(int k=0;k<nz;k++){
        if(zvec[k] > z){
            kloc = k-1;
            break;
        }
        else if (zvec[k] == z) {
            kloc = k;
        }
    }

    iloc = static_cast<int>(std::floor((x - xvec[0]) / dxvec));
    jloc = static_cast<int>(std::floor((y - yvec[0]) / dyvec));

    if(iloc > nx-1 or iloc < 0 or
       jloc > ny-1 or iloc < 0 or
       kloc > nz-1 or kloc < 0){
        printf("There is a problem inside here %d %d %d %d %d %d\n", iloc, jloc, kloc, nx-1, ny-1, nz-1);
        for(int k=0;k<nz;k++){
            printf("Value of z is %d %0.15g %0.15g\n", k, zvec[k], z);
       }
    }

    amrex::Real xlo = xvec[0] + iloc*dxvec;
    amrex::Real ylo = yvec[0] + jloc*dyvec;
    amrex::Real zlo = zvec[kloc];

    amrex::Real w_x = (x - xlo)/dxvec;
    amrex::Real w_y = (y - ylo)/dyvec;
    amrex::Real w_z = (z - zlo)/(zvec[kloc+1] - zvec[kloc]);

    int ind0 = get_single_index(iloc,jloc,kloc,nx,ny);
    int ind1 = get_single_index(iloc+1,jloc,kloc,nx,ny);
    int ind2 = get_single_index(iloc+1,jloc+1,kloc,nx,ny);
    int ind3 = get_single_index(iloc,jloc+1,kloc,nx,ny);
    int ind4 = get_single_index(iloc,jloc,kloc+1,nx,ny);
    int ind5 = get_single_index(iloc+1,jloc,kloc+1,nx,ny);
    int ind6 = get_single_index(iloc+1,jloc+1,kloc+1,nx,ny);
    int ind7 = get_single_index(iloc,jloc+1,kloc+1,nx,ny);

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


#endif
