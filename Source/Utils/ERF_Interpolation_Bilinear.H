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
        iloc = std::floor((x-xvec[0])/dxvec);
        jloc = std::floor((y-yvec[0])/dyvec);

        if(iloc > nx-1 or iloc < 0 or
           jloc > ny-1 or iloc < 0 or
           kloc > nz-1 or kloc < 0){
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

#endif
