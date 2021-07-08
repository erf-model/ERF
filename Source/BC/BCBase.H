#ifndef ERF_BCBASE_H
#define ERF_BCBASE_H

#include "ERF.H"
using namespace amrex;

namespace math_bcs {
  inline void foextrap_cell_lo(const amrex::Geometry& geom, const Box& b, int lo, Array4<Real> data, int idir, int ncomp = 0) {
    if (idir == 0) {
       amrex::ParallelFor(b, ncomp, [=] AMREX_GPU_DEVICE (int i, int j, int k, int n) noexcept
       {
         if (i < lo) {
           data(i,j,k,n) = 2.0*data(lo,j,k,n) - data(lo+1,j,k,n);
         }
       });
    } else if (idir == 1) {
       amrex::ParallelFor(b, ncomp, [=] AMREX_GPU_DEVICE (int i, int j, int k, int n) noexcept
       {
         if (j < lo) {
           data(i,j,k,n) = 2.0*data(i,lo,k,n) - data(i,lo+1,k,n);
         }
       });

    } else if (idir == 2) {
       amrex::ParallelFor(b, ncomp, [=] AMREX_GPU_DEVICE (int i, int j, int k, int n) noexcept
       {
         if (k < lo) {
           data(i,j,k,n) = 2.0*data(i,j,lo,n) - data(i,j,lo+1,n);
         }
       });
    }
  }

  inline void foextrap_cell_hi(const amrex::Geometry& geom, const Box& b, int hi, Array4<Real> data, int idir, int ncomp = 0) {
    if (idir == 0) {
       amrex::ParallelFor(b, ncomp, [=] AMREX_GPU_DEVICE (int i, int j, int k, int n) noexcept
       {
         if (i > hi) {
           data(i,j,k,n) = 2.*data(hi,j,k,n) - data(hi-1,j,k,n);
         }
       });
    } else if (idir == 1) {
       amrex::ParallelFor(b, ncomp, [=] AMREX_GPU_DEVICE (int i, int j, int k, int n) noexcept
       {
          if (j > hi) {
            data(i,j,k,n) = 2.*data(i,hi,k,n) - data(i,hi-1,k,n);
          }
       });
    } else if (idir == 2) {
       amrex::ParallelFor(b, ncomp, [=] AMREX_GPU_DEVICE (int i, int j, int k, int n) noexcept
       {
         if (k > hi) {
           data(i,j,k,n) = 2.*data(i,j,hi,n) - data(i,j,hi-1,n);
         }
       });
    }
  }

  inline void foextrap_face_lo(const amrex::Geometry& geom, const Box& b, int lo, Array4<Real> data, int idir, int ncomp = 0) {
    if (idir == 0) {
       amrex::ParallelFor(b, ncomp, [=] AMREX_GPU_DEVICE (int i, int j, int k, int n) noexcept
       {
        if (i < lo) {
          data(i,j,k,n) = 1.5*data(lo,j,k,n) - 0.5*data(lo+1,j,k,n);
        } else if (i == lo) {
          data(i,j,k,n) = 2.*data(lo,j,k,n) - data(lo+1,j,k,n);
        }
       });
    } else if (idir == 1) {
       amrex::ParallelFor(b, ncomp, [=] AMREX_GPU_DEVICE (int i, int j, int k, int n) noexcept
       {
        if (j < lo) {
          data(i,j,k,n) = 1.5*data(i,lo,k,n) - 0.5*data(i,lo+1,k,n);
        } else if (j == lo) {
          data(i,j,k,n) = 2.*data(i,lo,k,n) - data(i,lo+1,k,n);
        }
       });
    } else if (idir == 2) {
       amrex::ParallelFor(b, ncomp, [=] AMREX_GPU_DEVICE (int i, int j, int k, int n) noexcept
       {
        if (k < lo) {
          data(i,j,k,n) = 1.5*data(i,j,lo,n) - 0.5*data(i,j,lo+1,n);
        } else if (k == lo) {
          data(i,j,k,n) = 2.*data(i,j,lo,n) - data(i,j,lo+1,n);
        }
       });
    }
  }

  inline void foextrap_face_hi(const amrex::Geometry& geom, const Box& b, int hi, Array4<Real> data, int idir, int ncomp = 0) {
    if (idir == 0) {
       amrex::ParallelFor(b, ncomp, [=] AMREX_GPU_DEVICE (int i, int j, int k, int n) noexcept
       {
         if (i > hi) {
           data(i,j,k,n) = 1.5*data(hi,j,k,n) - 0.5*data(hi-1,j,k,n);
         } else if (i == hi) {
           data(i,j,k,n) = 2.*data(hi,j,k,n) - data(hi-1,j,k,n);
         }
       });
    } else if (idir == 1) {
       amrex::ParallelFor(b, ncomp, [=] AMREX_GPU_DEVICE (int i, int j, int k, int n) noexcept
       {
         if (j > hi) {
           data(i,j,k,n) = 1.5*data(i,hi,k,n) - 0.5*data(i,hi-1,k,n);
         } else if (j == hi) {
           data(i,j,k,n) = 2.*data(i,hi,k,n) - data(i,hi-1,k,n);
         }
       });
    } else if (idir == 2) {
       amrex::ParallelFor(b, ncomp, [=] AMREX_GPU_DEVICE (int i, int j, int k, int n) noexcept
       {
         if (k > hi) {
           data(i,j,k,n) = 1.5*data(i,j,hi,n) - 0.5*data(i,j,hi-1,n);
         } else if (k == hi) {
           data(i,j,k,n) = 2.*data(i,j,hi,n) - data(i,j,hi-1,n);
         }
       });
   }
  }


}

namespace phys_bcs {

class BCBase {
  public:
  virtual void applyBC (const amrex::Geometry geom, amrex::Vector<MultiFab*>& vars, int dim) = 0;
  virtual bool isInterior() { return false; }
};

}

#endif
