#ifndef ERF_BCBASE_H
#define ERF_BCBASE_H

#include "ERF.H"
using namespace amrex;

namespace math_bcs {

  enum BCBound { lower, upper };

  template<BCBound b>
  inline void foextrap(const amrex::Geometry& geom, const Box& b, int lo, Array4<Real> data, int idir, int ncomp = 1);

  template<BCBound b>
  inline void foextrap_edgebased(const amrex::Geometry& geom, const Box& b, int lo, Array4<Real> data, int idir);

  template<>
  inline void foextrap<BCBound::lower>(const amrex::Geometry& geom, const Box& b, int lo, Array4<Real> data, int idir, int ncomp = 1) {
    if(b.smallEnd(idir) < lo) {
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
  }

  template<>
  inline void foextrap<BCBound::upper>(const amrex::Geometry& geom, const Box& b, int hi, Array4<Real> data, int idir, int ncomp = 1) {
    if(b.bigEnd(idir) > hi) {
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
  }

  template<>
  inline void foextrap_edgebased<BCBound::lower>(const amrex::Geometry& geom, const Box& b, int lo, Array4<Real> data, int idir) {
    if(b.smallEnd(idir) < lo) {
      if (idir == 0) {
        amrex::ParallelFor(b, [=] AMREX_GPU_DEVICE (int i, int j, int k) noexcept
        {
          if (i <= lo) {
            data(i,j,k) = data(lo+1,j,k);
          }
        });
      } else if (idir == 1) {
        amrex::ParallelFor(b, [=] AMREX_GPU_DEVICE (int i, int j, int k) noexcept
        {
          if (j <= lo) {
            data(i,j,k) = data(i,lo+1,k);
          }
        });
      } else if (idir == 2) {
        amrex::ParallelFor(b, [=] AMREX_GPU_DEVICE (int i, int j, int k) noexcept
        {
          if (k <= lo) {
            data(i,j,k) = data(i,j,lo+1);
          }
        });
      }
    }
  }

  template<>
  inline void foextrap_edgebased<BCBound::upper>(const amrex::Geometry& geom, const Box& b, int hi, Array4<Real> data, int idir) {
    if(b.bigEnd(idir) > hi) {
      if (idir == 0) {
        amrex::ParallelFor(b, [=] AMREX_GPU_DEVICE (int i, int j, int k) noexcept
        {
          if (i >= hi) {
            data(i,j,k) = data(hi-1,j,k);
          }
        });
      } else if (idir == 1) {
        amrex::ParallelFor(b, [=] AMREX_GPU_DEVICE (int i, int j, int k) noexcept
        {
          if (j >= hi) {
            data(i,j,k) = data(i,hi-1,k);
          }
        });
      } else if (idir == 2) {
        amrex::ParallelFor(b, [=] AMREX_GPU_DEVICE (int i, int j, int k) noexcept
        {
          if (k >= hi) {
            data(i,j,k) = data(i,j,hi-1);
          }
        });
      }
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
