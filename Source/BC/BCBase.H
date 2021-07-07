#ifndef ERF_BCBASE_H
#define ERF_BCBASE_H

#include "ERF.H"
using namespace amrex;

namespace math_bcs {
  inline void foextrap_lo(const amrex::Geometry& geom, const Box& b, int lo, Array4<Real> data) {
    amrex::ParallelFor(b, [=] AMREX_GPU_DEVICE (int i, int j, int k) noexcept
    {
      if (i < lo) {
        data(i,j,k) = 2.0*data(lo,j,k) - data(lo+1,j,k);
      }
    });
  }

  inline void foextrap_hi(const amrex::Geometry& geom, const Box& b, int hi, Array4<Real> data) {
    amrex::ParallelFor(b, [=] AMREX_GPU_DEVICE (int i, int j, int k) noexcept
    {
      if (i > hi) {
        data(i,j,k) = 2.*data(hi,j,k) - data(hi-1,j,k);
      }
    });
  }

/*
  inline void DirichletBC(const amrex::Geometry& geom, const Box& b, int ncomp, Array4<Real> arr, ) {
    int nx = geom.Domain().bigEnd(0);
    amrex::ParallelFor(b, [=] AMREX_GPU_DEVICE (int i, int j, int k) noexcept
    {
        arr(-1,j,k,ncomp) = 0.0*2 - arr(0,j,k,ncomp);
        arr(nx+1,j,k,ncomp) = 1.0*2 - arr(nx,j,k,ncomp);
    });
  }

  inline void NuemannBC(const amrex::Geometry& geom, const Box& b, int ncomp, Array4<Real> arr) {
    Real dx = geom.CellSize();
    int nx = geom.Domain().bigEnd(0);
    amrex::ParallelFor(b, [=] AMREX_GPU_DEVICE (int i, int j, int k) noexcept
    {
        arr(-1,j,k,ncomp) = 0.0*dx + arr(0,j,k,ncomp);
        arr(nx+1,j,k,ncomp) = 1.0*dx + arr(nx,j,k,ncomp);
    });
  }
*/
}

namespace phys_bcs {

class BCBase {
  public:
  virtual void applyBC (const amrex::Geometry geom, amrex::Vector<MultiFab*>& vars, int dim) = 0;
  virtual bool isInterior() { return false; }
};

}

#endif
