#ifndef ERF_BCBASE_H
#define ERF_BCBASE_H

#include "ERF.H"
#include "IndexDefines.H"
using namespace amrex;

namespace math_bcs {

  enum BCBound { lower, upper };

  template<int IDIR, BCBound Bound> struct foextrap;

  template<int IDIR>
  struct foextrap<IDIR, BCBound::lower> {
    using self = foextrap<IDIR, BCBound::lower>;

    static inline void apply_cell_based(const amrex::Geometry& geom, const Box& b, Array4<Real> data, int ncomp=1) {
      int lo = geom.Domain().smallEnd(IDIR);
 
      if(b.smallEnd(IDIR) < lo) {
        if (IDIR == 0) {
          amrex::ParallelFor(b, ncomp, [=] AMREX_GPU_DEVICE (int i, int j, int k, int n) noexcept
          {
            if (i < lo) {
              data(i,j,k,n) = 2.0*data(lo,j,k,n) - data(lo+1,j,k,n);
            }
          });
        } else if (IDIR == 1) {
          amrex::ParallelFor(b, ncomp, [=] AMREX_GPU_DEVICE (int i, int j, int k, int n) noexcept
          {
            if (j < lo) {
              data(i,j,k,n) = 2.0*data(i,lo,k,n) - data(i,lo+1,k,n);
            }
          });
        } else if (IDIR == 2) {
          amrex::ParallelFor(b, ncomp, [=] AMREX_GPU_DEVICE (int i, int j, int k, int n) noexcept
          {
            if (k < lo) {
              data(i,j,k,n) = 2.0*data(i,j,lo,n) - data(i,j,lo+1,n);
            }
          });
        }
      }
    }

    template <bool val>
    static inline 
    std::enable_if_t<val>
    apply_face_based(const amrex::Geometry& geom, const Box& b, Array4<Real> data, int ncomp=1) {
      int lo = geom.Domain().smallEnd(IDIR);
 
      if(b.smallEnd(IDIR) < lo) {
        if (IDIR == 0) {
           amrex::ParallelFor(b, ncomp, [=] AMREX_GPU_DEVICE (int i, int j, int k, int n) noexcept
           {
            if (i < lo) {
              data(i,j,k,n) = 1.5*data(lo,j,k,n) - 0.5*data(lo+1,j,k,n);
            } else if (i == lo) {
              data(i,j,k,n) = 2.*data(lo,j,k,n) - data(lo+1,j,k,n);
            }
           });
        } else if (IDIR == 1) {
           amrex::ParallelFor(b, ncomp, [=] AMREX_GPU_DEVICE (int i, int j, int k, int n) noexcept
           {
            if (j < lo) {
              data(i,j,k,n) = 1.5*data(i,lo,k,n) - 0.5*data(i,lo+1,k,n);
            } else if (j == lo) {
              data(i,j,k,n) = 2.*data(i,lo,k,n) - data(i,lo+1,k,n);
            }
           });
        } else if (IDIR == 2) {
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
    }

    template <bool val>
    static inline 
    std::enable_if_t<!val>
    apply_face_based(const amrex::Geometry& geom, const Box& b, Array4<Real> data, int ncomp=1) {
      self::apply_cell_based(geom, b, data, ncomp);
    }

  };

  template<int IDIR>
  struct foextrap<IDIR, BCBound::upper> {
    using self = foextrap<IDIR, BCBound::upper>;

    static inline void apply_cell_based(const amrex::Geometry& geom, const Box& b, Array4<Real> data, int ncomp=1) {
      int hi = geom.Domain().bigEnd(IDIR);
 
      if(b.bigEnd(IDIR) > hi) {
        if (IDIR == 0) {
          amrex::ParallelFor(b, ncomp, [=] AMREX_GPU_DEVICE (int i, int j, int k, int n) noexcept
          {
            if (i > hi) {
              data(i,j,k,n) = 2.*data(hi,j,k,n) - data(hi-1,j,k,n);
            }
          });
        } else if (IDIR == 1) {
          amrex::ParallelFor(b, ncomp, [=] AMREX_GPU_DEVICE (int i, int j, int k, int n) noexcept
          {
            if (j > hi) {
              data(i,j,k,n) = 2.*data(i,hi,k,n) - data(i,hi-1,k,n);
            }
          });
        } else if (IDIR == 2) {
          amrex::ParallelFor(b, ncomp, [=] AMREX_GPU_DEVICE (int i, int j, int k, int n) noexcept
          {
            if (k > hi) {
              data(i,j,k,n) = 2.*data(i,j,hi,n) - data(i,j,hi-1,n);
            }
          });
        }
      }
    }

    template <bool val>
    static inline 
    std::enable_if_t<val>
    apply_face_based(const amrex::Geometry& geom, const Box& b, Array4<Real> data, int ncomp=1) {
      int hi = geom.Domain().bigEnd(IDIR);
 
      if(b.bigEnd(IDIR) > hi) {
        if (IDIR == 0) {
           amrex::ParallelFor(b, ncomp, [=] AMREX_GPU_DEVICE (int i, int j, int k, int n) noexcept
           {
             if (i > hi) {
               data(i,j,k,n) = 1.5*data(hi,j,k,n) - 0.5*data(hi-1,j,k,n);
             } else if (i == hi) {
               data(i,j,k,n) = 2.*data(hi,j,k,n) - data(hi-1,j,k,n);
             }
           });
        } else if (IDIR == 1) {
           amrex::ParallelFor(b, ncomp, [=] AMREX_GPU_DEVICE (int i, int j, int k, int n) noexcept
           {
             if (j > hi) {
               data(i,j,k,n) = 1.5*data(i,hi,k,n) - 0.5*data(i,hi-1,k,n);
             } else if (j == hi) {
               data(i,j,k,n) = 2.*data(i,hi,k,n) - data(i,hi-1,k,n);
             }
           });
        } else if (IDIR == 2) {
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

    template <bool val>
    static inline 
    std::enable_if_t<!val>
    apply_face_based(const amrex::Geometry& geom, const Box& b, Array4<Real> data, int ncomp=1) {
      self::apply_cell_based(geom, b, data, ncomp);
    }

  };

  template<>
  inline void extval<BCBound::lower>(const amrex::Geometry& geom, const Box& b, int lo, Array4<Real> data, int idir, int ncomp) {
    amrex::IntVect type = b.type();
    amrex::Vector bvalue(ncomp);
    if (type == amrex::IntVect::Zero) {
       amrex::Real rho0   = 0.0;
       amrex::Real theta0 = 0.0;
       amrex::Real press0 = 0.0;

       pp.query("rho0"  , rho0);
       pp.query("theta0", theta0);
       pp.query("press0", press0);

       bvalue[Rho_comp]       = rho0;
       bvalue[RhoTheta_comp]  = theta0;
       bvalue[RhoScalar_comp] = press0;
    } else if (type == amrex::IntVect(AMREX_D_DECL(1, 0, 0))) { // for X_VEL 
       amrex::Real velocity_x = 0.0;
       pp.query("velocity_x", velocity_x);
       bvalue[0] = velocity_x;
    } else if (type == amrex::IntVect(AMREX_D_DECL(0, 1, 0))) { // for Y_VEL 
       amrex::Real velocity_y = 0.0;
       pp.query("velocity_y", velocity_y);
       bvalue[0] = velocity_y;
    } else if (type == amrex::IntVect(AMREX_D_DECL(0, 1, 0))) { // for Z_VEL 
       amrex::Real velocity_y = 0.0;
       pp.query("velocity_y", velocity_y);
       bvalue[0] = velocity_y;
    }    

    if(b.smallEnd(idir) < lo) {
      if (idir == 0) {
         amrex::ParallelFor(b, ncomp, [=] AMREX_GPU_DEVICE (int i, int j, int k, int n) noexcept
         {
          if (i < lo) {
            data(i,j,k,n) = bvalue[n];
          }
         });
      } else if (idir == 1) {
         amrex::ParallelFor(b, ncomp, [=] AMREX_GPU_DEVICE (int i, int j, int k, int n) noexcept
         {
          if (j < lo) {
            data(i,j,k,n) = bvalue[n];
          }
         });
      } else if (idir == 2) {
         amrex::ParallelFor(b, ncomp, [=] AMREX_GPU_DEVICE (int i, int j, int k, int n) noexcept
         {
          if (k < lo) {
            data(i,j,k,n) = bvalue[n];
          }
         });
      }
    }
  }

  template<>
  inline void extval_face<BCBound::lower>(const amrex::Geometry& geom, const Box& b, int lo, Array4<Real> data, int idir, int ncomp) {
    amrex::IntVect type = b.type();
    amrex::Vector bvalue(ncomp);
    if (type == amrex::IntVect::Zero) {
       amrex::Real rho0   = 0.0;
       amrex::Real theta0 = 0.0;
       amrex::Real press0 = 0.0;

       pp.query("rho0"  , rho0);
       pp.query("theta0", theta0);
       pp.query("press0", press0);

       bvalue[Rho_comp]       = rho0;
       bvalue[RhoTheta_comp]  = theta0;
       bvalue[RhoScalar_comp] = press0;
    } else if (type == amrex::IntVect(AMREX_D_DECL(1, 0, 0))) { // for X_VEL 
       amrex::Real velocity_x = 0.0;
       pp.query("velocity_x", velocity_x);
       bvalue[0] = velocity_x;
    } else if (type == amrex::IntVect(AMREX_D_DECL(0, 1, 0))) { // for Y_VEL 
       amrex::Real velocity_y = 0.0;
       pp.query("velocity_y", velocity_y);
       bvalue[0] = velocity_y;
    } else if (type == amrex::IntVect(AMREX_D_DECL(0, 1, 0))) { // for Z_VEL 
       amrex::Real velocity_y = 0.0;
       pp.query("velocity_y", velocity_y);
       bvalue[0] = velocity_y;
    }

    if(b.smallEnd(idir) < lo) {
      if (idir == 0) {
         amrex::ParallelFor(b, ncomp, [=] AMREX_GPU_DEVICE (int i, int j, int k, int n) noexcept
         {
          if (i <= lo) {
            data(i,j,k,n) = bvalue[n];
          }
         });
      } else if (idir == 1) {
         amrex::ParallelFor(b, ncomp, [=] AMREX_GPU_DEVICE (int i, int j, int k, int n) noexcept
         {
          if (j <= lo) {
            data(i,j,k,n) = bvalue[n];
          }
         });
      } else if (idir == 2) {
         amrex::ParallelFor(b, ncomp, [=] AMREX_GPU_DEVICE (int i, int j, int k, int n) noexcept
         {
          if (k <= lo) {
            data(i,j,k,n) = bvalue[n];
          }
         });
      }
    }
  }

  template<>
  inline void extval<BCBound::upper>(const amrex::Geometry& geom, const Box& b, int hi, Array4<Real> data, int idir, int ncomp) {
    amrex::IntVect type = b.type();
    amrex::Vector bvalue(ncomp);
    if (type == amrex::IntVect::Zero) {
       amrex::Real rho0   = 0.0;
       amrex::Real theta0 = 0.0;
       amrex::Real press0 = 0.0;

       pp.query("rho0"  , rho0);
       pp.query("theta0", theta0);
       pp.query("press0", press0);

       bvalue[Rho_comp]       = rho0;
       bvalue[RhoTheta_comp]  = theta0;
       bvalue[RhoScalar_comp] = press0;
    } else if (type == amrex::IntVect(AMREX_D_DECL(1, 0, 0))) { // for X_VEL 
       amrex::Real velocity_x = 0.0;
       pp.query("velocity_x", velocity_x);
       bvalue[0] = velocity_x;
    } else if (type == amrex::IntVect(AMREX_D_DECL(0, 1, 0))) { // for Y_VEL 
       amrex::Real velocity_y = 0.0;
       pp.query("velocity_y", velocity_y);
       bvalue[0] = velocity_y;
    } else if (type == amrex::IntVect(AMREX_D_DECL(0, 1, 0))) { // for Z_VEL 
       amrex::Real velocity_y = 0.0;
       pp.query("velocity_y", velocity_y);
       bvalue[0] = velocity_y;
    }

    if(b.smallEnd(idir) < lo) {
      if (idir == 0) {
         amrex::ParallelFor(b, ncomp, [=] AMREX_GPU_DEVICE (int i, int j, int k, int n) noexcept
         {
          if (i > hi) {
            data(i,j,k,n) = bvalue[n];
          }
         });
      } else if (idir == 1) {
         amrex::ParallelFor(b, ncomp, [=] AMREX_GPU_DEVICE (int i, int j, int k, int n) noexcept
         {
          if (j > hi) {
            data(i,j,k,n) = bvalue[n];
          }
         });
      } else if (idir == 2) {
         amrex::ParallelFor(b, ncomp, [=] AMREX_GPU_DEVICE (int i, int j, int k, int n) noexcept
         {
          if (k > hi) {
            data(i,j,k,n) = bvalue[n];
          }
         });
      }
    }
  }

  template<>
  inline void extval<BCBound::upper>(const amrex::Geometry& geom, const Box& b, int hi, Array4<Real> data, int idir, int ncomp) {
    amrex::IntVect type = b.type();
    amrex::Vector bvalue(ncomp);
    if (type == amrex::IntVect::Zero) {
       amrex::Real rho0   = 0.0;
       amrex::Real theta0 = 0.0;
       amrex::Real press0 = 0.0;

       pp.query("rho0"  , rho0);
       pp.query("theta0", theta0);
       pp.query("press0", press0);

       bvalue[Rho_comp]       = rho0;
       bvalue[RhoTheta_comp]  = theta0;
       bvalue[RhoScalar_comp] = press0;
    } else if (type == amrex::IntVect(AMREX_D_DECL(1, 0, 0))) { // for X_VEL 
       amrex::Real velocity_x = 0.0;
       pp.query("velocity_x", velocity_x);
       bvalue[0] = velocity_x;
    } else if (type == amrex::IntVect(AMREX_D_DECL(0, 1, 0))) { // for Y_VEL 
       amrex::Real velocity_y = 0.0;
       pp.query("velocity_y", velocity_y);
       bvalue[0] = velocity_y;
    } else if (type == amrex::IntVect(AMREX_D_DECL(0, 1, 0))) { // for Z_VEL 
       amrex::Real velocity_y = 0.0;
       pp.query("velocity_y", velocity_y);
       bvalue[0] = velocity_y;
    }

    if(b.bigEnd(idir) > hi) {
      if (idir == 0) {
         amrex::ParallelFor(b, ncomp, [=] AMREX_GPU_DEVICE (int i, int j, int k, int n) noexcept
         {
          if (i >= hi) {
            data(i,j,k,n) = bvalue[n];
          }
         });
      } else if (idir == 1) {
         amrex::ParallelFor(b, ncomp, [=] AMREX_GPU_DEVICE (int i, int j, int k, int n) noexcept
         {
          if (j >= hi) {
            data(i,j,k,n) = bvalue[n];
          }
         });
      } else if (idir == 2) {
         amrex::ParallelFor(b, ncomp, [=] AMREX_GPU_DEVICE (int i, int j, int k, int n) noexcept
         {
          if (k >= hi) {
            data(i,j,k,n) = bvalue[n];
          }
         });
      }
    }
  }


}

namespace phys_bcs {

class BCBase {
  public:
  virtual void applyBC (const amrex::Geometry geom, amrex::Vector<MultiFab*>& vars) = 0;
  virtual bool isInterior() { return false; }
};

}

#endif
