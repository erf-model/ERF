#ifndef ERF_BCBASE_H
#define ERF_BCBASE_H

#include "ERF.H"
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
        if constexpr (IDIR == 0) {
          amrex::ParallelFor(b, ncomp, [=] AMREX_GPU_DEVICE (int i, int j, int k, int n) noexcept
          {
            if (i < lo) {
              data(i,j,k,n) = 2.0*data(lo,j,k,n) - data(lo+1,j,k,n);
            }
          });
        }
        if constexpr (IDIR == 1) {
          amrex::ParallelFor(b, ncomp, [=] AMREX_GPU_DEVICE (int i, int j, int k, int n) noexcept
          {
            if (j < lo) {
              data(i,j,k,n) = 2.0*data(i,lo,k,n) - data(i,lo+1,k,n);
            }
          });
        }
        if constexpr (IDIR == 2) {
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
        if constexpr (IDIR == 0) {
           amrex::ParallelFor(b, ncomp, [=] AMREX_GPU_DEVICE (int i, int j, int k, int n) noexcept
           {
            if (i < lo) {
              data(i,j,k,n) = 1.5*data(lo,j,k,n) - 0.5*data(lo+1,j,k,n);
            } else if (i == lo) {
              data(i,j,k,n) = 2.*data(lo,j,k,n) - data(lo+1,j,k,n);
            }
           });
        }
        if constexpr (IDIR == 1) {
           amrex::ParallelFor(b, ncomp, [=] AMREX_GPU_DEVICE (int i, int j, int k, int n) noexcept
           {
            if (j < lo) {
              data(i,j,k,n) = 1.5*data(i,lo,k,n) - 0.5*data(i,lo+1,k,n);
            } else if (j == lo) {
              data(i,j,k,n) = 2.*data(i,lo,k,n) - data(i,lo+1,k,n);
            }
           });
        }
        if constexpr (IDIR == 2) {
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
  }

  // TODO:
  // need to implement face variable extrapolation
  //
  inline void foextrap_face_lo(const amrex::Geometry& geom, const Box& b, int lo, Array4<Real> data, int idir, int ncomp = 0) {
    if (idir == 0) {
       amrex::ParallelFor(b, ncomp, [=] AMREX_GPU_DEVICE (int i, int j, int k, int n) noexcept
       {
        if (i < lo) {
          data(i,j,k,n) = 2.*data(lo,j,k,n) - data(lo+1,j,k,n);
        }
       });
    } else if (idir == 1) {
       amrex::ParallelFor(b, ncomp, [=] AMREX_GPU_DEVICE (int i, int j, int k, int n) noexcept
       {
        if (j < lo) {
          data(i,j,k,n) = 2.*data(lo,j,k,n) - data(i,lo+1,k,n);
        }
       });
    } else if (idir == 2) {
       amrex::ParallelFor(b, ncomp, [=] AMREX_GPU_DEVICE (int i, int j, int k, int n) noexcept
       {
        if (k < lo) {
          data(i,j,k,n) = 2.*data(lo,j,k,n) - data(i,j,lo+1,n);
        }
       });
  template<>
  inline void foextrap_face_based<BCBound::lower>(const amrex::Geometry& geom, const Box& b, int lo, Array4<Real> data, int idir, int ncomp) {
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

  template<>
  inline void foextrap_face_based<BCBound::lower>(const amrex::Geometry& geom, const Box& b, int lo, Array4<Real> data, int idir, int ncomp) {
    if(b.smallEnd(idir) < lo) {
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
        if constexpr (IDIR == 0) {
          amrex::ParallelFor(b, ncomp, [=] AMREX_GPU_DEVICE (int i, int j, int k, int n) noexcept
          {
            if (i > hi) {
              data(i,j,k,n) = 2.*data(hi,j,k,n) - data(hi-1,j,k,n);
            }
          });
        }
        if constexpr (IDIR == 1) {
          amrex::ParallelFor(b, ncomp, [=] AMREX_GPU_DEVICE (int i, int j, int k, int n) noexcept
          {
            if (j > hi) {
              data(i,j,k,n) = 2.*data(i,hi,k,n) - data(i,hi-1,k,n);
            }
          });
        }
        if constexpr (IDIR == 2) {
          amrex::ParallelFor(b, ncomp, [=] AMREX_GPU_DEVICE (int i, int j, int k, int n) noexcept
          {
            if (k > hi) {
              data(i,j,k,n) = 2.*data(i,j,hi,n) - data(i,j,hi-1,n);
            }
          });
        }
      }
    }
  }


  // TODO:
  // need to implement face variable extrapolation algorithm
  //
  inline void foextrap_face_hi(const amrex::Geometry& geom, const Box& b, int hi, Array4<Real> data, int idir, int ncomp = 0) {
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
           data(i,j,k,n) = 2.*data(hi,j,k,n) - data(i,hi-1,k,n);
         }
       });
    } else if (idir == 2) {
       amrex::ParallelFor(b, ncomp, [=] AMREX_GPU_DEVICE (int i, int j, int k, int n) noexcept
       {
         if (k > hi) {
           data(i,j,k,n) = 2.*data(hi,j,k,n) - data(i,j,hi-1,n);
         }
       });
   }
  template<>
  inline void foextrap_face_based<BCBound::upper>(const amrex::Geometry& geom, const Box& b, int hi, Array4<Real> data, int idir, int ncomp) {
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

  template<>
  inline void foextrap_face_based<BCBound::upper>(const amrex::Geometry& geom, const Box& b, int hi, Array4<Real> data, int idir, int ncomp) {
    if(b.bigEnd(idir) > hi) {
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

    template <bool val>
    static inline 
    std::enable_if_t<val>
    apply_face_based(const amrex::Geometry& geom, const Box& b, Array4<Real> data, int ncomp=1) {
      int hi = geom.Domain().bigEnd(IDIR);
 
      if(b.bigEnd(IDIR) > hi) {
        if constexpr (IDIR == 0) {
           amrex::ParallelFor(b, ncomp, [=] AMREX_GPU_DEVICE (int i, int j, int k, int n) noexcept
           {
             if (i > hi) {
               data(i,j,k,n) = 1.5*data(hi,j,k,n) - 0.5*data(hi-1,j,k,n);
             } else if (i == hi) {
               data(i,j,k,n) = 2.*data(hi,j,k,n) - data(hi-1,j,k,n);
             }
           });
        }
        if constexpr (IDIR == 1) {
           amrex::ParallelFor(b, ncomp, [=] AMREX_GPU_DEVICE (int i, int j, int k, int n) noexcept
           {
             if (j > hi) {
               data(i,j,k,n) = 1.5*data(i,hi,k,n) - 0.5*data(i,hi-1,k,n);
             } else if (j == hi) {
               data(i,j,k,n) = 2.*data(i,hi,k,n) - data(i,hi-1,k,n);
             }
           });
        }
        if constexpr (IDIR == 2) {
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

}

namespace phys_bcs {

class BCBase {
  public:
  virtual void applyBC (const amrex::Geometry geom, amrex::Vector<MultiFab*>& vars) = 0;
  virtual bool isInterior() { return false; }
};

}

#endif
