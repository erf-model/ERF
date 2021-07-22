#ifndef ERF_BCBASE_H
#define ERF_BCBASE_H

#include "ERF.H"
#include "IndexDefines.H"
#include "Constants.H"
#include "EOS.H"

using namespace amrex;

namespace math_bcs {

  enum BCBound { lower, upper };
  
  template<BCBound Bound>
  inline void brian_reflect_odd(const amrex::Geometry& geom, const Box& b, int extent, Array4<Real> data, int idir, int ncomp = 1);

  template<BCBound Bound>
  inline void brian_reflect_odd_face_based(const amrex::Geometry& geom, const Box& b, int extent, Array4<Real> data, int idir, int ncomp = 1);

  template<BCBound Bound>
  inline void wall_scalars(const amrex::Geometry& geom, const Box& b, int extent, Array4<Real> data, int idir, int ncomp = 1);


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
              amrex::Real grad = data(lo+1,j,k,n)-data(lo,j,k,n);
              data(i,j,k,n) = data(lo,j,k,n) - (lo-i)*grad;
            }
          });
        } else if (IDIR == 1) {
          amrex::ParallelFor(b, ncomp, [=] AMREX_GPU_DEVICE (int i, int j, int k, int n) noexcept
          {
            if (j < lo) {
              amrex::Real grad = data(i,lo+1,k,n) - data(i,lo,k,n);
              data(i,j,k,n) = data(i,lo,k,n) - (lo-j)*grad;
            }
          });
        } else if (IDIR == 2) {
          amrex::ParallelFor(b, ncomp, [=] AMREX_GPU_DEVICE (int i, int j, int k, int n) noexcept
          {
            if (k < lo) {
              amrex::Real grad = data(i,j,lo+1,n) - data(i,j,lo,n);
              data(i,j,k,n) = data(i,j,lo,n) - (lo-k)*grad;
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
              amrex::Real grad = data(lo+2,j,k,n) - data(lo+1,j,k,n);
              data(i,j,k,n) = data(lo+1,j,k,n) - (lo+1-i)*grad;
            } else if (i == lo) {
              amrex::Real grad = data(lo+2,j,k,n) - data(lo+1,j,k,n);
              data(i,j,k,n) = data(lo+1,j,k,n) - (lo+1-i)*grad;
            }
           });
        } else if (IDIR == 1) {
           amrex::ParallelFor(b, ncomp, [=] AMREX_GPU_DEVICE (int i, int j, int k, int n) noexcept
           {
            if (j < lo) {
              amrex::Real grad = data(i,lo+2,k,n) - data(i,lo+1,k,n);
              data(i,j,k,n) = data(i,lo+1,k,n) - (lo+1-j)*grad;
            } else if (j == lo) {
              amrex::Real grad = data(i,j,lo+2,n) - data(i,j,lo+1,n);
              data(i,j,k,n) = data(i,lo+1,k,n) - (lo+1-j)*grad;
            }
           });
        } else if (IDIR == 2) {
           amrex::ParallelFor(b, ncomp, [=] AMREX_GPU_DEVICE (int i, int j, int k, int n) noexcept
           {
            if (k < lo) {
              amrex::Real grad = data(i,j,lo+2,n) - data(i,j,lo+1,n);
              data(i,j,k,n) = data(i,j,lo+1,n) - (lo+1-k)*grad;
            } else if (k == lo) {
              amrex::Real grad = data(i,j,lo+2,n) - data(i,j,lo+1,n);
              data(i,j,k,n) = data(i,j,lo+1,n) - (lo+1-k)*grad;
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
              amrex::Real grad = data(hi,j,k,n) - data(hi-1,j,k,n);
              data(i,j,k,n) = data(hi,j,k,n) - (hi-i)*grad;
            }
          });
        } else if (IDIR == 1) {
          amrex::ParallelFor(b, ncomp, [=] AMREX_GPU_DEVICE (int i, int j, int k, int n) noexcept
          {
            if (j > hi) {
              amrex::Real grad = data(i,hi,k,n)-data(i,hi-1,k,n);
              data(i,j,k,n) = data(i,hi,k,n) - (hi-j)*grad;
            }
          });
        } else if (IDIR == 2) {
          amrex::ParallelFor(b, ncomp, [=] AMREX_GPU_DEVICE (int i, int j, int k, int n) noexcept
          {
            if (k > hi) {
              amrex::Real grad = data(i,j,hi,n) - data(i,j,hi-1,n);
              data(i,j,k,n) = data(i,j,hi,n) - (hi-k)*grad;
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
               amrex::Real grad = data(hi-1,j,k,n) - data(hi-2,j,k,n);
               data(i,j,k,n) = data(hi-1,j,k,n) - (hi-1-i)*grad;
             } else if (i == hi) {
               amrex::Real grad = data(hi-1,j,k,n) - data(hi-2,j,k,n);
               data(i,j,k,n) = data(hi-1,j,k,n) - (hi-1-i)*grad;
             }
           });
        } else if (IDIR == 1) {
           amrex::ParallelFor(b, ncomp, [=] AMREX_GPU_DEVICE (int i, int j, int k, int n) noexcept
           {
             if (j > hi) {
               amrex::Real grad = data(i,hi-1,k,n) - data(i,hi-2,k,n);
               data(i,j,k,n) = data(i,hi-1,k,n) - (hi-1-j)*grad; 
             } else if (j == hi) {
               amrex::Real grad = data(i,hi-1,k,n) - data(i,hi-2,k,n);
               data(i,j,k,n) = data(i,hi-1,k,n) - (hi-1-j)*grad;
             }
           });
        } else if (IDIR == 2) {
           amrex::ParallelFor(b, ncomp, [=] AMREX_GPU_DEVICE (int i, int j, int k, int n) noexcept
           {
             if (k > hi) {
               amrex::Real grad = data(i,j,hi-1,n) - data(i,j,hi-2,n);
               data(i,j,k,n) = data(i,j,hi-1,n) - (hi-1-k)*grad;
             } else if (k == hi) {
               amrex::Real grad = data(i,j,hi-1,n) - data(i,j,hi-2,n);
               data(i,j,k,n) = data(i,j,hi-1,n) - (hi-1-k)*grad; 
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

  template<int IDIR, BCBound Bound> struct reflect_even;

  template<int IDIR>
  struct reflect_even<IDIR, BCBound::lower> {
    using self = reflect_even<IDIR, BCBound::lower>;

    static inline void apply_cell_based(const amrex::Geometry& geom, const Box& b, Array4<Real> data, int ncomp=1) {
      int lo = geom.Domain().smallEnd(IDIR);
 
      if(b.smallEnd(IDIR) < lo) {
        if (IDIR == 0) {
          amrex::ParallelFor(b, ncomp, [=] AMREX_GPU_DEVICE (int i, int j, int k, int n) noexcept
          {
            // for zone center variables, if mirrored with zone (lo), then lo-1 is the mirror of lo+1 (see below)
            // |_______|_______|________|________|
            //   lo-2     lo-1    lo      lo+1
            //
            if (i < lo) {
              data(i,j,k,n) = data(2*lo-i,j,k,n);
            }
          });
        } else if (IDIR == 1) {
          amrex::ParallelFor(b, ncomp, [=] AMREX_GPU_DEVICE (int i, int j, int k, int n) noexcept
          {
            if (j < lo) {
              data(i,j,k,n) = data(i,2*lo-j,k,n);
            }
          });
        } else if (IDIR == 2) {
          amrex::ParallelFor(b, ncomp, [=] AMREX_GPU_DEVICE (int i, int j, int k, int n) noexcept
          {
            if (k < lo) {
              data(i,j,k,n) = data(i,j,2*lo-k,n);
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
            // for face center variables, if mirrored with zone (lo), then lo-1 is the mirror of lo+1 (see below)
            //  |_______|_______|________|________|
            // lo-2    lo-1    lo      lo+1     lo+2
            //      
            if (i < lo) {
              data(i,j,k,n) = data(2*lo-i,j,k,n);
            } else if (i == lo) {
              data(i,j,k,n) = 0.;
            }
           });
        } else if (IDIR == 1) {
           amrex::ParallelFor(b, ncomp, [=] AMREX_GPU_DEVICE (int i, int j, int k, int n) noexcept
           {
            if (j < lo) {
              data(i,j,k,n) = data(i,2*lo-j,k,n);
            } else if (j == lo) {
              data(i,j,k,n) = 0.;
            }
           });
        } else if (IDIR == 2) {
           amrex::ParallelFor(b, ncomp, [=] AMREX_GPU_DEVICE (int i, int j, int k, int n) noexcept
           {
            if (k < lo) {
              data(i,j,k,n) = data(i,j,2*lo-k,n);
            } else if (k == lo) {
              data(i,j,k,n) = 0.;
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
  struct reflect_even<IDIR, BCBound::upper> {
    using self = reflect_even<IDIR, BCBound::upper>;

    static inline void apply_cell_based(const amrex::Geometry& geom, const Box& b, Array4<Real> data, int ncomp=1) {
      int hi = geom.Domain().bigEnd(IDIR);
 
      if(b.bigEnd(IDIR) > hi) {
        if (IDIR == 0) {
            // for zone center variables, if mirrored with zone (hi), then hi-1 is the mirror of hi+1 (see below)
            // |_______|_______|________|________|________|
            //    hi-2    hi-1    hi       hi+1     hi+2
            //
          amrex::ParallelFor(b, ncomp, [=] AMREX_GPU_DEVICE (int i, int j, int k, int n) noexcept
          {
            if (i > hi) {
              data(i,j,k,n) = data(2*hi-i,j,k,n);
            }
          });
        } else if (IDIR == 1) {
          amrex::ParallelFor(b, ncomp, [=] AMREX_GPU_DEVICE (int i, int j, int k, int n) noexcept
          {
            if (j > hi) {
              data(i,j,k,n) = data(i,2*hi-j,k,n);
            }
          });
        } else if (IDIR == 2) {
          amrex::ParallelFor(b, ncomp, [=] AMREX_GPU_DEVICE (int i, int j, int k, int n) noexcept
          {
            if (k > hi) {
              data(i,j,k,n) = data(i,j,2*hi-k,n);
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
            // for face center variables, if mirrored with zone (hi), then hi-1 is the mirror of hi+1 (see below)
            //  |_______|_______|________|________|________|
            // hi-2    hi-1    hi       hi+1     hi+2    hi+3
            //
             if (i > hi) {
               data(i,j,k,n) = data(2*hi-i,j,k,n);
             } else if (i == hi) {
               data(i,j,k,n) = 0.;
             }
           });
        } else if (IDIR == 1) {
           amrex::ParallelFor(b, ncomp, [=] AMREX_GPU_DEVICE (int i, int j, int k, int n) noexcept
           {
             if (j > hi) {
               data(i,j,k,n) = data(i,2*hi-j,k,n);
             } else if (j == hi) {
               data(i,j,k,n) = 0.;
             }
           });
        } else if (IDIR == 2) {
           amrex::ParallelFor(b, ncomp, [=] AMREX_GPU_DEVICE (int i, int j, int k, int n) noexcept
           {
             if (k > hi) {
               data(i,j,k,n) = data(i,j,2*hi-k,n);
             } else if (k == hi) {
               data(i,j,k,n) = 0.;
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


  template<int IDIR, BCBound Bound> struct reflect_odd;

  template<int IDIR>
  struct reflect_odd<IDIR, BCBound::lower> {
    using self = reflect_odd<IDIR, BCBound::lower>;

    static inline void apply_cell_based(const amrex::Geometry& geom, const Box& b, Array4<Real> data, int ncomp=1) {
      int lo = geom.Domain().smallEnd(IDIR);
 
      if(b.smallEnd(IDIR) < lo) {
        if (IDIR == 0) {
          amrex::ParallelFor(b, ncomp, [=] AMREX_GPU_DEVICE (int i, int j, int k, int n) noexcept
          {
            // for zone center variables, if mirrored with zone (hi), then hi-1 is the mirror of hi+1 (see below)
            //   |_______|_______|________|________|________|
            //      lo-2    lo-1     lo       lo+1     lo+2
            //
            if (i < lo) {
              data(i,j,k,n) = -data(2*lo-i,j,k,n);
            }
          });
        } else if (IDIR == 1) {
          amrex::ParallelFor(b, ncomp, [=] AMREX_GPU_DEVICE (int i, int j, int k, int n) noexcept
          {
            if (j < lo) {
              data(i,j,k,n) = -data(i,2*lo-j,k,n);
            }
          });
        } else if (IDIR == 2) {
          amrex::ParallelFor(b, ncomp, [=] AMREX_GPU_DEVICE (int i, int j, int k, int n) noexcept
          {
            if (k < lo) {
              data(i,j,k,n) = -data(i,j,2*lo-k,n);
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
            // for face center variables, if mirrored with zone (lo), then lo-1 is the mirror of lo+1 (see below)
            //    |_______|_______|________|________|________|
            //   lo-2    lo-1    lo       lo+1     lo+2    lo+3
            //
            if (i < lo) {
              data(i,j,k,n) = -data(2*lo-i,j,k,n);
            } else if (i == lo) {
              data(i,j,k,n) = 0.;
            }
           });
        } else if (IDIR == 1) {
           amrex::ParallelFor(b, ncomp, [=] AMREX_GPU_DEVICE (int i, int j, int k, int n) noexcept
           {
            if (j < lo) {
              data(i,j,k,n) = -data(i,2*lo-j,k,n);
            } else if (j == lo) {
              data(i,j,k,n) = 0.;
            }
           });
        } else if (IDIR == 2) {
           amrex::ParallelFor(b, ncomp, [=] AMREX_GPU_DEVICE (int i, int j, int k, int n) noexcept
           {
            if (k < lo) {
              data(i,j,k,n) = -data(i,j,2*lo-k,n);
            } else if (k == lo) {
              data(i,j,k,n) = 0.;
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
  struct reflect_odd<IDIR, BCBound::upper> {
    using self = reflect_odd<IDIR, BCBound::upper>;

    static inline void apply_cell_based(const amrex::Geometry& geom, const Box& b, Array4<Real> data, int ncomp=1) {
      int hi = geom.Domain().bigEnd(IDIR);
 
      if(b.bigEnd(IDIR) > hi) {
        if (IDIR == 0) {
          amrex::ParallelFor(b, ncomp, [=] AMREX_GPU_DEVICE (int i, int j, int k, int n) noexcept
          {
            // for zone center variables, if mirrored with zone (hi), then hi-1 is the mirror of hi+1 (see below)
            //    |_______|_______|________|________|________|
            //       hi-2    hi-1    hi       hi+1     hi+2
            //
            if (i > hi) {
              data(i,j,k,n) = -data(2*hi-i,j,k,n);
            }
          });
        } else if (IDIR == 1) {
          amrex::ParallelFor(b, ncomp, [=] AMREX_GPU_DEVICE (int i, int j, int k, int n) noexcept
          {
            if (j > hi) {
              data(i,j,k,n) = -data(i,2*hi-j,k,n);
            }
          });
        } else if (IDIR == 2) {
          amrex::ParallelFor(b, ncomp, [=] AMREX_GPU_DEVICE (int i, int j, int k, int n) noexcept
          {
            if (k > hi) {
              data(i,j,k,n) = -data(i,j,2*hi-k,n);
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
            // for face center variables, if mirrored with zone (hi), then hi-1 is the mirror of hi+1 (see below)
            //   |_______|_______|________|________|________|
            //  hi-2    hi-1    hi       hi+1     hi+2
            //
             if (i > hi) {
               data(i,j,k,n) = -data(2*hi-i,j,k,n);
             } else if (i == hi) {
               data(i,j,k,n) = 0.;
             }
           });
        } else if (IDIR == 1) {
           amrex::ParallelFor(b, ncomp, [=] AMREX_GPU_DEVICE (int i, int j, int k, int n) noexcept
           {
             if (j > hi) {
               data(i,j,k,n) = -data(i,2*hi-j,k,n);
             } else if (j == hi) {
               data(i,j,k,n) = 0.;
             }
           });
        } else if (IDIR == 2) {
           amrex::ParallelFor(b, ncomp, [=] AMREX_GPU_DEVICE (int i, int j, int k, int n) noexcept
           {
             if (k > hi) {
               data(i,j,k,n) = -data(i,j,2*hi-k,n);
             } else if (k == hi) {
               data(i,j,k,n) = 0.;
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


  template<int IDIR, BCBound Bound> struct hoextrap;

  template<int IDIR>
  struct hoextrap<IDIR, BCBound::lower> {
    using self = hoextrap<IDIR, BCBound::lower>;

    static inline void apply_cell_based(const amrex::Geometry& geom, const Box& b, Array4<Real> data, int ncomp=1) {
      int lo = geom.Domain().smallEnd(IDIR);
 
      if(b.smallEnd(IDIR) < lo) {
        if (IDIR == 0) {
          amrex::ParallelFor(b, ncomp, [=] AMREX_GPU_DEVICE (int i, int j, int k, int n) noexcept
          {
            // NOTE: this is 2nd order of accuracy
            // using third order 3-points formula fprime = (-3f(x)+4f(x+h))-f(x+2h))/2h
            //
            if (i < lo) {
              amrex::Real grad = (-3.*data(lo,j,k,n)+4.*data(lo+1,j,k,n)-data(lo+2,j,k,n))/2.;
              data(i,j,k,n) = data(lo,j,k,n) - (lo-i)*grad;
            }
          });
        } else if (IDIR == 1) {
          amrex::ParallelFor(b, ncomp, [=] AMREX_GPU_DEVICE (int i, int j, int k, int n) noexcept
          {
            if (j < lo) {
              amrex::Real grad = (-3.*data(i,lo,k,n)+4.*data(i,lo+1,k,n)-data(i,lo+2,k,n))/2.;
              data(i,j,k,n) = data(i,lo,k,n) - (lo-j)*grad;
            }
          });
        } else if (IDIR == 2) {
          amrex::ParallelFor(b, ncomp, [=] AMREX_GPU_DEVICE (int i, int j, int k, int n) noexcept
          {
            if (k < lo) {
              amrex::Real grad = (-3.*data(i,j,lo,n)+4.*data(i,j,lo+1,n)-data(i,j,lo+2,n))/2.;
              data(i,j,k,n) = data(i,j,lo,n) - (lo-k)*grad;
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
            // NOTE: this is 2nd order of accuracy
            //  using third order 3-points formula fprime = (-3f(x)+4f(x+h))-f(x+2h))/2h
            //
            if (i < lo) {
              amrex::Real grad = (-3.*data(lo+1,j,k,n)+4.*data(lo+2,j,k,n)-data(lo+3,j,k,n))/2.;
              data(i,j,k,n) = data(lo+1,j,k,n) - (lo-i+1)*grad;
            } else if (i == lo) {
              amrex::Real grad = (-3.*data(lo+1,j,k,n)+4.*data(lo+2,j,k,n)-data(lo+3,j,k,n))/2.;
              data(i,j,k,n) = data(lo+1,j,k,n) - grad;
            }
           });
        } else if (IDIR == 1) {
           amrex::ParallelFor(b, ncomp, [=] AMREX_GPU_DEVICE (int i, int j, int k, int n) noexcept
           {
            if (j < lo) {
              amrex::Real grad = (-3.*data(i,lo+1,k,n)+4.*data(i,lo+2,k,n)-data(i,lo+3,k,n))/2.;
              data(i,j,k,n) = data(i,lo+1,k,n) - (lo-j+1)*grad;
            } else if (j == lo) {
              amrex::Real grad = (-3.*data(i,j,lo+1,n)+4.*data(i,j,lo+2,n)-data(i,j,lo+3,n))/2.;
              data(i,j,k,n) = data(i,lo+1,k,n)-grad;
            }
           });
        } else if (IDIR == 2) {
           amrex::ParallelFor(b, ncomp, [=] AMREX_GPU_DEVICE (int i, int j, int k, int n) noexcept
           {
            if (k < lo) {
              amrex::Real grad = (-3.*data(i,j,lo+1,n)+4.*data(i,j,lo+2,n)-data(i,j,lo+3,n))/2.;
              data(i,j,k,n) = data(i,j,lo+1,n) - (lo-k+1)*grad;
            } else if (k == lo) {
              amrex::Real grad = (-3.*data(i,j,lo+1,n)+4.*data(i,j,lo+2,n)-data(i,j,lo+3,n))/2.;
              data(i,j,k,n) = data(i,j,lo+1,n) - grad;
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
  struct hoextrap<IDIR, BCBound::upper> {
    using self = hoextrap<IDIR, BCBound::upper>;

    static inline void apply_cell_based(const amrex::Geometry& geom, const Box& b, Array4<Real> data, int ncomp=1) {
      int hi = geom.Domain().bigEnd(IDIR);
 
      if(b.bigEnd(IDIR) > hi) {
        if (IDIR == 0) {
          amrex::ParallelFor(b, ncomp, [=] AMREX_GPU_DEVICE (int i, int j, int k, int n) noexcept
          {
            // NOTE: this is 2nd order of accuracy
            //  using third order 3-points formula fprime = (-3f(x-2h)+4f(x-h))-f(x))/2h
            //
            if (i > hi) {
              amrex::Real grad = (-3.*data(hi-2,j,k,n)+4.*data(hi-1,j,k,n)-data(hi,j,k,n))/2.;
              data(i,j,k,n) = data(hi,j,k,n) - (hi-i)*grad;
            }
          });
        } else if (IDIR == 1) {
          amrex::ParallelFor(b, ncomp, [=] AMREX_GPU_DEVICE (int i, int j, int k, int n) noexcept
          {
            if (j > hi) {
              amrex::Real grad = (-3.*data(i,hi-2,k,n)+4.*data(i,hi-1,k,n)-data(i,hi,k,n))/2.;
              data(i,j,k,n) = data(i,hi,k,n) - (hi-j)*grad;
            }
          });
        } else if (IDIR == 2) {
          amrex::ParallelFor(b, ncomp, [=] AMREX_GPU_DEVICE (int i, int j, int k, int n) noexcept
          {
            if (k > hi) {
              amrex::Real grad = (-3.*data(i,j,hi-2,n)+4.*data(i,j,hi-1,n)-data(i,j,hi,n))/2.;
              data(i,j,k,n) = data(i,j,hi,n) - (hi-k)*grad;
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
            // NOTE: this is 2nd order of accuracy
            //  using third order 3-points formula fprime = (-3f(x-2h)+4f(x-h))-f(x))/2h
            //
             if (i > hi) {
               amrex::Real grad = (-3.*data(hi-3,j,k,n)+4.*data(hi-2,j,k,n)-data(hi-1,j,k,n))/2.;
               data(i,j,k,n) = data(hi-1,j,k,n) - (hi-i-1)*grad;
             } else if (i == hi) {
               amrex::Real grad = (-3.*data(hi-3,j,k,n)+4.*data(hi-2,j,k,n)-data(hi-1,j,k,n))/2.;
               data(i,j,k,n) = data(hi-1,j,k,n) - (hi-i-1)*grad;
             }
           });
        } else if (IDIR == 1) {
           amrex::ParallelFor(b, ncomp, [=] AMREX_GPU_DEVICE (int i, int j, int k, int n) noexcept
           {
             if (j > hi) {
               amrex::Real grad = (-3.*data(i,hi-3,k,n)+4.*data(i,hi-2,k,n)-data(i,hi-1,k,n))/2.;
               data(i,j,k,n) = data(i,hi-1,k,n) - (hi-j-1)*grad;
             } else if (j == hi) {
               amrex::Real grad = (-3.*data(i,hi-3,k,n)+4.*data(i,hi-2,k,n)-data(i,hi-1,k,n))/2.;
               data(i,j,k,n) = data(i,hi-1,k,n) - (hi-j-1)*grad;
             }
           });
        } else if (IDIR == 2) {
           amrex::ParallelFor(b, ncomp, [=] AMREX_GPU_DEVICE (int i, int j, int k, int n) noexcept
           {
             if (k > hi) {
               amrex::Real grad = (-3.*data(i,j,hi-3,n)+4.*data(i,j,hi-2,n)-data(i,j,hi-1))/2.;
               data(i,j,k,n) = data(i,j,hi-1,n) - (hi-k-1)*grad;
             } else if (k == hi) {
               amrex::Real grad = (-3.*data(i,j,hi-3,n)+4.*data(i,j,hi-2,n)-data(i,j,hi-1))/2.;
               data(i,j,k,n) = data(i,j,hi-1,n) - (hi-k-1)*grad;
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
  inline void brian_reflect_odd<BCBound::lower>(const amrex::Geometry& geom, const Box& b, int lo, Array4<Real> data, int idir, int ncomp) {
    if(b.smallEnd(idir) < lo) {
      if (idir == 0) {
        amrex::ParallelFor(b, ncomp, [=] AMREX_GPU_DEVICE (int i, int j, int k, int n) noexcept
        {
            // for zone center variables, if mirrored with zone interface (lo-1, lo), then lo-1 is the mirror of lo (see below)
            //     |_______|_______|________|________|________|
            //        lo-2    lo-1    lo       lo+1     lo+2
            // 
            //
          if (i < lo) {
            data(i,j,k,n) = -data(lo+(lo-i)-1,j,k,n);
          }
        });
      } else if (idir == 1) {
        amrex::ParallelFor(b, ncomp, [=] AMREX_GPU_DEVICE (int i, int j, int k, int n) noexcept
        {
          if (j < lo) {
            data(i,j,k,n) = -data(i,lo+(lo-j)-1,k,n);
          }
        });
      } else if (idir == 2) {
        amrex::ParallelFor(b, ncomp, [=] AMREX_GPU_DEVICE (int i, int j, int k, int n) noexcept
        {
          if (k < lo) {
            data(i,j,k,n) = -data(i,j,lo+(lo-k)-1,n);
          }
        });
      }
    }
  }

  template<>
  inline void brian_reflect_odd<BCBound::upper>(const amrex::Geometry& geom, const Box& b, int hi, Array4<Real> data, int idir, int ncomp) {
      if(b.bigEnd(idir) > hi) {
          if (idir == 0) {
              amrex::ParallelFor(b, ncomp, [=] AMREX_GPU_DEVICE (int i, int j, int k, int n) noexcept
              {
            // for zone center variables, if mirrored with zone interface (hi, hi+1), then hi-1 is the mirror of hi (see below)
            //    |_______|_______|________|________|________|
            //       hi-2    hi-1    hi       hi+1     hi+2
            //
                  if (i > hi) {
                      data(i,j,k,n) = -data(hi-(i-hi)+1,j,k,n);
                  }
              });
          } else if (idir == 1) {
              amrex::ParallelFor(b, ncomp, [=] AMREX_GPU_DEVICE (int i, int j, int k, int n) noexcept
              {
                  if (j > hi) {
                      data(i,j,k,n) =  -data(i,hi-(j-hi)+1,k,n);
                  }
              });
          } else if (idir == 2) {
              amrex::ParallelFor(b, ncomp, [=] AMREX_GPU_DEVICE (int i, int j, int k, int n) noexcept
              {
                  if (k > hi) {
                      data(i,j,k,n) = -data(i,j,hi-(k-hi)+1,n);
                  }
              });
          }
      }
  }

  template<>
  inline void brian_reflect_odd_face_based<BCBound::lower>(const amrex::Geometry& geom, const Box& b, int lo, Array4<Real> data, int idir, int ncomp) {
    if(b.smallEnd(idir) < lo) {
      if (idir == 0) {
        amrex::ParallelFor(b, ncomp, [=] AMREX_GPU_DEVICE (int i, int j, int k, int n) noexcept
        {
            // for zone center variables, if mirrored with zone (hi), then hi-1 is the mirror of hi+1 (see below)
            //   |_______|_______|________|________|________|
            //  lo-2    lo-1    lo       lo+1     lo+2     lo+3
            //
          if (i < lo) {
            data(i,j,k,n) = -data(lo+(lo-i),j,k,n);
          } else if (i == lo) {
            data(i,j,k,n) = 0.0;
	  }
        });
      } else if (idir == 1) {
        amrex::ParallelFor(b, ncomp, [=] AMREX_GPU_DEVICE (int i, int j, int k, int n) noexcept
        {
          if (j < lo) {
            data(i,j,k,n) = -data(i,lo+(lo-j),k,n);
          } else if (j == lo) {
            data(i,j,k,n) = 0.0;
          }
          
        });
      } else if (idir == 2) {
        amrex::ParallelFor(b, ncomp, [=] AMREX_GPU_DEVICE (int i, int j, int k, int n) noexcept
        {
          if (k < lo) {
            data(i,j,k,n) = -data(i,j,lo+(lo-k),n);
          } else if (k == lo) {
            data(i,j,k,n) = 0.0;
	  }
        });
      }
    }
  }

  template<>
  inline void brian_reflect_odd_face_based<BCBound::upper>(const amrex::Geometry& geom, const Box& b, int hi, Array4<Real> data, int idir, int ncomp) {
      if(b.bigEnd(idir) > hi) {
          if (idir == 0) {
              amrex::ParallelFor(b, ncomp, [=] AMREX_GPU_DEVICE (int i, int j, int k, int n) noexcept
              {
            // for zone center variables, if mirrored with zone (hi), then hi-1 is the mirror of hi+1 (see below)
            //     |_______|_______|________|________|________|
            //    hi-2    hi-1    hi       hi+1     hi+2     hi+3
            //
                  if (i > hi) {
                      data(i,j,k,n) = -data(hi-(i-hi),j,k,n);
                  } else if (i == hi) {
                      data(i,j,k,n) = 0.0;
                  }
              });
          } else if (idir == 1) {
              amrex::ParallelFor(b, ncomp, [=] AMREX_GPU_DEVICE (int i, int j, int k, int n) noexcept
              {
                  if (j > hi) {
                      data(i,j,k,n) =  -data(i,hi-(j-hi),k,n);
                  } else if (j == hi) {
                      data(i,j,k,n) = 0.0;
                  }
              });
          } else if (idir == 2) {
              amrex::ParallelFor(b, ncomp, [=] AMREX_GPU_DEVICE (int i, int j, int k, int n) noexcept
              {
                  if (k > hi) {
                      data(i,j,k,n) = -data(i,j,hi-(k-hi),n);
                  } else if (k == hi) {
                      data(i,j,k,n) = 0.0;
                  }
              });
          }
      }
  }

  template<>
  inline void wall_scalars<BCBound::lower>(const amrex::Geometry& geom, const Box& b, int lo, Array4<Real> data, int idir, int ncomp) {
    amrex::IntVect type = b.type();
    if (type == amrex::IntVect::Zero) {

      /*  Real gravity = use_gravity? CONST_GRAV: 0.0;
       const    Array<Real,AMREX_SPACEDIM> grav{0.0, 0.0, gravity};
    //  Note that CONST_GRAV is positive, but the vertical component of grav is negative //
     const GpuArray<Real,AMREX_SPACEDIM> grav_gpu{grav[0], grav[1], grav[2]};  */

    const auto dxarray = geom.CellSizeArray();


    if(b.smallEnd(idir) < lo) {
      if (idir == 0) {
         amrex::ParallelFor(b, ncomp, [=] AMREX_GPU_DEVICE (int i, int j, int k, int n) noexcept
         {
          if (i < lo) {
            if (n == Rho_comp) {
	      //  Assume Neumann-type BCs for density (zero gradient, and extrapolation of value beyond first ghost point //
              data(i,j,k,n) = data(lo,j,k,n);
	    }
            else if (n == RhoTheta_comp) {
	      // Ideally, the BCs for RhoTheta would be determined by the density BCs, assuming hydrostaticity must hold at wall.
	      // However, this would require passing the user-specified value of use_gravity into this scope.
              // Barring that, we will simply apply Neuman type BCs for RhoTheta as well.
              // This currently only impacts the vertical momentum equation at the extreme points.
              // So there should be no problems as long as the w = 0 is enforced on the wall immediately after the w tendency update.  BJG
               

              /* Real plo = getPgivenRTh(data(lo,j,k,n));
              Real plom1 = (-grav[idir] * dxarray[idir] / 2.0) * (data(lo,j,k,Rho_comp) + data(i,j,k,Rho_comp)) + plo;
              data(i,j,k,n) = getRhoThetagivenP(plom1);  */
              data(i,j,k,n) = data(lo,j,k,n);
            }
            else {
	      // For any passive scalars, also assume zero gradient //
              data(i,j,k,n) = data(lo,j,k,n);
            }
          }
         });
      } else if (idir == 1) {
         amrex::ParallelFor(b, ncomp, [=] AMREX_GPU_DEVICE (int i, int j, int k, int n) noexcept
         {
          if (j < lo) {
            if (n == Rho_comp) {
              data(i,j,k,n) = data(i,lo,k,n);
	    }
            else if (n == RhoTheta_comp) {
	      /*  Real plo = getPgivenRTh(data(i,lo,k,n));
              Real plom1 = (-grav[idir] * dxarray[idir] / 2.0) * (data(i,lo,k,Rho_comp) + data(i,j,k,Rho_comp)) + plo;
              data(i,j,k,n) = getRhoThetagivenP(plom1);  */
	      data(i,j,k,n) = data(i,lo,k,n); 
            }
            else {
              data(i,j,k,n) = data(i,lo,k,n);
            }
          }
         });
      } else if (idir == 2) {
         amrex::ParallelFor(b, ncomp, [=] AMREX_GPU_DEVICE (int i, int j, int k, int n) noexcept
         {
          if (k < lo) {
            if (n == Rho_comp) {
              data(i,j,k,n) = data(i,j,lo,n);
	    }
            else if (n == RhoTheta_comp) {
              /* Real plo = getPgivenRTh(data(i,j,lo,n));
              Real plom1 = (-grav[idir] * dxarray[idir] / 2.0) * (data(i,j,lo,Rho_comp) + data(i,j,k,Rho_comp)) + plo;
              data(i,j,k,n) = getRhoThetagivenP(plom1); */
	      data(i,j,k,n) = data(i,j,lo,n);
            }
            else {
              data(i,j,k,n) = data(i,j,lo,n);
            }
          }
         });
      }
    }

    }
   
    else {
      // print some message like "shouldn't be here, this condition applies to scalars only   //
    }
  }




  template<>
  inline void wall_scalars<BCBound::upper>(const amrex::Geometry& geom, const Box& b, int hi, Array4<Real> data, int idir, int ncomp) {
      amrex::IntVect type = b.type();
      if (type == amrex::IntVect::Zero) {

          /*    Real gravity = solverChoice.use_gravity? CONST_GRAV: 0.0;
        const    Array<Real,AMREX_SPACEDIM> grav{0.0, 0.0, gravity};
        const GpuArray<Real,AMREX_SPACEDIM> grav_gpu{grav[0], grav[1], grav[2]}; */

          const auto dxarray = geom.CellSizeArray();



          if(b.bigEnd(idir) > hi) {
              if (idir == 0) {
                  amrex::ParallelFor(b, ncomp, [=] AMREX_GPU_DEVICE (int i, int j, int k, int n) noexcept
                  {
                      if (i > hi) {
                          if (n == Rho_comp) {
                              //  Assume Neumann-type BCs for density (zero gradient, and extrapolation of value beyond first ghost point //
                              data(i,j,k,n) = data(hi,j,k,n);
                          }
                          else if (n == RhoTheta_comp) {
                              // Ideally, the BCs for RhoTheta would be determined by the density BCs, assuming hydrostaticity must hold at wall.
                              // However, this would require passing the user-specified value of use_gravity into this scope.
                              // Barring that, we will simply apply Neuman type BCs for RhoTheta as well.
                              // This currently only impacts the vertical momentum equation at the extreme points.
                              // So there should be no problems as long as the w = 0 is enforced on the wall immediately after the w tendency update.  BJG

                              /* Real phi = getPgivenRTh(data(hi,j,k,n));
                              Real phip1 = (grav[idir] * dxarray[idir] / 2.0) * (data(hi,j,k,Rho_comp) + data(i,j,k,Rho_comp)) + phi;
                              data(i,j,k,n) = getRhoThetagivenP(phip1); */
                              data(i,j,k,n) = data(hi,j,k,n);
                          }
                          else {
                              // For any passive scalars, also assume zero gradient //
                              data(i,j,k,n) = data(hi,j,k,n);
                          }
                      }
                  });
              } else if (idir == 1) {
                  amrex::ParallelFor(b, ncomp, [=] AMREX_GPU_DEVICE (int i, int j, int k, int n) noexcept
                  {
                      if (j > hi) {
                          if (n == Rho_comp) {
                              data(i,j,k,n) = data(i,hi,k,n);
                          }
                          else if (n == RhoTheta_comp) {
                              /*  Real phi = getPgivenRTh(data(i,hi,k,n));
                                  Real phip1 = (grav[idir] * dxarray[idir] / 2.0) * (data(i,hi,k,Rho_comp) + data(i,j,k,Rho_comp)) + phi;
                                  data(i,j,k,n) = getRhoThetagivenP(phip1);  */
                              data(i,j,k,n) = data(i,hi,k,n);
                          }
                          else {
                              data(i,j,k,n) = data(i,hi,k,n);
                          }
                      }
                  });
              } else if (idir == 2) {
                  amrex::ParallelFor(b, ncomp, [=] AMREX_GPU_DEVICE (int i, int j, int k, int n) noexcept
                  {
                      if (k > hi) {
                          if (n == Rho_comp) {
                              data(i,j,k,n) = data(i,j,hi,n);
                          }
                          else if (n == RhoTheta_comp) {
                              /* Real phi = getPgivenRTh(data(i,j,hi,n));
                              Real phip1 = (grav[idir] * dxarray[idir] / 2.0) * (data(i,j,hi,Rho_comp) + data(i,j,k,Rho_comp)) + phi;
                              data(i,j,k,n) = getRhoThetagivenP(phip1); */
                              data(i,j,k,n) = data(i,j,hi,n);
                          }
                          else {
                              data(i,j,k,n) = data(i,j,hi,n);
                          }
                      }
                  });
              }
          }

      }

      else {
          // print some message like "shouldn't be here, this is for scalars only   //
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
