#ifndef ERF_BCBASE_H
#define ERF_BCBASE_H

#include "ERF.H"
#include "ABLMost.H"

using namespace amrex;

namespace math_bcs {

  enum BCBound { lower, upper };

  /****************************************************************************************************************/

  template<int IDIR, BCBound Bound> struct pc_extrap;

  template<int IDIR>
  struct pc_extrap<IDIR, BCBound::lower> {
    using self = pc_extrap<IDIR, BCBound::lower>;

    static inline void apply_cell_based(const amrex::Geometry& geom, const Box& b, Array4<Real> data, int ncomp=1) {
      int lo = geom.Domain().smallEnd(IDIR);
      if(b.smallEnd(IDIR) < lo) {
        if (IDIR == 0) {

          amrex::ParallelFor(b, ncomp, [=] AMREX_GPU_DEVICE (int i, int j, int k, int n) noexcept
          {
            if (i < lo) {
              data(i,j,k,n) = data(lo,j,k,n);
            }
          });
        } else if (IDIR == 1) {
          amrex::ParallelFor(b, ncomp, [=] AMREX_GPU_DEVICE (int i, int j, int k, int n) noexcept
          {
            if (j < lo) {
              data(i,j,k,n) = data(i,lo,k,n);
            }
          });
        } else if (IDIR == 2) {
          amrex::ParallelFor(b, ncomp, [=] AMREX_GPU_DEVICE (int i, int j, int k, int n) noexcept
          {
            if (k < lo) {
              data(i,j,k,n) = data(i,j,lo,n);
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
              data(i,j,k,n) = data(lo+1,j,k,n);
            }
           });
        } else if (IDIR == 1) {
           amrex::ParallelFor(b, ncomp, [=] AMREX_GPU_DEVICE (int i, int j, int k, int n) noexcept
           {
            if (j < lo) {
              data(i,j,k,n) = data(i,lo+1,k,n);
            }
           });
        } else if (IDIR == 2) {
           amrex::ParallelFor(b, ncomp, [=] AMREX_GPU_DEVICE (int i, int j, int k, int n) noexcept
           {
            if (k < lo) {
              data(i,j,k,n) = data(i,j,lo+1,n);
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
  struct pc_extrap<IDIR, BCBound::upper> {
    using self = pc_extrap<IDIR, BCBound::upper>;

    static inline void apply_cell_based(const amrex::Geometry& geom, const Box& b, Array4<Real> data, int ncomp=1) {
      int hi = geom.Domain().bigEnd(IDIR);

      if(b.bigEnd(IDIR) > hi) {
        if (IDIR == 0) {
          amrex::ParallelFor(b, ncomp, [=] AMREX_GPU_DEVICE (int i, int j, int k, int n) noexcept
          {
            if (i > hi) {
              data(i,j,k,n) = data(hi,j,k,n);
            }
          });
        } else if (IDIR == 1) {
          amrex::ParallelFor(b, ncomp, [=] AMREX_GPU_DEVICE (int i, int j, int k, int n) noexcept
          {
            if (j > hi) {
              data(i,j,k,n) = data(i,hi,k,n);
            }
          });
        } else if (IDIR == 2) {
          amrex::ParallelFor(b, ncomp, [=] AMREX_GPU_DEVICE (int i, int j, int k, int n) noexcept
          {
            if (k > hi) {
              data(i,j,k,n) = data(i,j,hi,n);
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
             if (i > (hi+1)) {
               data(i,j,k,n) = data(hi+1,j,k,n);
             }
           });
        } else if (IDIR == 1) {
           amrex::ParallelFor(b, ncomp, [=] AMREX_GPU_DEVICE (int i, int j, int k, int n) noexcept
           {
             if (j > (hi+1)) {
               data(i,j,k,n) = data(i,hi+1,k,n);
             }
           });
        } else if (IDIR == 2) {
           amrex::ParallelFor(b, ncomp, [=] AMREX_GPU_DEVICE (int i, int j, int k, int n) noexcept
           {
             if (k > (hi+1)) {
               data(i,j,k,n) = data(i,j,hi+1,n);
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

  /****************************************************************************************************************/

  template<int IDIR, BCBound Bound> struct lin_extrap;

  template<int IDIR>
  struct lin_extrap<IDIR, BCBound::lower> {
    using self = lin_extrap<IDIR, BCBound::lower>;

    static inline void apply_cell_based(const amrex::Geometry& geom, const Box& b, Array4<Real> data, int ncomp=1) {
      int lo = geom.Domain().smallEnd(IDIR);
      if(b.smallEnd(IDIR) < lo) {
        if (IDIR == 0) {

          amrex::ParallelFor(b, ncomp, [=] AMREX_GPU_DEVICE (int i, int j, int k, int n) noexcept
          {
            if (i < lo) {
              data(i,j,k,n) = (lo-i+1)*data(lo,j,k,n) - (lo-i)*data(lo+1,j,k,n);
            }
          });
        } else if (IDIR == 1) {
          amrex::ParallelFor(b, ncomp, [=] AMREX_GPU_DEVICE (int i, int j, int k, int n) noexcept
          {
            if (j < lo) {
              data(i,j,k,n) = (lo-j+1)*data(i,lo,k,n) - (lo-j)*data(i,lo+1,k,n);
            }
          });
        } else if (IDIR == 2) {
          amrex::ParallelFor(b, ncomp, [=] AMREX_GPU_DEVICE (int i, int j, int k, int n) noexcept
          {
            if (k < lo) {
              data(i,j,k,n) = (lo-k+1)*data(i,j,lo,n) - (lo-k)*data(i,j,lo+1,n);
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
              data(i,j,k,n) = (lo-i+2)*data(lo+1,j,k,n) - (lo-i+1)*data(lo+2,j,k,n);
            }
           });
        } else if (IDIR == 1) {
           amrex::ParallelFor(b, ncomp, [=] AMREX_GPU_DEVICE (int i, int j, int k, int n) noexcept
           {
            if (j < lo) {
              data(i,j,k,n) = (lo-j+2)*data(i,lo+1,k,n) - (lo-j+1)*data(i,lo+2,k,n);
            }
           });
        } else if (IDIR == 2) {
           amrex::ParallelFor(b, ncomp, [=] AMREX_GPU_DEVICE (int i, int j, int k, int n) noexcept
           {
            if (k < lo) {
              data(i,j,k,n) = (lo-k+2)*data(i,j,lo+1,n) - (lo-k+1)*data(i,j,lo+2,n);
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
  struct lin_extrap<IDIR, BCBound::upper> {
    using self = lin_extrap<IDIR, BCBound::upper>;

    static inline void apply_cell_based(const amrex::Geometry& geom, const Box& b, Array4<Real> data, int ncomp=1) {
      int hi = geom.Domain().bigEnd(IDIR);

      if(b.bigEnd(IDIR) > hi) {
        if (IDIR == 0) {
          amrex::ParallelFor(b, ncomp, [=] AMREX_GPU_DEVICE (int i, int j, int k, int n) noexcept
          {
            if (i > hi) {
              data(i,j,k,n) = (i-hi+1)*data(hi,j,k,n) -(i-hi)*data(hi-1,j,k,n);
            }
          });
        } else if (IDIR == 1) {
          amrex::ParallelFor(b, ncomp, [=] AMREX_GPU_DEVICE (int i, int j, int k, int n) noexcept
          {
            if (j > hi) {
              data(i,j,k,n) = (j-hi+1)*data(i,hi,k,n) - (j-hi)*data(i,hi-1,k,n);
            }
          });
        } else if (IDIR == 2) {
          amrex::ParallelFor(b, ncomp, [=] AMREX_GPU_DEVICE (int i, int j, int k, int n) noexcept
          {
            if (k > hi) {
              data(i,j,k,n) = (k-hi+1)*data(i,j,hi,n) - (k-hi)*data(i,j,hi-1,n);
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
             if (i > (hi+1)) {
               data(i,j,k,n) = (i-hi)*data(hi+1,j,k,n) - (i-hi-1)*data(hi,j,k,n);
             }
           });
        } else if (IDIR == 1) {
           amrex::ParallelFor(b, ncomp, [=] AMREX_GPU_DEVICE (int i, int j, int k, int n) noexcept
           {
             if (j > (hi+1)) {
               data(i,j,k,n) = (j-hi)*data(i,hi+1,k,n) - (j-hi-1)*data(i,hi,k,n);
             }
           });
        } else if (IDIR == 2) {
           amrex::ParallelFor(b, ncomp, [=] AMREX_GPU_DEVICE (int i, int j, int k, int n) noexcept
           {
             if (k > (hi+1)) {
               data(i,j,k,n) = (k-hi)*data(i,j,hi+1,n) - (k-hi-1)*data(i,j,hi,n);
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

  /****************************************************************************************************************/

  template<int IDIR, BCBound Bound> struct reflect_even;

  template<int IDIR>
  struct reflect_even<IDIR, BCBound::lower> {
    using self = reflect_even<IDIR, BCBound::lower>;

    static inline void apply_cell_based(const amrex::Geometry& geom, const Box& b, Array4<Real> data, int ncomp=1) {
      int lo = geom.Domain().smallEnd(IDIR);
            //
            // for zone center variables, if mirrored with zone interface (lo-1/2), then lo-1 is the mirror of lo (see below)
            // |_______|_______|________|________|
            //   lo-2     lo-1    lo      lo+1
            //
      if(b.smallEnd(IDIR) < lo) {
        if (IDIR == 0) {
          amrex::ParallelFor(b, ncomp, [=] AMREX_GPU_DEVICE (int i, int j, int k, int n) noexcept
          {
            if (i < lo) {
              data(i,j,k,n) = data(2*lo-i-1,j,k,n);
            }
          });
        } else if (IDIR == 1) {
          amrex::ParallelFor(b, ncomp, [=] AMREX_GPU_DEVICE (int i, int j, int k, int n) noexcept
          {
            if (j < lo) {
              data(i,j,k,n) = data(i,2*lo-j-1,k,n);
            }
          });
        } else if (IDIR == 2) {
          amrex::ParallelFor(b, ncomp, [=] AMREX_GPU_DEVICE (int i, int j, int k, int n) noexcept
          {
            if (k < lo) {
              data(i,j,k,n) = data(i,j,2*lo-k-1,n);
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

            // for face center variables, if mirrored with zone (lo), then lo-1 is the mirror of lo+1 (see below)
            //   |_______|_______|________|________|
            //  lo-2    lo-1    lo      lo+1     lo+2
            //
      if(b.smallEnd(IDIR) < lo) {
        if (IDIR == 0) {
           amrex::ParallelFor(b, ncomp, [=] AMREX_GPU_DEVICE (int i, int j, int k, int n) noexcept
           {
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
            //
            // NOTE: we should mirror at the zone interface to match the face-based boundary configuration
            // for zone center variables, if mirrored with zone interface (hi+1/2), then hi is the mirror of hi+1 (see below)
            // |_______|_______|________|________|________|
            //    hi-2    hi-1    hi       hi+1     hi+2
            //
      if(b.bigEnd(IDIR) > hi) {
        if (IDIR == 0) {
          amrex::ParallelFor(b, ncomp, [=] AMREX_GPU_DEVICE (int i, int j, int k, int n) noexcept
          {
            if (i > hi) {
              data(i,j,k,n) = data(2*hi-i+1,j,k,n);
            }
          });
        } else if (IDIR == 1) {
          amrex::ParallelFor(b, ncomp, [=] AMREX_GPU_DEVICE (int i, int j, int k, int n) noexcept
          {
            if (j > hi) {
              data(i,j,k,n) = data(i,2*hi-j+1,k,n);
            }
          });
        } else if (IDIR == 2) {
          amrex::ParallelFor(b, ncomp, [=] AMREX_GPU_DEVICE (int i, int j, int k, int n) noexcept
          {
            if (k > hi) {
              data(i,j,k,n) = data(i,j,2*hi-k+1,n);
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
            //
            // for face center variables, if mirrored with zone (hi+1), then hi is the mirror of hi+2 (see below)
            //  |_______|_______|________|________|________|
            // hi-2    hi-1    hi       hi+1     hi+2    hi+3
            //
      if(b.bigEnd(IDIR) > hi) {
        if (IDIR == 0) {
           amrex::ParallelFor(b, ncomp, [=] AMREX_GPU_DEVICE (int i, int j, int k, int n) noexcept
           {
             if (i > hi+1) {
               data(i,j,k,n) = data(2*hi-i+2,j,k,n);
             } else if (i == hi+1) {
               data(i,j,k,n) = 0.;
             }
           });
        } else if (IDIR == 1) {
           amrex::ParallelFor(b, ncomp, [=] AMREX_GPU_DEVICE (int i, int j, int k, int n) noexcept
           {
             if (j > hi+1) {
               data(i,j,k,n) = data(i,2*hi-j+2,k,n);
             } else if (j == hi+1) {
               data(i,j,k,n) = 0.;
             }
           });
        } else if (IDIR == 2) {
           amrex::ParallelFor(b, ncomp, [=] AMREX_GPU_DEVICE (int i, int j, int k, int n) noexcept
           {
             if (k > hi) {
               data(i,j,k,n) = data(i,j,2*hi-k+2,n);
             } else if (k == hi+1) {
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

  /****************************************************************************************************************/

  template<int IDIR, BCBound Bound> struct reflect_odd;

  template<int IDIR>
  struct reflect_odd<IDIR, BCBound::lower> {
    using self = reflect_odd<IDIR, BCBound::lower>;

    static inline void apply_cell_based(const amrex::Geometry& geom, const Box& b, Array4<Real> data, int ncomp=1) {
      int lo = geom.Domain().smallEnd(IDIR);
             //
             // for zone center variables, if mirrored with zone interface (lo-1/2), then lo-1 is the mirror of lo (see below)
             //   |_______|_______|________|________|________|
             //      lo-2    lo-1     lo       lo+1     lo+2
             //
      if(b.smallEnd(IDIR) < lo) {
        if (IDIR == 0) {
          amrex::ParallelFor(b, ncomp, [=] AMREX_GPU_DEVICE (int i, int j, int k, int n) noexcept
          {
            if (i < lo) {
              data(i,j,k,n) = -data(lo+(lo-i)-1,j,k,n);
            }
          });
        } else if (IDIR == 1) {
          amrex::ParallelFor(b, ncomp, [=] AMREX_GPU_DEVICE (int i, int j, int k, int n) noexcept
          {
            if (j < lo) {
              data(i,j,k,n) = -data(i,lo+(lo-j)-1,k,n);
            }
          });
        } else if (IDIR == 2) {
          amrex::ParallelFor(b, ncomp, [=] AMREX_GPU_DEVICE (int i, int j, int k, int n) noexcept
          {
            if (k < lo) {
              data(i,j,k,n) = -data(i,j,lo+(lo-k)-1,n);
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

            // for face center variables, if mirrored with zone (lo), then lo-1 is the mirror of lo+1 (see below)
            //    |_______|_______|________|________|________|
            //   lo-2    lo-1    lo       lo+1     lo+2    lo+3
            //
      if(b.smallEnd(IDIR) <= lo) {
        if (IDIR == 0) {
           amrex::ParallelFor(b, ncomp, [=] AMREX_GPU_DEVICE (int i, int j, int k, int n) noexcept
           {
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
            // for zone center variables, if mirrored with zone (hi+1/2), then hi is the mirror of hi+1 (see below)
            //    |_______|_______|________|________|________|
            //       hi-2    hi-1    hi       hi+1     hi+2
            //
      if(b.bigEnd(IDIR) > hi) {
        if (IDIR == 0) {
          amrex::ParallelFor(b, ncomp, [=] AMREX_GPU_DEVICE (int i, int j, int k, int n) noexcept
          {
            if (i > hi) {
              data(i,j,k,n) = -data(hi-(i-hi)+1,j,k,n);
            }
          });
        } else if (IDIR == 1) {
          amrex::ParallelFor(b, ncomp, [=] AMREX_GPU_DEVICE (int i, int j, int k, int n) noexcept
          {
            if (j > hi) {
              data(i,j,k,n) =  -data(i,hi-(j-hi)+1,k,n);
            }
          });
        } else if (IDIR == 2) {
          amrex::ParallelFor(b, ncomp, [=] AMREX_GPU_DEVICE (int i, int j, int k, int n) noexcept
          {
            if (k > hi) {
              data(i,j,k,n) = -data(i,j,hi-(k-hi)+1,n);
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
            // for face center variables, if mirrored with zone (hi+1), then hi is the mirror of hi+2 (see below)
            //   |_______|_______|________|________|________|
            //  hi-2    hi-1    hi       hi+1     hi+2
            //
      if(b.bigEnd(IDIR) > hi) {
        if (IDIR == 0) {
           amrex::ParallelFor(b, ncomp, [=] AMREX_GPU_DEVICE (int i, int j, int k, int n) noexcept
           {
             if (i > (hi+1)) {
               data(i,j,k,n) = -data(2*(hi+1)-i,j,k,n);
             } else if (i == (hi+1)) {
               data(i,j,k,n) = 0.;
             }
           });
        } else if (IDIR == 1) {
           amrex::ParallelFor(b, ncomp, [=] AMREX_GPU_DEVICE (int i, int j, int k, int n) noexcept
           {
             if (j > (hi+1)) {
               data(i,j,k,n) = -data(i,2*(hi+1)-j,k,n);
             } else if (j == (hi+1)) {
               data(i,j,k,n) = 0.;
             }
           });
        } else if (IDIR == 2) {
           amrex::ParallelFor(b, ncomp, [=] AMREX_GPU_DEVICE (int i, int j, int k, int n) noexcept
           {
             if (k > (hi+1)) {
               data(i,j,k,n) = -data(i,j,2*(hi+1)-k,n);
             } else if (k == (hi+1)) {
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

  /****************************************************************************************************************/

  template<int IDIR, BCBound Bound> struct dirichlet;

  template<int IDIR>
  struct dirichlet<IDIR, BCBound::lower> {
    using self = dirichlet<IDIR, BCBound::lower>;

    static inline void apply_cell_based(
      const amrex::Geometry& geom,
      const Box& b,
      Array4<Real> data,
      int ncomp=1,
      Real bcval=0.0
    ) {
      int lo = geom.Domain().smallEnd(IDIR);
      // for i < lo, (data(i) + data(lo)) / 2 == bcval

      //amrex::Print() << "apply_cell_based (lo=" << lo << ") "
      //               << bcval << " for BC " << IDIR << std::endl;

      if(b.smallEnd(IDIR) < lo) {
        if (IDIR == 0) {
          amrex::ParallelFor(b, ncomp, [=] AMREX_GPU_DEVICE (int i, int j, int k, int n) noexcept
          {
             if (i < lo) {
                  data(i,j,k,n) = -data(lo,j,k,n) + 2*bcval;
            }
          });
        } else if (IDIR == 1) {
          amrex::ParallelFor(b, ncomp, [=] AMREX_GPU_DEVICE (int i, int j, int k, int n) noexcept
          {
            if (j < lo) {
                  data(i,j,k,n) = -data(i,lo,k,n) + 2*bcval;
            }
          });
        } else if (IDIR == 2) {
          amrex::ParallelFor(b, ncomp, [=] AMREX_GPU_DEVICE (int i, int j, int k, int n) noexcept
          {
            if (k < lo) {
                  data(i,j,k,n) = -data(i,j,lo,n) + 2*bcval;
            }
          });
        }
      }
    }

    static inline void apply_face_based(
      const amrex::Geometry& geom,
      const Box& b,
      Array4<Real> data,
      int ncomp=1,
      Real bcval=0.0
    ) {
      int lo = geom.Domain().smallEnd(IDIR);
      // for i = lo, data(i) == bcval
      // for i < lo, bcval - data(i) == data(lo-(i-lo)) - bcval

      //amrex::Print() << "apply_face_based (lo=" << lo << ") "
      //               << bcval << " for BC " << IDIR << std::endl;

      if(b.smallEnd(IDIR) <= lo) {
        if (IDIR == 0) {
           amrex::ParallelFor(b, ncomp, [=] AMREX_GPU_DEVICE (int i, int j, int k, int n) noexcept
           {
            if (i < lo) {
              data(i,j,k,n) = 2*bcval - data(2*lo-i,j,k,n);
            } else if (i == lo) {
              data(i,j,k,n) = bcval;
            }
           });
        } else if (IDIR == 1) {
           amrex::ParallelFor(b, ncomp, [=] AMREX_GPU_DEVICE (int i, int j, int k, int n) noexcept
           {
            if (j < lo) {
              data(i,j,k,n) = 2*bcval - data(i,2*lo-j,k,n);
            } else if (j == lo) {
              data(i,j,k,n) = bcval;
            }
           });
        } else if (IDIR == 2) {
           amrex::ParallelFor(b, ncomp, [=] AMREX_GPU_DEVICE (int i, int j, int k, int n) noexcept
           {
            if (k < lo) {
              data(i,j,k,n) = 2*bcval - data(i,j,2*lo-k,n);
            } else if (k == lo) {
              data(i,j,k,n) = bcval;
            }
           });
        }
      }
    }

    template <bool val>
    static inline
    std::enable_if_t<val>
    apply(const amrex::Geometry& geom, const Box& b, Array4<Real> data, int ncomp=1, Real bcval=0.0) {
      self::apply_face_based(geom, b, data, ncomp, bcval);
    }

    template <bool val>
    static inline
    std::enable_if_t<!val>
    apply(const amrex::Geometry& geom, const Box& b, Array4<Real> data, int ncomp=1, Real bcval=0.0) {
      self::apply_cell_based(geom, b, data, ncomp, bcval);
    }

  };

  template<int IDIR>
  struct dirichlet<IDIR, BCBound::upper> {
    using self = dirichlet<IDIR, BCBound::upper>;

    static inline void apply_cell_based(
      const amrex::Geometry& geom,
      const Box& b,
      Array4<Real> data,
      int ncomp=1,
      Real bcval=0.0
    ) {
      int hi = geom.Domain().bigEnd(IDIR);
      // hi: # of cells - 1
      // for i > hi, (data(i) + data(hi)) / 2 == bcval

      //amrex::Print() << "apply_cell_based (hi=" << hi << ") "
      //               << bcval << " for BC " << IDIR << std::endl;

      if(b.bigEnd(IDIR) > hi) {
        if (IDIR == 0) {
          amrex::ParallelFor(b, ncomp, [=] AMREX_GPU_DEVICE (int i, int j, int k, int n) noexcept
          {
             if (i > hi) {
                  data(i,j,k,n) = -data(hi,j,k,n) + 2*bcval;
            }
          });
        } else if (IDIR == 1) {
          amrex::ParallelFor(b, ncomp, [=] AMREX_GPU_DEVICE (int i, int j, int k, int n) noexcept
          {
            if (j > hi) {
                  data(i,j,k,n) = -data(i,hi,k,n) + 2*bcval;
            }
          });
        } else if (IDIR == 2) {
          amrex::ParallelFor(b, ncomp, [=] AMREX_GPU_DEVICE (int i, int j, int k, int n) noexcept
          {
            if (k > hi) {
                  data(i,j,k,n) = -data(i,j,hi,n) + 2*bcval;
            }
          });
        }
      }
    }

    static inline void apply_face_based(
      const amrex::Geometry& geom,
      const Box& b,
      Array4<Real> data,
      int ncomp=1,
      Real bcval=0.0
    ) {
      int hi = geom.Domain().bigEnd(IDIR);
      // hi+1: index of boundary face
      // for i = hi+1, data(i) == bcval
      // for i > hi+1, data(i) - bcval == bcval - data(hi+1-(i-(hi+1)))
      //                                               hi+1-i+hi+1 = 2*hi-i+2

      //amrex::Print() << "apply_face_based (hi=" << hi << ") "
      //               << bcval << " for BC " << IDIR << std::endl;

      if(b.bigEnd(IDIR) >= hi+1) {
        if (IDIR == 0) {
           amrex::ParallelFor(b, ncomp, [=] AMREX_GPU_DEVICE (int i, int j, int k, int n) noexcept
           {
            if (i > hi+1) {
              data(i,j,k,n) = 2*bcval - data(2*hi-i+2,j,k,n);
            } else if (i == hi+1) {
              data(i,j,k,n) = bcval;
            }
           });
        } else if (IDIR == 1) {
           amrex::ParallelFor(b, ncomp, [=] AMREX_GPU_DEVICE (int i, int j, int k, int n) noexcept
           {
            if (j > hi+1) {
              data(i,j,k,n) = 2*bcval - data(i,2*hi-j+2,k,n);
            } else if (j == hi+1) {
              data(i,j,k,n) = bcval;
            }
           });
        } else if (IDIR == 2) {
           amrex::ParallelFor(b, ncomp, [=] AMREX_GPU_DEVICE (int i, int j, int k, int n) noexcept
           {
            if (k > hi+1) {
              data(i,j,k,n) = 2*bcval - data(i,j,2*hi-k+2,n);
            } else if (k == hi+1) {
              data(i,j,k,n) = bcval;
            }
           });
        }
      }
    }

    template <bool val>
    static inline
    std::enable_if_t<val>
    apply(const amrex::Geometry& geom, const Box& b, Array4<Real> data, int ncomp=1, Real bcval=0.0) {
      self::apply_face_based(geom, b, data, ncomp, bcval);
    }

    template <bool val>
    static inline
    std::enable_if_t<!val>
    apply(const amrex::Geometry& geom, const Box& b, Array4<Real> data, int ncomp=1, Real bcval=0.0) {
      self::apply_cell_based(geom, b, data, ncomp, bcval);
    }

  };

  /****************************************************************************************************************/

  template<int IDIR, BCBound Bound> struct sim_bc;

  template<int IDIR>
  struct sim_bc<IDIR, BCBound::lower> {
    using self = sim_bc<IDIR, BCBound::lower>;

    static inline void apply_cell_based(const amrex::Geometry& geom, const Box& b, Array4<Real> data, int ncomp=1)
    {
      const auto dxarray = geom.CellSizeArray();
      int lo = geom.Domain().smallEnd(IDIR);
      const auto zrough = 0.01;
             //
             // For similarity BCs, the condition on the tangential velocities is found by
             // setting the ghost point value to that which causes the shear stress at the lower boundary
             // (S_(-1/2)) equal to ustar / (2 * vonkarman * z) * ( (u,v) / (u * u + v * v) ^ 0.5 )
             //
             // Should work in viscous case, may need to recast for Smagorinksy, as well as for > 1 ghost point
             // Also only uses local velocity fields for now -- may need to recast taking average tangential
             // velocity along lower boundary into account, as in Moeng, Schumann, etc.

      if(b.smallEnd(IDIR) < lo) {
        if (IDIR == 0) {
          amrex::ParallelFor(b, ncomp, [=] AMREX_GPU_DEVICE (int i, int j, int k, int n) noexcept
          {
             if (i < lo) {
           // Assume distance between lo and lo-1 is dxarray(IDIR), and height of lo
           // above boundary is dxarray(IDIR) / 2.0
           // zrough is set to constant value 0.1 m here.  Should be passed in through input file.

                  data(i,j,k,n) = data(lo,j,k,n) * (1.0 - dxarray[IDIR] / (0.5 * dxarray[IDIR] * log(0.5 * dxarray[IDIR]/zrough)));
            }
          });
        } else if (IDIR == 1) {
          amrex::ParallelFor(b, ncomp, [=] AMREX_GPU_DEVICE (int i, int j, int k, int n) noexcept
          {
            if (j < lo) {
          data(i,j,k,n) = data(i,lo,k,n) * (1.0 - dxarray[IDIR] / (0.5 * dxarray[IDIR] * log(0.5 * dxarray[IDIR]/zrough)));
            }
          });
        } else if (IDIR == 2) {
          amrex::ParallelFor(b, ncomp, [=] AMREX_GPU_DEVICE (int i, int j, int k, int n) noexcept
          {
            if (k < lo) {
          data(i,j,k,n) = data(i,j,lo,n) * (1.0 - dxarray[IDIR] / (0.5 * dxarray[IDIR] * log(0.5 * dxarray[IDIR]/zrough)));
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

      // For similarity BC, conditions on normal velocity are same as with no slip

            // for face center variables, if mirrored with zone (lo), then lo-1 is the mirror of lo+1 (see below)
            //    |_______|_______|________|________|________|
            //   lo-2    lo-1    lo       lo+1     lo+2    lo+3
            //
      if(b.smallEnd(IDIR) <= lo) {
        if (IDIR == 0) {
           amrex::ParallelFor(b, ncomp, [=] AMREX_GPU_DEVICE (int i, int j, int k, int n) noexcept
           {
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
  struct sim_bc<IDIR, BCBound::upper> {
    using self = sim_bc<IDIR, BCBound::upper>;

    static inline void apply_cell_based(const amrex::Geometry& geom, const Box& b, Array4<Real> data, int ncomp=1)
    {
      const auto dxarray = geom.CellSizeArray();
      int hi = geom.Domain().bigEnd(IDIR);
      const auto zrough = 0.01;
             //
             // For similarity BCs, the condition on the tangential velocities is found by
             // setting the ghost point value to that which causes the shear stress at the upper boundary
             // (S_(hi+1/2)) equal to ustar / (2 * vonkarman * z) * ( (u,v) / (u * u + v * v) ^ 0.5 )
             //
             // Should work in viscous case, may need to recast for Smagorinksy, as well as for > 1 ghost point
             // Also only uses local velocity fields for now -- may need to recast taking average tangential
             // velocity along lower boundary into account, as in Moeng, Schumann, etc.

            // for zone center variables, if mirrored with zone (hi+1/2), then hi is the mirror of hi+1 (see below)
            //    |_______|_______|________|________|________|
            //       hi-2    hi-1    hi       hi+1     hi+2
            //
      if(b.bigEnd(IDIR) > hi) {
        if (IDIR == 0) {
          amrex::ParallelFor(b, ncomp, [=] AMREX_GPU_DEVICE (int i, int j, int k, int n) noexcept
          {
            if (i > hi) {
           // Assume distance between hi and hi+1 is dxarray(IDIR), and height of hi
           // below boundary is dxarray(IDIR) / 2.0
               // zrough is set to constant value 0.1 m here.  Should be passed in through input file.

                 data(i,j,k,n) = data(hi,j,k,n) * (1.0 - dxarray[IDIR] / (0.5 * dxarray[IDIR] * log(0.5 * dxarray[IDIR]/zrough)));
            }
          });
        } else if (IDIR == 1) {
          amrex::ParallelFor(b, ncomp, [=] AMREX_GPU_DEVICE (int i, int j, int k, int n) noexcept
          {
            if (j > hi) {
               data(i,j,k,n) = data(i,hi,k,n) * (1.0 - dxarray[IDIR] / (0.5 * dxarray[IDIR] * log(0.5 * dxarray[IDIR]/zrough)));
            }
          });
        } else if (IDIR == 2) {
          amrex::ParallelFor(b, ncomp, [=] AMREX_GPU_DEVICE (int i, int j, int k, int n) noexcept
          {
            if (k > hi) {
               data(i,j,k,n) = data(i,j,hi,n) * (1.0 - dxarray[IDIR] / (0.5 * dxarray[IDIR] * log(0.5 * dxarray[IDIR]/zrough)));
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
      //  For similarity BCs, condition on normal velocity is same as with no slip

            // for face center variables, if mirrored with zone (hi+1), then hi is the mirror of hi+2 (see below)
            //   |_______|_______|________|________|________|
            //  hi-2    hi-1    hi       hi+1     hi+2
            //
      if(b.bigEnd(IDIR) > hi) {
        if (IDIR == 0) {
           amrex::ParallelFor(b, ncomp, [=] AMREX_GPU_DEVICE (int i, int j, int k, int n) noexcept
           {
             if (i > (hi+1)) {
               data(i,j,k,n) = -data(2*(hi+1)-i,j,k,n);
             } else if (i == (hi+1)) {
               data(i,j,k,n) = 0.;
             }
           });
        } else if (IDIR == 1) {
           amrex::ParallelFor(b, ncomp, [=] AMREX_GPU_DEVICE (int i, int j, int k, int n) noexcept
           {
             if (j > (hi+1)) {
               data(i,j,k,n) = -data(i,2*(hi+1)-j,k,n);
             } else if (j == (hi+1)) {
               data(i,j,k,n) = 0.;
             }
           });
        } else if (IDIR == 2) {
           amrex::ParallelFor(b, ncomp, [=] AMREX_GPU_DEVICE (int i, int j, int k, int n) noexcept
           {
             if (k > (hi+1)) {
               data(i,j,k,n) = -data(i,j,2*(hi+1)-k,n);
             } else if (k == (hi+1)) {
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

  /****************************************************************************************************************/

  template<int IDIR, BCBound Bound> struct quad_extrap;

  template<int IDIR>
  struct quad_extrap<IDIR, BCBound::lower> {
    using self = quad_extrap<IDIR, BCBound::lower>;

    static inline void apply_cell_based(const amrex::Geometry& geom, const Box& b, Array4<Real> data, int ncomp=1) {
      int lo = geom.Domain().smallEnd(IDIR);

            // NOTE: this is 2nd order of accuracy
            // using third order 3-points formula fprime = (-3f(x)+4f(x+h))-f(x+2h))/2h
            //
      if(b.smallEnd(IDIR) < lo) {
        if (IDIR == 0) {
          amrex::ParallelFor(b, ncomp, [=] AMREX_GPU_DEVICE (int i, int j, int k, int n) noexcept
          {
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
            if (i <= lo) {
              amrex::Real grad = (-3.*data(lo+1,j,k,n)+4.*data(lo+2,j,k,n)-data(lo+3,j,k,n))/2.;
              data(i,j,k,n) = data(lo+1,j,k,n) - (lo-i+1)*grad;
            }
           });
        } else if (IDIR == 1) {
           amrex::ParallelFor(b, ncomp, [=] AMREX_GPU_DEVICE (int i, int j, int k, int n) noexcept
           {
            if (j <= lo) {
              amrex::Real grad = (-3.*data(i,lo+1,k,n)+4.*data(i,lo+2,k,n)-data(i,lo+3,k,n))/2.;
              data(i,j,k,n) = data(i,lo+1,k,n) - (lo-j+1)*grad;
            }
           });
        } else if (IDIR == 2) {
           amrex::ParallelFor(b, ncomp, [=] AMREX_GPU_DEVICE (int i, int j, int k, int n) noexcept
           {
            if (k <= lo) {
              amrex::Real grad = (-3.*data(i,j,lo+1,n)+4.*data(i,j,lo+2,n)-data(i,j,lo+3,n))/2.;
              data(i,j,k,n) = data(i,j,lo+1,n) - (lo-k+1)*grad;
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
  struct quad_extrap<IDIR, BCBound::upper> {
    using self = quad_extrap<IDIR, BCBound::upper>;

    static inline void apply_cell_based(const amrex::Geometry& geom, const Box& b, Array4<Real> data, int ncomp=1) {
      int hi = geom.Domain().bigEnd(IDIR);

      if(b.bigEnd(IDIR) > hi) {
        if (IDIR == 0) {
          amrex::ParallelFor(b, ncomp, [=] AMREX_GPU_DEVICE (int i, int j, int k, int n) noexcept
          {
            if (i > hi) {
              amrex::Real grad = (-3.*data(hi,j,k,n)+4.*data(hi+1,j,k,n)-data(hi,j,k,n))/2.;
              data(i,j,k,n) = data(hi,j,k,n) - (i-hi)*grad;
            }
          });
        } else if (IDIR == 1) {
          amrex::ParallelFor(b, ncomp, [=] AMREX_GPU_DEVICE (int i, int j, int k, int n) noexcept
          {
            if (j > hi) {
              amrex::Real grad = (-3.*data(i,hi,k,n)+4.*data(i,hi+1,k,n)-data(i,hi,k,n))/2.;
              data(i,j,k,n) = data(i,hi,k,n) - (j-hi)*grad;
            }
          });
        } else if (IDIR == 2) {
          amrex::ParallelFor(b, ncomp, [=] AMREX_GPU_DEVICE (int i, int j, int k, int n) noexcept
          {
            if (k > hi) {
              amrex::Real grad = (-3.*data(i,j,hi,n)+4.*data(i,j,hi+1,n)-data(i,j,hi+2,n))/2.;
              data(i,j,k,n) = data(i,j,hi,n) - (k-hi)*grad;
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
             if (i >= hi+1) {
               amrex::Real grad = (-3.*data(hi-2,j,k,n)+4.*data(hi-1,j,k,n)-data(hi,j,k,n))/2.;
               data(i,j,k,n) = data(hi,j,k,n) - (i-hi)*grad;
             }
           });
        } else if (IDIR == 1) {
           amrex::ParallelFor(b, ncomp, [=] AMREX_GPU_DEVICE (int i, int j, int k, int n) noexcept
           {
             if (j >= hi+1) {
               amrex::Real grad = (-3.*data(i,hi-2,k,n)+4.*data(i,hi-1,k,n)-data(i,hi,k,n))/2.;
               data(i,j,k,n) = data(i,hi,k,n) - (j-hi)*grad;
             }
           });
        } else if (IDIR == 2) {
           amrex::ParallelFor(b, ncomp, [=] AMREX_GPU_DEVICE (int i, int j, int k, int n) noexcept
           {
             if (k >= hi+1) {
               amrex::Real grad = (-3.*data(i,j,hi-2,n)+4.*data(i,j,hi-1,n)-data(i,j,hi))/2.;
               data(i,j,k,n) = data(i,j,hi,n) - (k-hi)*grad;
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

  /****************************************************************************************************************/

  template<int IDIR, BCBound Bound> struct wall_scalars;

  template<int IDIR>
  struct wall_scalars<IDIR, BCBound::lower> {
    using self = wall_scalars<IDIR, BCBound::lower>;

    static inline void apply(const amrex::Geometry& geom, const Box& b, Array4<Real> data, int ncomp=1) {
      int lo = geom.Domain().smallEnd(IDIR);

      amrex::IntVect type = b.type();
      if (type == amrex::IntVect::Zero) {

        /*  Real gravity = use_gravity? CONST_GRAV: 0.0;
        const    Array<Real,AMREX_SPACEDIM> grav{0.0, 0.0, gravity};
        //  Note that CONST_GRAV is positive, but the vertical component of grav is negative //
        const GpuArray<Real,AMREX_SPACEDIM> grav_gpu{grav[0], grav[1], grav[2]};  */

        // const auto dxarray = geom.CellSizeArray();

        if(b.smallEnd(IDIR) < lo) {
          if (IDIR == 0) {
             amrex::ParallelFor(b, ncomp, [=] AMREX_GPU_DEVICE (int i, int j, int k, int n) noexcept
             {
               if (i < lo) {
                 if (n == Rho_comp) {
                   //  Assume Neumann-type BCs for density (zero gradient, and extrapolation of value beyond first ghost point //
                   data(i,j,k,n) = data(lo,j,k,n);
                 } else if (n == RhoTheta_comp) {
                   // Ideally, the BCs for RhoTheta would be determined by the density BCs, assuming hydrostaticity must hold at wall.
                   // However, this would require passing the user-specified value of use_gravity into this scope.
                   // Barring that, we will simply apply Neuman type BCs for RhoTheta as well.
                   // This currently only impacts the vertical momentum equation at the extreme points.
                   // So there should be no problems as long as the w = 0 is enforced on the wall immediately after the w tendency update.  BJG

                   /* Real plo = getPgivenRTh(data(lo,j,k,n));
                   Real plom1 = (-grav[IDIR] * dxarray[IDIR] / 2.0) * (data(lo,j,k,Rho_comp) + data(i,j,k,Rho_comp)) + plo;
                   data(i,j,k,n) = getRhoThetagivenP(plom1);  */
                   data(i,j,k,n) = data(lo,j,k,n);
                 }
                 else {
                   // For any passive scalars, also assume zero gradient //
                   data(i,j,k,n) = data(lo,j,k,n);
                 }
               }
             });
          } else if (IDIR == 1) {
             amrex::ParallelFor(b, ncomp, [=] AMREX_GPU_DEVICE (int i, int j, int k, int n) noexcept
             {
               if (j < lo) {
                 if (n == Rho_comp) {
                   data(i,j,k,n) = data(i,lo,k,n);
                 } else if (n == RhoTheta_comp) {
                   /*  Real plo = getPgivenRTh(data(i,lo,k,n));
                   Real plom1 = (-grav[IDIR] * dxarray[IDIR] / 2.0) * (data(i,lo,k,Rho_comp) + data(i,j,k,Rho_comp)) + plo;
                   data(i,j,k,n) = getRhoThetagivenP(plom1);  */
                   data(i,j,k,n) = data(i,lo,k,n);
                 }
                 else {
                   data(i,j,k,n) = data(i,lo,k,n);
                 }
               }
             });
          } else if (IDIR == 2) {
             amrex::ParallelFor(b, ncomp, [=] AMREX_GPU_DEVICE (int i, int j, int k, int n) noexcept
             {
               if (k < lo) {
                 if (n == Rho_comp) {
                   data(i,j,k,n) = data(i,j,lo,n);
                 }
                 else if (n == RhoTheta_comp) {
                   /* Real plo = getPgivenRTh(data(i,j,lo,n));
                   Real plom1 = (-grav[IDIR] * dxarray[IDIR] / 2.0) * (data(i,j,lo,Rho_comp) + data(i,j,k,Rho_comp)) + plo;
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
      } else {
        // print some message like "shouldn't be here, this condition applies to scalars only   //
      }
    }
  };

  template<int IDIR>
  struct wall_scalars<IDIR, BCBound::upper> {
    using self = wall_scalars<IDIR, BCBound::upper>;

    static inline void apply(const amrex::Geometry& geom, const Box& b, Array4<Real> data, int ncomp=1) {
      int hi = geom.Domain().bigEnd(IDIR);

      amrex::IntVect type = b.type();
      if (type == amrex::IntVect::Zero) {

        /*    Real gravity = solverChoice.use_gravity? CONST_GRAV: 0.0;
        const    Array<Real,AMREX_SPACEDIM> grav{0.0, 0.0, gravity};
        const GpuArray<Real,AMREX_SPACEDIM> grav_gpu{grav[0], grav[1], grav[2]}; */

        if(b.bigEnd(IDIR) > hi) {
          if (IDIR == 0) {
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
                    Real phip1 = (grav[IDIR] * dxarray[IDIR] / 2.0) * (data(hi,j,k,Rho_comp) + data(i,j,k,Rho_comp)) + phi;
                    data(i,j,k,n) = getRhoThetagivenP(phip1); */
                    data(i,j,k,n) = data(hi,j,k,n);
                }
                else {
                    // For any passive scalars, also assume zero gradient //
                    data(i,j,k,n) = data(hi,j,k,n);
                }
              }
            });
          } else if (IDIR == 1) {
            amrex::ParallelFor(b, ncomp, [=] AMREX_GPU_DEVICE (int i, int j, int k, int n) noexcept
            {
              if (j > hi) {
                if (n == Rho_comp) {
                  data(i,j,k,n) = data(i,hi,k,n);
                }
                else if (n == RhoTheta_comp) {
                  /*  Real phi = getPgivenRTh(data(i,hi,k,n));
                      Real phip1 = (grav[IDIR] * dxarray[IDIR] / 2.0) * (data(i,hi,k,Rho_comp) + data(i,j,k,Rho_comp)) + phi;
                      data(i,j,k,n) = getRhoThetagivenP(phip1);  */
                  data(i,j,k,n) = data(i,hi,k,n);
                }
                else {
                  data(i,j,k,n) = data(i,hi,k,n);
                }
              }
            });
          } else if (IDIR == 2) {
            amrex::ParallelFor(b, ncomp, [=] AMREX_GPU_DEVICE (int i, int j, int k, int n) noexcept
            {
              if (k > hi) {
                if (n == Rho_comp) {
                  data(i,j,k,n) = data(i,j,hi,n);
                }
                else if (n == RhoTheta_comp) {
                  /* Real phi = getPgivenRTh(data(i,j,hi,n));
                  Real phip1 = (grav[IDIR] * dxarray[IDIR] / 2.0) * (data(i,j,hi,Rho_comp) + data(i,j,k,Rho_comp)) + phi;
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
      } else {
        // print some message like "shouldn't be here, this is for scalars only   //
      }
    }
  };

  /****************************************************************************************************************/

  // MOST wall
  template<int IDIR, BCBound Bound> struct ablmostwall;

    //
    // variables vartype = 0 (cell center), 1 (x velocity), 2 (y velocity), 3 (z velocity)
    //
  template<int IDIR>
  struct ablmostwall<IDIR, BCBound::lower> {
    using self = ablmostwall<IDIR, BCBound::lower>;

    static inline void apply_cell_based(const amrex::Geometry& geom, const Box& b,
                                              amrex::Vector<Array4<Real>*>& states, Array4<Real>& data,
                                              const ABLMost& most, int vartype, int ncomp=1)
    {
      AMREX_ALWAYS_ASSERT(IDIR == 2);
      int zlo = geom.Domain().smallEnd()[IDIR];

      if (b.smallEnd(IDIR) >= zlo) return;

      const Array4<Real>& velx_array   = *states[0];
      const Array4<Real>& vely_array   = *states[1];
      const Array4<Real>& velz_array   = *states[2];
      const Array4<Real>& states_array = *states[3];
      const Array4<Real>& eta_array    = *states[4];

      const int xhi = geom.Domain().bigEnd(0);
      const int yhi = geom.Domain().bigEnd(1);

      amrex::ParallelFor(b, ncomp, [=] AMREX_GPU_DEVICE (int i, int j, int k, int n) noexcept
      {
        if (k < zlo) {
          Real velx, vely, rho, theta, eta;

              if (vartype == 1) {
//              if (i < xhi) {
                  velx  = velx_array(i,j,zlo);
                  vely  = vely_array(i,j,zlo);
                  rho   = states_array(i,j,zlo,Rho_comp);
                  theta = states_array(i,j,zlo,RhoTheta_comp)/rho;
                  eta   = eta_array(i,j,zlo);
//              } else {
//                velx  = velx_array(i,j,zlo);
//                vely  = 2.0*vely_array(xhi-1,j,zlo)-(vely_array)(xhi-2,j,zlo);
//                rho   = 2.0*states_array(xhi-1,j,zlo,Rho_comp)-states_array(xhi-2,j,zlo,Rho_comp);
//                theta = (2.0*states_array(xhi-1,j,zlo,RhoTheta_comp)-states_array(xhi-2,j,zlo,RhoTheta_comp))/rho;
//                eta   = 2.0*eta_array(xhi-1,j,zlo)-eta_array(xhi-2,j,zlo);
//              }
              } else if (vartype == 2) {
//              if (j < yhi) {
                  velx  = velx_array(i,j,zlo);
                  vely  = vely_array(i,j,zlo);
                  rho   = states_array(i,j,zlo,Rho_comp);
                  theta = states_array(i,j,zlo,RhoTheta_comp)/rho;
                  eta   = eta_array(i,j,zlo);
//              } else {
//                velx  = 2.0*velx_array(i,yhi-1,zlo)-velx_array(i,yhi-2,zlo);
//                vely  = vely_array(i,j,zlo);
//                rho   = 2.0*states_array(i,yhi-1,zlo,Rho_comp)-states_array(i,yhi-2,zlo,Rho_comp);
//                theta = (2.0*states_array(i,yhi-1,zlo,RhoTheta_comp)-states_array(i,yhi-2,zlo,RhoTheta_comp))/rho;
//                eta   = 2.0*eta_array(i,yhi-1,zlo)-eta_array(i,yhi-2,zlo);
//              }
              } else {
                  velx  = velx_array(i,j,zlo);
                  vely  = vely_array(i,j,zlo);
                  rho   = states_array(i,j,zlo,Rho_comp);
                  theta = states_array(i,j,zlo,RhoTheta_comp)/rho;
                  eta   = eta_array(i,j,zlo);
              }

              Real vmag  = sqrt(velx*velx+vely*vely);
              Real vgx = ((velx-most.vel_mean[0])*most.vmag_mean + vmag*most.vel_mean[0])/(most.vmag_mean*most.vmag_mean)*most.utau*most.utau;
              Real vgy = ((vely-most.vel_mean[1])*most.vmag_mean + vmag*most.vel_mean[1])/(most.vmag_mean*most.vmag_mean)*most.utau*most.utau;

              Real num1  = (theta-most.theta_mean)*most.vmag_mean;
              Real num2  = (most.theta_mean-most.surf_temp)*most.vmag_mean;
              Real motheta = (num1+num2)*most.utau*most.kappa/most.phi_h();

              //printf("i=%d, j=%d, k=%d, vx=%13.6e, vy=%13.6e, rho=%13.6e\n",i,j,k,vgx,vgy,rho);

              if (vartype == 0 && n == RhoTheta_comp) {
                 data(i,j,k,n) = rho*motheta;
              } else if (vartype == 1) { // for velx, and vely
                 data(i,j,k,n) = vgx*rho; //eta;
              } else if (vartype == 2) {
                 data(i,j,k,n) = vgy*rho; //eta;
              } else {
                 data(i,j,k,n) = (zlo-k+1)*data(i,j,zlo,n) - (zlo-k)*data(i,j,zlo+1,n);
              }
            }
          });
    }

    static inline void
    apply_face_based(const amrex::Geometry& geom, const Box& b, Array4<Real>& data, int ncomp=1) {
      int lo = geom.Domain().smallEnd(IDIR);
      if(b.smallEnd(IDIR) < lo) {
        if (IDIR == 0) {
           amrex::ParallelFor(b, ncomp, [=] AMREX_GPU_DEVICE (int i, int j, int k, int n) noexcept
           {
            if (i < lo) {
              data(i,j,k,n) = 0.0;
            }
           });
        } else if (IDIR == 1) {
           amrex::ParallelFor(b, ncomp, [=] AMREX_GPU_DEVICE (int i, int j, int k, int n) noexcept
           {
            if (j < lo) {
              data(i,j,k,n) = 0.0;
            }
           });
        } else if (IDIR == 2) {
           amrex::ParallelFor(b, ncomp, [=] AMREX_GPU_DEVICE (int i, int j, int k, int n) noexcept
           {
            if (k < lo) {
              data(i,j,k,n) = 0.0;
            }
           });
        }
      }
    }

 };

  template<int IDIR>
  struct ablmostwall<IDIR, BCBound::upper> {
    using self = ablmostwall<IDIR, BCBound::upper>;

    static inline void
    apply_cell_based(const amrex::Geometry& /*geom*/, const Box& /*b*/,
                     amrex::Vector<Array4<Real>*>& /*states*/, Array4<Real>& /*data*/,
                     const ABLMost& /*most*/, int /*vartype*/, int ncomp=1) {};

    static inline void
    apply_face_based(const amrex::Geometry& /*geom*/, const Box& /*b*/, Array4<Real>& /*data*/, int ncomp=1) {}
 };

  /****************************************************************************************************************/

}

namespace phys_bcs {

class BCBase {
  public:
  virtual void applyBC (const amrex::Geometry /*geom*/, amrex::Vector<MultiFab*>& /*vars*/, const SolverChoice& /*solverChoice*/) = 0;
  virtual bool isInterior() { return false; }
  virtual bool isNoSlipWall() { return false; }
};

}

#endif
