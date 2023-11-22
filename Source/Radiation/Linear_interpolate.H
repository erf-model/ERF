//
// interpolation of data in latitude, longitude, and time
//
#ifndef ERF_LINEAR_INTERP_H_
#define ERF_LINEAR_INTERP_H_

#include "rrtmgp_const.h"
using yakl::fortran::parallel_for;
using yakl::fortran::SimpleBounds;

class LinInterp
{
 public:
   enum InterpMethod{
     extrap_method_zero = 0,
     extrap_method_bndry = 1,
     extrap_method_cycle = 2,
   };

   struct InterpType {
     real1d wgts;
     real1d wgtn;
     int1d jjm;
     int1d jjp;
   };

  // Initialize a variable of type(interp_type) with weights for linear interpolation.
  //       this variable can then be used in calls to lininterp1d and lininterp2d.
  //   yin is a 1d array of length nin of locations to interpolate from - this array must
  //       be monotonic but can be increasing or decreasing
  //   yout is a 1d array of length nout of locations to interpolate to, this array need
  //       not be ordered
  //   extrap_method determines how to handle yout points beyond the bounds of yin
  //       if 0 set values outside output grid to 0
  //       if 1 set to boundary value
  //       if 2 set to cyclic boundaries
  //         optional values cyclicmin and cyclicmax can be used to set the bounds of the
  //         cyclic mapping - these default to 0 and 360.
  YAKL_INLINE
  static void init(const real1d& yin,
                   const int& nin,
                   const real1d& yout,
                   const int& nout,
                   const InterpMethod& extrap_method,
                   InterpType& interp_wgts,
                   real cyclicmin = 0.,
                   real cyclicmax= 0.) {
      real cmin, cmax;
      real dyinwrap;
      real avgdyin;

      //
      // Check validity of input coordinate arrays: must be monotonically increasing,
      // and have a total of at least 2 elements
      //
      if(cyclicmin != 0.) {
         cmin=cyclicmin;
      } else {
         cmin=0;
      }
      if(cyclicmax != 0.) {
         cmax=cyclicmax;
      } else {
         cmax=360.0;
      }

      bool increasing  = true;

      interp_wgts.jjm  = int1d("jjm", nout);
      interp_wgts.jjp  = int1d("jjp", nout);
      interp_wgts.wgts = real1d("wgts", nout);
      interp_wgts.wgtn = real1d("wgtn", nout);

      parallel_for(SimpleBounds<1>(nin), YAKL_LAMBDA (int j) {
         if(yin(j) > yin(j+1)) {
           printf("inputs are not monotonic!\n");
           return;
         }
      });

      //
      // Initialize index arrays for later checking
      //
      parallel_for(SimpleBounds<1>(nout), YAKL_LAMBDA (int j) {
        interp_wgts.jjm(j) = 0;
        interp_wgts.jjp(j) = 0;
      });

      switch (extrap_method) {
        case extrap_method_zero:
         //
         // For values which extend beyond boundaries, set weights
         // such that values will be 0.
         //
         if(increasing) {
            parallel_for(SimpleBounds<1>(nout), YAKL_LAMBDA (int j) {
               if (yout(j) < yin(1)) {
                  interp_wgts.jjm(j) = 1;
                  interp_wgts.jjp(j) = 1;
                  interp_wgts.wgts(j) = 0.;
                  interp_wgts.wgtn(j) = 0.;
               }
               else if (yout(j) > yin(nin)) {
                  interp_wgts.jjm(j) = nin;
                  interp_wgts.jjp(j) = nin;
                  interp_wgts.wgts(j) = 0.;
                  interp_wgts.wgtn(j) = 0.;
               }
            });
         }
         else {
            parallel_for(SimpleBounds<1>(nout), YAKL_LAMBDA (int j) {
               if (yout(j) > yin(1)) {
                  interp_wgts.jjm(j) = 1;
                  interp_wgts.jjp(j) = 1;
                  interp_wgts.wgts(j) = 0.;
                  interp_wgts.wgtn(j) = 0.;
               }
                else if (yout(j) < yin(nin)) {
                  interp_wgts.jjm(j) = nin;
                  interp_wgts.jjp(j) = nin;
                  interp_wgts.wgts(j) = 0.;
                  interp_wgts.wgtn(j) = 0.;
               }
           });
         }
         break;
        case extrap_method_bndry:
         //
         // For values which extend beyond boundaries, set weights
         // such that values will just be copied.
         //
         if(increasing) {
            parallel_for(SimpleBounds<1>(nout), YAKL_LAMBDA (int j) {
               if (yout(j) <= yin(1)) {
                  interp_wgts.jjm(j) = 1;
                  interp_wgts.jjp(j) = 1;
                  interp_wgts.wgts(j) = 1.;
                  interp_wgts.wgtn(j) = 0.;
               }
               else if (yout(j) > yin(nin)) {
                  interp_wgts.jjm(j) = nin;
                  interp_wgts.jjp(j) = nin;
                  interp_wgts.wgts(j) = 1.;
                  interp_wgts.wgtn(j) = 0.;
              }
          });
        }
        else {
           parallel_for(SimpleBounds<1>(nout), YAKL_LAMBDA (int j) {
               if (yout(j) > yin(1)) {
                  interp_wgts.jjm(j) = 1;
                  interp_wgts.jjp(j) = 1;
                  interp_wgts.wgts(j) = 1.;
                  interp_wgts.wgtn(j) = 0.;
               }
               else if (yout(j) <= yin(nin)) {
                  interp_wgts.jjm(j) = nin;
                  interp_wgts.jjp(j) = nin;
                  interp_wgts.wgts(j) = 1.;
                  interp_wgts.wgtn(j) = 0.;
              }
           });
        }
        break;
      case extrap_method_cycle:
         //
         // For values which extend beyond boundaries, set weights
         // for circular boundaries
         //
         dyinwrap = yin(1) + (cmax-cmin) - yin(nin);
         avgdyin = abs(yin(nin)-yin(1))/(nin-1.);
         auto ratio = dyinwrap/avgdyin;
         if (ratio < 0.9 || ratio > 1.1) {
            printf("ratio is too large, dyinwrap= %13.6e, avgdyin= %13.6e, yin(1) = %13.6e, yin(nin)= %13.6e\n", 
                  dyinwrap, avgdyin, yin(1), yin(nin));
         }

         if(increasing) {
            parallel_for(SimpleBounds<1>(nout), YAKL_LAMBDA (int j) {
               if (yout(j) <= yin(1)) {
                  interp_wgts.jjm(j) = nin;
                  interp_wgts.jjp(j) = 1;
                  interp_wgts.wgts(j) = (yin(1)-yout(j))/dyinwrap;
                  interp_wgts.wgtn(j) = (yout(j)+(cmax-cmin) - yin(nin))/dyinwrap;
               }
               else if (yout(j) > yin(nin)) {
                  interp_wgts.jjm(j) = nin;
                  interp_wgts.jjp(j) = 1;
                  interp_wgts.wgts(j) = (yin(1)+(cmax-cmin)-yout(j))/dyinwrap;
                  interp_wgts.wgtn(j) = (yout(j)-yin(nin))/dyinwrap;
               }
            });
          }
          else {
             parallel_for(SimpleBounds<1>(nout), YAKL_LAMBDA (int j) {
               if (yout(j) > yin(1)) {
                  interp_wgts.jjm(j) = nin;
                  interp_wgts.jjp(j) = 1;
                  interp_wgts.wgts(j) = (yin(1)-yout(j))/dyinwrap;
                  interp_wgts.wgtn(j) = (yout(j)+(cmax-cmin) - yin(nin))/dyinwrap;
               }
               else if (yout(j) <= yin(nin)) {
                  interp_wgts.jjm(j) = nin;
                  interp_wgts.jjp(j) = 1;
                  interp_wgts.wgts(j) = (yin(1)+(cmax-cmin)-yout(j))/dyinwrap;
                  interp_wgts.wgtn(j) = (yout(j)+(cmax-cmin)-yin(nin))/dyinwrap;
               }
             });
         }
         break;
      }
      //
      // Loop though output indices finding input indices and weights
      //
      if(increasing) {
         parallel_for(SimpleBounds<2>(nout, nin-1), YAKL_LAMBDA (int j, int jj) {
            if (yout(j) > yin(jj) && yout(j) <= yin(jj+1)) {
                  interp_wgts.jjm(j) = jj;
                  interp_wgts.jjp(j) = jj + 1;
                  interp_wgts.wgts(j) = (yin(jj+1)-yout(j))/(yin(jj+1)-yin(jj));
                  interp_wgts.wgtn(j) = (yout(j)-yin(jj))/(yin(jj+1)-yin(jj));
                  return;
            }
         });
      }
      else {
         parallel_for(SimpleBounds<2>(nout, nin-1), YAKL_LAMBDA (int j, int jj) {
             if (yout(j) <= yin(jj) && yout(j) > yin(jj+1)) {
                interp_wgts.jjm(j) = jj;
                interp_wgts.jjp(j) = jj + 1;
                interp_wgts.wgts(j) = (yin(jj+1)-yout(j))/(yin(jj+1)-yin(jj));
                interp_wgts.wgtn(j) = (yout(j)-yin(jj))/(yin(jj+1)-yin(jj));
                return;
             }
      });
    }

    //
    // Check that interp/extrap points have been found for all outputs
    //
    parallel_for(SimpleBounds<1>(nout), YAKL_LAMBDA (int j) {
       int count = 0;
       if (interp_wgts.jjm(j) == 0 || interp_wgts.jjp(j) == 0) count +=  1;
       auto frac = interp_wgts.wgts(j) + interp_wgts.wgtn(j);
       if ((frac < 0.9 || frac > 1.1) && extrap_method != 0)
         printf("interpolation error!\n");
    });

   }

   //
   // Purpose: Do a linear interpolation from input mesh to output
   //          mesh with weights as set in lininterp_init.
   YAKL_INLINE
   static void interp1d(const real1d& arrin,
                        const int& n1,
                        const real1d& arrout,
                        const int& m1,
                        const InterpType& interp_wgts) {
     //
     // Do the interpolation
     //
     parallel_for(SimpleBounds<1>(m1), YAKL_LAMBDA (int j) {
        arrout(j) = arrin(interp_wgts.jjm(j))*interp_wgts.wgts(j) +
                    arrin(interp_wgts.jjp(j))*interp_wgts.wgtn(j);
     });
   }

   YAKL_INLINE
   static void interp2d2d(const real2d& arrin,
                          const int& n1,
                          const int& n2,
                          const real2d& arrout,
                          const int& m1,
                          const int& m2,
                          const InterpType& wgt1,
                          const InterpType& wgt2) {

     real2d arrtmp("arrtmp",n1,m2);

     parallel_for(SimpleBounds<2>(n1, m2), YAKL_LAMBDA (int i, int j) {
        arrtmp(i,j) = arrin(i,wgt2.jjm(j))*wgt2.wgts(j) + arrin(i,wgt2.jjp(j))*wgt2.wgtn(j);
     });

     parallel_for(SimpleBounds<2>(n1, m2), YAKL_LAMBDA (int i, int j) {
        arrout(i,j) = arrtmp(wgt1.jjm(i),j)*wgt1.wgts(i) + arrtmp(wgt1.jjp(i),j)*wgt1.wgtn(i);
     });
   }

   YAKL_INLINE
   static void interp2d1d(const real2d& arrin,
                          const int& n1,
                          const int& n2,
                          const real1d& arrout,
                          const int& m1,
                          const InterpType& wgt1,
                          const InterpType& wgt2) {

     parallel_for(SimpleBounds<1>(m1), YAKL_LAMBDA (int i) {
        arrout(i) = arrin(wgt1.jjm(i),wgt2.jjm(i))*wgt1.wgts(i)*wgt2.wgts(i)+arrin(wgt1.jjp(i),wgt2.jjm(i))*wgt1.wgtn(i)*wgt2.wgts(i)
                  + arrin(wgt1.jjm(i),wgt2.jjp(i))*wgt1.wgtn(i)*wgt2.wgtn(i)+arrin(wgt1.jjp(i),wgt2.jjp(i))*wgt1.wgtn(i)*wgt2.wgtn(i);
     });
   }

   YAKL_INLINE
   static void interp3d2d(const real3d& arrin,
                          const int& n1,
                          const int& n2,
                          const int& n3,
                          const real2d& arrout,
                          const int& m1,
                          const int& len1,
                          const InterpType& wgt1,
                          const InterpType& wgt2) {

     parallel_for(SimpleBounds<2>(m1, n3), YAKL_LAMBDA (int i, int k) {
        arrout(i,k) = arrin(wgt1.jjm(i),wgt2.jjm(i),k)*wgt1.wgts(i)*wgt2.wgts(i)+arrin(wgt1.jjp(i),wgt2.jjm(i),k)*wgt1.wgtn(i)*wgt2.wgts(i) +
                      arrin(wgt1.jjm(i),wgt2.jjp(i),k)*wgt1.wgts(i)*wgt2.wgtn(i)+arrin(wgt1.jjp(i),wgt2.jjp(i),k)*wgt1.wgtn(i)*wgt2.wgtn(i);
    });
  }
};
#endif
