#ifndef NOAH_H
#define NOAH_H

#include <string>
#include <vector>
#include <memory>

#include <AMReX_FArrayBox.H>
#include <AMReX_Geometry.H>
#include <AMReX_MultiFabUtil.H>

#include <ERF_NullSurf.H>
#include <ERF_Constants.H>
#include <ERF_IndexDefines.H>
#include <ERF_DataStruct.H>

// External include from the noahmp library
#include <NoahmpIO.H>

class NOAH : public NullSurf {
public:
    // Constructor
    NOAH () {}

    // Destructor
    virtual ~NOAH () = default;

    // Initialize data structures
    void
    Init (const amrex::MultiFab& cons_in,
          const amrex::Geometry& geom,
          const amrex::Real& dt) override;


private:

   // C++ variable for NoahmpIO struct
   NoahmpIO_struct noahmpio;

};
#endif
