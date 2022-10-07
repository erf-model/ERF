#ifndef ERF_PhysBCFunct_H_
#define ERF_PhysBCFunct_H_
#include <AMReX_Config.H>

#include <AMReX_BCRec.H>
#include <AMReX_Geometry.H>
#include <AMReX_MultiFab.H>
#include <AMReX_ArrayLim.H>
#include <AMReX_FilCC_C.H>
#include <AMReX_FilND_C.H>
#include <AMReX_FilFC_C.H>
#include "AMReX_TypeTraits.H"
#include "AMReX_Orientation.H"

#include <ERF_ReadBndryPlanes.H>
#include <TimeInterpolatedData.H>
#include <IndexDefines.H>
#include <DataStruct.H>
#include <EddyViscosity.H>
#include <TerrainMetrics.H>

using PlaneVector = amrex::Vector<amrex::FArrayBox>;

class ERFPhysBCFunct
{
public:
    ERFPhysBCFunct (const int lev,
                    const amrex::Geometry& geom, const amrex::Vector<amrex::BCRec>& domain_bcs_type,
                    const amrex::Gpu::DeviceVector<amrex::BCRec>& domain_bcs_type_d,
                    const int& var_idx, const int& terrain_type, TimeInterpolatedData& data,
                    amrex::Array<amrex::Array<amrex::Real,AMREX_SPACEDIM*2>,AMREX_SPACEDIM+NVAR> bc_extdir_vals,
                    std::unique_ptr<amrex::MultiFab>& z_phys_nd,
                    std::unique_ptr<amrex::MultiFab>& detJ_cc)
        : m_lev(lev),
          m_geom(geom), m_domain_bcs_type(domain_bcs_type),
          m_domain_bcs_type_d(domain_bcs_type_d),
          m_var_idx(var_idx), m_terrain_type(terrain_type), m_data(data),
          m_bc_extdir_vals(bc_extdir_vals),
          m_z_phys_nd(z_phys_nd),
          m_detJ_cc(detJ_cc)
          {}

    ~ERFPhysBCFunct () {}

    //
    // mf is the multifab to be filled
    // icomp is the index into the MultiFab -- if cell-centered this can be any value
    //       from 0 to NVAR-1, if face-centered this must be 0
    // ncomp is the number of components -- if cell-centered (var_idx = 0) this can be any value
    //       from 1 to NVAR as long as icomp+ncomp <= NVAR-1.  If face-centered this
    //       must be 1
    // nghost is how many ghost cells to be filled
    // time is the time at which the data should be filled
    // bccomp is the index into both domain_bcs_type_bcr and bc_extdir_vals for icomp = 0  --
    //     so this follows the BCVars enum
    //
    void operator() (amrex::MultiFab& mf, int icomp, int ncomp, amrex::IntVect const& nghost,
                     amrex::Real time, int bccomp);

    void impose_xvel_bcs (const amrex::Array4<amrex::Real>& dest_arr, const amrex::Box& bx, const amrex::Box& domain,
                          const amrex::Array4<amrex::Real const>& z_nd,
                          const amrex::GpuArray<amrex::Real,AMREX_SPACEDIM> dxInv,
                          amrex::Real time, int bccomp);

    void impose_yvel_bcs (const amrex::Array4<amrex::Real>& dest_arr, const amrex::Box& bx, const amrex::Box& domain,
                          const amrex::Array4<amrex::Real const>& z_nd,
                          const amrex::GpuArray<amrex::Real,AMREX_SPACEDIM> dxInv,
                          amrex::Real time, int bccomp);


    void impose_zvel_bcs (const amrex::Array4<amrex::Real>& dest_arr,
                          const amrex::Box& bx, const amrex::Box& domain,
                          const amrex::Array4<amrex::Real const>& velx_arr,
                          const amrex::Array4<amrex::Real const>& vely_arr,
                          const amrex::Array4<amrex::Real const>& z_nd,
                          const amrex::GpuArray<amrex::Real,AMREX_SPACEDIM> dx,
                          const amrex::GpuArray<amrex::Real,AMREX_SPACEDIM> dxInv,
                          amrex::Real time, int bccomp, int terrain_type);

    void impose_cons_bcs (const amrex::Array4<amrex::Real>& mf, const amrex::Box& bx, const amrex::Box& domain,
                          const amrex::Array4<amrex::Real const>& z_nd,
                          const amrex::GpuArray<amrex::Real,AMREX_SPACEDIM> dxInv,
                          int icomp, int ncomp, amrex::Real time, int bccomp);

    // For backward compatibility
    // void FillBoundary (amrex::MultiFab& mf, int dcomp, int ncomp, amrex::IntVect const& nghost,
    //                    amrex::Real time, int bccomp) {
    //     this->operator()(mf,dcomp,ncomp,nghost,time,bccomp);
    // }

private:
    int                  m_lev;
    amrex::Geometry      m_geom;
    amrex::Vector<amrex::BCRec>            m_domain_bcs_type;
    amrex::Gpu::DeviceVector<amrex::BCRec> m_domain_bcs_type_d;
    int           m_var_idx;
    int           m_terrain_type;
    TimeInterpolatedData& m_data;
    amrex::Array<amrex::Array<amrex::Real, AMREX_SPACEDIM*2>,AMREX_SPACEDIM+NVAR> m_bc_extdir_vals;
    std::unique_ptr<amrex::MultiFab>& m_z_phys_nd;
    std::unique_ptr<amrex::MultiFab>&   m_detJ_cc;
};

#endif
