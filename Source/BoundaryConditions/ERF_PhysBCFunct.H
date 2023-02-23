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

class ERFPhysBCFunct
{
public:
    ERFPhysBCFunct (const int lev,
                    const amrex::Geometry& geom, const amrex::Vector<amrex::BCRec>& domain_bcs_type,
                    const amrex::Gpu::DeviceVector<amrex::BCRec>& domain_bcs_type_d,
                    const int& terrain_type,
                    amrex::Array<amrex::Array<amrex::Real,AMREX_SPACEDIM*2>,AMREX_SPACEDIM+NVAR> bc_extdir_vals,
                    amrex::Array<amrex::Array<amrex::Real,AMREX_SPACEDIM*2>,AMREX_SPACEDIM+NVAR> bc_neumann_vals,
                    std::unique_ptr<amrex::MultiFab>& z_phys_nd,
                    std::unique_ptr<amrex::MultiFab>& detJ_cc)
        : m_lev(lev),
          m_geom(geom), m_domain_bcs_type(domain_bcs_type),
          m_domain_bcs_type_d(domain_bcs_type_d),
          m_terrain_type(terrain_type),
          m_bc_extdir_vals(bc_extdir_vals),
          m_bc_neumann_vals(bc_neumann_vals),
          m_z_phys_nd(z_phys_nd),
          m_detJ_cc(detJ_cc)
          {}

    ~ERFPhysBCFunct () {}

    //
    // mfs is the vector of MultiFabs to be filled
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
    void operator() (const amrex::Vector<amrex::MultiFab*>& mfs, int icomp, int ncomp,
                     amrex::IntVect const& nghost_cons, amrex::IntVect const& nghost_vels,
                     amrex::Real time, std::string& init_type, bool cons_only);

    void impose_lateral_xvel_bcs (const amrex::Array4<amrex::Real>& dest_arr,
                                  const amrex::Box& bx, const amrex::Box& domain,
                                  amrex::Real time, int bccomp);
    void impose_vertical_xvel_bcs (const amrex::Array4<amrex::Real>& dest_arr,
                                  const amrex::Box& bx, const amrex::Box& domain,
                                   const amrex::Array4<amrex::Real const>& z_nd,
                                   const amrex::GpuArray<amrex::Real,AMREX_SPACEDIM> dxInv,
                                   amrex::Real time, int bccomp);

    void impose_lateral_yvel_bcs (const amrex::Array4<amrex::Real>& dest_arr,
                                  const amrex::Box& bx, const amrex::Box& domain,
                                  amrex::Real time, int bccomp);
    void impose_vertical_yvel_bcs (const amrex::Array4<amrex::Real>& dest_arr,
                                   const amrex::Box& bx, const amrex::Box& domain,
                                   const amrex::Array4<amrex::Real const>& z_nd,
                                   const amrex::GpuArray<amrex::Real,AMREX_SPACEDIM> dxInv,
                                   amrex::Real time, int bccomp);


    void impose_lateral_zvel_bcs (const amrex::Array4<amrex::Real>& dest_arr,
                                  const amrex::Box& bx, const amrex::Box& domain,
                                  const amrex::Array4<amrex::Real const>& z_nd,
                                  const amrex::GpuArray<amrex::Real,AMREX_SPACEDIM> dxInv,
                                  amrex::Real time, int bccomp_w);
    void impose_vertical_zvel_bcs (const amrex::Array4<amrex::Real>& dest_arr,
                                   const amrex::Box& bx, const amrex::Box& domain,
                                   const amrex::Array4<amrex::Real const>& z_nd,
                                   const amrex::GpuArray<amrex::Real,AMREX_SPACEDIM> dxInv,
                                   amrex::Real time, int bccomp_u, int bccomp_v, int bccomp_w,
                                   int terrain_type);

    void impose_lateral_cons_bcs (const amrex::Array4<amrex::Real>& mf,
                                  const amrex::Box& bx, const amrex::Box& domain,
                                  int icomp, int ncomp, amrex::Real time, int bccomp);
    void impose_vertical_cons_bcs (const amrex::Array4<amrex::Real>& mf,
                                   const amrex::Box& bx, const amrex::Box& domain,
                                   const amrex::Array4<amrex::Real const>& z_nd,
                                   const amrex::GpuArray<amrex::Real,AMREX_SPACEDIM> dxInv,
                                   int icomp, int ncomp, amrex::Real time, int bccomp);

private:
    int                  m_lev;
    amrex::Geometry      m_geom;
    amrex::Vector<amrex::BCRec>            m_domain_bcs_type;
    amrex::Gpu::DeviceVector<amrex::BCRec> m_domain_bcs_type_d;
    int           m_terrain_type;
    amrex::Array<amrex::Array<amrex::Real, AMREX_SPACEDIM*2>,AMREX_SPACEDIM+NVAR> m_bc_extdir_vals;
    amrex::Array<amrex::Array<amrex::Real, AMREX_SPACEDIM*2>,AMREX_SPACEDIM+NVAR> m_bc_neumann_vals;
    std::unique_ptr<amrex::MultiFab>& m_z_phys_nd;
    std::unique_ptr<amrex::MultiFab>&   m_detJ_cc;
};

#endif
