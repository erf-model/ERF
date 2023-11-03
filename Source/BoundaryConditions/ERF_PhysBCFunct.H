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
                    const TerrainType& terrain_type,
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
          m_z_phys_nd(z_phys_nd.get()),
          m_detJ_cc(detJ_cc.get())
    {}

    ~ERFPhysBCFunct () {}

   /*
    * Impose physical boundary conditions at domain boundaries
    *
    * @param[out] mfs        Vector of MultiFabs to be filled containing, in order: cons, xvel, yvel, and zvel data
    * @param[in] icomp_cons  starting component for conserved variables
    * @param[in] ncomp_cons  number of components for conserved variables
    * @param[in] nghost_cons number of ghost cells to be filled for conserved variables
    * @param[in] nghost_vels number of ghost cells to be filled for velocity components
    * @param[in] time        time at which the data should be filled
    * @param[in] init_type   if "real" then we fill boundary conditions for interior locations
    * @param[in] cons_only   if 1 then only fill conserved variables
    */
    void operator() (const amrex::Vector<amrex::MultiFab*>& mfs, int icomp, int ncomp,
                     amrex::IntVect const& nghost_cons, amrex::IntVect const& nghost_vels,
                     std::string& init_type, bool cons_only, int bccomp_cons, const amrex::Real time = 0.0);

    void impose_lateral_xvel_bcs (const amrex::Array4<amrex::Real>& dest_arr,
                                  const amrex::Box& bx, const amrex::Box& domain,
                                  int bccomp);
    void impose_vertical_xvel_bcs (const amrex::Array4<amrex::Real>& dest_arr,
                                   const amrex::Box& bx, const amrex::Box& domain,
                                   const amrex::Array4<amrex::Real const>& z_nd,
                                   const amrex::GpuArray<amrex::Real,AMREX_SPACEDIM> dxInv,
                                   int bccomp,
                                   const amrex::Real time);

    void impose_lateral_yvel_bcs (const amrex::Array4<amrex::Real>& dest_arr,
                                  const amrex::Box& bx, const amrex::Box& domain,
                                  int bccomp);
    void impose_vertical_yvel_bcs (const amrex::Array4<amrex::Real>& dest_arr,
                                   const amrex::Box& bx, const amrex::Box& domain,
                                   const amrex::Array4<amrex::Real const>& z_nd,
                                   const amrex::GpuArray<amrex::Real,AMREX_SPACEDIM> dxInv,
                                   int bccomp);


    void impose_lateral_zvel_bcs (const amrex::Array4<amrex::Real>& dest_arr,
                                  const amrex::Array4<amrex::Real const>& xvel_arr,
                                  const amrex::Array4<amrex::Real const>& yvel_arr,
                                  const amrex::Box& bx, const amrex::Box& domain,
                                  const amrex::Array4<amrex::Real const>& z_nd,
                                  const amrex::GpuArray<amrex::Real,AMREX_SPACEDIM> dxInv,
                                  int bccomp_w);
    void impose_vertical_zvel_bcs (const amrex::Array4<amrex::Real>& dest_arr,
                                   const amrex::Array4<amrex::Real const>& xvel_arr,
                                   const amrex::Array4<amrex::Real const>& yvel_arr,
                                   const amrex::Box& bx, const amrex::Box& domain,
                                   const amrex::Array4<amrex::Real const>& z_nd,
                                   const amrex::GpuArray<amrex::Real,AMREX_SPACEDIM> dxInv,
                                   int bccomp_u, int bccomp_v, int bccomp_w,
                                   TerrainType terrain_type);

    void impose_lateral_cons_bcs (const amrex::Array4<amrex::Real>& dest_arr,
                                  const amrex::Box& bx, const amrex::Box& domain,
                                  int icomp, int ncomp, int bccomp);
    void impose_vertical_cons_bcs (const amrex::Array4<amrex::Real>& dest_arr,
                                   const amrex::Box& bx, const amrex::Box& domain,
                                   const amrex::Array4<amrex::Real const>& z_nd,
                                   const amrex::GpuArray<amrex::Real,AMREX_SPACEDIM> dxInv,
                                   int icomp, int ncomp, int bccomp);

private:
    int                  m_lev;
    amrex::Geometry      m_geom;
    amrex::Vector<amrex::BCRec>            m_domain_bcs_type;
    amrex::Gpu::DeviceVector<amrex::BCRec> m_domain_bcs_type_d;
    TerrainType          m_terrain_type;
    amrex::Array<amrex::Array<amrex::Real, AMREX_SPACEDIM*2>,AMREX_SPACEDIM+NVAR> m_bc_extdir_vals;
    amrex::Array<amrex::Array<amrex::Real, AMREX_SPACEDIM*2>,AMREX_SPACEDIM+NVAR> m_bc_neumann_vals;
    amrex::MultiFab* m_z_phys_nd;
    amrex::MultiFab* m_detJ_cc;
};

#endif
