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

class ERFPhysBCFunct_cons
{
public:
    ERFPhysBCFunct_cons (const int lev,
                         const amrex::Geometry& geom, const amrex::Vector<amrex::BCRec>& domain_bcs_type,
                         const amrex::Gpu::DeviceVector<amrex::BCRec>& domain_bcs_type_d,
                         amrex::Array<amrex::Array<amrex::Real,AMREX_SPACEDIM*2>,AMREX_SPACEDIM+NBCVAR_max> bc_extdir_vals,
                         amrex::Array<amrex::Array<amrex::Real,AMREX_SPACEDIM*2>,AMREX_SPACEDIM+NBCVAR_max> bc_neumann_vals,
                         std::unique_ptr<amrex::MultiFab>& z_phys_nd,
                         const bool use_real_bcs)
        : m_lev(lev), m_geom(geom),
          m_domain_bcs_type(domain_bcs_type),
          m_domain_bcs_type_d(domain_bcs_type_d),
          m_bc_extdir_vals(bc_extdir_vals),
          m_bc_neumann_vals(bc_neumann_vals),
          m_z_phys_nd(z_phys_nd.get()),
          m_use_real_bcs(use_real_bcs)
    {}

    ~ERFPhysBCFunct_cons () {}

   /*
    * Impose physical boundary conditions at domain boundaries
    *
    * @param[out] mfs         Vector of MultiFabs to be filled containing, in order: cons, xvel, yvel, and zvel data
    * @param[in] icomp        starting component for conserved variables
    * @param[in] ncomp        number of components for conserved variables
    * @param[in] nghost       number of ghost cells to be filled for conserved variables
    * @param[in] time         time at which the data should be filled
    * @param[in] use_real_bcs if true then we fill boundary conditions for interior locations
    */
    void operator() (amrex::MultiFab& mf, int icomp, int ncomp,
                     amrex::IntVect const& nghost, const amrex::Real time, int bccomp_cons);

    void impose_lateral_cons_bcs (const amrex::Array4<amrex::Real>& dest_arr,
                                  const amrex::Box& bx, const amrex::Box& domain,
                                  int icomp, int ncomp, int ngz);
    void impose_vertical_cons_bcs (const amrex::Array4<amrex::Real>& dest_arr,
                                   const amrex::Box& bx, const amrex::Box& domain,
                                   const amrex::Array4<amrex::Real const>& z_nd,
                                   const amrex::GpuArray<amrex::Real,AMREX_SPACEDIM> dxInv,
                                   int icomp, int ncomp);

private:
    int                  m_lev;
    amrex::Geometry      m_geom;
    amrex::Vector<amrex::BCRec>            m_domain_bcs_type;
    amrex::Gpu::DeviceVector<amrex::BCRec> m_domain_bcs_type_d;
    amrex::Array<amrex::Array<amrex::Real, AMREX_SPACEDIM*2>,AMREX_SPACEDIM+NBCVAR_max> m_bc_extdir_vals;
    amrex::Array<amrex::Array<amrex::Real, AMREX_SPACEDIM*2>,AMREX_SPACEDIM+NBCVAR_max> m_bc_neumann_vals;
    amrex::MultiFab* m_z_phys_nd;
    bool                 m_use_real_bcs;
};

class ERFPhysBCFunct_u
{
public:
    ERFPhysBCFunct_u (const int lev, const amrex::Geometry& geom,
                      const amrex::Vector<amrex::BCRec>& domain_bcs_type,
                      const amrex::Gpu::DeviceVector<amrex::BCRec>& domain_bcs_type_d,
                      amrex::Array<amrex::Array<amrex::Real,AMREX_SPACEDIM*2>,AMREX_SPACEDIM+NBCVAR_max> bc_extdir_vals,
                      amrex::Array<amrex::Array<amrex::Real,AMREX_SPACEDIM*2>,AMREX_SPACEDIM+NBCVAR_max> bc_neumann_vals,
                      std::unique_ptr<amrex::MultiFab>& z_phys_nd,
                      const bool use_real_bcs,
                      amrex::Real* u_bc_data)
        : m_lev(lev), m_geom(geom),
          m_domain_bcs_type(domain_bcs_type),
          m_domain_bcs_type_d(domain_bcs_type_d),
          m_bc_extdir_vals(bc_extdir_vals),
          m_bc_neumann_vals(bc_neumann_vals),
          m_z_phys_nd(z_phys_nd.get()),
          m_use_real_bcs(use_real_bcs),
          m_u_bc_data(u_bc_data)
    { }

    ~ERFPhysBCFunct_u () {}

   /*
    * Impose physical boundary conditions at domain boundaries
    *
    * @param[out] mfs         Vector of MultiFabs to be filled containing, in order: cons, xvel, yvel, and zvel data
    * @param[in] nghost       number of ghost cells to be filled for velocity components
    * @param[in] time         time at which the data should be filled
    * @param[in] use_real_bcs if true then we fill boundary conditions for interior locations
    */
    void operator() (amrex::MultiFab& mf, int icomp, int ncomp,
                     amrex::IntVect const& nghost, const amrex::Real time, int bccomp);

    void impose_lateral_xvel_bcs (const amrex::Array4<amrex::Real>& dest_arr,
                                  const amrex::Box& bx, const amrex::Box& domain,
                                  int bccomp);
    void impose_vertical_xvel_bcs (const amrex::Array4<amrex::Real>& dest_arr,
                                   const amrex::Box& bx, const amrex::Box& domain,
                                   const amrex::Array4<amrex::Real const>& z_nd,
                                   const amrex::GpuArray<amrex::Real,AMREX_SPACEDIM> dxInv,
                                   int bccomp,
                                   const amrex::Real time);

private:
    int                  m_lev;
    amrex::Geometry      m_geom;
    amrex::Vector<amrex::BCRec>            m_domain_bcs_type;
    amrex::Gpu::DeviceVector<amrex::BCRec> m_domain_bcs_type_d;
    amrex::Array<amrex::Array<amrex::Real, AMREX_SPACEDIM*2>,AMREX_SPACEDIM+NBCVAR_max> m_bc_extdir_vals;
    amrex::Array<amrex::Array<amrex::Real, AMREX_SPACEDIM*2>,AMREX_SPACEDIM+NBCVAR_max> m_bc_neumann_vals;
    amrex::MultiFab* m_z_phys_nd;
    bool                 m_use_real_bcs;
    amrex::Real* m_u_bc_data;
};

class ERFPhysBCFunct_v
{
public:
    ERFPhysBCFunct_v (const int lev, const amrex::Geometry& geom,
                      const amrex::Vector<amrex::BCRec>& domain_bcs_type,
                      const amrex::Gpu::DeviceVector<amrex::BCRec>& domain_bcs_type_d,
                      amrex::Array<amrex::Array<amrex::Real,AMREX_SPACEDIM*2>,AMREX_SPACEDIM+NBCVAR_max> bc_extdir_vals,
                      amrex::Array<amrex::Array<amrex::Real,AMREX_SPACEDIM*2>,AMREX_SPACEDIM+NBCVAR_max> bc_neumann_vals,
                      std::unique_ptr<amrex::MultiFab>& z_phys_nd,
                      const bool use_real_bcs,
                      amrex::Real* v_bc_data)
        : m_lev(lev),
          m_geom(geom), m_domain_bcs_type(domain_bcs_type),
          m_domain_bcs_type_d(domain_bcs_type_d),
          m_bc_extdir_vals(bc_extdir_vals),
          m_bc_neumann_vals(bc_neumann_vals),
          m_z_phys_nd(z_phys_nd.get()),
          m_use_real_bcs(use_real_bcs),
          m_v_bc_data(v_bc_data)
    { }

    ~ERFPhysBCFunct_v () {}

   /*
    * Impose physical boundary conditions at domain boundaries
    *
    * @param[out] mfs         Vector of MultiFabs to be filled containing, in order: cons, xvel, yvel, and zvel data
    * @param[in] icomp_cons   starting component for conserved variables
    * @param[in] ncomp_cons   number of components for conserved variables
    * @param[in] nghost       number of ghost cells to be filled for velocity components
    * @param[in] time         time at which the data should be filled
    * @param[in] use_real_bcs if true then we fill boundary conditions for interior locations
    */
    void operator() (amrex::MultiFab& mf, int icomp, int ncomp,
                     amrex::IntVect const& nghost, const amrex::Real time, int bccomp);

    void impose_lateral_yvel_bcs (const amrex::Array4<amrex::Real>& dest_arr,
                                  const amrex::Box& bx, const amrex::Box& domain,
                                  int bccomp);
    void impose_vertical_yvel_bcs (const amrex::Array4<amrex::Real>& dest_arr,
                                   const amrex::Box& bx, const amrex::Box& domain,
                                   const amrex::Array4<amrex::Real const>& z_nd,
                                   const amrex::GpuArray<amrex::Real,AMREX_SPACEDIM> dxInv,
                                   int bccomp);

private:
    int                  m_lev;
    amrex::Geometry      m_geom;
    amrex::Vector<amrex::BCRec>            m_domain_bcs_type;
    amrex::Gpu::DeviceVector<amrex::BCRec> m_domain_bcs_type_d;
    amrex::Array<amrex::Array<amrex::Real, AMREX_SPACEDIM*2>,AMREX_SPACEDIM+NBCVAR_max> m_bc_extdir_vals;
    amrex::Array<amrex::Array<amrex::Real, AMREX_SPACEDIM*2>,AMREX_SPACEDIM+NBCVAR_max> m_bc_neumann_vals;
    amrex::MultiFab* m_z_phys_nd;
    bool                 m_use_real_bcs;
    amrex::Real* m_v_bc_data;
};

class ERFPhysBCFunct_w
{
public:
    ERFPhysBCFunct_w (const int lev, const amrex::Geometry& geom,
                      const amrex::Vector<amrex::BCRec>& domain_bcs_type,
                      const amrex::Gpu::DeviceVector<amrex::BCRec>& domain_bcs_type_d,
                      amrex::Array<amrex::Array<amrex::Real,AMREX_SPACEDIM*2>,AMREX_SPACEDIM+NBCVAR_max> bc_extdir_vals,
                      amrex::Array<amrex::Array<amrex::Real,AMREX_SPACEDIM*2>,AMREX_SPACEDIM+NBCVAR_max> bc_neumann_vals,
                      const TerrainType& terrain_type, std::unique_ptr<amrex::MultiFab>& z_phys_nd,
                      const bool use_real_bcs,
                      amrex::Real* w_bc_data)
        : m_lev(lev),
          m_geom(geom), m_domain_bcs_type(domain_bcs_type),
          m_domain_bcs_type_d(domain_bcs_type_d),
          m_bc_extdir_vals(bc_extdir_vals),
          m_bc_neumann_vals(bc_neumann_vals),
          m_terrain_type(terrain_type),
          m_z_phys_nd(z_phys_nd.get()),
          m_use_real_bcs(use_real_bcs),
          m_w_bc_data(w_bc_data)
    { }

    ~ERFPhysBCFunct_w () {}

   /*
    * Impose physical boundary conditions at domain boundaries
    *
    * @param[out] mfs         Vector of MultiFabs to be filled containing, in order: cons, xvel, yvel, and zvel data
    * @param[in] icomp_cons   starting component for conserved variables
    * @param[in] ncomp_cons   number of components for conserved variables
    * @param[in] nghost       number of ghost cells to be filled for velocity components
    * @param[in] time         time at which the data should be filled
    * @param[in] use_real_bcs if true then we fill boundary conditions for interior locations
    */
    void operator() (amrex::MultiFab& mf, amrex::MultiFab& xvel, amrex::MultiFab& yvel,
                     amrex::IntVect const& nghost, const amrex::Real time, int bccomp);

    void impose_lateral_zvel_bcs (const amrex::Array4<amrex::Real>& dest_arr,
                                  const amrex::Array4<amrex::Real const>& xvel_arr,
                                  const amrex::Array4<amrex::Real const>& yvel_arr,
                                  const amrex::Box& bx, const amrex::Box& domain,
                                  const amrex::Array4<amrex::Real const>& z_nd,
                                  const amrex::GpuArray<amrex::Real,AMREX_SPACEDIM> dxInv,
                                  int bccomp);
    void impose_vertical_zvel_bcs (const amrex::Array4<amrex::Real>& dest_arr,
                                   const amrex::Array4<amrex::Real const>& xvel_arr,
                                   const amrex::Array4<amrex::Real const>& yvel_arr,
                                   const amrex::Box& bx, const amrex::Box& domain,
                                   const amrex::Array4<amrex::Real const>& z_nd,
                                   const amrex::GpuArray<amrex::Real,AMREX_SPACEDIM> dxInv,
                                   int bccomp_u, int bccomp_v, int bccomp_w,
                                   TerrainType terrain_type);

private:
    int                  m_lev;
    amrex::Geometry      m_geom;
    amrex::Vector<amrex::BCRec>            m_domain_bcs_type;
    amrex::Gpu::DeviceVector<amrex::BCRec> m_domain_bcs_type_d;
    amrex::Array<amrex::Array<amrex::Real, AMREX_SPACEDIM*2>,AMREX_SPACEDIM+NBCVAR_max> m_bc_extdir_vals;
    amrex::Array<amrex::Array<amrex::Real, AMREX_SPACEDIM*2>,AMREX_SPACEDIM+NBCVAR_max> m_bc_neumann_vals;
    TerrainType          m_terrain_type;
    amrex::MultiFab* m_z_phys_nd;
    bool                 m_use_real_bcs;
    amrex::Real* m_w_bc_data;
};

class ERFPhysBCFunct_w_no_terrain
{
public:
    ERFPhysBCFunct_w_no_terrain (
                      const int lev, const amrex::Geometry& geom,
                      const amrex::Vector<amrex::BCRec>& domain_bcs_type,
                      const amrex::Gpu::DeviceVector<amrex::BCRec>& domain_bcs_type_d,
                      amrex::Array<amrex::Array<amrex::Real,AMREX_SPACEDIM*2>,AMREX_SPACEDIM+NBCVAR_max> bc_extdir_vals,
                      amrex::Array<amrex::Array<amrex::Real,AMREX_SPACEDIM*2>,AMREX_SPACEDIM+NBCVAR_max> bc_neumann_vals,
                      const bool use_real_bcs,
                      amrex::Real* w_bc_data)
        : m_lev(lev),
          m_geom(geom), m_domain_bcs_type(domain_bcs_type),
          m_domain_bcs_type_d(domain_bcs_type_d),
          m_bc_extdir_vals(bc_extdir_vals),
          m_bc_neumann_vals(bc_neumann_vals),
          m_use_real_bcs(use_real_bcs),
          m_w_bc_data(w_bc_data)
    { }

    ~ERFPhysBCFunct_w_no_terrain () {}

   /*
    * Impose physical boundary conditions at domain boundaries AS IF there is no terrain
    *
    * @param[out] mfs         Vector of MultiFabs to be filled containing, in order: cons, xvel, yvel, and zvel data
    * @param[in] icomp_cons   starting component for conserved variables
    * @param[in] ncomp_cons   number of components for conserved variables
    * @param[in] nghost       number of ghost cells to be filled for velocity components
    * @param[in] time         time at which the data should be filled
    * @param[in] use_real_bcs if true then we fill boundary conditions for interior locations
    */
    void operator() (amrex::MultiFab& mf, int icomp, int ncomp,
                     amrex::IntVect const& nghost, const amrex::Real time, int bccomp_w);

    void impose_lateral_zvel_bcs (const amrex::Array4<amrex::Real>& dest_arr,
                                  const amrex::Box& bx, const amrex::Box& domain,
                                  int bccomp_w);
    void impose_vertical_zvel_bcs (const amrex::Array4<amrex::Real>& dest_arr,
                                   const amrex::Box& bx, const amrex::Box& domain,
                                   int bccomp_w);

private:
    int                  m_lev;
    amrex::Geometry      m_geom;
    amrex::Vector<amrex::BCRec>            m_domain_bcs_type;
    amrex::Gpu::DeviceVector<amrex::BCRec> m_domain_bcs_type_d;
    amrex::Array<amrex::Array<amrex::Real, AMREX_SPACEDIM*2>,AMREX_SPACEDIM+NBCVAR_max> m_bc_extdir_vals;
    amrex::Array<amrex::Array<amrex::Real, AMREX_SPACEDIM*2>,AMREX_SPACEDIM+NBCVAR_max> m_bc_neumann_vals;
    bool                 m_use_real_bcs;
    amrex::Real* m_w_bc_data;
};

#endif
