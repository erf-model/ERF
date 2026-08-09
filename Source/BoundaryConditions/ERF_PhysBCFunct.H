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
#include <ERF_TimeInterpolatedData.H>
#include <ERF_IndexDefines.H>
#include <ERF_DataStruct.H>
#include <ERF_EddyViscosity.H>
#include <ERF_TerrainMetrics.H>

class ERFPhysBCFunct_cons
{
public:
    /**
     * Construct the conserved-variable physical boundary-condition functor.
     *
     * @param[in] lev level index
     * @param[in] geom geometry for the level
     * @param[in] domain_bcs_type host boundary-condition records
     * @param[in] domain_bcs_type_d device boundary-condition records
     * @param[in] bc_extdir_vals external Dirichlet values for boundary faces
     * @param[in] bc_neumann_vals Neumann values for boundary faces
     * @param[in] z_phys_nd nodal physical height MultiFab
     * @param[in] use_real_bcs whether real-data boundary conditions are active
     * @param[in] th_bc_data optional theta boundary data
     */
    ERFPhysBCFunct_cons (const int lev,
                         const amrex::Geometry& geom, const amrex::Vector<amrex::BCRec>& domain_bcs_type,
                         const amrex::Gpu::DeviceVector<amrex::BCRec>& domain_bcs_type_d,
                         amrex::Array<amrex::Array<amrex::Real,AMREX_SPACEDIM*2>,AMREX_SPACEDIM+NBCVAR_max> bc_extdir_vals,
                         amrex::Array<amrex::Array<amrex::Real,AMREX_SPACEDIM*2>,AMREX_SPACEDIM+NBCVAR_max> bc_neumann_vals,
                         std::unique_ptr<amrex::MultiFab>& z_phys_nd,
                         const bool use_real_bcs,
                         amrex::Real* th_bc_data)
        : m_lev(lev), m_geom(geom),
          m_domain_bcs_type(domain_bcs_type),
          m_domain_bcs_type_d(domain_bcs_type_d),
          m_bc_extdir_vals(bc_extdir_vals),
          m_bc_neumann_vals(bc_neumann_vals),
          m_z_phys_nd(z_phys_nd.get()),
          m_use_real_bcs(use_real_bcs),
          m_th_bc_data(th_bc_data)
    {}

    /**
     * Destroy the conserved-variable physical boundary-condition functor.
     */
    ~ERFPhysBCFunct_cons () {}

   /**
    * Impose physical boundary conditions at domain boundaries
    *
    * @param[in,out] mf                    conserved-variable MultiFab to be filled
    * @param[in]     xvel                  x-velocity MultiFab used by upwind boundary conditions
    * @param[in]     yvel                  y-velocity MultiFab used by upwind boundary conditions
    * @param[in]     icomp                 starting component for conserved variables
    * @param[in]     ncomp                 number of components for conserved variables
    * @param[in]     nghost                number of ghost cells to be filled for conserved variables
    * @param[in]     time                  time at which the data should be filled
    * @param[in]     bccomp_cons           index into m_domain_bcs_type for conserved variables
    * @param[in]     do_fb                 whether to fill periodic/interior ghost cells first
    * @param[in]     do_terrain_adjustment whether to apply terrain-aware Neumann adjustments
    */
    void operator() (amrex::MultiFab& mf, amrex::MultiFab& xvel, amrex::MultiFab& yvel,
                     int icomp, int ncomp, amrex::IntVect const& nghost,
                     const double time, int bccomp_cons, bool do_fb = true, bool do_terrain_adjustment = true);

    /**
     * Impose lateral boundary conditions on conserved variables.
     *
     * @param[in,out] dest_arr cell-centered data to be filled
     * @param[in]     xvel_arr x-velocity used to determine upwind inflow
     * @param[in]     yvel_arr y-velocity used to determine upwind inflow
     * @param[in]     bx       box holding data to be filled
     * @param[in]     domain   simulation domain
     * @param[in]     icomp    starting component
     * @param[in]     ncomp    number of components
     * @param[in]     ng       number of ghost cells
     * @param[in]     time     time at which data should be filled
     */
    void impose_lateral_cons_bcs (const amrex::Array4<amrex::Real>& dest_arr,
                                  const amrex::Array4<amrex::Real const>& xvel_arr,
                                  const amrex::Array4<amrex::Real const>& yvel_arr,
                                  const amrex::Box& bx, const amrex::Box& domain,
                                  int icomp, int ncomp, amrex::IntVect ng,
                                  const double time);

    /**
     * Impose vertical boundary conditions on conserved variables.
     *
     * @param[in,out] dest_arr              cell-centered data to be filled
     * @param[in]     bx                    box holding data to be filled
     * @param[in]     domain                simulation domain
     * @param[in]     z_nd                  nodal physical height
     * @param[in]     dxInv                 inverse cell size
     * @param[in]     icomp                 starting component
     * @param[in]     ncomp                 number of components
     * @param[in]     time                  time at which data should be filled
     * @param[in]     do_terrain_adjustment whether to apply terrain-aware Neumann adjustments
     */
    void impose_vertical_cons_bcs (const amrex::Array4<amrex::Real>& dest_arr,
                                   const amrex::Box& bx, const amrex::Box& domain,
                                   const amrex::Array4<amrex::Real const>& z_nd,
                                   const amrex::GpuArray<amrex::Real,AMREX_SPACEDIM> dxInv,
                                   int icomp, int ncomp, const double time,
                                   bool do_terrain_adjustment = true);

private:
    int                  m_lev;
    amrex::Geometry      m_geom;
    amrex::Vector<amrex::BCRec>            m_domain_bcs_type;
    amrex::Gpu::DeviceVector<amrex::BCRec> m_domain_bcs_type_d;
    amrex::Array<amrex::Array<amrex::Real, AMREX_SPACEDIM*2>,AMREX_SPACEDIM+NBCVAR_max> m_bc_extdir_vals;
    amrex::Array<amrex::Array<amrex::Real, AMREX_SPACEDIM*2>,AMREX_SPACEDIM+NBCVAR_max> m_bc_neumann_vals;
    amrex::MultiFab* m_z_phys_nd;
    bool             m_use_real_bcs;
    amrex::Real*     m_th_bc_data;
};

class ERFPhysBCFunct_u
{
public:
    /**
     * Construct the x-velocity physical boundary-condition functor.
     *
     * @param[in] lev level index
     * @param[in] geom geometry for the level
     * @param[in] domain_bcs_type host boundary-condition records
     * @param[in] domain_bcs_type_d device boundary-condition records
     * @param[in] bc_extdir_vals external Dirichlet values for boundary faces
     * @param[in] bc_neumann_vals Neumann values for boundary faces
     * @param[in] z_phys_nd nodal physical height MultiFab
     * @param[in] use_real_bcs whether real-data boundary conditions are active
     * @param[in] u_bc_data optional x-velocity boundary data
     */
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

    /**
     * Destroy the x-velocity physical boundary-condition functor.
     */
    ~ERFPhysBCFunct_u () {}

   /**
    * Impose physical boundary conditions at domain boundaries
    *
    * @param[in,out] mf      x-velocity MultiFab to be filled
    * @param[in]     xvel    x-velocity MultiFab used by upwind boundary conditions
    * @param[in]     yvel    y-velocity MultiFab used by upwind boundary conditions
    * @param[in]     nghost  number of ghost cells to be filled for velocity components
    * @param[in]     time    time at which the data should be filled
    * @param[in]     bccomp  index into m_domain_bcs_type
    * @param[in]     do_fb   whether to fill periodic/interior ghost cells first
    */
    void operator() (amrex::MultiFab& mf, amrex::MultiFab& xvel, amrex::MultiFab& yvel,
                     amrex::IntVect const& nghost, const double time, int bccomp, bool do_fb);

    /**
     * Impose lateral boundary conditions on x velocity.
     *
     * @param[in,out] dest_arr x-velocity data to be filled
     * @param[in]     xvel_arr x-velocity used to determine upwind inflow
     * @param[in]     yvel_arr y-velocity used to determine upwind inflow
     * @param[in]     bx       box holding data to be filled
     * @param[in]     domain   simulation domain
     * @param[in]     bccomp   boundary-condition component
     * @param[in]     time     time at which data should be filled
     */
    void impose_lateral_xvel_bcs (const amrex::Array4<amrex::Real>& dest_arr,
                                  const amrex::Array4<amrex::Real const>& xvel_arr,
                                  const amrex::Array4<amrex::Real const>& yvel_arr,
                                  const amrex::Box& bx, const amrex::Box& domain,
                                  int bccomp,
                                  const double time);

    /**
     * Impose vertical boundary conditions on x velocity.
     *
     * @param[in,out] dest_arr x-velocity data to be filled
     * @param[in]     bx       box holding data to be filled
     * @param[in]     domain   simulation domain
     * @param[in]     z_nd     nodal physical height
     * @param[in]     dxInv    inverse cell size
     * @param[in]     bccomp   boundary-condition component
     * @param[in]     time     time at which data should be filled
     */
    void impose_vertical_xvel_bcs (const amrex::Array4<amrex::Real>& dest_arr,
                                   const amrex::Box& bx, const amrex::Box& domain,
                                   const amrex::Array4<amrex::Real const>& z_nd,
                                   const amrex::GpuArray<amrex::Real,AMREX_SPACEDIM> dxInv,
                                   int bccomp,
                                   const double time);

private:
    int                  m_lev;
    amrex::Geometry      m_geom;
    amrex::Vector<amrex::BCRec>            m_domain_bcs_type;
    amrex::Gpu::DeviceVector<amrex::BCRec> m_domain_bcs_type_d;
    amrex::Array<amrex::Array<amrex::Real, AMREX_SPACEDIM*2>,AMREX_SPACEDIM+NBCVAR_max> m_bc_extdir_vals;
    amrex::Array<amrex::Array<amrex::Real, AMREX_SPACEDIM*2>,AMREX_SPACEDIM+NBCVAR_max> m_bc_neumann_vals;
    amrex::MultiFab* m_z_phys_nd;
    bool             m_use_real_bcs;
    amrex::Real*     m_u_bc_data;
};

class ERFPhysBCFunct_v
{
public:
    /**
     * Construct the y-velocity physical boundary-condition functor.
     *
     * @param[in] lev level index
     * @param[in] geom geometry for the level
     * @param[in] domain_bcs_type host boundary-condition records
     * @param[in] domain_bcs_type_d device boundary-condition records
     * @param[in] bc_extdir_vals external Dirichlet values for boundary faces
     * @param[in] bc_neumann_vals Neumann values for boundary faces
     * @param[in] z_phys_nd nodal physical height MultiFab
     * @param[in] use_real_bcs whether real-data boundary conditions are active
     * @param[in] v_bc_data optional y-velocity boundary data
     */
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

    /**
     * Destroy the y-velocity physical boundary-condition functor.
     */
    ~ERFPhysBCFunct_v () {}

   /**
    * Impose physical boundary conditions at domain boundaries
    *
    * @param[in,out] mf      y-velocity MultiFab to be filled
    * @param[in]     xvel    x-velocity MultiFab used by upwind boundary conditions
    * @param[in]     yvel    y-velocity MultiFab used by upwind boundary conditions
    * @param[in]     nghost  number of ghost cells to be filled for velocity components
    * @param[in]     time    time at which the data should be filled
    * @param[in]     bccomp  index into m_domain_bcs_type
    * @param[in]     do_fb   whether to fill periodic/interior ghost cells first
    */
    void operator() (amrex::MultiFab& mf, amrex::MultiFab& xvel, amrex::MultiFab& yvel,
                     amrex::IntVect const& nghost, const double time, int bccomp, bool do_fb);

    /**
     * Impose lateral boundary conditions on y velocity.
     *
     * @param[in,out] dest_arr y-velocity data to be filled
     * @param[in]     xvel_arr x-velocity used to determine upwind inflow
     * @param[in]     yvel_arr y-velocity used to determine upwind inflow
     * @param[in]     bx       box holding data to be filled
     * @param[in]     domain   simulation domain
     * @param[in]     bccomp   boundary-condition component
     * @param[in]     time     time at which data should be filled
     */
    void impose_lateral_yvel_bcs (const amrex::Array4<amrex::Real>& dest_arr,
                                  const amrex::Array4<amrex::Real const>& xvel_arr,
                                  const amrex::Array4<amrex::Real const>& yvel_arr,
                                  const amrex::Box& bx, const amrex::Box& domain,
                                  int bccomp,
                                  const double time);

    /**
     * Impose vertical boundary conditions on y velocity.
     *
     * @param[in,out] dest_arr y-velocity data to be filled
     * @param[in]     bx       box holding data to be filled
     * @param[in]     domain   simulation domain
     * @param[in]     z_nd     nodal physical height
     * @param[in]     dxInv    inverse cell size
     * @param[in]     bccomp   boundary-condition component
     * @param[in]     time     time at which data should be filled
     */
    void impose_vertical_yvel_bcs (const amrex::Array4<amrex::Real>& dest_arr,
                                   const amrex::Box& bx, const amrex::Box& domain,
                                   const amrex::Array4<amrex::Real const>& z_nd,
                                   const amrex::GpuArray<amrex::Real,AMREX_SPACEDIM> dxInv,
                                   int bccomp,
                                   const double time);

private:
    int                  m_lev;
    amrex::Geometry      m_geom;
    amrex::Vector<amrex::BCRec>            m_domain_bcs_type;
    amrex::Gpu::DeviceVector<amrex::BCRec> m_domain_bcs_type_d;
    amrex::Array<amrex::Array<amrex::Real, AMREX_SPACEDIM*2>,AMREX_SPACEDIM+NBCVAR_max> m_bc_extdir_vals;
    amrex::Array<amrex::Array<amrex::Real, AMREX_SPACEDIM*2>,AMREX_SPACEDIM+NBCVAR_max> m_bc_neumann_vals;
    amrex::MultiFab* m_z_phys_nd;
    bool             m_use_real_bcs;
    amrex::Real*     m_v_bc_data;
};

class ERFPhysBCFunct_w
{
public:
    /**
     * Construct the z-velocity physical boundary-condition functor.
     *
     * @param[in] lev level index
     * @param[in] geom geometry for the level
     * @param[in] domain_bcs_type host boundary-condition records
     * @param[in] domain_bcs_type_d device boundary-condition records
     * @param[in] bc_extdir_vals external Dirichlet values for boundary faces
     * @param[in] bc_neumann_vals Neumann values for boundary faces
     * @param[in] terrain_type terrain type for the level
     * @param[in] mapfac_lev map-factor MultiFabs for the level
     * @param[in] z_phys_nd nodal physical height MultiFab
     * @param[in] use_real_bcs whether real-data boundary conditions are active
     * @param[in] w_bc_data optional z-velocity boundary data
     */
    ERFPhysBCFunct_w (const int lev, const amrex::Geometry& geom,
                      const amrex::Vector<amrex::BCRec>& domain_bcs_type,
                      const amrex::Gpu::DeviceVector<amrex::BCRec>& domain_bcs_type_d,
                      amrex::Array<amrex::Array<amrex::Real,AMREX_SPACEDIM*2>,AMREX_SPACEDIM+NBCVAR_max> bc_extdir_vals,
                      amrex::Array<amrex::Array<amrex::Real,AMREX_SPACEDIM*2>,AMREX_SPACEDIM+NBCVAR_max> bc_neumann_vals,
                      const TerrainType& terrain_type,
                      amrex::Vector<std::unique_ptr<amrex::MultiFab>>& mapfac_lev,
                      std::unique_ptr<amrex::MultiFab>& z_phys_nd,
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
    {
          m_mapfac_u = mapfac_lev[MapFacType::u_x].get();
          m_mapfac_v = mapfac_lev[MapFacType::v_y].get();
    }

    /**
     * Destroy the z-velocity physical boundary-condition functor.
     */
    ~ERFPhysBCFunct_w () {}

   /**
    * Impose physical boundary conditions at domain boundaries
    *
    * @param[in,out] mf       z-velocity MultiFab to be filled
    * @param[in]     xvel     x-velocity MultiFab used by upwind boundary conditions
    * @param[in]     yvel     y-velocity MultiFab used by upwind boundary conditions
    * @param[in]     nghost   number of ghost cells to be filled for velocity components
    * @param[in]     time     time at which the data should be filled
    * @param[in]     bccomp   index into m_domain_bcs_type for z velocity
    * @param[in]     do_fb    whether to fill periodic/interior ghost cells first
    */
    void operator() (amrex::MultiFab& mf, amrex::MultiFab& xvel, amrex::MultiFab& yvel,
                     amrex::IntVect const& nghost, const double time, int bccomp,
                     bool do_fb);

    /**
     * Impose lateral boundary conditions on z velocity.
     *
     * @param[in,out] dest_arr     z-velocity data to be filled
     * @param[in]     xvel_arr     x-velocity used to determine upwind inflow
     * @param[in]     yvel_arr     y-velocity used to determine upwind inflow
     * @param[in]     bx           box holding data to be filled
     * @param[in]     domain       simulation domain
     * @param[in]     mf_u         x map-factor field
     * @param[in]     mf_v         y map-factor field
     * @param[in]     z_nd         nodal physical height
     * @param[in]     dxInv        inverse cell size
     * @param[in]     terrain_type terrain representation
     * @param[in]     bccomp       boundary-condition component
     * @param[in]     time         time at which data should be filled
     */
    void impose_lateral_zvel_bcs (const amrex::Array4<amrex::Real>& dest_arr,
                                  const amrex::Array4<amrex::Real const>& xvel_arr,
                                  const amrex::Array4<amrex::Real const>& yvel_arr,
                                  const amrex::Box& bx, const amrex::Box& domain,
                                  const amrex::Array4<amrex::Real const>& mf_u,
                                  const amrex::Array4<amrex::Real const>& mf_v,
                                  const amrex::Array4<amrex::Real const>& z_nd,
                                  const amrex::GpuArray<amrex::Real,AMREX_SPACEDIM> dxInv,
                                  TerrainType terrain_type, int bccomp,
                                  const double time);

    /**
     * Impose vertical boundary conditions on z velocity.
     *
     * @param[in,out] dest_arr     z-velocity data to be filled
     * @param[in]     xvel_arr     x-velocity used by terrain-aware vertical BCs
     * @param[in]     yvel_arr     y-velocity used by terrain-aware vertical BCs
     * @param[in]     bx           box holding data to be filled
     * @param[in]     domain       simulation domain
     * @param[in]     mf_u         x map-factor field
     * @param[in]     mf_v         y map-factor field
     * @param[in]     z_nd         nodal physical height
     * @param[in]     dxInv        inverse cell size
     * @param[in]     bccomp_u     x-velocity boundary-condition component
     * @param[in]     bccomp_v     y-velocity boundary-condition component
     * @param[in]     bccomp_w     z-velocity boundary-condition component
     * @param[in]     terrain_type terrain representation
     * @param[in]     time         time at which data should be filled
     */
    void impose_vertical_zvel_bcs (const amrex::Array4<amrex::Real>& dest_arr,
                                   const amrex::Array4<amrex::Real const>& xvel_arr,
                                   const amrex::Array4<amrex::Real const>& yvel_arr,
                                   const amrex::Box& bx, const amrex::Box& domain,
                                   const amrex::Array4<amrex::Real const>& mf_u,
                                   const amrex::Array4<amrex::Real const>& mf_v,
                                   const amrex::Array4<amrex::Real const>& z_nd,
                                   const amrex::GpuArray<amrex::Real,AMREX_SPACEDIM> dxInv,
                                   int bccomp_u, int bccomp_v, int bccomp_w,
                                   TerrainType terrain_type,
                                   const double time);

private:
    int                  m_lev;
    amrex::Geometry      m_geom;
    amrex::Vector<amrex::BCRec>            m_domain_bcs_type;
    amrex::Gpu::DeviceVector<amrex::BCRec> m_domain_bcs_type_d;
    amrex::Array<amrex::Array<amrex::Real, AMREX_SPACEDIM*2>,AMREX_SPACEDIM+NBCVAR_max> m_bc_extdir_vals;
    amrex::Array<amrex::Array<amrex::Real, AMREX_SPACEDIM*2>,AMREX_SPACEDIM+NBCVAR_max> m_bc_neumann_vals;
    TerrainType      m_terrain_type;
    amrex::MultiFab* m_mapfac_u;
    amrex::MultiFab* m_mapfac_v;
    amrex::MultiFab* m_z_phys_nd;
    bool             m_use_real_bcs;
    amrex::Real*     m_w_bc_data;
};

class ERFPhysBCFunct_base
{
public:
    /**
     * Construct the base-state physical boundary-condition functor.
     *
     * @param[in] lev level index
     * @param[in] geom geometry for the level
     * @param[in] domain_bcs_type host boundary-condition records
     * @param[in] domain_bcs_type_d device boundary-condition records
     * @param[in] z_phys_nd nodal physical height MultiFab
     * @param[in] moving_terrain whether moving terrain is active
     */
    ERFPhysBCFunct_base (const int lev,
                         const amrex::Geometry& geom, const amrex::Vector<amrex::BCRec>& domain_bcs_type,
                         const amrex::Gpu::DeviceVector<amrex::BCRec>& domain_bcs_type_d,
                         std::unique_ptr<amrex::MultiFab>& z_phys_nd,
                         bool moving_terrain)
        : m_lev(lev), m_geom(geom),
          m_moving_terrain(moving_terrain),
          m_domain_bcs_type(domain_bcs_type),
          m_domain_bcs_type_d(domain_bcs_type_d),
          m_z_phys_nd(z_phys_nd.get())
    {}

    /**
     * Destroy the base-state physical boundary-condition functor.
     */
    ~ERFPhysBCFunct_base () {}

   /**
    * Impose physical boundary conditions at domain boundaries
    *
    * @param[in,out] mf     base-state MultiFab to be filled
    * @param[in]     icomp  starting component for base-state variables
    * @param[in]     ncomp  number of components for base-state variables
    * @param[in]     nghost number of ghost cells to be filled
    */
    void operator() (amrex::MultiFab& mf, int icomp, int ncomp, amrex::IntVect const& nghost);

    /**
     * Impose lateral boundary conditions on base-state variables.
     *
     * @param[in,out] dest_arr base-state data to be filled
     * @param[in]     bx       box holding data to be filled
     * @param[in]     domain   simulation domain
     * @param[in]     ncomp    number of components
     * @param[in]     nghost   number of ghost cells to be filled
     */
    void impose_lateral_basestate_bcs (const amrex::Array4<amrex::Real>& dest_arr,
                                       const amrex::Box& bx, const amrex::Box& domain,
                                       int ncomp, const amrex::IntVect& nghost);

    /**
     * Impose vertical boundary conditions on base-state variables.
     *
     * @param[in,out] dest_arr base-state data to be filled
     * @param[in]     z_nd     nodal physical height
     * @param[in]     bx       box holding data to be filled
     * @param[in]     domain   simulation domain
     * @param[in]     ncomp    number of components
     * @param[in]     nghost   number of ghost cells to be filled
     */
    void impose_vertical_basestate_bcs (const amrex::Array4<amrex::Real>& dest_arr,
                                        const amrex::Array4<amrex::Real const>& z_nd,
                                        const amrex::Box& bx, const amrex::Box& domain,
                                        int ncomp, const amrex::IntVect& nghost);

private:
    int                  m_lev;
    amrex::Geometry      m_geom;
    bool                 m_moving_terrain;
    amrex::Vector<amrex::BCRec>            m_domain_bcs_type;
    amrex::Gpu::DeviceVector<amrex::BCRec> m_domain_bcs_type_d;
    amrex::MultiFab*     m_z_phys_nd;
};

#endif
