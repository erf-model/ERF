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
#include <PBLModels.H>
#include <EddyViscosity.H>

#include <TerrainMetrics.H>

using namespace amrex;

using PlaneVector = amrex::Vector<amrex::FArrayBox>;

struct NullFill
{
    AMREX_GPU_DEVICE
    void operator() (const amrex::IntVect& /*iv*/, amrex::Array4<amrex::Real> const& /*dest*/,
                     const int /*dcomp*/, const int /*numcomp*/,
                     amrex::GeometryData const& /*geom*/, const amrex::Real /*time*/,
                     const amrex::BCRec* /*bcr*/, const int /*bcomp*/,
                     const int /*orig_comp*/) const
        {
        }
};

class ERFPhysBCFunct
{
public:
    ERFPhysBCFunct (const int lev,
                    const amrex::Geometry& geom, const amrex::Vector<amrex::BCRec>& domain_bcs_type,
                    const amrex::Gpu::DeviceVector<amrex::BCRec>& domain_bcs_type_d,
                    const int& var_idx, const int& terrain_type, TimeInterpolatedData& data,
                    amrex::Array<amrex::Array<amrex::Real,AMREX_SPACEDIM*2>,AMREX_SPACEDIM+NVAR> bc_extdir_vals,
                    std::unique_ptr<amrex::MultiFab>& z_phys_nd,
                    std::unique_ptr<amrex::MultiFab>& detJ_cc,
#ifdef ERF_USE_NETCDF
                    const std::string& init_type,
                    const amrex::Vector<amrex::Vector<amrex::FArrayBox>>& bdy_data_xlo,
                    const amrex::Vector<amrex::Vector<amrex::FArrayBox>>& bdy_data_xhi,
                    const amrex::Vector<amrex::Vector<amrex::FArrayBox>>& bdy_data_ylo,
                    const amrex::Vector<amrex::Vector<amrex::FArrayBox>>& bdy_data_yhi,
                    const amrex::Real bdy_time_interval,
#endif
                    std::unique_ptr<ReadBndryPlanes>& r2d)
        : m_lev(lev), m_geom(geom), m_domain_bcs_type(domain_bcs_type), m_domain_bcs_type_d(domain_bcs_type_d),
          m_var_idx(var_idx), m_terrain_type(terrain_type), m_data(data), m_bc_extdir_vals(bc_extdir_vals),
          m_z_phys_nd(z_phys_nd), m_detJ_cc(detJ_cc),

#ifdef ERF_USE_NETCDF
          m_init_type(init_type),
          m_bdy_data_xlo(bdy_data_xlo),
          m_bdy_data_xhi(bdy_data_xhi),
          m_bdy_data_ylo(bdy_data_ylo),
          m_bdy_data_yhi(bdy_data_yhi),
          m_bdy_time_interval(bdy_time_interval),
#endif
          m_r2d(r2d)
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

    void impose_xvel_bcs (const amrex::Array4<Real>& dest_arr, const amrex::Box& bx, const amrex::Box& domain,
                          const amrex::Array4<Real const>& z_nd,
                          const GpuArray<Real,AMREX_SPACEDIM> dxInv,
                          amrex::Real time, int bccomp);

    void impose_yvel_bcs (const amrex::Array4<Real>& dest_arr, const amrex::Box& bx, const amrex::Box& domain,
                          const amrex::Array4<Real const>& z_nd,
                          const GpuArray<Real,AMREX_SPACEDIM> dxInv,
                          amrex::Real time, int bccomp);

    void impose_zvel_bcs (const amrex::Array4<Real>& dest_arr, const amrex::Box& bx, const amrex::Box& domain,
                          const amrex::Array4<Real const>& velx_arr,
                          const amrex::Array4<Real const>& vely_arr,
                          const amrex::Array4<Real const>& z_nd,
                          const GpuArray<Real,AMREX_SPACEDIM> dx,
                          const GpuArray<Real,AMREX_SPACEDIM> dxInv,
                          amrex::Real time, int bccomp, int terrain_type);

    void impose_cons_bcs (const amrex::Array4<Real>& mf, const amrex::Box& bx, const amrex::Box& domain,
                          const amrex::Array4<Real const>& z_nd,
                          const GpuArray<Real,AMREX_SPACEDIM> dxInv,
                          int icomp, int ncomp, amrex::Real time, int bccomp);

    void fill_from_bndryregs (const int lev,   const amrex::Box& bx, const amrex::Array4<Real>& dest_arr,
                              const int icomp, const int bccomp, const int ncomp,
                              const amrex::Box& domain, const amrex::BCRec* bc_ptr,
                              const amrex::Real time);

#ifdef ERF_USE_NETCDF
    void fill_from_wrfbdy (const int lev,   const amrex::Box& bx, const amrex::Array4<Real>& dest_arr,
                           const int icomp, const int bccomp, const int ncomp,
                           const amrex::Box& domain, const amrex::BCRec* bc_ptr,
                           const amrex::Real time, const amrex::Real bdy_time_interval);
#endif

    // For backward compatibility
    void FillBoundary (amrex::MultiFab& mf, int dcomp, int ncomp, amrex::IntVect const& nghost,
                       amrex::Real time, int bccomp) {
        this->operator()(mf,dcomp,ncomp,nghost,time,bccomp);
    }

private:
    int                  m_lev;
    amrex::Geometry      m_geom;
    amrex::Vector<amrex::BCRec> m_domain_bcs_type;
    amrex::Gpu::DeviceVector<amrex::BCRec> m_domain_bcs_type_d;
    int           m_var_idx;
    int           m_terrain_type;
    TimeInterpolatedData& m_data;
    amrex::Array<amrex::Array<amrex::Real, AMREX_SPACEDIM*2>,AMREX_SPACEDIM+NVAR> m_bc_extdir_vals;
    std::unique_ptr<amrex::MultiFab>& m_z_phys_nd;
    std::unique_ptr<amrex::MultiFab>&   m_detJ_cc;
#ifdef ERF_USE_NETCDF
    const std::string m_init_type;
    const amrex::Vector<amrex::Vector<amrex::FArrayBox>>& m_bdy_data_xlo;
    const amrex::Vector<amrex::Vector<amrex::FArrayBox>>& m_bdy_data_xhi;
    const amrex::Vector<amrex::Vector<amrex::FArrayBox>>& m_bdy_data_ylo;
    const amrex::Vector<amrex::Vector<amrex::FArrayBox>>& m_bdy_data_yhi;
    const amrex::Real m_bdy_time_interval;
#endif
    std::unique_ptr<ReadBndryPlanes>& m_r2d;
};

#endif
