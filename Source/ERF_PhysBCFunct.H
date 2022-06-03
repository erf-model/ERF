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

#ifdef ERF_USE_TERRAIN
#include <TerrainMetrics.H>
#endif

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

#ifdef ERF_USE_TERRAIN
struct TerrainBcFill
{
    const Array4<const Real>& m_z_nd;

    TerrainBcFill(const Array4<const Real>& z_nd): m_z_nd(z_nd){}

    AMREX_GPU_DEVICE
    void operator() (const amrex::IntVect& iv, amrex::Array4<amrex::Real> const& dest,
                     const int dcomp, const int numcomp,
                     amrex::GeometryData const& geom, const amrex::Real time,
                     const amrex::BCRec* bcr, const int bcomp,
                     const int orig_comp) const
    {
        // This needs re-work for multicore. The box owned by the process
        // needs sent to this routine to avoid segfault with derivatives.
        return;

        // Indices of parallelfor
        const int i = iv[0]; const int j = iv[1]; const int k = iv[2];

        // Inverse cell size
        amrex::GpuArray<amrex::Real, AMREX_SPACEDIM> dx,dxInv;
        dx[0] = geom.CellSize(0);
        dx[1] = geom.CellSize(1);
        dx[2] = geom.CellSize(2);
        dxInv[0] = 1. / dx[0];
        dxInv[1] = 1. / dx[1];
        dxInv[2] = 1. / dx[2];

        // Domain box indices
        const amrex::Box& domain = geom.Domain();
        int ilo = domain.smallEnd(0);
        int jlo = domain.smallEnd(1);
        int klo = domain.smallEnd(2);
        int ihi = domain.bigEnd(0);
        int jhi = domain.bigEnd(1);
        int khi = domain.bigEnd(2);

        // Metric terms
        amrex::Real met_h_xi,met_h_eta,met_h_zeta;

        // Modify BC on scalar and x/y-velocity
        if(orig_comp <= BCVars::yvel_bc) {

            // Set nodal bounds
            if(orig_comp == BCVars::yvel_bc) {
                jhi += 1;
            } else if (orig_comp == BCVars::xvel_bc) {
                ihi += 1;
            }

            // Indices for metrics
            int im = (i  < ilo) ? ilo : i;
                im = (im > ihi) ? ihi : im;
            int jm = (j  < jlo) ? jlo : j;
                jm = (jm > jhi) ? jhi : jm;
            int km = (k  < klo) ? klo : k;
                km = (km > khi) ? khi : km;

            // Only modify bottom wall
            if(k < klo) {
                // Get metrics
                if (orig_comp == BCVars::yvel_bc) {
                    ComputeMetricAtJface(im,jm,km,met_h_xi,met_h_eta,met_h_zeta,dxInv,m_z_nd,TerrainMet::all);
                } else if (orig_comp == BCVars::xvel_bc) {
                    ComputeMetricAtIface(im,jm,km,met_h_xi,met_h_eta,met_h_zeta,dxInv,m_z_nd,TerrainMet::all);
                } else {
                    ComputeMetricAtCellCenter(im,jm,km,met_h_xi,met_h_eta,met_h_zeta,dxInv,m_z_nd,TerrainMet::all);
                }

                // GradX at IJK location inside domain
                amrex::Real GradVarx = 0.5 * dxInv[0] * (dest(im+1,jm,km) - dest(im-1,jm,km));

                // GradY at IJK location inside domain
                amrex::Real GradVary = 0.5 * dxInv[1] * (dest(im,jm+1,km) - dest(im,jm-1,km));

                // Prefactor
                amrex::Real met_fac =  met_h_zeta / ( met_h_xi*met_h_xi + met_h_eta*met_h_eta + 1. );

                // Accumulate in bottom ghost cell (Neumann condition along surface normal)
                dest(i,j,k) -= dx[2] * met_fac * ( met_h_xi * GradVarx + met_h_eta * GradVary );
            } // klo
        } // y_vel
    } // operator
};
#endif

class ERFPhysBCFunct
{
public:
    ERFPhysBCFunct (const int lev,
                    const amrex::Geometry& geom, const amrex::Vector<amrex::BCRec>& domain_bcs_type,
                    const amrex::Gpu::DeviceVector<amrex::BCRec>& domain_bcs_type_d,
                    const int& var_idx, TimeInterpolatedData& data,
                    amrex::Array<amrex::Array<amrex::Real,AMREX_SPACEDIM*2>,AMREX_SPACEDIM+NVAR> bc_extdir_vals,
#ifdef ERF_USE_TERRAIN
                    const amrex::MultiFab& z_phys_nd,
                    const amrex::MultiFab& detJ_cc,
#endif
#ifdef ERF_USE_NETCDF
                    const amrex::Vector<amrex::FArrayBox>& bdy_data_xlo,
                    const amrex::Vector<amrex::FArrayBox>& bdy_data_xhi,
                    const amrex::Vector<amrex::FArrayBox>& bdy_data_ylo,
                    const amrex::Vector<amrex::FArrayBox>& bdy_data_yhi,
#endif
                    std::unique_ptr<ReadBndryPlanes>& r2d)
        : m_lev(lev), m_geom(geom), m_domain_bcs_type(domain_bcs_type), m_domain_bcs_type_d(domain_bcs_type_d),
          m_var_idx(var_idx), m_data(data), m_bc_extdir_vals(bc_extdir_vals),
#ifdef ERF_USE_TERRAIN
          m_z_phys_nd(z_phys_nd), m_detJ_cc(detJ_cc),
#endif
#ifdef ERF_USE_NETCDF
          m_bdy_data_xlo(bdy_data_xlo),
          m_bdy_data_xhi(bdy_data_xhi),
          m_bdy_data_ylo(bdy_data_ylo),
          m_bdy_data_yhi(bdy_data_yhi),
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
    TimeInterpolatedData& m_data;
    amrex::Array<amrex::Array<amrex::Real, AMREX_SPACEDIM*2>,AMREX_SPACEDIM+NVAR> m_bc_extdir_vals;
#ifdef ERF_USE_TERRAIN
    const amrex::MultiFab& m_z_phys_nd;
    const amrex::MultiFab& m_detJ_cc;
#endif
    std::unique_ptr<ReadBndryPlanes>& m_r2d;
#ifdef ERF_USE_NETCDF
    const amrex::Vector<amrex::FArrayBox>& m_bdy_data_xlo;
    const amrex::Vector<amrex::FArrayBox>& m_bdy_data_xhi;
    const amrex::Vector<amrex::FArrayBox>& m_bdy_data_ylo;
    const amrex::Vector<amrex::FArrayBox>& m_bdy_data_yhi;
#endif
};

#endif
