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
#include <PBLModels.H>
#include <ABLMost.H>
#ifdef ERF_USE_TERRAIN
#include <SpatialStencils.H>
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

class ERFPhysBCFunct
{
public:
    ERFPhysBCFunct (const amrex::Geometry& geom, const amrex::Vector<amrex::BCRec>& domain_bcs_type,
                    const amrex::Gpu::DeviceVector<amrex::BCRec>& domain_bcs_type_d,
                    const int& var_idx, TimeInterpolatedData& data,
                    amrex::Array<amrex::Array<amrex::Real,AMREX_SPACEDIM*2>,AMREX_SPACEDIM+NVAR> bc_extdir_vals,
                    const SolverChoice& solverChoice,
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
                    ABLMost& most, std::unique_ptr<ReadBndryPlanes>& r2d)
        : m_geom(geom), m_domain_bcs_type(domain_bcs_type), m_domain_bcs_type_d(domain_bcs_type_d),
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
          m_most(most), m_r2d(r2d)
        {
            m_has_most_bcs = false;
            for (const auto& bcr : domain_bcs_type) {
                for (amrex::OrientationIter oiter; oiter && !m_has_most_bcs; ++oiter) {
                    int side = (int) oiter();
                    if (bcr.data()[side] == static_cast<int>(ERFBCType::MOST)) {
                        m_has_most_bcs = true;
                    }
                }
                if (m_has_most_bcs) break;
            }

            if (m_has_most_bcs) {
                const auto& cons_mf = m_data.get_var(Vars::cons);
                m_viscosity.define(cons_mf.boxArray(), cons_mf.DistributionMap(), EddyDiff::NumDiffs, cons_mf.nGrowVect());
                AMREX_ALWAYS_ASSERT_WITH_MESSAGE(solverChoice.les_type == LESType::Smagorinsky ||
                                                 solverChoice.les_type == LESType::Deardorff ||,
                                                 solverChoice.pbl_type != PBLType::None,
                                                 "Must use an LES or PBL model to compute turbulent viscosity for MOST boundaries");
		// PBL or PBL+LES - get vertical viscosity component with PBL model. LES - get it with LES model
		if (solverChoice.pbl_type != PBLType::None) {
		  ComputeTurbulentViscosityPBL(m_data.get_var(Vars::xvel),
					       m_data.get_var(Vars::yvel),
					       m_data.get_var(Vars::zvel),
					       cons_mf, m_viscosity, geom, solverChoice, m_most);
		} else {
		  ComputeTurbulentViscosity(m_data.get_var(Vars::xvel),
					    m_data.get_var(Vars::yvel),
					    m_data.get_var(Vars::zvel),
					    cons_mf, m_viscosity, geom, solverChoice,
					    domain_bcs_type_d);
		}
            }
        }

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
    amrex::Geometry      m_geom;
    amrex::Vector<amrex::BCRec> m_domain_bcs_type;
    amrex::Gpu::DeviceVector<amrex::BCRec> m_domain_bcs_type_d;
    int           m_var_idx;
    TimeInterpolatedData& m_data;
    amrex::Array<amrex::Array<amrex::Real, AMREX_SPACEDIM*2>,AMREX_SPACEDIM+NVAR> m_bc_extdir_vals;
    amrex::MultiFab m_viscosity;
    bool m_has_most_bcs;
#ifdef ERF_USE_TERRAIN
    const amrex::MultiFab& m_z_phys_nd;
    const amrex::MultiFab& m_detJ_cc;
#endif
    ABLMost m_most;
    std::unique_ptr<ReadBndryPlanes>& m_r2d;
#ifdef ERF_USE_NETCDF
    const amrex::Vector<amrex::FArrayBox>& m_bdy_data_xlo;
    const amrex::Vector<amrex::FArrayBox>& m_bdy_data_xhi;
    const amrex::Vector<amrex::FArrayBox>& m_bdy_data_ylo;
    const amrex::Vector<amrex::FArrayBox>& m_bdy_data_yhi;
#endif
};

#endif
