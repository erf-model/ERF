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
#include <AMReX_TypeTraits.H>

#include <TimeInterpolatedData.H>
#include <IndexDefines.H>

struct NullFill
{
    AMREX_GPU_DEVICE
    void operator() (const IntVect& /*iv*/, Array4<Real> const& /*dest*/,
                     const int /*dcomp*/, const int /*numcomp*/,
                     GeometryData const& /*geom*/, const Real /*time*/,
                     const BCRec* /*bcr*/, const int /*bcomp*/,
                     const int /*orig_comp*/) const
        {
        }
};

class ERFPhysBCFunct
{
public:
    ERFPhysBCFunct (const Geometry& geom, const Vector<BCRec>& bcr, const int var_idx, TimeInterpolatedData& data,
                    GpuArray<GpuArray<Real, NVAR          >,AMREX_SPACEDIM*2 >bc_cons_extdir, 
                    GpuArray<GpuArray<Real, AMREX_SPACEDIM>,AMREX_SPACEDIM*2 >bc_vels_extdir)
        : m_geom(geom), m_var_idx(var_idx), m_bcr(bcr), m_data(data), m_bc_cons(bc_cons_extdir), m_bc_vels(bc_vels_extdir)
        {}

    ~ERFPhysBCFunct () {}

    void operator() (MultiFab& mf, int icomp, int ncomp, IntVect const& nghost,
                     Real time, int bccomp)
    {
        if (m_geom.isAllPeriodic()) return;

        BL_PROFILE("ERFPhysBCFunct::()");

        const auto& domain = m_geom.Domain();
        const auto& dom_lo = amrex::lbound(domain);
        const auto& dom_hi = amrex::ubound(domain);

        //! create a grown domain box containing valid + periodic cells
        Box gdomain = amrex::convert(domain, mf.boxArray().ixType());
        for (int i = 0; i < AMREX_SPACEDIM; ++i) {
            if (m_geom.isPeriodic(i)) {
                gdomain.grow(i, nghost[i]);
            }
        }

#ifdef AMREX_USE_OMP
#pragma omp parallel if (Gpu::notInLaunchRegion())
#endif
        {
            Vector<BCRec> bcrs(ncomp);
            for (MFIter mfi(mf); mfi.isValid(); ++mfi)
            {
                FArrayBox& dest = mf[mfi];
                auto dest_array = dest.array();
                const Box& bx = mfi.fabbox();

                // make a cell-centered box for our state data
                const Box& ccbx = m_data.get_var(Vars::cons)[mfi].box();

                // make a vector of Array4's for our cell-centered & velocity data
                Vector<Array4<const Real> > vars_arrays_h;
                for (int i = 0; i < m_data.num_vars(); ++i) {
                    vars_arrays_h.push_back(m_data.get_var(i)[mfi].array());
                }
                amrex::Gpu::DeviceVector<Array4<const Real> > vars_arrays_d;
                vars_arrays_d.resize(vars_arrays_h.size());
                amrex::Gpu::copy(amrex::Gpu::hostToDevice, vars_arrays_h.begin(), vars_arrays_h.end(), vars_arrays_d.begin());
                Array4<const Real>* vars_arrays_p = vars_arrays_d.data();

                //! if there are cells not in the valid + periodic grown box
                //! we need to fill them here
                //!
                if (!gdomain.contains(bx))
                {
                    //! Based on BCRec for the domain, we need to make BCRec for this Box
                    amrex::setBC(bx, domain, bccomp, 0, ncomp, m_bcr, bcrs);

                    // Call the default fill functions
                    //! Note that we pass 0 as starting component of bcrs.
                    GpuBndryFuncFab<NullFill> bndry_fill_cc_fc_nd(NullFill{});

                    // Calls routines to fill all the foextrap, hoextrap, etc types of bc's
                    bndry_fill_cc_fc_nd(bx, dest, icomp, ncomp, m_geom, time, bcrs, 0, bccomp);

                    // xlo: ori = 0
                    // ylo: ori = 1
                    // zlo: ori = 2
                    // xhi: ori = 3
                    // yhi: ori = 4
                    // zhi: ori = 5

                    // Fill here all the "generic" ext_dir bc's
                    ParallelFor(bx, ncomp, [=] AMREX_GPU_DEVICE (int i, int j, int k, int n) 
                    {
                        if (m_var_idx == Vars::cons) {
                            if (i < dom_lo.x && bcrs[n].lo(0) == ERFBCType::ext_dir)
                                dest_array(i,j,k,n) = m_bc_cons[0][n];
                            if (j < dom_lo.y && bcrs[n].lo(1) == ERFBCType::ext_dir)
                                dest_array(i,j,k,n) = m_bc_cons[1][n];
                            if (k < dom_lo.z && bcrs[n].lo(2) == ERFBCType::ext_dir)
                                dest_array(i,j,k,n) = m_bc_cons[2][n];
                            if (i > dom_hi.x && bcrs[n].hi(0) == ERFBCType::ext_dir)
                                dest_array(i,j,k,n) = m_bc_cons[3][n];
                            if (j > dom_hi.y && bcrs[n].hi(1) == ERFBCType::ext_dir) 
                                dest_array(i,j,k,n) = m_bc_cons[4][n];
                            if (k > dom_hi.z && bcrs[n].hi(2) == ERFBCType::ext_dir)
                                dest_array(i,j,k,n) = m_bc_cons[5][n];

                        } else if (m_var_idx == Vars::xvel) {
                            if (i <= dom_lo.x   && bcrs[n].lo(0) == ERFBCType::ext_dir)
                                dest_array(i,j,k,n) = m_bc_vels[0][0];
                            if (j < dom_lo.y   && bcrs[n].lo(1) == ERFBCType::ext_dir)
                                dest_array(i,j,k,n) = m_bc_vels[1][0];
                            if (k < dom_lo.z   && bcrs[n].lo(2) == ERFBCType::ext_dir)
                                dest_array(i,j,k,n) = m_bc_vels[2][0];
                            if (i > dom_hi.x   && bcrs[n].hi(0) == ERFBCType::ext_dir)
                                dest_array(i,j,k,n) = m_bc_vels[3][0];
                            if (j > dom_hi.y   && bcrs[n].hi(1) == ERFBCType::ext_dir) 
                                dest_array(i,j,k,n) = m_bc_vels[4][0];
                            if (k > dom_hi.z   && bcrs[n].hi(2) == ERFBCType::ext_dir)
                                dest_array(i,j,k,n) = m_bc_vels[5][0];

                        } else if (m_var_idx == Vars::yvel) {
                            if (i < dom_lo.x   && bcrs[n].lo(0) == ERFBCType::ext_dir)
                                dest_array(i,j,k,n) = m_bc_vels[0][1];
                            if (j <= dom_lo.y   && bcrs[n].lo(1) == ERFBCType::ext_dir)
                                dest_array(i,j,k,n) = m_bc_vels[1][1];
                            if (k < dom_lo.z   && bcrs[n].lo(2) == ERFBCType::ext_dir)
                                dest_array(i,j,k,n) = m_bc_vels[2][1];
                            if (i > dom_hi.x   && bcrs[n].hi(0) == ERFBCType::ext_dir)
                                dest_array(i,j,k,n) = m_bc_vels[3][1];
                            if (j > dom_hi.y   && bcrs[n].hi(1) == ERFBCType::ext_dir) 
                                dest_array(i,j,k,n) = m_bc_vels[4][1];
                            if (k > dom_hi.z   && bcrs[n].hi(2) == ERFBCType::ext_dir)
                                dest_array(i,j,k,n) = m_bc_vels[5][1];


                        } else if (m_var_idx == Vars::zvel) {
                            if (i < dom_lo.x   && bcrs[n].lo(0) == ERFBCType::ext_dir)
                                dest_array(i,j,k,n) = m_bc_vels[0][2];
                            if (j < dom_lo.y   && bcrs[n].lo(1) == ERFBCType::ext_dir)
                                dest_array(i,j,k,n) = m_bc_vels[1][2];
                            if (k <= dom_lo.z   && bcrs[n].lo(2) == ERFBCType::ext_dir)
                                dest_array(i,j,k,n) = m_bc_vels[2][2];
                            if (i > dom_hi.x   && bcrs[n].hi(0) == ERFBCType::ext_dir)
                                dest_array(i,j,k,n) = m_bc_vels[3][2];
                            if (j > dom_hi.y   && bcrs[n].hi(1) == ERFBCType::ext_dir) 
                                dest_array(i,j,k,n) = m_bc_vels[4][2];
                            if (k > dom_hi.z && bcrs[n].hi(2) == ERFBCType::ext_dir)
                                dest_array(i,j,k,n) = m_bc_vels[5][2];
                        }
                    });

/* EXAMPLE CUSTOM PHYSICAL BC FUNCTION USING ALL VARIABLES
                    // Call our custom BC fill functions here using ccbx and vars_arrays to pass
                    // the cell-centered box and all the state variables. Note that var_idx tells us which
                    // of Vars::cons, Vars::xvel, Vars::yvel, or Vars::zvel we are filling in case we need logic.
                    ParallelFor(bx, [=] AMREX_GPU_DEVICE (int i, int j, int k) {
                        // we probably want some logic to compare (i,j,k) to the domain bounds here ...
                        if (var_idx == Vars::cons) { // depending on which variable type we are filling
                            dest_array(i,j,k,0) = vars_arrays_p[Vars::xvel](i,j,k); // (use other variable types, for example)
                        }
                });
*/
                }
            }
        }
    }

    // For backward compatibility
    void FillBoundary (MultiFab& mf, int dcomp, int ncomp, IntVect const& nghost,
                       Real time, int bccomp) {
        this->operator()(mf,dcomp,ncomp,nghost,time,bccomp);
    }

private:
    Geometry      m_geom;
    Vector<BCRec> m_bcr;
    int           m_var_idx;
    TimeInterpolatedData& m_data;
    GpuArray<GpuArray<Real, NVAR          >,AMREX_SPACEDIM*2 >m_bc_cons;
    GpuArray<GpuArray<Real, AMREX_SPACEDIM>,AMREX_SPACEDIM*2 >m_bc_vels;
};

#endif
