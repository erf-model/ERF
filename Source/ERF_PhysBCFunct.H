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
    ERFPhysBCFunct (const Geometry& geom, const Vector<BCRec>& domain_bcs_type_bcr,
                    const int& var_idx, TimeInterpolatedData& data,
                    Array<Array<Real,AMREX_SPACEDIM*2>,AMREX_SPACEDIM+NVAR> bc_extdir_vals)
        : m_geom(geom), m_var_idx(var_idx), m_domain_bcs_type_bcr(domain_bcs_type_bcr),
          m_data(data), m_bc_extdir_vals(bc_extdir_vals)
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
    void operator() (MultiFab& mf, int icomp, int ncomp, IntVect const& nghost,
                     Real time, int bccomp)
    {
        if (m_geom.isAllPeriodic()) return;

        BL_PROFILE("ERFPhysBCFunct::()");

        const auto& domain = m_geom.Domain();
        const auto& dom_lo = amrex::lbound(domain);
        const auto& dom_hi = amrex::ubound(domain);

        // Create a grown domain box containing valid + periodic cells
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
                const Array4<Real>& dest_array = mf.array(mfi);
                const Box& bx = mfi.fabbox();

                // make a cell-centered box for our state data
                //const Box& ccbx = m_data.get_var(Vars::cons)[mfi].box();

/*
                // make a vector of Array4's for our cell-centered & velocity data
                Gpu::DeviceVector<Array4<const Real> > vars_arrays_h;
                for (int i = 0; i < m_data.num_vars(); ++i) {
                    vars_arrays_h.push_back(m_data.get_var(i)[mfi].array());
                }
                Gpu::DeviceVector<Array4<const Real> > vars_arrays_d;
                vars_arrays_d.resize(vars_arrays_h.size());
                Gpu::copy(amrex::Gpu::hostToDevice, vars_arrays_h.begin(), vars_arrays_h.end(), vars_arrays_d.begin());
                Array4<const Real>* vars_arrays_p = vars_arrays_d.data();
*/

                //! if there are cells not in the valid + periodic grown box
                //! we need to fill them here
                //!
                if (!gdomain.contains(bx))
                {
                    //! Based on BCRec for the domain, we need to make BCRec for this Box
                    // bccomp is used as starting index for m_domain_bcs_type_bcr
                    //      0 is used as starting index for bcrs
                    amrex::setBC(bx, domain, bccomp, 0, ncomp, m_domain_bcs_type_bcr, bcrs);

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

                    amrex::Gpu::DeviceVector<BCRec> bcrs_d(ncomp);
#ifdef AMREX_USE_GPU
                    Gpu::htod_memcpy
                        (bcrs_d.data(), bcrs.data(), sizeof(BCRec)*ncomp);
#else
                    std::memcpy
                        (bcrs_d.data(), bcrs.data(), sizeof(BCRec)*ncomp);
#endif

                    if (icomp == BCVars::xvel_bc || icomp == BCVars::yvel_bc || icomp == BCVars::zvel_bc) {
                        AMREX_ALWAYS_ASSERT(ncomp == 1);
                    } else {
                        AMREX_ALWAYS_ASSERT(icomp+ncomp <= NVAR);
                    }

                    amrex::GpuArray<amrex::GpuArray<amrex::Real, AMREX_SPACEDIM*2>,
                                                                 AMREX_SPACEDIM+NVAR> l_bc_extdir_vals_d;

                    for (int i = 0; i < ncomp; i++)
                        for (int ori = 0; ori < 2*AMREX_SPACEDIM; ori++)
                            l_bc_extdir_vals_d[i][ori] = m_bc_extdir_vals[bccomp+i][ori];

                    const amrex::BCRec* bc_ptr = bcrs_d.data();

                    // Fill here all the "generic" ext_dir bc's
                    ParallelFor(bx, ncomp, [=] AMREX_GPU_DEVICE (int i, int j, int k, int n)
                    {
                        if (i < dom_lo.x && bc_ptr[n].lo(0) == ERFBCType::ext_dir)
                            dest_array(i,j,k,icomp+n) = l_bc_extdir_vals_d[n][0];
                        if (j < dom_lo.y && bc_ptr[n].lo(1) == ERFBCType::ext_dir)
                            dest_array(i,j,k,icomp+n) = l_bc_extdir_vals_d[n][1];
                        if (k < dom_lo.z && bc_ptr[n].lo(2) == ERFBCType::ext_dir)
                            dest_array(i,j,k,icomp+n) = l_bc_extdir_vals_d[n][2];
                        if (i > dom_hi.x && bc_ptr[n].hi(0) == ERFBCType::ext_dir)
                            dest_array(i,j,k,icomp+n) = l_bc_extdir_vals_d[n][3];
                        if (j > dom_hi.y && bc_ptr[n].hi(1) == ERFBCType::ext_dir)
                            dest_array(i,j,k,icomp+n) = l_bc_extdir_vals_d[n][4];
                        if (k > dom_hi.z && bc_ptr[n].hi(2) == ERFBCType::ext_dir)
                            dest_array(i,j,k,icomp+n) = l_bc_extdir_vals_d[n][5];

                        if (bccomp == BCVars::xvel_bc)
                        {
                            if (i == dom_lo.x && bc_ptr[n].lo(0) == ERFBCType::ext_dir)
                                dest_array(i,j,k,n) = l_bc_extdir_vals_d[n][0];
                        }
                        if (bccomp == BCVars::yvel_bc)
                        {
                            if (j == dom_lo.y && bc_ptr[n].lo(1) == ERFBCType::ext_dir)
                                dest_array(i,j,k,n) = l_bc_extdir_vals_d[n][1];
                        }
                        if (bccomp == BCVars::zvel_bc)
                        {
                            if (k == dom_lo.z && bc_ptr[n].lo(2) == ERFBCType::ext_dir)
                                dest_array(i,j,k,n) = l_bc_extdir_vals_d[n][2];
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
    int           m_var_idx;
    Vector<BCRec> m_domain_bcs_type_bcr;
    TimeInterpolatedData& m_data;
    Array<Array<Real, AMREX_SPACEDIM*2>,AMREX_SPACEDIM+NVAR> m_bc_extdir_vals;
};

#endif
