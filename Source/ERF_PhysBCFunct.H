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
#include <AMReX_Orientation.H>

#include <TimeInterpolatedData.H>
#include <IndexDefines.H>
#include <DataStruct.H>
#include <EddyViscosity.H>
#include <ABLMost.H>

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
    ERFPhysBCFunct (const Geometry& geom, const Vector<BCRec>& domain_bcs_type,
                    const Gpu::DeviceVector<BCRec>& domain_bcs_type_d,
                    const int& var_idx, TimeInterpolatedData& data,
                    Array<Array<Real,AMREX_SPACEDIM*2>,AMREX_SPACEDIM+NVAR> bc_extdir_vals,
                    const SolverChoice& solverChoice,
                    ABLMost& most)
        : m_geom(geom), m_domain_bcs_type(domain_bcs_type), m_domain_bcs_type_d(domain_bcs_type_d),
          m_var_idx(var_idx), m_data(data), m_bc_extdir_vals(bc_extdir_vals), m_most(most)
        {
            m_has_most_bcs = false;
            for (const auto& bcr : domain_bcs_type) {
                for (OrientationIter oiter; oiter && !m_has_most_bcs; ++oiter) {
                    int side = (int) oiter();
                    if (bcr.data()[side] == static_cast<int>(ERFBCType::MOST)) {
                        m_has_most_bcs = true;
                    }
                }
                if (m_has_most_bcs) break;
            }

            if (m_has_most_bcs) {
                const auto& cons_mf = m_data.get_var(Vars::cons);
                m_viscosity.define(cons_mf.boxArray(), cons_mf.DistributionMap(), 1, cons_mf.nGrowVect());
                ComputeTurbulentViscosity(m_data.get_var(Vars::xvel),
                                          m_data.get_var(Vars::yvel),
                                          m_data.get_var(Vars::zvel),
                                          cons_mf, m_viscosity, geom, solverChoice,
                                          domain_bcs_type_d);
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

                // make Array4's for our data
                const auto cons_arr = m_data.get_var(Vars::cons)[mfi].array();
                const auto velx_arr = m_data.get_var(Vars::xvel)[mfi].array();
                const auto vely_arr = m_data.get_var(Vars::yvel)[mfi].array();
                const auto velz_arr = m_data.get_var(Vars::zvel)[mfi].array();
                const auto eta_arr = m_has_most_bcs ? m_viscosity[mfi].array() : Array4<Real>();

                //! if there are cells not in the valid + periodic grown box
                //! we need to fill them here
                //!
                if (!gdomain.contains(bx))
                {
                    //! Based on BCRec for the domain, we need to make BCRec for this Box
                    // bccomp is used as starting index for m_domain_bcs_type
                    //      0 is used as starting index for bcrs
                    amrex::setBC(bx, domain, bccomp, 0, ncomp, m_domain_bcs_type, bcrs);

                    // in the case of MOST BCs, use reflect_odd in special cases depending
                    // on which variables we are working with and in which directions we fill BCs.
                    for (int i = 0; i < bcrs.size(); i++) {
                        for (int iface = 0; iface <= 2*AMREX_SPACEDIM; iface++) {
                            if (bcrs[i].data()[iface] == static_cast<int>(ERFBCType::MOST)) {
                                const int idir = iface % AMREX_SPACEDIM;
                                if (((m_var_idx == Vars::xvel || m_var_idx == Vars::xmom) && idir == 0) ||
                                    ((m_var_idx == Vars::yvel || m_var_idx == Vars::ymom) && idir == 1) ||
                                    ((m_var_idx == Vars::zvel || m_var_idx == Vars::zmom) && idir == 2)) {
                                    bcrs[i].data()[iface] == ERFBCType::reflect_odd;
                                }
                            }
                        }
                    }

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

                    if (m_var_idx == Vars::xvel || m_var_idx == Vars::xmom ||
                        m_var_idx == Vars::yvel || m_var_idx == Vars::ymom ||
                        m_var_idx == Vars::zvel || m_var_idx == Vars::zmom) {
                        AMREX_ALWAYS_ASSERT(ncomp == 1 && icomp == 0);
                    } else {
                        AMREX_ALWAYS_ASSERT(icomp+ncomp <= NVAR);
                    }

                    amrex::GpuArray<amrex::GpuArray<amrex::Real, AMREX_SPACEDIM*2>,
                                                                 AMREX_SPACEDIM+NVAR> l_bc_extdir_vals_d;

                    for (int i = 0; i < ncomp; i++)
                        for (int ori = 0; ori < 2*AMREX_SPACEDIM; ori++)
                            l_bc_extdir_vals_d[i][ori] = m_bc_extdir_vals[bccomp+i][ori];

                    const amrex::BCRec* bc_ptr = bcrs_d.data();

                    int zlo = m_geom.Domain().smallEnd(2);
                    int zhi = m_geom.Domain().bigEnd(2);

                    // Fill here all the "generic" ext_dir bc's
                    ParallelFor(bx, ncomp, [=,m_var_idx=m_var_idx,m_most=m_most] AMREX_GPU_DEVICE (int i, int j, int k, int n)
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

                        // check m_var_idx to distinguish between Vars::Dvel and Vars::Dmom (in Legacy this was controlled by the is_derived flag)
                        // check the bc_ptr for ERFBCType::MOST
                        if (bccomp == BCVars::xvel_bc)
                        {
                            if (i == dom_lo.x && bc_ptr[n].lo(0) == ERFBCType::ext_dir)
                                dest_array(i,j,k,icomp+n) = l_bc_extdir_vals_d[n][0];
                        }
                        if (bccomp == BCVars::yvel_bc)
                        {
                            if (j == dom_lo.y && bc_ptr[n].lo(1) == ERFBCType::ext_dir)
                                dest_array(i,j,k,icomp+n) = l_bc_extdir_vals_d[n][1];
                        }
                        if (bccomp == BCVars::zvel_bc)
                        {
                            if (k == dom_lo.z && bc_ptr[n].lo(2) == ERFBCType::ext_dir)
                                dest_array(i,j,k,icomp+n) = l_bc_extdir_vals_d[n][2];
                        }

                        bool var_is_derived = false;
                        if (m_var_idx == Vars::xvel || m_var_idx == Vars::yvel || m_var_idx == Vars::zvel) {
                            var_is_derived = true;
                        }

                        if (bc_ptr[n].lo(2) == ERFBCType::MOST) {
                            if (k < zlo) {
                            Real velx, vely, rho, theta, eta;
                            int ix, jx, iy, jy, ie, je;
                            if (m_var_idx == Vars::cons && n == Cons::RhoTheta) {
                                /**
                                * for zone centered variables, the face centered velocity need to be interpolated to
                                * zone center.
                                * NOTE: the number of ghost zone for state variables are different from face centered
                                *       variables in the new version.
                                */
                                ix = i < lbound(velx_arr).x ? lbound(velx_arr).x : i;
                                jx = j < lbound(velx_arr).y ? lbound(velx_arr).y : j;
                                ix = i > ubound(velx_arr).x ? ubound(velx_arr).x-1 : i;
                                jx = j > ubound(velx_arr).y ? ubound(velx_arr).y : j;

                                iy = i < lbound(vely_arr).x ? lbound(vely_arr).x : i;
                                jy = j < lbound(vely_arr).y ? lbound(vely_arr).y : j;
                                iy = i > ubound(vely_arr).x ? ubound(vely_arr).x : i;
                                jy = j > ubound(vely_arr).y ? ubound(vely_arr).y-1 : j;

                                ie = i < lbound(eta_arr).x ? lbound(eta_arr).x : i;
                                je = j < lbound(eta_arr).y ? lbound(eta_arr).y : j;
                                ie = i > ubound(eta_arr).x ? ubound(eta_arr).x : i;
                                je = j > ubound(eta_arr).y ? ubound(eta_arr).y : j;

                                velx  = 0.5*(velx_arr(ix,jx,zlo)+velx_arr(ix+1,jx,zlo));
                                vely  = 0.5*(vely_arr(iy,jy,zlo)+velx_arr(iy,jy+1,zlo));
                                rho   = cons_arr(i,j,zlo,Rho_comp);
                                theta = cons_arr(i,j,zlo,RhoTheta_comp)/rho;
                                eta   = eta_arr(ie,je,zlo);

                                Real vmag    = sqrt(velx*velx+vely*vely);
                                Real num1    = (theta-m_most.theta_mean)*m_most.vmag_mean;
                                Real num2    = (m_most.theta_mean-m_most.surf_temp)*vmag;
                                Real motheta = (num1+num2)*m_most.utau*m_most.kappa/m_most.phi_h();

                                //printf("i=%d,j=%d,k=%d,theta=%13.6e,rho=%13.6e,eta=%13.6e,tau=%13.6e,val=%13.6e\n",i,j,k,motheta,rho,(eta_arr)(60,10,1),m_most.utau,motheta*rho/eta);
                                if (!var_is_derived) {
                                    dest_array(i,j,k,icomp+n) = rho*(m_most.surf_temp + motheta*rho/eta);
                                } else {
                                    dest_array(i,j,k,icomp+n) = m_most.surf_temp + motheta/eta;
                                }

                            } else if (m_var_idx == Vars::xvel || m_var_idx == Vars::xmom) { //for velx
                                iy = i < lbound(vely_arr).x ? lbound(vely_arr).x : i;
                                jy = j < lbound(vely_arr).y ? lbound(vely_arr).y : j;
                                iy = i > ubound(vely_arr).x ? ubound(vely_arr).x : i;
                                jy = j > ubound(vely_arr).y ? ubound(vely_arr).y-1 : j;

                                ie = i < lbound(eta_arr).x ? lbound(eta_arr).x : i;
                                je = j < lbound(eta_arr).y ? lbound(eta_arr).y : j;
                                ie = i > ubound(eta_arr).x ? ubound(eta_arr).x-1 : i;
                                je = j > ubound(eta_arr).y ? ubound(eta_arr).y : j;

                                velx  = velx_arr(i,j,zlo);
                                vely  = 0.5*(vely_arr(iy,jy,zlo)+vely_arr(iy,jy+1,zlo));
                                rho   = 0.5*(cons_arr(i,j,zlo,Rho_comp)+cons_arr(i+1,j,zlo,Rho_comp));
                                eta   = 0.5*(eta_arr(ie,je,zlo)+eta_arr(ie+1,je,zlo));
                                Real vmag  = sqrt(velx*velx+vely*vely);
                                Real vgx   = ((velx-m_most.vel_mean[0])*m_most.vmag_mean + vmag*m_most.vel_mean[0])/(m_most.vmag_mean*m_most.vmag_mean)*m_most.utau*m_most.utau;
                                Real vgy   = ((vely-m_most.vel_mean[1])*m_most.vmag_mean + vmag*m_most.vel_mean[1])/(m_most.vmag_mean*m_most.vmag_mean)*m_most.utau*m_most.utau;

                                //printf("i=%d,j=%d,k=%d,vx=%13.6e,rho=%13.6e,eta=%13.6e,tau=%13.6e,val=%13.6e\n",i,j,k,vgx,rho,eta,m_most.utau,vgx*rho/eta);
                                if (!var_is_derived) {
                                    dest_array(i,j,k,icomp+n) = dest_array(i,j,zlo,icomp+n) - vgx*rho/eta;
                                } else {
                                    dest_array(i,j,k,icomp+n) = dest_array(i,j,zlo,icomp+n) - vgx/eta;
                                }

                            } else if (m_var_idx == Vars::yvel || m_var_idx == Vars::ymom) { //for vely
                                ix = i < lbound(velx_arr).x ? lbound(velx_arr).x : i;
                                jx = j < lbound(velx_arr).y ? lbound(velx_arr).y : j;
                                ix = i > ubound(velx_arr).x ? ubound(velx_arr).x : i;
                                jx = j > ubound(velx_arr).y ? ubound(velx_arr).y-1 : j;

                                ie = i < lbound(eta_arr).x ? lbound(eta_arr).x : i;
                                je = j < lbound(eta_arr).y ? lbound(eta_arr).y : j;
                                ie = i > ubound(eta_arr).x ? ubound(eta_arr).x : i;
                                je = j > ubound(eta_arr).y ? ubound(eta_arr).y-1 : j;

                                velx  = 0.5*(velx_arr(ix,jx,zlo)+velx_arr(ix,jx+1,zlo));
                                vely  = vely_arr(i,j,zlo);
                                rho   = 0.5*(cons_arr(i,j,zlo,Rho_comp)+cons_arr(i,j+1,zlo,Rho_comp));
                                eta   = 0.5*(eta_arr(ie,je,zlo)+eta_arr(ie,je+1,zlo));
                                Real vmag  = sqrt(velx*velx+vely*vely);
                                Real vgx   = ((velx-m_most.vel_mean[0])*m_most.vmag_mean + vmag*m_most.vel_mean[0])/(m_most.vmag_mean*m_most.vmag_mean)*m_most.utau*m_most.utau;
                                Real vgy   = ((vely-m_most.vel_mean[1])*m_most.vmag_mean + vmag*m_most.vel_mean[1])/(m_most.vmag_mean*m_most.vmag_mean)*m_most.utau*m_most.utau;

                                //printf("i=%d,j=%d,k=%d,vy=%13.6e,rho=%13.6e,eta=%13.6e,tau=%13.6e,val=%13.6e\n",i,j,k,vgy,rho,eta,m_most.utau,vgy*rho/eta);
                                if (!var_is_derived) {
                                    dest_array(i,j,k,icomp+n) = dest_array(i,j,zlo,icomp+n) - vgy*rho/eta;
                                } else {
                                    dest_array(i,j,k,icomp+n) = dest_array(i,j,zlo,icomp+n) - vgy/eta;
                                }
                            } else if (m_var_idx == Vars::zvel || m_var_idx == Vars::zmom) { //for velz
                                dest_array(i,j,k,icomp+n) = (zlo-k+1)*dest_array(i,j,zlo,icomp+n) - (zlo-k)*dest_array(i,j,zlo+1,icomp+n);
                            }
                            }
                        }

                        if (bc_ptr[n].hi(2) == ERFBCType::MOST) {
                            // Not tested yet, FIX ME if you find bugs
                            if (k > zhi) {
                            Real velx, vely, rho, theta, eta;
                            int ix, jx, iy, jy, ie, je;
                            if (m_var_idx == Vars::cons && n == Cons::RhoTheta) {
                                /**
                                * for zone centered variables, the face centered velocity need to be interpolated to
                                * zone center.
                                * NOTE: the number of ghost zone for state variables are different from face centered
                                *       variables in the new version.
                                */
                                ix = i < lbound(velx_arr).x ? lbound(velx_arr).x : i;
                                jx = j < lbound(velx_arr).y ? lbound(velx_arr).y : j;
                                ix = i > ubound(velx_arr).x ? ubound(velx_arr).x-1 : i;
                                jx = j > ubound(velx_arr).y ? ubound(velx_arr).y : j;

                                iy = i < lbound(vely_arr).x ? lbound(vely_arr).x : i;
                                jy = j < lbound(vely_arr).y ? lbound(vely_arr).y : j;
                                iy = i > ubound(vely_arr).x ? ubound(vely_arr).x : i;
                                jy = j > ubound(vely_arr).y ? ubound(vely_arr).y-1 : j;

                                ie = i < lbound(eta_arr).x ? lbound(eta_arr).x : i;
                                je = j < lbound(eta_arr).y ? lbound(eta_arr).y : j;
                                ie = i > ubound(eta_arr).x ? ubound(eta_arr).x : i;
                                je = j > ubound(eta_arr).y ? ubound(eta_arr).y : j;

                                velx  = 0.5*(velx_arr(ix,jx,zhi)+velx_arr(ix+1,jx,zhi));
                                vely  = 0.5*(vely_arr(iy,jy,zhi)+velx_arr(iy,jy+1,zhi));
                                rho   = cons_arr(i,j,zhi,Rho_comp);
                                theta = cons_arr(i,j,zhi,RhoTheta_comp)/rho;
                                eta   = eta_arr(ie,je,zhi);

                                Real vmag    = sqrt(velx*velx+vely*vely);
                                Real num1    = (theta-m_most.theta_mean)*m_most.vmag_mean;
                                Real num2    = (m_most.theta_mean-m_most.surf_temp)*vmag;
                                Real motheta = (num1+num2)*m_most.utau*m_most.kappa/m_most.phi_h();

                                //printf("i=%d,j=%d,k=%d,theta=%13.6e,rho=%13.6e,eta=%13.6e,tau=%13.6e,val=%13.6e\n",i,j,k,motheta,rho,(eta_arr)(60,10,1),m_most.utau,motheta*rho/eta);
                                if (!var_is_derived) {
                                    dest_array(i,j,k,icomp+n) = rho*(dest_array(i,j,zhi,icomp+n) - motheta*rho/eta);
                                } else {
                                    dest_array(i,j,k,icomp+n) = dest_array(i,j,zhi,icomp+n) - motheta/eta;
                                }

                            } else if (m_var_idx == Vars::xvel || m_var_idx == Vars::xmom) { //for velx
                                iy = i < lbound(vely_arr).x ? lbound(vely_arr).x : i;
                                jy = j < lbound(vely_arr).y ? lbound(vely_arr).y : j;
                                iy = i > ubound(vely_arr).x ? ubound(vely_arr).x : i;
                                jy = j > ubound(vely_arr).y ? ubound(vely_arr).y-1 : j;

                                ie = i < lbound(eta_arr).x ? lbound(eta_arr).x : i;
                                je = j < lbound(eta_arr).y ? lbound(eta_arr).y : j;
                                ie = i > ubound(eta_arr).x ? ubound(eta_arr).x-1 : i;
                                je = j > ubound(eta_arr).y ? ubound(eta_arr).y : j;

                                velx  = velx_arr(i,j,zhi);
                                vely  = 0.5*(vely_arr(iy,jy,zhi)+vely_arr(iy,jy+1,zhi));
                                rho   = 0.5*(cons_arr(i,j,zhi,Rho_comp)+cons_arr(i+1,j,zhi,Rho_comp));
                                eta   = 0.5*(eta_arr(ie,je,zhi)+eta_arr(ie+1,je,zhi));
                                Real vmag  = sqrt(velx*velx+vely*vely);
                                Real vgx   = ((velx-m_most.vel_mean[0])*m_most.vmag_mean + vmag*m_most.vel_mean[0])/(m_most.vmag_mean*m_most.vmag_mean)*m_most.utau*m_most.utau;
                                Real vgy   = ((vely-m_most.vel_mean[1])*m_most.vmag_mean + vmag*m_most.vel_mean[1])/(m_most.vmag_mean*m_most.vmag_mean)*m_most.utau*m_most.utau;

                                //printf("i=%d,j=%d,k=%d,vx=%13.6e,rho=%13.6e,eta=%13.6e,tau=%13.6e,val=%13.6e\n",i,j,k,vgx,rho,eta,m_most.utau,vgx*rho/eta);

                                if (!var_is_derived) {
                                    dest_array(i,j,k,icomp+n) = dest_array(i,j,zhi,icomp+n) + vgx*rho/eta;
                                } else {
                                    dest_array(i,j,k,icomp+n) = dest_array(i,j,zhi,icomp+n) + vgx/eta;
                                }

                            } else if (m_var_idx == Vars::yvel || m_var_idx == Vars::ymom) { //for vely
                                ix = i < lbound(velx_arr).x ? lbound(velx_arr).x : i;
                                jx = j < lbound(velx_arr).y ? lbound(velx_arr).y : j;
                                ix = i > ubound(velx_arr).x ? ubound(velx_arr).x : i;
                                jx = j > ubound(velx_arr).y ? ubound(velx_arr).y-1 : j;

                                ie = i < lbound(eta_arr).x ? lbound(eta_arr).x : i;
                                je = j < lbound(eta_arr).y ? lbound(eta_arr).y : j;
                                ie = i > ubound(eta_arr).x ? ubound(eta_arr).x : i;
                                je = j > ubound(eta_arr).y ? ubound(eta_arr).y-1 : j;

                                velx  = 0.5*(velx_arr(ix,jx,zhi)+velx_arr(ix,jx+1,zhi));
                                vely  = vely_arr(i,j,zhi);
                                rho   = 0.5*(cons_arr(i,j,zhi,Rho_comp)+cons_arr(i,j+1,zhi,Rho_comp));
                                eta   = 0.5*(eta_arr(ie,je,zhi)+eta_arr(ie,je+1,zhi));
                                Real vmag  = sqrt(velx*velx+vely*vely);
                                Real vgx   = ((velx-m_most.vel_mean[0])*m_most.vmag_mean + vmag*m_most.vel_mean[0])/(m_most.vmag_mean*m_most.vmag_mean)*m_most.utau*m_most.utau;
                                Real vgy   = ((vely-m_most.vel_mean[1])*m_most.vmag_mean + vmag*m_most.vel_mean[1])/(m_most.vmag_mean*m_most.vmag_mean)*m_most.utau*m_most.utau;

                                //printf("i=%d,j=%d,k=%d,vy=%13.6e,rho=%13.6e,eta=%13.6e,tau=%13.6e,val=%13.6e\n",i,j,k,vgy,rho,eta,m_most.utau,vgy*rho/eta);
                                if (!var_is_derived) {
                                    dest_array(i,j,k,icomp+n) = dest_array(i,j,zhi,icomp+n) + vgy*rho/eta;
                                } else {
                                    dest_array(i,j,k,icomp+n) = dest_array(i,j,zhi,icomp+n) + vgy/eta;
                                }
                            } else if (m_var_idx == Vars::zvel || m_var_idx == Vars::zmom) { //for velz
                                dest_array(i,j,k,icomp+n) = (k-zhi+1)*dest_array(i,j,zhi,icomp+n) - (k-zhi)*dest_array(i,j,zhi-1,icomp+n);
                            }
                            }
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
    Vector<BCRec> m_domain_bcs_type;
    Gpu::DeviceVector<BCRec> m_domain_bcs_type_d;
    int           m_var_idx;
    TimeInterpolatedData& m_data;
    Array<Array<Real, AMREX_SPACEDIM*2>,AMREX_SPACEDIM+NVAR> m_bc_extdir_vals;
    MultiFab m_viscosity;
    bool m_has_most_bcs;
    ABLMost m_most;
};

#endif
