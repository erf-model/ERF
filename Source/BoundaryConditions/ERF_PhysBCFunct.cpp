#include "AMReX_PhysBCFunct.H"
#include "IndexDefines.H"
#include <ERF_PhysBCFunct.H>

using namespace amrex;

/*
 * Impose physical boundary conditions at domain boundaries
 *
 * @param[inout] mf  MultiFab of cell-centered quantities to be filled
 * @param[in] icomps starting component
 * @param[in] ncomps number of components
 * @param[in] nghost number of ghost cells to be filled
 */

void ERFPhysBCFunct_cons::operator() (MultiFab& mf, int icomp, int ncomp,
                                      IntVect const& nghost, const Real /*time*/, int bccomp)
{
    BL_PROFILE("ERFPhysBCFunct_cons::()");

    if (m_geom.isAllPeriodic()) return;

    const auto& domain = m_geom.Domain();
    const auto dxInv   = m_geom.InvCellSizeArray();

    // Create a grown domain box containing valid + periodic cells
    Box gdomain  = domain;
    for (int i = 0; i < AMREX_SPACEDIM; ++i) {
        if (m_geom.isPeriodic(i)) {
            gdomain.grow(i, nghost[i]);
        }
    }

#ifdef AMREX_USE_OMP
#pragma omp parallel if (Gpu::notInLaunchRegion())
#endif
    {
        for (MFIter mfi(mf,false); mfi.isValid(); ++mfi)
        {
            //
            // This is the box we pass to the different routines
            // NOTE -- this is the full grid box NOT the tile box
            //
            Box bx  = mfi.validbox();

            //
            // These are the boxes we use to test on relative to the domain
            //
            Box cbx1 = bx; cbx1.grow(IntVect(nghost[0],nghost[1],0));
            Box cbx2 = bx; cbx2.grow(nghost);

            Array4<const Real> z_nd_arr;

            if (m_z_phys_nd)
            {
                z_nd_arr = m_z_phys_nd->const_array(mfi);
            }

            if (!gdomain.contains(cbx2))
            {
                const Array4<Real> cons_arr = mf.array(mfi);;

                if (!m_use_real_bcs)
                {
                    impose_lateral_cons_bcs(cons_arr,cbx1,domain,icomp,ncomp,bccomp);
                }

                impose_vertical_cons_bcs(cons_arr,cbx2,domain,z_nd_arr,dxInv,icomp,ncomp,bccomp);
            }

        } // MFIter
    } // OpenMP
} // operator()

void ERFPhysBCFunct_u::operator() (MultiFab& mf, int /*icomp*/, int /*ncomp*/,
                                   IntVect const& nghost, const Real time, int bccomp)
{
    BL_PROFILE("ERFPhysBCFunct_u::()");

    if (m_geom.isAllPeriodic()) return;

    const auto& domain = m_geom.Domain();
    const auto dxInv   = m_geom.InvCellSizeArray();

    Box gdomainx = surroundingNodes(domain,0);
    for (int i = 0; i < AMREX_SPACEDIM; ++i) {
        if (m_geom.isPeriodic(i)) {
            gdomainx.grow(i, nghost[i]);
        }
    }

#ifdef AMREX_USE_OMP
#pragma omp parallel if (Gpu::notInLaunchRegion())
#endif
    {
        for (MFIter mfi(mf,false); mfi.isValid(); ++mfi)
        {
            //
            // This is the box we pass to the different routines
            // NOTE -- this is the full grid NOT the tile box
            //
            Box bx  = mfi.validbox();

            //
            // These are the boxes we use to test on relative to the domain
            //
            Box xbx1 = surroundingNodes(bx,0); xbx1.grow(IntVect(nghost[0],nghost[1],0));
            Box xbx2 = surroundingNodes(bx,0); xbx2.grow(nghost); 

            Array4<const Real> z_nd_arr;

            if (m_z_phys_nd)
            {
                z_nd_arr = m_z_phys_nd->const_array(mfi);
            }

            if (!gdomainx.contains(xbx2))
            {
                const Array4<Real> velx_arr = mf.array(mfi);;

                if (!m_use_real_bcs)
                {
                    if (!gdomainx.contains(xbx1))
                    {
                        impose_lateral_xvel_bcs(velx_arr,xbx1,domain,bccomp);
                    }
                }

                impose_vertical_xvel_bcs(velx_arr,xbx2,domain,z_nd_arr,dxInv,bccomp,time);
            }

        } // MFIter
    } // OpenMP
} // operator()

void ERFPhysBCFunct_v::operator() (MultiFab& mf, int /*icomp*/, int /*ncomp*/,
                                   IntVect const& nghost, const Real /*time*/, int bccomp)
{
    BL_PROFILE("ERFPhysBCFunct_v::()");

    if (m_geom.isAllPeriodic()) return;

    const auto& domain = m_geom.Domain();
    const auto dxInv   = m_geom.InvCellSizeArray();

    Box gdomainy = surroundingNodes(domain,1);
    for (int i = 0; i < AMREX_SPACEDIM; ++i) {
        if (m_geom.isPeriodic(i)) {
            gdomainy.grow(i, nghost[i]);
        }
    }

#ifdef AMREX_USE_OMP
#pragma omp parallel if (Gpu::notInLaunchRegion())
#endif
    {
        for (MFIter mfi(mf,false); mfi.isValid(); ++mfi)
        {
            //
            // This is the box we pass to the different routines
            // NOTE -- this is the full grid NOT the tile box
            //
            Box bx  = mfi.validbox();

            //
            // These are the boxes we use to test on relative to the domain
            //
            Box ybx1 = surroundingNodes(bx,1); ybx1.grow(IntVect(nghost[0],nghost[1],0));
            Box ybx2 = surroundingNodes(bx,1); ybx2.grow(nghost); 

            Array4<const Real> z_nd_arr;

            if (m_z_phys_nd)
            {
                z_nd_arr = m_z_phys_nd->const_array(mfi);
            }

            if (!gdomainy.contains(ybx2))
            {
                const Array4<Real> vely_arr = mf.array(mfi);;

                if (!m_use_real_bcs)
                {
                    impose_lateral_yvel_bcs(vely_arr,ybx1,domain,bccomp);
                }

                impose_vertical_yvel_bcs(vely_arr,ybx2,domain,z_nd_arr,dxInv,bccomp);
            }

        } // MFIter
    } // OpenMP
} // operator()

void ERFPhysBCFunct_w::operator() (MultiFab& mf, MultiFab& xvel, MultiFab& yvel,
                                   IntVect const& nghost, const Real /*time*/,
                                   const int bccomp_u, const int bccomp_v, const int bccomp_w)
{
    BL_PROFILE("ERFPhysBCFunct_w::()");

    if (m_geom.isAllPeriodic()) return;

    const auto& domain = m_geom.Domain();
    const auto dxInv   = m_geom.InvCellSizeArray();

    Box gdomainz = surroundingNodes(domain,2);
    for (int i = 0; i < AMREX_SPACEDIM; ++i) {
        if (m_geom.isPeriodic(i)) {
            gdomainz.grow(i, nghost[i]);
        }
    }
    // We want to make sure we impose the z-vels at k=0  if the box includes k=0
    if (gdomainz.smallEnd(2) == 0) gdomainz.setSmall(2,1);

#ifdef AMREX_USE_OMP
#pragma omp parallel if (Gpu::notInLaunchRegion())
#endif
    {
        for (MFIter mfi(mf,false); mfi.isValid(); ++mfi)
        {
            //
            // This is the box we pass to the different routines
            // NOTE -- this is the full grid NOT the tile box
            //
            Box bx  = mfi.validbox();

            //
            // These are the boxes we use to test on relative to the domain
            //
            Box zbx = surroundingNodes(bx,2); zbx.grow(0,nghost[0]);
                                              zbx.grow(1,nghost[1]);
            Array4<const Real> z_nd_arr;

            if (m_z_phys_nd)
            {
                z_nd_arr = m_z_phys_nd->const_array(mfi);
            }

            Array4<const Real> const& velx_arr = xvel.const_array(mfi);;
            Array4<const Real> const& vely_arr = yvel.const_array(mfi);;
            Array4<      Real> const& velz_arr = mf.array(mfi);;

            if (!m_use_real_bcs)
            {
                if (!gdomainz.contains(zbx))
                {
                    impose_lateral_zvel_bcs(velz_arr,velx_arr,vely_arr,zbx,domain,z_nd_arr,dxInv,bccomp_w);
                }
            } // m_use_real_bcs

            if (!gdomainz.contains(zbx)) {
                impose_vertical_zvel_bcs(velz_arr,velx_arr,vely_arr,zbx,domain,z_nd_arr,dxInv,
                                         bccomp_u, bccomp_v, bccomp_w, m_terrain_type);
            }
        } // MFIter
    } // OpenMP
} // operator()

void ERFPhysBCFunct_w_no_terrain::operator() (MultiFab& mf, int /*icomp*/, int /*ncomp*/,
                                              IntVect const& nghost, const Real /*time*/, int bccomp)
{
    BL_PROFILE("ERFPhysBCFunct_w::()");

    if (m_geom.isAllPeriodic()) return;

    const auto& domain = m_geom.Domain();

    Box gdomainz = surroundingNodes(domain,2);
    for (int i = 0; i < AMREX_SPACEDIM; ++i) {
        if (m_geom.isPeriodic(i)) {
            gdomainz.grow(i, nghost[i]);
        }
    }

#ifdef AMREX_USE_OMP
#pragma omp parallel if (Gpu::notInLaunchRegion())
#endif
    {
        for (MFIter mfi(mf,false); mfi.isValid(); ++mfi)
        {
            //
            // This is the box we pass to the different routines
            // NOTE -- this is the full grid NOT the tile box
            //
            Box bx  = mfi.validbox();

            //
            // These are the boxes we use to test on relative to the domain
            //
            Box zbx = surroundingNodes(bx,2); zbx.grow(0,nghost[0]);
                                              zbx.grow(1,nghost[1]);

            if (!m_use_real_bcs)
            {
                Array4<      Real> const& velz_arr = mf.array(mfi);;
                if (!gdomainz.contains(zbx))
                {
                    impose_lateral_zvel_bcs(velz_arr,zbx,domain,bccomp);
                }
            } // m_use_real_bcs

            const Array4<      Real> velz_arr = mf.array(mfi);;
            if (!gdomainz.contains(zbx)) {
                impose_vertical_zvel_bcs(velz_arr,zbx,domain,bccomp);
            }

        } // MFIter
    } // OpenMP
} // operator()
