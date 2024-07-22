#ifdef AMREX_USE_EB
#include <AMReX_EBMultiFabUtil.H>
#endif

#include <AMReX_MultiFabUtil.H>
#include <AMReX_ParmParse.H>
#include <AMReX_BC_TYPES.H>

using namespace amrex;

void set_masks (Array< MultiFab*, AMREX_SPACEDIM> a_umac,
                Array<iMultiFab , AMREX_SPACEDIM>& inflow_masks,
                Array<iMultiFab , AMREX_SPACEDIM>& outflow_masks,
                const Box& domain)
{
    for (int dir = 0; dir < AMREX_SPACEDIM; dir++)
    {
        for (MFIter mfi(*a_umac[dir], false); mfi.isValid(); ++mfi)
        {
            auto&  inflow_mask = inflow_masks[dir];
            auto& outflow_mask = outflow_masks[dir];

            auto  mac_vel = a_umac[dir]->array(mfi);
            auto  in_mask = inflow_mask.array(mfi);
            auto out_mask = outflow_mask.array(mfi);

            const Box vbx = mfi.validbox();

            // Low side
            const Box box2d_lo = makeSlab(vbx,dir,vbx.smallEnd(dir));
            ParallelFor(box2d_lo, [=] AMREX_GPU_DEVICE (int i, int j, int k)
            {
                if (mac_vel(i,j,k) >= 0) {
                      in_mask(i,j,k) = 1;
                 } else {
                     out_mask(i,j,k) = 1;
                 }
            });

            // High side
            const Box box2d_hi = makeSlab(vbx,dir,vbx.bigEnd(dir));
            ParallelFor(box2d_hi, [=] AMREX_GPU_DEVICE (int i, int j, int k)
            {
                if (mac_vel(i,j,k) <= 0) {
                     in_mask(i,j,k) = 1;
                } else {
                    out_mask(i,j,k) = 1;
                }
            });

        } // mfi
    } // dir
}


Real
compute_influx_outflux(Array< MultiFab*, AMREX_SPACEDIM> a_umac,
                       const Array<iMultiFab, AMREX_SPACEDIM>& inflow_masks,
                       const Array<iMultiFab, AMREX_SPACEDIM>& outflow_masks,
                       const Real* a_dx)
{
    Real  influx = 0.0;
    Real outflux = 0.0;

    // Loop over the three dimensions
    for (int dir = 0; dir < AMREX_SPACEDIM; dir++)
    {
        // normal face area
        const Real ds = a_dx[(dir+1) % AMREX_SPACEDIM] * a_dx[(dir+2) % AMREX_SPACEDIM];

        const MultiFab*  mac_vel_mf = a_umac[dir];

        auto&  inflow_mask = inflow_masks[dir];
        auto& outflow_mask = outflow_masks[dir];

        auto const&      mac_vel_ma = mac_vel_mf->const_arrays();
        auto const&  inflow_mask_ma = inflow_mask.const_arrays();
        auto const& outflow_mask_ma = outflow_mask.const_arrays();

        Real contrib_in =
            ParReduce(TypeList<ReduceOpSum>{},
                     TypeList<Real>{},
                     *mac_vel_mf, IntVect(0), // zero ghost cells
           [=] AMREX_GPU_DEVICE (int box_no, int i, int j, int k)
               noexcept -> GpuTuple<Real>
           {
               if (inflow_mask_ma[box_no](i,j,k)) {
                   return { std::abs(mac_vel_ma[box_no](i,j,k)) };
               } else {
                   return { 0. };
               }
           });
        influx += ds * contrib_in;

        Real contrib_out =
            ParReduce(TypeList<ReduceOpSum>{},
                     TypeList<Real>{},
                     *mac_vel_mf, IntVect(0), // zero ghost cells
           [=] AMREX_GPU_DEVICE (int box_no, int i, int j, int k)
               noexcept -> GpuTuple<Real>
           {
               if (outflow_mask_ma[box_no](i,j,k)) {
                   return { std::abs(mac_vel_ma[box_no](i,j,k)) };
               } else {
                   return { 0. };
               }
           });
        outflux += ds * contrib_out;
    }

    // amrex::Print() << "INFLUX OUTFLUX " << influx << " " << outflux << std::endl;
    Real alpha = 0.;
    if (outflux > 0.) alpha = influx/outflux;
    return alpha;
}

void correct_outflow (      Array<MultiFab*, AMREX_SPACEDIM>& a_umac,
                      const Array<iMultiFab, AMREX_SPACEDIM>& outflow_masks,
                      const Real alpha)
{
    for (int dir = 0; dir < AMREX_SPACEDIM; dir++)
    {
        for (MFIter mfi(*a_umac[dir], false); mfi.isValid(); ++mfi)
        {
            auto  mac_vel = a_umac[dir]->array(mfi);
            auto out_mask = outflow_masks[dir].array(mfi);

            const Box vbx = mfi.validbox();

            // Low side
            const Box box2d_lo = makeSlab(vbx,dir,vbx.smallEnd(dir));
            ParallelFor(box2d_lo, [=] AMREX_GPU_DEVICE (int i, int j, int k)
            {
                if (out_mask(i,j,k)) {
                    mac_vel(i,j,k) *= alpha;
                }
            });

            // Low side
            const Box box2d_hi = makeSlab(vbx,dir,vbx.bigEnd(dir));
            ParallelFor(box2d_hi, [=] AMREX_GPU_DEVICE (int i, int j, int k)
            {
                if (out_mask(i,j,k)) {
                    mac_vel(i,j,k) *= alpha;
                }
            });
        } // mfi
    } // dir
}

void
enforceInflowOutflowSolvability (Array<MultiFab*, AMREX_SPACEDIM>& a_umac,
                                 const Geometry& geom)
{
    const Box domain = geom.Domain();

    // masks to tag in/out flow at in-out boundaries
    // separate iMultifab for each velocity direction
    Array<iMultiFab, AMREX_SPACEDIM> inflow_masks;
    Array<iMultiFab, AMREX_SPACEDIM> outflow_masks;

    for (int dir = 0; dir < AMREX_SPACEDIM; dir++)
    {
        auto& mac_vel_mf = a_umac[dir];    // normal mac velocity multifab

        inflow_masks[dir].define(mac_vel_mf->boxArray(), mac_vel_mf->DistributionMap(), 1, 0);
        inflow_masks[dir].setVal(0);

        outflow_masks[dir].define(mac_vel_mf->boxArray(), mac_vel_mf->DistributionMap(), 1, 0);
        outflow_masks[dir].setVal(0);
    }

    set_masks(a_umac, inflow_masks, outflow_masks, domain);

    const Real* a_dx = geom.CellSize();
    Real influx = 0.0, outflux = 0.0;

    // Calculate the ratio of influx to outflux needed to correct for solvability
    Real alpha = compute_influx_outflux(a_umac, inflow_masks, outflow_masks, a_dx);

    // Apply correction factor to outflow
    correct_outflow(a_umac, outflow_masks, alpha);

    // Calculate the ratio of influx to outflux needed to correct for solvability
    Real beta = compute_influx_outflux(a_umac, inflow_masks, outflow_masks, a_dx);
}
