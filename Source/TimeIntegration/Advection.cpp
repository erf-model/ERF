#include <TimeIntegration.H>

using namespace amrex;

AMREX_GPU_DEVICE
Real
AdvectionContributionForXMom(const int &i, const int &j, const int &k,
                            const Array4<const Real>& rho_u, const Array4<const Real>& rho_v, const Array4<const Real>& rho_w,
                            const Array4<Real>& u,
                            const GpuArray<Real, AMREX_SPACEDIM>& cellSizeInv,
                            const int& spatial_order)
{
    auto dxInv = cellSizeInv[0], dyInv = cellSizeInv[1], dzInv = cellSizeInv[2];
    Real rho_u_avg, rho_v_avg, rho_w_avg;

    rho_u_avg = 0.5 * (rho_u(i+1, j, k) + rho_u(i, j, k));
    Real centFluxXXNext = rho_u_avg *
                          InterpolateFromCellOrFace(i+1, j, k, u, 0, rho_u_avg, Coord::x, spatial_order);

    rho_u_avg = 0.5 * (rho_u(i-1, j, k) + rho_u(i, j, k));
    Real centFluxXXPrev = rho_u_avg *
                          InterpolateFromCellOrFace(i  , j, k, u, 0, rho_u_avg, Coord::x, spatial_order);

    rho_v_avg = 0.5 * (rho_v(i, j+1, k) + rho_v(i-1, j+1, k));
    Real edgeFluxXYNext = rho_v_avg *
                          InterpolateFromCellOrFace(i, j+1, k, u, 0, rho_v_avg, Coord::y, spatial_order);

    rho_v_avg = 0.5 * (rho_v(i, j  , k) + rho_v(i-1, j  , k));
    Real edgeFluxXYPrev = rho_v_avg *
                          InterpolateFromCellOrFace(i, j  , k, u, 0, rho_v_avg, Coord::y, spatial_order);

    rho_w_avg = 0.5 * (rho_w(i, j, k+1) + rho_w(i-1, j, k+1));
    Real edgeFluxXZNext = rho_w_avg *
                          InterpolateFromCellOrFace(i, j, k+1, u, 0, rho_w_avg, Coord::z, spatial_order);

    rho_w_avg = 0.5 * (rho_w(i, j, k) + rho_w(i-1, j, k));
    Real edgeFluxXZPrev = rho_w_avg *
                          InterpolateFromCellOrFace(i, j, k  , u, 0, rho_w_avg, Coord::z, spatial_order);

    Real advectionContribution = (centFluxXXNext - centFluxXXPrev) * dxInv
                               + (edgeFluxXYNext - edgeFluxXYPrev) * dyInv
                               + (edgeFluxXZNext - edgeFluxXZPrev) * dzInv;

    return advectionContribution;
}

AMREX_GPU_DEVICE
Real
AdvectionContributionForYMom(const int &i, const int &j, const int &k,
                            const Array4<const Real>& rho_u, const Array4<const Real>& rho_v, const Array4<const Real>& rho_w,
                            const Array4<Real>& v,
                            const GpuArray<Real, AMREX_SPACEDIM>& cellSizeInv,
                            const int& spatial_order)
{
    auto dxInv = cellSizeInv[0], dyInv = cellSizeInv[1], dzInv = cellSizeInv[2];
    Real rho_u_avg, rho_v_avg, rho_w_avg;

    rho_u_avg = 0.5*(rho_u(i+1, j, k) + rho_u(i+1, j-1, k));
    Real edgeFluxYXNext = rho_u_avg *
                          InterpolateFromCellOrFace(i+1, j, k, v, 0, rho_u_avg, Coord::x, spatial_order);

    rho_u_avg = 0.5*(rho_u(i  , j, k) + rho_u(i  , j-1, k));
    Real edgeFluxYXPrev = rho_u_avg *
                          InterpolateFromCellOrFace(i  , j, k, v, 0, rho_u_avg, Coord::x, spatial_order);

    rho_v_avg = 0.5*(rho_v(i, j, k) + rho_v(i, j+1, k));
    Real centFluxYYNext = rho_v_avg *
                          InterpolateFromCellOrFace(i, j+1, k, v, 0, rho_v_avg, Coord::y, spatial_order);

    rho_v_avg = 0.5*(rho_v(i, j, k) + rho_v(i, j-1, k));
    Real centFluxYYPrev = rho_v_avg *
                          InterpolateFromCellOrFace(i, j  , k, v, 0, rho_v_avg, Coord::y, spatial_order);

    rho_w_avg = 0.5*(rho_w(i, j, k+1) + rho_w(i, j-1, k+1));
    Real edgeFluxYZNext = rho_w_avg *
                          InterpolateFromCellOrFace(i, j, k+1, v, 0, rho_w_avg, Coord::z, spatial_order);

    rho_w_avg = 0.5*(rho_w(i, j, k) + rho_w(i, j-1, k));
    Real edgeFluxYZPrev = rho_w_avg *
                          InterpolateFromCellOrFace(i, j, k  , v, 0, rho_w_avg, Coord::z, spatial_order);

    Real advectionContribution = (edgeFluxYXNext - edgeFluxYXPrev) * dxInv
                               + (centFluxYYNext - centFluxYYPrev) * dyInv
                               + (edgeFluxYZNext - edgeFluxYZPrev) * dzInv;

    return advectionContribution;
}

AMREX_GPU_DEVICE
Real
AdvectionContributionForZMom(const int &i, const int &j, const int &k,
                            const Array4<const Real>& rho_u, const Array4<const Real>& rho_v, const Array4<const Real>& rho_w,
                            const Array4<Real>& w,
                            const GpuArray<Real, AMREX_SPACEDIM>& cellSizeInv,
                            const int& spatial_order)
{
    auto dxInv = cellSizeInv[0], dyInv = cellSizeInv[1], dzInv = cellSizeInv[2];
    Real rho_u_avg, rho_v_avg, rho_w_avg;

    rho_u_avg = 0.5*(rho_u(i+1, j, k) + rho_u(i+1, j, k-1));
    Real edgeFluxZXNext = rho_u_avg *
                          InterpolateFromCellOrFace(i+1, j, k, w, 0, rho_u_avg, Coord::x, spatial_order);

    rho_u_avg = 0.5*(rho_u(i  , j, k) + rho_u(i  , j, k-1));
    Real edgeFluxZXPrev = rho_u_avg *
                          InterpolateFromCellOrFace(i  , j, k, w, 0, rho_u_avg, Coord::x, spatial_order);

    rho_v_avg = 0.5*(rho_v(i, j+1, k) + rho_v(i, j+1, k-1));
    Real edgeFluxZYNext = rho_v_avg *
                          InterpolateFromCellOrFace(i, j+1, k, w, 0, rho_v_avg, Coord::y, spatial_order);

    rho_v_avg = 0.5*(rho_v(i, j  , k) + rho_v(i, j  , k-1));
    Real edgeFluxZYPrev = rho_v_avg *
                          InterpolateFromCellOrFace(i, j  , k, w, 0, rho_v_avg, Coord::y, spatial_order);

    rho_w_avg = 0.5*(rho_w(i, j  , k+1) + rho_w(i, j, k));
    Real centFluxZZNext = rho_w_avg *
                          InterpolateFromCellOrFace(i, j, k+1, w, 0, rho_w_avg, Coord::z, spatial_order);

    rho_w_avg = 0.5*(rho_w(i, j  , k-1) + rho_w(i, j, k));
    Real centFluxZZPrev = rho_w_avg *
                          InterpolateFromCellOrFace(i, j, k  , w, 0, rho_w_avg, Coord::z, spatial_order);

    Real advectionContribution = (edgeFluxZXNext - edgeFluxZXPrev) * dxInv
                               + (edgeFluxZYNext - edgeFluxZYPrev) * dyInv
                               + (centFluxZZNext - centFluxZZPrev) * dzInv;

    return advectionContribution;
}

AMREX_GPU_DEVICE
Real
AdvectionContributionForState(const int &i, const int &j, const int &k,
                              const Array4<const Real>& rho_u, const Array4<const Real>& rho_v, const Array4<const Real>& rho_w,
                              const Array4<const Real>& cell_data, const int &qty_index,
                              const Array4<Real>& xflux, const Array4<Real>& yflux, const Array4<Real>& zflux,
                              const GpuArray<Real, AMREX_SPACEDIM>& cellSizeInv,
                              const int &spatial_order) {

    auto dxInv = cellSizeInv[0], dyInv = cellSizeInv[1], dzInv = cellSizeInv[2];

    Real advectionContribution;

    if (qty_index == Rho_comp)
    {
        xflux(i+1,j,k,qty_index) = rho_u(i+1,j,k);
        xflux(i  ,j,k,qty_index) = rho_u(i  ,j,k);
        yflux(i,j+1,k,qty_index) = rho_v(i,j+1,k);
        yflux(i,j  ,k,qty_index) = rho_v(i,j  ,k);
        zflux(i,j,k+1,qty_index) = rho_w(i,j,k+1);
        zflux(i,j,k  ,qty_index) = rho_w(i,j,k  );

    } else {

        Real uadv_hi = rho_u(i+1,j,k);
        Real uadv_lo = rho_u(i  ,j,k);
        Real vadv_hi = rho_v(i,j+1,k);
        Real vadv_lo = rho_v(i,j  ,k);
        Real wadv_hi = rho_w(i,j,k+1);
        Real wadv_lo = rho_w(i,j,k  );
        xflux(i+1,j,k,qty_index) = rho_u(i+1,j,k) *
            InterpolateFromCellOrFace(i+1, j, k, cell_data, qty_index, uadv_hi, Coord::x, spatial_order) /
            InterpolateFromCellOrFace(i+1, j, k, cell_data, Rho_comp , uadv_hi, Coord::x, spatial_order);

        xflux(i  ,j,k,qty_index) = rho_u(i  ,j,k) *
            InterpolateFromCellOrFace(i  , j, k, cell_data, qty_index, uadv_lo, Coord::x, spatial_order) /
            InterpolateFromCellOrFace(i  , j, k, cell_data, Rho_comp , uadv_lo, Coord::x, spatial_order);

        yflux(i,j+1,k,qty_index) = rho_v(i,j+1,k) *
            InterpolateFromCellOrFace(i, j+1, k, cell_data, qty_index, vadv_hi, Coord::y, spatial_order) /
            InterpolateFromCellOrFace(i, j+1, k, cell_data, Rho_comp , vadv_hi, Coord::y, spatial_order);
        yflux(i,j  ,k,qty_index) = rho_v(i,j  ,k) *
            InterpolateFromCellOrFace(i, j  , k, cell_data, qty_index, vadv_lo, Coord::y, spatial_order) /
            InterpolateFromCellOrFace(i, j  , k, cell_data, Rho_comp , vadv_lo, Coord::y, spatial_order);

        zflux(i,j,k+1,qty_index) = rho_w(i,j,k+1)*
            InterpolateFromCellOrFace(i, j, k+1, cell_data, qty_index, wadv_hi, Coord::z, spatial_order) /
            InterpolateFromCellOrFace(i, j, k+1, cell_data, Rho_comp , wadv_hi, Coord::z, spatial_order);
        zflux(i,j,k  ,qty_index) = rho_w(i,j,k  )*
            InterpolateFromCellOrFace(i, j, k  , cell_data, qty_index, wadv_lo, Coord::z, spatial_order) /
            InterpolateFromCellOrFace(i, j, k  , cell_data, Rho_comp , wadv_lo, Coord::z, spatial_order);
    }

    advectionContribution = (xflux(i+1,j,k,qty_index) - xflux(i  ,j,k,qty_index)) * dxInv
                          + (yflux(i,j+1,k,qty_index) - yflux(i,j  ,k,qty_index)) * dyInv
                          + (zflux(i,j,k+1,qty_index) - zflux(i,j,k  ,qty_index)) * dzInv;

    return advectionContribution;
}
