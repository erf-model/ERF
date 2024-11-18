#ifndef ERF_TERRAINPOISSON_H_
#define ERF_TERRAINPOISSON_H_

#include <AMReX_Config.H>
#include <AMReX_MLCellABecLap.H>

#include "ERF_TerrainPoisson_3D_K.H"

namespace amrex {

// del dot grad phi

template <typename MF>
class MLTerrainPoissonT
    : public MLCellABecLapT<MF>
{
public:

    using MLCellABecLapT<MF>::getFluxes;

    using FAB = typename MF::fab_type;
    using RT  = typename MF::value_type;

    using BCType = LinOpBCType;
    using Location  = typename MLLinOpT<MF>::Location;

    MLTerrainPoissonT () = default;
    MLTerrainPoissonT (const Vector<Geometry>& a_geom,
                const Vector<BoxArray>& a_grids,
                const Vector<DistributionMapping>& a_dmap,
                const LPInfo& a_info = LPInfo(),
                const Vector<FabFactory<FAB> const*>& a_factory = {});
    ~MLTerrainPoissonT () override;

    MLTerrainPoissonT (const MLTerrainPoissonT<MF>&) = delete;
    MLTerrainPoissonT (MLTerrainPoissonT<MF>&&) = delete;
    MLTerrainPoissonT<MF>& operator= (const MLTerrainPoissonT<MF>&) = delete;
    MLTerrainPoissonT<MF>& operator= (MLTerrainPoissonT<MF>&&) = delete;

    void define (const Vector<Geometry>& a_geom,
                 const Vector<BoxArray>& a_grids,
                 const Vector<DistributionMapping>& a_dmap,
                 const LPInfo& a_info = LPInfo(),
                 const Vector<FabFactory<FAB> const*>& a_factory = {});

    /**
     * Define the height coordinate
     *
     * \param [in] amrlev The level of the multifab for the solver, with
     *                    \p amrlev = 0 always being the lowest level in the
     *                    AMR hierarchy represented in the solve.
     * \param [in] zphys  Multifab of z_phys_nd
     */
    template <typename AMF,
              std::enable_if_t<IsFabArray<AMF>::value &&
                               std::is_convertible_v<typename AMF::value_type,
                                                     typename  MF::value_type>,
                               int> = 0>
    void setZPhysNd (int amrlev, const AMF& zphys);

    void prepareForSolve () final;
    [[nodiscard]] bool isSingular (int amrlev) const final { return m_is_singular[amrlev]; }
    [[nodiscard]] bool isBottomSingular () const final { return m_is_singular[0]; }

    void Fapply (int amrlev, int mglev, MF& out, const MF& in) const final;
    void Fsmooth (int amrlev, int mglev, MF& sol, const MF& rhs, int redblack) const final;
    void FFlux (int amrlev, const MFIter& mfi,
                const Array<FAB*,AMREX_SPACEDIM>& flux,
                const FAB& sol, Location loc, int face_only=0) const final;

    void getFluxes (const Vector<Array<MF*,AMREX_SPACEDIM> >& a_flux,
                    const Vector<MF*>& a_sol); // const final;

    Vector<MF> m_a_zphys_nd;

    void normalize (int amrlev, int mglev, MF& mf) const final;

    [[nodiscard]] RT getAScalar () const final { return RT(0.0); }
    [[nodiscard]] RT getBScalar () const final { return RT(-1.0); }
    [[nodiscard]] MF const* getACoeffs (int /*amrlev*/, int /*mglev*/) const final { return nullptr; }
    [[nodiscard]] Array<MF const*,AMREX_SPACEDIM> getBCoeffs (int /*amrlev*/, int /*mglev*/) const final
        { return {{ AMREX_D_DECL(nullptr,nullptr,nullptr)}}; }

private:

    Vector<int> m_is_singular;
};

template <typename MF>
MLTerrainPoissonT<MF>::MLTerrainPoissonT (const Vector<Geometry>& a_geom,
                                      const Vector<BoxArray>& a_grids,
                                      const Vector<DistributionMapping>& a_dmap,
                                      const LPInfo& a_info,
                                      const Vector<FabFactory<FAB> const*>& a_factory)
{
    define(a_geom, a_grids, a_dmap, a_info, a_factory);
}

template <typename MF>
void
MLTerrainPoissonT<MF>::define (const Vector<Geometry>& a_geom,
                             const Vector<BoxArray>& a_grids,
                             const Vector<DistributionMapping>& a_dmap,
                             const LPInfo& a_info,
                             const Vector<FabFactory<FAB> const*>& a_factory)
{
    BL_PROFILE("MLTerrainPoisson::define()");
    MLCellABecLapT<MF>::define(a_geom, a_grids, a_dmap, a_info, a_factory);
}

template <typename MF>
MLTerrainPoissonT<MF>::~MLTerrainPoissonT () = default;

template <typename MF>
void
MLTerrainPoissonT<MF>::prepareForSolve ()
{
    BL_PROFILE("MLTerrainPoisson::prepareForSolve()");

    MLCellABecLapT<MF>::prepareForSolve();

    m_is_singular.clear();
    m_is_singular.resize(this->m_num_amr_levels, false);
    auto itlo = std::find(this->m_lobc[0].begin(), this->m_lobc[0].end(), BCType::Dirichlet);
    auto ithi = std::find(this->m_hibc[0].begin(), this->m_hibc[0].end(), BCType::Dirichlet);
    if (itlo == this->m_lobc[0].end() && ithi == this->m_hibc[0].end())
    if (!m_is_singular[0] && this->m_needs_coarse_data_for_bc &&
        this->m_coarse_fine_bc_type == LinOpBCType::Neumann)
    {
        auto bbox = this->m_grids[0][0].minimalBox();
        for (int idim = 0; idim < AMREX_SPACEDIM; ++idim) {
            if (this->m_lobc[0][idim] == LinOpBCType::Dirichlet) {
                bbox.growLo(idim,1);
            }
            if (this->m_hibc[0][idim] == LinOpBCType::Dirichlet) {
                bbox.growHi(idim,1);
            }
        }
        if (this->m_geom[0][0].Domain().contains(bbox)) {
            m_is_singular[0] = true;
        }
    }
}

template <typename MF>
template <typename AMF,
          std::enable_if_t<IsFabArray<AMF>::value &&
                           std::is_convertible_v<typename AMF::value_type,
                                                 typename  MF::value_type>, int>>

void
MLTerrainPoissonT<MF>::setZPhysNd (int amrlev, const AMF& zphys)
{
    AMREX_ASSERT_WITH_MESSAGE(zphys.nComp() == 1,
                              "MLTerrainPoisson::setZphysNd: zphys_nd is supposed to be single component.");
    m_a_zphys_nd.resize(this->m_num_amr_levels);
    m_a_zphys_nd[amrlev].define
        (convert(this->m_grids[amrlev][0],IntVect(0,0,1)),
                 this->m_dmap[amrlev][0], 1,
                 zphys.nGrowVect(),
                 MFInfo(), *(this->m_factory[amrlev][0]));
    m_a_zphys_nd[amrlev].LocalCopy(zphys, 0, 0, 1, zphys.nGrowVect());
}

template <typename MF>
void
MLTerrainPoissonT<MF>::Fapply (int amrlev, int mglev, MF& out, const MF& in) const
{
    BL_PROFILE("MLTerrainPoisson::Fapply()");

    const Real* dxinv = this->m_geom[amrlev][mglev].InvCellSize();

    AMREX_D_TERM(const RT dx = RT(dxinv[0]);,
                 const RT dy = RT(dxinv[1]);,
                 const RT dz = RT(dxinv[2]););

    AMREX_D_TERM(const RT dhx = RT(dxinv[0]*dxinv[0]);,
                 const RT dhy = RT(dxinv[1]*dxinv[1]);,
                 const RT dhz = RT(dxinv[2]*dxinv[2]););

#if (AMREX_SPACEDIM == 3)
    RT dh0 = this->get_d0(dhx, dhy, dhz);
    RT dh1 = this->get_d1(dhx, dhy, dhz);
#endif

#if (AMREX_SPACEDIM < 3)
    const RT dx = RT(this->m_geom[amrlev][mglev].CellSize(0));
    const RT probxlo = RT(this->m_geom[amrlev][mglev].ProbLo(0));
#endif

#ifdef AMREX_USE_GPU
    if (Gpu::inLaunchRegion() && out.isFusingCandidate() && !this->hasHiddenDimension()) {
        auto const& xma = in.const_arrays();
        auto const& yma = out.arrays();
        auto const& zpa = m_a_zphys_nd[amrlev].arrays();
        {
            ParallelFor(out,
            [=] AMREX_GPU_DEVICE (int box_no, int i, int j, int k) noexcept
            {
                amrex::ignore_unused(j,k);
                terrpoisson_adotx(AMREX_D_DECL(i,j,k), yma[box_no], xma[box_no], zpa[box_no],
                                  AMREX_D_DECL(dx,dy,dz));
            });
        }
        Gpu::streamSynchronize();
    } else
#endif
    {
#ifdef AMREX_USE_OMP
#pragma omp parallel if (Gpu::notInLaunchRegion())
#endif
        for (MFIter mfi(out, TilingIfNotGPU()); mfi.isValid(); ++mfi)
        {
            Box  const& bx = mfi.tilebox();
            auto const& xfab = in.array(mfi);
            auto const& yfab = out.array(mfi);
            auto const& zpa = m_a_zphys_nd[amrlev].array(mfi);
            {
                if (this->hasHiddenDimension()) {
                    Box const& bx2d = this->compactify(bx);
                    const auto& xfab2d = this->compactify(xfab);
                    const auto& yfab2d = this->compactify(yfab);
                    AMREX_HOST_DEVICE_PARALLEL_FOR_3D(bx2d, i, j, k,
                    {
                        amrex::ignore_unused(k);
                        // TwoD::terrpoisson_adotx(i, j, yfab2d, xfab2d, dh0, dh1);
                    });
                } else {
                    AMREX_HOST_DEVICE_PARALLEL_FOR_3D(bx, i, j, k,
                    {
                        terrpoisson_adotx(i, j, k, yfab, xfab, zpa, dx, dy, dz);
                    });
                }
            }
        } // mfi
    }
}

template <typename MF>
void
MLTerrainPoissonT<MF>::normalize (int /*amrlev*/, int /*mglev*/, MF& /*mf*/) const
{
    BL_PROFILE("MLTerrainPoisson::normalize()");
    amrex::Abort("What am I doing in MLTerrainPoisson::normalize??");
}

template <typename MF>
void
MLTerrainPoissonT<MF>::Fsmooth (int amrlev, int /*mglev*/, MF& /*sol*/, const MF& /*rhs*/, int /*redblack*/) const
{
    BL_PROFILE("MLTerrainPoisson::Fsmooth()");
    amrex::Abort("What am I doing in MLTerrainPoisson::Fsmooth??");
}

template <typename MF>
void
MLTerrainPoissonT<MF>::getFluxes (const Vector<Array<MF*,AMREX_SPACEDIM> >& a_flux,
                                  const Vector<MF*>& a_sol)
{
    BL_PROFILE("MLTerrainPoisson::getFluxes()");
    amrex::Abort("What am I doing in MLTerrainPoisson::getFluxes??");
}

template <typename MF>
void
MLTerrainPoissonT<MF>::FFlux (int amrlev, const MFIter& mfi,
                       const Array<FAB*,AMREX_SPACEDIM>& flux,
                       const FAB& sol, Location, const int /*face_only*/) const
{
    BL_PROFILE("MLTerrainPoisson::FFlux()");

    const int mglev = 0;
    const Box& box = mfi.tilebox();
    const Real* dxinv = this->m_geom[amrlev][mglev].InvCellSize();

    AMREX_D_TERM(const auto& fxarr = flux[0]->array();,
                 const auto& fyarr = flux[1]->array();,
                 const auto& fzarr = flux[2]->array(););

    const auto& solarr = sol.array();
    const auto& zpa    = m_a_zphys_nd[amrlev].array(mfi);

    RT dxinv = RT(dxinv[0]);
    RT dyinv = RT(dxinv[1]);
    RT dzinv = RT(dxinv[2]);

    if (this->hiddenDirection() != 0) {
        Box bflux = amrex::surroundingNodes(box, 0);
        AMREX_LAUNCH_HOST_DEVICE_LAMBDA ( bflux, tbox,
        {
            terrpoisson_flux_x(tbox, fxarr, solarr, zpa, invdx, invdz);
        });
    } else {
        flux[0]->template setVal<RunOn::Device>(RT(0.0));
    }
    if (this->hiddenDirection() != 1) {
        Box bflux = amrex::surroundingNodes(box, 1);
        AMREX_LAUNCH_HOST_DEVICE_LAMBDA ( bflux, tbox,
        {
            terrpoisson_flux_y(tbox, fyarr, solarr, zpa, invdy, invdz);
        });
    } else {
        flux[1]->template setVal<RunOn::Device>(RT(0.0));
    }
    if (this->hiddenDirection() != 2) {
        Box bflux = amrex::surroundingNodes(box, 2);
        AMREX_LAUNCH_HOST_DEVICE_LAMBDA ( bflux, tbox,
        {
            terrpoisson_flux_z(tbox, fzarr, solarr, zpa, invdx, invdy, invdz);
        });
    } else {
        flux[2]->template setVal<RunOn::Device>(RT(0.0));
    }
}

extern template class MLTerrainPoissonT<MultiFab>;

using MLTerrainPoisson = MLTerrainPoissonT<MultiFab>;

}

#endif
