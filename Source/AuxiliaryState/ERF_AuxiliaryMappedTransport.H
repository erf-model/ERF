#ifndef ERF_AUXILIARY_MAPPED_TRANSPORT_H_
#define ERF_AUXILIARY_MAPPED_TRANSPORT_H_

#include <AMReX_Array4.H>
#include <AMReX.H>
#include <AMReX_GpuQualifiers.H>
#include <AMReX_MultiFab.H>

#include <array>
#include <memory>
#include <string>

namespace erf_auxiliary {

struct MappedFaceFluxRateTag {};
struct IntegratedMappedFaceFluxTag {};

/** A directionally face-centered mapped-coordinate field with a semantic tag. */
template <typename Tag>
class MappedFaceField {
public:
    MappedFaceField() = default;
    MappedFaceField(const MappedFaceField&) = delete;
    MappedFaceField& operator=(const MappedFaceField&) = delete;
    MappedFaceField(MappedFaceField&&) noexcept = default;
    MappedFaceField& operator=(MappedFaceField&&) noexcept = default;

    void define(const amrex::BoxArray& cell_ba,
                const amrex::DistributionMapping& dm,
                const int ncomp,
                const int ngrow = 0)
    {
        AMREX_ALWAYS_ASSERT(ncomp > 0);
        AMREX_ALWAYS_ASSERT(ngrow >= 0);
        m_ncomp = ncomp;
        for (int dir = 0; dir < AMREX_SPACEDIM; ++dir) {
            m_faces[static_cast<std::size_t>(dir)] = std::make_unique<amrex::MultiFab>(
                amrex::convert(cell_ba, amrex::IntVect::TheDimensionVector(dir)),
                dm, ncomp, ngrow);
        }
    }

    [[nodiscard]] bool is_defined() const noexcept { return m_ncomp > 0; }
    [[nodiscard]] int nComp() const noexcept { return m_ncomp; }

    amrex::MultiFab& dir(const int direction)
    {
        AMREX_ALWAYS_ASSERT(direction >= 0 && direction < AMREX_SPACEDIM);
        AMREX_ALWAYS_ASSERT(m_faces[static_cast<std::size_t>(direction)] != nullptr);
        return *m_faces[static_cast<std::size_t>(direction)];
    }

    const amrex::MultiFab& dir(const int direction) const
    {
        AMREX_ALWAYS_ASSERT(direction >= 0 && direction < AMREX_SPACEDIM);
        AMREX_ALWAYS_ASSERT(m_faces[static_cast<std::size_t>(direction)] != nullptr);
        return *m_faces[static_cast<std::size_t>(direction)];
    }

    void setVal(const amrex::Real value)
    {
        AMREX_ALWAYS_ASSERT(is_defined());
        for (auto& face : m_faces) { face->setVal(value); }
    }

private:
    std::array<std::unique_ptr<amrex::MultiFab>, AMREX_SPACEDIM> m_faces{};
    int m_ncomp{0};
};

using MappedFaceFluxRate = MappedFaceField<MappedFaceFluxRateTag>;
using IntegratedMappedFaceFlux = MappedFaceField<IntegratedMappedFaceFluxTag>;

/** Computational-coordinate divergence only: no metric or physical face area. */
AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE amrex::Real
ComputationalMappedDivergence(const amrex::Real fx_hi,
                              const amrex::Real fx_lo,
                              const amrex::Real fy_hi,
                              const amrex::Real fy_lo,
                              const amrex::Real fz_hi,
                              const amrex::Real fz_lo,
                              const amrex::Real dx_inv,
                              const amrex::Real dy_inv,
                              const amrex::Real dzeta_inv) noexcept
{
    return (fx_hi - fx_lo) * dx_inv +
           (fy_hi - fy_lo) * dy_inv +
           (fz_hi - fz_lo) * dzeta_inv;
}

/** Build omega=detJ/(mx*my), returning false and a diagnostic on invalid cells. */
bool BuildMappedCellMeasure(amrex::MultiFab& omega,
                            const amrex::MultiFab& detJ,
                            const amrex::MultiFab& mx,
                            const amrex::MultiFab& my,
                            std::string& diagnostic);

/** Check that a selected cell-centered component is finite and strictly positive. */
bool ValidatePositiveFiniteComponent(const amrex::MultiFab& field,
                                     int component,
                                     std::string& diagnostic);

/** Check that a selected cell-centered component is finite. */
bool ValidateFiniteComponent(const amrex::MultiFab& field,
                             int component,
                             std::string& diagnostic);

/**
 * Apply -D_xi(flux)/omega into one selected output component.  The supplied
 * face values are already in ERF's mapped coordinate convention.
 */
template <typename Tag>
void ApplyMappedFluxTendency(
    const MappedFaceField<Tag>& flux,
    const int flux_comp,
    const amrex::MultiFab& omega,
    const int omega_comp,
    amrex::MultiFab& tendency,
    const int tendency_comp,
    const amrex::GpuArray<amrex::Real, AMREX_SPACEDIM>& dx_inv)
{
    AMREX_ALWAYS_ASSERT(flux.is_defined());
    AMREX_ALWAYS_ASSERT(flux_comp >= 0 && flux_comp < flux.nComp());
    AMREX_ALWAYS_ASSERT(omega_comp >= 0 && omega_comp < omega.nComp());
    AMREX_ALWAYS_ASSERT(tendency_comp >= 0 && tendency_comp < tendency.nComp());

    for (amrex::MFIter mfi(tendency, amrex::TilingIfNotGPU()); mfi.isValid(); ++mfi) {
        const amrex::Box bx = mfi.tilebox();
        const auto fx = flux.dir(0).const_array(mfi);
        const auto fy = flux.dir(1).const_array(mfi);
        const auto fz = flux.dir(2).const_array(mfi);
        const auto measure = omega.const_array(mfi);
        const auto out = tendency.array(mfi);
        const auto dx = dx_inv[0];
        const auto dy = dx_inv[1];
        const auto dz = dx_inv[2];

        amrex::ParallelFor(bx, [=] AMREX_GPU_DEVICE(int i, int j, int k) noexcept {
            const amrex::Real div = ComputationalMappedDivergence(
                fx(i + 1, j, k, flux_comp), fx(i, j, k, flux_comp),
                fy(i, j + 1, k, flux_comp), fy(i, j, k, flux_comp),
                fz(i, j, k + 1, flux_comp), fz(i, j, k, flux_comp),
                dx, dy, dz);
            out(i, j, k, tendency_comp) = -div / measure(i, j, k, omega_comp);
        });
    }
}

/** Return the max absolute difference between two component-selected face fields. */
amrex::Real MaxFaceFieldDifference(const MappedFaceFluxRate& lhs,
                                  int lhs_comp,
                                  const MappedFaceFluxRate& rhs,
                                  int rhs_comp);

/** Scale every component and direction of a flux rate into an integrated ledger. */
void AccumulateIntegratedFaceFlux(IntegratedMappedFaceFlux& ledger,
                                 const MappedFaceFluxRate& rate,
                                 amrex::Real weight);

} // namespace erf_auxiliary

#endif
