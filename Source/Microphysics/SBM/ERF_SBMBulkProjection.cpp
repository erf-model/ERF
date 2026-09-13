#include "ERF_SBMBulkProjection.H"

#include "ERF_IndexDefines.H"

#include <AMReX_MFParallelFor.H>

#include <algorithm>
#include <stdexcept>

namespace erf_sbm {

SBMBulkProjection::SBMBulkProjection(const SBMLayout& layout)
{
    const auto& projection = layout.liquid_projection();
    const auto p = std::find_if(layout.populations().begin(), layout.populations().end(),
        [&](const PopulationLayout& candidate) {
            return candidate.population_id == projection.population_id;
        });
    m_cloud_offset = p->mass_offset;
    m_cloud_count = projection.cloud_rain_split;
    m_rain_offset = m_cloud_offset + m_cloud_count;
    m_rain_count = p->grid.nbins() - m_cloud_count;
}

void SBMBulkProjection::apply_to_face_transfer(
    const ::erf_auxiliary::AuxiliaryFaceTransfer& spectral,
    ::erf_auxiliary::AuxiliaryFaceTransfer& bulk) const
{
    if (bulk.ncomp() != 2 || spectral.ncomp() < m_rain_offset + m_rain_count) {
        throw std::invalid_argument("face-transfer projection has incompatible component counts");
    }
    for (int dir = 0; dir < AMREX_SPACEDIM; ++dir) {
        const auto& source = spectral.direction(dir);
        auto& target = bulk.direction(dir);
        for (amrex::MFIter mfi(target); mfi.isValid(); ++mfi) {
            const amrex::Box box = mfi.validbox();
            const auto in = source.const_array(mfi);
            const auto out = target.array(mfi);
            const int cloud_offset = m_cloud_offset;
            const int cloud_count = m_cloud_count;
            const int rain_offset = m_rain_offset;
            const int rain_count = m_rain_count;
            amrex::ParallelFor(box, [=] AMREX_GPU_DEVICE (int i, int j, int k) noexcept {
                amrex::Real qc = amrex::Real(0.0);
                amrex::Real qr = amrex::Real(0.0);
                for (int b = 0; b < cloud_count; ++b) qc += in(i,j,k,cloud_offset+b);
                for (int b = 0; b < rain_count; ++b) qr += in(i,j,k,rain_offset+b);
                out(i,j,k,0) = qc;
                out(i,j,k,1) = qr;
            });
        }
    }
}

amrex::Real SBMBulkProjection::max_face_projection_error(
    const ::erf_auxiliary::AuxiliaryFaceTransfer& spectral,
    const ::erf_auxiliary::AuxiliaryFaceTransfer& bulk) const
{
    if (bulk.ncomp() != 2 || spectral.ncomp() < m_rain_offset + m_rain_count) {
        throw std::invalid_argument("face-transfer error check has incompatible component counts");
    }
    amrex::Real maximum = amrex::Real(0.0);
    for (int dir = 0; dir < AMREX_SPACEDIM; ++dir) {
        const auto& source = spectral.direction(dir);
        const auto& target = bulk.direction(dir);
        amrex::MultiFab error(target.boxArray(), target.DistributionMap(), 1, target.nGrowVect());
        error.setVal(amrex::Real(0.0));
        for (amrex::MFIter mfi(target); mfi.isValid(); ++mfi) {
            const amrex::Box box = mfi.validbox();
            const auto in = source.const_array(mfi);
            const auto out = target.const_array(mfi);
            const auto err = error.array(mfi);
            const int cloud_offset = m_cloud_offset;
            const int cloud_count = m_cloud_count;
            const int rain_offset = m_rain_offset;
            const int rain_count = m_rain_count;
            amrex::ParallelFor(box, [=] AMREX_GPU_DEVICE (int i, int j, int k) noexcept {
                amrex::Real qc = amrex::Real(0.0);
                amrex::Real qr = amrex::Real(0.0);
                for (int b = 0; b < cloud_count; ++b) qc += in(i,j,k,cloud_offset+b);
                for (int b = 0; b < rain_count; ++b) qr += in(i,j,k,rain_offset+b);
                err(i,j,k) = amrex::max(amrex::Math::abs(out(i,j,k,0) - qc),
                                         amrex::Math::abs(out(i,j,k,1) - qr));
            });
        }
        maximum = amrex::max(maximum, error.max(0));
    }
    return maximum;
}

BulkProjection SBMBulkProjection::apply(const std::vector<amrex::Real>& auxiliary) const
{
    if (m_rain_offset + m_rain_count > static_cast<int>(auxiliary.size())) {
        throw std::invalid_argument("auxiliary state is smaller than SBM projection");
    }
    BulkProjection result;
    for (int b = 0; b < m_cloud_count; ++b) result.qc += auxiliary[static_cast<std::size_t>(m_cloud_offset + b)];
    for (int b = 0; b < m_rain_count; ++b) result.qr += auxiliary[static_cast<std::size_t>(m_rain_offset + b)];
    return result;
}

void SBMBulkProjection::apply_to_core(const amrex::Box& box,
                                      const amrex::Array4<const amrex::Real>& auxiliary,
                                      const amrex::Array4<amrex::Real>& core) const
{
    const int cloud_offset = m_cloud_offset;
    const int cloud_count = m_cloud_count;
    const int rain_offset = m_rain_offset;
    const int rain_count = m_rain_count;
    amrex::ParallelFor(box, [=] AMREX_GPU_DEVICE (int i, int j, int k) noexcept {
        amrex::Real qc = 0.0;
        amrex::Real qr = 0.0;
        for (int b = 0; b < cloud_count; ++b) qc += auxiliary(i,j,k,cloud_offset+b);
        for (int b = 0; b < rain_count; ++b) qr += auxiliary(i,j,k,rain_offset+b);
        core(i,j,k,RhoQ2_comp) = qc;
        core(i,j,k,RhoQ3_comp) = qr;
    });
}

} // namespace erf_sbm
