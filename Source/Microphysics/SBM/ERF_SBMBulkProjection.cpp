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
    if (p == layout.populations().end()) {
        throw std::invalid_argument("SBM liquid projection references an unknown population");
    }
    m_cloud_offset = p->mass_offset;
    m_cloud_count = projection.cloud_rain_split;
    m_rain_offset = m_cloud_offset + m_cloud_count;
    m_rain_count = p->grid.nbins() - m_cloud_count;
}

BulkProjection SBMBulkProjection::apply(const std::vector<amrex::Real>& auxiliary) const
{
    if (m_rain_offset + m_rain_count > static_cast<int>(auxiliary.size())) {
        throw std::invalid_argument("auxiliary state is smaller than SBM projection");
    }
    BulkProjection result;
    for (int b = 0; b < m_cloud_count; ++b) {
        result.qc += auxiliary[static_cast<std::size_t>(m_cloud_offset + b)];
    }
    for (int b = 0; b < m_rain_count; ++b) {
        result.qr += auxiliary[static_cast<std::size_t>(m_rain_offset + b)];
    }
    return result;
}

void SBMBulkProjection::apply_to_core(const amrex::Box& box,
                                      const amrex::Array4<const amrex::Real>& auxiliary,
                                      const amrex::Array4<amrex::Real>& core,
                                      const int qc_component,
                                      const int qr_component) const
{
    const int cloud_offset = m_cloud_offset;
    const int cloud_count = m_cloud_count;
    const int rain_offset = m_rain_offset;
    const int rain_count = m_rain_count;
    amrex::ParallelFor(box, [=] AMREX_GPU_DEVICE (int i, int j, int k) noexcept {
        amrex::Real qc = amrex::Real(0.0);
        amrex::Real qr = amrex::Real(0.0);
        for (int b = 0; b < cloud_count; ++b) qc += auxiliary(i,j,k,cloud_offset+b);
        for (int b = 0; b < rain_count; ++b) qr += auxiliary(i,j,k,rain_offset+b);
        core(i,j,k,qc_component) = qc;
        core(i,j,k,qr_component) = qr;
    });
}

} // namespace erf_sbm
