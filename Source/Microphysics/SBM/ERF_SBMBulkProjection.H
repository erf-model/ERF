#ifndef ERF_SBM_BULK_PROJECTION_H_
#define ERF_SBM_BULK_PROJECTION_H_

#include "ERF_SBMLayout.H"

#include <AMReX_MultiFab.H>

#include <vector>

namespace erf_sbm {

struct BulkProjection {
    amrex::Real qc{0.0};
    amrex::Real qr{0.0};
};

/** Fixed cloud/rain projection from the liquid population's mass bins. */
class SBMBulkProjection {
public:
    explicit SBMBulkProjection(const SBMLayout& layout);

    [[nodiscard]] BulkProjection apply(const std::vector<amrex::Real>& auxiliary) const;
    void apply_to_core(const amrex::Box& box,
                       const amrex::Array4<const amrex::Real>& auxiliary,
                       const amrex::Array4<amrex::Real>& core,
                       int qc_component, int qr_component) const;

    [[nodiscard]] int cloud_offset() const noexcept { return m_cloud_offset; }
    [[nodiscard]] int cloud_count() const noexcept { return m_cloud_count; }
    [[nodiscard]] int rain_offset() const noexcept { return m_rain_offset; }
    [[nodiscard]] int rain_count() const noexcept { return m_rain_count; }

private:
    int m_cloud_offset{0};
    int m_cloud_count{0};
    int m_rain_offset{0};
    int m_rain_count{0};
};

} // namespace erf_sbm

#endif
