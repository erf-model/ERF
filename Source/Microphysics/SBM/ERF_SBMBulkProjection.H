#ifndef ERF_SBM_BULK_PROJECTION_H_
#define ERF_SBM_BULK_PROJECTION_H_

#include "ERF_SBMLayout.H"
#include "ERF_AuxiliaryFaceTransfer.H"

#include <AMReX_MultiFab.H>

#include <vector>

namespace erf_sbm {

struct BulkProjection {
    amrex::Real qc{0.0};
    amrex::Real qr{0.0};
};

/** Fixed cloud/rain projection from the first population's liquid mass only. */
class SBMBulkProjection {
public:
    explicit SBMBulkProjection(const SBMLayout& layout);

    [[nodiscard]] BulkProjection apply(const std::vector<amrex::Real>& auxiliary) const;
    void apply_to_core(const amrex::Box& box, const amrex::Array4<const amrex::Real>& auxiliary,
                       const amrex::Array4<amrex::Real>& core) const;
    void apply_to_face_transfer(const ::erf_auxiliary::AuxiliaryFaceTransfer& spectral,
                                ::erf_auxiliary::AuxiliaryFaceTransfer& bulk) const;
    [[nodiscard]] amrex::Real max_face_projection_error(
        const ::erf_auxiliary::AuxiliaryFaceTransfer& spectral,
        const ::erf_auxiliary::AuxiliaryFaceTransfer& bulk) const;

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
