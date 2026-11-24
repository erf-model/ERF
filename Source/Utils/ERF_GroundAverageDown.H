#ifndef ERF_GROUND_AVERAGE_DOWN_H_
#define ERF_GROUND_AVERAGE_DOWN_H_

#include <AMReX_MultiFab.H>

/**
 * Restriction of a quantity held in the ground layer of cells (a surface accumulation) from a
 * fine level to the coarse level below it.  Each coarse ground cell under a fine box that stands
 * on the ground gets a scale times the sum of the rr[0]*rr[1] fine ground cells above it: the
 * inverse area for a depth, 1 for an amount per cell.  The coarse ground boxes, distributed with
 * the fine boxes they lie under, are built once and serve every quantity restricted with them.
 */
class GroundAverageDown
{
public:
    /**
     * @param a_fba, a_fdm BoxArray and DistributionMapping of the fine data
     * @param a_rr         refinement ratio between the two levels
     * @param a_kf, a_kc   index of the ground layer on the fine and on the coarse level
     */
    GroundAverageDown (const amrex::BoxArray& a_fba, const amrex::DistributionMapping& a_fdm,
                       const amrex::IntVect& a_rr, int a_kf, int a_kc)
        : m_fba(a_fba), m_fdm(a_fdm), m_rr(a_rr), m_kf(a_kf)
    {
        amrex::BoxList bl;
        amrex::Vector<int> pmap;
        for (int ib = 0; ib < a_fba.size(); ++ib) {
            if (a_fba[ib].smallEnd(2) != a_kf) { continue; }
            amrex::Box s = amrex::coarsen(a_fba[ib], a_rr);
            s.setRange(2, a_kc);
            bl.push_back(s);
            m_fine_box.push_back(ib);
            pmap.push_back(a_fdm[ib]);
        }
        if (!bl.isEmpty()) {
            m_tmp.define(amrex::BoxArray(std::move(bl)), amrex::DistributionMapping(std::move(pmap)), 1, 0);
        }
    }

    /**
     * Component a_ccomp of a_crse = a_scale times the sum of the fine ground cells of component
     * a_fcomp of a_fine, in the coarse ground cells under the fine level.  a_fine must be defined
     * on the fine BoxArray and DistributionMapping given to the constructor.
     */
    void apply (const amrex::MultiFab& a_fine, amrex::MultiFab& a_crse, const amrex::Real a_scale,
                const int a_fcomp = 0, const int a_ccomp = 0)
    {
        AMREX_ASSERT(a_fine.boxArray() == m_fba && a_fine.DistributionMap() == m_fdm);
        if (m_fine_box.empty()) { return; }
        const amrex::IntVect rr = m_rr;
        const int kf = m_kf;
        for (amrex::MFIter mfi(m_tmp); mfi.isValid(); ++mfi) {
            const auto f = a_fine.const_array(m_fine_box[mfi.index()], a_fcomp);
            const auto c = m_tmp.array(mfi);
            amrex::ParallelFor(mfi.validbox(), [=] AMREX_GPU_DEVICE (int ic, int jc, int kc) noexcept
            {
                amrex::Real sum = 0.0;
                for (int j = jc*rr[1]; j < (jc+1)*rr[1]; ++j) {
                    for (int i = ic*rr[0]; i < (ic+1)*rr[0]; ++i) {
                        sum += f(i,j,kf);
                    }
                }
                c(ic,jc,kc) = sum * a_scale;
            });
        }
        a_crse.ParallelCopy(m_tmp, 0, a_ccomp, 1);
    }

private:
    amrex::BoxArray m_fba;
    amrex::DistributionMapping m_fdm;
    amrex::IntVect m_rr;
    int m_kf;
    amrex::Vector<int> m_fine_box;
    amrex::MultiFab m_tmp;
};

/**
 * The reverse of GroundAverageDown: each ground cell of a_fine (component a_fcomp) gets a_scale
 * times the coarse ground cell of a_crse (component a_ccomp) it lies over -- 1/(rr[0]*rr[1]) for
 * an amount per cell, which keeps the amount under each coarse cell.  a_crse must cover the fine
 * ground boxes, coarsened.
 */
inline void
groundFillFromCoarse (amrex::MultiFab& a_fine, const amrex::MultiFab& a_crse, const amrex::IntVect& a_rr,
                      const int a_kf, const int a_kc, const amrex::Real a_scale,
                      const int a_fcomp = 0, const int a_ccomp = 0)
{
    const amrex::BoxArray& fba = a_fine.boxArray();
    amrex::BoxList bl;
    amrex::Vector<int> fine_box, pmap;
    for (int ib = 0; ib < fba.size(); ++ib) {
        if (fba[ib].smallEnd(2) != a_kf) { continue; }
        amrex::Box s = amrex::coarsen(fba[ib], a_rr);
        s.setRange(2, a_kc);
        bl.push_back(s);
        fine_box.push_back(ib);
        pmap.push_back(a_fine.DistributionMap()[ib]);
    }
    if (bl.isEmpty()) { return; }
    amrex::MultiFab tmp(amrex::BoxArray(std::move(bl)), amrex::DistributionMapping(std::move(pmap)), 1, 0);
    tmp.setVal(0.0);
    tmp.ParallelCopy(a_crse, a_ccomp, 0, 1);
    const amrex::IntVect rr = a_rr;
    for (amrex::MFIter mfi(tmp); mfi.isValid(); ++mfi) {
        const int ib = fine_box[mfi.index()];
        const auto c = tmp.const_array(mfi);
        const auto f = a_fine.array(ib, a_fcomp);
        amrex::Box gbx = fba[ib];
        gbx.setRange(2, a_kf);
        amrex::ParallelFor(gbx, [=] AMREX_GPU_DEVICE (int i, int j, int k) noexcept
        {
            // floor division: the coarse cell over fine cell (i, j)
            const int ic = (i >= 0) ? i / rr[0] : -((-i - 1) / rr[0]) - 1;
            const int jc = (j >= 0) ? j / rr[1] : -((-j - 1) / rr[1]) - 1;
            f(i,j,k) = a_scale * c(ic,jc,a_kc);
        });
    }
}

#endif
