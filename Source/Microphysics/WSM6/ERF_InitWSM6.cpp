#include "ERF_WSM6.H"

#include <AMReX_Gpu.H>
#include "ERF_EOS.H"

using namespace amrex;

void
WSM6::Init(const MultiFab& cons_in,
           const BoxArray&,
           const Geometry& geom,
           const Real& dt_advance,
           std::unique_ptr<MultiFab>& z_phys_nd,
           std::unique_ptr<MultiFab>& detJ_cc)
{
    dt = dt_advance;
    m_geom = geom;

    m_z_phys_nd = z_phys_nd.get();
    m_detJ_cc = detJ_cc.get();

    MicVarMap.resize(m_qmoist_size);
    MicVarMap = {MicVar_WSM6::rain_accum, MicVar_WSM6::snow_accum, MicVar_WSM6::graup_accum};

    for (int ivar = 0; ivar < MicVar_WSM6::NumVars; ++ivar) {
        mic_fab_vars[ivar] = std::make_shared<MultiFab>(cons_in.boxArray(), cons_in.DistributionMap(),
                                                        1, cons_in.nGrowVect());
        mic_fab_vars[ivar]->setVal(0.0);
    }

    nlev = m_geom.Domain().length(2);
    zlo = m_geom.Domain().smallEnd(2);
    zhi = m_geom.Domain().bigEnd(2);
}

void
WSM6::Copy_State_to_Micro(const MultiFab& cons_in)
{
    for (MFIter mfi(cons_in); mfi.isValid(); ++mfi) {
        // Match Morrison behavior: refresh microphysics ghost zones from state.
        // WSM6 Fortran reads the full (ims:ime, jms:jme, kms:kme) slab.
        const auto& box3d = mfi.growntilebox();
        auto states = cons_in.array(mfi);

        auto rho = mic_fab_vars[MicVar_WSM6::rho]->array(mfi);
        auto theta = mic_fab_vars[MicVar_WSM6::theta]->array(mfi);
        auto tabs = mic_fab_vars[MicVar_WSM6::tabs]->array(mfi);
        auto pres = mic_fab_vars[MicVar_WSM6::pres]->array(mfi);

        auto qv = mic_fab_vars[MicVar_WSM6::qv]->array(mfi);
        auto qc = mic_fab_vars[MicVar_WSM6::qc]->array(mfi);
        auto qi = mic_fab_vars[MicVar_WSM6::qi]->array(mfi);
        auto qr = mic_fab_vars[MicVar_WSM6::qr]->array(mfi);
        auto qs = mic_fab_vars[MicVar_WSM6::qs]->array(mfi);
        auto qg = mic_fab_vars[MicVar_WSM6::qg]->array(mfi);

        ParallelFor(box3d, [=] AMREX_GPU_DEVICE(int i, int j, int k) {
            rho(i,j,k) = states(i,j,k,Rho_comp);
            theta(i,j,k) = states(i,j,k,RhoTheta_comp) / states(i,j,k,Rho_comp);

            qv(i,j,k) = amrex::max(Real(0.0), states(i,j,k,RhoQ1_comp) / states(i,j,k,Rho_comp));
            qc(i,j,k) = amrex::max(Real(0.0), states(i,j,k,RhoQ2_comp) / states(i,j,k,Rho_comp));
            qi(i,j,k) = amrex::max(Real(0.0), states(i,j,k,RhoQ3_comp) / states(i,j,k,Rho_comp));
            qr(i,j,k) = amrex::max(Real(0.0), states(i,j,k,RhoQ4_comp) / states(i,j,k,Rho_comp));
            qs(i,j,k) = amrex::max(Real(0.0), states(i,j,k,RhoQ5_comp) / states(i,j,k,Rho_comp));
            qg(i,j,k) = amrex::max(Real(0.0), states(i,j,k,RhoQ6_comp) / states(i,j,k,Rho_comp));

            tabs(i,j,k) = getTgivenRandRTh(states(i,j,k,Rho_comp),
                                           states(i,j,k,RhoTheta_comp),
                                           qv(i,j,k));
            pres(i,j,k) = getPgivenRTh(states(i,j,k,RhoTheta_comp), qv(i,j,k));
        });
    }
}
