/**
 * @file ERF_DustModelExample.cpp
 * @brief A minimal dust model, to check ERF's dust coupling (ERF_DustModel.H)
 *
 * It emits a constant mass flux (dust.example.flux [kg/m2/s]) from the ground
 * into the first dust scalar over the whole of level 0, and keeps one number
 * of its own, the mass it has emitted, which it writes to a CSV file
 * (dust.example.output) and to the checkpoint. It checks the coupling, not
 * dust physics: the registry, the sources at every Runge-Kutta stage, the
 * calls after each step and at output time, and the checkpoint and restart.
 *
 * With walls at the top and bottom and periodic sides nothing leaves the
 * domain, so the dust ERF reports (RHO DUST) equals the mass emitted.
 */

#include <ERF.H>
#include <ERF_DustModel.H>
#include <ERF_DustModelRegistry.H>
#include <ERF_IndexDefines.H>

#include <AMReX_ParmParse.H>

#include <fstream>
#include <iomanip>
#include <limits>
#include <string>

namespace {

class DustModelExample : public DustModel
{
public:
    const char* name () const override { return "example"; }

    void initialize (const DustInit& init) override
    {
        amrex::ParmParse pp("dust.example");
        pp.query("flux", m_flux);
        pp.query("output", m_output);
        if (!(m_flux >= 0.0)) {
            amrex::Abort("dust.example.flux must be >= 0 kg/m2/s, got " + std::to_string(m_flux));
        }

        const amrex::Geometry& g = init.erf.Geom(0);
        m_area = (g.ProbHi(0) - g.ProbLo(0)) * (g.ProbHi(1) - g.ProbLo(1));

        if (init.restart_chkfile.empty()) {
            m_emitted = 0.0;
            if (amrex::ParallelDescriptor::IOProcessor()) {
                std::ofstream csv(m_output, std::ios::trunc);
                csv << "step,time,emitted\n";
            }
        } else {
            // every rank reads the one number the checkpoint holds
            std::ifstream in(init.restart_chkfile + "/DustExample");
            if (!(in >> m_emitted)) {
                amrex::Abort("dust model example: no state in " + init.restart_chkfile + "/DustExample");
            }
        }
    }

    void add_sources (int lev, amrex::MultiFab& cc_src, const amrex::MultiFab& /*state*/,
                      const amrex::MultiFab& detJ_cc, const amrex::Geometry& geom,
                      amrex::Real /*dt*/) override
    {
        if (lev != 0) { return; }
        const int klo = geom.Domain().smallEnd(2);
        const amrex::Real dz = geom.CellSize(2);
        const amrex::Real flux = m_flux;
        for (amrex::MFIter mfi(cc_src); mfi.isValid(); ++mfi) {
            amrex::Box bx = mfi.validbox();
            if (bx.smallEnd(2) != klo) { continue; }
            bx.setBig(2, klo);
            const auto& src = cc_src.array(mfi);
            const auto& dJ  = detJ_cc.const_array(mfi);
            amrex::ParallelFor(bx, [=] AMREX_GPU_DEVICE (int i, int j, int k) noexcept
            {
                // the flux into the bottom cell, per unit volume: its height is detJ dz
                src(i, j, k, RhoDust_comp) += flux / (dJ(i, j, k) * dz);
            });
        }
    }

    void advance (const DustStep& step) override
    {
        m_emitted += m_flux * m_area * step.dt;
    }

    void write_output (int step, double time, bool /*is_final*/) override
    {
        if (amrex::ParallelDescriptor::IOProcessor()) {
            std::ofstream csv(m_output, std::ios::app);
            csv << step << "," << std::setprecision(std::numeric_limits<amrex::Real>::max_digits10)
                << time << "," << m_emitted << "\n";
        }
    }

    void write_checkpoint (const std::string& checkpointname) const override
    {
        if (amrex::ParallelDescriptor::IOProcessor()) {
            std::ofstream out(checkpointname + "/DustExample");
            out << std::setprecision(std::numeric_limits<amrex::Real>::max_digits10) << m_emitted << "\n";
        }
    }

private:
    amrex::Real m_flux    = 1.0e-4;   ///< dust.example.flux [kg/m2/s]
    std::string m_output  = "dust_example.csv";
    amrex::Real m_area    = 0.0;      ///< area of the domain [m2]
    amrex::Real m_emitted = 0.0;      ///< mass emitted so far [kg]
};

} // namespace

ERF_REGISTER_DUST_MODEL(DustModelExample, "example")
