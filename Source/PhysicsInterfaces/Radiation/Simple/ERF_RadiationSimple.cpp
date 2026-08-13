
#include <ERF_RadiationSimple.H>

#include <AMReX_MFIter.H>

using namespace amrex;

void RadiationSimple::Init(const amrex::Geometry& geom,
                           const amrex::BoxArray &ba,
                           amrex::MultiFab* cons_in)
{
    m_geom = geom;
    m_ba = ba;

    DistributionMapping dm = cons_in->DistributionMap();
    deltaq.define(ba, dm, 1, 0);
    flux.define(ba, dm, 1, IntVect(0, 0, 1));
    radlwdn = std::make_unique<amrex::MultiFab>(ba, dm, 1, 0);
    radqrlw = std::make_unique<amrex::MultiFab>(ba, dm, 1, 0);

    deltaq.setVal(0.0);
    flux.setVal(0.0);
    radlwdn->setVal(0.0);
    radqrlw->setVal(0.0);
}

void RadiationSimple::Run(int& /*level*/,
                          int& /*step*/,
                          double& /*time*/,
                          const double& /*dt*/,
                          const amrex::BoxArray& /*ba*/,
                          amrex::Geometry& geom,
                          amrex::MultiFab* cons_in,
                          amrex::iMultiFab* /*lmask*/,
                          amrex::MultiFab* /*t_surf*/,
                          amrex::Vector<amrex::MultiFab*>& /*lsm_input_ptrs*/,
                          amrex::Vector<amrex::MultiFab*>& /*lsm_output_ptrs*/,
                          amrex::MultiFab* qheating_rates,
                          amrex::MultiFab* rad_fluxes,
                          amrex::MultiFab* z_phys,
                          amrex::MultiFab* /*lat*/,
                          amrex::MultiFab* /*lon*/,
                          const bool /*updated_lsm*/)
{

    constexpr amrex::Real cp_spec = 1015.0;
    constexpr amrex::Real coef = 70.0;
    constexpr amrex::Real coef1 = 22.0;
    constexpr amrex::Real f0=3.75e-6;
    constexpr amrex::Real xk = 85.0;

    const Real fixed_dz = geom.CellSize(2);
    const int nz = geom.Domain().length(2);

    flux.setVal(0.0);
    deltaq.setVal(0.0);

    auto moist = m_moist;
    auto ice = m_ice;

    for (MFIter mfi(*cons_in, TileNoZ()); mfi.isValid(); ++mfi)
    {
        Box box = mfi.validbox();
        box.makeSlab(2, 0);

        const Array4<const Real>& cons_arr = cons_in->const_array(mfi);
        const Array4<const Real>& z_nd_arr = (z_phys) ? z_phys->const_array(mfi) : Array4<const Real>{};
        const Array4<Real>& deltaq_arr = deltaq.array(mfi);
        const Array4<Real>& flux_arr = flux.array(mfi);
        const Array4<Real>& radlwdn_arr = radlwdn->array(mfi);
        const Array4<Real>& radqrlw_arr = radqrlw->array(mfi);
        const Array4<Real>& qheating_arr = qheating_rates->array(mfi);

        const Array4<Real>& radfluxes_arr = rad_fluxes->array(mfi);

        ParallelFor(box, [=] AMREX_GPU_DEVICE (int i, int j, int)
        {
            int itop = 0;
            Real qzinf = 0.0;  // holds accumulated optical depth between z and domain top
            Real qzeroz = 0.0; // holds accumulated optical depth between surface and z

            for (int k = 0; k < nz; k++)
            {
                Real rho = cons_arr(i, j, k, Rho_comp);
                Real qv = (moist) ? cons_arr(i, j, k, RhoQ1_comp) / rho : Real(0);
                Real qc = (moist) ? cons_arr(i, j, k, RhoQ2_comp) / rho : Real(0);
                Real qi = (ice)   ? cons_arr(i, j, k, RhoQ3_comp) / rho : Real(0);

                Real dz = (z_nd_arr) ? Real(0.25) * ( (z_nd_arr(i  ,j  ,k+1) - z_nd_arr(i  ,j  ,k))
                                            + (z_nd_arr(i+1,j  ,k+1) - z_nd_arr(i+1,j  ,k))
                                            + (z_nd_arr(i  ,j+1,k+1) - z_nd_arr(i  ,j+1,k))
                                            + (z_nd_arr(i+1,j+1,k+1) - z_nd_arr(i+1,j+1,k)) ) : fixed_dz;

                // optical depth only includes that due to liquid water
                if (qc + qi > 0.0)
                {
                    deltaq_arr(i, j, k) = xk * rho * (qc + qi) * dz;
                }

                // inversion height is top of highest layer w/q>8g/kg
                if (qv + qc + qi > 0.008)
                {
                    itop = max(itop, k+1); // note zi(k+1) is inversion height
                }

                // accumulate optical depth in qzinf
                qzinf = qzinf + deltaq_arr(i, j, k);
            }

            // note: qzinf now holds total optical depth of column (due to cloud)
            //       qzeroz is initialized to zero since first level is at surface

            // compute net upward longwave flux up to inversion height
            for (int k = 0; k <= itop; k++)
            {
                flux_arr(i, j, k) = coef*exp(-qzinf) + coef1*exp(-qzeroz);
                qzinf = qzinf - deltaq_arr(i, j, k);
                qzeroz = qzeroz + deltaq_arr(i, j, k);
            }

            // compute net upward longwave flux above inversion height
            // this includes correction for clearsky fluxes which balances
            // the prescribed subsidence heating above the inversion
            Real z_itop = (z_nd_arr) ? z_nd_arr(i, j, itop) : itop * fixed_dz;
            for (int k = itop + 1; k < nz; k++) {
                Real rho_nd = 0.5 * (cons_arr(i, j, k - 1, Rho_comp) + cons_arr(i, j, k, Rho_comp)); // rho at interface
                Real zi = (z_nd_arr) ? z_nd_arr(i, j, k) : k * fixed_dz;
                flux_arr(i, j, k) = coef * exp(-qzinf) + coef1 * exp(-qzeroz)
                                    + cp_spec * rho_nd * f0 * (0.25 * zi + 0.75 * z_itop)
                                    * pow(zi - z_itop, 1.0 / 3.0);
                qzinf = qzinf - deltaq_arr(i, j, k);
                qzeroz = qzeroz + deltaq_arr(i, j, k);
            }
            Real rho_nd = 0.5 * (cons_arr(i, j, nz - 1, Rho_comp) + cons_arr(i, j, nz, Rho_comp));
            Real zi = (z_nd_arr) ? z_nd_arr(i, j, nz) : nz * fixed_dz;
            flux_arr(i, j, nz) = coef * exp(-qzinf) + coef1 * exp(-qzeroz)
                                 + cp_spec * rho_nd * f0 * (0.25 * zi + 0.75 * z_itop)
                                 * pow(zi - z_itop, 1.0 / 3.0);

            // note that our flux differs from the specification in that it
            // uses the local density (rather than the interface density) in
            // computing the correction to the flux above the inversion.
            // Formulating things in this way ensures a rough balance between
            // the subsidence heating and radiative cooling above the inversion.

            // compute radiative heating as divergence of net upward lw flux
            for (int k = 0; k < nz; k++)
            {
                Real dz = (z_nd_arr) ? 0.25 * ( (z_nd_arr(i  ,j  ,k+1) - z_nd_arr(i  ,j  ,k))
                                            + (z_nd_arr(i+1,j  ,k+1) - z_nd_arr(i+1,j  ,k))
                                            + (z_nd_arr(i  ,j+1,k+1) - z_nd_arr(i  ,j+1,k))
                                            + (z_nd_arr(i+1,j+1,k+1) - z_nd_arr(i+1,j+1,k)) ) : fixed_dz;

                Real cpmassl = cp_spec * cons_arr(i, j, k, Rho_comp) * dz; // thermal mass
                Real FTHRL = -(flux_arr(i, j, k+1) - flux_arr(i, j, k)) / cpmassl;
                // convert heating rate FTHRL for theta_d
                Real pres = getPgivenRTh(cons_arr(i,j,k,RhoTheta_comp),cons_arr(i,j,k,RhoQ1_comp)/cons_arr(i,j,k,Rho_comp));
                Real exner = getExnergivenP(pres, R_d/Cp_d);
                FTHRL *= exner;
                qheating_arr(i, j, k, 1) = FTHRL;          // radiative heating for source term
                radlwdn_arr(i, j, k) = flux_arr(i, j, k);  // net lw flux
                radqrlw_arr(i, j, k) = FTHRL;              // net lw heating

                // write lw dn to ERF fluxes
                radfluxes_arr(i, j, k, 3) = flux_arr(i, j, k);  // net lw flux
            }
        });
    }
}

void RadiationSimple::WriteDataLog(const double &time)
{
    constexpr int datwidth = 14;
    constexpr int datprecision = 9;
    constexpr int timeprecision = 13;

    Gpu::HostVector<Real> h_avg_radlwdn, h_avg_radqrlw;

    auto domain = m_geom.Domain();
    h_avg_radlwdn = sumToLine(*radlwdn, 0, 1, domain, 2);
    h_avg_radqrlw = sumToLine(*radqrlw, 0, 1, domain, 2);

    Real area_z = static_cast<Real>(domain.length(0)*domain.length(1));
    int nz = domain.length(2);
    for (int k = 0; k < nz; k++) {
        h_avg_radlwdn[k] /= area_z;
        h_avg_radqrlw[k] /= area_z;
    }

    if (ParallelDescriptor::IOProcessor()) {
        std::ostream& log = *datalog;

        if (log.good()) {
            for (int k = 0; k < nz; k++)
            {
                Real z = k * m_geom.CellSize(2);
                log << std::setw(datwidth) << std::setprecision(timeprecision) << time << " "
                    << std::setw(datwidth) << std::setprecision(datprecision) << z << " "
                    << h_avg_radlwdn[k] << " "
                    << h_avg_radqrlw[k] << std::endl;
            }
            // Write top face values
            Real z = nz * m_geom.CellSize(2);
            log << std::setw(datwidth) << std::setprecision(timeprecision) << time << " "
                << std::setw(datwidth) << std::setprecision(datprecision) << z << " "
                << 0.0 << " "
                << 0.0 << std::endl;
        }
    }
}
