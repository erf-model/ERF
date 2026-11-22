
#include <ERF_RadiationSimple.H>
#include "ERF_Constants.H"
#include "ERF_TerrainMetrics.H"

#include <AMReX_MFIter.H>

#include <limits>

using namespace amrex;

void RadiationSimple::Init (const amrex::Geometry& geom,
                            const amrex::BoxArray &ba,
                            amrex::MultiFab* cons_in)
{
    m_geom = geom;
    m_ba = ba;

    // Set (or cleared) on every Init, so a level that spans the column again after a
    // regrid resumes the sweep.
    m_partial_columns = !rad_boxes_span_column(ba, geom.Domain());

    DistributionMapping dm = cons_in->DistributionMap();
    deltaq.define(ba, dm, 1, 0);
    flux.define(ba, dm, 1, IntVect(0, 0, 1));
    deltaq.setVal(zero);
    flux.setVal(zero);

    // Run() skips a level made of partial columns (and aborts on level 0), so the datalog
    // is never filled there; drop any MultiFabs left over from an earlier grid.
    if (datalog_int > 0 && !m_partial_columns) {
        datalog_mf.define(ba, dm, n_datalog_comp, 0);
        datalog_mf.setVal(zero);

        // One cell thick in z; each box spans the full column (m_partial_columns is false)
        BoxList bl_bdy = ba.boxList();
        for (auto& b : bl_bdy) { b.setRange(2, 0, 1); }
        datalog_bdy_mf.define(BoxArray(std::move(bl_bdy)), dm, n_datalog_bdy_comp, 0);
        datalog_bdy_mf.setVal(zero);
    } else {
        datalog_mf.clear();
        datalog_bdy_mf.clear();
    }
}

void RadiationSimple::Run (int& level,
                           int& /*step*/,
                           double& /*time*/,
                           const double& /*dt*/,
                           const amrex::BoxArray& /*ba*/,
                           amrex::Geometry& geom,
                           amrex::MultiFab* cons_in,
                           amrex::iMultiFab* /*lmask*/,
                           amrex::MultiFab* /*t_surf*/,
                           const amrex::Vector<const amrex::MultiFab*>& /*lsm_input_ptrs*/,
                           const amrex::Vector<amrex::MultiFab*>& /*lsm_output_ptrs*/,
                           amrex::MultiFab* qheating_rates,
                           amrex::MultiFab* /*rad_fluxes*/,
                           amrex::MultiFab* z_phys,
                           amrex::MultiFab* z_sfc,
                           amrex::MultiFab* /*lat*/,
                           amrex::MultiFab* /*lon*/,
                           const bool /*force_update*/,
                           amrex::Real& /*declin*/,
                           amrex::Real& /*calday*/)
{
    if (m_partial_columns) {
        // A level 0 split in z is rejected when its grids are built (ERF::check_rad_level_grids)
        AMREX_ASSERT(level > 0);
        amrex::ignore_unused(level);
        return;
    }

    constexpr amrex::Real cp_spec = 1015.0;
    constexpr amrex::Real coef = 70.0;
    constexpr amrex::Real coef1 = 22.0;
    constexpr amrex::Real f0=3.75e-6;
    constexpr amrex::Real xk = 85.0;

    const Real fixed_dz = geom.CellSize(2);
    const int nz = geom.Domain().length(2);

    flux.setVal(zero);
    deltaq.setVal(zero);

    const bool have_qv = m_have_qv;
    const bool have_qc = m_have_qc;
    const bool ice     = m_ice;
    const int  qv_comp = m_qv_comp;
    const int  qc_comp = m_qc_comp;
    const int  qi_comp = m_qi_comp;

    const bool do_datalog = (datalog_int > 0);

    // Heights above ground need both the node heights and the surface heights
    const bool has_z = (z_phys != nullptr) && (z_sfc != nullptr);

    for (MFIter mfi(*cons_in, TileNoZ()); mfi.isValid(); ++mfi)
    {
        Box box = mfi.tilebox();
        box.makeSlab(2, 0);

        const Array4<const Real>& cons_arr = cons_in->const_array(mfi);
        const Array4<const Real>& z_nd_arr = (z_phys) ? z_phys->const_array(mfi) : Array4<const Real>{};
        const Array4<Real>& deltaq_arr = deltaq.array(mfi);
        const Array4<Real>& flux_arr = flux.array(mfi);
        const Array4<Real>& qheating_arr = qheating_rates->array(mfi);
        const Array4<const Real>& zsfc_arr = (has_z) ? z_sfc->const_array(mfi) : Array4<const Real>{};
        const Array4<Real>& dl_arr  = (do_datalog) ? datalog_mf.array(mfi)     : Array4<Real>{};
        const Array4<Real>& bdy_arr = (do_datalog) ? datalog_bdy_mf.array(mfi) : Array4<Real>{};

        ParallelFor(box, [=] AMREX_GPU_DEVICE (int i, int j, int)
        {
            int itop = 0;
            Real qzinf = 0.0;  // holds accumulated optical depth between z and domain top
            Real qzeroz = 0.0; // holds accumulated optical depth between surface and z

            for (int k = 0; k < nz; k++)
            {
                Real rho = cons_arr(i, j, k, Rho_comp);
                Real qv = (have_qv) ? cons_arr(i, j, k, qv_comp) / rho : zero;
                Real qc = (have_qc) ? cons_arr(i, j, k, qc_comp) / rho : zero;
                Real qi = (ice)     ? cons_arr(i, j, k, qi_comp) / rho : zero;

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

            // compute net upward longwave flux up to inversion height.  When the
            // moist layer reaches the top cell, itop = nz: the face at nz is then
            // written below, where its above-inversion term vanishes (zi = z_itop),
            // and deltaq has no cell at nz to read.
            for (int k = 0; k <= amrex::min(itop, nz - 1); k++)
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
                Real qv_k = (have_qv) ? cons_arr(i,j,k,qv_comp)/cons_arr(i,j,k,Rho_comp) : zero;
                Real pres = getPgivenRTh(cons_arr(i,j,k,RhoTheta_comp), qv_k);
                Real exner = getExnergivenP(pres, R_d/Cp_d);
                FTHRL /= exner;
                qheating_arr(i, j, k, 1) = FTHRL;          // radiative heating for source term

                // Datalog profiles at the cell center; the net flux is averaged from the faces
                if (do_datalog) {
                    dl_arr(i, j, k, dl_z_comp)       = (has_z) ? Compute_Zrel_AtCellCenter(i,j,k,zsfc_arr,z_nd_arr)
                                                               : (k + myhalf)*fixed_dz;
                    dl_arr(i, j, k, dl_flux_comp)    = myhalf * (flux_arr(i, j, k) + flux_arr(i, j, k+1));
                    dl_arr(i, j, k, dl_heating_comp) = FTHRL;
                }

                // The net flux is not written to rad_fluxes: that array holds
                // (SW up, SW down, LW up, LW down), and this model computes
                // none of those. It stays zero, and HasRadiation() is false.
            }

            // Net flux on the surface and top-of-atmosphere faces, and the top-face height
            if (do_datalog) {
                bdy_arr(i, j, 0, bdy_sfc_comp)  = flux_arr(i, j, 0);
                bdy_arr(i, j, 0, bdy_toa_comp)  = flux_arr(i, j, nz);
                bdy_arr(i, j, 0, bdy_ztop_comp) = (has_z) ?
                    fourth*( z_nd_arr(i,j,nz) + z_nd_arr(i+1,j,nz) + z_nd_arr(i,j+1,nz) + z_nd_arr(i+1,j+1,nz) )
                  - fourth*( zsfc_arr(i,j,0) + zsfc_arr(i+1,j,0) + zsfc_arr(i,j+1,0) + zsfc_arr(i+1,j+1,0) )
                  : nz*fixed_dz;
            }
        });
    }
}

void RadiationSimple::WriteDataLog (const double &time)
{
    auto domain = m_geom.Domain();
    const Real area_z = static_cast<Real>(domain.length(0)*domain.length(1));
    const int nz = domain.length(2);

    // Planar averages of the cell-centered profiles
    Vector<Gpu::HostVector<Real>> cc_avg(n_datalog_comp);
    for (int n = 0; n < n_datalog_comp; ++n) {
        cc_avg[n] = sumToLine(datalog_mf, n, 1, domain, 2);
        for (auto& val : cc_avg[n]) { val /= area_z; }
    }

    // Planar averages of the surface and top-of-atmosphere values
    Vector<Real> bdy_avg(n_datalog_bdy_comp);
    for (int n = 0; n < n_datalog_bdy_comp; ++n) { bdy_avg[n] = datalog_bdy_mf.sum(n, true); }
    ParallelDescriptor::ReduceRealSum(bdy_avg.data(), n_datalog_bdy_comp);
    for (auto& val : bdy_avg) { val /= area_z; }

    if (ParallelDescriptor::IOProcessor()) {
        std::ostream& log = *datalog;
        if (log.good()) {
            // The log is opened in append mode, so only a new file gets the header
            log.seekp(0, std::ios::end);
            if (log.tellp() == 0) {
                log << "# Planar averages; each output is a block of nz+2 rows: the surface face,\n"
                    << "# the nz cell centers, then the top-of-atmosphere face. z_agl is the height\n"
                    << "# above ground. The heating rate exists only at cell centers and is nan on\n"
                    << "# the face rows.\n"
                    << "# time z_agl lw_flux_net lw_heating" << std::endl;
            }

            const Real nan = std::numeric_limits<Real>::quiet_NaN();

            // Row -1 is the surface face, rows 0..nz-1 the cell centers, row nz the top face
            for (int row = -1; row <= nz; ++row)
            {
                Real z_agl, lw_flux, lw_heating;
                if (row < 0) {
                    z_agl      = zero;
                    lw_flux    = bdy_avg[bdy_sfc_comp];
                    lw_heating = nan;
                } else if (row == nz) {
                    z_agl      = bdy_avg[bdy_ztop_comp];
                    lw_flux    = bdy_avg[bdy_toa_comp];
                    lw_heating = nan;
                } else {
                    z_agl      = cc_avg[dl_z_comp][row];
                    lw_flux    = cc_avg[dl_flux_comp][row];
                    lw_heating = cc_avg[dl_heating_comp][row];
                }

                log << std::setw(datalog_width) << std::setprecision(datalog_time_precision) << time << " "
                    << std::setw(datalog_width) << std::setprecision(datalog_precision) << z_agl << " "
                    << lw_flux << " "
                    << lw_heating << std::endl;
            } // for row
        } // log is good
    } // IOProc
}
