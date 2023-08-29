#ifndef ABLMOST_H
#define ABLMOST_H

#include <AMReX_Geometry.H>
#include <AMReX_ParmParse.H>
#include <AMReX_FArrayBox.H>
#include <AMReX_MultiFab.H>
#include <AMReX_iMultiFab.H>

#include <IndexDefines.H>
#include <ERF_Constants.H>
#include <MOSTAverage.H>

/** Monin-Obukhov surface layer profile
 *
 * van der Laan, P., Kelly, M. C., & Sørensen, N. N. (2017). A new k-epsilon
 * model consistent with Monin-Obukhov similarity theory. Wind Energy,
 * 20(3), 479–489. https://doi.org/10.1002/we.2017
 *
 * Consistent with Dyer (1974) formulation from page 57, Chapter 2, Modeling
 * the vertical ABL structure in Modelling of Atmospheric Flow Fields,
 * Demetri P Lalas and Corrado F Ratto, January 1996,
 * https://doi.org/10.1142/2975.
 */

class ABLMostData
{
public:

    amrex::Real z0_const{0.1};       ///< Roughness height -- default constant value(m)
    amrex::Real kappa{KAPPA};        ///< von Karman constant
    amrex::Real gravity{CONST_GRAV}; ///< Acceleration due to gravity (m/s^2)
    amrex::Real surf_temp_flux{0.0}; ///< Heat flux
    amrex::Real surf_temp;           ///< Surface temperature
    amrex::Real beta_m{5.0};         ///< Constants from Dyer, BLM, 1974
    amrex::Real beta_h{5.0};         ///< https://doi.org/10.1007/BF00240838
    amrex::Real Cnk_a{0.0185};       ///< Standard Charnock constant https://doi.org/10.1175/JAMC-D-17-0137.1
    amrex::Real Cnk_b1{1.0/30.0};    ///< Modified Charnock Eq (4) https://doi.org/10.1175/JAMC-D-17-0137.1
    amrex::Real Cnk_b2{1260.0};      ///< Modified Charnock Eq (4) https://doi.org/10.1175/JAMC-D-17-0137.1
    amrex::Real Cnk_d{30.0};         ///< Modified Charnock Eq (4) https://doi.org/10.1175/JAMC-D-17-0137.1
    amrex::Real gamma_m{16.0};
    amrex::Real gamma_h{16.0};

    /**
     * Function to compute psi_m.
     *
     * @param[in] zeta Ratio of the query heigh to the Obukhov length scale
     */
    AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE
    amrex::Real calc_psi_m(amrex::Real zeta) const
    {
        if (zeta > 0) {
            return -beta_m * zeta;
        } else {
            amrex::Real x = std::sqrt(std::sqrt(1.0 - gamma_m * zeta));
            return 2.0 * std::log(0.5 * (1.0 + x)) + log(0.5 * (1.0 + x * x)) -
                   2.0 * std::atan(x) + PIoTwo;
        }
    }

    /**
     * Function to compute psi_h.
     *
     * @param[in] zeta Ratio of the query heigh to the Obukhov length scale
     */
    AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE
    amrex::Real calc_psi_h(amrex::Real zeta) const
    {
        if (zeta > 0) {
            return -beta_h * zeta;
        } else {
            amrex::Real x = std::sqrt(1.0 - gamma_h * zeta);
            return 2.0 * std::log(0.5 * (1.0 + x));
        }
    }
};

class ABLMost : public ABLMostData
{

public:

    // Constructor
    explicit ABLMost(const amrex::Vector<amrex::Geometry>& geom,
                     amrex::Vector<amrex::Vector<amrex::MultiFab>>& vars_old,
                     amrex::Vector<std::unique_ptr<amrex::MultiFab>>& Theta_prim,
                     amrex::Vector<std::unique_ptr<amrex::MultiFab>>& z_phys_nd)
      : m_geom(geom), m_ma(geom,vars_old,Theta_prim,z_phys_nd)
    {
        amrex::ParmParse pp("erf");
        pp.query("most.z0" , z0_const);

        // Specify the momentum calc type
        std::string mom_type_in;
        auto erf_mom = pp.query("most.roughness_type", mom_type_in);
        if (erf_mom) {
            if (mom_type_in == "constant") {
                mom_type = MomCalcType::DEFAULT;
            } else if (mom_type_in == "charnock") {
                mom_type = MomCalcType::CHARNOCK;
                pp.query("most.charnock_constant", Cnk_a);
            } else if (mom_type_in == "modified_charnock") {
                mom_type = MomCalcType::MODIFIED_CHARNOCK;
                pp.query("most.modified_charnock_depth", Cnk_d);
            } else {
                amrex::Abort("MOST momentum calc type not recognized!");
            }
        } else {
            mom_type = MomCalcType::DEFAULT;
        }

        // Specify surface temperature or surface flux
        auto erf_st = pp.query("most.surf_temp", surf_temp);
        if (erf_st) {
            theta_type = ThetaCalcType::SURFACE_TEMPERATURE;
        } else {
            pp.query("most.surf_temp_flux", surf_temp_flux);
            theta_type = ThetaCalcType::HEAT_FLUX;
        }

        int nlevs = m_geom.size();
        z_0.resize(nlevs);
        u_star.resize(nlevs);
        t_star.resize(nlevs);
        t_surf.resize(nlevs);
        olen.resize(nlevs);

        for (int lev = 0; lev < nlevs; lev++)
        {
            // Z0 heights
            //--------------------------------------------------------
            // TODO: generalize the "3" for the number of ghost cells
            amrex::Box bx = amrex::grow(m_geom[lev].Domain(),3);
            bx.setSmall(2,0);
            bx.setBig(2,0);
            z_0[lev].resize(bx,1);
            z_0[lev].setVal<amrex::RunOn::Device>(z0_const);

            // 2D MFs for U*, T*, T_surf
            //--------------------------------------------------------
            { // CC vars
                auto& mf = vars_old[lev][Vars::cons];
                // Create a 2D ba, dm, & ghost cells
                amrex::BoxArray ba  = mf.boxArray();
                amrex::BoxList bl2d = ba.boxList();
                for (auto& b : bl2d) {
                    b.setRange(2,0);
                }
                amrex::BoxArray ba2d(std::move(bl2d));
                const amrex::DistributionMapping& dm = mf.DistributionMap();
                const int ncomp   = 1;
                amrex::IntVect ng = mf.nGrowVect(); ng[2]=0;

                u_star[lev] = new amrex::MultiFab(ba2d,dm,ncomp,ng);
                u_star[lev]->setVal(1.E34);

                t_star[lev] = new amrex::MultiFab(ba2d,dm,ncomp,ng);
                t_star[lev]->setVal(1.E34);

                olen[lev] = new amrex::MultiFab(ba2d,dm,ncomp,ng);
                olen[lev]->setVal(1.E34);

                t_surf[lev] = new amrex::MultiFab(ba2d,dm,ncomp,ng);
                if (erf_st) {
                   t_surf[lev]->setVal(surf_temp);
                } else {
                    t_surf[lev]->setVal(0.0);
                }
            }
        }// lev
    }

    // Destructor
    ~ABLMost()
    {
        for (int lev(0); lev<m_geom.size(); ++lev)
        {
          delete u_star[lev];
          delete t_star[lev];
          delete olen[lev];
          delete t_surf[lev];
        }
    }

    ABLMostData
    get_most_data ()
    {
        return *this;
    }

    void
    impose_most_bcs(const int lev,
                    const amrex::Vector<amrex::MultiFab*>& mfs,
                    amrex::MultiFab* eddyDiffs);

    void
    update_fluxes(int lev, int max_iters = 25);

    void
    update_mac_ptrs(int lev,
                    amrex::Vector<amrex::Vector<amrex::MultiFab>>& vars_old,
                    amrex::Vector<std::unique_ptr<amrex::MultiFab>>& Theta_prim)
    { m_ma.update_field_ptrs(lev,vars_old,Theta_prim); }

    const amrex::MultiFab*
    get_u_star(int lev) { return u_star[lev]; }

    const amrex::MultiFab*
    get_t_star(int lev) { return t_star[lev]; }

    const amrex::MultiFab*
    get_olen(int lev) { return olen[lev]; }

    const amrex::MultiFab*
    get_mac_avg(int lev, int comp) { return m_ma.get_average(lev,comp); }

    enum MomCalcType {
        DEFAULT = 0,
        CHARNOCK,
        MODIFIED_CHARNOCK
    };

    enum ThetaCalcType {
        HEAT_FLUX = 0,      ///< Heat-flux specified
        SURFACE_TEMPERATURE ///< Surface temperature specified
    };

    MomCalcType   mom_type;
    ThetaCalcType theta_type;

    void print() const
    {
        amrex::Print() << "ABLMost:\n";
        amrex::Print() << " kappa: "          << kappa          << "\n";
        amrex::Print() << " gravity: "        << gravity        << "\n";
        amrex::Print() << " surf_temp_flux: " << surf_temp_flux << "\n";
        amrex::Print() << " surf_temp: "      << surf_temp      << "\n";
        amrex::Print() << " beta_m: "         << beta_m         << "\n";
        amrex::Print() << " beta_h: "         << beta_h         << "\n";
        amrex::Print() << " gamma_m: "        << gamma_m        << "\n";
        amrex::Print() << " gamma_h: "        << gamma_h        << "\n";
    }

    private:
        amrex::Vector<amrex::Geometry>  m_geom;
        amrex::Vector<amrex::FArrayBox> z_0;

        MOSTAverage m_ma;
        amrex::Vector<amrex::MultiFab*> u_star;
        amrex::Vector<amrex::MultiFab*> t_star;
        amrex::Vector<amrex::MultiFab*> olen;
        amrex::Vector<amrex::MultiFab*> t_surf;
};

#endif /* ABLMOST_H */
