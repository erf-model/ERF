#ifndef ERF_TURB_STRUCT_H_
#define ERF_TURB_STRUCT_H_

#include <ERF_MYNNStruct.H>

AMREX_ENUM(LESType, None, Smagorinsky, Deardorff);

AMREX_ENUM(PBLType, None, MYNN25, YSU);

template <typename T>
void query_one_or_per_level(const amrex::ParmParse& pp, const char* query_string, T& query_var, const int lev, const int maxlev)
{
    int count = pp.countval(query_string);
    if (count == 0) {
        return; // nothing to do
    } else if (count == 1) {
        pp.query(query_string, query_var);
    } else if (count == maxlev + 1) {
        pp.query(query_string, query_var, lev);
    } else {
        amrex::Error("For parmparse variable " + pp.prefixedName(query_string)
                     + ": if specified, specify once total or once for each level");
    }
}

/**
 * Container holding quantities related to turbulence parametrizations
 */
struct TurbChoice {
  public:
    void init_params(int lev, int max_level)
    {
        amrex::ParmParse pp(pp_prefix);

        // Which LES closure?
        std::string les_type_string = "None";
        query_one_or_per_level(pp, "les_type", les_type, lev, max_level);

        // Which PBL Closure
        static std::string pbl_type_string = "None";
        query_one_or_per_level(pp, "pbl_type", pbl_type, lev, max_level);

        // Do some more stuff for PBL Modeling
        if (pbl_type != PBLType::None) {
            // Check for compatibility between PBL, LES, Molec Transport
            if (les_type != LESType::None) {
                amrex::Print() << "Selected a PBL model and an LES model: " <<
                    "Using PBL for vertical transport, LES for horizontal" << std::endl;
            } else if (les_type == LESType::Deardorff) {
                amrex::Error("It is not recommended to use Deardorff LES and a PBL model");
            }

            if (pbl_type == PBLType::MYNN25) {
                query_one_or_per_level(pp, "pbl_mynn_A1", pbl_mynn.A1, lev, max_level);
                query_one_or_per_level(pp, "pbl_mynn_A2", pbl_mynn.A2, lev, max_level);
                query_one_or_per_level(pp, "pbl_mynn_B1", pbl_mynn.B1, lev, max_level);
                query_one_or_per_level(pp, "pbl_mynn_B2", pbl_mynn.B2, lev, max_level);
                query_one_or_per_level(pp, "pbl_mynn_C1", pbl_mynn.C1, lev, max_level);
                query_one_or_per_level(pp, "pbl_mynn_C2", pbl_mynn.C2, lev, max_level);
                query_one_or_per_level(pp, "pbl_mynn_C3", pbl_mynn.C3, lev, max_level);
                query_one_or_per_level(pp, "pbl_mynn_C4", pbl_mynn.C4, lev, max_level);
                query_one_or_per_level(pp, "pbl_mynn_C5", pbl_mynn.C5, lev, max_level);
                pbl_mynn_level2.init_coeffs(pbl_mynn.A1, pbl_mynn.A2, pbl_mynn.B1, pbl_mynn.B2,
                                            pbl_mynn.C1, pbl_mynn.C2, pbl_mynn.C3, pbl_mynn.C4, pbl_mynn.C5);
                query_one_or_per_level(pp, "pbl_mynn_diffuse_moistvars", pbl_mynn.diffuse_moistvars, lev, max_level);
                query_one_or_per_level(pp, "pbl_mynn_SMmin", pbl_mynn.SMmin, lev, max_level);
                query_one_or_per_level(pp, "pbl_mynn_SMmax", pbl_mynn.SMmax, lev, max_level);
                query_one_or_per_level(pp, "pbl_mynn_SHmin", pbl_mynn.SHmin, lev, max_level);
                query_one_or_per_level(pp, "pbl_mynn_SHmax", pbl_mynn.SHmax, lev, max_level);
            } else if (pbl_type == PBLType::YSU) {
                query_one_or_per_level(pp, "pbl_ysu_coriolis_freq", pbl_ysu_coriolis_freq, lev, max_level);
                query_one_or_per_level(pp, "pbl_ysu_use_consistent_coriolis", pbl_ysu_use_consistent_coriolis, lev, max_level);
                query_one_or_per_level(pp, "pbl_ysu_force_over_water", pbl_ysu_force_over_water, lev, max_level);
                query_one_or_per_level(pp, "pbl_ysu_land_Ribcr", pbl_ysu_land_Ribcr, lev, max_level);
                query_one_or_per_level(pp, "pbl_ysu_unst_Ribcr", pbl_ysu_unst_Ribcr, lev, max_level);
            }
        }

        // Right now, solving the QKE equation is only supported when MYNN PBL is turned on
        if (pbl_type == PBLType::MYNN25) {
            use_QKE = true;
            query_one_or_per_level(pp, "advect_QKE"    , advect_QKE, lev, max_level);
            query_one_or_per_level(pp, "diffuse_QKE_3D", diffuse_QKE_3D, lev, max_level);
        }

        // LES constants...
        query_one_or_per_level(pp, "Cs"  ,Cs, lev, max_level);
        query_one_or_per_level(pp, "CI"  ,CI, lev, max_level);
        query_one_or_per_level(pp, "Pr_t",Pr_t, lev, max_level);
        query_one_or_per_level(pp, "Sc_t",Sc_t, lev, max_level);

        // Compute relevant forms of diffusion parameters
        Pr_t_inv = amrex::Real(1.0) / Pr_t;
        Sc_t_inv = amrex::Real(1.0) / Sc_t;

        query_one_or_per_level(pp, "Ce"      , Ce, lev, max_level);
        query_one_or_per_level(pp, "Ce_wall" , Ce_wall, lev, max_level);
        query_one_or_per_level(pp, "sigma_k" , sigma_k, lev, max_level);

        if (les_type == LESType::Deardorff) {
            query_one_or_per_level(pp, "Ck" , Ck, lev, max_level);
        }

        query_one_or_per_level(pp, "theta_ref", theta_ref, lev, max_level);

        // Validate inputs
        if (les_type == LESType::Smagorinsky) {
            if (Cs == 0) {
                amrex::Error("Need to specify Cs for Smagorsinky LES");
            }
        }
    }

    void display(int lev)
    {
        amrex::Print() << "Turbulence Settings at level " << lev << std::endl;

        if (les_type == LESType::None && pbl_type == PBLType::None) {
            amrex::Print() << "    Using DNS model at level " << lev << std::endl;
        } else if  (les_type == LESType::Smagorinsky) {
            amrex::Print() << "    Using Smagorinsky LES model at level " << lev << std::endl;
        } else if  (les_type == LESType::Deardorff) {
            amrex::Print() << "    Using Deardorff LES model at level " << lev << std::endl;
        } else if  (pbl_type == PBLType::MYNN25) {
            amrex::Print() << "    Using MYNN2.5 PBL model at level " << lev << std::endl;
        } else if  (pbl_type == PBLType::YSU) {
            amrex::Print() << "    Using YSU PBL model at level " << lev << std::endl;
        }

        if (les_type != LESType::None) {
            if (les_type == LESType::Smagorinsky) {
                amrex::Print() << "Cs                          : " << Cs << std::endl;
                amrex::Print() << "CI                          : " << CI << std::endl;
                amrex::Print() << "Pr_t                        : " << Pr_t << std::endl;
                amrex::Print() << "Sc_t                        : " << Sc_t << std::endl;
            }
            if (les_type == LESType::Deardorff) {
                amrex::Print() << "Ce                          : " << Ce << std::endl;
                amrex::Print() << "Ce at wall                  : " << Ce_wall << std::endl;
                amrex::Print() << "Ck                          : " << Ck << std::endl;
                amrex::Print() << "sigma_k                     : " << sigma_k << std::endl;
            }
            amrex::Print() << "reference theta             : " << theta_ref << std::endl;
        }

        if (pbl_type == PBLType::MYNN25) {
            amrex::Print() << "pbl_mynn_A1                      : " << pbl_mynn.A1 << std::endl;
            amrex::Print() << "pbl_mynn_A2                      : " << pbl_mynn.A2 << std::endl;
            amrex::Print() << "pbl_mynn_B1                      : " << pbl_mynn.B1 << std::endl;
            amrex::Print() << "pbl_mynn_B2                      : " << pbl_mynn.B2 << std::endl;
            amrex::Print() << "pbl_mynn_C1                      : " << pbl_mynn.C1 << std::endl;
            amrex::Print() << "pbl_mynn_C2                      : " << pbl_mynn.C2 << std::endl;
            amrex::Print() << "pbl_mynn_C3                      : " << pbl_mynn.C3 << std::endl;
            amrex::Print() << "pbl_mynn_C4                      : " << pbl_mynn.C4 << std::endl;
            amrex::Print() << "pbl_mynn_C5                      : " << pbl_mynn.C5 << std::endl;
        } else if (pbl_type == PBLType::YSU) {
            amrex::Print() << "pbl_ysu_coriolis_freq            : " << pbl_ysu_coriolis_freq << std::endl;
            amrex::Print() << "pbl_ysu_use_consistent_coriolis  : " << pbl_ysu_use_consistent_coriolis << std::endl;
            amrex::Print() << "pbl_ysu_force_over_water         : " << pbl_ysu_force_over_water << std::endl;
            amrex::Print() << "pbl_ysu_land_Ribcr               : " << pbl_ysu_land_Ribcr << std::endl;
            amrex::Print() << "pbl_ysu_unst_Ribcr               : " << pbl_ysu_unst_Ribcr << std::endl;
        }
    }

    // Default prefix
    std::string pp_prefix {"erf"};

    // LES model
    LESType les_type;
    // Smagorinsky Cs coefficient
    amrex::Real Cs = 0.0;
    // Smagorinsky CI coefficient
    amrex::Real CI = 0.0;
    // Smagorinsky Turbulent Prandtl Number
    amrex::Real Pr_t = amrex::Real(1.0) / amrex::Real(3.0);
    amrex::Real Pr_t_inv = amrex::Real(3.0);
    // Smagorinsky Turbulent Schmidt Number
    amrex::Real Sc_t = 1.0;
    amrex::Real Sc_t_inv = 1.0;

    // Deardorff Ce coefficient
    amrex::Real Ce = 0.93;
    amrex::Real Ce_wall = 0.0; // if > 0, then set Ce to this at k=0
    // Deardorff Ck coefficient
    amrex::Real Ck = 0.1;
    // Deardorff sigma_k coefficient
    amrex::Real sigma_k = 0.5;

    amrex::Real theta_ref = 300.0;

    // PBL model
    PBLType pbl_type;

    MYNNLevel25 pbl_mynn;
    MYNNLevel2 pbl_mynn_level2; // for filtering

    // Model coefficients - YSU
    // TODO: Add parmparse for all of these above
    amrex::Real pbl_ysu_coriolis_freq = 1.0e-4; // 1e-4 is hardcoded in WRF, we let the user specify or tske the value from ERF coriolis forcing
    bool        pbl_ysu_use_consistent_coriolis = false; // ignore input pbl_ysu_coriolis_freq, take value from ERF coriolis forcing instead
    bool        pbl_ysu_force_over_water = false; // Force YSU to act as if it is over water regardless of other inputs (for testing)
    amrex::Real pbl_ysu_land_Ribcr = 0.25; // Critical Bulk Richardson number of Land for stable conditions
    amrex::Real pbl_ysu_unst_Ribcr = 0.0;  // Critical Bulk Richardson number for unstable conditions

    // QKE stuff - default is not to use it, if MYNN2.5 PBL is used default is turb transport in Z-direction only
    bool use_QKE = false;
    bool diffuse_QKE_3D = false;
    bool advect_QKE = true;
};
#endif
