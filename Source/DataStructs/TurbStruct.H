#ifndef _TURB_STRUCT_H_
#define _TURB_STRUCT_H_

#include <MYNNStruct.H>

enum struct LESType {
    None, Smagorinsky, Deardorff
};

enum struct PBLType {
    None, MYNN25, YSU
};


/**
 * Container holding quantities related to turbulence parametrizations
 */
struct TurbChoice {
  public:
    void init_params(int lev, int max_level)
    {
        amrex::ParmParse pp(pp_prefix);

        int nvals_les = pp.countval("les_type");
        int nvals_pbl = pp.countval("pbl_type");

        // If nvals_les and nvals_pbl are both zero then we skip everything below
        if ( (nvals_les == 0) && (nvals_pbl == 0) ) {
            les_type = LESType::None;
            pbl_type = PBLType::None;
            return;
        }

        // If nvals_les and nvals_pbl are both > 0 then they must be the same, and either 1 or max_levels
        if ( (nvals_les > 0) && (nvals_pbl > 0) ) {
            if (nvals_les != nvals_pbl) {
                amrex::Error("If specifying both, we must specify the same number of values for les_type and pbl_type");
            }
            if (nvals_les != 1 && nvals_les != max_level+1) {
                amrex::Error("If specifying both, we must either input one value for all levels, or one value per level");
            }
        }

        // Here we cover the case where either one is 1 and the other is 0, or they both are = 1
        if (nvals_les == 1 || nvals_pbl == 1) {

            // Which LES closure?
            std::string les_type_string = "None";
            pp.query("les_type",les_type_string);

            if (!les_type_string.compare("Smagorinsky")) {
                les_type = LESType::Smagorinsky;
            } else if (!les_type_string.compare("Deardorff")) {
                les_type = LESType::Deardorff;
            } else if (!les_type_string.compare("None")) {
                les_type = LESType::None; // Means DNS
            } else {
                amrex::Error("Don't know this les_type");
            }

            // Which PBL Closure
            static std::string pbl_type_string = "None";
            pp.query("pbl_type",pbl_type_string);
            if (pbl_type_string == "MYNN2.5") {
                pbl_type = PBLType::MYNN25;
            } else if (pbl_type_string == "YSU") {
                pbl_type = PBLType::YSU;
            } else if (pbl_type_string == "None") {
                pbl_type = PBLType::None;
            } else {
                amrex::Error("Don't know this pbl_type");
            }

            // Do some more stuff for PBL Modeling
            if (pbl_type != PBLType::None) {
                // Check for compatibility between PBL, LES, Molec Transport
                if (les_type == LESType::Deardorff) {
                    amrex::Error("It is not recommended to use Deardorff LES and a PBL model");
                } else if (les_type != LESType::None) {
                    amrex::Print() << "Selected a PBL model and an LES model: " <<
                    "Using PBL for vertical transport, LES for horizontal" << std::endl;
                }

                if (pbl_type == PBLType::MYNN25) {
                    pp.query("pbl_mynn_A1", pbl_mynn_A1);
                    pp.query("pbl_mynn_A2", pbl_mynn_A2);
                    pp.query("pbl_mynn_B1", pbl_mynn_B1);
                    pp.query("pbl_mynn_B2", pbl_mynn_B2);
                    pp.query("pbl_mynn_C1", pbl_mynn_C1);
                    pp.query("pbl_mynn_C2", pbl_mynn_C2);
                    pp.query("pbl_mynn_C3", pbl_mynn_C3);
                    pp.query("pbl_mynn_C4", pbl_mynn_C4);
                    pp.query("pbl_mynn_C5", pbl_mynn_C5);
                    pbl_mynn_level2.init_coeffs(pbl_mynn_A1, pbl_mynn_A2, pbl_mynn_B1, pbl_mynn_B2,
                                                pbl_mynn_C1, pbl_mynn_C2, pbl_mynn_C3, pbl_mynn_C4, pbl_mynn_C5);
                }
            }

            // Right now, solving the QKE equation is only supported when MYNN PBL is turned on
            if (pbl_type == PBLType::MYNN25) {
                use_QKE = true;
                pp.query("diffuse_QKE_3D", diffuse_QKE_3D);
                pp.query("advect_QKE", advect_QKE);
            }

            // LES constants...
            pp.query("Cs" , Cs);
            pp.query("CI" , CI);
            pp.query("Pr_t", Pr_t);
            pp.query("Sc_t", Sc_t);

            // Compute relevant forms of diffusion parameters
            Pr_t_inv = amrex::Real(1.0) / Pr_t;
            Sc_t_inv = amrex::Real(1.0) / Sc_t;

            pp.query("Ce" , Ce);
            pp.query("Ce_wall" , Ce_wall);
            pp.query("sigma_k" , sigma_k);

            if (les_type == LESType::Deardorff) {
                pp.query("Ck" , Ck);
            }

            pp.query("theta_ref", theta_ref);

            // Validate inputs
            if (les_type == LESType::Smagorinsky) {
                if (Cs == 0) {
                    amrex::Error("Need to specify Cs for Smagorsinky LES");
                }
            }

        // Here we cover the case where either one is > 1 and the other is 0, or they both are the same and > 1
        } else {

            // Which LES closure?
            std::string les_type_string = "None";
            pp.get("les_type", les_type_string, lev);

            if (!les_type_string.compare("Smagorinsky")) {
                les_type = LESType::Smagorinsky;
            } else if (!les_type_string.compare("Deardorff")) {
                les_type = LESType::Deardorff;
            } else if (!les_type_string.compare("None")) {
                les_type = LESType::None; // Means DNS
            } else {
                amrex::Error("Don't know this les_type");
            }

            // Which PBL Closure
            static std::string pbl_type_string = "None";
            pp.get("pbl_type", pbl_type_string, lev);
            if (pbl_type_string == "MYNN2.5") {
                pbl_type = PBLType::MYNN25;
            } else if (pbl_type_string == "YSU") {
                pbl_type = PBLType::YSU;
            } else if (pbl_type_string == "None") {
                pbl_type = PBLType::None;
            } else {
                amrex::Error("Don't know this pbl_type");
            }

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
                    pp.query("pbl_mynn_A1", pbl_mynn_A1, lev);
                    pp.query("pbl_mynn_A2", pbl_mynn_A2, lev);
                    pp.query("pbl_mynn_B1", pbl_mynn_B1, lev);
                    pp.query("pbl_mynn_B2", pbl_mynn_B2, lev);
                    pp.query("pbl_mynn_C1", pbl_mynn_C1, lev);
                    pp.query("pbl_mynn_C2", pbl_mynn_C2, lev);
                    pp.query("pbl_mynn_C3", pbl_mynn_C3, lev);
                    pp.query("pbl_mynn_C4", pbl_mynn_C4, lev);
                    pp.query("pbl_mynn_C5", pbl_mynn_C5, lev);
                    pbl_mynn_level2.init_coeffs(pbl_mynn_A1, pbl_mynn_A2, pbl_mynn_B1, pbl_mynn_B2,
                                                pbl_mynn_C1, pbl_mynn_C2, pbl_mynn_C3, pbl_mynn_C4, pbl_mynn_C5);
                }
            }

            // Right now, solving the QKE equation is only supported when MYNN PBL is turned on
            if (pbl_type == PBLType::MYNN25) {
                use_QKE = true;
                pp.query("diffuse_QKE_3D", diffuse_QKE_3D, lev);
                pp.query("advect_QKE"    , advect_QKE, lev);
            }

            // LES constants...
            pp.query("Cs"  ,Cs, lev);
            pp.query("CI"  ,CI, lev);
            pp.query("Pr_t",Pr_t, lev);
            pp.query("Sc_t",Sc_t, lev);

            // Compute relevant forms of diffusion parameters
            Pr_t_inv = amrex::Real(1.0) / Pr_t;
            Sc_t_inv = amrex::Real(1.0) / Sc_t;

            pp.query("Ce"      , Ce, lev);
            pp.query("Ce_wall" , Ce_wall, lev);
            pp.query("sigma_k" , sigma_k, lev);

            if (les_type == LESType::Deardorff) {
                pp.query("Ck" , Ck, lev);
            }

            pp.query("theta_ref", theta_ref, lev);
        }
    }

    void display(int lev)
    {
        amrex::Print() << "Turbulence Settings at level " << lev << std::endl;

        if (les_type == LESType::None && pbl_type == PBLType::None) {
            amrex::Print() << "Using DNS model at level " << lev << std::endl;
        } else if  (les_type == LESType::Smagorinsky) {
            amrex::Print() << "Using Smagorinsky LES model at level " << lev << std::endl;
        } else if  (les_type == LESType::Deardorff) {
            amrex::Print() << "Using Deardorff LES model at level " << lev << std::endl;
        } else if  (pbl_type == PBLType::MYNN25) {
            amrex::Print() << "Using MYNN2.5 PBL model at level " << lev << std::endl;
        } else if  (pbl_type == PBLType::YSU) {
            amrex::Print() << "Using YSU PBL model at level " << lev << std::endl;
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
            amrex::Print() << "pbl_mynn_A1                      : " << pbl_mynn_A1 << std::endl;
            amrex::Print() << "pbl_mynn_A2                      : " << pbl_mynn_A2 << std::endl;
            amrex::Print() << "pbl_mynn_B1                      : " << pbl_mynn_B1 << std::endl;
            amrex::Print() << "pbl_mynn_B2                      : " << pbl_mynn_B2 << std::endl;
            amrex::Print() << "pbl_mynn_C1                      : " << pbl_mynn_C1 << std::endl;
            amrex::Print() << "pbl_mynn_C2                      : " << pbl_mynn_C2 << std::endl;
            amrex::Print() << "pbl_mynn_C3                      : " << pbl_mynn_C3 << std::endl;
            amrex::Print() << "pbl_mynn_C4                      : " << pbl_mynn_C4 << std::endl;
            amrex::Print() << "pbl_mynn_C5                      : " << pbl_mynn_C5 << std::endl;
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

    // Model coefficients - MYNN2.5 (from Nakanishi & Niino 2009 [NN09])
    // TODO: Move to MYNNStruct.H
    amrex::Real pbl_mynn_A1 = 1.18;
    amrex::Real pbl_mynn_A2 = 0.665;
    amrex::Real pbl_mynn_B1 = 24.0;
    amrex::Real pbl_mynn_B2 = 15.0;
    amrex::Real pbl_mynn_C1 = 0.137;
    amrex::Real pbl_mynn_C2 = 0.75;
    amrex::Real pbl_mynn_C3 = 0.352;
    amrex::Real pbl_mynn_C4 = 0.0;
    amrex::Real pbl_mynn_C5 = 0.2;
    MYNNLevel2 pbl_mynn_level2;

    // Model coefficients - YSU
    // TODO: Add parmparse for all of these above
    amrex::Real pbl_ysu_coriolis_freq = 1.0e-4; // TODO: make this consistent with coriolis forcing (but note WRF always uses 1e-4)
    bool        pbl_ysu_over_land = false; // TODO: pull from other inputs and make local
    amrex::Real pbl_ysu_land_Ribcr = 0.25; // Critical Bulk Richardson number of Land for stable conditions
    amrex::Real pbl_ysu_unst_Ribcr = 0.0;  // Critical Bulk Richardson number for unstable conditions

    // QKE stuff - default is not to use it, if MYNN2.5 PBL is used default is turb transport in Z-direction only
    bool use_QKE = false;
    bool diffuse_QKE_3D = false;
    bool advect_QKE = true;
};
#endif
