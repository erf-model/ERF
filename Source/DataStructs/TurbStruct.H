#ifndef _TURB_STRUCT_H_
#define _TURB_STRUCT_H_

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
    void init_params()
    {
        amrex::ParmParse pp(pp_prefix);

        // Which LES closure?
        static std::string les_type_string = "None";
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
            if (les_type != LESType::None) {
                amrex::Print() << "Selected a PBL model and an LES model: " <<
                "Using PBL for vertical transport, LES for horizontal" << std::endl;
            } else if (les_type == LESType::Deardorff) {
                amrex::Error("It is not recommended to use Deardorff LES and a PBL model");
            }
            pp.query("pbl_A1", pbl_A1);
            pp.query("pbl_A2", pbl_A2);
            pp.query("pbl_B1", pbl_B1);
            pp.query("pbl_B2", pbl_B2);
            pp.query("pbl_C1", pbl_C1);
            pp.query("pbl_C2", pbl_C2);
            pp.query("pbl_C3", pbl_C3);
            pp.query("pbl_C4", pbl_C4);
            pp.query("pbl_C5", pbl_C5);
        }

        // Right now, solving the QKE equation is only supported when MYNN PBL is turned on
        if (pbl_type == PBLType::MYNN25 ||
            pbl_type == PBLType::YSU       ) {
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
        Pr_t_inv = 1.0 / Pr_t;
        Sc_t_inv = 1.0 / Sc_t;

        pp.query("Ce" , Ce);
        pp.query("Ce_wall" , Ce_wall);
        pp.query("sigma_k" , sigma_k);

        if (les_type == LESType::Deardorff) {
            pp.query("Ck" , Ck);
        }

        pp.query("theta_ref", theta_ref);
    }

    void display()
    {
        amrex::Print() << "Turbulence Settings         : " << std::endl;

        if (les_type == LESType::None) {
            amrex::Print() << "Using DNS model " << std::endl;
        } else if  (les_type == LESType::Smagorinsky) {
            amrex::Print() << "Using Smagorinsky LES model " << std::endl;
        } else if  (les_type == LESType::Deardorff) {
            amrex::Print() << "Using Deardorff LES model " << std::endl;
        }

        amrex::Print() << "Cs                          : " << Cs << std::endl;
        amrex::Print() << "CI                          : " << CI << std::endl;
        amrex::Print() << "Ce                          : " << Ce << std::endl;
        amrex::Print() << "Ce at wall                  : " << Ce_wall << std::endl;
        amrex::Print() << "Ck                          : " << Ck << std::endl;
        amrex::Print() << "reference theta             : " << theta_ref << std::endl;
        amrex::Print() << "sigma_k                     : " << sigma_k << std::endl;
        amrex::Print() << "Pr_t                        : " << Pr_t << std::endl;
        amrex::Print() << "Sc_t                        : " << Sc_t << std::endl;
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
    amrex::Real Pr_t = 1. / 3.;
    amrex::Real Pr_t_inv = 3.0;
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
    // Model coefficients
    amrex::Real pbl_A1 = 1.18;
    amrex::Real pbl_A2 = 0.665;
    amrex::Real pbl_B1 = 24.0;
    amrex::Real pbl_B2 = 15.0;
    amrex::Real pbl_C1 = 0.137;
    amrex::Real pbl_C2 = 0.75;
    amrex::Real pbl_C3 = 0.352;
    amrex::Real pbl_C4 = 0.0;
    amrex::Real pbl_C5 = 0.2;

    // QKE stuff - default is not to use it, if MYNN2.5 or YSU PBL is used default is turb transport in Z-direction only
    bool use_QKE = false;
    bool diffuse_QKE_3D = false;
    bool advect_QKE = true;
};
#endif
