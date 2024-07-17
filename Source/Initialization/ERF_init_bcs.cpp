#include <AMReX_Vector.H>
#include <AMReX_BC_TYPES.H>
#include <AMReX_ParmParse.H>
#include <ERF.H>
#include <Utils/ParFunctions.H>

using namespace amrex;

/**
 * Initializes data structures in the ERF class that specify
 * which boundary conditions we are implementing on each face
 * of the domain.
 *
 * This function also maps the selected boundary condition types
 * (e.g. Outflow, Inflow, Periodic, Dirichlet, ...) to the
 * specific implementation needed for each variable.
 *
 * Stores this information in both host and device vectors
 * so it is available for GPU kernels.
 */
void ERF::init_bcs ()
{
    auto f = [this] (std::string const& bcid, Orientation ori)
    {
        // These are simply defaults for Dirichlet faces -- they should be over-written below
        m_bc_extdir_vals[BCVars::Rho_bc_comp][ori]       =  1.0;
        m_bc_extdir_vals[BCVars::RhoTheta_bc_comp][ori]  = -1.0; // It is important to set this negative
                                                                 // because the sign is tested on below
        m_bc_extdir_vals[BCVars::RhoKE_bc_comp][ori]     = 0.0;
        m_bc_extdir_vals[BCVars::RhoQKE_bc_comp][ori]    = 0.0;
        m_bc_extdir_vals[BCVars::RhoScalar_bc_comp][ori] = 0.0;
        m_bc_extdir_vals[BCVars::RhoQ1_bc_comp][ori]     = 0.0;
        m_bc_extdir_vals[BCVars::RhoQ2_bc_comp][ori]     = 0.0;

        m_bc_extdir_vals[BCVars::xvel_bc][ori] = 0.0; // default
        m_bc_extdir_vals[BCVars::yvel_bc][ori] = 0.0;
        m_bc_extdir_vals[BCVars::zvel_bc][ori] = 0.0;

        // These are simply defaults for Neumann gradients -- they should be over-written below
        m_bc_neumann_vals[BCVars::Rho_bc_comp][ori]       = 0.0;
        m_bc_neumann_vals[BCVars::RhoTheta_bc_comp][ori]  = 0.0;
        m_bc_neumann_vals[BCVars::RhoKE_bc_comp][ori]     = 0.0;
        m_bc_neumann_vals[BCVars::RhoQKE_bc_comp][ori]    = 0.0;
        m_bc_neumann_vals[BCVars::RhoScalar_bc_comp][ori] = 0.0;
        m_bc_neumann_vals[BCVars::RhoQ1_bc_comp][ori]     = 0.0;
        m_bc_neumann_vals[BCVars::RhoQ2_bc_comp][ori]     = 0.0;
        m_bc_neumann_vals[BCVars::xvel_bc][ori] = 0.0;
        m_bc_neumann_vals[BCVars::yvel_bc][ori] = 0.0;
        m_bc_neumann_vals[BCVars::zvel_bc][ori] = 0.0;

        std::string pp_text;
        if (pp_prefix == "erf") {
          pp_text = bcid;
        } else {
          pp_text = pp_prefix + "." + bcid;
        }
        ParmParse pp(pp_text);
        std::string bc_type_in = "null";
        pp.query("type", bc_type_in);
        //if (pp.query("type", bc_type_in))
        //   Print() << "INPUT BC TYPE " << bcid << " " << bc_type_in << std::endl;
        std::string bc_type = amrex::toLower(bc_type_in);

        if (bc_type == "symmetry")
        {
            // Print() << bcid << " set to symmetry.\n";
              phys_bc_type[ori] = ERF_BC::symmetry;
            domain_bc_type[ori] = "Symmetry";
        }
        else if (bc_type == "outflow")
        {
            // Print() << bcid << " set to outflow.\n";
              phys_bc_type[ori] = ERF_BC::outflow;
            domain_bc_type[ori] = "Outflow";
        }
        else if (bc_type == "open")
        {
            // Print() << bcid << " set to open.\n";
            AMREX_ASSERT_WITH_MESSAGE((ori.coordDir() != 2), "Open boundary not valid on zlo or zhi!");
              phys_bc_type[ori] = ERF_BC::open;
            domain_bc_type[ori] = "Open";
        }
        else if (bc_type == "ho_outflow")
        {
            phys_bc_type[ori] = ERF_BC::ho_outflow;
            domain_bc_type[ori] = "HO_Outflow";
        }

        else if (bc_type == "inflow")
        {
            // Print() << bcid << " set to inflow.\n";
              phys_bc_type[ori] = ERF_BC::inflow;
            domain_bc_type[ori] = "Inflow";

            std::vector<Real> v;
            if (input_bndry_planes && m_r2d->ingested_velocity()) {
                m_bc_extdir_vals[BCVars::xvel_bc][ori] = 0.;
                m_bc_extdir_vals[BCVars::yvel_bc][ori] = 0.;
                m_bc_extdir_vals[BCVars::zvel_bc][ori] = 0.; 
            } else {
                // Test for input data file if at xlo face
                std::string dirichlet_file;
                auto file_exists = pp.query("dirichlet_file", dirichlet_file);
                if (ori.isLow() && (ori.coordDir() == 0) && file_exists) {
                    init_Dirichlet_bc_data(ori, dirichlet_file);
                } else {
                    pp.getarr("velocity", v, 0, AMREX_SPACEDIM);
                    m_bc_extdir_vals[BCVars::xvel_bc][ori] = v[0];
                    m_bc_extdir_vals[BCVars::yvel_bc][ori] = v[1];
                    m_bc_extdir_vals[BCVars::zvel_bc][ori] = v[2];
                }
            }

            Real rho_in;
            if (input_bndry_planes && m_r2d->ingested_density()) {
                m_bc_extdir_vals[BCVars::Rho_bc_comp][ori] = 0.;
            } else {
                pp.get("density", rho_in);
                m_bc_extdir_vals[BCVars::Rho_bc_comp][ori] = rho_in;
            }
            Real theta_in;
            if (input_bndry_planes && m_r2d->ingested_theta()) {
                m_bc_extdir_vals[BCVars::RhoTheta_bc_comp][ori] = 0.;
            } else {
                pp.get("theta", theta_in);
                m_bc_extdir_vals[BCVars::RhoTheta_bc_comp][ori] = rho_in*theta_in;
            }
            Real scalar_in = 0.;
            if (input_bndry_planes && m_r2d->ingested_scalar()) {
                m_bc_extdir_vals[BCVars::RhoScalar_bc_comp][ori] = 0.;
            } else {
                if (pp.query("scalar", scalar_in))
                m_bc_extdir_vals[BCVars::RhoScalar_bc_comp][ori] = rho_in*scalar_in;
            }

            if (solverChoice.moisture_type != MoistureType::None) {
                Real qv_in = 0.;
                if (input_bndry_planes && m_r2d->ingested_q1()) {
                    m_bc_extdir_vals[BCVars::RhoQ1_bc_comp][ori] = 0.;
                } else {
                    if (pp.query("qv", qv_in))
                    m_bc_extdir_vals[BCVars::RhoQ1_bc_comp][ori] = rho_in*qv_in;
                }
                Real qc_in = 0.;
                if (input_bndry_planes && m_r2d->ingested_q2()) {
                    m_bc_extdir_vals[BCVars::RhoQ2_bc_comp][ori] = 0.;
                } else {
                    if (pp.query("qc", qc_in))
                    m_bc_extdir_vals[BCVars::RhoQ2_bc_comp][ori] = rho_in*qc_in;
                }
            }

            Real KE_in = 0.;
            if (input_bndry_planes && m_r2d->ingested_KE()) {
                m_bc_extdir_vals[BCVars::RhoKE_bc_comp][ori] = 0.;
            } else {
                if (pp.query("KE", KE_in))
                m_bc_extdir_vals[BCVars::RhoKE_bc_comp][ori] = rho_in*KE_in;
            }
            Real QKE_in = 0.;
            if (input_bndry_planes && m_r2d->ingested_QKE()) {
                m_bc_extdir_vals[BCVars::RhoQKE_bc_comp][ori] = 0.;
            } else {
                if (pp.query("QKE", QKE_in))
                m_bc_extdir_vals[BCVars::RhoQKE_bc_comp][ori] = rho_in*QKE_in;
            }

        }
        else if (bc_type == "noslipwall")
        {
            // Print() << bcid <<" set to no-slip wall.\n";
              phys_bc_type[ori] = ERF_BC::no_slip_wall;
            domain_bc_type[ori] = "NoSlipWall";

            std::vector<Real> v;

            // The values of m_bc_extdir_vals default to 0.
            // But if we find "velocity" in the inputs file, use those values instead.
            if (pp.queryarr("velocity", v, 0, AMREX_SPACEDIM))
            {
                v[ori.coordDir()] = 0.0;
                m_bc_extdir_vals[BCVars::xvel_bc][ori] = v[0];
                m_bc_extdir_vals[BCVars::yvel_bc][ori] = v[1];
                m_bc_extdir_vals[BCVars::zvel_bc][ori] = v[2];
            }

            Real rho_in;
            if (pp.query("density", rho_in))
            {
                m_bc_extdir_vals[BCVars::Rho_bc_comp][ori] = rho_in;
            }

            Real theta_in;
            if (pp.query("theta", theta_in))
            {
               m_bc_extdir_vals[BCVars::RhoTheta_bc_comp][ori] = theta_in*m_bc_extdir_vals[BCVars::Rho_bc_comp][ori];
            }

            Real theta_grad_in;
            if (pp.query("theta_grad", theta_grad_in))
            {
                m_bc_neumann_vals[BCVars::RhoTheta_bc_comp][ori] = theta_grad_in;
            }
        }
        else if (bc_type == "slipwall")
        {
            // Print() << bcid <<" set to slip wall.\n";

              phys_bc_type[ori] = ERF_BC::slip_wall;
            domain_bc_type[ori] = "SlipWall";

            Real rho_in;
            if (pp.query("density", rho_in))
            {
                m_bc_extdir_vals[BCVars::Rho_bc_comp][ori] = rho_in;
            }

            Real theta_in;
            if (pp.query("theta", theta_in))
            {
               m_bc_extdir_vals[BCVars::RhoTheta_bc_comp][ori] = theta_in*m_bc_extdir_vals[BCVars::Rho_bc_comp][ori];
            }

            Real rho_grad_in;
            if (pp.query("density_grad", rho_grad_in))
            {
               m_bc_neumann_vals[BCVars::Rho_bc_comp][ori] = rho_grad_in;
            }

            Real theta_grad_in;
            if (pp.query("theta_grad", theta_grad_in))
            {
               m_bc_neumann_vals[BCVars::RhoTheta_bc_comp][ori] = theta_grad_in;
            }
        }
        else if (bc_type == "most")
        {
              phys_bc_type[ori] = ERF_BC::MOST;
            domain_bc_type[ori] = "MOST";
        }
        else
        {
            phys_bc_type[ori] = ERF_BC::undefined;
        }

        if (geom[0].isPeriodic(ori.coordDir())) {
            domain_bc_type[ori] = "Periodic";
            if (phys_bc_type[ori] == ERF_BC::undefined)
            {
                phys_bc_type[ori] = ERF_BC::periodic;
            } else {
                Abort("Wrong BC type for periodic boundary");
            }
        }

        if (phys_bc_type[ori] == ERF_BC::undefined)
        {
             Print() << "BC Type specified for face " << bcid << " is " << bc_type_in << std::endl;
             Abort("This BC type is unknown");
        }
    };

    f("xlo", Orientation(Direction::x,Orientation::low));
    f("xhi", Orientation(Direction::x,Orientation::high));
    f("ylo", Orientation(Direction::y,Orientation::low));
    f("yhi", Orientation(Direction::y,Orientation::high));
    f("zlo", Orientation(Direction::z,Orientation::low));
    f("zhi", Orientation(Direction::z,Orientation::high));

    // *****************************************************************************
    //
    // Here we translate the physical boundary conditions -- one type per face --
    //     into logical boundary conditions for each velocity component
    //
    // *****************************************************************************
    {
        domain_bcs_type.resize(AMREX_SPACEDIM+NVAR_max);
        domain_bcs_type_d.resize(AMREX_SPACEDIM+NVAR_max);

        for (OrientationIter oit; oit; ++oit) {
            Orientation ori = oit();
            int dir = ori.coordDir();
            Orientation::Side side = ori.faceDir();
            auto const bct = phys_bc_type[ori];
            if ( bct == ERF_BC::symmetry )
            {
                if (side == Orientation::low) {
                    for (int i = 0; i < AMREX_SPACEDIM; i++) {
                        domain_bcs_type[BCVars::xvel_bc+i].setLo(dir, ERFBCType::reflect_even);
                    }
                    domain_bcs_type[BCVars::xvel_bc+dir].setLo(dir, ERFBCType::reflect_odd);
                } else {
                    for (int i = 0; i < AMREX_SPACEDIM; i++) {
                        domain_bcs_type[BCVars::xvel_bc+i].setHi(dir, ERFBCType::reflect_even);
                    }
                    domain_bcs_type[BCVars::xvel_bc+dir].setHi(dir, ERFBCType::reflect_odd);
                }
            }
            else if (bct == ERF_BC::outflow or bct == ERF_BC::ho_outflow )
            {
                if (side == Orientation::low) {
                    for (int i = 0; i < AMREX_SPACEDIM; i++) {
                        domain_bcs_type[BCVars::xvel_bc+i].setLo(dir, ERFBCType::foextrap);
                    }
                    domain_bcs_type[BCVars::xvel_bc+dir].setLo(dir, ERFBCType::neumann_int);
                } else {
                    for (int i = 0; i < AMREX_SPACEDIM; i++) {
                        domain_bcs_type[BCVars::xvel_bc+i].setHi(dir, ERFBCType::foextrap);
                    }
                    domain_bcs_type[BCVars::xvel_bc+dir].setHi(dir, ERFBCType::neumann_int);
                }
            }
            else if (bct == ERF_BC::open)
            {
                if (side == Orientation::low) {
                    for (int i = 0; i < AMREX_SPACEDIM; i++)
                        domain_bcs_type[BCVars::xvel_bc+i].setLo(dir, ERFBCType::open);
                } else {
                    for (int i = 0; i < AMREX_SPACEDIM; i++)
                        domain_bcs_type[BCVars::xvel_bc+i].setHi(dir, ERFBCType::open);
                }
            }
            else if (bct == ERF_BC::inflow)
            {
                if (side == Orientation::low) {
                    for (int i = 0; i < AMREX_SPACEDIM; i++) {
                        domain_bcs_type[BCVars::xvel_bc+i].setLo(dir, ERFBCType::ext_dir);
                        if (input_bndry_planes && dir < 2 && m_r2d->ingested_velocity()) {
                            domain_bcs_type[BCVars::xvel_bc+i].setLo(dir, ERFBCType::ext_dir_ingested);
                        }
                    }
                } else {
                    for (int i = 0; i < AMREX_SPACEDIM; i++) {
                        domain_bcs_type[BCVars::xvel_bc+i].setHi(dir, ERFBCType::ext_dir);
                        if (input_bndry_planes && dir < 2 && m_r2d->ingested_velocity()) {
                            domain_bcs_type[BCVars::xvel_bc+i].setHi(dir, ERFBCType::ext_dir_ingested);
                        }
                    }
                }
            }
            else if (bct == ERF_BC::no_slip_wall)
            {
                if (side == Orientation::low) {
                    for (int i = 0; i < AMREX_SPACEDIM; i++) {
                        domain_bcs_type[BCVars::xvel_bc+i].setLo(dir, ERFBCType::ext_dir);
                    }
                } else {
                    for (int i = 0; i < AMREX_SPACEDIM; i++) {
                        domain_bcs_type[BCVars::xvel_bc+i].setHi(dir, ERFBCType::ext_dir);
                    }
                }
            }
            else if (bct == ERF_BC::slip_wall)
            {
                if (side == Orientation::low) {
                    for (int i = 0; i < AMREX_SPACEDIM; i++) {
                        domain_bcs_type[BCVars::xvel_bc+i].setLo(dir, ERFBCType::foextrap);
                    }
                    // Only normal direction has ext_dir
                    domain_bcs_type[BCVars::xvel_bc+dir].setLo(dir, ERFBCType::ext_dir);

                } else {
                    for (int i = 0; i < AMREX_SPACEDIM; i++) {
                        domain_bcs_type[BCVars::xvel_bc+i].setHi(dir, ERFBCType::foextrap);
                    }
                    // Only normal direction has ext_dir
                    domain_bcs_type[BCVars::xvel_bc+dir].setHi(dir, ERFBCType::ext_dir);
                }
            }
            else if (bct == ERF_BC::periodic)
            {
                if (side == Orientation::low) {
                    for (int i = 0; i < AMREX_SPACEDIM; i++) {
                        domain_bcs_type[BCVars::xvel_bc+i].setLo(dir, ERFBCType::int_dir);
                    }
                } else {
                    for (int i = 0; i < AMREX_SPACEDIM; i++) {
                        domain_bcs_type[BCVars::xvel_bc+i].setHi(dir, ERFBCType::int_dir);
                    }
                }
            }
            else if ( bct == ERF_BC::MOST )
            {
                AMREX_ALWAYS_ASSERT(dir == 2 && side == Orientation::low);
                domain_bcs_type[BCVars::xvel_bc+0].setLo(dir, ERFBCType::reflect_odd);
                domain_bcs_type[BCVars::xvel_bc+1].setLo(dir, ERFBCType::reflect_odd);
                domain_bcs_type[BCVars::xvel_bc+2].setLo(dir, ERFBCType::ext_dir);
            }
        }
    }

    // *****************************************************************************
    //
    // Here we translate the physical boundary conditions -- one type per face --
    //     into logical boundary conditions for each cell-centered variable
    //
    // *****************************************************************************
    {
        for (OrientationIter oit; oit; ++oit) {
            Orientation ori = oit();
            int dir = ori.coordDir();
            Orientation::Side side = ori.faceDir();
            auto const bct = phys_bc_type[ori];
            if ( bct == ERF_BC::symmetry )
            {
                if (side == Orientation::low) {
                    for (int i = 0; i < NVAR_max; i++) {
                        domain_bcs_type[BCVars::cons_bc+i].setLo(dir, ERFBCType::reflect_even);
                    }
                } else {
                    for (int i = 0; i < NVAR_max; i++) {
                        domain_bcs_type[BCVars::cons_bc+i].setHi(dir, ERFBCType::reflect_even);
                    }
                }
            }
            else if ( bct == ERF_BC::outflow )
            {
                if (side == Orientation::low) {
                    for (int i = 0; i < NVAR_max; i++) {
                        domain_bcs_type[BCVars::cons_bc+i].setLo(dir, ERFBCType::foextrap);
                    }
                } else {
                    for (int i = 0; i < NVAR_max; i++) {
                        domain_bcs_type[BCVars::cons_bc+i].setHi(dir, ERFBCType::foextrap);
                    }
                }
            }
            else if ( bct == ERF_BC::ho_outflow )
            {
                if (side == Orientation::low) {
                    for (int i = 0; i < NVAR_max; i++) {
                        domain_bcs_type[BCVars::cons_bc+i].setLo(dir, ERFBCType::hoextrapcc);
                    }
                } else {
                    for (int i = 0; i < NVAR_max; i++) {
                        domain_bcs_type[BCVars::cons_bc+i].setHi(dir, ERFBCType::hoextrapcc);
                    }
                }
            }
            else if ( bct == ERF_BC::open )
            {
                if (side == Orientation::low) {
                    for (int i = 0; i < NVAR_max; i++)
                        domain_bcs_type[BCVars::cons_bc+i].setLo(dir, ERFBCType::open);
                } else {
                    for (int i = 0; i < NVAR_max; i++)
                        domain_bcs_type[BCVars::cons_bc+i].setHi(dir, ERFBCType::open);
                }
            }
            else if ( bct == ERF_BC::no_slip_wall)
            {
                if (side == Orientation::low) {
                    for (int i = 0; i < NVAR_max; i++) {
                        domain_bcs_type[BCVars::cons_bc+i].setLo(dir, ERFBCType::foextrap);
                    }
                    if (m_bc_extdir_vals[BCVars::RhoTheta_bc_comp][ori] > 0.) {
                        domain_bcs_type[BCVars::RhoTheta_bc_comp].setLo(dir, ERFBCType::ext_dir);
                    }
                    if (std::abs(m_bc_neumann_vals[BCVars::RhoTheta_bc_comp][ori]) > 0.) {
                        domain_bcs_type[BCVars::RhoTheta_bc_comp].setLo(dir, ERFBCType::neumann);
                    }
                } else {
                    for (int i = 0; i < NVAR_max; i++) {
                        domain_bcs_type[BCVars::cons_bc+i].setHi(dir, ERFBCType::foextrap);
                    }
                    if (m_bc_extdir_vals[BCVars::RhoTheta_bc_comp][ori] > 0.) {
                        domain_bcs_type[BCVars::RhoTheta_bc_comp].setHi(dir, ERFBCType::ext_dir);
                    }
                    if (std::abs(m_bc_neumann_vals[BCVars::RhoTheta_bc_comp][ori]) > 0.) {
                        domain_bcs_type[BCVars::RhoTheta_bc_comp].setHi(dir, ERFBCType::neumann);
                    }
                }
            }
            else if (bct == ERF_BC::slip_wall)
            {
                if (side == Orientation::low) {
                    for (int i = 0; i < NVAR_max; i++) {
                        domain_bcs_type[BCVars::cons_bc+i].setLo(dir, ERFBCType::foextrap);
                    }
                    if (m_bc_extdir_vals[BCVars::RhoTheta_bc_comp][ori] > 0.) {
                        domain_bcs_type[BCVars::RhoTheta_bc_comp].setLo(dir, ERFBCType::ext_dir);
                    }
                    if (std::abs(m_bc_neumann_vals[BCVars::RhoTheta_bc_comp][ori]) > 0.) {
                        domain_bcs_type[BCVars::RhoTheta_bc_comp].setLo(dir, ERFBCType::neumann);
                    }
                    if (std::abs(m_bc_neumann_vals[BCVars::Rho_bc_comp][ori]) > 0.) {
                        domain_bcs_type[BCVars::Rho_bc_comp].setLo(dir, ERFBCType::neumann);
                    }
                } else {
                    for (int i = 0; i < NVAR_max; i++) {
                        domain_bcs_type[BCVars::cons_bc+i].setHi(dir, ERFBCType::foextrap);
                    }
                    if (m_bc_extdir_vals[BCVars::RhoTheta_bc_comp][ori] > 0.) {
                        domain_bcs_type[BCVars::RhoTheta_bc_comp].setHi(dir, ERFBCType::ext_dir);
                    }
                    if (std::abs(m_bc_neumann_vals[BCVars::RhoTheta_bc_comp][ori]) > 0.) {
                        domain_bcs_type[BCVars::RhoTheta_bc_comp].setHi(dir, ERFBCType::neumann);
                    }
                    if (std::abs(m_bc_neumann_vals[BCVars::Rho_bc_comp][ori]) > 0.) {
                        domain_bcs_type[BCVars::Rho_bc_comp].setHi(dir, ERFBCType::neumann);
                    }
                }
            }
            else if (bct == ERF_BC::inflow)
            {
                if (side == Orientation::low) {
                    for (int i = 0; i < NVAR_max; i++) {
                        domain_bcs_type[BCVars::cons_bc+i].setLo(dir, ERFBCType::ext_dir);
                        if (input_bndry_planes && dir < 2 && (
                           ( (BCVars::cons_bc+i == BCVars::Rho_bc_comp)       && m_r2d->ingested_density()) ||
                           ( (BCVars::cons_bc+i == BCVars::RhoTheta_bc_comp)  && m_r2d->ingested_theta()  ) ||
                           ( (BCVars::cons_bc+i == BCVars::RhoKE_bc_comp)     && m_r2d->ingested_KE()     ) ||
                           ( (BCVars::cons_bc+i == BCVars::RhoQKE_bc_comp)    && m_r2d->ingested_QKE()    ) ||
                           ( (BCVars::cons_bc+i == BCVars::RhoScalar_bc_comp) && m_r2d->ingested_scalar() ) ||
                           ( (BCVars::cons_bc+i == BCVars::RhoQ1_bc_comp)     && m_r2d->ingested_q1()     ) ||
                           ( (BCVars::cons_bc+i == BCVars::RhoQ2_bc_comp)     && m_r2d->ingested_q2()     )) ) {
                            domain_bcs_type[BCVars::cons_bc+i].setLo(dir, ERFBCType::ext_dir_ingested);
                           }
                    }
                } else {
                    for (int i = 0; i < NVAR_max; i++) {
                        domain_bcs_type[BCVars::cons_bc+i].setHi(dir, ERFBCType::ext_dir);
                        if (input_bndry_planes && dir < 2 && (
                           ( (BCVars::cons_bc+i == BCVars::Rho_bc_comp)       && m_r2d->ingested_density()) ||
                           ( (BCVars::cons_bc+i == BCVars::RhoTheta_bc_comp)  && m_r2d->ingested_theta()  ) ||
                           ( (BCVars::cons_bc+i == BCVars::RhoKE_bc_comp)     && m_r2d->ingested_KE()     ) ||
                           ( (BCVars::cons_bc+i == BCVars::RhoQKE_bc_comp)    && m_r2d->ingested_QKE()    ) ||
                           ( (BCVars::cons_bc+i == BCVars::RhoScalar_bc_comp) && m_r2d->ingested_scalar() ) ||
                           ( (BCVars::cons_bc+i == BCVars::RhoQ1_bc_comp)     && m_r2d->ingested_q1()     ) ||
                           ( (BCVars::cons_bc+i == BCVars::RhoQ2_bc_comp)     && m_r2d->ingested_q2()     )
                           ) ) {
                            domain_bcs_type[BCVars::cons_bc+i].setHi(dir, ERFBCType::ext_dir_ingested);
                           }
                    }
                }
            }
            else if (bct == ERF_BC::periodic)
            {
                if (side == Orientation::low) {
                    for (int i = 0; i < NVAR_max; i++) {
                        domain_bcs_type[BCVars::cons_bc+i].setLo(dir, ERFBCType::int_dir);
                    }
                } else {
                    for (int i = 0; i < NVAR_max; i++) {
                       domain_bcs_type[BCVars::cons_bc+i].setHi(dir, ERFBCType::int_dir);
                    }
                }
            }
            else if ( bct == ERF_BC::MOST )
            {
                AMREX_ALWAYS_ASSERT(dir == 2 && side == Orientation::low);
                for (int i = 0; i < NVAR_max; i++) {
                    domain_bcs_type[BCVars::cons_bc+i].setLo(dir, ERFBCType::foextrap);
                }
            }
        }
    }

    // NOTE: Gpu:copy is a wrapper to htod_memcpy (GPU) or memcpy (CPU) and is a blocking comm
    Gpu::copy(Gpu::hostToDevice, domain_bcs_type.begin(), domain_bcs_type.end(), domain_bcs_type_d.begin());
}

void ERF::init_Dirichlet_bc_data (Orientation ori,
                                  const std::string input_file)
{
    Print() << "START\n";
    // Only on coarsest level!
    int lev = 0;
  
    const int klo = 0;
    const int khi = geom[lev].Domain().bigEnd()[2];
    const int Nz  = geom[lev].Domain().size()[2];
    const Real dz = geom[lev].CellSize()[2];

    const bool use_terrain = (zlevels_stag.size() > 0);
    const Real zbot = (use_terrain) ? zlevels_stag[klo]   : geom[lev].ProbLo(2);
    const Real ztop = (use_terrain) ? zlevels_stag[khi+1] : geom[lev].ProbHi(2);

    Vector<Real> z_inp(Nz+2);
    Vector<Real> u_inp(Nz+2);
    Vector<Real> v_inp(Nz+2);
    Vector<Real> w_inp(Nz+2);

    Vector<Real> z_inp_tmp, u_inp_tmp, v_inp_tmp, w_inp_tmp;

    // Read the dirichlet_input file
    Print() << "dirichlet_input file location : " << input_file << std::endl;
    std::ifstream input_reader(input_file);
    if(!input_reader.is_open()) {
        Error("Error opening the dirichlet_input file.\n");
    } else {
        // Read the contents of the input file
        Print() << "Successfully opened the dirichlet_input file. Now reading... " << std::endl;
        std::string line;

        // First, read the input data into temp vectors; then, interpolate vectors to the
        // domain lo/hi and cell centers (from level 0)

        // Add surface
        z_inp_tmp.push_back(zbot); // height above sea level [m]
        u_inp_tmp.push_back(0.);
        v_inp_tmp.push_back(0.);
        w_inp_tmp.push_back(0.);

        // Read the vertical profile at each given height
        Real z, u, v, w;
        while(std::getline(input_reader, line)) {
            std::istringstream iss_z(line);
            iss_z >> z >> u >> v >> w;
            if (z == zbot) {
                u_inp_tmp[0] = u;
                v_inp_tmp[0] = v;
                w_inp_tmp[0] = w;
            } else {
                AMREX_ALWAYS_ASSERT(z > z_inp_tmp[z_inp_tmp.size()-1]); // sounding is increasing in height
                z_inp_tmp.push_back(z);
                u_inp_tmp.push_back(u);
                v_inp_tmp.push_back(v);
                w_inp_tmp.push_back(w);
                if (z >= ztop) break;
            }
        }

        // At this point, we have an input from zbot up to
        // z_inp_tmp[N-1] >= ztop. Now, interpolate to grid level 0 heights
        const int Ninp = z_inp_tmp.size();
        z_inp[0] = zbot;
        u_inp[0] = u_inp_tmp[0];
        v_inp[0] = v_inp_tmp[0];
        w_inp[0] = w_inp_tmp[0];
        for (int k=0; k < Nz; ++k) {
            z_inp[k+1] = (use_terrain) ? 0.5 * (zlevels_stag[k] + zlevels_stag[k+1])
                                       : zbot + (k + 0.5) * dz;
            u_inp[k+1] = interpolate_1d(z_inp_tmp.dataPtr(), u_inp_tmp.dataPtr(), z_inp[k+1], Ninp);
            v_inp[k+1] = interpolate_1d(z_inp_tmp.dataPtr(), v_inp_tmp.dataPtr(), z_inp[k+1], Ninp);
            w_inp[k+1] = interpolate_1d(z_inp_tmp.dataPtr(), w_inp_tmp.dataPtr(), z_inp[k+1], Ninp);
        }
        z_inp[Nz+1] = ztop;
        u_inp[Nz+1] = interpolate_1d(z_inp_tmp.dataPtr(), u_inp_tmp.dataPtr(), ztop, Ninp);
        v_inp[Nz+1] = interpolate_1d(z_inp_tmp.dataPtr(), v_inp_tmp.dataPtr(), ztop, Ninp);
        w_inp[Nz+1] = interpolate_1d(z_inp_tmp.dataPtr(), w_inp_tmp.dataPtr(), ztop, Ninp);
    }

    amrex::Print() << "Successfully read the dirichlet_input file..." << std::endl;
    input_reader.close();
    
    // Populate host vectors at their respective positions (only on lev==0)
    Vector<Real> zcc(khi+1);
    Vector<Real> znd(khi+2);

    if (z_phys_cc[lev]) {
      // use_terrain=1
      // calculate the damping strength based on the max height at each k
      reduce_to_max_per_level(zcc, z_phys_cc[lev]);
      reduce_to_max_per_level(znd, z_phys_nd[lev]);
    } else {
      const auto *const prob_lo = geom[lev].ProbLo();
      const auto *const dx = geom[lev].CellSize();
      for (int k = 0; k <= khi; k++) {
          zcc[k] = prob_lo[2] + (k+0.5) * dx[2];
          znd[k] = prob_lo[2] + (k    ) * dx[2];
      }
      znd[khi+1] = prob_lo[2] + (khi+1) * dx[2];
    }

    Vector<Real> u_dir_h(Nz  ,0.); xvel_bc_data[ori].resize(Nz  ,0.0);
    Vector<Real> v_dir_h(Nz  ,0.); yvel_bc_data[ori].resize(Nz  ,0.0);
    Vector<Real> w_dir_h(Nz+1,0.); zvel_bc_data[ori].resize(Nz+1,0.0);
    for (int k = 0; k <= khi; k++) {
        u_dir_h[k] = interpolate_1d(z_inp.dataPtr(), u_inp.dataPtr(), zcc[k], Nz  );
        v_dir_h[k] = interpolate_1d(z_inp.dataPtr(), v_inp.dataPtr(), zcc[k], Nz  );
        w_dir_h[k] = interpolate_1d(z_inp.dataPtr(), w_inp.dataPtr(), znd[k], Nz+1);
    }
    w_dir_h[khi+1] = interpolate_1d(z_inp.dataPtr(), w_inp.dataPtr(), znd[khi+1], Nz+1);

    
    // Copy from host version to device version
    Gpu::copy(Gpu::hostToDevice, u_dir_h.begin(), u_dir_h.end(), xvel_bc_data[ori].begin());
    Gpu::copy(Gpu::hostToDevice, v_dir_h.begin(), v_dir_h.end(), yvel_bc_data[ori].begin());
    Gpu::copy(Gpu::hostToDevice, w_dir_h.begin(), w_dir_h.end(), zvel_bc_data[ori].begin());
}
