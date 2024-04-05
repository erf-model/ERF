#include <AMReX_Vector.H>
#include <AMReX_BC_TYPES.H>
#include <AMReX_ParmParse.H>

#include <ERF.H>

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
                pp.getarr("velocity", v, 0, AMREX_SPACEDIM);
                m_bc_extdir_vals[BCVars::xvel_bc][ori] = v[0];
                m_bc_extdir_vals[BCVars::yvel_bc][ori] = v[1];
                m_bc_extdir_vals[BCVars::zvel_bc][ori] = v[2];
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
                Real qt_in = 0.;
                if (input_bndry_planes && m_r2d->ingested_q1()) {
                    m_bc_extdir_vals[BCVars::RhoQ1_bc_comp][ori] = 0.;
                } else {
                    if (pp.query("qt", qt_in))
                    m_bc_extdir_vals[BCVars::RhoQ1_bc_comp][ori] = rho_in*qt_in;
                }
                Real qp_in = 0.;
                if (input_bndry_planes && m_r2d->ingested_q2()) {
                    m_bc_extdir_vals[BCVars::RhoQ2_bc_comp][ori] = 0.;
                } else {
                    if (pp.query("qp", qp_in))
                    m_bc_extdir_vals[BCVars::RhoQ2_bc_comp][ori] = rho_in*qp_in;
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
            else if (bct == ERF_BC::outflow)
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

