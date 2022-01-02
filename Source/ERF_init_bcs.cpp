#include <AMReX_Vector.H>
#include <AMReX_BC_TYPES.H>
#include <AMReX_ParmParse.H>

#include <ERF.H>

using namespace amrex;

void ERF::init_bcs ()
{
    auto f = [this] (std::string const& bcid, Orientation ori)
    {
        // These are simply defaults for Dirichlet faces -- they should be over-written below
        m_bc_extdir_vals[BCVars::Rho_bc_comp][ori]       =  1.0;
        m_bc_extdir_vals[BCVars::RhoTheta_bc_comp][ori] = -1.0; // It is important to set this negative
                                               // because the sign is tested on below
        m_bc_extdir_vals[BCVars::RhoKE_bc_comp][ori]     = 0.0;
        m_bc_extdir_vals[BCVars::RhoScalar_bc_comp][ori] = 0.0;

        m_bc_extdir_vals[BCVars::xvel_bc][ori] = 0.0; // default
        m_bc_extdir_vals[BCVars::yvel_bc][ori] = 0.0;
        m_bc_extdir_vals[BCVars::zvel_bc][ori] = 0.0;

        ParmParse pp(bcid);
        std::string bc_type_in = "null";
        pp.query("type", bc_type_in);
        //if (pp.query("type", bc_type_in))
        //   amrex::Print() << "INPUT BC TYPE " << bcid << " " << bc_type_in << std::endl;
        std::string bc_type = amrex::toLower(bc_type_in);

        if (bc_type == "symmetry")
        {
            // amrex::Print() << bcid << " set to symmetry.\n";

           phys_bc_type[ori] = BC::symmetry;
        }
        else if (bc_type == "outflow")
        {
            // amrex::Print() << bcid << " set to outflow.\n";

           phys_bc_type[ori] = BC::outflow;
        }
        else if (bc_type == "inflow")
        {
            // amrex::Print() << bcid << " set to inflow.\n";

            phys_bc_type[ori] = BC::inflow;

            std::vector<Real> v;
            pp.getarr("velocity", v, 0, AMREX_SPACEDIM);
            m_bc_extdir_vals[BCVars::xvel_bc][ori] = v[0];
            m_bc_extdir_vals[BCVars::yvel_bc][ori] = v[1];
            m_bc_extdir_vals[BCVars::zvel_bc][ori] = v[2];

            Real rho_in;
            pp.get("density", rho_in);
            {
               m_bc_extdir_vals[BCVars::Rho_bc_comp][ori] = rho_in;
            }
            Real theta_in;
            pp.get("theta", theta_in);
            {
               m_bc_extdir_vals[BCVars::RhoTheta_bc_comp][ori] = rho_in*theta_in;
            }
            Real KE_in = 0.;
            if (pp.query("KE", KE_in))
            {
               m_bc_extdir_vals[BCVars::RhoKE_bc_comp][ori] = rho_in*KE_in;
            }
            Real scalar_in = 0.;
            if (pp.query("scalar", scalar_in))
            {
               m_bc_extdir_vals[BCVars::RhoScalar_bc_comp][ori] = rho_in*scalar_in;
            }

        }
        else if (bc_type == "noslipwall")
        {
            // amrex::Print() << bcid <<" set to no-slip wall.\n";

            phys_bc_type[ori] = BC::no_slip_wall;

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

            Real theta_in;
            if (pp.query("theta", theta_in))
            {
               m_bc_extdir_vals[BCVars::RhoTheta_bc_comp][ori] = theta_in*m_bc_extdir_vals[BCVars::Rho_bc_comp][ori];
            }
        }
        else if (bc_type == "slipwall")
        {
            // amrex::Print() << bcid <<" set to slip wall.\n";

            phys_bc_type[ori] = BC::slip_wall;

            Real theta_in;
            if (pp.query("theta", theta_in))
            {
               m_bc_extdir_vals[BCVars::RhoTheta_bc_comp][ori] = theta_in*m_bc_extdir_vals[BCVars::Rho_bc_comp][ori];
            }

        }
        else if (bc_type == "most")
        {
            // amrex::Print() << bcid <<" set to MOST.\n";

            phys_bc_type[ori] = BC::MOST;

        }
        else
        {
            phys_bc_type[ori] = BC::undefined;
        }

        if (geom[0].isPeriodic(ori.coordDir())) {
            if (phys_bc_type[ori] == BC::undefined)
            {
                phys_bc_type[ori] = BC::periodic;
            } else {
                amrex::Abort("Wrong BC type for periodic boundary");
            }
        }

        if (phys_bc_type[ori] == BC::undefined)
        {
             // amrex::Print() << "BC Type specified for face " << bcid << " is " << bc_type_in << std::endl;
             amrex::Abort("This type is unknown");
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
        domain_bcs_type.resize(AMREX_SPACEDIM+NVAR);
        domain_bcs_type_d.resize(AMREX_SPACEDIM+NVAR);

        for (OrientationIter oit; oit; ++oit) {
            Orientation ori = oit();
            int dir = ori.coordDir();
            Orientation::Side side = ori.faceDir();
            auto const bct = phys_bc_type[ori];
            if ( bct == BC::symmetry )
            {
                if (side == Orientation::low) {
                    for (int i = 0; i < AMREX_SPACEDIM; i++)
                        domain_bcs_type[BCVars::xvel_bc+i].setLo(dir, ERFBCType::reflect_even);
                    domain_bcs_type[BCVars::xvel_bc+dir].setLo(dir, ERFBCType::reflect_odd);
                } else {
                    for (int i = 0; i < AMREX_SPACEDIM; i++)
                        domain_bcs_type[BCVars::xvel_bc+i].setHi(dir, ERFBCType::reflect_even);
                    domain_bcs_type[BCVars::xvel_bc+dir].setHi(dir, ERFBCType::reflect_odd);
                }
            }
            else if (bct == BC::outflow)
            {
                if (side == Orientation::low) {
                    for (int i = 0; i < AMREX_SPACEDIM; i++)
                        domain_bcs_type[BCVars::xvel_bc+i].setLo(dir, ERFBCType::foextrap);
                } else {
                    for (int i = 0; i < AMREX_SPACEDIM; i++)
                        domain_bcs_type[BCVars::xvel_bc+i].setHi(dir, ERFBCType::foextrap);
                }
            }
            else if (bct == BC::inflow || bct == BC::no_slip_wall)
            {
                if (side == Orientation::low) {
                    for (int i = 0; i < AMREX_SPACEDIM; i++)
                        domain_bcs_type[BCVars::xvel_bc+i].setLo(dir, ERFBCType::ext_dir);
                    if (dir == 2) lo_z_is_dirichlet = true;
                } else {
                    for (int i = 0; i < AMREX_SPACEDIM; i++)
                        domain_bcs_type[BCVars::xvel_bc+i].setHi(dir, ERFBCType::ext_dir);
                    if (dir == 2) hi_z_is_dirichlet = true;
                }
            }
            else if (bct == BC::slip_wall)
            {
                if (side == Orientation::low) {
                    // Tangential directions have hoextrap
                    for (int i = 0; i < AMREX_SPACEDIM; i++)
                        domain_bcs_type[BCVars::xvel_bc+i].setLo(dir, ERFBCType::foextrap);
                    // Only normal direction has ext_dir
                    domain_bcs_type[BCVars::xvel_bc+dir].setLo(dir, ERFBCType::ext_dir);

                } else {
                    // Tangential directions have hoextrap
                    for (int i = 0; i < AMREX_SPACEDIM; i++)
                        domain_bcs_type[BCVars::xvel_bc+i].setHi(dir, ERFBCType::foextrap);
                    // Only normal direction has ext_dir
                    domain_bcs_type[BCVars::xvel_bc+dir].setHi(dir, ERFBCType::ext_dir);
                }
            }
            else if (bct == BC::periodic)
            {
                if (side == Orientation::low) {
                    for (int i = 0; i < AMREX_SPACEDIM; i++)
                        domain_bcs_type[BCVars::xvel_bc+i].setLo(dir, ERFBCType::int_dir);
                } else {
                    for (int i = 0; i < AMREX_SPACEDIM; i++)
                        domain_bcs_type[BCVars::xvel_bc+i].setHi(dir, ERFBCType::int_dir);
                }
            }
            else if ( bct == BC::MOST )
            {
                if (dir == 2 && side == Orientation::low) {
                    for (int i = 0; i < AMREX_SPACEDIM; i++)
                        domain_bcs_type[BCVars::xvel_bc+i].setLo(dir, ERFBCType::MOST);
                } else {
                    amrex::Error("MOST bc can only be applied on low z-face");
                }
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
            if ( bct == BC::symmetry )
            {
                if (side == Orientation::low) {
                    for (int i = 0; i < NVAR; i++)
                        domain_bcs_type[BCVars::cons_bc+i].setLo(dir, ERFBCType::reflect_even);
                } else {
                    for (int i = 0; i < NVAR; i++)
                        domain_bcs_type[BCVars::cons_bc+i].setHi(dir, ERFBCType::reflect_even);
                }
            }
            else if ( bct == BC::outflow )
            {
                if (side == Orientation::low) {
                    for (int i = 0; i < NVAR; i++)
                        domain_bcs_type[BCVars::cons_bc+i].setLo(dir, ERFBCType::foextrap);
                } else {
                    for (int i = 0; i < NVAR; i++)
                        domain_bcs_type[BCVars::cons_bc+i].setHi(dir, ERFBCType::foextrap);
                }
            }
            else if ( bct == BC::no_slip_wall)
            {
                if (side == Orientation::low) {
                    for (int i = 0; i < NVAR; i++)
                        domain_bcs_type[BCVars::cons_bc+i].setLo(dir, ERFBCType::foextrap);
                    if (m_bc_extdir_vals[BCVars::RhoTheta_bc_comp][ori] > 0.)
                        domain_bcs_type[BCVars::RhoTheta_bc_comp].setLo(dir, ERFBCType::ext_dir);
                } else {
                    for (int i = 0; i < NVAR; i++)
                        domain_bcs_type[BCVars::cons_bc+i].setHi(dir, ERFBCType::foextrap);
                    if (m_bc_extdir_vals[BCVars::RhoTheta_bc_comp][ori] > 0.)
                        domain_bcs_type[BCVars::RhoTheta_bc_comp].setHi(dir, ERFBCType::ext_dir);
                }
            }
            else if (bct == BC::slip_wall)
            {
                if (side == Orientation::low) {
                    for (int i = 0; i < NVAR; i++)
                        domain_bcs_type[BCVars::cons_bc+i].setLo(dir, ERFBCType::foextrap);
                    if (m_bc_extdir_vals[BCVars::RhoTheta_bc_comp][ori] > 0.)
                        domain_bcs_type[BCVars::RhoTheta_bc_comp].setLo(dir, ERFBCType::ext_dir);
                } else {
                    for (int i = 0; i < NVAR; i++)
                        domain_bcs_type[BCVars::cons_bc+i].setHi(dir, ERFBCType::foextrap);
                    if (m_bc_extdir_vals[BCVars::RhoTheta_bc_comp][ori] > 0.)
                        domain_bcs_type[BCVars::RhoTheta_bc_comp].setHi(dir, ERFBCType::ext_dir);
                }
            }
            else if (bct == BC::inflow)
            {
                if (side == Orientation::low) {
                    for (int i = 0; i < NVAR; i++)
                        domain_bcs_type[BCVars::cons_bc+i].setLo(dir, ERFBCType::ext_dir);
                } else {
                    for (int i = 0; i < NVAR; i++)
                        domain_bcs_type[BCVars::cons_bc+i].setHi(dir, ERFBCType::ext_dir);
                }
            }
            else if (bct == BC::periodic)
            {
                if (side == Orientation::low) {
                    for (int i = 0; i < NVAR; i++)
                        domain_bcs_type[BCVars::cons_bc+i].setLo(dir, ERFBCType::int_dir);
                } else {
                    for (int i = 0; i < NVAR; i++)
                       domain_bcs_type[BCVars::cons_bc+i].setHi(dir, ERFBCType::int_dir);
                }
            }
            else if ( bct == BC::MOST )
            {
                if (dir == 2 && side == Orientation::low) {
                    for (int i = 0; i < NVAR; i++)
                        domain_bcs_type[BCVars::cons_bc+i].setLo(dir, ERFBCType::MOST);
                } else {
                    amrex::Error("MOST bc can only be applied on low z-face");
                }
            }
        }
    }

#ifdef AMREX_USE_GPU
    Gpu::htod_memcpy
        (domain_bcs_type_d.data(), domain_bcs_type.data(),
         sizeof(amrex::BCRec)*(NVAR+AMREX_SPACEDIM));
#else
    std::memcpy
        (domain_bcs_type_d.data(), domain_bcs_type.data(),
         sizeof(amrex::BCRec)*(NVAR+AMREX_SPACEDIM));
#endif
}

