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
        m_bc_cons[ori][Rho_comp]       =  1.0;
        m_bc_cons[ori][RhoTheta_comp]  = -1.0; // It is important to set this negative
                                               // because the sign is tested on below
        m_bc_cons[ori][RhoKE_comp]     = 0.0;
        m_bc_cons[ori][RhoScalar_comp] = 0.0;
        AMREX_D_TERM(m_bc_vels[ori][0] = 0.0;, // default
                     m_bc_vels[ori][1] = 0.0;,
                     m_bc_vels[ori][2] = 0.0;);

        ParmParse pp(bcid);
        std::string bc_type_in = "null";
        pp.query("type", bc_type_in);
        if (pp.query("type", bc_type_in))
           amrex::Print() << "INPUT BC TYPE " << bcid << " " << bc_type_in << std::endl;
        std::string bc_type = amrex::toLower(bc_type_in);

        if (bc_type == "symmetry")
        {
            amrex::Print() << bcid << " set to symmetry.\n";

            m_bc_type[ori] = BC::symmetry;
        }
        else if (bc_type == "outflow")
        {
            amrex::Print() << bcid << " set to outflow.\n";

            m_bc_type[ori] = BC::outflow;
        }
        else if (bc_type == "inflow")
        {
            amrex::Print() << bcid << " set to inflow.\n";

            m_bc_type[ori] = BC::inflow;

            std::vector<Real> v;
            pp.getarr("velocity", v, 0, AMREX_SPACEDIM);
            for (int i=0; i<AMREX_SPACEDIM; i++)
            {
                m_bc_vels[ori][i] = v[i];
            }

            Real rho_in;
            pp.get("density", rho_in);
            {
               m_bc_cons[ori][Rho_comp] = rho_in;
            }
            Real theta_in;
            pp.get("theta", theta_in);
            {
               m_bc_cons[ori][RhoTheta_comp] = rho_in*theta_in;
            }
            Real KE_in = 0.;
            if (pp.query("KE", KE_in))
            {
               m_bc_cons[ori][RhoKE_comp] = rho_in*KE_in;
            }
            Real scalar_in = 0.;
            if (pp.query("scalar", scalar_in))
            {
               m_bc_cons[ori][RhoScalar_comp] = rho_in*scalar_in;
            }

        }
        else if (bc_type == "noslipwall")
        {
            amrex::Print() << bcid <<" set to no-slip wall.\n";

            m_bc_type[ori] = BC::no_slip_wall;

            // Note that m_bc_vels defaults to 0 above so we are ok if
            //      queryarr finds nothing
            // Here we make sure that we only use the tangential components
            //      of a specified velocity field -- the wall is not allowed
            //      to move in the normal direction
            std::vector<Real> v;
            if (pp.queryarr("velocity", v, 0, AMREX_SPACEDIM)) {
                v[ori.coordDir()] = 0.0;
                for (int i=0; i<AMREX_SPACEDIM; i++){
                    m_bc_vels[ori][i] = v[i];
                }
            }

            Real theta_in;
            if (pp.query("theta", theta_in))
            {
               m_bc_cons[ori][RhoTheta_comp] = theta_in*m_bc_cons[ori][Rho_comp];
            }
        }
        else if (bc_type == "slipwall")
        {
            amrex::Print() << bcid <<" set to slip wall.\n";

            m_bc_type[ori] = BC::slip_wall;

            Real theta_in;
            if (pp.query("theta", theta_in))
            {
               m_bc_cons[ori][RhoTheta_comp] = theta_in*m_bc_cons[ori][Rho_comp];
            }

        }
        else if (bc_type == "most")
        {
            amrex::Print() << bcid <<" set to MOST.\n";

            m_bc_type[ori] = BC::MOST;

        }
        else
        {
            m_bc_type[ori] = BC::undefined;
        }

        if (geom[0].isPeriodic(ori.coordDir())) {
            if (m_bc_type[ori] == BC::undefined)
            {
                m_bc_type[ori] = BC::periodic;
            } else {
                amrex::Abort("Wrong BC type for periodic boundary");
            }
        }

        if (m_bc_type[ori] == BC::undefined)
        {
             amrex::Print() << "BC Type specified for face " << bcid << " is " << bc_type_in << std::endl;
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
        int nvar_cc  = Cons::NumVars;
        int nvar_vel = AMREX_SPACEDIM;

        bcs.resize(2);
        bcs[BCVars::cons].resize(nvar_cc);
        bcs[BCVars::vels].resize(nvar_vel);

        bcs_d.resize(2);
        bcs_d[BCVars::cons].resize(nvar_cc);
        bcs_d[BCVars::vels].resize(nvar_vel);

        for (OrientationIter oit; oit; ++oit) {
            Orientation ori = oit();
            int dir = ori.coordDir();
            Orientation::Side side = ori.faceDir();
            auto const bct = m_bc_type[ori];
            if ( bct == BC::symmetry )
            {
                if (side == Orientation::low) {
                    for (int i = 0; i < AMREX_SPACEDIM; i++)
                        bcs[BCVars::vels][i].setLo(dir, ERFBCType::reflect_even);
                    bcs[BCVars::vels][dir].setLo(dir, ERFBCType::reflect_odd);
                } else {
                    for (int i = 0; i < AMREX_SPACEDIM; i++)
                        bcs[BCVars::vels][i].setHi(dir, ERFBCType::reflect_even);
                    bcs[BCVars::vels][dir].setHi(dir, ERFBCType::reflect_odd);
                }
            }
            else if (bct == BC::outflow)
            {
                if (side == Orientation::low) {
                    for (int i = 0; i < AMREX_SPACEDIM; i++)
                        bcs[BCVars::vels][i].setLo(dir, ERFBCType::foextrap);
                } else {
                    for (int i = 0; i < AMREX_SPACEDIM; i++)
                        bcs[BCVars::vels][i].setHi(dir, ERFBCType::foextrap);
                }
            }
            else if (bct == BC::inflow || bct == BC::no_slip_wall)
            {
                if (side == Orientation::low) {
                    for (int i = 0; i < AMREX_SPACEDIM; i++)
                        bcs[BCVars::vels][i].setLo(dir, ERFBCType::ext_dir);
                    if (dir == 2) lo_z_is_dirichlet = true;
                } else {
                    for (int i = 0; i < AMREX_SPACEDIM; i++)
                        bcs[BCVars::vels][i].setHi(dir, ERFBCType::ext_dir);
                    if (dir == 2) hi_z_is_dirichlet = true;
                }
            }
            else if (bct == BC::slip_wall)
            {
                if (side == Orientation::low) {
                    // Tangential directions have hoextrap
                    for (int i = 0; i < AMREX_SPACEDIM; i++)
                        bcs[BCVars::vels][i].setLo(dir, ERFBCType::foextrap);
                    // Only normal direction has ext_dir
                    bcs[BCVars::vels][dir].setLo(dir, ERFBCType::ext_dir);

                } else {
                    // Tangential directions have hoextrap
                    for (int i = 0; i < AMREX_SPACEDIM; i++)
                        bcs[BCVars::vels][i].setHi(dir, ERFBCType::foextrap);
                    // Only normal direction has ext_dir
                    bcs[BCVars::vels][dir].setHi(dir, ERFBCType::ext_dir);
                }
            }
            else if (bct == BC::periodic)
            {
                if (side == Orientation::low) {
                    for (int i = 0; i < AMREX_SPACEDIM; i++)
                        bcs[BCVars::vels][i].setLo(dir, ERFBCType::int_dir);
                } else {
                    for (int i = 0; i < AMREX_SPACEDIM; i++)
                        bcs[BCVars::vels][i].setHi(dir, ERFBCType::int_dir);
                }
            }
            else if ( bct == BC::MOST )
            {
                if (dir == 2 && side == Orientation::low) {
                    for (int i = 0; i < AMREX_SPACEDIM; i++)
                        bcs[BCVars::vels][i].setLo(dir, ERFBCType::MOST);
                } else {
                    amrex::Error("MOST bc can only be applied on low z-face");
                }
            }
        }

#ifdef AMREX_USE_GPU
        Gpu::htod_memcpy
#else
        std::memcpy
#endif
        (bcs_d[BCVars::vels].data(), bcs[BCVars::vels].data(), sizeof(BCRec)*AMREX_SPACEDIM);
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
            auto const bct = m_bc_type[ori];
            if ( bct == BC::symmetry )
            {
                if (side == Orientation::low) {
                    for (int i = 0; i < NVAR; i++)
                        bcs[BCVars::cons][i].setLo(dir, ERFBCType::reflect_even);
                } else {
                    for (int i = 0; i < NVAR; i++)
                        bcs[BCVars::cons][i].setHi(dir, ERFBCType::reflect_even);
                }
            }
            else if ( bct == BC::outflow )
            {
                if (side == Orientation::low) {
                    for (int i = 0; i < NVAR; i++)
                        bcs[BCVars::cons][i].setLo(dir, ERFBCType::foextrap);
                } else {
                    for (int i = 0; i < NVAR; i++)
                        bcs[BCVars::cons][i].setHi(dir, ERFBCType::foextrap);
                }
            }
            else if ( bct == BC::no_slip_wall)
            {
                if (side == Orientation::low) {
                    for (int i = 0; i < NVAR; i++)
                        bcs[BCVars::cons][i].setLo(dir, ERFBCType::foextrap);
                    if (m_bc_cons[ori][RhoTheta_comp] > 0.)
                        bcs[BCVars::cons][RhoTheta_comp].setLo(dir, ERFBCType::ext_dir);
                } else {
                    for (int i = 0; i < NVAR; i++)
                        bcs[BCVars::cons][i].setHi(dir, ERFBCType::foextrap);
                    if (m_bc_cons[ori][RhoTheta_comp] > 0.)
                        bcs[BCVars::cons][RhoTheta_comp].setHi(dir, ERFBCType::ext_dir);
                }
            }
            else if (bct == BC::slip_wall)
            {
                if (side == Orientation::low) {
                    for (int i = 0; i < NVAR; i++)
                        bcs[BCVars::cons][i].setLo(dir, ERFBCType::foextrap);
                    if (m_bc_cons[ori][RhoTheta_comp] > 0.)
                        bcs[BCVars::cons][RhoTheta_comp].setLo(dir, ERFBCType::ext_dir);
                } else {
                    for (int i = 0; i < NVAR; i++)
                        bcs[BCVars::cons][i].setHi(dir, ERFBCType::foextrap);
                    if (m_bc_cons[ori][RhoTheta_comp] > 0.)
                        bcs[BCVars::cons][RhoTheta_comp].setHi(dir, ERFBCType::ext_dir);
                }
            }
            else if (bct == BC::inflow)
            {
                if (side == Orientation::low) {
                    for (int i = 0; i < NVAR; i++)
                        bcs[BCVars::cons][i].setLo(dir, ERFBCType::ext_dir);
                } else {
                    for (int i = 0; i < NVAR; i++)
                        bcs[BCVars::cons][i].setHi(dir, ERFBCType::ext_dir);
                }
            }
            else if (bct == BC::periodic)
            {
                if (side == Orientation::low) {
                    for (int i = 0; i < NVAR; i++)
                        bcs[BCVars::cons][i].setLo(dir, ERFBCType::int_dir);
                } else {
                    for (int i = 0; i < NVAR; i++)
                       bcs[BCVars::cons][i].setHi(dir, ERFBCType::int_dir);
                }
            }
            else if ( bct == BC::MOST )
            {
                if (dir == 2 && side == Orientation::low) {
                    for (int i = 0; i < NVAR; i++)
                        bcs[BCVars::cons][i].setLo(dir, ERFBCType::MOST);
                } else {
                    amrex::Error("MOST bc can only be applied on low z-face");
                }
            }
        }

#ifdef AMREX_USE_GPU
        Gpu::htod_memcpy
#else
        std::memcpy
#endif
        (bcs_d[BCVars::cons].data(), bcs[BCVars::cons].data(), sizeof(BCRec)*Cons::NumVars);
    }
}

