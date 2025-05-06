#include "ERF.H"
#include "ERF_Utils.H"

#include "AMReX_MLMG.H"
#include "AMReX_MLNodeLaplacian.H"

using namespace amrex;

/**
 * Calculate wall distances using the Poisson equation
 *
 * The zlo boundary is assumed to correspond to the land surface. If there are
 * no boundary walls, then the other use case is to calculate wall distances
 * for immersed boundaries (embedded or thin body).
 *
 * See Tucker, P. G. (2003). Differential equation-based wall distance
 * computation for DES and RANS. Journal of Computational Physics,
 * 190(1), 229–248. https://doi.org/10.1016/S0021-9991(03)00272-9
 */
void ERF::poisson_wall_dist (int lev)
{
    BL_PROFILE("ERF::poisson_wall_dist()");

    bool havewall{false};
    Orientation zlo(Direction::z, Orientation::low);
    if ( ( phys_bc_type[zlo] == ERF_BC::surface_layer                      ) ||
         ( phys_bc_type[zlo] == ERF_BC::no_slip_wall                       ) )/*||
         ((phys_bc_type[zlo] == ERF_BC::slip_wall) && (dom_hi.z > dom_lo.z)) )*/
    {
        havewall = true;
    }

    auto const& geomdata = geom[lev];

    if (havewall) {
        if (solverChoice.mesh_type == MeshType::ConstantDz) {
// Comment this out to test the wall dist calc in the trivial case:
//#if 0
            Print() << "Directly calculating direct wall distance for constant dz" << std::endl;
            const Real* prob_lo = geomdata.ProbLo();
            const Real* dx = geomdata.CellSize();
            for (MFIter mfi(*walldist[lev]); mfi.isValid(); ++mfi) {
                const Box& bx = mfi.validbox();
                auto dist_arr = walldist[lev]->array(mfi);
                ParallelFor(bx, [=] AMREX_GPU_DEVICE(int i, int j, int k) {
                    dist_arr(i, j, k) = prob_lo[2] + (k + 0.5) * dx[2];
                });
            }
            return;
//#endif
        } else if (solverChoice.mesh_type == MeshType::StretchedDz) {
            // TODO: Handle this trivial case
            Error("Wall dist calc not implemented with grid stretching yet");
        } else {
            // TODO
            Error("Wall dist calc not implemented over terrain yet");
        }
    }

    Print() << "Calculating Poisson wall distance" << std::endl;

    // Make sure the solver only sees the levels over which we are solving
    BoxArray nba = walldist[lev]->boxArray();
    nba.surroundingNodes();
    Vector<Geometry>          geom_tmp; geom_tmp.push_back(geom[lev]);
    Vector<BoxArray>            ba_tmp;   ba_tmp.push_back(nba);
    Vector<DistributionMapping> dm_tmp;   dm_tmp.push_back(walldist[lev]->DistributionMap());

    Vector<MultiFab> rhs;
    Vector<MultiFab> phi;

    if (solverChoice.terrain_type == TerrainType::EB) {
        amrex::Error("Wall dist calc not implemented for EB");
    } else {
        rhs.resize(1);   rhs[0].define(ba_tmp[0], dm_tmp[0], 1, 0);
        phi.resize(1);   phi[0].define(ba_tmp[0], dm_tmp[0], 1, 1);
    }

    rhs[0].setVal(-1.0);

    // Define an overset mask to set dirichlet nodes on walls
    iMultiFab mask(ba_tmp[0], dm_tmp[0], 1, 0);
    Vector<const iMultiFab*> overset_mask = {&mask};

    auto const dom_lo = lbound(geom[lev].Domain());
    auto const dom_hi = ubound(geom[lev].Domain());

    // ****************************************************************************
    // Initialize phi
    // (It is essential that we do this in order to fill the corners; this is
    // used if we include blanking.)
    // ****************************************************************************
    phi[0].setVal(0.0);

    // ****************************************************************************
    // Interior boundaries are marked with phi=0
    // ****************************************************************************
    // Overset mask is 0/1: 1 means the node is an unknown. 0 means it's known.
    mask.setVal(1);
    if (solverChoice.advChoice.have_zero_flux_faces) {
        Warning("Poisson distance is inaccurate for bodies in open domains that are small compared to the domain size, skipping");
        return;
#if 0
        Gpu::DeviceVector<IntVect> xfacelist, yfacelist, zfacelist;

        xfacelist.resize(solverChoice.advChoice.zero_xflux.size());
        yfacelist.resize(solverChoice.advChoice.zero_yflux.size());
        zfacelist.resize(solverChoice.advChoice.zero_zflux.size());

        if (xfacelist.size() > 0) {
            Gpu::copy(amrex::Gpu::hostToDevice,
                      solverChoice.advChoice.zero_xflux.begin(),
                      solverChoice.advChoice.zero_xflux.end(),
                      xfacelist.begin());
            Print() << "  masking interior xfaces" << std::endl;
        }
        if (yfacelist.size() > 0) {
            Gpu::copy(amrex::Gpu::hostToDevice,
                      solverChoice.advChoice.zero_yflux.begin(),
                      solverChoice.advChoice.zero_yflux.end(),
                      yfacelist.begin());
            Print() << "  masking interior yfaces" << std::endl;
        }
        if (zfacelist.size() > 0) {
            Gpu::copy(amrex::Gpu::hostToDevice,
                      solverChoice.advChoice.zero_zflux.begin(),
                      solverChoice.advChoice.zero_zflux.end(),
                      zfacelist.begin());
            Print() << "  masking interior zfaces" << std::endl;
        }

        for (MFIter mfi(phi[0]); mfi.isValid(); ++mfi) {
            const Box& bx = mfi.validbox();

            auto phi_arr  = phi[0].array(mfi);
            auto mask_arr = mask.array(mfi);

            if (xfacelist.size() > 0) {
                ParallelFor(bx, [=] AMREX_GPU_DEVICE(int i, int j, int k) {
                    for (int iface=0; iface < xfacelist.size(); ++iface) {
                        if ((i == xfacelist[iface][0]) &&
                            (j == xfacelist[iface][1]) &&
                            (k == xfacelist[iface][2]))
                        {
                            mask_arr(i, j  , k  ) = 0;
                            mask_arr(i, j  , k+1) = 0;
                            mask_arr(i, j+1, k  ) = 0;
                            mask_arr(i, j+1, k+1) = 0;
                        }
                    }
                });
            }

            if (yfacelist.size() > 0) {
                ParallelFor(bx, [=] AMREX_GPU_DEVICE(int i, int j, int k) {
                    for (int iface=0; iface < yfacelist.size(); ++iface) {
                        if ((i == yfacelist[iface][0]) &&
                            (j == yfacelist[iface][1]) &&
                            (k == yfacelist[iface][2]))
                        {
                            mask_arr(i  , j, k  ) = 0;
                            mask_arr(i  , j, k+1) = 0;
                            mask_arr(i+1, j, k  ) = 0;
                            mask_arr(i+1, j, k+1) = 0;
                        }
                    }
                });
            }

            if (zfacelist.size() > 0) {
                ParallelFor(bx, [=] AMREX_GPU_DEVICE(int i, int j, int k) {
                    for (int iface=0; iface < zfacelist.size(); ++iface) {
                        if ((i == xfacelist[iface][0]) &&
                            (j == xfacelist[iface][1]) &&
                            (k == xfacelist[iface][2]))
                        {
                            mask_arr(i  , j  , k) = 0;
                            mask_arr(i  , j+1, k) = 0;
                            mask_arr(i+1, j  , k) = 0;
                            mask_arr(i+1, j+1, k) = 0;
                        }
                    }
                });
            }
        }
#endif
    }

    // ****************************************************************************
    // Setup BCs, with solid domain boundaries being dirichlet
    // ****************************************************************************
    amrex::Array<amrex::LinOpBCType,AMREX_SPACEDIM> bc3d_lo, bc3d_hi;
    for (int dir = 0; dir < AMREX_SPACEDIM; ++dir) {
        if (geom[0].isPeriodic(dir)) {
            bc3d_lo[dir] = LinOpBCType::Periodic;
            bc3d_hi[dir] = LinOpBCType::Periodic;
        } else {
            bc3d_lo[dir] = LinOpBCType::Neumann;
            bc3d_hi[dir] = LinOpBCType::Neumann;
        }
    }
    if (havewall) {
        Print() << "  Poisson zlo BC is dirichlet" << std::endl;
        bc3d_lo[2] = LinOpBCType::Dirichlet;
    }
    Print() << "  bc lo : " << bc3d_lo << std::endl;
    Print() << "  bc hi : " << bc3d_hi << std::endl;

    if (!solverChoice.advChoice.have_zero_flux_faces && !havewall) {
        Error("No solid boundaries in the computational domain");
    }

    LPInfo info;
/* Nodal solver cannot have hidden dimensions */
#if 0
    // Allow a hidden direction if the domain is one cell wide
    if (dom_lo.x == dom_hi.x) {
        info.setHiddenDirection(0);
        Print() << "  domain is 2D in yz" << std::endl;
    } else if (dom_lo.y == dom_hi.y) {
        info.setHiddenDirection(1);
        Print() << "  domain is 2D in xz" << std::endl;
    } else if (dom_lo.z == dom_hi.z) {
        info.setHiddenDirection(2);
        Print() << "  domain is 2D in xy" << std::endl;
    }
#endif

    // ****************************************************************************
    // Solve nodal masked Poisson problem with MLMG
    // TODO: different solver for terrain?
    // ****************************************************************************
    const Real reltol = solverChoice.poisson_reltol;
    const Real abstol = solverChoice.poisson_abstol;

    Real sigma = 1.0;
    Vector<EBFArrayBoxFactory const*> factory_vec = { &EBFactory(lev) };
    MLNodeLaplacian mlpoisson(geom_tmp, ba_tmp, dm_tmp, info, factory_vec, sigma);

    mlpoisson.setDomainBC(bc3d_lo, bc3d_hi);

    if (lev > 0) {
        mlpoisson.setCoarseFineBC(nullptr, ref_ratio[lev-1], LinOpBCType::Neumann);
    }

    mlpoisson.setLevelBC(0, nullptr);

    mlpoisson.setOversetMask(0, mask);

    // Solve
    MLMG mlmg(mlpoisson);
    int max_iter = 100;
    mlmg.setMaxIter(max_iter);

    mlmg.setVerbose(mg_verbose);
    mlmg.setBottomVerbose(0);

    mlmg.solve(GetVecOfPtrs(phi),
               GetVecOfConstPtrs(rhs),
               reltol, abstol);

    // Now overwrite with periodic fill outside domain and fine-fine fill inside
    phi[0].FillBoundary(geom[lev].periodicity());

    // ****************************************************************************
    // Compute grad(phi) to get distances
    // - Note that phi is nodal and walldist is cell-centered
    // - TODO: include terrain metrics for dphi/dz
    // ****************************************************************************
    for (MFIter mfi(*walldist[lev]); mfi.isValid(); ++mfi) {
        const Box& bx = mfi.validbox();

        const auto invCellSize = geomdata.InvCellSizeArray();

        auto const& phi_arr = phi[0].const_array(mfi);
        auto dist_arr = walldist[lev]->array(mfi);

        ParallelFor(bx, [=] AMREX_GPU_DEVICE(int i, int j, int k) {
            Real dpdx{0}, dpdy{0}, dpdz{0};

            // dphi/dx
            if (dom_lo.x != dom_hi.x) {
                dpdx = 0.25 * invCellSize[0] * (
                        (phi_arr(i+1, j  , k  ) - phi_arr(i, j  , k  ))
                      + (phi_arr(i+1, j  , k+1) - phi_arr(i, j  , k+1))
                      + (phi_arr(i+1, j+1, k  ) - phi_arr(i, j+1, k  ))
                      + (phi_arr(i+1, j+1, k+1) - phi_arr(i, j+1, k+1)) );
            }

            // dphi/dy
            if (dom_lo.y != dom_hi.y) {
                dpdy = 0.25 * invCellSize[1] * (
                        (phi_arr(i  , j+1, k  ) - phi_arr(i  , j, k  ))
                      + (phi_arr(i  , j+1, k+1) - phi_arr(i  , j, k+1))
                      + (phi_arr(i+1, j+1, k  ) - phi_arr(i+1, j, k  ))
                      + (phi_arr(i+1, j+1, k+1) - phi_arr(i+1, j, k+1)) );
            }

            // dphi/dz
            if (dom_lo.z != dom_hi.z) {
                dpdz = 0.25 * invCellSize[2] * (
                        (phi_arr(i  , j  , k+1) - phi_arr(i  , j  , k))
                      + (phi_arr(i  , j+1, k+1) - phi_arr(i  , j+1, k))
                      + (phi_arr(i+1, j  , k+1) - phi_arr(i+1, j  , k))
                      + (phi_arr(i+1, j+1, k+1) - phi_arr(i+1, j+1, k)) );
            }

            Real dp_dot_dp = dpdx*dpdx + dpdy*dpdy + dpdz*dpdz;
            Real phi_avg = 0.125 * (
                    phi_arr(i  , j  , k  ) + phi_arr(i  , j  , k+1) + phi_arr(i  , j+1, k  ) + phi_arr(i  , j+1, k+1)
                  + phi_arr(i+1, j  , k  ) + phi_arr(i+1, j  , k+1) + phi_arr(i+1, j+1, k  ) + phi_arr(i+1, j+1, k+1) );
            dist_arr(i, j, k) = -std::sqrt(dp_dot_dp) + std::sqrt(dp_dot_dp + 2*phi_avg);

            // DEBUG: output phi instead
            //dist_arr(i, j, k) = phi_arr(i, j, k);
        });
    }
}
