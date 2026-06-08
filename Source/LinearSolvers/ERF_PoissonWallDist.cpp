#include "ERF.H"
#include "ERF_Utils.H"
#include "ERF_TerrainPoisson_3D_K.H"
#include "ERF_TerrainMetrics.H"

#include "AMReX_MLMG.H"
#include "AMReX_MLABecLaplacian.H"

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
 * 190(1), 229–Real(248.) https://doi.org/Real(10.1016)/S0021-9991(03)00272-9
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
    auto const& dxinv    = geomdata.InvCellSizeArray();

    auto const& zphys_arr = z_phys_nd[lev]->const_arrays();

    if (havewall) {
#if 1
        // Bypass wall dist calc in the trivial cases

        if (solverChoice.mesh_type == MeshType::ConstantDz) {
            Print() << "Directly calculating direct wall distance for constant dz" << std::endl;
            const Real* prob_lo = geomdata.ProbLo();
            const Real* dx = geomdata.CellSize();
            for (MFIter mfi(*walldist[lev]); mfi.isValid(); ++mfi) {
                const Box& bx = mfi.validbox();
                auto dist_arr = walldist[lev]->array(mfi);
                ParallelFor(bx, [=] AMREX_GPU_DEVICE(int i, int j, int k) {
                    dist_arr(i, j, k) = prob_lo[2] + (k + myhalf) * dx[2];
                });
            }
            return;
        }

        if (solverChoice.mesh_type == MeshType::StretchedDz) {
            Print() << "Directly calculating direct wall distance for stretched dz" << std::endl;
            for (MFIter mfi(*walldist[lev],TileNoZ()); mfi.isValid(); ++mfi) {
                const Box& bx = mfi.validbox();
                auto dist_arr = walldist[lev]->array(mfi);
                const auto zcc_arr = z_phys_cc[lev]->const_array(mfi);
                const auto znd_arr = z_phys_nd[lev]->const_array(mfi);
                ParallelFor(bx, [=] AMREX_GPU_DEVICE(int i, int j, int k) {
                    dist_arr(i, j, k) = zcc_arr(i, j, k) - znd_arr(i, j, 0);
                });
            }
            return;
        }
#endif
    }
    else
    {
        Error("No solid boundaries in the computational domain");
    }

    Print() << "Calculating Poisson wall distance for general terrain" << std::endl;

    // Make sure the solver only sees the levels over which we are solving
    Vector<Geometry>          geom_tmp; geom_tmp.push_back(geom[lev]);
    Vector<BoxArray>            ba_tmp;   ba_tmp.push_back(walldist[lev]->boxArray());
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
#if 0
    // Define an overset mask (0 or 1) to set dirichlet nodes on walls
    // 1 means the node is an unknown. 0 means it's known.
    iMultiFab mask(ba_tmp[0], dm_tmp[0], 1, 0);
    Vector<const iMultiFab*> overset_mask = {&mask};

    mask.setVal(1);
    if (solverChoice.advChoice.have_zero_flux_faces) {
        Warning("Poisson distance is inaccurate for bodies in open domains that are small compared to the domain size, skipping");
        return;

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
    }
#endif

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

    LPInfo info; // defaults

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

#if 0
    Vector<EBFArrayBoxFactory const*> factory_vec;
    factory_vec.push_back(static_cast<FabFactory<FArrayBox> const*>(&EBFactory(lev));
#endif

    // ****************************************************************************
    // Setup Poisson problem
    // (A \alpha - B \nabla \cdot \beta \nabla ) \phi = f
    //
    // In physical space:
    //   \nabla \cdot \nabla \phi = -1
    //
    // In computational space:
    //   grad(phi) = T^T \nabla \phi
    // and
    //   \nabla \cdot (h_zeta T (T^T \nabla \phi)) = -h_zeta
    // where T = inv(J), T^T is the transpose of inv(J)
    // ****************************************************************************
    constexpr Real constA = zero;
    constexpr Real constB = -one;

    MLABecLaplacian mlabec(geom_tmp, ba_tmp, dm_tmp, info);

    mlabec.setScalars(constA, constB);
    mlabec.setACoeffs(0, zero);
#if 1
    // Set beta coefficients at faces
    Array<MultiFab, AMREX_SPACEDIM> beta;

    for (int idim = 0; idim < AMREX_SPACEDIM; ++idim) {
        BoxArray ba_face = ba_tmp[0];
        ba_face.surroundingNodes(idim);  // Convert to face-centered in direction idim
        beta[idim].define(ba_face, dm_tmp[0], 1, 0);
    }

    auto beta0_arr = beta[0].arrays();
    auto beta1_arr = beta[1].arrays();
    auto beta2_arr = beta[2].arrays();

    // Note: This ignores the off-diagonal components of (h_zeta T T^T), which
    //       is equivalent to assuming that h_xi and h_eta are small.

    ParallelFor(beta[0], [=] AMREX_GPU_DEVICE(int b, int i, int j, int k) {
        beta0_arr[b](i, j, k) = Compute_h_zeta_AtIface(i, j, k, dxinv, zphys_arr[b]);;
    });
    ParallelFor(beta[1], [=] AMREX_GPU_DEVICE(int b, int i, int j, int k) {
        beta1_arr[b](i, j, k) = Compute_h_zeta_AtJface(i, j, k, dxinv, zphys_arr[b]);;
    });
    ParallelFor(beta[2], [=] AMREX_GPU_DEVICE(int b, int i, int j, int k) {
        Real inv_h_zeta = one / Compute_h_zeta_AtKface(i, j, k, dxinv, zphys_arr[b]);
        Real h_xi = Compute_h_xi_AtKface(i, j, k, dxinv, zphys_arr[b]);
        Real h_eta = Compute_h_eta_AtKface(i, j, k, dxinv, zphys_arr[b]);
        beta2_arr[b](i, j, k) = inv_h_zeta * (1 + h_xi*h_xi + h_eta*h_eta);
    });

    mlabec.setBCoeffs(0, GetArrOfConstPtrs(beta));

    // Set RHS := -h_zeta
    auto rhs_arr = rhs[0].arrays();
    ParallelFor(rhs[0], [=] AMREX_GPU_DEVICE(int b, int i, int j, int k) {
        rhs_arr[b](i, j, k) = -Compute_h_zeta_AtCellCenter(i, j, k, dxinv, zphys_arr[b]);
    });
#else
    mlabec.setBCoeffs(0, one);
#endif

    mlabec.setDomainBC(bc3d_lo, bc3d_hi);

    if (lev > 0) {
        mlabec.setCoarseFineBC(nullptr, ref_ratio[lev-1], LinOpBCType::Neumann);
    }

    // If we have inhomogeneous BCs -- do this after setCoarseFineBC
    mlabec.setLevelBC(0, nullptr);

    // ****************************************************************************
    // Solve Poisson problem with MLMG
    // ****************************************************************************
    const Real reltol = solverChoice.poisson_reltol;
    const Real abstol = solverChoice.poisson_abstol;
    const int n_corr = solverChoice.ncorr;
    constexpr int max_iter = 100;

    MLMG mlmg(mlabec);
    mlmg.setMaxIter(max_iter);
    mlmg.setVerbose(mg_verbose);
    mlmg.setBottomVerbose(0);

    for (int icorr=0; icorr <= n_corr; ++icorr) {
        Print()<< "Solving wall distance poisson, icorr=" << icorr << std::endl;

        mlmg.solve(GetVecOfPtrs(phi),
                   GetVecOfConstPtrs(rhs),
                   reltol, abstol);

        // ****************************************************************************
        // Apply BCs: dirichlet (odd) on zlo, neumann (even) / periodic elsewhere
        // ****************************************************************************

        // Overwrite with periodic fill outside domain and fine-fine fill inside
        phi[0].FillBoundary(geom[lev].periodicity());

        if (!geom[lev].isPeriodic(0)) {
            for (MFIter mfi(phi[0],true); mfi.isValid(); ++mfi)
            {
                Box bx = mfi.tilebox();
                const Array4<Real>& phi_arr = phi[0].array(mfi);
                if (bx.smallEnd(0) <= dom_lo.x) {
                    ParallelFor(makeSlab(bx,0,dom_lo.x),
                    [=] AMREX_GPU_DEVICE (int i, int j, int k)
                    {
                        phi_arr(i-1,j,k) =  phi_arr(i,j,k); // even BC
                    });
                } // lo x
                if (bx.bigEnd(0) >= dom_hi.x) {
                    ParallelFor(makeSlab(bx,0,dom_hi.x),
                    [=] AMREX_GPU_DEVICE (int i, int j, int k)
                    {
                        phi_arr(i+1,j,k) =  phi_arr(i,j,k); // even BC
                    });
                } // hi x
            } // mfi
        } // not periodic in x

        if (!geom[lev].isPeriodic(1)) {
            for (MFIter mfi(phi[0],true); mfi.isValid(); ++mfi)
            {
                Box bx = mfi.tilebox();
                Box bx2(bx); bx2.grow(0,1);
                const Array4<Real>& phi_arr = phi[0].array(mfi);
                if (bx.smallEnd(1) <= dom_lo.y) {
                    ParallelFor(makeSlab(bx2,1,dom_lo.y),
                    [=] AMREX_GPU_DEVICE (int i, int j, int k)
                    {
                        phi_arr(i,j-1,k) =  phi_arr(i,j,k); // even BC
                    });
                } // lo y
                if (bx.bigEnd(1) >= dom_hi.y) {
                    ParallelFor(makeSlab(bx2,1,dom_hi.y),
                    [=] AMREX_GPU_DEVICE (int i, int j, int k)
                    {
                        phi_arr(i,j+1,k) =  phi_arr(i,j,k); // even BC
                    });
                } // hi y

            } // mfi
        } // not periodic in y

        for (MFIter mfi(phi[0],true); mfi.isValid(); ++mfi)
        {
            Box bx = mfi.tilebox();
            Box bx3(bx); bx3.grow(0,1); bx3.grow(1,1);
            const Array4<Real>& phi_arr = phi[0].array(mfi);
            if (bx.smallEnd(2) <= dom_lo.z) {
                ParallelFor(makeSlab(bx3,2,dom_lo.z),
                [=] AMREX_GPU_DEVICE (int i, int j, int k)
                {
                    phi_arr(i,j,k-1) = -phi_arr(i,j,k); // ODD BC
                });
            } // lo z
            if (bx.bigEnd(2) >= dom_hi.z) {
                ParallelFor(makeSlab(bx3,2,dom_hi.z),
                [=] AMREX_GPU_DEVICE (int i, int j, int k)
                {
                    phi_arr(i,j,k+1) =  phi_arr(i,j,k); // even BC
                });
            } // hi z
        } // mfi

        // ****************************************************************************
        // Compute grad(phi) to get distances
        // ****************************************************************************
        auto const& phi_arr = phi[0].const_arrays();
        //auto rhs_arr = rhs[0].arrays();
        auto dist_arr = walldist[lev]->arrays();

        ParallelFor(*walldist[lev], [=] AMREX_GPU_DEVICE(int b, int i, int j, int k) {
            Real dpdx{0}, dpdy{0}, dpdz{0};

            dpdx = terrpoisson_flux_x(i, j, k, phi_arr[b], zphys_arr[b], dxinv[0]);
            dpdy = terrpoisson_flux_y(i, j, k, phi_arr[b], zphys_arr[b], dxinv[1]);
            if (k == dom_lo.z) {
                dpdz = terrpoisson_flux_zlo_dir(i, j, k, phi_arr[b], zphys_arr[b], dxinv[0], dxinv[1]);
            } else {
                // This returns 0 at the wall, hence the need for the separate calc above
                dpdz = terrpoisson_flux_z(i, j, k, phi_arr[b], zphys_arr[b], dxinv[0], dxinv[1]);
            }

            Real magsqr_dphi = dpdx*dpdx + dpdy*dpdy + dpdz*dpdz;
            Real mag_dphi = std::sqrt(magsqr_dphi);
#if 1
            // Tucker 2003 Eqn 2
            dist_arr[b](i, j, k) = -mag_dphi + std::sqrt(magsqr_dphi + 2*phi_arr[b](i, j, k));
#else
            // DEBUG: output phi instead
            if (i==0 && j==0) AllPrint() << "walldist"<<IntVect(i,j,k) << " = " << dist_arr[b](i,j,k) << std::endl;
            dist_arr[b](i, j, k) = phi_arr[b](i, j, k);
#endif
            // Update RHS source term to explicitly include cross-terms
            if (n_corr > 0) {
                // d/dxi ( h_xi * dphi/dzeta )
                Real phi_zeta_xlo = fourth * dxinv[2] * ( phi_arr[b](i  , j, k+1) - phi_arr[b](i  , j, k-1)
                                                      + phi_arr[b](i-1, j, k+1) - phi_arr[b](i-1, j, k-1) );
                Real phi_zeta_xhi = fourth * dxinv[2] * ( phi_arr[b](i  , j, k+1) - phi_arr[b](i  , j, k-1)
                                                      + phi_arr[b](i+1, j, k+1) - phi_arr[b](i+1, j, k-1) );
                Real h_xi_xlo = Compute_h_xi_AtIface(i  , j, k, dxinv, zphys_arr[b]);
                Real h_xi_xhi = Compute_h_xi_AtIface(i+1, j, k, dxinv, zphys_arr[b]);

                // d/deta ( h_eta * dphi/dzeta )
                Real phi_zeta_ylo = fourth * dxinv[2] * ( phi_arr[b](i, j  , k+1) - phi_arr[b](i, j  , k-1)
                                                      + phi_arr[b](i, j-1, k+1) - phi_arr[b](i, j-1, k-1) );
                Real phi_zeta_yhi = fourth * dxinv[2] * ( phi_arr[b](i, j  , k+1) - phi_arr[b](i, j  , k-1)
                                                      + phi_arr[b](i, j+1, k+1) - phi_arr[b](i, j+1, k-1) );
                Real h_eta_ylo = Compute_h_eta_AtJface(i, j  , k, dxinv, zphys_arr[b]);
                Real h_eta_yhi = Compute_h_eta_AtJface(i, j+1, k, dxinv, zphys_arr[b]);

                // d/dzeta ( h_xi * dphi/dxi )
                Real phi_xi_zlo = fourth * dxinv[0] * ( phi_arr[b](i+1, j, k  ) - phi_arr[b](i-1, j, k  )
                                                    + phi_arr[b](i+1, j, k-1) - phi_arr[b](i-1, j, k-1) );
                Real phi_xi_zhi = fourth * dxinv[0] * ( phi_arr[b](i+1, j, k  ) - phi_arr[b](i-1, j, k  )
                                                    + phi_arr[b](i+1, j, k+1) - phi_arr[b](i-1, j, k+1) );
                Real h_xi_zlo = Compute_h_xi_AtKface(i, j, k  , dxinv, zphys_arr[b]);
                Real h_xi_zhi = Compute_h_xi_AtKface(i, j, k+1, dxinv, zphys_arr[b]);

                // d/dzeta ( h_eta * dphi/deta )
                Real phi_eta_zlo = fourth * dxinv[1] * ( phi_arr[b](i, j+1, k  ) - phi_arr[b](i, j-1, k  )
                                                     + phi_arr[b](i, j+1, k-1) - phi_arr[b](i, j-1, k-1) );
                Real phi_eta_zhi = fourth * dxinv[1] * ( phi_arr[b](i, j+1, k  ) - phi_arr[b](i, j-1, k  )
                                                     + phi_arr[b](i, j+1, k+1) - phi_arr[b](i, j-1, k+1) );
                Real h_eta_zlo = Compute_h_eta_AtKface(i, j, k  , dxinv, zphys_arr[b]);
                Real h_eta_zhi = Compute_h_eta_AtKface(i, j, k+1, dxinv, zphys_arr[b]);

                Real detJ = Compute_h_zeta_AtCellCenter(i, j, k, dxinv, zphys_arr[b]);

                rhs_arr[b](i, j, k) = -detJ
                                    + dxinv[0] * ( h_xi_xhi * phi_zeta_xhi - h_xi_xlo * phi_zeta_xlo)
                                    + dxinv[1] * ( h_eta_yhi * phi_zeta_yhi - h_eta_ylo * phi_zeta_ylo)
                                    + dxinv[2] * ( h_xi_zhi * phi_xi_zhi - h_xi_zlo * phi_xi_zlo
                                                 + h_eta_zhi * phi_eta_zhi - h_eta_zlo * phi_eta_zlo);
            }
        });
    } // corrector loop
}
