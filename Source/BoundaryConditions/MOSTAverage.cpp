#include <MOSTAverage.H>

// Constructor
MOSTAverage::MOSTAverage (amrex::Vector<const amrex::MultiFab*> fields,
                          amrex::Vector<amrex::MultiFab*> averages,
                          amrex::Vector<amrex::iMultiFab*> k_indx,
                          amrex::Vector<amrex::Real> z_ref,
                          const amrex::MultiFab* z_nd,
                          amrex::Geometry geom,
                          int policy,
                          bool update_k)
    : m_fields(fields), m_averages(averages), m_k_indx(k_indx), m_z_ref(z_ref)  ,
      m_z_nd(z_nd)    , m_geom(geom)        , m_policy(policy), m_update_k(update_k)
{
    AMREX_ALWAYS_ASSERT(fields.size() >= 2);

    // Domain bottom and cell-size
    m_zlo = m_geom.ProbLo  (2);
    m_dz  = m_geom.CellSize(2);

    // Num components, plane avg, cells per plane, interp coeff
    amrex::Box domain = m_geom.Domain();
    amrex::IntVect dom_lo(domain.loVect());
    amrex::IntVect dom_hi(domain.hiVect());
    int asize = m_averages.size();
    m_ncomps.resize( asize );
    m_plane_average.resize( asize );
    m_ncell_plane.resize( asize );
    c_interp.resize( asize );
    for (int i(0); i<asize; ++i) {
        m_ncomps[i] = m_averages[i]->nComp();
        m_plane_average[i].resize(2,0.0);

        m_ncell_plane[i] = 1;
        amrex::IndexType ixt = m_averages[i]->boxArray().ixType();
        for (int j = 0; j < AMREX_SPACEDIM; ++j) {
            if (j != 2) {
                if (ixt.nodeCentered(j)) {
                    m_ncell_plane[i] *= (dom_hi[j] - dom_lo[j] + 2);
                } else {
                    m_ncell_plane[i] *= (dom_hi[j] - dom_lo[j] + 1);
                }
            }
        }
    }

    // Select a policy for k-index
    switch(m_policy) {
    case 0: // Standard plane average
        if (m_update_k) {
            compute_plane_k_indices();
        } else {
            for (int i(0); i<m_averages.size(); ++i) {
                amrex::Real z = m_z_ref[i];
                int k_lo = static_cast<int>(floor((z - m_zlo) / m_dz - 0.5));
                const amrex::Real z_lo = m_zlo + (k_lo + 0.5) * m_dz;
                c_interp[i] = (z - z_lo) / m_dz;
            }
        }
        break;

    case 1: // Local region/point
        if (m_update_k) compute_point_k_indices();
        break;

    case 2: // Fixed height above the terrain surface
        break;

    case 3: // Along a normal vector
        break;

    default:
        AMREX_ASSERT_WITH_MESSAGE(false, "Unknown policy for MOSTAverage!");
    }
}

// Populate a 2D iMF with the k indices for averaging
void
MOSTAverage::compute_plane_k_indices()
{
    for (int iavg(0); iavg<m_averages.size(); ++iavg) {
        amrex::Real z = m_z_ref[iavg];
        AMREX_ALWAYS_ASSERT(z >= m_zlo + 0.5 * m_dz);

        int k_lo = static_cast<int>(floor((z - m_zlo) / m_dz - 0.5));
        int k_hi = k_lo + 1;
        const amrex::Real z_lo = m_zlo + (k_lo + 0.5) * m_dz;
        c_interp[iavg] = (z - z_lo) / m_dz;

        amrex::IntVect ng = m_k_indx[iavg]->nGrowVect(); ng[2]=0;

        for (amrex::MFIter mfi(*m_k_indx[iavg], amrex::TilingIfNotGPU()); mfi.isValid(); ++mfi) {
            amrex::Box gbx = mfi.growntilebox(ng);

            auto k_arr     = m_k_indx[iavg]->array(mfi);

            ParallelFor(gbx,gbx,
            [=] AMREX_GPU_DEVICE (int i, int j, int k) noexcept {
                k_arr(i,j,k,0) = k_lo;
            },
            [=] AMREX_GPU_DEVICE (int i, int j, int k) noexcept {
                k_arr(i,j,k,1) = k_hi;
            });
       }
    }
}

// Populate a 2D iMF with the k indices for averaging
void
MOSTAverage::compute_point_k_indices()
{
    amrex::ParmParse pp("erf");
    auto read_k = pp.query("most.k_indx", m_k_ind);
    pp.query("most.radius", m_radius);

    // Specified k index
    if (read_k) {
        m_k_ind = std::max(m_radius,m_k_ind);
        for (int iavg(0); iavg<m_averages.size(); ++iavg) m_k_indx[iavg]->setVal(m_k_ind);
    // Set k from zref
    } else {
        for (int iavg(0); iavg<m_averages.size(); ++iavg) {
            amrex::Real z = m_z_ref[iavg];
            AMREX_ALWAYS_ASSERT(z >= m_zlo + 0.5 * m_dz);

            int k_lo = static_cast<int>(floor((z - m_zlo) / m_dz - 0.5));
            k_lo = std::max(m_radius,k_lo);

            amrex::IntVect ng = m_k_indx[iavg]->nGrowVect(); ng[2]=0;

            for (amrex::MFIter mfi(*m_k_indx[iavg], amrex::TilingIfNotGPU()); mfi.isValid(); ++mfi) {
                amrex::Box gbx = mfi.growntilebox(ng);
                auto k_arr     = m_k_indx[iavg]->array(mfi);

                ParallelFor(gbx, [=] AMREX_GPU_DEVICE (int i, int j, int k) noexcept
                {
                    k_arr(i,j,k,0) = k_lo;
                });
            } // MFiter
        } // iavg
    } // read_k
}


// Driver to call appropriate average member function
void
MOSTAverage::compute_averages()
{
    switch(m_policy) {
    case 0: // Standard plane average
        compute_plane_averages();
        break;

    case 1: // Local region/point
        compute_point_averages();
        break;

    case 2: // Fixed height above the terrain surface
        break;

    case 3: // Along a normal vector
        break;

    default:
        AMREX_ASSERT_WITH_MESSAGE(false, "Unknown policy for MOSTAverage!");
    }
}

// Fill plane storage with averages
void
MOSTAverage::compute_plane_averages()
{
    // Averages over all the fields
    //----------------------------------------------------------
    for (int imf(0); imf<m_fields.size(); ++imf) {
        const amrex::Real denom = 1.0 / (amrex::Real)m_ncell_plane[imf];

        amrex::AsyncArray<amrex::Real> pavg(m_plane_average[imf].data(), m_plane_average[imf].size());
        amrex::Real* plane_avg = pavg.data();

#ifdef _OPENMP
#pragma omp parallel if (amrex::Gpu::notInLaunchRegion())
#endif
        for (amrex::MFIter mfi(*m_fields[imf], amrex::TilingIfNotGPU()); mfi.isValid(); ++mfi) {
            amrex::Box bx  = mfi.tilebox();
            amrex::Box pbx = bx; pbx.setSmall(2,0); pbx.setBig(2,1);

            auto mf_arr = m_fields[imf]->const_array(mfi);
            auto k_arr  = m_k_indx[imf]->const_array(mfi);

            ParallelFor(amrex::Gpu::KernelInfo().setReduction(true), pbx, [=]
            AMREX_GPU_DEVICE(int i, int j, int n, amrex::Gpu::Handler const& handler) noexcept
            {
                int k = k_arr(i,j,0,n);
                amrex::Gpu::deviceReduceSum(&plane_avg[n],mf_arr(i, j, k) * denom, handler);
            });
        }

        pavg.copyToHost(m_plane_average[imf].data(), m_plane_average[imf].size());
        amrex::ParallelDescriptor::ReduceRealSum(m_plane_average[imf].data(), m_plane_average[imf].size());

        // No spatial variation with plane averages
        amrex::Real c_val = c_interp[imf];
        amrex::Real a_val = m_plane_average[imf][0]*(1.0 - c_val) + m_plane_average[imf][1]*c_val;
        m_averages[imf]->setVal(a_val);
    }

    // Averages for the tangential velocity magnitude
    //----------------------------------------------------------
    int imf  = 0;
    int iavg = m_fields.size();
    const amrex::Real denom = 1.0 / (amrex::Real)m_ncell_plane[iavg];

    amrex::AsyncArray<amrex::Real> pavg(m_plane_average[iavg].data(), m_plane_average[iavg].size());
    amrex::Real* plane_avg = pavg.data();

#ifdef _OPENMP
#pragma omp parallel if (amrex::Gpu::notInLaunchRegion())
#endif
    for (amrex::MFIter mfi(*m_fields[imf], amrex::TilingIfNotGPU()); mfi.isValid(); ++mfi) {
        amrex::Box bx  = amrex::enclosedCells(mfi.tilebox());
        amrex::Box pbx = bx; pbx.setSmall(2,0); pbx.setBig(2,1);

        auto u_mf_arr = m_fields[imf]->const_array(mfi);
        auto v_mf_arr = m_fields[imf+1]->const_array(mfi);

        auto k_u_arr  = m_k_indx[imf]->const_array(mfi);
        auto k_v_arr  = m_k_indx[imf+1]->const_array(mfi);

        ParallelFor(amrex::Gpu::KernelInfo().setReduction(true), pbx, [=]
        AMREX_GPU_DEVICE(int i, int j, int n, amrex::Gpu::Handler const& handler) noexcept
        {
            int k_u = k_u_arr(i,j,0,n);
            int k_v = k_v_arr(i,j,0,n);

            const amrex::Real u_val = 0.5 * (u_mf_arr(i,j,k_u) + u_mf_arr(i+1,j  ,k_u));
            const amrex::Real v_val = 0.5 * (v_mf_arr(i,j,k_v) + v_mf_arr(i  ,j+1,k_v));
            const amrex::Real mag   = std::sqrt(u_val*u_val + v_val*v_val);
            amrex::Gpu::deviceReduceSum(&plane_avg[n],mag * denom, handler);
        });
    }

    pavg.copyToHost(m_plane_average[iavg].data(), m_plane_average[iavg].size());
    amrex::ParallelDescriptor::ReduceRealSum(m_plane_average[iavg].data(), m_plane_average[iavg].size());

    // No spatial variation with plane averages
    amrex::Real c_val = c_interp[iavg];
    amrex::Real a_val = m_plane_average[iavg][0]*(1.0 - c_val) + m_plane_average[iavg][1]*c_val;
    m_averages[iavg]->setVal(a_val);
}

// Fill 2D MF with local averages
void
MOSTAverage::compute_point_averages()
{
    // Number of cells contained in the local average
    int ncell_avg = (2 * m_radius + 1) * (2 * m_radius + 1) * (2 * m_radius + 1);
    const amrex::Real denom = 1.0 / (amrex::Real) ncell_avg;

    // Capture radius for device
    int d_radius = m_radius;

    // Averages over all the fields
    //----------------------------------------------------------
    for (int imf(0); imf<m_fields.size(); ++imf) {
#ifdef _OPENMP
#pragma omp parallel if (amrex::Gpu::notInLaunchRegion())
#endif
        for (amrex::MFIter mfi(*m_fields[imf], amrex::TilingIfNotGPU()); mfi.isValid(); ++mfi) {
            amrex::Box bx  = mfi.tilebox();
            amrex::Box pbx = bx; pbx.setSmall(2,0); pbx.setBig(2,0);

            auto mf_arr = m_fields[imf]->const_array(mfi);
            auto k_arr  = m_k_indx[imf]->const_array(mfi);
            auto ma_arr = m_averages[imf]->array(mfi);

            ParallelFor(pbx, [=] AMREX_GPU_DEVICE(int i, int j, int) noexcept
            {
                ma_arr(i,j,0) = 0.0;

                int k = k_arr(i,j,0,0);
                for (int lk(k-d_radius); lk<=(k+d_radius); ++lk) {
                    for (int lj(j-d_radius); lj<=(j+d_radius); ++lj) {
                        for (int li(i-d_radius); li<=(i+d_radius); ++li) {
                            ma_arr(i,j,0) += mf_arr(li, lj, lk) * denom;
                        }
                    }
                }
            });
        }

        // Fill interior ghost cells and any ghost cells outside a periodic domain
        //***********************************************************************************
        m_averages[imf]->FillBoundary(m_geom.periodicity());
    }

    // Averages for the tangential velocity magnitude
    //----------------------------------------------------------
    int imf  = 0;
    int iavg = m_fields.size();

#ifdef _OPENMP
#pragma omp parallel if (amrex::Gpu::notInLaunchRegion())
#endif
    for (amrex::MFIter mfi(*m_averages[iavg], amrex::TilingIfNotGPU()); mfi.isValid(); ++mfi) {
        amrex::Box bx  = mfi.tilebox();
        amrex::Box pbx = bx; pbx.setSmall(2,0); pbx.setBig(2,0);

        auto u_mf_arr = m_fields[imf]->const_array(mfi);
        auto v_mf_arr = m_fields[imf+1]->const_array(mfi);
        auto k_arr    = m_k_indx[imf]->const_array(mfi);
        auto ma_arr   = m_averages[iavg]->array(mfi);

        ParallelFor(pbx, [=] AMREX_GPU_DEVICE(int i, int j, int) noexcept
        {
            ma_arr(i,j,0) = 0.0;

            int k = k_arr(i,j,0,0);
            for (int lk(k-d_radius); lk<=(k+d_radius); ++lk) {
                for (int lj(j-d_radius); lj<=(j+d_radius); ++lj) {
                    for (int li(i-d_radius); li<=(i+d_radius); ++li) {
                        const amrex::Real u_val = 0.5 * (u_mf_arr(li,lj,lk) + u_mf_arr(li+1,lj  ,lk));
                        const amrex::Real v_val = 0.5 * (v_mf_arr(li,lj,lk) + v_mf_arr(li  ,lj+1,lk));
                        const amrex::Real mag   = std::sqrt(u_val*u_val + v_val*v_val);
                        ma_arr(i,j,0) += mag * denom;
                    }
                }
            }
        });

        // Fill interior ghost cells and any ghost cells outside a periodic domain
        //***********************************************************************************
        m_averages[iavg]->FillBoundary(m_geom.periodicity());
    }


    // Need to fill ghost cells outside the domain if not periodic
    bool not_per_x = !(m_geom.periodicity().isPeriodic(0));
    bool not_per_y = !(m_geom.periodicity().isPeriodic(1));
    if (not_per_x || not_per_y) {
        amrex::Box domain = m_geom.Domain();
        for (int ima(0); ima<m_averages.size(); ++ima) {
#ifdef _OPENMP
#pragma omp parallel if (amrex::Gpu::notInLaunchRegion())
#endif
            amrex::IndexType ixt     = m_averages[ima]->boxArray().ixType();
            amrex::Box ldomain       = domain; ldomain.convert(ixt);
            amrex::IntVect ng = m_averages[ima]->nGrowVect(); ng[2]=0;
            for (amrex::MFIter mfi(*m_averages[ima], amrex::TilingIfNotGPU()); mfi.isValid(); ++mfi) {
                amrex::Box gbx = mfi.growntilebox(ng);
                gbx.setSmall(2,0); gbx.setBig(2,0);

                if (ldomain.contains(gbx)) continue;

                auto ma_arr = m_averages[ima]->array(mfi);

                int i_lo = ldomain.smallEnd(0); int i_hi = ldomain.bigEnd(0);
                int j_lo = ldomain.smallEnd(1); int j_hi = ldomain.bigEnd(1);
                ParallelFor(gbx, [=] AMREX_GPU_DEVICE(int i, int j, int) noexcept
                {
                    int li, lj;
                    li = i  < i_lo ? i_lo : i;
                    li = li > i_hi ? i_hi : li;
                    lj = j  < j_lo ? j_lo : j;
                    lj = lj > j_hi ? j_hi : lj;

                    ma_arr(i,j,0) = ma_arr(li,lj,0);
                });
            } // MFiter
        } // ima
    } // Not periodic
}

// Write k indices
void
MOSTAverage::write_k_indices()
{
    // Last component u_mag_mean is cc
    int navg = m_averages.size() - 1;

    std::ofstream ofile;
    ofile.open ("MOST_k_indices.txt");
    ofile << "K indices used to compute averages via MOSTAverages class:\n";

    for (amrex::MFIter mfi(*m_averages[navg], amrex::TilingIfNotGPU()); mfi.isValid(); ++mfi) {
        amrex::Box bx  = mfi.tilebox(); bx.setBig(2,0);
        int il = bx.smallEnd(0); int iu = bx.bigEnd(0);
        int jl = bx.smallEnd(1); int ju = bx.bigEnd(1);

        for (int j(jl); j<=ju; ++j) {
            for (int i(il); i<=iu; ++i) {
                ofile << "(I,J): " << "(" << i << "," << j << ")" << "\n";
                int k = 0;
                for (int iavg(0); iavg<=navg; ++iavg) {
                    auto k_arr = m_k_indx[iavg]->array(mfi);
                    ofile << "iavg K_lo K_hi c_interp: "
                          << iavg << ' '
                          << k_arr(i,j,k,0) << ' '
                          << k_arr(i,j,k,1) << ' '
                          << c_interp[iavg] << "\n";
                }
                ofile << "\n";
            }
        }
    }
    ofile.close();
}

// Write averages
void
MOSTAverage::write_averages()
{
    // Last component u_mag_mean is cc
    int navg = m_averages.size() - 1;

    std::ofstream ofile;
    ofile.open ("MOST_averages.txt");
    ofile << "Averages computed via MOSTAverages class:\n";

    for (amrex::MFIter mfi(*m_averages[navg], amrex::TilingIfNotGPU()); mfi.isValid(); ++mfi) {
        amrex::Box bx  = mfi.tilebox(); bx.setBig(2,0);
        int il = bx.smallEnd(0); int iu = bx.bigEnd(0);
        int jl = bx.smallEnd(1); int ju = bx.bigEnd(1);

        for (int j(jl); j<=ju; ++j) {
            for (int i(il); i<=iu; ++i) {
                ofile << "(I,J): " << "(" << i << "," << j << ")" << "\n";
                int k = 0;
                for (int iavg(0); iavg<=navg; ++iavg) {
                    auto mf_arr = m_averages[iavg]->array(mfi);
                    ofile << "iavg val: "
                          << iavg << ' '
                          << mf_arr(i,j,k) << "\n";
                }
                ofile << "\n";
            }
        }
    }
    ofile.close();
}
