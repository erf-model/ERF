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
                int z = m_z_ref[i];
                int k_lo = static_cast<int>(floor((z - m_zlo) / m_dz - 0.5));
                const amrex::Real z_lo = m_zlo + (k_lo + 0.5) * m_dz;
                c_interp[i] = (z - z_lo) / m_dz;
            }
        }
        break;

    case 1: // Fixed height above the terrain surface
        break;

    case 2: // Along a normal vector
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
        int z = m_z_ref[iavg];
        AMREX_ALWAYS_ASSERT(z >= m_zlo + 0.5 * m_dz);

        int k_lo = static_cast<int>(floor((z - m_zlo) / m_dz - 0.5));
        int k_hi = k_lo + 1;
        const amrex::Real z_lo = m_zlo + (k_lo + 0.5) * m_dz;
        c_interp[iavg] = (z - z_lo) / m_dz;

        const amrex::IntVect& ng = m_k_indx[iavg]->nGrowVect();

        for (amrex::MFIter mfi(*m_k_indx[iavg], amrex::TilingIfNotGPU()); mfi.isValid(); ++mfi) {
            amrex::Box gbx = mfi.growntilebox({ng[0], ng[1], 0});
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

// Driver to call appropriate average member function
void
MOSTAverage::compute_averages()
{
    switch(m_policy) {
    case 0: // Standard plane average
        compute_plane_averages();
        break;

    case 1: // Fixed height above the terrain surface
        break;

    case 2: // Along a normal vector
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
        const int ncomp = m_ncomps[imf];
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
    const int ncomp = m_ncomps[iavg];
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

// Write k indices
void
MOSTAverage::write_k_indices()
{
    // Last component u_mag_mean is cc
    int navg = m_averages.size() - 1;

    std::ofstream ofile;
    ofile.open ("MOST_K_indices.txt");
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
