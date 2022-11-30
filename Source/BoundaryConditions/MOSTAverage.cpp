#include <MOSTAverage.H>

// Constructor
MOSTAverage::MOSTAverage (const amrex::Vector<amrex::Geometry>& geom,
                          amrex::Vector<amrex::Vector<amrex::MultiFab>>& vars_old,
                          amrex::Vector<std::unique_ptr<amrex::MultiFab>>& Theta_prim,
                          amrex::Vector<std::unique_ptr<amrex::MultiFab>>& z_phys_nd)
  : m_geom(geom)
{
    // Get basic info
    //--------------------------------------------------------
    amrex::ParmParse pp(m_pp_prefix);
    pp.query("most.radius",m_radius);
    pp.query("most.time_average",m_t_avg);
    pp.query("most.average_policy",m_policy);
    pp.query("most.use_normal_vector",m_norm_vec);

    // Set up fields and 2D MF/iMFs for averages
    //--------------------------------------------------------
    m_maxlev = m_geom.size();
    m_fields.resize(m_maxlev);
    
    m_k_in.resize(m_maxlev);
    
    m_k_indx.resize(m_maxlev);
    m_j_indx.resize(m_maxlev);
    m_i_indx.resize(m_maxlev);
    
    m_averages.resize(m_maxlev);
    
    m_z_phys_nd.resize(m_maxlev);
    for (int lev(0); lev < m_maxlev; lev++) {
      m_fields[lev].resize(m_nvar);
      m_averages[lev].resize(m_navg);
      m_z_phys_nd[lev] = z_phys_nd[lev].get();
      { // Nodal in x
        auto& mf  = vars_old[lev][Vars::xvel];
        amrex::MultiFab* mfp = &vars_old[lev][Vars::xvel];
        // Create a 2D ba, dm, & ghost cells
        amrex::BoxArray ba  = mf.boxArray();
        amrex::BoxList bl2d = ba.boxList();
        for (auto& b : bl2d) b.setRange(2,0);
        amrex::BoxArray ba2d(std::move(bl2d));
        const amrex::DistributionMapping& dm = mf.DistributionMap();
        const int ncomp   = 1;
        amrex::IntVect ng = mf.nGrowVect(); ng[2]=0;

          m_fields[lev][0] = mfp;
        m_averages[lev][0] = new amrex::MultiFab(ba2d,dm,ncomp,ng);
        m_averages[lev][0]->setVal(1.E34);
      }
      { // Nodal in y
        auto& mf  = vars_old[lev][Vars::yvel];
        amrex::MultiFab* mfp = &vars_old[lev][Vars::yvel];
        // Create a 2D ba, dm, & ghost cells
        amrex::BoxArray ba  = mf.boxArray();
        amrex::BoxList bl2d = ba.boxList();
        for (auto& b : bl2d) b.setRange(2,0);
        amrex::BoxArray ba2d(std::move(bl2d));
        const amrex::DistributionMapping& dm = mf.DistributionMap();
        const int ncomp   = 1;
        amrex::IntVect ng = mf.nGrowVect(); ng[2]=0;

          m_fields[lev][1] = mfp;
        m_averages[lev][1] = new amrex::MultiFab(ba2d,dm,ncomp,ng);
        m_averages[lev][1]->setVal(1.E34);
      }
      { // CC vars
        auto& mf  = *Theta_prim[lev];
        amrex::MultiFab* mfp = Theta_prim[lev].get();
        // Create a 2D ba, dm, & ghost cells
        amrex::BoxArray ba  = mf.boxArray();
        amrex::BoxList bl2d = ba.boxList();
        for (auto& b : bl2d) b.setRange(2,0);
        amrex::BoxArray ba2d(std::move(bl2d));
        const amrex::DistributionMapping& dm = mf.DistributionMap();
        const int ncomp   = 1;
        const int incomp  = 1;
        amrex::IntVect ng = mf.nGrowVect(); ng[2]=0;

          m_fields[lev][2] = mfp;
        m_averages[lev][2] = new amrex::MultiFab(ba2d,dm,ncomp,ng);
        m_averages[lev][2]->setVal(1.E34);

        m_averages[lev][3] = new amrex::MultiFab(ba2d,dm,ncomp,ng);
        m_averages[lev][3]->setVal(1.E34);

        m_k_indx[lev] = new amrex::iMultiFab(ba2d,dm,incomp,ng);
        
        if (m_norm_vec) {
            m_j_indx[lev] = new amrex::iMultiFab(ba2d,dm,incomp,ng);
            m_i_indx[lev] = new amrex::iMultiFab(ba2d,dm,incomp,ng);
        } else {
            m_j_indx[lev] = nullptr;
            m_i_indx[lev] = nullptr;
        }
      }
    } // lev

    // Setup 2d iMF for spatial configuration
    //--------------------------------------------------------
    if (m_norm_vec && m_z_phys_nd[0]) {
        set_ijk_indices_T();
    } else if (m_z_phys_nd[0]) {
        set_k_indices_T();
    } else {
        set_k_indices_N();
    }

    // Setup auxiliary data for the chosen policy
    //--------------------------------------------------------
    switch(m_policy) {
    case 0: // Standard plane average
        set_plane_normalization();
        break;

    case 1: // Local region/point
        set_region_normalization();
        break;

    default:
        AMREX_ASSERT_WITH_MESSAGE(false, "Unknown policy for MOSTAverage!");
    }

    // Set up the exponential time filtering
    //--------------------------------------------------------
    if (m_t_avg) {
        // m_time_window is normalized by the time-step "dt"
        pp.query("most.time_window", m_time_window);

        // Exponential filter function
        m_fact_old = std::exp(-1.0 / m_time_window);

        // Enforce discrete normalization: (mfn*val_new + mfo*val_old)
        m_fact_new = 1.0 - m_fact_old;

        // None of the averages are initialized
        m_t_init.resize(m_maxlev,0);
    }
}


// Reset the pointers to field MFs
void
MOSTAverage::update_field_ptrs(int lev,
                               amrex::Vector<amrex::Vector<amrex::MultiFab>>& vars_old,
                               amrex::Vector<std::unique_ptr<amrex::MultiFab>>& Theta_prim)
{
    m_fields[lev][0] = &vars_old[lev][Vars::xvel];
    m_fields[lev][1] = &vars_old[lev][Vars::yvel];
    m_fields[lev][2] = Theta_prim[lev].get();
}

// Compute ncells per plane
void
MOSTAverage::set_plane_normalization()
{
    // Cells per plane and temp avg storage
    m_ncell_plane.resize(m_maxlev);
    m_plane_average.resize(m_maxlev);

    for (int lev(0); lev < m_maxlev; lev++) {
        // Num components, plane avg, cells per plane
        amrex::Box domain = m_geom[lev].Domain();
        amrex::IntVect dom_lo(domain.loVect());
        amrex::IntVect dom_hi(domain.hiVect());
        m_ncell_plane[lev].resize(m_navg);
        m_plane_average[lev].resize(m_navg);
        for (int iavg(0); iavg < m_navg; ++iavg) {
            m_plane_average[lev][iavg] = 0.0;

            m_ncell_plane[lev][iavg] = 1;
            amrex::IndexType ixt = m_averages[lev][iavg]->boxArray().ixType();
            for (int idim(0); idim < AMREX_SPACEDIM; ++idim) {
                if (idim != 2) {
                    if (ixt.nodeCentered(idim)) {
                        m_ncell_plane[lev][iavg] *= (dom_hi[idim] - dom_lo[idim] + 2);
                    } else {
                        m_ncell_plane[lev][iavg] *= (dom_hi[idim] - dom_lo[idim] + 1);
                    }
                }
            } // idim
        } // iavg
    } // lev
}


// Populate a 2D iMF with the k indices for averaging
void
MOSTAverage::set_k_indices_N()
{
    amrex::ParmParse pp(m_pp_prefix);
    auto read_z = pp.query("most.zref",m_zref);
    auto read_k = pp.queryarr("most.k_arr_in",m_k_in);

    // Specify z_ref & compute k_indx (z_ref takes precedence)
    if (read_z) {
        for (int lev(0); lev < m_maxlev; lev++) {
            amrex::Real m_zlo = m_geom[lev].ProbLo(2);
            amrex::Real m_dz  = m_geom[lev].CellSize(2);
            
            AMREX_ALWAYS_ASSERT(m_zref >= m_zlo + 0.5 * m_dz);
            
            int lk = static_cast<int>(floor((m_zref - m_zlo) / m_dz - 0.5));
            lk = std::max(m_radius,lk);
            
            for (amrex::MFIter mfi(*m_k_indx[lev], amrex::TilingIfNotGPU()); mfi.isValid(); ++mfi) {
                amrex::Box gpbx = mfi.tilebox();
                auto k_arr      = m_k_indx[lev]->array(mfi);
                ParallelFor(gpbx, [=] AMREX_GPU_DEVICE (int i, int j, int k) noexcept
                {
                    k_arr(i,j,k) = lk;
                });
            }
        }
    // Specified k_indx & compute z_ref
    } else if (read_k) {
        for (int lev(0); lev < m_maxlev; lev++){
            m_k_in[lev] = std::max(m_radius,m_k_in[lev]);
            m_k_indx[lev]->setVal(m_k_in[lev]);
        }
        
        // TODO: check that z_ref is constant across levels
        amrex::Real m_zlo = m_geom[0].ProbLo(2);
        amrex::Real m_dz  = m_geom[0].CellSize(2);
        m_zref = ((amrex::Real)m_k_in[0] + 0.5) * m_dz + m_zlo;
    }
}


// Populate a 2D iMF with the k indices for averaging
void
MOSTAverage::set_k_indices_T()
{
    amrex::ParmParse pp(m_pp_prefix);
    auto read_z = pp.query("most.zref",m_zref);
    auto read_k = pp.queryarr("most.k_arr_in",m_k_in);

    // Specify z_ref & compute k_indx (z_ref takes precedence)
    if (read_z) {
        for (int lev(0); lev < m_maxlev; lev++) {
            int kmax = m_geom[lev].Domain().bigEnd(2);
            for (amrex::MFIter mfi(*m_k_indx[lev], amrex::TilingIfNotGPU()); mfi.isValid(); ++mfi) {
                amrex::Box gpbx = mfi.tilebox();
                auto z_arr      = m_z_phys_nd[lev]->array(mfi);
                auto k_arr      = m_k_indx[lev]->array(mfi);
                ParallelFor(gpbx, [=] AMREX_GPU_DEVICE (int i, int j, int k) noexcept
                {
                    k_arr(i,j,k) = 0;
                    
                    for (int lk(0); lk<=kmax; ++lk) {
                        if (m_zref > z_arr(i,j,lk) && m_zref < z_arr(i,j,lk+1)){
                            k_arr(i,j,k) = lk;
                            break;
                        }
                    }
                });
            }
        }   
    // Specified k_indx & compute z_ref
    } else if (read_k) {
        AMREX_ASSERT_WITH_MESSAGE(false, "Specified k-indx with terrain not implemented!");
    }
}


// Populate all 2D iMFs for averaging
void
MOSTAverage::set_ijk_indices_T()
{

}


// Driver to call appropriate average member function
void
MOSTAverage::compute_averages(int lev)
{
    switch(m_policy) {
    case 0: // Standard plane average
        compute_plane_averages(lev);
        break;

    case 1: // Local region/point
        compute_region_averages(lev);
        break;

    default:
        AMREX_ASSERT_WITH_MESSAGE(false, "Unknown policy for MOSTAverage!");
    }

    // We have initialized the averages
    if (m_t_avg) m_t_init[lev] = 1;
}


// Fill plane storage with averages
void
MOSTAverage::compute_plane_averages(int lev)
{
    // Peel back the level
    auto& fields        = m_fields[lev];
    auto& averages      = m_averages[lev];
    auto& k_indx        = m_k_indx[lev];
    auto& j_indx        = m_j_indx[lev];
    auto& i_indx        = m_i_indx[lev];
    auto& ncell_plane   = m_ncell_plane[lev];
    auto& plane_average = m_plane_average[lev];

    // Set factors for time averaging
    amrex::Real d_fact_new, d_fact_old;
    if (m_t_avg && m_t_init[lev]) {
        d_fact_new = m_fact_new;
        d_fact_old = m_fact_old;
    } else {
        d_fact_new = 1.0;
        d_fact_old = 0.0;
    }

    // GPU array to accumulate averages into
    amrex::Gpu::DeviceVector<amrex::Real> pavg(plane_average.size(), 0.0);
    amrex::Real* plane_avg = pavg.data();

    // Averages over all the fields
    //----------------------------------------------------------
    for (int imf(0); imf < m_nvar; ++imf) {
        const amrex::Real denom = 1.0 / (amrex::Real)ncell_plane[imf];

#ifdef _OPENMP
#pragma omp parallel if (amrex::Gpu::notInLaunchRegion())
#endif
        for (amrex::MFIter mfi(*fields[imf], amrex::TilingIfNotGPU()); mfi.isValid(); ++mfi) {
            amrex::Box pbx = mfi.tilebox(); pbx.setSmall(2,0); pbx.setBig(2,0);

            auto mf_arr = fields[imf]->const_array(mfi);
            auto k_arr  = k_indx->const_array(mfi);
            auto j_arr  = j_indx ? j_indx->const_array(mfi) : amrex::Array4<const int> {};
            auto i_arr  = i_indx ? i_indx->const_array(mfi) : amrex::Array4<const int> {};
            
            amrex::Real d_val_old = plane_average[imf]*d_fact_old;

            ParallelFor(amrex::Gpu::KernelInfo().setReduction(true), pbx, [=]
            AMREX_GPU_DEVICE(int i, int j, int k, amrex::Gpu::Handler const& handler) noexcept
            {
                int mk = k_arr(i,j,k);
                int mj = j_arr ? j_arr(i,j,k) : j;
                int mi = i_arr ? i_arr(i,j,k) : i;
                amrex::Real val = denom * ( mf_arr(mi,mj,mk)*d_fact_new + d_val_old );
                amrex::Gpu::deviceReduceSum(&plane_avg[imf], val, handler);
            });
        }
    }

    // Averages for the tangential velocity magnitude
    //----------------------------------------------------------
    {
        int imf  = 0;
        int iavg = m_navg - 1;
        const amrex::Real denom = 1.0 / (amrex::Real)ncell_plane[iavg];

#ifdef _OPENMP
#pragma omp parallel if (amrex::Gpu::notInLaunchRegion())
#endif
        for (amrex::MFIter mfi(*averages[iavg], amrex::TilingIfNotGPU()); mfi.isValid(); ++mfi) {
            amrex::Box pbx = mfi.tilebox(); pbx.setSmall(2,0); pbx.setBig(2,0);

            auto u_mf_arr = fields[imf  ]->const_array(mfi);
            auto v_mf_arr = fields[imf+1]->const_array(mfi);
            auto k_arr    = k_indx->const_array(mfi);
            auto j_arr    = j_indx ? j_indx->const_array(mfi) : amrex::Array4<const int> {};
            auto i_arr    = i_indx ? i_indx->const_array(mfi) : amrex::Array4<const int> {};

            amrex::Real d_val_old = plane_average[iavg]*d_fact_old;

            ParallelFor(amrex::Gpu::KernelInfo().setReduction(true), pbx, [=]
            AMREX_GPU_DEVICE(int i, int j, int k, amrex::Gpu::Handler const& handler) noexcept
            {
                int mk = k_arr(i,j,k);
                int mj = j_arr ? j_arr(i,j,k) : j;
                int mi = i_arr ? i_arr(i,j,k) : i;
                
                const amrex::Real u_val = 0.5 * (u_mf_arr(mi,mj,mk) + u_mf_arr(mi+1,mj  ,mk));
                const amrex::Real v_val = 0.5 * (v_mf_arr(mi,mj,mk) + v_mf_arr(mi  ,mj+1,mk));
                const amrex::Real mag   = std::sqrt(u_val*u_val + v_val*v_val);
                amrex::Real val = denom * ( mag*d_fact_new + d_val_old);
                amrex::Gpu::deviceReduceSum(&plane_avg[iavg], val, handler);
            });
        }
    }

    // Copy to host and sum across procs
    amrex::Gpu::copy(amrex::Gpu::deviceToHost, pavg.begin(), pavg.end(), plane_average.begin());
    amrex::ParallelDescriptor::ReduceRealSum(plane_average.data(), plane_average.size());

    // No spatial variation with plane averages
    for (int iavg(0); iavg < m_navg; ++iavg) averages[iavg]->setVal(plane_average[iavg]);
}


// Fill 2D MF with local averages
void
MOSTAverage::compute_region_averages(int lev)
{
    // Peel back the level
    auto& fields   = m_fields[lev];
    auto& averages = m_averages[lev];
    auto& k_indx   = m_k_indx[lev];
    auto& j_indx   = m_j_indx[lev];
    auto& i_indx   = m_i_indx[lev];
    auto& geom     = m_geom[lev];

    // Set factors for time averaging
    amrex::Real d_fact_new, d_fact_old;
    if (m_t_avg && m_t_init[lev]) {
        d_fact_new = m_fact_new;
        d_fact_old = m_fact_old;
    } else {
        d_fact_new = 1.0;
        d_fact_old = 0.0;
    }

    // Number of cells contained in the local average
    const amrex::Real denom = 1.0 / (amrex::Real) m_ncell_region;

    // Capture radius for device
    int d_radius = m_radius;

    // Averages over all the fields
    //----------------------------------------------------------
    for (int imf(0); imf < m_nvar; ++imf) {
#ifdef _OPENMP
#pragma omp parallel if (amrex::Gpu::notInLaunchRegion())
#endif
        for (amrex::MFIter mfi(*fields[imf], amrex::TilingIfNotGPU()); mfi.isValid(); ++mfi) {
            amrex::Box pbx = mfi.tilebox(); pbx.setSmall(2,0); pbx.setBig(2,0);

            auto mf_arr = fields[imf]->const_array(mfi);
            auto ma_arr = averages[imf]->array(mfi);
            auto k_arr  = k_indx->const_array(mfi);
            auto j_arr  = j_indx ? j_indx->const_array(mfi) : amrex::Array4<const int> {};
            auto i_arr  = i_indx ? i_indx->const_array(mfi) : amrex::Array4<const int> {};

            ParallelFor(pbx, [=] AMREX_GPU_DEVICE(int i, int j, int k) noexcept
            {
                ma_arr(i,j,k) *= d_fact_old;

                int mk = k_arr(i,j,k);
                int mj = j_arr ? j_arr(i,j,k) : j;
                int mi = i_arr ? i_arr(i,j,k) : i;
                for (int lk(mk-d_radius); lk <= (mk+d_radius); ++lk) {
                    for (int lj(mj-d_radius); lj <= (mj+d_radius); ++lj) {
                        for (int li(mi-d_radius); li <= (mi+d_radius); ++li) {
                            amrex::Real val = denom * mf_arr(li, lj, lk) * d_fact_new;
                            ma_arr(i,j,k) += val;
                        }
                    }
                }
            });
        }

        // Fill interior ghost cells and any ghost cells outside a periodic domain
        //***********************************************************************************
        averages[imf]->FillBoundary(geom.periodicity());
    }

    // Averages for the tangential velocity magnitude
    //----------------------------------------------------------
    {
        int imf  = 0;
        int iavg = m_navg - 1;

#ifdef _OPENMP
#pragma omp parallel if (amrex::Gpu::notInLaunchRegion())
#endif
        for (amrex::MFIter mfi(*averages[iavg], amrex::TilingIfNotGPU()); mfi.isValid(); ++mfi) {
            amrex::Box pbx = mfi.tilebox(); pbx.setSmall(2,0); pbx.setBig(2,0);

            auto u_mf_arr = fields[imf]->const_array(mfi);
            auto v_mf_arr = fields[imf+1]->const_array(mfi);
            auto ma_arr   = averages[iavg]->array(mfi);
            auto k_arr    = k_indx->const_array(mfi);
            auto j_arr    = j_indx ? j_indx->const_array(mfi) : amrex::Array4<const int> {};
            auto i_arr    = i_indx ? i_indx->const_array(mfi) : amrex::Array4<const int> {};
            
            ParallelFor(pbx, [=] AMREX_GPU_DEVICE(int i, int j, int k) noexcept
            {
                ma_arr(i,j,k) *= d_fact_old;

                int mk = k_arr(i,j,k);
                int mj = j_arr ? j_arr(i,j,k) : j;
                int mi = i_arr ? i_arr(i,j,k) : i;
                for (int lk(mk-d_radius); lk <= (mk+d_radius); ++lk) {
                    for (int lj(mj-d_radius); lj <= (mj+d_radius); ++lj) {
                        for (int li(mi-d_radius); li <= (mi+d_radius); ++li) {
                            const amrex::Real u_val = 0.5 * (u_mf_arr(li,lj,lk) + u_mf_arr(li+1,lj  ,lk));
                            const amrex::Real v_val = 0.5 * (v_mf_arr(li,lj,lk) + v_mf_arr(li  ,lj+1,lk));
                            const amrex::Real mag   = std::sqrt(u_val*u_val + v_val*v_val);
                            amrex::Real val = denom * mag * d_fact_new;
                            ma_arr(i,j,k) += val;
                        }
                    }
                }
            });

            // Fill interior ghost cells and any ghost cells outside a periodic domain
            //***********************************************************************************
            averages[iavg]->FillBoundary(geom.periodicity());
        }
    }


    // Need to fill ghost cells outside the domain if not periodic
    bool not_per_x = !(geom.periodicity().isPeriodic(0));
    bool not_per_y = !(geom.periodicity().isPeriodic(1));
    if (not_per_x || not_per_y) {
        amrex::Box domain = geom.Domain();
        for (int iavg(0); iavg < m_navg; ++iavg) {
#ifdef _OPENMP
#pragma omp parallel if (amrex::Gpu::notInLaunchRegion())
#endif
            amrex::IndexType ixt = averages[iavg]->boxArray().ixType();
            amrex::Box ldomain   = domain; ldomain.convert(ixt);
            amrex::IntVect ng    = averages[iavg]->nGrowVect(); ng[2]=0;
            for (amrex::MFIter mfi(*averages[iavg], amrex::TilingIfNotGPU()); mfi.isValid(); ++mfi) {
                amrex::Box gpbx = mfi.growntilebox(ng); gpbx.setSmall(2,0); gpbx.setBig(2,0);

                if (ldomain.contains(gpbx)) continue;

                auto ma_arr = averages[iavg]->array(mfi);

                int i_lo = ldomain.smallEnd(0); int i_hi = ldomain.bigEnd(0);
                int j_lo = ldomain.smallEnd(1); int j_hi = ldomain.bigEnd(1);
                ParallelFor(gpbx, [=] AMREX_GPU_DEVICE(int i, int j, int k) noexcept
                {
                    int li, lj;
                    li = i  < i_lo ? i_lo : i;
                    li = li > i_hi ? i_hi : li;
                    lj = j  < j_lo ? j_lo : j;
                    lj = lj > j_hi ? j_hi : lj;

                    ma_arr(i,j,k) = ma_arr(li,lj,k);
                });
            } // MFiter
        } // iavg
    } // Not periodic
}


// Write k indices
void
MOSTAverage::write_k_indices(int lev)
{
    // Peel back the level
    auto& averages = m_averages[lev];
    auto& k_indx   = m_k_indx[lev];

    int navg = m_navg - 1;

    std::ofstream ofile;
    ofile.open ("MOST_k_indices.txt");
    ofile << "K indices used to compute averages via MOSTAverages class:\n";

    for (amrex::MFIter mfi(*averages[navg], amrex::TilingIfNotGPU()); mfi.isValid(); ++mfi) {
        amrex::Box bx  = mfi.tilebox(); bx.setBig(2,0);
        int il = bx.smallEnd(0); int iu = bx.bigEnd(0);
        int jl = bx.smallEnd(1); int ju = bx.bigEnd(1);

        auto k_arr = k_indx->array(mfi);

        for (int j(jl); j <= ju; ++j) {
            for (int i(il); i <= iu; ++i) {
                ofile << "(I,J): " << "(" << i << "," << j << ")" << "\n";
                int k = 0;
                ofile << "K_ind: "
                      << k_arr(i,j,k) << "\n";
                ofile << "\n";
            }
        }
    }
    ofile.close();
}

// Write averages
void
MOSTAverage::write_averages(int lev)
{
    // Peel back the level
    auto& averages = m_averages[lev];

    int navg = m_navg - 1;

    std::ofstream ofile;
    ofile.open ("MOST_averages.txt");
    ofile << "Averages computed via MOSTAverages class:\n";

    for (amrex::MFIter mfi(*averages[navg], amrex::TilingIfNotGPU()); mfi.isValid(); ++mfi) {
        amrex::Box bx  = mfi.tilebox(); bx.setBig(2,0);
        int il = bx.smallEnd(0); int iu = bx.bigEnd(0);
        int jl = bx.smallEnd(1); int ju = bx.bigEnd(1);

        for (int j(jl); j <= ju; ++j) {
            for (int i(il); i <= iu; ++i) {
                ofile << "(I,J): " << "(" << i << "," << j << ")" << "\n";
                int k = 0;
                for (int iavg(0); iavg <= navg; ++iavg) {
                    auto mf_arr = averages[iavg]->array(mfi);
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
