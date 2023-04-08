#include <MOSTAverage.H>
#include <utility>
#include <TileNoZ.H>

using namespace amrex;

// Constructor
MOSTAverage::MOSTAverage (Vector<Geometry>  geom,
                          Vector<Vector<MultiFab>>& vars_old,
                          Vector<std::unique_ptr<MultiFab>>& Theta_prim,
                          Vector<std::unique_ptr<MultiFab>>& z_phys_nd)
  : m_geom(std::move(geom))
{
    // Get basic info
    //--------------------------------------------------------
    ParmParse pp(m_pp_prefix);
    pp.query("most.radius",m_radius);
    pp.query("most.time_average",m_t_avg);
    pp.query("most.average_policy",m_policy);
    pp.query("most.use_interpolation",m_interp);
    pp.query("most.use_normal_vector",m_norm_vec);

    AMREX_ASSERT_WITH_MESSAGE(m_radius<=2, "Radius must be less than nGhost=3!");
    if (m_interp) AMREX_ASSERT_WITH_MESSAGE((z_phys_nd[0].get()), "Interpolation only implemented with terrain!");

    // Set up fields and 2D MF/iMFs for averages
    //--------------------------------------------------------
    m_maxlev = m_geom.size();

    m_fields.resize(m_maxlev);
    m_averages.resize(m_maxlev);
    m_z_phys_nd.resize(m_maxlev);

    m_k_in.resize(m_maxlev);

    m_x_pos.resize(m_maxlev);
    m_y_pos.resize(m_maxlev);
    m_z_pos.resize(m_maxlev);

    m_i_indx.resize(m_maxlev);
    m_j_indx.resize(m_maxlev);
    m_k_indx.resize(m_maxlev);


    for (int lev(0); lev < m_maxlev; lev++) {
      m_fields[lev].resize(m_nvar);
      m_averages[lev].resize(m_navg);
      m_z_phys_nd[lev] = z_phys_nd[lev].get();
      { // Nodal in x
        auto& mf  = vars_old[lev][Vars::xvel];
        MultiFab* mfp = &vars_old[lev][Vars::xvel];
        // Create a 2D ba, dm, & ghost cells
        const BoxArray& ba  = mf.boxArray();
        BoxList bl2d = ba.boxList();
        for (auto& b : bl2d) b.setRange(2,0);
        BoxArray ba2d(std::move(bl2d));
        const DistributionMapping& dm = mf.DistributionMap();
        const int ncomp   = 1;
        IntVect ng = mf.nGrowVect(); ng[2]=0;

          m_fields[lev][0] = mfp;
        m_averages[lev][0] = new MultiFab(ba2d,dm,ncomp,ng);
        m_averages[lev][0]->setVal(1.E34);
      }
      { // Nodal in y
        auto& mf  = vars_old[lev][Vars::yvel];
        MultiFab* mfp = &vars_old[lev][Vars::yvel];
        // Create a 2D ba, dm, & ghost cells
        const BoxArray& ba  = mf.boxArray();
        BoxList bl2d = ba.boxList();
        for (auto& b : bl2d) b.setRange(2,0);
        BoxArray ba2d(std::move(bl2d));
        const DistributionMapping& dm = mf.DistributionMap();
        const int ncomp   = 1;
        IntVect ng = mf.nGrowVect(); ng[2]=0;

          m_fields[lev][1] = mfp;
        m_averages[lev][1] = new MultiFab(ba2d,dm,ncomp,ng);
        m_averages[lev][1]->setVal(1.E34);
      }
      { // CC vars
        auto& mf  = *Theta_prim[lev];
        MultiFab* mfp = Theta_prim[lev].get();
        // Create a 2D ba, dm, & ghost cells
        const BoxArray& ba  = mf.boxArray();
        BoxList bl2d = ba.boxList();
        for (auto& b : bl2d) b.setRange(2,0);
        BoxArray ba2d(std::move(bl2d));
        const DistributionMapping& dm = mf.DistributionMap();
        const int ncomp   = 1;
        const int incomp  = 1;
        IntVect ng = mf.nGrowVect(); ng[2]=0;

          m_fields[lev][2] = mfp;
        m_averages[lev][2] = new MultiFab(ba2d,dm,ncomp,ng);
        m_averages[lev][2]->setVal(1.E34);

        m_averages[lev][3] = new MultiFab(ba2d,dm,ncomp,ng);
        m_averages[lev][3]->setVal(1.E34);

        if (m_z_phys_nd[0] && m_norm_vec && m_interp) {
            m_x_pos[lev] = new MultiFab(ba2d,dm,ncomp,ng);
            m_y_pos[lev] = new MultiFab(ba2d,dm,ncomp,ng);
            m_z_pos[lev] = new MultiFab(ba2d,dm,ncomp,ng);
        } else if (m_z_phys_nd[0] && m_interp) {
            m_x_pos[lev] = new MultiFab(ba2d,dm,ncomp,ng);
            m_y_pos[lev] = new MultiFab(ba2d,dm,ncomp,ng);
            m_z_pos[lev] = new MultiFab(ba2d,dm,ncomp,ng);
        } else if (m_z_phys_nd[0] && m_norm_vec) {
            m_i_indx[lev] = new iMultiFab(ba2d,dm,incomp,ng);
            m_j_indx[lev] = new iMultiFab(ba2d,dm,incomp,ng);
            m_k_indx[lev] = new iMultiFab(ba2d,dm,incomp,ng);
        } else {
            m_k_indx[lev] = new iMultiFab(ba2d,dm,incomp,ng);
        }
      }
    } // lev

    // Setup auxiliary data for spatial configuration & policy
    //--------------------------------------------------------
    if (m_z_phys_nd[0] && m_norm_vec && m_interp) { // Terrain w/ norm & w/ interpolation
        set_norm_positions_T();
    } else if (m_z_phys_nd[0] && m_interp) {        // Terrain w/ interpolation
        set_z_positions_T();
    } else if (m_z_phys_nd[0] && m_norm_vec) {      // Terrain w/ norm & w/o interpolation
        set_norm_indices_T();
    } else if (m_z_phys_nd[0]) {                    // Terrain
        set_k_indices_T();
    } else {                                        // No Terrain
        set_k_indices_N();
    }

    // Setup normalization data for the chosen policy
    //--------------------------------------------------------
    switch(m_policy) {
    case 0: // Plane average
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
                               Vector<Vector<MultiFab>>& vars_old,
                               Vector<std::unique_ptr<MultiFab>>& Theta_prim)
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
        Array<int,AMREX_SPACEDIM> is_per = {0,0,0};
        for (int idim(0); idim < AMREX_SPACEDIM-1; ++idim) {
            if (m_geom[lev].isPeriodic(idim)) is_per[idim] = 1;
        }
        Box domain = m_geom[lev].Domain();
        m_ncell_plane[lev].resize(m_navg);
        m_plane_average[lev].resize(m_navg);
        for (int iavg(0); iavg < m_navg; ++iavg) {
            // Convert domain to current index type
            IndexType ixt = m_averages[lev][iavg]->boxArray().ixType();
            domain.convert(ixt);
            IntVect dom_lo(domain.loVect());
            IntVect dom_hi(domain.hiVect());

            m_plane_average[lev][iavg] = 0.0;

            m_ncell_plane[lev][iavg] = 1;
            for (int idim(0); idim < AMREX_SPACEDIM; ++idim) {
                if (idim != 2) {
                    if (ixt.nodeCentered(idim) && is_per[idim]) {
                        m_ncell_plane[lev][iavg] *= (dom_hi[idim] - dom_lo[idim]);
                    } else {
                        m_ncell_plane[lev][iavg] *= (dom_hi[idim] - dom_lo[idim] + 1);
                    }
                }
            } // idim
        } // iavg
    } // lev

}


// Populate a 2D iMF with the k indices for averaging (w/o terrain)
void
MOSTAverage::set_k_indices_N()
{
    ParmParse pp(m_pp_prefix);
    auto read_z = pp.query("most.zref",m_zref);
    auto read_k = pp.queryarr("most.k_arr_in",m_k_in);

    // Specify z_ref & compute k_indx (z_ref takes precedence)
    if (read_z) {
        for (int lev(0); lev < m_maxlev; lev++) {
            Real m_zlo = m_geom[lev].ProbLo(2);
            Real m_dz  = m_geom[lev].CellSize(2);

            AMREX_ASSERT_WITH_MESSAGE(m_zref >= m_zlo + 0.5 * m_dz,
                                      "Query point must be past first z-cell!");

            int lk = static_cast<int>(floor((m_zref - m_zlo) / m_dz - 0.5));

            AMREX_ALWAYS_ASSERT(lk >= m_radius);

            m_k_indx[lev]->setVal(lk);
        }
    // Specified k_indx & compute z_ref
    } else if (read_k) {
        for (int lev(0); lev < m_maxlev; lev++){
            AMREX_ASSERT_WITH_MESSAGE(m_k_in[lev] >= m_radius,
                                      "K index must be larger than averaging radius!");
            m_k_indx[lev]->setVal(m_k_in[lev]);
        }

        // TODO: check that z_ref is constant across levels
        Real m_zlo = m_geom[0].ProbLo(2);
        Real m_dz  = m_geom[0].CellSize(2);
        m_zref = ((Real)m_k_in[0] + 0.5) * m_dz + m_zlo;
    }
}


// Populate a 2D iMF with the k indices for averaging (w/ terrain)
void
MOSTAverage::set_k_indices_T()
{
    ParmParse pp(m_pp_prefix);
    auto read_z = pp.query("most.zref",m_zref);
    auto read_k = pp.queryarr("most.k_arr_in",m_k_in);

    // Capture for device
    Real d_zref   = m_zref;
    Real d_radius = m_radius;

    // Specify z_ref & compute k_indx (z_ref takes precedence)
    if (read_z) {
        for (int lev(0); lev < m_maxlev; lev++) {
            int kmax = m_geom[lev].Domain().bigEnd(2);
            for (MFIter mfi(*m_k_indx[lev], TileNoZ()); mfi.isValid(); ++mfi) {
                Box npbx  = mfi.tilebox(); npbx.convert({1,1,0});
                const auto z_phys_arr = m_z_phys_nd[lev]->const_array(mfi);
                auto k_arr = m_k_indx[lev]->array(mfi);
                ParallelFor(npbx, [=] AMREX_GPU_DEVICE (int i, int j, int k) noexcept
                {
                    Real z_target = d_zref + z_phys_arr(i,j,k);
                    for (int lk(0); lk<=kmax; ++lk) {
                        Real z_lo = 0.25 * ( z_phys_arr(i,j  ,lk  ) + z_phys_arr(i+1,j  ,lk  )
                                                  + z_phys_arr(i,j+1,lk  ) + z_phys_arr(i+1,j+1,lk  ) );
                        Real z_hi = 0.25 * ( z_phys_arr(i,j  ,lk+1) + z_phys_arr(i+1,j  ,lk+1)
                                                  + z_phys_arr(i,j+1,lk+1) + z_phys_arr(i+1,j+1,lk+1) );
                        if (z_target > z_lo && z_target < z_hi){
                            AMREX_ASSERT_WITH_MESSAGE(lk >= d_radius,
                                                      "K index must be larger than averaging radius!");
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


// Populate all 2D iMFs for averaging (w/ terrain & norm vector)
void
MOSTAverage::set_norm_indices_T()
{
    ParmParse pp(m_pp_prefix);
    pp.query("most.zref",m_zref);

    // Capture for device
    Real d_zref   = m_zref;
    Real d_radius = m_radius;

    for (int lev(0); lev < m_maxlev; lev++) {
        int kmax = m_geom[lev].Domain().bigEnd(2);
        const auto dxInv  = m_geom[lev].InvCellSizeArray();
        IntVect ng = m_k_indx[lev]->nGrowVect(); ng[2]=0;
        for (MFIter mfi(*m_k_indx[lev], TileNoZ()); mfi.isValid(); ++mfi) {
            Box npbx  = mfi.tilebox(); npbx.convert({1,1,0});
            Box gpbx  = mfi.growntilebox(ng);
            const auto z_phys_arr = m_z_phys_nd[lev]->const_array(mfi);
            auto i_arr = m_i_indx[lev]->array(mfi);
            auto j_arr = m_j_indx[lev]->array(mfi);
            auto k_arr = m_k_indx[lev]->array(mfi);
            ParallelFor(npbx, [=] AMREX_GPU_DEVICE (int i, int j, int k) noexcept
            {
                // Elements of normal vector
                Real met_h_xi  = Compute_h_xi_AtCellCenter (i,j,k,dxInv,z_phys_arr);
                Real met_h_eta = Compute_h_eta_AtCellCenter(i,j,k,dxInv,z_phys_arr);
                Real mag = std::sqrt(met_h_xi*met_h_xi + met_h_eta*met_h_eta + 1.0);

                // Unit-normal vector scaled by z_ref
                Real delta_x = -met_h_xi/mag  * d_zref;
                Real delta_y = -met_h_eta/mag * d_zref;
                Real delta_z = 1.0/mag * d_zref;

                // Compute i & j as displacements (no grid stretching)
                int delta_i  = static_cast<int>(std::round(delta_x*dxInv[0]));
                int delta_j  = static_cast<int>(std::round(delta_y*dxInv[1]));
                int i_new    = i + delta_i;
                int j_new    = j + delta_j;
                i_arr(i,j,k) = i_new;
                j_arr(i,j,k) = j_new;

                // Search for k (grid is stretched in z)
                Real z_target = delta_z + z_phys_arr(i,j,k);
                for (int lk(0); lk<=kmax; ++lk) {
                    Real z_lo = 0.25 * ( z_phys_arr(i_new,j_new  ,lk  ) + z_phys_arr(i_new+1,j_new  ,lk  )
                                              + z_phys_arr(i_new,j_new+1,lk  ) + z_phys_arr(i_new+1,j_new+1,lk  ) );
                    Real z_hi = 0.25 * ( z_phys_arr(i_new,j_new  ,lk+1) + z_phys_arr(i_new+1,j_new  ,lk+1)
                                              + z_phys_arr(i_new,j_new+1,lk+1) + z_phys_arr(i_new+1,j_new+1,lk+1) );
                    if (z_target > z_lo && z_target < z_hi){
                        AMREX_ASSERT_WITH_MESSAGE(lk >= d_radius,
                                                  "K index must be larger than averaging radius!");
                        k_arr(i,j,k) = lk;
                        break;
                    }
                }

                // Destination cell must be contained on the current process!
                AMREX_ASSERT_WITH_MESSAGE(gpbx.contains(i_arr(i,j,k),j_arr(i,j,k),k_arr(i,j,k)),
                                          "Query index outside of proc domain!");
            });
        }
    }
}


// Populate positions (w/ terrain & interpolation)
void
MOSTAverage::set_z_positions_T()
{
    ParmParse pp(m_pp_prefix);
    pp.query("most.zref",m_zref);

    // Capture for device
    Real d_zref = m_zref;

    for (int lev(0); lev < m_maxlev; lev++) {
        RealVect base;
        const auto dx = m_geom[lev].CellSizeArray();
        IntVect ng = m_x_pos[lev]->nGrowVect(); ng[2]=0;
        for (MFIter mfi(*m_x_pos[lev], TileNoZ()); mfi.isValid(); ++mfi) {
            Box npbx  = mfi.tilebox(); npbx.convert({1,1,0});
            Box gpbx  = mfi.growntilebox(ng);
            RealBox grb{gpbx,dx.data(),base.dataPtr()};

            const auto z_phys_arr = m_z_phys_nd[lev]->const_array(mfi);
            auto x_pos_arr   = m_x_pos[lev]->array(mfi);
            auto y_pos_arr   = m_y_pos[lev]->array(mfi);
            auto z_pos_arr   = m_z_pos[lev]->array(mfi);
            ParallelFor(npbx, [=] AMREX_GPU_DEVICE (int i, int j, int k) noexcept
            {
                // Final position at end of vector
                x_pos_arr(i,j,k) = ((Real) i + 0.5) * dx[0];
                y_pos_arr(i,j,k) = ((Real) j + 0.5) * dx[1];
                z_pos_arr(i,j,k) = z_phys_arr(i,j,k) + d_zref;

                // Destination position must be contained on the current process!
                Real pos[] = {x_pos_arr(i,j,k),y_pos_arr(i,j,k),0.5*dx[2]};
                AMREX_ASSERT_WITH_MESSAGE( grb.contains(&pos[0]),
                                           "Query point outside of proc domain!");
            });
        }
    }
}


// Populate positions (w/ terrain & norm vector & interpolation)
void
MOSTAverage::set_norm_positions_T()
{
    ParmParse pp(m_pp_prefix);
    pp.query("most.zref",m_zref);

    // Capture for device
    Real d_zref = m_zref;

    for (int lev(0); lev < m_maxlev; lev++) {
        RealVect base;
        const auto dx = m_geom[lev].CellSizeArray();
        const auto dxInv  = m_geom[lev].InvCellSizeArray();
        IntVect ng = m_x_pos[lev]->nGrowVect(); ng[2]=0;
        for (MFIter mfi(*m_x_pos[lev], TileNoZ()); mfi.isValid(); ++mfi) {
            Box npbx  = mfi.tilebox(); npbx.convert({1,1,0});
            Box gpbx  = mfi.growntilebox(ng);
            RealBox grb{gpbx,dx.data(),base.dataPtr()};

            const auto z_phys_arr = m_z_phys_nd[lev]->const_array(mfi);
            auto x_pos_arr   = m_x_pos[lev]->array(mfi);
            auto y_pos_arr   = m_y_pos[lev]->array(mfi);
            auto z_pos_arr   = m_z_pos[lev]->array(mfi);
            ParallelFor(npbx, [=] AMREX_GPU_DEVICE (int i, int j, int k) noexcept
            {
                // Elements of normal vector
                Real met_h_xi  = Compute_h_xi_AtCellCenter (i,j,k,dxInv,z_phys_arr);
                Real met_h_eta = Compute_h_eta_AtCellCenter(i,j,k,dxInv,z_phys_arr);
                Real mag = std::sqrt(met_h_xi*met_h_xi + met_h_eta*met_h_eta + 1.0);

                // Unit-normal vector scaled by z_ref
                Real delta_x = -met_h_xi/mag  * d_zref;
                Real delta_y = -met_h_eta/mag * d_zref;
                Real delta_z = 1.0/mag * d_zref;

                // Position of the current node (indx:0,0,1)
                Real x0 = ((Real) i + 0.5) * dx[0];
                Real y0 = ((Real) j + 0.5) * dx[1];

                // Final position at end of vector
                x_pos_arr(i,j,k) = x0 + delta_x;
                y_pos_arr(i,j,k) = y0 + delta_y;
                z_pos_arr(i,j,k) = z_phys_arr(i,j,k) + delta_z;

                // Destination position must be contained on the current process!
                Real pos[] = {x_pos_arr(i,j,k),y_pos_arr(i,j,k),0.5*dx[2]};
                AMREX_ASSERT_WITH_MESSAGE( grb.contains(&pos[0]),
                                           "Query point outside of proc domain!");
            });
        }
    }
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
    auto& fields   = m_fields[lev];
    auto& averages = m_averages[lev];
    const auto & geom     = m_geom[lev];

    auto& z_phys   = m_z_phys_nd[lev];
    auto& x_pos    = m_x_pos[lev];
    auto& y_pos    = m_y_pos[lev];
    auto& z_pos    = m_z_pos[lev];

    auto& i_indx   = m_i_indx[lev];
    auto& j_indx   = m_j_indx[lev];
    auto& k_indx   = m_k_indx[lev];

    auto& ncell_plane   = m_ncell_plane[lev];
    auto& plane_average = m_plane_average[lev];

    // Set factors for time averaging
    Real d_fact_new, d_fact_old;
    if (m_t_avg && m_t_init[lev]) {
        d_fact_new = m_fact_new;
        d_fact_old = m_fact_old;
    } else {
        d_fact_new = 1.0;
        d_fact_old = 0.0;
    }

    // GPU array to accumulate averages into
    Gpu::DeviceVector<Real> pavg(plane_average.size(), 0.0);
    Real* plane_avg = pavg.data();

    // Vectors for normalization and buffer storage
    Vector<Real> denom(plane_average.size(),0.0);
    Vector<Real> val_old(plane_average.size(),0.0);

    // Averages over all the fields
    //----------------------------------------------------------
    Box domain = geom.Domain();

    Array<int,AMREX_SPACEDIM> is_per = {0,0,0};
    for (int idim(0); idim < AMREX_SPACEDIM-1; ++idim) {
        if (geom.isPeriodic(idim)) is_per[idim] = 1;
    }

    for (int imf(0); imf < m_nvar; ++imf) {
        denom[imf]   = 1.0 / (Real)ncell_plane[imf];
        val_old[imf] = plane_average[imf]*d_fact_old;

#ifdef _OPENMP
#pragma omp parallel if (Gpu::notInLaunchRegion())
#endif
        for (MFIter mfi(*fields[imf], TileNoZ()); mfi.isValid(); ++mfi) {
            Box vbx = mfi.validbox(); // This is the grid (not tile)
            Box pbx = mfi.tilebox();  // This is the tile (not grid)
            pbx.setSmall(2,0); pbx.setBig(2,0);

            // Avoid double counting nodal data by changing the high end when we are
            //     at the high side of the grid (not just of the tile)
            IndexType ixt = averages[imf]->boxArray().ixType();
            for (int idim(0); idim < AMREX_SPACEDIM-1; ++idim) {
                if ( ixt.nodeCentered(idim)  && (pbx.bigEnd(idim) == vbx.bigEnd(idim)) ) {
                    int dom_hi = domain.bigEnd(idim)+1;
                    if (pbx.bigEnd(idim) < dom_hi || is_per[idim]) {
                        pbx.growHi(idim,-1);
                    }
                }
            }

            auto mf_arr = fields[imf]->const_array(mfi);

            if (m_interp) {
                const auto plo   = geom.ProbLoArray();
                const auto dxInv = geom.InvCellSizeArray();
                const auto z_phys_arr = z_phys->const_array(mfi);
                auto x_pos_arr = x_pos->array(mfi);
                auto y_pos_arr = y_pos->array(mfi);
                auto z_pos_arr = z_pos->array(mfi);
                ParallelFor(Gpu::KernelInfo().setReduction(true), pbx, [=]
                AMREX_GPU_DEVICE(int i, int j, int k, Gpu::Handler const& handler) noexcept
                {
                    Real interp{0};
                    trilinear_interp_T(x_pos_arr(i,j,k), y_pos_arr(i,j,k), z_pos_arr(i,j,k),
                                       &interp, mf_arr, z_phys_arr, plo, dxInv, 1);
                    Real val = interp;
                    Gpu::deviceReduceSum(&plane_avg[imf], val, handler);
                });
            } else {
                auto k_arr  = k_indx->const_array(mfi);
                auto j_arr  = j_indx ? j_indx->const_array(mfi) : Array4<const int> {};
                auto i_arr  = i_indx ? i_indx->const_array(mfi) : Array4<const int> {};
                ParallelFor(Gpu::KernelInfo().setReduction(true), pbx, [=]
                AMREX_GPU_DEVICE(int i, int j, int k, Gpu::Handler const& handler) noexcept
                {
                    int mk = k_arr(i,j,k);
                    int mj = j_arr ? j_arr(i,j,k) : j;
                    int mi = i_arr ? i_arr(i,j,k) : i;
                    Real val = mf_arr(mi,mj,mk);
                    Gpu::deviceReduceSum(&plane_avg[imf], val, handler);
                });
            }
        }
    }

    //
    //------------------------------------------------------------------------
    // Averages for the tangential velocity magnitude
    // (This is cell-centered so we don't need to worry about double-counting)
    //------------------------------------------------------------------------
    //
    {
        int imf  = 0;
        int iavg = m_navg - 1;
        denom[iavg]   = 1.0 / (Real)ncell_plane[iavg];
        val_old[iavg] = plane_average[iavg]*d_fact_old;

#ifdef _OPENMP
#pragma omp parallel if (Gpu::notInLaunchRegion())
#endif
        for (MFIter mfi(*averages[iavg], TileNoZ()); mfi.isValid(); ++mfi)
        {
            Box pbx = mfi.tilebox();
            pbx.setSmall(2,0); pbx.setBig(2,0);

            // Last element is Umag and always cell centered
            auto u_mf_arr = fields[imf  ]->const_array(mfi);
            auto v_mf_arr = fields[imf+1]->const_array(mfi);

            if (m_interp) {
                const auto plo   = m_geom[lev].ProbLoArray();
                const auto dxInv = m_geom[lev].InvCellSizeArray();
                const auto z_phys_arr = z_phys->const_array(mfi);
                auto x_pos_arr = x_pos->array(mfi);
                auto y_pos_arr = y_pos->array(mfi);
                auto z_pos_arr = z_pos->array(mfi);
                ParallelFor(Gpu::KernelInfo().setReduction(true), pbx, [=]
                AMREX_GPU_DEVICE(int i, int j, int k, Gpu::Handler const& handler) noexcept
                {
                    Real u_interp{0};
                    Real v_interp{0};
                    trilinear_interp_T(x_pos_arr(i,j,k), y_pos_arr(i,j,k), z_pos_arr(i,j,k),
                                            &u_interp, u_mf_arr, z_phys_arr, plo, dxInv, 1);
                    trilinear_interp_T(x_pos_arr(i,j,k), y_pos_arr(i,j,k), z_pos_arr(i,j,k),
                                            &v_interp, v_mf_arr, z_phys_arr, plo, dxInv, 1);
                    const Real val   = std::sqrt(u_interp*u_interp + v_interp*v_interp);
                    Gpu::deviceReduceSum(&plane_avg[iavg], val, handler);
                });
            } else {
                auto k_arr = k_indx->const_array(mfi);
                auto j_arr = j_indx ? j_indx->const_array(mfi) : Array4<const int> {};
                auto i_arr = i_indx ? i_indx->const_array(mfi) : Array4<const int> {};
                ParallelFor(Gpu::KernelInfo().setReduction(true), pbx, [=]
                AMREX_GPU_DEVICE(int i, int j, int k, Gpu::Handler const& handler) noexcept
                {
                    int mk = k_arr(i,j,k);
                    int mj = j_arr ? j_arr(i,j,k) : j;
                    int mi = i_arr ? i_arr(i,j,k) : i;
                    const Real u_val = 0.5 * (u_mf_arr(mi,mj,mk) + u_mf_arr(mi+1,mj  ,mk));
                    const Real v_val = 0.5 * (v_mf_arr(mi,mj,mk) + v_mf_arr(mi  ,mj+1,mk));
                    const Real val   = std::sqrt(u_val*u_val + v_val*v_val);
                    Gpu::deviceReduceSum(&plane_avg[iavg], val, handler);
                });
            }
        }
    }

    // Copy to host and sum across procs
    Gpu::copy(Gpu::deviceToHost, pavg.begin(), pavg.end(), plane_average.begin());
    ParallelDescriptor::ReduceRealSum(plane_average.data(), plane_average.size());

    // No spatial variation with plane averages
    for (int iavg(0); iavg < m_navg; ++iavg){
        plane_average[iavg] *= denom[iavg]*d_fact_new;
        plane_average[iavg] += val_old[iavg];
        averages[iavg]->setVal(plane_average[iavg]);
    }
}


// Fill 2D MF with local averages
void
MOSTAverage::compute_region_averages(int lev)
{
    // Peel back the level
    auto& fields   = m_fields[lev];
    auto& averages = m_averages[lev];
    const auto & geom     = m_geom[lev];

    auto& z_phys   = m_z_phys_nd[lev];
    auto& x_pos    = m_x_pos[lev];
    auto& y_pos    = m_y_pos[lev];
    auto& z_pos    = m_z_pos[lev];

    auto& i_indx   = m_i_indx[lev];
    auto& j_indx   = m_j_indx[lev];
    auto& k_indx   = m_k_indx[lev];

    // Set factors for time averaging
    Real d_fact_new, d_fact_old;
    if (m_t_avg && m_t_init[lev]) {
        d_fact_new = m_fact_new;
        d_fact_old = m_fact_old;
    } else {
        d_fact_new = 1.0;
        d_fact_old = 0.0;
    }

    // Number of cells contained in the local average
    const Real denom = 1.0 / (Real) m_ncell_region;

    // Capture radius for device
    int d_radius = m_radius;

    // Averages over all the fields
    //----------------------------------------------------------
    for (int imf(0); imf < m_nvar; ++imf) {
#ifdef _OPENMP
#pragma omp parallel if (Gpu::notInLaunchRegion())
#endif
        for (MFIter mfi(*fields[imf], TileNoZ()); mfi.isValid(); ++mfi) {
            Box pbx = mfi.tilebox(); pbx.setSmall(2,0); pbx.setBig(2,0);

            auto mf_arr = fields[imf]->const_array(mfi);
            auto ma_arr = averages[imf]->array(mfi);

            if (m_interp) {
                const auto plo   = geom.ProbLoArray();
                const auto dx    = geom.CellSizeArray();
                const auto dxInv = geom.InvCellSizeArray();
                const auto z_phys_arr = z_phys->const_array(mfi);
                auto x_pos_arr = x_pos->array(mfi);
                auto y_pos_arr = y_pos->array(mfi);
                auto z_pos_arr = z_pos->array(mfi);
                ParallelFor(pbx, [=] AMREX_GPU_DEVICE(int i, int j, int k) noexcept
                {
                    ma_arr(i,j,k) *= d_fact_old;

                    Real met_h_zeta = Compute_h_zeta_AtCellCenter(i,j,k,dxInv,z_phys_arr);
                    for (int lk(-d_radius); lk <= (d_radius); ++lk) {
                      for (int lj(-d_radius); lj <= (d_radius); ++lj) {
                        for (int li(-d_radius); li <= (d_radius); ++li) {
                            Real interp{0};
                            Real xp = x_pos_arr(i,j,k) + li*dx[0];
                            Real yp = y_pos_arr(i,j,k) + lj*dx[1];
                            Real zp = z_pos_arr(i,j,k) + met_h_zeta*lk*dx[2];
                            trilinear_interp_T(xp, yp, zp, &interp, mf_arr, z_phys_arr, plo, dxInv, 1);
                            Real val = denom * interp * d_fact_new;
                            ma_arr(i,j,k) += val;
                        }
                      }
                    }
                });
            } else {
                auto k_arr = k_indx->const_array(mfi);
                auto j_arr = j_indx ? j_indx->const_array(mfi) : Array4<const int> {};
                auto i_arr = i_indx ? i_indx->const_array(mfi) : Array4<const int> {};
                ParallelFor(pbx, [=] AMREX_GPU_DEVICE(int i, int j, int k) noexcept
                {
                    ma_arr(i,j,k) *= d_fact_old;

                    int mk = k_arr(i,j,k);
                    int mj = j_arr ? j_arr(i,j,k) : j;
                    int mi = i_arr ? i_arr(i,j,k) : i;
                    for (int lk(mk-d_radius); lk <= (mk+d_radius); ++lk) {
                      for (int lj(mj-d_radius); lj <= (mj+d_radius); ++lj) {
                        for (int li(mi-d_radius); li <= (mi+d_radius); ++li) {
                            Real val = denom * mf_arr(li, lj, lk) * d_fact_new;
                            ma_arr(i,j,k) += val;
                        }
                      }
                    }
                });
            }
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
#pragma omp parallel if (Gpu::notInLaunchRegion())
#endif
        for (MFIter mfi(*averages[iavg], TileNoZ()); mfi.isValid(); ++mfi) {
            Box pbx = mfi.tilebox(); pbx.setSmall(2,0); pbx.setBig(2,0);

            auto u_mf_arr = fields[imf]->const_array(mfi);
            auto v_mf_arr = fields[imf+1]->const_array(mfi);
            auto ma_arr   = averages[iavg]->array(mfi);

            if (m_interp) {
                const auto plo   = geom.ProbLoArray();
                const auto dx    = geom.CellSizeArray();
                const auto dxInv = geom.InvCellSizeArray();
                const auto z_phys_arr = z_phys->const_array(mfi);
                auto x_pos_arr = x_pos->array(mfi);
                auto y_pos_arr = y_pos->array(mfi);
                auto z_pos_arr = z_pos->array(mfi);
                ParallelFor(pbx, [=] AMREX_GPU_DEVICE(int i, int j, int k) noexcept
                {
                    ma_arr(i,j,k) *= d_fact_old;

                    Real met_h_zeta = Compute_h_zeta_AtCellCenter(i,j,k,dxInv,z_phys_arr);
                    for (int lk(-d_radius); lk <= (d_radius); ++lk) {
                      for (int lj(-d_radius); lj <= (d_radius); ++lj) {
                        for (int li(-d_radius); li <= (d_radius); ++li) {
                            Real u_interp{0};
                            Real v_interp{0};
                            Real xp = x_pos_arr(i,j,k) + li*dx[0];
                            Real yp = y_pos_arr(i,j,k) + lj*dx[1];
                            Real zp = z_pos_arr(i,j,k) + met_h_zeta*lk*dx[2];
                            trilinear_interp_T(xp, yp, zp, &u_interp, u_mf_arr, z_phys_arr, plo, dxInv, 1);
                            trilinear_interp_T(xp, yp, zp, &v_interp, v_mf_arr, z_phys_arr, plo, dxInv, 1);
                            Real mag = std::sqrt(u_interp*u_interp + v_interp*v_interp);
                            Real val = denom * mag * d_fact_new;
                            ma_arr(i,j,k) += val;
                        }
                      }
                    }
                });
            } else {
                auto k_arr = k_indx->const_array(mfi);
                auto j_arr = j_indx ? j_indx->const_array(mfi) : Array4<const int> {};
                auto i_arr = i_indx ? i_indx->const_array(mfi) : Array4<const int> {};
                ParallelFor(pbx, [=] AMREX_GPU_DEVICE(int i, int j, int k) noexcept
                {
                    ma_arr(i,j,k) *= d_fact_old;

                    int mk = k_arr(i,j,k);
                    int mj = j_arr ? j_arr(i,j,k) : j;
                    int mi = i_arr ? i_arr(i,j,k) : i;
                    for (int lk(mk-d_radius); lk <= (mk+d_radius); ++lk) {
                      for (int lj(mj-d_radius); lj <= (mj+d_radius); ++lj) {
                        for (int li(mi-d_radius); li <= (mi+d_radius); ++li) {
                            const Real u_val = 0.5 * (u_mf_arr(li,lj,lk) + u_mf_arr(li+1,lj  ,lk));
                            const Real v_val = 0.5 * (v_mf_arr(li,lj,lk) + v_mf_arr(li  ,lj+1,lk));
                            const Real mag   = std::sqrt(u_val*u_val + v_val*v_val);
                            Real val = denom * mag * d_fact_new;
                            ma_arr(i,j,k) += val;
                        }
                      }
                    }
                });
            }

            // Fill interior ghost cells and any ghost cells outside a periodic domain
            //***********************************************************************************
            averages[iavg]->FillBoundary(geom.periodicity());
        }
    }


    // Need to fill ghost cells outside the domain if not periodic
    bool not_per_x = !(geom.periodicity().isPeriodic(0));
    bool not_per_y = !(geom.periodicity().isPeriodic(1));
    if (not_per_x || not_per_y) {
        Box domain = geom.Domain();
        for (int iavg(0); iavg < m_navg; ++iavg) {
            IndexType ixt = averages[iavg]->boxArray().ixType();
            Box ldomain   = domain; ldomain.convert(ixt);
            IntVect ng    = averages[iavg]->nGrowVect(); ng[2]=0;
#ifdef _OPENMP
#pragma omp parallel if (Gpu::notInLaunchRegion())
#endif
            for (MFIter mfi(*averages[iavg], TileNoZ()); mfi.isValid(); ++mfi) {
                Box gpbx = mfi.growntilebox(ng); gpbx.setSmall(2,0); gpbx.setBig(2,0);

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

    for (MFIter mfi(*averages[navg], TileNoZ()); mfi.isValid(); ++mfi) {
        Box bx  = mfi.tilebox(); bx.setBig(2,0);
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


// Write ijk indices
void
MOSTAverage::write_norm_indices(int lev)
{
    // Peel back the level
    auto& averages = m_averages[lev];
    auto& k_indx   = m_k_indx[lev];
    auto& j_indx   = m_j_indx[lev];
    auto& i_indx   = m_i_indx[lev];

    int navg = m_navg - 1;

    std::ofstream ofile;
    ofile.open ("MOST_ijk_indices.txt");
    ofile << "IJK indices used to compute averages via MOSTAverages class:\n";

    for (MFIter mfi(*averages[navg], TileNoZ()); mfi.isValid(); ++mfi) {
        Box bx  = mfi.tilebox(); bx.setBig(2,0);
        int il = bx.smallEnd(0); int iu = bx.bigEnd(0);
        int jl = bx.smallEnd(1); int ju = bx.bigEnd(1);

        auto k_arr = k_indx->array(mfi);
        auto j_arr = j_indx ? j_indx->array(mfi) : Array4<int> {};
        auto i_arr = i_indx ? i_indx->array(mfi) : Array4<int> {};

        for (int j(jl); j <= ju; ++j) {
            for (int i(il); i <= iu; ++i) {
                ofile << "(I1,J1,K1): " << "(" << i << "," << j << "," << 0 << ")" << "\n";

                int k = 0;
                int km = k_arr(i,j,k);
                int jm = j_arr ? j_arr(i,j,k) : j;
                int im = i_arr ? i_arr(i,j,k) : i;

                ofile << "(I2,J2,K2): "
                      << "(" << im << "," << jm << "," << km << ")" << "\n";
                ofile << "\n";
            }
        }
    }
    ofile.close();
}


// Write position of XZ plane
void
MOSTAverage::write_xz_positions(int lev, int j)
{
    // Peel back the level
    auto& x_pos_mf  = m_x_pos[lev];
    auto& z_pos_mf  = m_z_pos[lev];

    std::ofstream ofile;
    ofile.open ("MOST_xz_positions.txt");

    for (MFIter mfi(*x_pos_mf, TileNoZ()); mfi.isValid(); ++mfi) {
        Box bx  = mfi.tilebox(); bx.setBig(2,0);
        int il = bx.smallEnd(0); int iu = bx.bigEnd(0);

        auto x_pos_arr  = x_pos_mf->array(mfi);
        auto z_pos_arr  = z_pos_mf->array(mfi);

        int k  = 0;
        for (int i(il); i <= iu; ++i)
            ofile << x_pos_arr(i,j,k) << ' ' << z_pos_arr(i,j,k) << "\n";
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

    for (MFIter mfi(*averages[navg], TileNoZ()); mfi.isValid(); ++mfi) {
        Box bx  = mfi.tilebox(); bx.setBig(2,0);
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
