#ifndef ERF_MOSTAverage_H
#define ERF_MOSTAverage_H

#include "AMReX_Gpu.H"
#include "AMReX_FArrayBox.H"
#include "AMReX_MultiFab.H"
#include "AMReX_iMultiFab.H"
#include "AMReX_ParmParse.H"

#include "ERF_IndexDefines.H"
#include "ERF_TerrainMetrics.H"
#include "ERF_DataStruct.H"
#include "ERF_EB.H"

class MOSTAverage {
public:
    /**
     * Construct the MOST averaging helper.
     *
     * @param[in] geom           geometry for all AMR levels
     * @param[in] has_zphys      whether nodal physical-height data are available
     * @param[in] a_pp_prefix    ParmParse prefix
     * @param[in] m_mesh_type    mesh type
     * @param[in] m_terrain_type terrain representation
     * @param[in] eb_vec         optional embedded-boundary geometry data
     */
    explicit MOSTAverage (amrex::Vector<amrex::Geometry>  geom,
                          const bool& has_zphys,
                          std::string a_pp_prefix,
                          const MeshType& m_mesh_type,
                          const TerrainType& m_terrain_type,
                          const amrex::Vector<const eb_*>& eb_vec = {});

    // MOSTAverage() = default;

    /**
     * Destroy the MOST averaging helper.
     */
    ~MOSTAverage ()
    {}

    /**
     * Default move constructor.
     */
    MOSTAverage (MOSTAverage&&)  noexcept = default;

    /**
     * Deleted move-assignment operator.
     *
     * @param[in] other source object
     */
    MOSTAverage& operator=(MOSTAverage&& other)  noexcept = delete;

    /**
     * Deleted copy constructor.
     *
     * @param[in] other source object
     */
    MOSTAverage (const MOSTAverage& other) = delete;

    /**
     * Deleted copy-assignment operator.
     *
     * @param[in] other source object
     */
    MOSTAverage& operator=(const MOSTAverage& other) = delete;

    /**
     * Make MOST-average data structures at one level.
     *
     * @param[in] lev        level index
     * @param[in] vars_old   old-time state and velocity fields
     * @param[in] Theta_prim primitive potential-temperature field
     * @param[in] Qv_prim    primitive water-vapor field
     * @param[in] Qr_prim    primitive rain-water field
     * @param[in] z_phys_nd  nodal physical-height field
     */
    void make_MOSTAverage_at_level (const int& lev,
                                    const amrex::Vector<amrex::MultiFab*>& vars_old,
                                    std::unique_ptr<amrex::MultiFab>& Theta_prim,
                                    std::unique_ptr<amrex::MultiFab>& Qv_prim,
                                    std::unique_ptr<amrex::MultiFab>& Qr_prim,
                                    std::unique_ptr<amrex::MultiFab>& z_phys_nd);

    /**
     * Reset pointers to field MultiFabs.
     *
     * @param[in] lev        level index
     * @param[in] vars_old   old-time state and velocity fields
     * @param[in] Theta_prim primitive potential-temperature fields by level
     * @param[in] Qv_prim    primitive water-vapor fields by level
     * @param[in] Qr_prim    primitive rain-water fields by level
     */
    void update_field_ptrs (const int& lev,
                            amrex::Vector<amrex::Vector<amrex::MultiFab>>& vars_old,
                            amrex::Vector<std::unique_ptr<amrex::MultiFab>>& Theta_prim,
                            amrex::Vector<std::unique_ptr<amrex::MultiFab>>& Qv_prim,
                            amrex::Vector<std::unique_ptr<amrex::MultiFab>>& Qr_prim);

    /**
     * Update the rotated fields.
     *
     * @param[in] lev level index
     */
    void set_rotated_fields (const int& lev);

    /**
     * Compute number of cells per averaging plane.
     *
     * @param[in] lev level index
     */
    void set_plane_normalization (const int& lev);

    /**
     * Compute total embedded-boundary surface area.
     *
     * @param[in] lev level index
     */
    void set_eb_normalization (const int& lev);

    /**
     * Compute number of cells in the region average.
     *
     * @param[in] lev level index
     */
    void set_region_normalization (const int& /*lev*/)
    {m_ncell_region = (2 * m_radius + 1) * (2 * m_radius + 1) * (2 * m_radius + 1);}

    /**
     * Populate the 2D k-index iMultiFab without terrain.
     *
     * @param[in] lev level index
     */
    void set_k_indices_N (const int& lev);

    /**
     * Populate the 2D k-index iMultiFab with terrain.
     *
     * @param[in] lev level index
     */
    void set_k_indices_T (const int& lev);

    /**
     * Populate all 2D normal-index iMultiFabs with terrain.
     *
     * @param[in] lev level index
     */
    void set_norm_indices_T (const int& lev);

    /**
     * Populate terrain-aware z positions.
     *
     * @param[in] lev level index
     */
    void set_z_positions_T (const int& lev);

    /**
     * Populate embedded-boundary z positions.
     *
     * @param[in] lev level index
     */
    void set_z_positions_EB (const int& lev);

    /**
     * Populate terrain-aware normal positions.
     *
     * @param[in] lev level index
     */
    void set_norm_positions_T (const int& lev);

    /**
     * Driver for the selected average policy.
     *
     * @param[in] lev level index
     */
    void compute_averages (const int& lev);

    /**
     * Fill averages for the plane policy.
     *
     * @param[in] lev level index
     */
    void compute_plane_averages (const int& lev);

    /**
     * Fill averages for the point or region policy.
     *
     * @param[in] lev level index
     */
    void compute_region_averages (const int& lev);

    /**
     * Fill averages for the embedded-boundary policy.
     *
     * @param[in] lev level index
     */
    void compute_eb_averages (const int& lev);

    /**
     * Write k-index data.
     *
     * @param[in] lev level index
     */
    void write_k_indices (const int& lev);

    /**
     * Write normal-index data.
     *
     * @param[in] lev level index
     */
    void write_norm_indices (const int& lev);

    /**
     * Write XZ planar positions.
     *
     * @param[in] lev level index
     * @param[in] j   y-index of the plane to write
     */
    void write_xz_positions (const int& lev,
                             const int& j);

    /**
     * Write averages on the 2D MultiFab.
     *
     * @param[in] lev level index
     */
    void write_averages (const int& lev);

    /**
     * Return one 2D average MultiFab.
     *
     * @param[in] lev  level index
     * @param[in] comp average component index
     */
    [[nodiscard]] const amrex::MultiFab* get_average (const int& lev, const int& comp) const { return m_averages[lev][comp].get(); }

    /**
     * Return z_ref, which may be computed from a specified k-index.
     *
     * @param[in] lev level index
     */
    [[nodiscard]] amrex::MultiFab* get_zref (const int& lev) const { return m_zref[lev].get(); }

    /**
     * Return the k-index iMultiFab.
     *
     * @param[in] lev level index
     */
    [[nodiscard]] const amrex::iMultiFab* get_k_indices (const int& lev) const { return m_k_indx[lev].get(); }

    /**
     * Function to compute trilinear interpolation with terrain.
     *
     * @param[in] xp X-position
     * @param[in] yp Y-position
     * @param[in] zp Z-position
     * @param[out] interp_vals Values interpolated
     * @param[in] interp_array Array to interpolate on
     * @param[in] z_arr Physical heights
     * @param[in] plo Problem lower bounds
     * @param[in] dxi Inverse cell size array
     * @param[in] interp_comp Number of components to interpolate
     */
    AMREX_GPU_HOST_DEVICE AMREX_INLINE
    static void trilinear_interp_T (const amrex::Real& xp,
                                    const amrex::Real& yp,
                                    const amrex::Real& zp,
                                    amrex::Real* interp_vals,
                                    amrex::Array4<amrex::Real const> const& interp_array,
                                    amrex::Array4<amrex::Real const> const& z_arr,
                                    const amrex::GpuArray<amrex::Real, AMREX_SPACEDIM>& plo,
                                    const amrex::GpuArray<amrex::Real, AMREX_SPACEDIM>& dxi,
                                    const int interp_comp)
    {
        // Search to get z/k
        bool found = false;
        int kmax   = ubound(z_arr).z;
        amrex::Real zval     = zero;
        amrex::Real z_target = zp;

        // Map position to i,j (must be same mapping in cpp file)
        amrex::Real ireal = (xp - plo[0]) * dxi[0];
        amrex::Real jreal = (yp - plo[1]) * dxi[1];
        int i_new = (int) (ireal - myhalf);
        int j_new = (int) (jreal - myhalf);

        for (int lk(0); lk<kmax; ++lk) {
            amrex::Real z_lo = fourth * ( z_arr(i_new,j_new  ,lk  ) + z_arr(i_new+1,j_new  ,lk  )
                                      + z_arr(i_new,j_new+1,lk  ) + z_arr(i_new+1,j_new+1,lk  ) );
            amrex::Real z_hi = fourth * ( z_arr(i_new,j_new  ,lk+1) + z_arr(i_new+1,j_new  ,lk+1)
                                      + z_arr(i_new,j_new+1,lk+1) + z_arr(i_new+1,j_new+1,lk+1) );
            if (z_target > z_lo && z_target < z_hi){
                found = true;
                zval  = (amrex::Real) lk + ((z_target - z_lo) / (z_hi - z_lo)) + myhalf;
                break;
            }
        }

        amrex::ignore_unused(found);
        AMREX_ALWAYS_ASSERT_WITH_MESSAGE(found, "MOSTAverage: Height above terrain not found, try increasing z_ref!");

        // NOTE: This is the point ahead of the current i,j (e.g. i/j_new + 1)
        const amrex::RealVect lx(ireal + myhalf, jreal + myhalf, zval);

        const amrex::IntVect ijk = lx.floor();

        int i = ijk[0]; int j = ijk[1]; int k = ijk[2];

        // Convert ijk (IntVect) to a RealVect explicitly
        amrex::RealVect ijk_r(static_cast<amrex::Real>(ijk[0]),
                              static_cast<amrex::Real>(ijk[1]),
                              static_cast<amrex::Real>(ijk[2]));

        // Weights
        const amrex::RealVect sx_hi = lx  - ijk_r;
        const amrex::RealVect sx_lo = one - sx_hi;

        for (int n = 0; n < interp_comp; n++) {
            interp_vals[n] = sx_lo[0]*sx_lo[1]*sx_lo[2]*interp_array(i-1, j-1, k-1,n) +
                             sx_lo[0]*sx_lo[1]*sx_hi[2]*interp_array(i-1, j-1, k  ,n) +
                             sx_lo[0]*sx_hi[1]*sx_lo[2]*interp_array(i-1, j  , k-1,n) +
                             sx_lo[0]*sx_hi[1]*sx_hi[2]*interp_array(i-1, j  , k  ,n) +
                             sx_hi[0]*sx_lo[1]*sx_lo[2]*interp_array(i  , j-1, k-1,n) +
                             sx_hi[0]*sx_lo[1]*sx_hi[2]*interp_array(i  , j-1, k  ,n) +
                             sx_hi[0]*sx_hi[1]*sx_lo[2]*interp_array(i  , j  , k-1,n) +
                             sx_hi[0]*sx_hi[1]*sx_hi[2]*interp_array(i  , j  , k  ,n);
        }
    }

protected:

    // Passed through constructor
    //--------------------------------------------
    const amrex::Vector<amrex::Geometry> m_geom;                     // Geometry at each level
    amrex::Vector<amrex::Vector<amrex::MultiFab*>> m_fields;         // Ptr to fields to be averaged
    amrex::Vector<amrex::MultiFab*> m_z_phys_nd;                     // Ptr to terrain height coords
    std::string m_pp_prefix;                                         // ParmParse prefix
    MeshType m_mesh_type;                                            // Mesh type
    TerrainType m_terrain_type;                                      // Terrain type

    // General vars for multiple or all policies
    //--------------------------------------------
    int  m_nvar{6};                                                  // 6 fields for U/V/T/Qv/Qr/W
    int  m_navg{6};                                                  // 6 averages for U/V/T/Qv/Tv/Umag
    int  m_maxlev{0};                                                // Total number of levels
    int  m_policy{0};                                                // Policy for type of averaging
    bool m_rotate{false};                                            // Do vector rotations for terrain?
    amrex::Vector<std::unique_ptr<amrex::MultiFab>> m_zref;          // Height above surface for MOST BC
    amrex::Vector<std::unique_ptr<amrex::MultiFab>> m_x_pos;         // Ptr to 2D mf to hold x position (maxlev)
    amrex::Vector<std::unique_ptr<amrex::MultiFab>> m_y_pos;         // Ptr to 2D mf to hold y position (maxlev)
    amrex::Vector<std::unique_ptr<amrex::MultiFab>> m_z_pos;         // Ptr to 2D mf to hold z position (maxlev)
    amrex::Vector<std::unique_ptr<amrex::iMultiFab>> m_i_indx;       // Ptr to 2D imf to hold i indices (maxlev)
    amrex::Vector<std::unique_ptr<amrex::iMultiFab>> m_j_indx;       // Ptr to 2D imf to hold j indices (maxlev)
    amrex::Vector<std::unique_ptr<amrex::iMultiFab>> m_k_indx;       // Ptr to 2D imf to hold k indices (maxlev)
    amrex::Vector<amrex::Vector<std::unique_ptr<amrex::MultiFab>>> m_averages; // Ptr to 2D mf to hold averages (maxlev,navg)
    amrex::Vector<amrex::Vector<std::unique_ptr<amrex::MultiFab>>> m_rot_fields; // Rotated field data

    // Vars for planar average policy
    //--------------------------------------------
    amrex::Vector<amrex::Vector<int>> m_ncell_plane;                 // Number of cells in plane (maxlev,navg)
    amrex::Vector<amrex::Vector<amrex::Real>> m_plane_average;       // Plane avgs (maxlev,navg)

    // Vars for point/region average policy
    //--------------------------------------------
    int m_radius{0};                                                 // Radius around k_index
    int m_ncell_region{1};                                           // Number of cells in local region
    amrex::Vector<int> m_k_in;                                       // Specified k_index for region avg (maxlev)

    // Vars for normal vector policy
    //--------------------------------------------
    bool m_interp{false};                                            // Do interpolation on destination?
    bool m_norm_vec{false};                                          // Use normal vector to find IJK?

    // Vars for EB averaging policy
    //--------------------------------------------
    amrex::Vector<const eb_*> m_eb_vec;                              // Vector of pointers to EB objects (one per level)
    amrex::Vector<amrex::Vector<amrex::Real>> m_total_bndry_area;   // Total EB surface areas (maxlev, navg)

    // Time average w/ exponential filter fun
    //--------------------------------------------
    bool m_t_avg{false};                                             // Flag to do moving average in time
    amrex::Vector<int> m_t_init;                                     // Flag to specify if averages are initialized
    double m_time_window{1.0e-16};                                                // Width of the exp filter function
    amrex::Real m_fact_new, m_fact_old;                              // Time average factors for new and old means

    // Surface velocity correction
    //--------------------------------------------
    bool include_subgrid_vel = false;
    amrex::Vector<amrex::Real> m_Vsg;                                // Subgrid velocity scale (Mahrt & Sun 1995 MWR)

    // Default values
    //--------------------------------------------
    const amrex::Real zref_default = amrex::Real(10.0);
};
#endif
