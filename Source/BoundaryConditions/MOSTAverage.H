#ifndef MOSTAverage_H
#define MOSTAverage_H

#include <AMReX_Gpu.H>
#include <AMReX_FArrayBox.H>
#include <AMReX_MultiFab.H>
#include <AMReX_iMultiFab.H>
#include <AMReX_ParmParse.H>
#include <IndexDefines.H>
#include <TerrainMetrics.H>

class MOSTAverage {
public:
    explicit MOSTAverage (amrex::Vector<amrex::Geometry>  geom,
                          amrex::Vector<amrex::Vector<amrex::MultiFab>>& vars_old,
                          amrex::Vector<std::unique_ptr<amrex::MultiFab>>& Theta_prim,
                          amrex::Vector<std::unique_ptr<amrex::MultiFab>>& Qv_prim,
                          amrex::Vector<std::unique_ptr<amrex::MultiFab>>& Qr_prim,
                          amrex::Vector<std::unique_ptr<amrex::MultiFab>>& z_phys_nd);

    // MOSTAverage() = default;
    ~MOSTAverage ()
    {}

    // Declare a default move constructor so we ensure the destructor is
    // not called when we return an object of this class by value
    MOSTAverage (MOSTAverage&&)  noexcept = default;

    // Delete the move assignment operator
    MOSTAverage& operator=(MOSTAverage&& other)  noexcept = delete;

    // Delete the copy constructor
    MOSTAverage (const MOSTAverage& other) = delete;
    //
    // Delete the copy assignment operator
    MOSTAverage& operator=(const MOSTAverage& other) = delete;

    // Reset the pointers to field MFs
    void update_field_ptrs (int lev,
                            amrex::Vector<amrex::Vector<amrex::MultiFab>>& vars_old,
                            amrex::Vector<std::unique_ptr<amrex::MultiFab>>& Theta_prim,
                            amrex::Vector<std::unique_ptr<amrex::MultiFab>>& Qv_prim,
                            amrex::Vector<std::unique_ptr<amrex::MultiFab>>& Qr_prim);

    // Compute ncells per plane
    void set_plane_normalization ();

    // Compute ncells per plane
    void set_region_normalization ()
    {m_ncell_region = (2 * m_radius + 1) * (2 * m_radius + 1) * (2 * m_radius + 1);}

    // Populate a 2D iMF k_indx (w/o terrain)
    void set_k_indices_N ();

    // Populate a 2D iMF k_indx (w/ terrain)
    void set_k_indices_T ();

    // Populate all 2D iMFs ijk_indx (w/ terrain)
    void set_norm_indices_T ();

    // Populate positions (w/ terrain & norm vector & interpolation)
    void set_z_positions_T ();

    // Populate positions (w/ terrain & norm vector & interpolation)
    void set_norm_positions_T ();

    // Driver for the different average policies
    void compute_averages (int lev);

    // Fill averages for policy::plane
    void compute_plane_averages (int lev);

    // Fill averages for policy::point
    void compute_region_averages (int lev);

    // Write k indices
    void write_k_indices (int lev);

    // Write ijk indices
    void write_norm_indices (int lev);

    // Write XZ planar positions
    void write_xz_positions (int lev, int j);

    // Write averages on 2D mf
    void write_averages (int lev);

    // Get pointer to the 2D mf of averages
    [[nodiscard]] const amrex::MultiFab* get_average (int lev, int comp) const { return m_averages[lev][comp].get(); }

    // Get z_ref (may be computed from specified k_indx)
    [[nodiscard]] amrex::Real get_zref () const { return m_zref; }

    /**
     * Function to compute trilinear interpolation with terrain.
     *
     * @param[in] xp X-position
     * @param[in] yp Y-position
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
        amrex::Real zval= 0.0;
        int kmax  = ubound(z_arr).z;
        int i_new = (int) (xp * dxi[0] - 0.5);
        int j_new = (int) (yp * dxi[1] - 0.5);
        amrex::Real z_target = zp;
        for (int lk(0); lk<kmax; ++lk) {
            amrex::Real z_lo = 0.25 * ( z_arr(i_new,j_new  ,lk  ) + z_arr(i_new+1,j_new  ,lk  )
                                      + z_arr(i_new,j_new+1,lk  ) + z_arr(i_new+1,j_new+1,lk  ) );
            amrex::Real z_hi = 0.25 * ( z_arr(i_new,j_new  ,lk+1) + z_arr(i_new+1,j_new  ,lk+1)
                                      + z_arr(i_new,j_new+1,lk+1) + z_arr(i_new+1,j_new+1,lk+1) );
            if (z_target > z_lo && z_target < z_hi){
                found = true;
                zval = (amrex::Real) lk + ((z_target - z_lo) / (z_hi - z_lo)) + 0.5;
                break;
            }
        }

        amrex::ignore_unused(found);
        AMREX_ASSERT_WITH_MESSAGE(found, "MOSTAverage: Height above terrain not found, try increasing z_ref!");

        const amrex::RealVect lx((xp - plo[0])*dxi[0] + 0.5,
                                 (yp - plo[1])*dxi[1] + 0.5,
                                  zval);

        const amrex::IntVect ijk = lx.floor();

        int i = ijk[0]; int j = ijk[1]; int k = ijk[2];

        // Weights
        const amrex::RealVect sx_hi = lx  - ijk;
        const amrex::RealVect sx_lo = 1.0 - sx_hi;

        for (int n = 0; n < interp_comp; n++)
            interp_vals[n] = sx_lo[0]*sx_lo[1]*sx_lo[2]*interp_array(i-1, j-1, k-1,n) +
                             sx_lo[0]*sx_lo[1]*sx_hi[2]*interp_array(i-1, j-1, k  ,n) +
                             sx_lo[0]*sx_hi[1]*sx_lo[2]*interp_array(i-1, j  , k-1,n) +
                             sx_lo[0]*sx_hi[1]*sx_hi[2]*interp_array(i-1, j  , k  ,n) +
                             sx_hi[0]*sx_lo[1]*sx_lo[2]*interp_array(i  , j-1, k-1,n) +
                             sx_hi[0]*sx_lo[1]*sx_hi[2]*interp_array(i  , j-1, k  ,n) +
                             sx_hi[0]*sx_hi[1]*sx_lo[2]*interp_array(i  , j  , k-1,n) +
                             sx_hi[0]*sx_hi[1]*sx_hi[2]*interp_array(i  , j  , k  ,n);
    }

protected:

    // Passed through constructor
    //--------------------------------------------
    const amrex::Vector<amrex::Geometry> m_geom;                     // Geometry at each level
    amrex::Vector<amrex::Vector<amrex::MultiFab*>> m_fields;         // Ptr to fields to be averaged
    amrex::Vector<amrex::MultiFab*> m_z_phys_nd;                     // Ptr to terrain height coords

    // General vars for multiple or all policies
    //--------------------------------------------
    int m_nvar{5};                                                   // 5 fields for U/V/T/Qv/Qr
    int m_navg{6};                                                   // 6 averages for U/V/T/Qv/Tv/Umag
    int m_maxlev{0};                                                 // Total number of levels
    int m_policy{0};                                                 // Policy for type of averaging
    amrex::Real m_zref{10.0};                                        // Height above surface for MOST BC
    std::string m_pp_prefix {"erf"};                                 // ParmParse prefix
    amrex::Vector<std::unique_ptr<amrex::MultiFab>> m_x_pos;         // Ptr to 2D mf to hold x position (maxlev)
    amrex::Vector<std::unique_ptr<amrex::MultiFab>> m_y_pos;         // Ptr to 2D mf to hold y position (maxlev)
    amrex::Vector<std::unique_ptr<amrex::MultiFab>> m_z_pos;         // Ptr to 2D mf to hold z position (maxlev)
    amrex::Vector<std::unique_ptr<amrex::iMultiFab>> m_i_indx;       // Ptr to 2D imf to hold i indices (maxlev)
    amrex::Vector<std::unique_ptr<amrex::iMultiFab>> m_j_indx;       // Ptr to 2D imf to hold j indices (maxlev)
    amrex::Vector<std::unique_ptr<amrex::iMultiFab>> m_k_indx;       // Ptr to 2D imf to hold k indices (maxlev)
    amrex::Vector<amrex::Vector<std::unique_ptr<amrex::MultiFab>>> m_averages; // Ptr to 2D mf to hold averages (maxlev,navg)

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

    // Time average w/ exponential filter fun
    //--------------------------------------------
    bool m_t_avg{false};                                             // Flag to do moving average in time
    amrex::Vector<int> m_t_init;                                     // Flag to specify if averages are initialized
    amrex::Real m_time_window{1.0e-16};                              // Width of the exp filter function
    amrex::Real m_fact_new, m_fact_old;                              // Time average factors for new and old means
};
#endif
