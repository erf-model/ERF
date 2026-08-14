/**
 * \file ERF_EB.H
 * \brief Declares the embedded-boundary factory manager used by ERF levels.
 */
#ifndef ERF_EB_H_
#define ERF_EB_H_

#include <AMReX_Geometry.H>
#include <AMReX_DistributionMapping.H>
#include <AMReX_BoxArray.H>

#include <AMReX_EB2.H>
#include <AMReX_EBFabFactory.H>

#include <ERF_EBAux.H>

/**
 * \brief Owns and exposes cell-centered and face-centered EB factories.
 *
 * The class builds AMReX EB factories for cell-centered state data and, depending
 * on USE_FC_FACTORY, either native AMReX or ERF auxiliary factories for staggered
 * velocity components.
 */
class eb_ {

  public:

    //! Destroy the EB factory manager.
    ~eb_ ();

    /**
     * \brief Construct an EB description from terrain data.
     * \param a_geom Geometry used to place the terrain EB.
     * \param terrain_fab Terrain height data.
     * \param a_dz_stretched Device vector containing stretched vertical spacing.
     * \param is_anelastic Whether the owning level is using anelastic dynamics.
     */
    eb_ (amrex::Geometry const& a_geom,
         amrex::FArrayBox const& terrain_fab,
         amrex::Gpu::DeviceVector<amrex::Real>& a_dz_stretched,
         bool is_anelastic);
    //! Construct an empty EB factory manager.
    eb_ ();

    /**
     * \brief Define the EB level metadata used by this factory manager.
     * \param level AMR level index.
     * \param a_geom Geometry for the level.
     * \param a_eb_level AMReX EB level owned by the global EB index space.
     * \param is_anelastic Whether the owning level is using anelastic dynamics.
     */
    void define (int level,
                 amrex::Geometry const& a_geom,
                 amrex::EB2::Level const* a_eb_level,
                 bool is_anelastic);

    /**
     * \brief Build cell-centered and staggered EB factories for the level.
     * \param level AMR level index.
     * \param a_geom Geometry for the level.
     * \param ba BoxArray over which factories are defined.
     * \param dm DistributionMapping for the BoxArray.
     * \param a_eb_level AMReX EB level used to populate factory data.
     */
    void make_all_factories ( int level,
                        amrex::Geometry   const& a_geom,
                        amrex::BoxArray            const& ba,
                        amrex::DistributionMapping const& dm,
                        amrex::EB2::Level const& a_eb_level);

    /**
     * \brief Build only the cell-centered EB factory for the level.
     * \param level AMR level index.
     * \param a_geom Geometry for the level.
     * \param ba BoxArray over which the factory is defined.
     * \param dm DistributionMapping for the BoxArray.
     * \param a_eb_level AMReX EB level used to populate factory data.
     */
    void make_cc_factory ( int level,
                        amrex::Geometry   const& a_geom,
                        amrex::BoxArray            const& ba,
                        amrex::DistributionMapping const& dm,
                        amrex::EB2::Level const& a_eb_level);

    //! Number of ghost cells needed for basic EB flags.
    int nghost_basic () const  { return 5; } // nghost_eb_basic ()
    //! Number of ghost cells needed for EB volume data.
    int nghost_volume () const { return 5; } // nghost_eb_volume ()
    //! Number of ghost cells needed for full EB geometric data.
    int nghost_full () const   { return 4; } // nghost_eb_full ()

    // Toggle between eb_aux_ and native AMReX FC factories
    #ifndef USE_FC_FACTORY
    /**
     * \def USE_FC_FACTORY
     * \brief Select native AMReX face-centered factories instead of ERF auxiliary factories.
     */
    #define USE_FC_FACTORY 0 // Default to 0: use eb_aux_ (original)
    #endif

    //! Return the cell-centered EB factory.
    const std::unique_ptr<amrex::EBFArrayBoxFactory>& get_const_factory () const noexcept { return m_factory; }

    /**
     * \brief Disconnect covered neighbors in an EB factory's cell flags.
     * \param factory Factory whose non-const flag data are updated.
     */
    void set_connection_flags (amrex::EBFArrayBoxFactory* factory);

    #if USE_FC_FACTORY
    //! Return the native AMReX x-face EB factory.
    amrex::EBFArrayBoxFactory const* get_u_const_factory() const noexcept { return m_u_factory_fc.get(); }
    //! Return the native AMReX y-face EB factory.
    amrex::EBFArrayBoxFactory const* get_v_const_factory() const noexcept { return m_v_factory_fc.get(); }
    //! Return the native AMReX z-face EB factory.
    amrex::EBFArrayBoxFactory const* get_w_const_factory() const noexcept { return m_w_factory_fc.get(); }
    #else
    //! Return the ERF auxiliary x-face EB factory.
    eb_aux_ const* get_u_const_factory() const noexcept { return &m_u_factory; }
    //! Return the ERF auxiliary y-face EB factory.
    eb_aux_ const* get_v_const_factory() const noexcept { return &m_v_factory; }
    //! Return the ERF auxiliary z-face EB factory.
    eb_aux_ const* get_w_const_factory() const noexcept { return &m_w_factory; }
    #endif

    //! Forward declaration for EB surface output.
    class EBToPVD;

  private:

    int m_has_eb;
    std::string m_type;

    amrex::EBSupport m_support_level;

    int m_write_eb_surface;

    amrex::FabArray<amrex::EBCellFlagFab>* m_cellflags = nullptr;

    //! EB level constructed from building GeometryShop
    amrex::EB2::Level const* m_eb_level;

    std::unique_ptr<amrex::EBFArrayBoxFactory> m_factory = nullptr;

    // Original eb_aux_ factories (USE_FC_FACTORY=0)
    eb_aux_ m_u_factory;
    eb_aux_ m_v_factory;
    eb_aux_ m_w_factory;

    // New native AMReX FC factories (USE_FC_FACTORY=1)
    std::unique_ptr<amrex::EBFArrayBoxFactory> m_u_factory_fc = nullptr;
    std::unique_ptr<amrex::EBFArrayBoxFactory> m_v_factory_fc = nullptr;
    std::unique_ptr<amrex::EBFArrayBoxFactory> m_w_factory_fc = nullptr;

    /**
     * \brief Build terrain EB geometry for the supplied Geometry.
     * \param a_geom Geometry used to define the EB index space.
     */
    void make_terrain (amrex::Geometry const& a_geom);

    /**
     * \brief Construct EB levels from an AMReX GeometryShop.
     * \tparam F Implicit-function type stored by the GeometryShop.
     * \param a_geom Geometry used to build the EB index space.
     * \param a_gshop GeometryShop describing the implicit EB.
     */
    template<class F>
    void build_level (amrex::Geometry const& a_geom,
                      amrex::EB2::GeometryShop<F> a_gshop)
    {
      int const req_lev(0);
      int const max_lev(2);

      amrex::EB2::Build(a_gshop, a_geom, req_lev, max_lev);
      const amrex::EB2::IndexSpace& ebis = amrex::EB2::IndexSpace::top();
      m_eb_level = &(ebis.getLevel(a_geom));
    }

    //! Return mutable EB cell flags from an AMReX factory.
    inline amrex::FabArray<amrex::EBCellFlagFab>&
    getNonConstEBCellFlags(const amrex::EBFArrayBoxFactory& ebfact)
    {
      const amrex::FabArray<amrex::EBCellFlagFab>& flags_const = ebfact.getMultiEBCellFlagFab();
      return const_cast<amrex::FabArray<amrex::EBCellFlagFab>&>(flags_const);
    }

    //! Return mutable volume fractions from an AMReX factory.
    inline amrex::MultiFab&
    getNonConstVolFrac(const amrex::EBFArrayBoxFactory& ebfact)
    {
      const amrex::MultiFab& vfrac_const = ebfact.getVolFrac();
      return const_cast<amrex::MultiFab&>(vfrac_const);
    }

    //! Return mutable cell centroids from an AMReX factory.
    inline amrex::MultiCutFab&
    getNonConstCentroid(const amrex::EBFArrayBoxFactory& ebfact)
    {
      const amrex::MultiCutFab& vcent_const = ebfact.getCentroid();
      return const_cast<amrex::MultiCutFab&>(vcent_const);
    }

    //! Return mutable area-fraction arrays from an AMReX factory.
    inline amrex::Array<amrex::MultiCutFab*, AMREX_SPACEDIM>
    getNonConstAreaFrac(const amrex::EBFArrayBoxFactory& ebfact)
    {
      auto afrac_const = ebfact.getAreaFrac();
      amrex::Array<amrex::MultiCutFab*, AMREX_SPACEDIM> afrac;
      for (int dir = 0; dir < AMREX_SPACEDIM; ++dir) {
          afrac[dir] = const_cast<amrex::MultiCutFab*>(afrac_const[dir]);
      }
      return afrac;
    }

    //! Return mutable face-centroid arrays from an AMReX factory.
    inline amrex::Array<amrex::MultiCutFab*, AMREX_SPACEDIM>
    getNonConstFaceCent(const amrex::EBFArrayBoxFactory& ebfact)
    {
      auto fcent_const = ebfact.getFaceCent();
      amrex::Array<amrex::MultiCutFab*, AMREX_SPACEDIM> fcent;
      for (int dir = 0; dir < AMREX_SPACEDIM; ++dir) {
          fcent[dir] = const_cast<amrex::MultiCutFab*>(fcent_const[dir]);
      }
      return fcent;
    }

};
#endif
