#ifndef ERF_EB_H_
#define ERF_EB_H_

#include <AMReX_Geometry.H>
#include <AMReX_DistributionMapping.H>
#include <AMReX_BoxArray.H>

#include <AMReX_EB2.H>
#include <AMReX_EBFabFactory.H>

#include <ERF_EBAux.H>

class eb_ {

  public:

    ~eb_ ();

    eb_ (amrex::Geometry const& a_geom,
         amrex::FArrayBox const& terrain_fab,
         amrex::Gpu::DeviceVector<amrex::Real>& a_dz_stretched,
         bool is_anelastic);
    eb_ ();

    void define (int level,
                 amrex::Geometry const& a_geom,
                 amrex::EB2::Level const* a_eb_level,
                 bool is_anelastic);

    void make_all_factories ( int level,
                        amrex::Geometry   const& a_geom,
                        amrex::BoxArray            const& ba,
                        amrex::DistributionMapping const& dm,
                        amrex::EB2::Level const& a_eb_level);

    void make_cc_factory ( int level,
                        amrex::Geometry   const& a_geom,
                        amrex::BoxArray            const& ba,
                        amrex::DistributionMapping const& dm,
                        amrex::EB2::Level const& a_eb_level);

    int nghost_basic () const  { return 5; } // nghost_eb_basic ()
    int nghost_volume () const { return 5; } // nghost_eb_volume ()
    int nghost_full () const   { return 4; } // nghost_eb_full ()

    const std::unique_ptr<amrex::EBFArrayBoxFactory>& get_const_factory () const noexcept { return m_factory; }

    void set_connection_flags ();

    eb_aux_ const* get_u_const_factory() const noexcept { return &m_u_factory; }
    eb_aux_ const* get_v_const_factory() const noexcept { return &m_v_factory; }
    eb_aux_ const* get_w_const_factory() const noexcept { return &m_w_factory; }

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

    eb_aux_ m_u_factory;
    eb_aux_ m_v_factory;
    eb_aux_ m_w_factory;

    void make_terrain (amrex::Geometry const& a_geom);

    //! Construct EB levels from Geometry shop.
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

    // Wrapper to get non-const version of EBCellFlagFab using const_cast
    inline amrex::FabArray<amrex::EBCellFlagFab>&
    getNonConstEBCellFlags(const amrex::EBFArrayBoxFactory& ebfact)
    {
      const amrex::FabArray<amrex::EBCellFlagFab>& flags_const = ebfact.getMultiEBCellFlagFab();
      return const_cast<amrex::FabArray<amrex::EBCellFlagFab>&>(flags_const);
    }

    inline amrex::MultiFab&
    getNonConstVolFrac(const amrex::EBFArrayBoxFactory& ebfact)
    {
      const amrex::MultiFab& vfrac_const = ebfact.getVolFrac();
      return const_cast<amrex::MultiFab&>(vfrac_const);
    }

    inline amrex::MultiCutFab&
    getNonConstCentroid(const amrex::EBFArrayBoxFactory& ebfact)
    {
      const amrex::MultiCutFab& vcent_const = ebfact.getCentroid();
      return const_cast<amrex::MultiCutFab&>(vcent_const);
    }

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
