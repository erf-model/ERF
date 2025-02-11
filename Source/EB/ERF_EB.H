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

    void make_factory ( amrex::Geometry            const& a_geom,
                        amrex::DistributionMapping const& a_dmap,
                        amrex::BoxArray            const& a_grids);

    int nghost_basic () const  { return 2; }
    int nghost_volume () const { return 2; }
    int nghost_full () const   { return 2; }

    amrex::EBFArrayBoxFactory const* get_const_factory () noexcept
    { return m_factory.get(); }

    amrex::EB2::Level const* get_level () const noexcept
    { return m_eb_level; }

    void WriteEBSurface (const amrex::BoxArray & ba,
                        const amrex::DistributionMapping & dmap,
                        const amrex::Geometry & geom,
                        const amrex::EBFArrayBoxFactory * ebf);

    eb_aux_ const* get_u_const_factory() { return &m_u_factory; }
    eb_aux_ const* get_v_const_factory() { return &m_v_factory; }
    eb_aux_ const* get_w_const_factory() { return &m_w_factory; }

    class EBToPVD;

  private:

    int m_has_eb;
    std::string m_type;

    amrex::EBSupport m_support_level;

    int m_write_eb_surface;

    amrex::FabArray<amrex::EBCellFlagFab>* m_cellflags = nullptr;

    //! EB level constructed from building GeometryShop
    amrex::EB2::Level const* m_eb_level;

    std::unique_ptr<amrex::EBFArrayBoxFactory> m_factory;

    eb_aux_ m_u_factory;
    eb_aux_ m_v_factory;
    eb_aux_ m_w_factory;

    void make_box     (amrex::Geometry const& a_geom);
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

};
#endif
