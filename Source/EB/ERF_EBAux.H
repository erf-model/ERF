#ifndef ERF_EB_AUX_H_
#define ERF_EB_AUX_H_

#include <AMReX_Geometry.H>
#include <AMReX_DistributionMapping.H>
#include <AMReX_BoxArray.H>

#include <AMReX_EB2.H>
#include <AMReX_EBFabFactory.H>

class eb_aux_
{
  public:

    eb_aux_ ();
    ~eb_aux_ ();

    void define ( int const& a_idim,
                  amrex::Geometry            const& a_geom,
                  amrex::BoxArray            const& a_grids,
                  amrex::DistributionMapping const& a_dmap,
                  amrex::Vector<int>         const& a_ngrow,
                  amrex::EBFArrayBoxFactory  const* a_factory);

    void set_verbose ( ) { m_verbose = 1; }

    [[nodiscard]] const amrex::MultiFab& getVolFrac () const;
    [[nodiscard]] const amrex::MultiCutFab& getVolCent () const;
    [[nodiscard]] const amrex::MultiCutFab& getBndryArea () const;
    [[nodiscard]] const amrex::MultiCutFab& getBndryCent () const;
    [[nodiscard]] const amrex::MultiCutFab& getBndryNorm () const;

    [[nodiscard]] amrex::Array<const amrex::MultiCutFab*, AMREX_SPACEDIM> getAreaFrac () const;
    [[nodiscard]] amrex::Array<const amrex::MultiCutFab*, AMREX_SPACEDIM> getFaceCent () const;

  private:

    int m_verbose;

    // int m_defined;

    amrex::FabArray<amrex::EBCellFlagFab>* m_cellflags = nullptr;

    amrex::MultiFab* m_volfrac = nullptr;

    amrex::MultiCutFab* m_volcent = nullptr;
    amrex::MultiCutFab* m_bndryarea = nullptr;
    amrex::MultiCutFab* m_bndrycent = nullptr;
    amrex::MultiCutFab* m_bndrynorm = nullptr;

    amrex::Array<amrex::MultiCutFab*,AMREX_SPACEDIM> m_areafrac {{AMREX_D_DECL(nullptr, nullptr, nullptr)}};
    amrex::Array<amrex::MultiCutFab*,AMREX_SPACEDIM> m_facecent {{AMREX_D_DECL(nullptr, nullptr, nullptr)}};
    // amrex::Array<amrex::MultiCutFab*,AMREX_SPACEDIM> m_edgecent {{AMREX_D_DECL(nullptr, nullptr, nullptr)}};

};

#endif
