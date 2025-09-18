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

    void define ( int const& a_level,
                  int const& a_idim,
                  amrex::Geometry            const& a_geom,
                  amrex::BoxArray            const& a_grids,
                  amrex::DistributionMapping const& a_dmap,
                  amrex::Vector<int>         const& a_ngrow,
                  amrex::EBFArrayBoxFactory  const* a_factory);

    void set_verbose ( ) { m_verbose = 1; }

    [[nodiscard]] const amrex::FabArray<amrex::EBCellFlagFab>& getMultiEBCellFlagFab () const;
    [[nodiscard]] const amrex::MultiFab& getVolFrac () const;
    [[nodiscard]] const amrex::MultiFab& getCentroid () const;
    [[nodiscard]] const amrex::MultiFab& getBndryArea () const;
    [[nodiscard]] const amrex::MultiFab& getBndryCent () const;
    [[nodiscard]] const amrex::MultiFab& getBndryNorm () const;

    [[nodiscard]] amrex::Array<const amrex::MultiFab*, AMREX_SPACEDIM> getAreaFrac () const;
    [[nodiscard]] amrex::Array<const amrex::MultiFab*, AMREX_SPACEDIM> getFaceCent () const;

  private:

    int m_verbose;

    // int m_defined;

    amrex::FabArray<amrex::EBCellFlagFab>* m_cellflags = nullptr;
    amrex::MultiFab* m_volfrac = nullptr;
    amrex::MultiFab* m_volcent = nullptr;
    amrex::MultiFab* m_bndryarea = nullptr;
    amrex::MultiFab* m_bndrycent = nullptr;
    amrex::MultiFab* m_bndrynorm = nullptr;

    amrex::Array<amrex::MultiFab*,AMREX_SPACEDIM> m_areafrac {{AMREX_D_DECL(nullptr, nullptr, nullptr)}};
    amrex::Array<amrex::MultiFab*,AMREX_SPACEDIM> m_facecent {{AMREX_D_DECL(nullptr, nullptr, nullptr)}};
    // amrex::Array<amrex::MultiCutFab*,AMREX_SPACEDIM> m_edgecent {{AMREX_D_DECL(nullptr, nullptr, nullptr)}};

};

#endif
