/**
 * \file ERF_EBAux.H
 * \brief Declares auxiliary EB factories for face-centered data.
 */
#ifndef ERF_EB_AUX_H_
#define ERF_EB_AUX_H_

#include <memory>

#include <AMReX_Geometry.H>
#include <AMReX_DistributionMapping.H>
#include <AMReX_BoxArray.H>

#include <AMReX_EB2.H>
#include <AMReX_EBFabFactory.H>

/**
 * \brief Stores face-centered EB geometry reconstructed from a cell-centered factory.
 *
 * This class mirrors the subset of the AMReX EB factory interface needed by ERF
 * staggered momentum fields when native face-centered EB factories are disabled.
 */
class eb_aux_
{
  public:

    //! Construct an undefined auxiliary EB factory.
    eb_aux_ ();
    //! Destroy the auxiliary EB factory; the owned geometry is released by its unique_ptrs.
    ~eb_aux_ () = default;

    /**
     * \brief Define face-centered EB geometry for one staggered direction.
     *
     * May be called more than once on the same object -- ERF::RemakeLevel and
     * MakeNewLevelFromCoarse both redefine the factories of an existing level.
     * Each call replaces the previously owned geometry, which is freed on
     * reassignment of the owning unique_ptr.
     * \param a_level AMR level index.
     * \param a_idim Staggered direction for the target velocity grid.
     * \param a_geom Cell-centered level geometry.
     * \param a_grids Cell-centered grid boxes.
     * \param a_dmap Distribution mapping for the grids.
     * \param a_ngrow Ghost-cell counts for flags, volumes, and centroids.
     * \param a_factory Source cell-centered EB factory.
     */
    void define ( int const& a_level,
                  int const& a_idim,
                  amrex::Geometry            const& a_geom,
                  amrex::BoxArray            const& a_grids,
                  amrex::DistributionMapping const& a_dmap,
                  amrex::Vector<int>         const& a_ngrow,
                  amrex::EBFArrayBoxFactory  const* a_factory);

    //! Enable verbose diagnostic output during factory construction.
    void set_verbose ( ) { m_verbose = 1; }

    //! Return the reconstructed EB cell flags.
    [[nodiscard]] const amrex::FabArray<amrex::EBCellFlagFab>& getMultiEBCellFlagFab () const;
    //! Return the reconstructed volume fractions.
    [[nodiscard]] const amrex::MultiFab& getVolFrac () const;
    //! Return the reconstructed volume centroids.
    [[nodiscard]] const amrex::MultiFab& getCentroid () const;
    //! Return the reconstructed EB boundary areas.
    [[nodiscard]] const amrex::MultiFab& getBndryArea () const;
    //! Return the reconstructed EB boundary centroids.
    [[nodiscard]] const amrex::MultiFab& getBndryCent () const;
    //! Return the reconstructed EB boundary normals.
    [[nodiscard]] const amrex::MultiFab& getBndryNormal () const;

    //! Return reconstructed area fractions for all coordinate directions.
    [[nodiscard]] amrex::Array<const amrex::MultiFab*, AMREX_SPACEDIM> getAreaFrac () const;
    //! Return reconstructed face centroids for all coordinate directions.
    [[nodiscard]] amrex::Array<const amrex::MultiFab*, AMREX_SPACEDIM> getFaceCent () const;

  private:

    int m_verbose;

    // int m_defined;

    // These own the geometry they point at.  define() may be called repeatedly on the
    // same object, so ownership has to be explicit: assigning a new unique_ptr frees
    // whatever the previous define() allocated.
    std::unique_ptr<amrex::FabArray<amrex::EBCellFlagFab>> m_cellflags;
    std::unique_ptr<amrex::MultiFab> m_volfrac;
    std::unique_ptr<amrex::MultiFab> m_volcent;
    std::unique_ptr<amrex::MultiFab> m_bndryarea;
    std::unique_ptr<amrex::MultiFab> m_bndrycent;
    std::unique_ptr<amrex::MultiFab> m_bndrynorm;

    amrex::Array<std::unique_ptr<amrex::MultiFab>,AMREX_SPACEDIM> m_areafrac;
    amrex::Array<std::unique_ptr<amrex::MultiFab>,AMREX_SPACEDIM> m_facecent;
    // amrex::Array<std::unique_ptr<amrex::MultiCutFab>,AMREX_SPACEDIM> m_edgecent;

};

#endif
