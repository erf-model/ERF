/**
 * \file ERF_MLTerrainPoisson.H
 */
#ifndef ERF_MLTERRAINPOISSON_H_
#define ERF_MLTERRAINPOISSON_H_

#include <AMReX_MLCellLinOp.H>
#include <AMReX_MultiFab.H>
#include <AMReX_iMultiFab.H>
#include <AMReX_LayoutData.H>

#include "ERF_MLTerrainPoisson_K.H"

/**
 * Geometric multigrid form of the terrain-fitted anelastic Poisson operator.
 *
 * The operator is the 15-point stencil of ERF_TerrainPoisson_3D_K.H, the same one
 * the FFT-preconditioned GMRES solver applies, with the full metric terms (h_xi and
 * h_eta cross derivatives included).  Each multigrid level carries its own nodal
 * terrain surface, obtained by keeping every ratio-th node of the finer level, and
 * the face areas and Jacobian recomputed from it with the ERF metric routines, so
 * the coarse operators discretize the same terrain rather than an averaged one.
 *
 * The smoother is a red-black (by column parity) z-line Gauss-Seidel: every column
 * segment of a box is relaxed with the exact tridiagonal block of the operator,
 * whose entries are read off the stencil itself with a ghost-fill-aware unit vector.
 * The operator is scaled by the inverse Jacobian, so its left null vector is dJ and
 * the solvability offset is the dJ-weighted mean of the right-hand side.
 *
 * Boundary conditions: periodic, homogeneous Neumann and homogeneous Dirichlet, the
 * last filled at second order (ghost = -interior) like the GMRES operator; the box
 * faces inside the domain (coarse/fine boundary of a fine level) use the
 * coarse/fine type set with setCoarseFineBC.  The ghost cells are filled by this
 * class's own applyBC because the stencil needs the x-z and y-z edge ghosts.
 */
class MLTerrainPoisson : public amrex::MLCellLinOp
{
public:

    using MF  = amrex::MultiFab;
    using FAB = amrex::FArrayBox;
    using RT  = amrex::Real;

    /**
     * Construct the operator on one AMR level.
     *
     * @param a_geom Geometry of the level (one entry)
     * @param a_grids Boxes to solve on (one entry)
     * @param a_dmap Distribution of the boxes (one entry)
     * @param a_info Multigrid options (hidden direction, agglomeration, ...)
     * @param z_phys_nd Node-centered terrain height on a_grids, at least one ghost node
     * @param ax Metric coefficient on x-faces (map-factor scaled as the projection uses it)
     * @param ay Metric coefficient on y-faces
     * @param az Metric coefficient on z-faces
     * @param dJ Cell-centered Jacobian determinant
     */
    MLTerrainPoisson (const amrex::Vector<amrex::Geometry>& a_geom,
                      const amrex::Vector<amrex::BoxArray>& a_grids,
                      const amrex::Vector<amrex::DistributionMapping>& a_dmap,
                      const amrex::LPInfo& a_info,
                      const amrex::MultiFab& z_phys_nd,
                      const amrex::MultiFab& ax,
                      const amrex::MultiFab& ay,
                      const amrex::MultiFab& az,
                      const amrex::MultiFab& dJ);

    ~MLTerrainPoisson () override = default;

    MLTerrainPoisson (const MLTerrainPoisson&) = delete;
    MLTerrainPoisson (MLTerrainPoisson&&) = delete;
    MLTerrainPoisson& operator= (const MLTerrainPoisson&) = delete;
    MLTerrainPoisson& operator= (MLTerrainPoisson&&) = delete;

    [[nodiscard]] std::string name () const override { return std::string("MLTerrainPoisson"); }

    //! The stencil reads edge ghost cells, so the base-class cross-stencil fill does not apply.
    [[nodiscard]] bool isCrossStencil () const override { return false; }

    /**
     * Restrict the residual with the fine Jacobian as weight (the conservative
     * choice for an operator scaled by 1/dJ) instead of the plain cell average.
     *
     * @param flag True for the weighted restriction
     */
    void setWeightedRestriction (bool flag) noexcept { m_weighted_restriction = flag; }

    void prepareForSolve () override;

    [[nodiscard]] bool isSingular (int amrlev) const override { return m_is_singular[amrlev]; }
    [[nodiscard]] bool isBottomSingular () const override { return m_is_singular[0]; }

    void Fapply (int amrlev, int mglev, MF& out, const MF& in) const override;
    void Fsmooth (int amrlev, int mglev, MF& sol, const MF& rhs, int redblack) const override;
    void FFlux (int amrlev, const amrex::MFIter& mfi,
                const amrex::Array<FAB*,AMREX_SPACEDIM>& flux,
                const FAB& sol, Location loc, int face_only=0) const override;

    void applyBC (int amrlev, int mglev, MF& in, BCMode bc_mode, StateMode s_mode,
                  const amrex::MLMGBndry* bndry=nullptr, bool skip_fillboundary=false) const override;

    void restriction (int amrlev, int cmglev, MF& crse, MF& fine) const override;

    //! Face fluxes of a solution (what MLMG::getFluxes returns): the stencil's -grad terms
    void getFluxes (const amrex::Vector<amrex::Array<MF*,AMREX_SPACEDIM>>& a_flux,
                    const amrex::Vector<MF*>& a_sol, Location a_loc) const override;

    [[nodiscard]] amrex::Vector<RT> getSolvabilityOffset (int amrlev, int mglev, MF const& rhs) const override;
    void fixSolvabilityByOffset (int amrlev, int mglev, MF& rhs, amrex::Vector<RT> const& offset) const override;

    //! Number of multigrid levels built on the AMR level
    [[nodiscard]] int numMGLevels () const noexcept { return NMGLevels(0); }

    //! Nodal terrain height of a multigrid level (level 0 is the one passed in)
    [[nodiscard]] const amrex::MultiFab& zPhysNd (int mglev) const { return m_zphys[mglev]; }
    //! x-face metric coefficient of a multigrid level
    [[nodiscard]] const amrex::MultiFab& axCoef (int mglev) const { return m_ax[mglev]; }
    //! y-face metric coefficient of a multigrid level
    [[nodiscard]] const amrex::MultiFab& ayCoef (int mglev) const { return m_ay[mglev]; }
    //! z-face metric coefficient of a multigrid level
    [[nodiscard]] const amrex::MultiFab& azCoef (int mglev) const { return m_az[mglev]; }
    //! Jacobian determinant of a multigrid level
    [[nodiscard]] const amrex::MultiFab& detJ (int mglev) const { return m_dJ[mglev]; }
    //! Column tridiagonal coefficients (sub, main, super) of a multigrid level; built by prepareForSolve
    [[nodiscard]] const amrex::MultiFab& triDiag (int mglev) const { return m_tri[mglev]; }
    //! Uncovered-ghost-cell mask of a multigrid level (1 where applyBC fills)
    [[nodiscard]] const amrex::iMultiFab& ghostMask (int mglev) const { return m_mask[mglev]; }
    //! Source cell (3 components) each uncovered ghost cell of a level is filled from
    [[nodiscard]] const amrex::iMultiFab& ghostSource (int mglev) const { return m_src[mglev]; }
    //! Sign each uncovered ghost cell of a level carries relative to its source
    [[nodiscard]] const amrex::MultiFab& ghostSign (int mglev) const { return m_sgn[mglev]; }

    /**
     * Number of times a box can be coarsened by 2 in each direction while keeping at
     * least min_width cells there (0 in the hidden direction).  The multigrid coarsens
     * every direction together as long as all allow it, then the remaining ones alone
     * (semicoarsening), so the hierarchy is as deep as the deepest direction.
     *
     * @param b Box to test
     * @param hidden_dir Hidden direction (-1 for none), never coarsened
     * @param min_width Smallest extent allowed after coarsening
     * @return Coarsening depth per direction
     */
    static amrex::IntVect coarsening_depth (const amrex::Box& b, int hidden_dir, int min_width = 2);

    /**
     * Number of times every box of a BoxArray can be coarsened by 2 in each direction
     * (see the Box version).
     *
     * @param ba Boxes to test
     * @param hidden_dir Hidden direction (-1 for none)
     * @param min_width Smallest extent allowed after coarsening
     * @return Coarsening depth per direction
     */
    static amrex::IntVect coarsening_depth (const amrex::BoxArray& ba, int hidden_dir, int min_width = 2);

    /**
     * Boxes covering a domain whose corners and extents are multiples of 2^depth[d] in
     * each direction d, so that the multigrid can coarsen them depth[d] times, and no
     * longer than about max_size cells in a direction that is coarsened.
     *
     * @param domain Box to cover
     * @param hidden_dir Hidden direction (-1 for none), never split
     * @param depth Coarsenings the boxes must allow per direction (at most coarsening_depth(domain))
     * @param max_size Target box extent in cells
     * @return The BoxArray
     */
    static amrex::BoxArray multigrid_grids (const amrex::Box& domain, int hidden_dir,
                                            const amrex::IntVect& depth, int max_size);

    /**
     * Copy every ratio-th cell or node of fine into crse (valid region only).
     *
     * crse(I,J,K) = fine(ratio*I, ratio*J, ratio*K); works for any index type as long as
     * crse.boxArray() covers the same region as fine.boxArray() coarsened by ratio.
     *
     * @param fine Source field
     * @param crse Destination field
     * @param ratio Coarsening ratio per direction
     */
    static void subsample (const amrex::MultiFab& fine, amrex::MultiFab& crse, const amrex::IntVect& ratio);

    /**
     * Fill the ghost nodes of a nodal terrain field that no box of its BoxArray
     * covers, the way ERF fills them on the finest level: copy the nearest node
     * inside the box in a lateral direction, extrapolate linearly below the lowest
     * and above the highest node.  Covered ghost nodes are filled by FillBoundary.
     *
     * @param znd Nodal terrain height, ghost nodes updated
     * @param geom Geometry giving the periodicity
     */
    static void fill_uncovered_nodes (amrex::MultiFab& znd, const amrex::Geometry& geom);

    // The next three are internal, but public because they launch device lambdas, which
    // nvcc does not accept inside private or protected member functions.

    //! Build the ghost mask and the ghost source map of one level (needs the domain BCs)
    void define_ghost_fill (int mglev);

    //! Build the coarse-level terrain, metrics and ghost masks
    void define_levels (const amrex::MultiFab& z_phys_nd,
                        const amrex::MultiFab& ax, const amrex::MultiFab& ay,
                        const amrex::MultiFab& az, const amrex::MultiFab& dJ);

    //! Probe the stencil for the column tridiagonal coefficients of one level
    void compute_tridiagonal (int mglev);

private:

    //! Signs the ghost fill applies on each face of a level
    [[nodiscard]] TerrainGhostRule ghost_rule (int mglev) const;

    //! Index value in the hidden direction (0 if there is none)
    [[nodiscard]] int hidden_index (int mglev) const;

    //! Domain length in each periodic direction, 0 in the others
    [[nodiscard]] amrex::Dim3 periodic_lengths (int mglev) const;

    amrex::Vector<amrex::MultiFab>  m_zphys;   //!< nodal terrain height per multigrid level
    amrex::Vector<amrex::MultiFab>  m_ax;      //!< x-face metric coefficient per level
    amrex::Vector<amrex::MultiFab>  m_ay;      //!< y-face metric coefficient per level
    amrex::Vector<amrex::MultiFab>  m_az;      //!< z-face metric coefficient per level
    amrex::Vector<amrex::MultiFab>  m_dJ;      //!< Jacobian per level
    amrex::Vector<amrex::MultiFab>  m_tri;     //!< column tridiagonal (sub, main, super) per level
    amrex::Vector<amrex::iMultiFab> m_mask;    //!< 1 on uncovered ghost cells, per level
    amrex::Vector<amrex::iMultiFab> m_src;     //!< source cell of each uncovered ghost cell (3 comps), per level
    amrex::Vector<amrex::MultiFab>  m_sgn;     //!< sign of each uncovered ghost cell, per level
    //! Per box: 1 if it touches a coarse/fine face (the flux through those faces is zero)
    amrex::Vector<amrex::LayoutData<int>> m_box_has_cf;

    //! Map-factor scale of ax, ay, az relative to the plain metric areas, per level
    amrex::Vector<amrex::MultiFab>  m_sx, m_sy, m_sz;
    bool m_has_scale = false;

    amrex::Vector<int> m_is_singular;
    bool m_coeffs_ready = false;
    bool m_weighted_restriction = false;
};

#endif
