#ifndef ERF_NODAL_RECONSTRUCTION_H_
#define ERF_NODAL_RECONSTRUCTION_H_

#include <algorithm>
#include <cmath>
#include <limits>
#include <stdexcept>
#include <utility>

#include "AMReX_Array4.H"
#include "AMReX_Box.H"
#include "AMReX_FArrayBox.H"
#include "AMReX_Geometry.H"
#include "AMReX_REAL.H"

#include "ERF_Constants.H"

enum class VariationOperator
{
    FirstDeriv,
    Laplacian
};

/**
 * @brief Convergence and diagnostic information for the reconstruction solve.
 */
struct SolveInfo {
    int  iterations           = 0;     //!< total CG iterations, summed over attempts
    amrex::Real final_residual = zero; //!< CG residual of the accepted solve
    bool converged            = false;
    int  refinements          = 0;     //!< times the regularization had to be raised
    amrex::Real regularization = zero; //!< accepted value of mu

    // Diagnostics on the reconstructed field
    amrex::Real deviation            = zero; //!< max |S - S_ref|
    amrex::Real deviation_cap        = zero; //!< bound imposed on the above
    amrex::Real max_average_error    = zero; //!< max |avg4(S) - T|
    amrex::Real interp_average_error = zero; //!< max |avg4(S_ref) - T|, for comparison
    amrex::Real min_value            = zero; //!< min over nodes of S
    amrex::Real max_value            = zero; //!< max over nodes of S
};

/**
 * Reconstruct a nodal field S from data T that lives at cell centers in (x,y),
 * such that the four-node average of S matches T as closely as a *bounded*,
 * smooth nodal field can.
 *
 * ---------------------------------------------------------------------------
 * Why we do not simply invert the averaging operator
 * ---------------------------------------------------------------------------
 * The four-node averaging operator A,
 *
 *     (A S)(i,j) = 1/4 [ S(i,j) + S(i+1,j) + S(i,j+1) + S(i+1,j+1) ] ,
 *
 * has symbol cos(kx/2) cos(ky/2), which vanishes at the grid Nyquist mode.
 * Inverting A exactly -- e.g. with the back substitution
 *
 *     S(i+1,j+1) = 4 T(i,j) - S(i,j) - S(i+1,j) - S(i,j+1)
 *
 * -- therefore amplifies grid-scale content of T without bound.  The Green's
 * function of that recursion is 4 / [(1+x)(1+y)], a doubly alternating
 * cumulative sum, so a checkerboard component of T of amplitude a produces a
 * nodal response of amplitude a*i*j.  On a 128^2 grid an O(10 m) grid-scale
 * roughness in the WRF terrain thus yields nodal heights that are kilometers
 * away from the WRF terrain, including physically impossible negative
 * elevations (ERF issue #3709).  Adding a free homogeneous solution cannot
 * repair this: the null space of A is spanned by (-1)^i f(j) + (-1)^j g(i),
 * which cannot represent the (-1)^(i+j) i j growth of the particular solution.
 * The failure is silent, because checking that avg4(S) reproduces T is then
 * satisfied by construction.
 *
 * ---------------------------------------------------------------------------
 * What we do instead
 * ---------------------------------------------------------------------------
 * Let S_ref be the direct interpolation of T to the nodes (makeReference) and
 * R = T - A S_ref the part of T that interpolation fails to reproduce.  We
 * solve for a *correction* E to the interpolant,
 *
 *     min_E  || A E - R ||^2 + lambda || V E ||^2 + mu || E ||^2 ,   S = S_ref + E,
 *
 * where V is a roughness operator (first derivative or Laplacian).  The normal
 * equations
 *
 *     ( A^T A + lambda V^T V + mu I ) E = A^T R
 *
 * are symmetric positive definite for any mu > 0, so E exists, is unique, and
 * -- unlike the exact inverse -- is bounded.  Specifically, writing
 * ( A^T A + lambda V^T V + mu I ) E = A^T R and pairing with E gives
 *
 *     || A E ||^2 + lambda || V E ||^2 + mu || E ||^2 <= || A E || || R || ,
 *
 * and since 2 sqrt(mu) ||A E|| ||E|| <= ||A E||^2 + mu ||E||^2, we get the
 * a priori bound
 *
 *     || E || <= || R || / ( 2 sqrt(mu) ) .
 *
 * The correction is therefore controlled by how badly interpolation already
 * fails, never by the conditioning of A.
 *
 * ---------------------------------------------------------------------------
 * Choosing mu
 * ---------------------------------------------------------------------------
 * A single fixed mu cannot serve both smooth and rough terrain: small mu
 * recovers well-resolved features almost exactly but lets rough fields drift
 * far from the interpolant, and large mu does the reverse.  We therefore start
 * from a small mu and raise it until the correction satisfies the pointwise cap
 *
 *     max |E| <= min( 2 max|R| , 1/4 (max T - min T) ) ,
 *
 * i.e. the nodes may not move away from the direct interpolation by more than
 * twice the interpolation misfit, nor by more than a quarter of the relief of
 * the data.  Because ||E|| -> 0 as mu -> infinity the loop always terminates,
 * and on failure we fall back to the interpolant itself.  Smooth terrain keeps
 * the smallest mu and is reproduced to well below the interpolation error;
 * rough terrain is pulled back toward the interpolant instead of blowing up.
 *
 * The averaging error max |avg4(S) - T| returned in SolveInfo is now a genuine
 * diagnostic -- it is no longer satisfied by construction -- and is reported
 * alongside the same quantity for the plain interpolant and the range of S, so
 * callers can sanity check the result.
 */
class NodalReconstruction
{
public:
    using Real = amrex::Real;

    /**
     * @brief Construct a nodal reconstruction object for a given face box.
     * @param face_box Box specifying the cell-centered data region.
     * @param geom Geometry used for grid spacing and anisotropy weights.
     */
    explicit NodalReconstruction (amrex::Box const& face_box,
                                  amrex::Geometry const& geom)
        : m_face_box(face_box),
          m_ilo(face_box.smallEnd(0)),
          m_ihi(face_box.bigEnd(0)),
          m_jlo(face_box.smallEnd(1)),
          m_jhi(face_box.bigEnd(1))
    {
        if (face_box.length(2) != 1) {
            throw std::invalid_argument(
                "face_box must contain exactly one z cell.");
        }
        if (face_box.smallEnd(2) != 0) {
            throw std::invalid_argument(
                "face_box must be a slab at k = 0; all fabs are indexed at k = 0.");
        }

        m_nodal_box = amrex::surroundingNodes(face_box);
        m_nodal_box.setSmall(2, 0);
        m_nodal_box.setBig(2, 0);

        Real const dx = geom.CellSizeArray()[0];
        Real const dy = geom.CellSizeArray()[1];
        if (!(dx > zero) || !(dy > zero)) {
            throw std::invalid_argument("Cell sizes must be positive.");
        }

        // Non-dimensional, O(1) stencil weights that respect grid anisotropy.
        // They are unity on an isotropic grid, so the smoothing weight has the
        // same meaning regardless of the mesh spacing.
        Real const h = std::min(dx,dy);
        m_gx  = h / dx;
        m_gy  = h / dy;
        m_gxx = m_gx * m_gx;
        m_gyy = m_gy * m_gy;
    }

    [[nodiscard]] amrex::Box const& nodalBox () const noexcept {
        return m_nodal_box;
    }

    /**
     * Direct interpolation of the cell-centered data to the nodes: the average
     * of the (up to four) cells that touch each node, with edge and corner
     * nodes taking the average of the cells that exist.  This is the field the
     * solve corrects, and the field it falls back to.
     */
    [[nodiscard]] amrex::FArrayBox makeReference (amrex::FArrayBox const& T) const
    {
        /**
         * @brief Local storage for the reference nodal field.
         */
        amrex::FArrayBox R(m_nodal_box,1,amrex::The_Managed_Arena());
        auto const T_arr = T.const_array();
        auto const R_arr = R.array();

        for (int j=m_jlo; j<=m_jhi+1; ++j) {
            int const jm = std::max(m_jlo, std::min(j-1, m_jhi));
            int const jp = std::max(m_jlo, std::min(j  , m_jhi));
            for (int i=m_ilo; i<=m_ihi+1; ++i) {
                int const im = std::max(m_ilo, std::min(i-1, m_ihi));
                int const ip = std::max(m_ilo, std::min(i  , m_ihi));
                R_arr(i,j,0) = fourth * ( T_arr(im,jm,0) + T_arr(ip,jm,0)
                                        + T_arr(im,jp,0) + T_arr(ip,jp,0) );
            }
        }
        return R;
    }

    /**
     * Solve the regularized least-squares problem described above.
     *
     * \param T                  cell-centered data to be matched
     * \param reference          the interpolant to correct (see makeReference)
     * \param var_op             roughness operator used as the smoothness prior
     * \param relative_tolerance CG stopping tolerance
     * \param smoothing          lambda; weight of the roughness penalty
     * \param min_regularization smallest mu tried; raised by factors of four
     *                           until the correction respects its cap
     * \param max_iterations     CG iteration cap per attempt
     */
    [[nodiscard]] std::pair<amrex::FArrayBox,SolveInfo>
    solve (amrex::FArrayBox const& T,
           amrex::FArrayBox const& reference,
           VariationOperator var_op,
           Real relative_tolerance  = Real(1.e-10),
           Real smoothing           = Real(1.e-2),
           Real min_regularization  = Real(2.5e-4),
           int  max_iterations      = 1000)
    {
        if (!(smoothing > zero)) {
            throw std::invalid_argument("Smoothing weight must be positive.");
        }
        if (!(min_regularization > zero)) {
            // With mu == 0 the operator can be singular: the roughness penalty
            // alone does not control every mode (the Laplacian penalty, for
            // instance, is blind to bilinear and discrete-harmonic fields).
            throw std::invalid_argument("Regularization must be positive.");
        }

        m_var_op = var_op;
        m_lambda = smoothing;

        m_cc_scratch.resize(m_face_box ,1,amrex::The_Managed_Arena());
        m_nd_scratch.resize(m_nodal_box,1,amrex::The_Managed_Arena());

        SolveInfo info;

        // Misfit of the plain interpolation, and the cap it earns the solve.
        /**
         * @brief Local storage for the interpolation misfit.
         */
        amrex::FArrayBox misfit(m_face_box,1,amrex::The_Managed_Arena());
        Real misfit_max = zero;
        Real t_min      =  std::numeric_limits<Real>::max();
        Real t_max      = -std::numeric_limits<Real>::max();
        {
            applyA(reference,misfit);
            auto const T_arr = T.const_array();
            auto const m_arr = misfit.array();
            for (int j=m_jlo; j<=m_jhi; ++j) {
                for (int i=m_ilo; i<=m_ihi; ++i) {
                    m_arr(i,j,0) = T_arr(i,j,0) - m_arr(i,j,0);
                    misfit_max = std::max(misfit_max, std::abs(m_arr(i,j,0)));
                    t_min      = std::min(t_min, T_arr(i,j,0));
                    t_max      = std::max(t_max, T_arr(i,j,0));
                }
            }
        }
        info.interp_average_error = misfit_max;
        info.deviation_cap = std::min(deviation_factor * misfit_max,
                                      relief_factor * (t_max - t_min));

        // rhs = A^T (T - A S_ref); the correction solves M E = rhs from E = 0.
        /**
         * @brief Local storage for the right-hand side of the normal equations.
         */
        amrex::FArrayBox rhs(m_nodal_box,1,amrex::The_Managed_Arena());
        applyAT(misfit,rhs);

        /**
         * @brief Local storage for the nodal correction.
         */
        amrex::FArrayBox E(m_nodal_box,1,amrex::The_Managed_Arena());
        /**
         * @brief Local storage for the reconstructed nodal field.
         */
        amrex::FArrayBox S(m_nodal_box,1,amrex::The_Managed_Arena());

        m_mu = min_regularization;
        bool accepted = false;
        for (int attempt=0; attempt<max_attempts; ++attempt) {
            SolveInfo attempt_info = conjugateGradient(rhs,E,relative_tolerance,max_iterations);
            info.iterations += attempt_info.iterations;
            info.final_residual = attempt_info.final_residual;
            info.converged      = attempt_info.converged;

            info.deviation = maxAbs(E);
            if (info.deviation <= info.deviation_cap) { accepted = true; break; }

            m_mu *= mu_growth;
            ++info.refinements;
        }

        if (accepted) {
            auto const S_arr = S.array();
            auto const r_arr = reference.const_array();
            auto const E_arr = E.const_array();
            for (int j=m_jlo; j<=m_jhi+1; ++j) {
                for (int i=m_ilo; i<=m_ihi+1; ++i) {
                    S_arr(i,j,0) = r_arr(i,j,0) + E_arr(i,j,0);
                }
            }
        } else {
            // Should not happen -- ||E|| -> 0 as mu -> infinity -- but if it
            // does, the interpolant is a safe, bounded answer.
            S.copy<amrex::RunOn::Host>(reference,0,0,1);
            info.deviation = zero;
        }
        info.regularization = m_mu;

        computeDiagnostics(T,S,info);
        return {std::move(S),info};
    }

    //! max_ij | avg4(S) - T |
    [[nodiscard]] Real maxAverageError (amrex::FArrayBox const& T,
                                        amrex::FArrayBox const& S) const
    {
        auto const T_arr = T.const_array();
        auto const S_arr = S.const_array();
        Real err = zero;
        for (int j=m_jlo; j<=m_jhi; ++j) {
            for (int i=m_ilo; i<=m_ihi; ++i) {
                Real const avg = fourth * ( S_arr(i  ,j  ,0) + S_arr(i+1,j  ,0)
                                          + S_arr(i  ,j+1,0) + S_arr(i+1,j+1,0) );
                err = std::max(err, std::abs(avg - T_arr(i,j,0)));
            }
        }
        return err;
    }

    [[nodiscard]] Real totalSquaredVariation (amrex::FArrayBox const& S) const
    {
        auto const S_arr = S.const_array();
        Real value = zero;

        for (int j=m_jlo; j<=m_jhi+1; ++j) {
            for (int i=m_ilo; i<=m_ihi; ++i) {
                Real const d = S_arr(i+1,j,0) - S_arr(i,j,0);
                value += d*d;
            }
        }
        for (int j=m_jlo; j<=m_jhi; ++j) {
            for (int i=m_ilo; i<=m_ihi+1; ++i) {
                Real const d = S_arr(i,j+1,0) - S_arr(i,j,0);
                value += d*d;
            }
        }
        return value;
    }

private:
    //! The correction may not move a node further from the direct interpolation
    //! than this multiple of the interpolation misfit ...
    static constexpr amrex::Real deviation_factor = amrex::Real(2.0);
    //! ... nor than this fraction of the relief of the data.
    static constexpr amrex::Real relief_factor    = amrex::Real(0.25);
    //! Factor by which mu is raised when the cap is violated, and the number of
    //! times we are willing to do so.
    static constexpr amrex::Real mu_growth        = amrex::Real(4.0);
    static constexpr int         max_attempts     = 24;

    amrex::Box m_face_box;
    amrex::Box m_nodal_box;
    int m_ilo=0;
    int m_ihi=-1;
    int m_jlo=0;
    int m_jhi=-1;
    amrex::Real m_gx  = one;
    amrex::Real m_gy  = one;
    amrex::Real m_gxx = one;
    amrex::Real m_gyy = one;
    amrex::Real m_lambda = zero;
    amrex::Real m_mu     = zero;
    VariationOperator m_var_op = VariationOperator::FirstDeriv;

    mutable amrex::FArrayBox m_cc_scratch;
    mutable amrex::FArrayBox m_nd_scratch;

    //! (A S)(i,j) = 1/4 [ S(i,j) + S(i+1,j) + S(i,j+1) + S(i+1,j+1) ]
    /**
     * @brief Apply the four-node averaging operator A.
     * @param[in] S Nodal field to be averaged.
     * @param[out] T Resulting cell-centered field.
     */
    void applyA (amrex::FArrayBox const& S,
                 amrex::FArrayBox& T) const
    {
        auto const S_arr = S.const_array();
        auto const T_arr = T.array();
        for (int j=m_jlo; j<=m_jhi; ++j) {
            for (int i=m_ilo; i<=m_ihi; ++i) {
                T_arr(i,j,0) = fourth * ( S_arr(i  ,j  ,0) + S_arr(i+1,j  ,0)
                                        + S_arr(i  ,j+1,0) + S_arr(i+1,j+1,0) );
            }
        }
    }

    //! Adjoint of applyA: scatter each cell value to its four nodes.
    /**
     * @brief Apply the adjoint of the averaging operator A.
     * @param[in] T Cell-centered field to be scattered.
     * @param[out] R Resulting nodal field.
     */
    void applyAT (amrex::FArrayBox const& T,
                  amrex::FArrayBox& R) const
    {
        R.setVal<amrex::RunOn::Host>(zero);
        auto const T_arr = T.const_array();
        auto const R_arr = R.array();
        for (int j=m_jlo; j<=m_jhi; ++j) {
            for (int i=m_ilo; i<=m_ihi; ++i) {
                Real const q = fourth * T_arr(i,j,0);
                R_arr(i  ,j  ,0) += q;
                R_arr(i+1,j  ,0) += q;
                R_arr(i  ,j+1,0) += q;
                R_arr(i+1,j+1,0) += q;
            }
        }
    }

    /**
     * @brief Compute the first derivatives of the nodal field.
     * @param[in] S Nodal field.
     * @param[out] difx Resulting x-derivatives.
     * @param[out] dify Resulting y-derivatives.
     */
    void applyD (amrex::FArrayBox const& S,
                 amrex::FArrayBox& difx,
                 amrex::FArrayBox& dify) const
    {
        difx.setVal<amrex::RunOn::Host>(zero);
        dify.setVal<amrex::RunOn::Host>(zero);

        auto const S_arr    = S.const_array();
        auto const difx_arr = difx.array();
        auto const dify_arr = dify.array();

        for (int j = m_jlo; j <= m_jhi + 1; ++j) {
            for (int i = m_ilo; i <= m_ihi; ++i) {
                difx_arr(i, j, 0) = m_gx * ( S_arr(i+1, j, 0) - S_arr(i, j, 0) );
            }
        }

        for (int j = m_jlo; j <= m_jhi; ++j) {
            for (int i = m_ilo; i <= m_ihi + 1; ++i) {
                dify_arr(i, j, 0) = m_gy * ( S_arr(i, j+1, 0) - S_arr(i, j, 0) );
            }
        }
    }

    /**
     * @brief Apply the adjoint of the derivative operator.
     * @param[in] difx X-derivatives.
     * @param[in] dify Y-derivatives.
     * @param[out] R Resulting nodal field.
     */
    void applyDT (amrex::FArrayBox const& difx,
                  amrex::FArrayBox const& dify,
                  amrex::FArrayBox& R) const
    {
        R.setVal<amrex::RunOn::Host>(zero);

        auto const qx_arr = difx.const_array();
        auto const qy_arr = dify.const_array();
        auto const R_arr  = R.array();

        for (int j = m_jlo; j <= m_jhi + 1; ++j) {
            for (int i = m_ilo; i <= m_ihi; ++i) {
                R_arr(i  , j, 0) -= m_gx * qx_arr(i, j, 0);
                R_arr(i+1, j, 0) += m_gx * qx_arr(i, j, 0);
            }
        }

        for (int j = m_jlo; j <= m_jhi; ++j) {
            for (int i = m_ilo; i <= m_ihi + 1; ++i) {
                R_arr(i, j  , 0) -= m_gy * qy_arr(i, j, 0);
                R_arr(i, j+1, 0) += m_gy * qy_arr(i, j, 0);
            }
        }
    }

    /**
     * @brief Apply the product of the derivative operator and its adjoint.
     * @param[in] S Nodal field.
     * @param[out] R Resulting nodal field.
     */
    void applyDTD (amrex::FArrayBox const& S,
                   amrex::FArrayBox& R) const
    {
        amrex::Box x_edge_box = m_nodal_box;
        amrex::Box y_edge_box = m_nodal_box;
        x_edge_box.growHi(0,-1);
        y_edge_box.growHi(1,-1);
        amrex::FArrayBox difx(x_edge_box,1,amrex::The_Managed_Arena());
        amrex::FArrayBox dify(y_edge_box,1,amrex::The_Managed_Arena());
        applyD(S, difx, dify);
        applyDT(difx, dify, R);
    }

    /**
     * @brief Apply the Laplacian operator to a nodal field.
     * @param[in] S Nodal field.
     * @param[out] lap Resulting Laplacian field.
     */
    void applyL (amrex::FArrayBox const& S,
                 amrex::FArrayBox& lap) const
    {
        lap.setVal<amrex::RunOn::Host>(zero);

        auto const S_arr   = S.const_array();
        auto const lap_arr = lap.array();

        for (int j = m_jlo + 1; j <= m_jhi; ++j) {
            for (int i = m_ilo + 1; i <= m_ihi; ++i) {
                lap_arr(i, j, 0) = m_gxx * (         S_arr(i+1, j  , 0)
                                            - two *  S_arr(i  , j  , 0)
                                                   + S_arr(i-1, j  , 0) )
                                 + m_gyy * (         S_arr(i  , j+1, 0)
                                            - two *  S_arr(i  , j  , 0)
                                                   + S_arr(i  , j-1, 0) );
            }
        }
    }

    /**
     * @brief Apply the adjoint of the Laplacian operator.
     * @param[in] lap Laplacian field.
     * @param[out] R Resulting nodal field.
     */
    void applyLT (amrex::FArrayBox const& lap,
                  amrex::FArrayBox& R) const
    {
        R.setVal<amrex::RunOn::Host>(zero);

        auto const q_arr = lap.const_array();
        auto const R_arr = R.array();

        for (int j = m_jlo + 1; j <= m_jhi; ++j) {
            for (int i = m_ilo + 1; i <= m_ihi; ++i) {
                const amrex::Real val = q_arr(i, j, 0);
                R_arr(i-1, j  , 0) += m_gxx * val;
                R_arr(i+1, j  , 0) += m_gxx * val;
                R_arr(i  , j-1, 0) += m_gyy * val;
                R_arr(i  , j+1, 0) += m_gyy * val;
                R_arr(i  , j  , 0) -= two * (m_gxx + m_gyy) * val;
            }
        }
    }

    /**
     * @brief Apply the product of the Laplacian operator and its adjoint.
     * @param[in] S Nodal field.
     * @param[out] R Resulting nodal field.
     */
    void applyLTL (amrex::FArrayBox const& S,
                   amrex::FArrayBox& R) const
    {
        amrex::Box interior_nodal_box = m_nodal_box;
        interior_nodal_box.grow(0,-1);
        interior_nodal_box.grow(1,-1);
        if (!interior_nodal_box.ok()) {
            // Fewer than two cells in a direction: there is no interior node to
            // evaluate the Laplacian at.  mu keeps the operator definite.
            R.setVal<amrex::RunOn::Host>(zero);
            return;
        }
        amrex::FArrayBox lap(interior_nodal_box,1,amrex::The_Managed_Arena());
        applyL(S, lap);
        applyLT(lap, R);
    }

    /**
     * @brief Apply the selected variation operator (first derivative or Laplacian).
     * @param[in] S Nodal field.
     * @param[out] R Resulting nodal field.
     */
    void applyVariationOperator (amrex::FArrayBox const& S,
                                 amrex::FArrayBox& R) const
    {
        switch (m_var_op) {
        case VariationOperator::FirstDeriv:
            applyDTD(S, R);
            break;

        case VariationOperator::Laplacian:
            applyLTL(S, R);
            break;

        default:
            throw std::runtime_error("Unsupported variation operator.");
        }
    }

    //! R = ( A^T A + lambda V^T V + mu I ) S
    /**
     * @brief Apply the full system operator M = A^T A + lambda V^T V + mu I.
     * @param[in] S Nodal field.
     * @param[out] R Resulting nodal field.
     */
    void applyM (amrex::FArrayBox const& S,
                 amrex::FArrayBox& R) const
    {
        applyA(S,m_cc_scratch);
        applyAT(m_cc_scratch,R);
        applyVariationOperator(S,m_nd_scratch);

        auto const R_arr = R.array();
        auto const v_arr = m_nd_scratch.const_array();
        auto const S_arr = S.const_array();
        for (int j=m_jlo; j<=m_jhi+1; ++j) {
            for (int i=m_ilo; i<=m_ihi+1; ++i) {
                R_arr(i,j,0) += m_lambda * v_arr(i,j,0) + m_mu * S_arr(i,j,0);
            }
        }
    }

    /**
     * Diagonal of M, assembled by accumulating the square of every stencil
     * coefficient that touches a node.  Used as a Jacobi preconditioner.
     */
    void computeDiagonal (amrex::FArrayBox& diag) const
    {
        diag.setVal<amrex::RunOn::Host>(m_mu);
        auto const d_arr = diag.array();

        // A^T A: each incident cell contributes (1/4)^2
        Real const qa = fourth * fourth;
        for (int j=m_jlo; j<=m_jhi; ++j) {
            for (int i=m_ilo; i<=m_ihi; ++i) {
                d_arr(i  ,j  ,0) += qa;
                d_arr(i+1,j  ,0) += qa;
                d_arr(i  ,j+1,0) += qa;
                d_arr(i+1,j+1,0) += qa;
            }
        }

        // lambda V^T V
        if (m_var_op == VariationOperator::FirstDeriv) {
            Real const qx = m_lambda * m_gx * m_gx;
            Real const qy = m_lambda * m_gy * m_gy;
            for (int j=m_jlo; j<=m_jhi+1; ++j) {
                for (int i=m_ilo; i<=m_ihi; ++i) {
                    d_arr(i  ,j,0) += qx;
                    d_arr(i+1,j,0) += qx;
                }
            }
            for (int j=m_jlo; j<=m_jhi; ++j) {
                for (int i=m_ilo; i<=m_ihi+1; ++i) {
                    d_arr(i,j  ,0) += qy;
                    d_arr(i,j+1,0) += qy;
                }
            }
        } else {
            Real const qx = m_lambda * m_gxx * m_gxx;
            Real const qy = m_lambda * m_gyy * m_gyy;
            Real const qc = m_lambda * four * (m_gxx + m_gyy) * (m_gxx + m_gyy);
            for (int j=m_jlo+1; j<=m_jhi; ++j) {
                for (int i=m_ilo+1; i<=m_ihi; ++i) {
                    d_arr(i-1,j  ,0) += qx;
                    d_arr(i+1,j  ,0) += qx;
                    d_arr(i  ,j-1,0) += qy;
                    d_arr(i  ,j+1,0) += qy;
                    d_arr(i  ,j  ,0) += qc;
                }
            }
        }
    }

    /**
     * @brief Compute diagnostic metrics for the reconstructed nodal field.
     * @param[in] T Cell-centered data.
     * @param[in] S Reconstructed nodal field.
     * @param[in,out] info Structure to store diagnostics.
     */
    void computeDiagnostics (amrex::FArrayBox const& T,
                             amrex::FArrayBox const& S,
                             SolveInfo& info) const
    {
        auto const S_arr = S.const_array();
        Real smin =  std::numeric_limits<Real>::max();
        Real smax = -std::numeric_limits<Real>::max();
        for (int j=m_jlo; j<=m_jhi+1; ++j) {
            for (int i=m_ilo; i<=m_ihi+1; ++i) {
                smin = std::min(smin, S_arr(i,j,0));
                smax = std::max(smax, S_arr(i,j,0));
            }
        }
        info.min_value = smin;
        info.max_value = smax;
        info.max_average_error = maxAverageError(T,S);
    }

    [[nodiscard]] Real maxAbs (amrex::FArrayBox const& F) const
    {
        auto const F_arr = F.const_array();
        Real v = zero;
        for (int j=m_jlo; j<=m_jhi+1; ++j) {
            for (int i=m_ilo; i<=m_ihi+1; ++i) {
                v = std::max(v, std::abs(F_arr(i,j,0)));
            }
        }
        return v;
    }

    //! Jacobi-preconditioned conjugate gradient on the SPD operator M.
    //! x is both the initial guess and the result; it is reset to zero here
    //! because every attempt restarts the correction from the interpolant.
    [[nodiscard]] SolveInfo conjugateGradient (amrex::FArrayBox const& rhs,
                                               amrex::FArrayBox& x,
                                               Real relative_tolerance,
                                               int max_iterations) const
    {
        /**
         * @brief Local storage for the operator diagonal.
         */
        amrex::FArrayBox diag(m_nodal_box,1,amrex::The_Managed_Arena());
        computeDiagonal(diag);

        /**
         * @brief Local storage for the residual.
         */
        amrex::FArrayBox r (m_nodal_box,1,amrex::The_Managed_Arena());
        /**
         * @brief Local storage for the preconditioned residual.
         */
        amrex::FArrayBox z (m_nodal_box,1,amrex::The_Managed_Arena());
        /**
         * @brief Local storage for the search direction.
         */
        amrex::FArrayBox p (m_nodal_box,1,amrex::The_Managed_Arena());
        /**
         * @brief Local storage for the product of the operator and the search direction.
         */
        amrex::FArrayBox Mp(m_nodal_box,1,amrex::The_Managed_Arena());

        x.setVal<amrex::RunOn::Host>(zero);
        r.copy<amrex::RunOn::Host>(rhs,0,0,1);   // r = rhs - M*0
        applyJacobi(diag,r,z);
        p.copy<amrex::RunOn::Host>(z,0,0,1);

        Real rz = dot(r,z);
        Real const initial = std::sqrt(dot(r,r));
        Real const tol = relative_tolerance * std::max({initial,one});

        SolveInfo info;
        info.final_residual = initial;
        if (initial <= tol) {
            info.converged = true;
            return info;
        }

        for (int it=0; it<max_iterations; ++it) {
            applyM(p,Mp);
            Real const pMp = dot(p,Mp);
            if (!(pMp > zero) || !std::isfinite(pMp)) {
                throw std::runtime_error("CG curvature failure.");
            }

            Real const alpha = rz/pMp;
            {
                auto const x_arr  = x.array();
                auto const r_arr  = r.array();
                auto const p_arr  = p.const_array();
                auto const Mp_arr = Mp.const_array();
                for (int j=m_jlo; j<=m_jhi+1; ++j) {
                    for (int i=m_ilo; i<=m_ihi+1; ++i) {
                        x_arr(i,j,0) += alpha * p_arr(i,j,0);
                        r_arr(i,j,0) -= alpha * Mp_arr(i,j,0);
                    }
                }
            }

            info.iterations = it+1;
            info.final_residual = std::sqrt(dot(r,r));
            if (info.final_residual <= tol) {
                info.converged = true;
                return info;
            }

            applyJacobi(diag,r,z);
            Real const rz_new = dot(r,z);
            Real const beta = rz_new/rz;
            {
                auto const p_arr = p.array();
                auto const z_arr = z.const_array();
                for (int j=m_jlo; j<=m_jhi+1; ++j) {
                    for (int i=m_ilo; i<=m_ihi+1; ++i) {
                        p_arr(i,j,0) = z_arr(i,j,0) + beta * p_arr(i,j,0);
                    }
                }
            }
            rz = rz_new;
        }
        return info;
    }

    /**
     * @brief Apply the Jacobi preconditioner.
     * @param[in] diag Operator diagonal.
     * @param[in] r Residual field.
     * @param[out] z Preconditioned residual field.
     */
    void applyJacobi (amrex::FArrayBox const& diag,
                      amrex::FArrayBox const& r,
                      amrex::FArrayBox& z) const
    {
        auto const d_arr = diag.const_array();
        auto const r_arr = r.const_array();
        auto const z_arr = z.array();
        for (int j=m_jlo; j<=m_jhi+1; ++j) {
            for (int i=m_ilo; i<=m_ihi+1; ++i) {
                z_arr(i,j,0) = r_arr(i,j,0) / d_arr(i,j,0);
            }
        }
    }

    [[nodiscard]] Real dot (amrex::FArrayBox const& a,
                            amrex::FArrayBox const& b) const
    {
        auto const a_arr = a.const_array();
        auto const b_arr = b.const_array();
        Real v = zero;
        for (int j=m_jlo; j<=m_jhi+1; ++j) {
            for (int i=m_ilo; i<=m_ihi+1; ++i) {
                v += a_arr(i,j,0) * b_arr(i,j,0);
            }
        }
        return v;
    }
};

#endif
