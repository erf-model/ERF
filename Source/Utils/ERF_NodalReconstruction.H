#ifndef ERF_NODAL_RECONSTRUCTION_H_
#define ERF_NODAL_RECONSTRUCTION_H_

#include <algorithm>
#include <cmath>
#include <cstddef>
#include <limits>
#include <stdexcept>
#include <utility>

#include "AMReX_Array4.H"
#include "AMReX_Box.H"
#include "AMReX_FArrayBox.H"
#include "AMReX_REAL.H"
#include "AMReX_Vector.H"

#include "ERF_Constants.H"

enum class VariationOperator
{
    FirstDeriv,
    Laplacian
};

struct SolveInfo {
    int iterations = 0;
    amrex::Real final_residual = zero;
    bool converged = false;
};

class NodalReconstruction
{
public:
    using Real   = amrex::Real;
    using Vector = amrex::Vector<Real>;

    explicit NodalReconstruction (amrex::Box const& face_box,
                                  amrex::Geometry const& geom)
        : m_face_box(face_box),
          m_ilo(face_box.smallEnd(0)),
          m_ihi(face_box.bigEnd(0)),
          m_jlo(face_box.smallEnd(1)),
          m_jhi(face_box.bigEnd(1)),
          m_nx(face_box.length(0)),
          m_ny(face_box.length(1)),
          m_boundary_count(m_nx + m_ny + 1)
    {
        if (face_box.length(2) != 1) {
            throw std::invalid_argument(
                "face_box must contain exactly one z cell.");
        }

        m_nodal_box = amrex::surroundingNodes(face_box);
        m_nodal_box.setSmall(2, 0);
        m_nodal_box.setBig(2, 0);

        m_dx = geom.CellSizeArray()[0];
        m_dy = geom.CellSizeArray()[1];
    }

    [[nodiscard]] int boundaryCount () const noexcept {
        return m_boundary_count;
    }

    [[nodiscard]] amrex::Box const& nodalBox () const noexcept {
        return m_nodal_box;
    }

    [[nodiscard]] Vector makeBoundaryFromT (amrex::FArrayBox const& T) const
    {
        auto const T_arr = T.const_array();
        Vector b(static_cast<std::size_t>(m_boundary_count), zero);

        b[bottomIndex(0)]         = T_arr(m_ilo, m_jlo, 0);
        for (int ii = 1; ii < m_nx; ++ii) {
            int const i = m_ilo + ii;
            b[bottomIndex(ii)] = myhalf * ( T_arr(i-1, m_jlo, 0) + T_arr(i,  m_jlo, 0) );
        }
        b[bottomIndex(m_nx)] = T_arr(m_ihi, m_jlo, 0);

        for (int jj = 1; jj < m_ny; ++jj) {
            int const j = m_jlo + jj;
            b[leftIndex(jj)] = myhalf * ( T_arr(m_ilo, j-1, 0) + T_arr(m_ilo, j, 0) );
        }
        b[leftIndex(m_ny)] = T_arr(m_ilo, m_jhi, 0);
        return b;
    }

    [[nodiscard]] std::pair<amrex::FArrayBox,SolveInfo> solve (amrex::FArrayBox const& T,
                                                               Vector const& boundary_reference,
                                                               VariationOperator var_op,
                                                               Real epsilon = Real(1.e-3),
                                                               Real relative_tolerance = Real(1.e-8),
                                                               int max_iterations = 100)
    {
        m_var_op = var_op;

        amrex::FArrayBox c(m_nodal_box,1);
        computeParticularSolution(T,c);

        amrex::FArrayBox vtv_c(m_nodal_box,1);
        applyVariationOperator(c,vtv_c);

        Vector rhs;
        applyGT(vtv_c,rhs);
        for (int m=0; m<m_boundary_count; ++m) {
            auto const k = static_cast<std::size_t>(m);
            rhs[k] = -rhs[k] + epsilon*boundary_reference[k];
        }

        Vector b = boundary_reference;
        SolveInfo info = conjugateGradient(rhs, b, epsilon, relative_tolerance, max_iterations);

        amrex::FArrayBox S(m_nodal_box,1);
        reconstruct(T,b,S);
        return {std::move(S),info};
    }

    void reconstruct (amrex::FArrayBox const& T,
                      Vector const& boundary,
                      amrex::FArrayBox& S) const
    {
        S.setVal(zero);
        auto const T_arr = T.const_array();
        auto const S_arr = S.array();
        injectBoundary(boundary,S_arr);

        for (int j=m_jlo; j<=m_jhi; ++j) {
            for (int i=m_ilo; i<=m_ihi; ++i) {
                S_arr(i+1,j+1,0) = four * T_arr(i,j,0)
                                 - S_arr(i,  j  , 0)
                                 - S_arr(i+1,j  , 0)
                                 - S_arr(i,  j+1, 0);
            }
        }
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
    amrex::Box m_face_box;
    amrex::Box m_nodal_box;
    int m_ilo=0;
    int m_ihi=-1;
    int m_jlo=0;
    int m_jhi=-1;
    int m_nx=0;
    int m_ny=0;
    int m_boundary_count=0;
    amrex::Real m_dx;
    amrex::Real m_dy;
    VariationOperator m_var_op;

    [[nodiscard]] std::size_t bottomIndex (int ii) const noexcept
    {
        return static_cast<std::size_t>(ii);
    }

    [[nodiscard]] std::size_t leftIndex (int jj) const noexcept
    {
        return static_cast<std::size_t>(m_nx+jj);
    }

    void injectBoundary (Vector const& b,
                         amrex::Array4<Real> const& u) const
    {
        for (int ii=0; ii<=m_nx; ++ii) {
            u(m_ilo+ii, m_jlo , 0) = b[bottomIndex(ii)];
        }
        for (int jj=1; jj<=m_ny; ++jj) {
            u(m_ilo, m_jlo+jj, 0) = b[leftIndex(jj)];
        }
    }

    void extractBoundary (amrex::Array4<Real const> const& u,
                          Vector& b) const
    {
        b.assign(static_cast<std::size_t>(m_boundary_count), zero);
        for (int ii=0; ii<=m_nx; ++ii) {
            b[bottomIndex(ii)]=u(m_ilo+ii, m_jlo, 0);
        }
        for (int jj=1; jj<=m_ny; ++jj) {
            b[leftIndex(jj)]=u(m_ilo, m_jlo+jj, 0);
        }
    }

    void computeParticularSolution (amrex::FArrayBox const& T,
                                    amrex::FArrayBox& c) const
    {
        c.setVal(zero);
        auto const T_arr=T.const_array();
        auto const c_arr=c.array();
        for (int j=m_jlo; j<=m_jhi; ++j) {
            for (int i=m_ilo; i<=m_ihi; ++i) {
                c_arr(i+1,j+1,0) = four * T_arr(i,j,0)
                                 - c_arr(i,  j  , 0)
                                 - c_arr(i+1,j  , 0)
                                 - c_arr(i,  j+1, 0);
            }
        }
    }

    void applyG (Vector const& b,
                 amrex::FArrayBox& U) const
    {
        U.setVal(zero);
        auto const U_arr = U.array();
        injectBoundary(b,U_arr);
        for (int j=m_jlo; j<=m_jhi; ++j) {
            for (int i=m_ilo; i<=m_ihi; ++i) {
                U_arr(i+1,j+1,0) = -U_arr(i  ,j  , 0)
                                   -U_arr(i+1,j  , 0)
                                   -U_arr(i  ,j+1, 0);
            }
        }
    }

    void applyD (amrex::FArrayBox const& S,
                 amrex::FArrayBox& difx,
                 amrex::FArrayBox& dify) const
    {
        difx.setVal(zero);
        dify.setVal(zero);

        auto const S_arr    = S.const_array();
        auto const difx_arr = difx.array();
        auto const dify_arr = dify.array();

        const amrex::Real idx = one / m_dx;
        const amrex::Real idy = one / m_dy;

        for (int j = m_jlo; j <= m_jhi + 1; ++j) {
            for (int i = m_ilo; i <= m_ihi; ++i) {
                difx_arr(i, j, 0) = idx * ( S_arr(i+1, j, 0) - S_arr(i, j, 0) );
            }
        }

        for (int j = m_jlo; j <= m_jhi; ++j) {
            for (int i = m_ilo; i <= m_ihi + 1; ++i) {
                dify_arr(i, j, 0) = idy * ( S_arr(i, j+1, 0) - S_arr(i, j, 0) );
            }
        }
    }

    void applyDT (amrex::FArrayBox const& difx,
                  amrex::FArrayBox const& dify,
                  amrex::FArrayBox& R) const
    {
        R.setVal(zero);

        auto const qx_arr = difx.const_array();
        auto const qy_arr = dify.const_array();
        auto const R_arr  = R.array();

        const amrex::Real idx = one / m_dx;
        const amrex::Real idy = one / m_dy;

        for (int j = m_jlo; j <= m_jhi + 1; ++j) {
            for (int i = m_ilo; i <= m_ihi; ++i) {
                R_arr(i  , j, 0) -= idx * qx_arr(i, j, 0);
                R_arr(i+1, j, 0) += idx * qx_arr(i, j, 0);
            }
        }

        for (int j = m_jlo; j <= m_jhi; ++j) {
            for (int i = m_ilo; i <= m_ihi + 1; ++i) {
                R_arr(i, j  , 0) -= idy * qy_arr(i, j, 0);
                R_arr(i, j+1, 0) += idy * qy_arr(i, j, 0);
            }
        }

    }

    void applyDTD (amrex::FArrayBox const& S,
                   amrex::FArrayBox& R) const
    {
        amrex::Box x_edge_box = m_nodal_box;
        amrex::Box y_edge_box = m_nodal_box;
        x_edge_box.growHi(0,-1);
        y_edge_box.growHi(1,-1);
        amrex::FArrayBox difx(x_edge_box, 1);
        amrex::FArrayBox dify(y_edge_box, 1);
        applyD(S, difx, dify);
        applyDT(difx, dify, R);
    }

    void applyL (amrex::FArrayBox const& S,
                 amrex::FArrayBox& lap) const
    {
        lap.setVal(zero);

        auto const S_arr   = S.const_array();
        auto const lap_arr = lap.array();

        const amrex::Real idx2 = one / (m_dx * m_dx);
        const amrex::Real idy2 = one / (m_dy * m_dy);

        for (int j = m_jlo; j <= m_jhi; ++j) {
            for (int i = m_ilo; i <= m_ihi; ++i) {
                amrex::Real im_val = (i==m_ilo) ? S_arr(i, j, 0) : S_arr(i-1, j  , 0);
                amrex::Real ip_val = (i==m_ihi) ? S_arr(i, j, 0) : S_arr(i+1, j  , 0);
                amrex::Real jm_val = (j==m_jlo) ? S_arr(i, j, 0) : S_arr(i  , j-1, 0);
                amrex::Real jp_val = (j==m_jhi) ? S_arr(i, j, 0) : S_arr(i  , j+1, 0);
                lap_arr(i, j, 0) = idx2 * ( ip_val - two * S_arr(i, j, 0) + im_val )
                                 + idy2 * ( jp_val - two * S_arr(i, j, 0) + jm_val );
            }
        }
    }

    void applyLT (amrex::FArrayBox const& lap,
                  amrex::FArrayBox& R) const
    {
        R.setVal(zero);

        auto const q_arr = lap.const_array();
        auto const R_arr = R.array();

        const amrex::Real idx2 = one / (m_dx * m_dx);
        const amrex::Real idy2 = one / (m_dy * m_dy);

        for (int j = m_jlo + 1; j <= m_jhi; ++j) {
            for (int i = m_ilo + 1; i <= m_ihi; ++i) {
                const amrex::Real val = q_arr(i, j, 0);
                R_arr(i-1, j  , 0) += idx2 * val;
                R_arr(i+1, j  , 0) += idx2 * val;
                R_arr(i  , j-1, 0) += idy2 * val;
                R_arr(i  , j+1, 0) += idy2 * val;
                R_arr(i  , j  , 0) -= two * (idx2 + idy2) * val;
            }
        }
    }

    void applyLTL (amrex::FArrayBox const& S,
                   amrex::FArrayBox& R) const
    {
        amrex::FArrayBox lap(m_nodal_box, 1);
        applyL(S, lap);
        applyLT(lap, R);
    }

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

    void applyGT (amrex::FArrayBox const& R,
                  Vector& result) const
    {
        amrex::FArrayBox adj(m_nodal_box,1);
        adj.copy(R,0,0,1);
        auto const a_arr = adj.array();

        for (int j=m_jhi; j>=m_jlo; --j) {
            for (int i=m_ihi; i>=m_ilo; --i) {
                Real const child= a_arr(i+1,j+1,0);
                a_arr(i,  j  ,0) -= child;
                a_arr(i+1,j  ,0) -= child;
                a_arr(i,  j+1,0) -= child;
            }
        }
        extractBoundary(adj.const_array(),result);
    }

    void applyH (Vector const& v,
                 Vector& Hv,
                 Real epsilon) const
    {
        amrex::FArrayBox Gv(m_nodal_box,1);
        amrex::FArrayBox q(m_nodal_box,1);
        applyG(v,Gv);
        applyVariationOperator(Gv,q);
        applyGT(q,Hv);
        for (int m=0; m<m_boundary_count; ++m) {
            auto const k=static_cast<std::size_t>(m);
            Hv[k]+=epsilon*v[k];
        }
    }

    [[nodiscard]] SolveInfo conjugateGradient (Vector const& rhs,
                                               Vector& x,
                                               Real epsilon,
                                               Real relative_tolerance,
                                               int max_iterations) const
    {
        Vector Hx;
        applyH(x,Hx,epsilon);

        Vector r(rhs.size()),p(rhs.size()),Hp;
        for (int k=0; k<rhs.size(); ++k) {
            r[k]=rhs[k]-Hx[k];
            p[k]=r[k];
        }

        Real rr=dot(r,r);
        Real const initial=std::sqrt(rr);
        Real const tol=relative_tolerance*std::max({norm(rhs),initial,one});

        SolveInfo info;
        info.final_residual=initial;
        if (initial<=tol) {
            info.converged=true;
            return info;
        }

        for (int it=0; it<max_iterations; ++it) {
            applyH(p,Hp,epsilon);
            Real const pHp=dot(p,Hp);
            if (!(pHp>zero) || !std::isfinite(pHp)) {
                throw std::runtime_error("CG curvature failure.");
            }

            Real const alpha=rr/pHp;
            for (int k=0; k<x.size(); ++k) {
                x[k]+=alpha*p[k];
                r[k]-=alpha*Hp[k];
            }

            Real const rr_new=dot(r,r);
            info.iterations=it+1;
            info.final_residual=std::sqrt(rr_new);
            if (info.final_residual<=tol) {
                info.converged=true;
                return info;
            }

            Real const beta=rr_new/rr;
            for (int k=0; k<p.size(); ++k) {
                p[k]=r[k]+beta*p[k];
            }
            rr=rr_new;
        }
        return info;
    }

    [[nodiscard]] static Real dot (Vector const& a,
                                   Vector const& b)
    {
        if (a.size()!=b.size()) {
            throw std::invalid_argument("Vector size mismatch.");
        }
        Real v = zero;
        for (int k=0; k<a.size(); ++k) v+=a[k]*b[k];
        return v;
    }

    [[nodiscard]] static Real norm (Vector const& a) {
        return std::sqrt(dot(a,a));
    }
};

#endif
