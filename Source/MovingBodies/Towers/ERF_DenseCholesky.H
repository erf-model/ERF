#ifndef ERF_DENSE_CHOLESKY_H
#define ERF_DENSE_CHOLESKY_H

/**
 * A dense Cholesky factorisation A = L L^T of a symmetric positive-definite matrix, in double
 * precision, for the frame model's stiffness matrix: a lattice tower has a few hundred to a few
 * thousand degrees of freedom, factored once and solved many times on every rank alike.
 */

#include <cmath>
#include <cstddef>
#include <vector>

namespace erf_towers {

/** One factorisation A = L L^T, L lower triangular, stored row by row in packed form. */
class DenseCholesky
{
public:
    DenseCholesky () = default;

    /**
     * Factor the n x n symmetric matrix a (row-major; only the lower triangle a[i*n+j], j <= i, is
     * read). Returns -1 on success, or the 0-based row whose pivot is not positive: that row's
     * unknown is not held by the rows before it (a mechanism, for a stiffness matrix). A pivot
     * counts as not positive when it falls to relative_tolerance times its diagonal entry or below.
     */
    long factor (const std::vector<double>& a, std::size_t n, double relative_tolerance = 1.0e-12)
    {
        m_n = n;
        m_l.assign(n * (n + 1) / 2, 0.0);
        for (std::size_t i = 0; i < n; ++i) {
            for (std::size_t j = 0; j <= i; ++j) {
                double s = a[i * n + j];
                const double* li = &m_l[row(i)];
                const double* lj = &m_l[row(j)];
                for (std::size_t k = 0; k < j; ++k) { s -= li[k] * lj[k]; }
                if (i == j) {
                    const double d = a[i * n + i];
                    if (!(s > relative_tolerance * std::abs(d)) || !(d > 0.0)) { m_n = 0; m_l.clear(); return static_cast<long>(i); }
                    m_l[row(i) + i] = std::sqrt(s);
                } else {
                    m_l[row(i) + j] = s / m_l[row(j) + j];
                }
            }
        }
        return -1;
    }

    /** Whether factor() succeeded. */
    bool factored () const { return m_n > 0; }
    /** The order n of the factored matrix. */
    std::size_t size () const { return m_n; }

    /** Solve A x = b in place (b holds n values on entry, x on return). Precondition: factored(). */
    void solve (std::vector<double>& b) const
    {
        const std::size_t n = m_n;
        for (std::size_t i = 0; i < n; ++i) {            // L y = b
            const double* li = &m_l[row(i)];
            double s = b[i];
            for (std::size_t k = 0; k < i; ++k) { s -= li[k] * b[k]; }
            b[i] = s / li[i];
        }
        for (std::size_t i = n; i-- > 0;) {               // L^T x = y
            double s = b[i];
            for (std::size_t k = i + 1; k < n; ++k) { s -= m_l[row(k) + i] * b[k]; }
            b[i] = s / m_l[row(i) + i];
        }
    }

private:
    static std::size_t row (std::size_t i) { return i * (i + 1) / 2; }
    std::size_t m_n{0};
    std::vector<double> m_l;   //!< L, row i at [i(i+1)/2, i(i+1)/2 + i]
};

/**
 * The number of negative pivots of the LDL^T factorisation of the n x n symmetric matrix a
 * (row-major, lower triangle read), without pivoting: by Sylvester's law of inertia, the number of
 * negative eigenvalues of a. For a = K - sigma M with K and M positive definite, it is the number of
 * eigenvalues of K phi = lambda M phi below sigma (the Sturm sequence count). A zero pivot counts as
 * negative; sigma should not sit on an eigenvalue.
 */
inline std::size_t negative_pivots (std::vector<double> a, std::size_t n)
{
    std::size_t count = 0;
    for (std::size_t k = 0; k < n; ++k) {
        const double d = a[k * n + k];
        if (!(d > 0.0)) { ++count; }
        if (d == 0.0) { continue; }
        for (std::size_t i = k + 1; i < n; ++i) {
            const double l = a[i * n + k] / d;
            if (l == 0.0) { continue; }
            for (std::size_t j = k + 1; j <= i; ++j) { a[i * n + j] -= l * a[j * n + k]; }
        }
    }
    return count;
}

} // namespace erf_towers

#endif
