// The frame's natural modes (subspace iteration with a Sturm check) and its Newmark time response.

#include "ERF_FrameDynamics.H"

#include <algorithm>
#include <cmath>
#include <limits>
#include <numeric>
#include <iterator>
#include <tuple>
#include <memory>
#include <map>

#include <AMReX.H>
#include <AMReX_BLassert.H>

namespace erf_towers {

namespace {

/**
 * K_eff = K (1 + 2 a1/h) + M (4/h^2 + 2 a0/h) of frame f factored, shared by every FrameDynamics of the
 * same frame, damping and step: the towers of one type step with one factorisation, not one each. The
 * cache holds weak pointers, so a factorisation lives as long as a FrameDynamics uses it. Not thread-safe:
 * the towers are stepped one after another on each rank.
 */
std::shared_ptr<const DenseCholesky> effective_stiffness (const Frame& f, double a0, double a1, double h)
{
    AMREX_ALWAYS_ASSERT_WITH_MESSAGE(std::isfinite(a0) && std::isfinite(a1) && std::isfinite(h),
                                     "FrameDynamics: the Rayleigh coefficients and the step must be finite");
    using Key = std::tuple<const Frame*, double, double, double>;
    static std::map<Key, std::weak_ptr<const DenseCholesky>> cache;
    const Key key{&f, a0, a1, h};
    const auto found = cache.find(key);
    if (found != cache.end()) {
        if (auto k = found->second.lock()) { return k; }
    }
    for (auto it = cache.begin(); it != cache.end();) { it = it->second.expired() ? cache.erase(it) : std::next(it); }
    const double c0 = 4.0 / (h * h), c1 = 2.0 / h;
    auto k = std::make_shared<DenseCholesky>();
    const long bad = k->factor(f.assemble_free(1.0 + c1 * a1, c0 + c1 * a0), f.num_free_dofs());
    AMREX_ALWAYS_ASSERT_WITH_MESSAGE(bad < 0, "FrameDynamics::step: the effective stiffness is not positive definite");
    cache[key] = k;
    return k;
}

/**
 * The free mass matrix of frame f factored, shared by every FrameDynamics of the same frame like
 * effective_stiffness(). Not thread-safe, for the same reason. A singular mass matrix gives null,
 * with bad the free degree of freedom where the factorisation failed.
 */
std::shared_ptr<const DenseCholesky> mass_factor (const Frame& f, long& bad)
{
    static std::map<const Frame*, std::weak_ptr<const DenseCholesky>> cache;
    bad = -1;
    const auto found = cache.find(&f);
    if (found != cache.end()) {
        if (auto m = found->second.lock()) { return m; }
    }
    for (auto it = cache.begin(); it != cache.end();) { it = it->second.expired() ? cache.erase(it) : std::next(it); }
    auto m = std::make_shared<DenseCholesky>();
    bad = m->factor(f.assemble_free(0.0, 1.0), f.num_free_dofs());
    if (bad >= 0) { return nullptr; }
    cache[&f] = m;
    return m;
}

constexpr double two_pi = 2.0 * 3.14159265358979323846;

/**
 * The eigenvalues w and eigenvectors (columns of v, row-major n x n) of the symmetric n x n matrix a
 * (row-major) by cyclic Jacobi rotations, to an off-diagonal norm of 1e-15 of the matrix norm.
 */
void jacobi (std::vector<double> a, std::size_t n, std::vector<double>& w, std::vector<double>& v)
{
    v.assign(n * n, 0.0);
    for (std::size_t i = 0; i < n; ++i) { v[i * n + i] = 1.0; }
    double norm = 0.0;
    for (const double x : a) { norm += x * x; }
    norm = std::sqrt(norm);
    for (int sweep = 0; sweep < 100; ++sweep) {
        double off = 0.0;
        for (std::size_t i = 0; i < n; ++i) { for (std::size_t j = i + 1; j < n; ++j) { off += a[i * n + j] * a[i * n + j]; } }
        if (std::sqrt(2.0 * off) <= 1.0e-15 * norm) { break; }
        for (std::size_t p = 0; p < n; ++p) {
            for (std::size_t q = p + 1; q < n; ++q) {
                const double apq = a[p * n + q];
                if (apq == 0.0) { continue; }
                const double theta = (a[q * n + q] - a[p * n + p]) / (2.0 * apq);
                const double t = (theta >= 0.0 ? 1.0 : -1.0) / (std::abs(theta) + std::sqrt(theta * theta + 1.0));
                const double c = 1.0 / std::sqrt(t * t + 1.0), s = t * c;
                for (std::size_t k = 0; k < n; ++k) {                 // A <- A J
                    const double akp = a[k * n + p], akq = a[k * n + q];
                    a[k * n + p] = c * akp - s * akq;
                    a[k * n + q] = s * akp + c * akq;
                }
                for (std::size_t k = 0; k < n; ++k) {                 // A <- J^T A
                    const double apk = a[p * n + k], aqk = a[q * n + k];
                    a[p * n + k] = c * apk - s * aqk;
                    a[q * n + k] = s * apk + c * aqk;
                }
                for (std::size_t k = 0; k < n; ++k) {                 // V <- V J
                    const double vkp = v[k * n + p], vkq = v[k * n + q];
                    v[k * n + p] = c * vkp - s * vkq;
                    v[k * n + q] = s * vkp + c * vkq;
                }
            }
        }
    }
    w.resize(n);
    for (std::size_t i = 0; i < n; ++i) { w[i] = a[i * n + i]; }
}

double dot (const std::vector<double>& x, const std::vector<double>& y)
{
    double s = 0.0;
    for (std::size_t i = 0; i < x.size(); ++i) { s += x[i] * y[i]; }
    return s;
}

} // namespace

std::string frame_modes (const Frame& f, std::size_t count, FrameModes& modes)
{
    const std::size_t nf = f.num_free_dofs();
    const std::vector<std::size_t>& free = f.free_dofs();
    if (count == 0 || count > nf) { return "frame_modes: ask for 1 to " + std::to_string(nf) + " modes"; }
    const std::size_t q = std::min(nf, std::max(2 * count, count + 8));
    // free <-> full vectors and the mass operator on free vectors
    auto full = [&] (const std::vector<double>& x) {
        std::vector<double> y(f.num_dofs(), 0.0);
        for (std::size_t i = 0; i < nf; ++i) { y[free[i]] = x[i]; }
        return y;
    };
    auto mass_times = [&] (const std::vector<double>& x) {
        const std::vector<double> y = f.apply_mass(full(x));
        std::vector<double> r(nf);
        for (std::size_t i = 0; i < nf; ++i) { r[i] = y[free[i]]; }
        return r;
    };
    // starting vectors (Bathe): the mass diagonal, unit vectors where M_ii/K_ii is largest, and a spread vector
    std::vector<double> kd(nf), md(nf);
    {
        const std::vector<double> k = f.assemble_free(1.0, 0.0), m = f.assemble_free(0.0, 1.0);
        for (std::size_t i = 0; i < nf; ++i) { kd[i] = k[i * nf + i]; md[i] = m[i * nf + i]; }
    }
    std::vector<std::size_t> order(nf);
    std::iota(order.begin(), order.end(), std::size_t(0));
    std::stable_sort(order.begin(), order.end(), [&] (std::size_t a, std::size_t b) { return md[a] / kd[a] > md[b] / kd[b]; });
    std::vector<std::vector<double>> x(q, std::vector<double>(nf, 0.0));
    x[0] = md;
    for (std::size_t j = 1; j + 1 < q; ++j) { x[j][order[j - 1]] = 1.0; }
    if (q > 1) { for (std::size_t i = 0; i < nf; ++i) { x[q - 1][i] = std::sin(1.3 * static_cast<double>(i) + 0.7); } }
    if (q == nf) {   // the whole space: unit vectors
        for (std::size_t j = 0; j < q; ++j) { x[j].assign(nf, 0.0); x[j][j] = 1.0; }
    }
    // the eigenvalues that must converge: those asked for and the next one, to place the Sturm shift
    // between them; more when the next ones are too close to separate from the last one asked for
    std::size_t pc = std::min(q, count + 1);
    std::vector<double> lambda(q, 0.0), previous(q, 0.0);
    bool converged = false;
    for (int it = 0; it < 500 && !converged; ++it) {
        std::vector<std::vector<double>> mx(q), y(q), my(q);
        for (std::size_t j = 0; j < q; ++j) {
            mx[j] = mass_times(x[j]);
            y[j] = mx[j];
            f.solve_free(y[j]);
            my[j] = mass_times(y[j]);
        }
        std::vector<double> kr(q * q), mr(q * q);
        for (std::size_t i = 0; i < q; ++i) {
            for (std::size_t j = i; j < q; ++j) {
                kr[i * q + j] = kr[j * q + i] = 0.5 * (dot(y[i], mx[j]) + dot(y[j], mx[i]));
                mr[i * q + j] = mr[j * q + i] = 0.5 * (dot(y[i], my[j]) + dot(y[j], my[i]));
            }
        }
        // Kr Q = Mr Q Lambda through Kr, positive definite: Kr = W D W^T, T = W D^(-1/2), then
        // T^T Mr T V = V mu with mu = 1/lambda; a degree of freedom without mass gives mu = 0
        std::vector<double> dk, w;
        jacobi(kr, q, dk, w);
        for (std::size_t j = 0; j < q; ++j) {
            if (!(dk[j] > 0.0)) { return "frame_modes: the stiffness of the trial vectors is not positive (they lost their rank)"; }
        }
        std::vector<double> t(q * q), mt(q * q, 0.0), b(q * q, 0.0);
        for (std::size_t i = 0; i < q; ++i) { for (std::size_t j = 0; j < q; ++j) { t[i * q + j] = w[i * q + j] / std::sqrt(dk[j]); } }
        for (std::size_t i = 0; i < q; ++i) {           // Mr T
            for (std::size_t l = 0; l < q; ++l) {
                const double m = mr[i * q + l];
                for (std::size_t j = 0; j < q; ++j) { mt[i * q + j] += m * t[l * q + j]; }
            }
        }
        for (std::size_t k = 0; k < q; ++k) {           // T^T (Mr T)
            for (std::size_t i = 0; i < q; ++i) {
                const double tki = t[k * q + i];
                for (std::size_t j = 0; j < q; ++j) { b[i * q + j] += tki * mt[k * q + j]; }
            }
        }
        std::vector<double> mu, v;
        jacobi(b, q, mu, v);
        std::vector<std::size_t> idx(q);
        std::iota(idx.begin(), idx.end(), std::size_t(0));
        std::sort(idx.begin(), idx.end(), [&] (std::size_t a, std::size_t c) { return mu[a] > mu[c]; });
        double mu_max = 0.0;
        for (const double m : mu) { mu_max = std::max(mu_max, m); }
        // the new trial vectors X = Y T V, each scaled to unit mass (V^T B V = diag(mu))
        std::vector<std::vector<double>> xn(q, std::vector<double>(nf, 0.0));
        for (std::size_t jj = 0; jj < q; ++jj) {
            const std::size_t j = idx[jj];
            const bool massless = !(mu[j] > 1.0e-14 * mu_max);
            if (massless && jj < pc) {
                return "frame_modes: mode " + std::to_string(jj + 1) + " has no mass (a degree of freedom without mass or inertia?)";
            }
            lambda[jj] = massless ? std::numeric_limits<double>::max() : 1.0 / mu[j];
            const double scale = massless ? 1.0 : 1.0 / std::sqrt(mu[j]);
            for (std::size_t k = 0; k < q; ++k) {
                double c = 0.0;
                for (std::size_t l = 0; l < q; ++l) { c += t[k * q + l] * v[l * q + j]; }
                c *= scale;
                if (c == 0.0) { continue; }
                for (std::size_t i = 0; i < nf; ++i) { xn[jj][i] += c * y[k][i]; }
            }
        }
        x.swap(xn);
        converged = (it > 0);
        for (std::size_t j = 0; j < pc; ++j) {
            if (!(std::abs(lambda[j] - previous[j]) <= 1.0e-12 * std::abs(lambda[j]))) { converged = false; }
        }
        previous = lambda;
        if (converged && q < nf && !(lambda[pc - 1] > lambda[pc - 2] * (1.0 + 1.0e-8))) {
            if (pc + 1 >= q) { return "frame_modes: the modes from " + std::to_string(count) + " on are too close to place a Sturm check"; }
            ++pc;
            converged = false;
        }
    }
    if (!converged) { return "frame_modes: the subspace iteration did not converge in 500 iterations"; }
    // Sturm check: the number of eigenvalues below a shift between two converged, separated ones
    if (q < nf) {
        const std::size_t boundary = pc - 1;
        const double sigma = 0.5 * (lambda[boundary - 1] + lambda[boundary]);
        const std::size_t below = negative_pivots(f.assemble_free(1.0, -sigma), nf);
        if (below != boundary) {
            return "frame_modes: the Sturm check counts " + std::to_string(below) + " modes below " +
                   std::to_string(std::sqrt(sigma) / two_pi) + " Hz, the subspace iteration found " + std::to_string(boundary);
        }
    }
    modes.frequency.resize(count);
    modes.shape.resize(count);
    for (std::size_t j = 0; j < count; ++j) {
        modes.frequency[j] = std::sqrt(std::max(lambda[j], 0.0)) / two_pi;
        modes.shape[j] = full(x[j]);
    }
    return std::string();
}

void rayleigh_coefficients (double f1, double zeta1, double f2, double zeta2, double& a0, double& a1)
{
    AMREX_ALWAYS_ASSERT_WITH_MESSAGE(f1 > 0.0 && f2 > 0.0 && f1 != f2,
                                     "rayleigh_coefficients: two distinct positive frequencies are needed");
    const double w1 = two_pi * f1, w2 = two_pi * f2, d = w2 * w2 - w1 * w1;
    a0 = 2.0 * w1 * w2 * (zeta1 * w2 - zeta2 * w1) / d;
    a1 = 2.0 * (zeta2 * w2 - zeta1 * w1) / d;
}

FrameDynamics::FrameDynamics (const Frame& frame, double a0, double a1)
    : m_frame(frame), m_a0(a0), m_a1(a1)
{
    AMREX_ALWAYS_ASSERT_WITH_MESSAGE(std::isfinite(a0) && a0 >= 0.0 && std::isfinite(a1) && a1 >= 0.0,
                                     "FrameDynamics: the Rayleigh coefficients must be finite and >= 0");
    m_u.assign(frame.num_dofs(), 0.0);
    m_v = m_u;
    m_a = m_u;
}

std::vector<double> FrameDynamics::to_free (const std::vector<double>& full) const
{
    const auto& free = m_frame.free_dofs();
    std::vector<double> r(free.size());
    for (std::size_t i = 0; i < free.size(); ++i) { r[i] = full[free[i]]; }
    return r;
}

std::vector<double> FrameDynamics::to_full (const std::vector<double>& x) const
{
    const auto& free = m_frame.free_dofs();
    std::vector<double> r(m_frame.num_dofs(), 0.0);
    for (std::size_t i = 0; i < free.size(); ++i) { r[free[i]] = x[i]; }
    return r;
}

void FrameDynamics::start_static (const std::vector<double>& node_loads, double gravity, double t)
{
    std::vector<double> b = to_free(m_frame.load_vector(node_loads, gravity));
    m_frame.solve_free(b);
    m_u = to_full(b);
    m_v.assign(m_u.size(), 0.0);
    m_a.assign(m_u.size(), 0.0);
    m_t = t;
}

std::string FrameDynamics::start (const std::vector<double>& u, const std::vector<double>& v,
                                  const std::vector<double>& node_loads, double gravity, double t)
{
    const std::size_t n = m_frame.num_dofs();
    if (u.size() != n || v.size() != n) {
        return "FrameDynamics::start: one displacement and one velocity per degree of freedom are needed";
    }
    for (std::size_t i = 0; i < n; ++i) {
        if (!std::isfinite(u[i]) || !std::isfinite(v[i])) { return "FrameDynamics::start: the displacement and velocity must be finite"; }
    }
    // M a = f - C v - K u, C = a0 M + a1 K
    const std::vector<double> f = m_frame.load_vector(node_loads, gravity);
    std::vector<double> mv(n, 0.0);
    std::vector<double> w(n);
    for (std::size_t i = 0; i < n; ++i) { w[i] = u[i] + m_a1 * v[i]; }
    const std::vector<double> kw = m_frame.apply_stiffness(w);
    if (m_a0 != 0.0) { mv = m_frame.apply_mass(v); }
    std::vector<double> r(n);
    for (std::size_t i = 0; i < n; ++i) { r[i] = f[i] - kw[i] - m_a0 * mv[i]; }
    long bad = -1;
    std::shared_ptr<const DenseCholesky> mass = mass_factor(m_frame, bad);
    if (!mass) {
        return "FrameDynamics::start: the mass matrix is singular at free degree of freedom " + std::to_string(bad) +
               " (a node without mass or rotary inertia)";
    }
    m_mass = mass;
    std::vector<double> b = to_free(r);
    mass->solve(b);
    m_u = to_full(to_free(u));
    m_v = to_full(to_free(v));
    m_a = to_full(b);
    m_t = t;
    return std::string();
}

void FrameDynamics::step (double h, const std::vector<double>& node_loads, double gravity)
{
    AMREX_ALWAYS_ASSERT_WITH_MESSAGE(std::isfinite(h) && h > 0.0, "FrameDynamics::step: the step must be finite and positive (s)");
    const double c0 = 4.0 / (h * h), c1 = 2.0 / h, c2 = 4.0 / h;
    if (!m_keff || !(std::abs(h - m_h) <= 1.0e-14 * h)) {
        m_keff = effective_stiffness(m_frame, m_a0, m_a1, h);
        m_h = h;
    }
    const std::size_t n = m_u.size();
    // K_eff u' = f' + M (c0 u + c2 v + a + a0 (c1 u + v)) + a1 K (c1 u + v)
    std::vector<double> wm(n), wk(n);
    for (std::size_t i = 0; i < n; ++i) {
        const double cv = c1 * m_u[i] + m_v[i];
        wm[i] = c0 * m_u[i] + c2 * m_v[i] + m_a[i] + m_a0 * cv;
        wk[i] = m_a1 * cv;
    }
    const std::vector<double> f = m_frame.load_vector(node_loads, gravity);
    const std::vector<double> mw = m_frame.apply_mass(wm);
    const std::vector<double> kw = (m_a1 != 0.0) ? m_frame.apply_stiffness(wk) : std::vector<double>(n, 0.0);
    std::vector<double> r(n);
    for (std::size_t i = 0; i < n; ++i) { r[i] = f[i] + mw[i] + kw[i]; }
    std::vector<double> b = to_free(r);
    m_keff->solve(b);
    const std::vector<double> un = to_full(b);
    for (std::size_t i = 0; i < n; ++i) {
        const double an = c0 * (un[i] - m_u[i]) - c2 * m_v[i] - m_a[i];
        m_v[i] += 0.5 * h * (m_a[i] + an);
        m_a[i] = an;
        m_u[i] = un[i];
    }
    m_t += h;
}

void FrameDynamics::step_held (double h, const std::vector<double>& node_loads, double gravity)
{
    AMREX_ALWAYS_ASSERT_WITH_MESSAGE(std::isfinite(h) && h > 0.0, "FrameDynamics::step_held: the step must be finite and positive (s)");
    if (!m_keff || !(std::abs(h - m_h) <= 1.0e-14 * h)) {
        m_keff = effective_stiffness(m_frame, m_a0, m_a1, h);
        m_h = h;
    }
    if (!m_mass) {
        long bad = -1;
        m_mass = mass_factor(m_frame, bad);
        AMREX_ALWAYS_ASSERT_WITH_MESSAGE(m_mass, "FrameDynamics::step_held: the mass matrix is singular (a node without mass or rotary inertia)");
    }
    const std::size_t n = m_u.size();
    // Newmark's K_eff u' = f + M (c0 u + c2 v + a) + C (c1 u + v) with M a = f - C v - K u gives
    // K_eff (u' - u) = 2 (f - K u) + (4/h) M v; then v' = 2 (u' - u)/h - v, and M a' = f - C v' - K u'
    const std::vector<double> f = m_frame.load_vector(node_loads, gravity);
    const std::vector<double> ku = m_frame.apply_stiffness(m_u);
    const std::vector<double> mv = m_frame.apply_mass(m_v);
    std::vector<double> r(n);
    for (std::size_t i = 0; i < n; ++i) { r[i] = 2.0 * (f[i] - ku[i]) + (4.0 / h) * mv[i]; }
    std::vector<double> b = to_free(r);
    m_keff->solve(b);
    const std::vector<double> du = to_full(b);
    for (std::size_t i = 0; i < n; ++i) {
        m_u[i] += du[i];
        m_v[i] = 2.0 * du[i] / h - m_v[i];
    }
    // the acceleration at the end of the step, M a' = f - (a0 M + a1 K) v' - K u'
    std::vector<double> w(n);
    for (std::size_t i = 0; i < n; ++i) { w[i] = m_u[i] + m_a1 * m_v[i]; }
    const std::vector<double> kw = m_frame.apply_stiffness(w);
    const std::vector<double> mv1 = (m_a0 != 0.0) ? m_frame.apply_mass(m_v) : std::vector<double>(n, 0.0);
    for (std::size_t i = 0; i < n; ++i) { r[i] = f[i] - kw[i] - m_a0 * mv1[i]; }
    b = to_free(r);
    m_mass->solve(b);
    m_a = to_full(b);
    m_t += h;
}

double FrameDynamics::kinetic_energy () const { return 0.5 * dot(m_v, m_frame.apply_mass(m_v)); }

double FrameDynamics::strain_energy () const { return 0.5 * dot(m_u, m_frame.apply_stiffness(m_u)); }

std::vector<double> FrameDynamics::state () const
{
    std::vector<double> s;
    s.reserve(1 + 3 * m_u.size());
    s.push_back(m_t);
    s.insert(s.end(), m_u.begin(), m_u.end());
    s.insert(s.end(), m_v.begin(), m_v.end());
    s.insert(s.end(), m_a.begin(), m_a.end());
    return s;
}

bool FrameDynamics::set_state (const std::vector<double>& s)
{
    const std::size_t n = m_frame.num_dofs();
    if (s.size() != 1 + 3 * n) { return false; }
    for (const double x : s) { if (!std::isfinite(x)) { return false; } }
    m_t = s[0];
    m_u.assign(s.begin() + 1, s.begin() + 1 + static_cast<long>(n));
    m_v.assign(s.begin() + 1 + static_cast<long>(n), s.begin() + 1 + 2 * static_cast<long>(n));
    m_a.assign(s.begin() + 1 + 2 * static_cast<long>(n), s.end());
    return true;
}

} // namespace erf_towers
