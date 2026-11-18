#ifndef ERF_DA_EnKFSRF_H
#define ERF_DA_EnKFSRF_H

#include <ERF.H>
#include <ERF_DataStruct.H>
#include <AMReX_Geometry.H>
#include <AMReX_MultiFab.H>
#include <AMReX_PlotFileUtil.H>

extern "C" {
    void dgetrf_(int* M, int* N, double* A, int* lda, int* ipiv, int* info);
    void dgetri_(int* N, double* A, int* lda, int* ipiv,
                 double* work, int* lwork, int* info);
    void dpotrf_(char* UPLO, int* N, double* A, int* LDA, int* INFO);
    void dsyev_(char* JOBZ, char* UPLO, int* N, double* A, int* LDA,
                double* W, double* WORK, int* LWORK, int* INFO);
}

class Matrix
{
public:
    Matrix(int n)
        : n_(n), data_(n * n, 0.0)
    {}

    Matrix(int n, const std::vector<double>& in)
        : n_(n), data_(in)
    {
        if (data_.size() != n * n)
            throw std::runtime_error("Matrix size mismatch");
    }

    // Fortran-style access: A(i,j)
    double& operator()(int i, int j)
    {
        return data_[j * n_ + i];   // column-major
    }

    const double& operator()(int i, int j) const
    {
        return data_[j * n_ + i];
    }

    int size() const { return n_; }

    const double* data() const { return data_.data(); }
    double* data() { return data_.data(); }

    // Inverse (returns new Matrix)
    Matrix inverse() const
    {
        int n = n_;
        Matrix Ainv = *this;  // copy

        std::vector<int> ipiv(n);
        int lda = n;
        int info;

        // LU factorization
        dgetrf_(&n, &n, Ainv.data(), &lda, ipiv.data(), &info);
        if (info != 0)
            throw std::runtime_error("LU factorization failed");

        // workspace query
        int lwork = -1;
        double work_query;

        dgetri_(&n, Ainv.data(), &lda, ipiv.data(),
                &work_query, &lwork, &info);

        lwork = static_cast<int>(work_query);
        std::vector<double> work(lwork);

        // inversion
        dgetri_(&n, Ainv.data(), &lda, ipiv.data(),
                work.data(), &lwork, &info);

        if (info != 0)
            throw std::runtime_error("Matrix inversion failed");

        return Ainv;
    }

    Matrix cholesky_lower() const
    {
        int n = n_;
        Matrix L = *this;  // copy (LAPACK overwrites)

        int lda = n;
        int info;
        char uplo = 'L';   // lower triangular

        dpotrf_(&uplo, &n, L.data(), &lda, &info);

        if (info != 0) {
            throw std::runtime_error("Cholesky failed (matrix not SPD?)");
        }

        // Zero out upper part for clarity
        for (int i = 0; i < n; ++i) {
            for (int j = i+1; j < n; ++j) {
                L(i,j) = 0.0;
            }
        }

        return L;
    }

    // Symmetric inverse square root A^{-1/2} = V Lambda^{-1/2} V^T of a
    // symmetric positive definite matrix, via eigendecomposition (dsyev).
    // Unlike a Cholesky factor, the result is symmetric, so for the ETKF
    // matrix S (which has S*1 = 1) it satisfies T*1 = 1 and keeps the
    // analysis perturbations zero-mean.
    Matrix inverse_sqrt_symmetric() const
    {
        int n = n_;

        // Symmetrize to remove round-off asymmetry; dsyev reads the lower triangle
        Matrix V(n);
        for (int i = 0; i < n; ++i) {
            for (int j = 0; j < n; ++j) {
                V(i,j) = 0.5 * ((*this)(i,j) + (*this)(j,i));
            }
        }

        std::vector<double> w(n);
        int lda = n;
        int info;
        char jobz = 'V';   // eigenvalues and eigenvectors
        char uplo = 'L';

        // workspace query
        int lwork = -1;
        double work_query;
        dsyev_(&jobz, &uplo, &n, V.data(), &lda, w.data(),
               &work_query, &lwork, &info);

        lwork = static_cast<int>(work_query);
        std::vector<double> work(lwork);

        // V is overwritten with orthonormal eigenvectors (columns)
        dsyev_(&jobz, &uplo, &n, V.data(), &lda, w.data(),
               work.data(), &lwork, &info);

        if (info != 0) {
            throw std::runtime_error("Symmetric eigendecomposition failed");
        }

        for (int k = 0; k < n; ++k) {
            if (w[k] <= 0.0) {
                throw std::runtime_error("Matrix not positive definite in inverse_sqrt_symmetric");
            }
        }

        // T = V * diag(w^{-1/2}) * V^T
        Matrix T(n);
        for (int i = 0; i < n; ++i) {
            for (int j = 0; j < n; ++j) {
                double sum = 0.0;
                for (int k = 0; k < n; ++k) {
                    sum += V(i,k) * V(j,k) / std::sqrt(w[k]);
                }
                T(i,j) = sum;
            }
        }
        return T;
    }

    Matrix transpose() const
    {
        int n = n_;
        Matrix At(n);

        for (int i = 0; i < n; ++i) {
            for (int j = 0; j < n; ++j) {
                At(i,j) = (*this)(j,i);
            }
        }
        return At;
    }

private:
    int n_;
    std::vector<double> data_;
};


void
read_plot_file(amrex::PlotFileData& pf,
               const std::vector<std::string> varnames,
               amrex::MultiFab& mf);

amrex::MultiFab
read_member_multifab(const int da_iter,
                     const int n,
                     const std::string& pf_name,
                     const amrex::Vector<std::string>& varnames);

amrex::MultiFab
compute_ensemble_mean(const int da_iter,
                      const int Nens,
                      const std::string& pf_name,
                      const amrex::Vector<std::string>& varnames);

void Apply_H(const amrex::MultiFab& x_mf, amrex::MultiFab& y_mf);

void read_in_observations(const int& da_iter, amrex::MultiFab& y_obs);

std::vector<std::string> get_plotfile_list(const int da_iter);
void read_in_observations(const int& da_iter,
                          const amrex::Vector<std::string>& varnames,
                          amrex::MultiFab& y_obs);
void compute_d_vec(const amrex::MultiFab& mf1,
                   const amrex::MultiFab& mf2,
                   amrex::MultiFab& d_vec);

void compute_R_diag_vals(amrex::Vector<amrex::Real>& R_diag);

void compute_d_prime_vec(amrex::MultiFab& d_prime_vec,
                    const amrex::MultiFab& d_vec,
                    const amrex::Vector<amrex::Real>& R_diag);

void compute_mean_H_xf(const int da_iter,
                       amrex::MultiFab& mean_H_xf,
                       const int Nens,
                       const std::string& last_pf_name,
                       const amrex::Vector<std::string>& varnames);

void compute_r_vec (const int da_iter,
                    int Nens,
                    const std::string& last_pf_name,
                    const amrex::Vector<std::string>& varnames,
                    const amrex::MultiFab& mean_H_xf,
                    const amrex::MultiFab& d_prime_vec,
                    amrex::Vector<amrex::Real>& r_vec);

void compute_alpha_vec (const int& Nens,
                        const Matrix& S,
                        const amrex::Vector<amrex::Real>& r_vec,
                        amrex::Vector<amrex::Real>& alpha_vec);

void compute_Xf_prime_times_vector (const int da_iter,
                                    const int Nens,
                                    const std::string& last_pf_name,
                                    const amrex::Vector<std::string>& varnames,
                                    const amrex::MultiFab& xf_bar,
                                    const amrex::Vector<amrex::Real>& vec_in,
                                    amrex::MultiFab& result);

void update_ensemble (const int da_iter,
                      const int Nens,
                      const std::string& last_pf_name,
                      const amrex::Vector<std::string>& varnames,
                      const amrex::MultiFab& xf_bar,
                      const Matrix& T,
                      const int& n,
                      amrex::MultiFab& result);

void add_multifabs (const amrex::MultiFab& xf_bar,
                    const amrex::MultiFab& Xf_prime_alpha,
                    amrex::MultiFab& result);

void compute_T_matrix (const Matrix& S_mat,
                        Matrix& T_mat);

// Simple matrix multiplication
Matrix matrix_multiply(const Matrix& A, const Matrix& B);
void matrix_print(const Matrix& A);
void lapack_testing();

void compute_S_matrix(const int da_iter,
                      Matrix& S,
                      const int& Nens,
                      const amrex::MultiFab& mean_H_xf,
                      const amrex::Vector<amrex::Real>& R_diag,
                      const std::string& last_pf_name,
                      const amrex::Vector<std::string>& varnames);
std::string
MakeEnsembleCheckpointName (int da_iter, int ens_no);

std::string
MakeFullPath (const std::string& filename, const int da_iter);

void
WriteInnovationRMS(const amrex::MultiFab& dvec,
                   const std::string& filename,
                   int da_iter);
#endif
