#ifndef ERF_DA_EnKFSRF_H
#define ERF_DA_EnKFSRF_H

#include <ERF.H>
#include <ERF_DataStruct.H>
#include <AMReX_Geometry.H>
#include <AMReX_MultiFab.H>
#include <AMReX_PlotFileUtil.H>
void
read_plot_file(amrex::PlotFileData& pf,
               const std::vector<std::string> varnames,
               amrex::MultiFab& mf);

amrex::MultiFab
read_member_multifab(int n,
                     const std::string& pf_name,
                     const amrex::Vector<std::string>& varnames);

amrex::MultiFab
compute_ensemble_mean(const std::string& pf_name,
                      int Nens,
                      const amrex::Vector<std::string>& varnames);


extern "C" {
    void dgetrf_(int* M, int* N, double* A, int* lda, int* ipiv, int* info);
    void dgetri_(int* N, double* A, int* lda, int* ipiv,
                 double* work, int* lwork, int* info);
    void dpotrf_(char* UPLO, int* N, double* A, int* LDA, int* INFO);
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
// Simple matrix multiplication
Matrix matrix_multiply(const Matrix& A, const Matrix& B);
void matrix_print(const Matrix& A);

void lapack_testing();
#endif
