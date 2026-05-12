#include <ERF.H>
#include <ERF_DataStruct.H>
#include <AMReX_Geometry.H>
#include <AMReX_MultiFab.H>
#include <AMReX_PlotFileUtil.H>
#include "ERF_DA_EnKFSRF.H"

#include <vector>
#include <stdexcept>
#include <cassert>
#include <filesystem>

using namespace amrex;
namespace fs = std::filesystem;


// Simple matrix multiplication
Matrix
matrix_multiply(const Matrix& A, const Matrix& B)
{
    int n = A.size();
    Matrix C(n);

    for (int i = 0; i < n; ++i) {
        for (int j = 0; j < n; ++j) {

            double sum = 0.0;
            for (int k = 0; k < n; ++k) {
                sum += A(i,k) * B(k,j);
            }

            C(i,j) = sum;
        }
    }

    return C;
}

// Print matrix
void matrix_print(const Matrix& A)
{
    int n = A.size();

    for (int i = 0; i < n; ++i) {
        for (int j = 0; j < n; ++j) {
            std::cout << A(i,j) << " ";
        }
        std::cout << "\n";
    }
}

void
lapack_testing()
{
 Matrix S(3);
    S(0,0) = 4.0;  S(0,1) = 1.0;  S(0,2) = 1.0;
    S(1,0) = 1.0;  S(1,1) = 3.0;  S(1,2) = 0.5;
    S(2,0) = 1.0;  S(2,1) = 0.5;  S(2,2) = 2.5;

    Matrix Sinv = S.inverse();
    // fill Sinv ...

    Matrix T = Sinv.cholesky_lower();
    Matrix T_trans = T.transpose();

    Matrix T_T_trans = matrix_multiply(T, T_trans);
    if (ParallelDescriptor::IOProcessor()) {
        matrix_print(T_T_trans);
        std::cout << "T*T_trans done" << std::endl;
    }

    if (ParallelDescriptor::IOProcessor()) {
        std::cout  << "Sinv is " << std::endl;
        matrix_print(Sinv);
    }

    Matrix I_mat = matrix_multiply(S, Sinv);
    if (ParallelDescriptor::IOProcessor()) {
        std::cout << "Checking S*Sinv " << std::endl;
        matrix_print(I_mat);
    }
}


std::vector<std::string>
get_plotfile_list()
{
    std::vector<std::string> pltfiles;
    const std::string member_prefix = "member_";

    std::string pf_dir = member_prefix + "00/plotfiles";

    if (!fs::exists(pf_dir)) {
        amrex::Abort("Plotfile directory not found: " + pf_dir);
    }

    for (const auto& entry : fs::directory_iterator(pf_dir)) {
        if (entry.is_directory()) {
            std::string name = entry.path().filename().string();
            if (name.rfind("plt", 0) == 0) {  // starts with "plt"
                pltfiles.push_back(name);
            }
        }
    }

    std::sort(pltfiles.begin(), pltfiles.end());
    return pltfiles;
}

// Reads the plotfile data into cell cenetred multifab
// Does not fill ghost cells
void
read_plot_file(PlotFileData& pf,
               const std::vector<std::string> varnames,
               MultiFab& mf)
{
    // ------------------------------------------------------------
    // Open plotfile
    // ------------------------------------------------------------
    const std::vector<std::string>& var_names_pf = pf.varNames();

    // ------------------------------------------------------------
    // Validate requested variables
    // ------------------------------------------------------------
    for (auto const& v : varnames) {
        bool found = false;
        for (auto const& vpf : var_names_pf) {
            if (v == vpf) {
                found = true;
                break;
            }
        }
        if (!found) {
            Abort("read_plot_file: invalid variable name: " + v);
        }
    }

    // ------------------------------------------------------------
    // Define destination MultiFab (single level only)
    // ------------------------------------------------------------
    const int level = 0;

    BoxArray ba = pf.boxArray(level);
    DistributionMapping dm(ba);

    int ncomp = varnames.size();

    mf.define(ba, dm, ncomp, 0);

    // ------------------------------------------------------------
    // Copy plotfile data → mf
    // ------------------------------------------------------------
    for (int comp = 0; comp < ncomp; ++comp)
    {
        const MultiFab& src = pf.get(level, varnames[comp]);
        MultiFab::Copy(mf, src, 0, comp, 1, 0);
    }
}

MultiFab
read_member_multifab(int n,
                     const std::string& pf_name,
                     const Vector<std::string>& varnames)
{
    const std::string member_prefix = "member_";
    std::string member_dir = member_prefix + amrex::Concatenate("", n, 2);
    std::string pf_path = member_dir + "/plotfiles/" + pf_name;

    PlotFileData pf(pf_path);

    const BoxArray& ba = pf.boxArray(0);
    const DistributionMapping& dm = pf.DistributionMap(0);
    int ncomp = varnames.size();

    MultiFab mf(ba, dm, ncomp, 0);

    read_plot_file(pf, varnames, mf);

    return mf;
}

MultiFab
compute_ensemble_mean(int Nens,
                      const std::string& pf_name,
                      const Vector<std::string>& varnames)
{
    MultiFab mf_mean;
    bool initialized = false;

    for (int n = 0; n < Nens; ++n)
    {
        MultiFab mf_tmp = read_member_multifab(n, pf_name, varnames);

        if (!initialized) {
            mf_mean.define(mf_tmp.boxArray(),
                           mf_tmp.DistributionMap(),
                           mf_tmp.nComp(),
                           mf_tmp.nGrow());
            mf_mean.setVal(0.0);
            initialized = true;
        }

        MultiFab::Add(mf_mean, mf_tmp, 0, 0, mf_tmp.nComp(), mf_tmp.nGrow());
    }

    mf_mean.mult(1.0 / Real(Nens));
    return mf_mean;
}

void
compute_d_vec(const MultiFab& mf1,
              const MultiFab& mf2,
              MultiFab& d_vec)
{
    AMREX_ALWAYS_ASSERT(mf1.boxArray() == mf2.boxArray());
    AMREX_ALWAYS_ASSERT(mf1.DistributionMap() == mf2.DistributionMap());
    AMREX_ALWAYS_ASSERT(mf1.nComp() == mf2.nComp());
    AMREX_ALWAYS_ASSERT(mf1.nGrowVect() == mf2.nGrowVect());

    d_vec.define(mf1.boxArray(),
                 mf1.DistributionMap(),
                 mf1.nComp(),
                 mf1.nGrowVect());

    MultiFab::Copy(d_vec,
                   mf1,
                   0, 0,
                   mf1.nComp(),
                   mf1.nGrowVect());

    MultiFab::Subtract(d_vec,
                       mf2,
                       0, 0,
                       mf1.nComp(),
                       mf1.nGrowVect());
}

void
compute_d_prime_vec(MultiFab& d_prime_vec,
                    const MultiFab& d_vec,
                    const Vector<Real>& R_diag)
{
    AMREX_ALWAYS_ASSERT(d_vec.nComp() == R_diag.size());

    const int ncomp = d_vec.nComp();

    d_prime_vec.define(d_vec.boxArray(),
                       d_vec.DistributionMap(),
                       ncomp,
                       d_vec.nGrowVect());

    amrex::Gpu::DeviceVector<Real> R_diag_d(R_diag.size());
    amrex::Gpu::copy(amrex::Gpu::hostToDevice, R_diag.begin(), R_diag.end(), R_diag_d.begin());
    const Real* R_diag_d_ptr = R_diag_d.data();

    for (MFIter mfi(d_vec, TilingIfNotGPU()); mfi.isValid(); ++mfi)
    {
        const Box& bx = mfi.tilebox();
        const Array4<const Real> d_arr = d_vec.const_array(mfi);
        const Array4<Real> dp_arr = d_prime_vec.array(mfi);

        amrex::ParallelFor(bx, ncomp,
                            [=] AMREX_GPU_DEVICE (int i, int j, int k, int n)
        {
            dp_arr(i,j,k,n) = d_arr(i,j,k,n) / R_diag_d_ptr[n];
        });
    }
}

void
compute_mean_H_xf(MultiFab& mean_H_xf,
                  const int Nens,
                  const std::string& last_pf_name,
                  const Vector<std::string>& varnames)
{

    for (int n = 0; n < Nens; ++n)
    {
        MultiFab xf_i = read_member_multifab(n, last_pf_name, varnames);
        if (n==0) {
            mean_H_xf.define(xf_i.boxArray(),
                             xf_i.DistributionMap(),
                             xf_i.nComp(),
                             xf_i.nGrow());
            mean_H_xf.setVal(0.0);
        }

        MultiFab H_xf_i;

        Apply_H(xf_i, H_xf_i);

        if(n==0){
            MultiFab sum_H_xf_i(H_xf_i.boxArray(), H_xf_i.DistributionMap(), H_xf_i.nComp(), H_xf_i.nGrow());
        }

        MultiFab::Add(mean_H_xf, H_xf_i, 0, 0, H_xf_i.nComp(), H_xf_i.nGrow());
    }

    mean_H_xf.mult(1.0 / Nens);
}

Real
Compute_yf_prime_i_T_Rinv_yf_prime_j(const MultiFab& Yi,
                                   const MultiFab& Yj,
                                   const Vector<Real>& R_diag)
{
    const int ncomp = Yi.nComp();

    ReduceOps<ReduceOpSum> reduce_op;
    ReduceData<Real> reduce_data(reduce_op);

    using ReduceTuple = typename decltype(reduce_data)::Type;

    amrex::Gpu::DeviceVector<Real> R_diag_d(R_diag.size());
    amrex::Gpu::copy(amrex::Gpu::hostToDevice, R_diag.begin(), R_diag.end(), R_diag_d.begin());
    const Real* R_diag_d_ptr = R_diag_d.data();

    for (MFIter mfi(Yi, TilingIfNotGPU()); mfi.isValid(); ++mfi)
    {
        const Box& bx = mfi.tilebox();

        auto const& yi = Yi.const_array(mfi);
        auto const& yj = Yj.const_array(mfi);

        reduce_op.eval(bx, reduce_data,
        [=] AMREX_GPU_DEVICE (int i, int j, int k) -> ReduceTuple
        {
            Real local = 0.0;

            for (int n = 0; n < ncomp; ++n) {
                local += yi(i,j,k,n) * yj(i,j,k,n)/R_diag_d_ptr[n];
            }

            return {local};
        });
    }

    ReduceTuple hv = reduce_data.value();
    return amrex::get<0>(hv);
}

void
compute_yf_prime (int i,
                const std::string& last_pf_name,
                const Vector<std::string>& varnames,
                const MultiFab& mean_H_xf,
                MultiFab& yf_prime_i)
{
    // Read ensemble member
    MultiFab xf_i = read_member_multifab(i, last_pf_name, varnames);

    // Apply observation operator
    MultiFab H_xf_i;
    Apply_H(xf_i, H_xf_i);

    // Allocate output
    yf_prime_i.define(H_xf_i.boxArray(),
                      H_xf_i.DistributionMap(),
                      H_xf_i.nComp(),
                      H_xf_i.nGrow());

    // y'_i = H(xf_i) - mean(H(xf))
    MultiFab::Copy(yf_prime_i,
                   H_xf_i,
                   0, 0,
                   H_xf_i.nComp(),
                   H_xf_i.nGrow());

    MultiFab::Subtract(yf_prime_i,
                       mean_H_xf,
                       0, 0,
                       H_xf_i.nComp(),
                       H_xf_i.nGrow());
}

void
compute_S_matrix(Matrix& S,
                 const int& Nens,
                 const MultiFab& mean_H_xf,
                 const Vector<Real>& R_diag,
                 const std::string& last_pf_name,
                 const Vector<std::string>& varnames)
{
    for (int i = 0; i < Nens; ++i) {
        MultiFab yf_prime_i;
        compute_yf_prime(i, last_pf_name, varnames, mean_H_xf, yf_prime_i);

        for (int j = 0; j < Nens; ++j) {
            MultiFab yf_prime_j;
            compute_yf_prime(j, last_pf_name, varnames, mean_H_xf, yf_prime_j);
            Real val = Compute_yf_prime_i_T_Rinv_yf_prime_j(yf_prime_i, yf_prime_j, R_diag);
            S(i,j) = val/(Nens-1);
            if(i==j) {
                S(i,j) = 1.0 + S(i,j);
            }
        }
    }
}

Real
compute_yf_prime_T_d_prime_vec (const MultiFab& yf_prime,
                                const MultiFab& d_prime_vec)
{
    AMREX_ALWAYS_ASSERT(yf_prime.nComp() ==
                        d_prime_vec.nComp());

    const int ncomp = yf_prime.nComp();

    ReduceOps<ReduceOpSum> reduce_op;
    ReduceData<Real> reduce_data(reduce_op);

    using ReduceTuple =
        typename decltype(reduce_data)::Type;

    for (MFIter mfi(yf_prime, TilingIfNotGPU());
         mfi.isValid(); ++mfi)
    {
        const Box& bx = mfi.tilebox();

        auto const& yf_arr =
            yf_prime.const_array(mfi);

        auto const& d_arr =
            d_prime_vec.const_array(mfi);

        reduce_op.eval(
        bx,
        reduce_data,
        [=] AMREX_GPU_DEVICE (int i, int j, int k)
            -> ReduceTuple
        {
            Real local = 0.0;

            for (int n = 0; n < ncomp; ++n)
            {
                local +=
                    yf_arr(i,j,k,n) *
                    d_arr(i,j,k,n);
            }

            return {local};
        });
    }

    ReduceTuple hv = reduce_data.value();

    return amrex::get<0>(hv);
}

void
compute_r_vec (int Nens,
               const std::string& last_pf_name,
               const Vector<std::string>& varnames,
               const MultiFab& mean_H_xf,
               const MultiFab& d_prime_vec,
               Vector<Real>& r_vec)
{
    r_vec.resize(Nens);

    for (int i = 0; i < Nens; ++i)
    {
        MultiFab yf_prime_i;

        compute_yf_prime(i, last_pf_name, varnames, mean_H_xf, yf_prime_i);

        r_vec[i] = compute_yf_prime_T_d_prime_vec(yf_prime_i, d_prime_vec);
    }
}


void
compute_alpha_vec (const int& Nens,
                   const Matrix& S,
                   const Vector<Real>& r_vec,
                   Vector<Real>& alpha_vec)
{
    AMREX_ALWAYS_ASSERT(r_vec.size() == Nens);

    Matrix Sinv = S.inverse();

    alpha_vec.resize(Nens, 0.0);

    const Real fac = 1.0_rt / static_cast<Real>(Nens - 1);

    for (int i = 0; i < Nens; ++i)
    {
        Real sum = 0.0_rt;

        for (int j = 0; j < Nens; ++j)
        {
            sum += Sinv(i,j) * r_vec[j];
        }

        alpha_vec[i] = fac * sum;
    }
}

// Perform a matrix time a small vector multiply as a summation of
// column-vector element multiply

void
compute_Xf_prime_times_vector (const int Nens,
                            const std::string& last_pf_name,
                            const Vector<std::string>& varnames,
                            const MultiFab& xf_bar,
                            const Vector<Real>& vec_in,
                            MultiFab& result)
{
    AMREX_ALWAYS_ASSERT(vec_in.size() == Nens);

    // Read first member to define result
    MultiFab xf_0 =
        read_member_multifab(0,
                             last_pf_name,
                             varnames);

    result.define(xf_0.boxArray(),
                  xf_0.DistributionMap(),
                  xf_0.nComp(),
                  xf_0.nGrow());

    result.setVal(0.0);

    for (int n = 0; n < Nens; ++n)
    {
        MultiFab xf_n =
            read_member_multifab(n, last_pf_name, varnames);

        // xf_n_prime = xf_n - xf_bar
        MultiFab xf_n_prime(xf_n.boxArray(),
                            xf_n.DistributionMap(),
                            xf_n.nComp(),
                            xf_n.nGrow());

        MultiFab::Copy(xf_n_prime,
                       xf_n,
                       0, 0,
                       xf_n.nComp(),
                       xf_n.nGrow());

        MultiFab::Subtract(xf_n_prime,
                           xf_bar,
                           0, 0,
                           xf_n.nComp(),
                           xf_n.nGrow());

        // result += vec_in[n] * xf_n_prime
        MultiFab::Saxpy(result,
                        vec_in[n],
                        xf_n_prime,
                        0, 0,
                        xf_n.nComp(),
                        xf_n.nGrow());
    }
}

void
add_multifabs (const MultiFab& xf_bar,
              const MultiFab& Xf_prime_alpha,
              MultiFab& result)
{
    result.define(
        xf_bar.boxArray(),
        xf_bar.DistributionMap(),
        xf_bar.nComp(),
        xf_bar.nGrowVect()
    );
    result.setVal(0.0);
    MultiFab::Copy(result, xf_bar, 0, 0, xf_bar.nComp(), 0);

    MultiFab::Add(result, Xf_prime_alpha,
                  0, 0, Xf_prime_alpha.nComp(), 0);
}

void
compute_T_matrix (const Matrix& S_mat,
                  Matrix& T_mat)
{
    Matrix Sinv = S_mat.inverse();
    // fill Sinv ...

    T_mat = Sinv.cholesky_lower();
}


void
update_ensemble (const int Nens,
                 const std::string& last_pf_name,
                 const Vector<std::string>& varnames,
                 const MultiFab& xf_bar,
                 const Matrix& T,
                 const int& n,
                 MultiFab& result)
{
    Vector<Real> T_colvec(Nens);
    for (int i = 0; i < Nens; ++i) {
        T_colvec[i] = T(i,n);
    }
    compute_Xf_prime_times_vector(Nens, last_pf_name, varnames, xf_bar, T_colvec, result);
}
