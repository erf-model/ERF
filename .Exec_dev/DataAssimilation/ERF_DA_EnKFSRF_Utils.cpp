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

// Simple matrix multiplication
Matrix matrix_multiply(const Matrix& A, const Matrix& B)
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
        MultiFab H_xf_i;
        Apply_H(xf_i, H_xf_i);

        if(n==0){
            MultiFab sum_H_xf_i(H_xf_i.boxArray(), H_xf_i.DistributionMap(), H_xf_i.nComp(), H_xf_i.nGrow());
        }

        MultiFab::Add(mean_H_xf, xf_i, 0, 0, H_xf_i.nComp(), H_xf_i.nGrow());
    }

    mean_H_xf.mult(1.0 / Nens);
}

Real
Compute_y_prime_i_T_Rinv_y_prime_j(const MultiFab& Yi,
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
compute_S_matrix(Matrix& S,
                 const int& Nens,
                 const MultiFab& mean_H_xf,
                 const Vector<Real>& R_diag,
                 const std::string& last_pf_name,
                 const Vector<std::string>& varnames)
{
    for (int i = 0; i < Nens; ++i) {
        MultiFab yf_prime_i;
        MultiFab xf_i = read_member_multifab(i, last_pf_name, varnames);
        MultiFab H_xf_i;
        Apply_H(xf_i, H_xf_i);
        MultiFab y_prime_i(xf_i.boxArray(), xf_i.DistributionMap(), 2, xf_i.nGrow());
        MultiFab::Subtract(y_prime_i, mean_H_xf, 0, 0, 2, mean_H_xf.nGrow());

        for (int j = 0; j < Nens; ++j) {
            MultiFab yf_prime_j;
            MultiFab xf_j = read_member_multifab(j, last_pf_name, varnames);
            MultiFab H_xf_j;
            Apply_H(xf_j, H_xf_j);

            MultiFab y_prime_j(mean_H_xf.boxArray(), mean_H_xf.DistributionMap(), 2, mean_H_xf.nGrow());
            // y_prime = H_xf_i
            MultiFab::Copy(y_prime_i, H_xf_j, 0, 0, 2, H_xf_i.nGrow());

           // y_prime -= mean_H_xf
            MultiFab::Subtract(y_prime_j, mean_H_xf, 0, 0, 2, H_xf_j.nGrow());
            Real val = Compute_y_prime_i_T_Rinv_y_prime_j(yf_prime_i, y_prime_j, R_diag);
            S(i,j) = val;
        }
    }
}

