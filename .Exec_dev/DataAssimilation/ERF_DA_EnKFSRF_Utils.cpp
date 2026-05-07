#include <ERF.H>
#include <ERF_DataStruct.H>
#include <AMReX_Geometry.H>
#include <AMReX_MultiFab.H>
#include <AMReX_PlotFileUtil.H>
#include "ERF_DA_EnKFSRF.H"

#include <vector>
#include <stdexcept>
#include <cassert>


using namespace amrex;
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
