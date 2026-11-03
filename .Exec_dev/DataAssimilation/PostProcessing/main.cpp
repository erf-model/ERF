#include <AMReX.H>
#include <AMReX_MultiFab.H>
#include <AMReX_Geometry.H>
#include <AMReX_PlotFileUtil.H>
#include <AMReX_Reduce.H>
#include <AMReX_ParallelDescriptor.H>

#include "ReadPlotFile.H"

#include <fstream>
#include <iostream>
#include <iomanip>
#include <string>
#include <vector>
#include <cmath>
#include <stdexcept>

using namespace amrex;

// Read plotfile paths from a text file, one path per line.
static std::vector<std::string>
ReadPlotFileList(const std::string& filename)
{
    std::ifstream in(filename);

    if (!in.is_open()) {
        Abort("Cannot open plotfile list: " + filename);
    }

    std::vector<std::string> files;
    std::string line;

    while (std::getline(in, line)) {
        // Ignore empty lines.
        if (!line.empty()) {
            files.push_back(line);
        }
    }

    return files;
}

// Compute RMS error between plotfiles.
static Real
ComputeRMS(const MultiFab& mf_1,
           const MultiFab& mf_2)
{
    if (mf_1.nComp() != mf_2.nComp()) {
        Abort("The plot files in list 1 and list 2 have different component counts.");
    }

    if (mf_1.boxArray().size() == 0 ||
        mf_2.boxArray().size() == 0) {
        Abort("Empty MultiFab.");
    }

    // Copy onto BoxArray and
    // DistributionMapping. This permits different box decompositions,
    // provided both MultiFabs represent the same mesh resolution/domain.
    MultiFab mf_2_aligned(
        mf_1.boxArray(),
        mf_1.DistributionMap(),
        mf_1.nComp(),
        0);

    mf_2_aligned.ParallelCopy(
        mf_2, 0, 0, mf_1.nComp(), 0, 0);

    MultiFab dvec(
        mf_1.boxArray(),
        mf_1.DistributionMap(),
        mf_1.nComp(),
        0);

    MultiFab::Copy(dvec, mf_1, 0, 0, mf_1.nComp(), 0);
    MultiFab::Subtract(dvec, mf_2_aligned, 0, 0,
                       mf_1.nComp(), 0);

    Real sum_sq = 0.0;
    long count = 0;

    for (MFIter mfi(dvec); mfi.isValid(); ++mfi)
    {
        const auto& arr = dvec.const_array(mfi);
        const Box& box = mfi.validbox();

        ReduceOps<ReduceOpSum> reduce_op;
        ReduceData<Real> reduce_data(reduce_op);
        using ReduceTuple = typename decltype(reduce_data)::Type;

        const int ncomp = dvec.nComp();

        reduce_op.eval(
            box, reduce_data,
            [=] AMREX_GPU_DEVICE (int i, int j, int k) noexcept
                -> ReduceTuple
            {
                Real local_sum = 0.0;

                for (int n = 0; n < ncomp; ++n) {
                    const Real val = arr(i, j, k, n);
                    local_sum += val * val;
                }

                return {local_sum};
            });

        sum_sq += get<0>(reduce_data.value());

        count += static_cast<long>(box.numPts()) * ncomp;
    }

    // Global reductions across MPI ranks.
    ParallelDescriptor::ReduceRealSum(sum_sq);
    ParallelDescriptor::ReduceLongSum(count);

    return (count > 0)
        ? std::sqrt(sum_sq / static_cast<Real>(count))
        : 0.0;
}

int main(int argc, char* argv[])
{
    amrex::Initialize(argc, argv);

    int exit_code = 0;

    {

        std::string a_file_list_1;
        std::string a_file_list_2;

        for (int i = 1; i < argc; ++i) {
            std::string arg(argv[i]);

            if (arg.find("--file-list-1=") == 0) {
                a_file_list_1 = arg.substr(14);
            }
            else if (arg.find("--file-list-2=") == 0) {
                a_file_list_2 = arg.substr(14);
            }
        }

        if (a_file_list_1.empty() || a_file_list_2.empty()) {
            amrex::Print()
                << "Usage: " << argv[0]
                << " --f1ile-list-1=<list1>"
                << " --f1ile-list-2=<list2>\n";
            amrex::Finalize();
            return 1;
        }

        const std::string vars_file = "vars.txt";
        const std::string output_file = "rms_error.txt";

        // Read both lists.
        const std::vector<std::string> file_list_1 =
            ReadPlotFileList(a_file_list_1);

        const std::vector<std::string> file_list_2 =
            ReadPlotFileList(a_file_list_2);

        if (file_list_1.size() != file_list_2.size()) {
            Abort("file_list_1 and file_list_2 have different "
                  "numbers of plotfiles.");
        }

        if (file_list_1.empty()) {
            Abort("No plotfiles found in the input lists.");
        }

        if (ParallelDescriptor::IOProcessor()) {
            Print() << "Number of plotfile pairs: "
                    << file_list_1.size() << '\n';
            Print() << "Writing RMS errors to: "
                    << output_file << '\n';
        }

        // Open output file on the I/O rank.
        std::ofstream out;

        if (ParallelDescriptor::IOProcessor()) {
            out.open(output_file);

            if (!out.is_open()) {
                Abort("Cannot open output file: " + output_file);
            }

            out << "# iter rms_error\n";
            out << std::setprecision(16);
        }

        // Read corresponding plotfiles and compute RMS error.
        for (std::size_t iter = 0;
             iter < file_list_1.size(); ++iter)
        {
            if (ParallelDescriptor::IOProcessor()) {
                Print() << "\nPair " << iter << '\n';
                Print() << "  File list 1: "
                        << file_list_1[iter] << '\n';
                Print() << "  File list 2:     "
                        << file_list_2[iter] << '\n';
            }

            MultiFab mf_cc_1;

            PlotFileData pf_1(file_list_1[iter]);

            ReadPlotFile(
                vars_file, pf_1, mf_cc_1);

            // Read simulation plotfile.
            MultiFab mf_cc_2;

            PlotFileData pf_2(file_list_2[iter]);

            ReadPlotFile(
                vars_file, pf_2, mf_cc_2);

            // Compute RMS error for this time instant.
            const Real rms =
                ComputeRMS(mf_cc_1, mf_cc_2);

            // Write: iteration, RMS error.
            if (ParallelDescriptor::IOProcessor()) {
                out << iter << " "
                    << std::setprecision(16)
                    << rms << '\n';

                out.flush();

                Print() << "  RMS error = "
                        << std::setprecision(16)
                        << rms << '\n';
            }
        }

        if (ParallelDescriptor::IOProcessor()) {
            out.close();
            Print() << "\nFinished computing RMS errors.\n";
        }
    }

    amrex::Finalize();

    return exit_code;
}
