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

// Compute RMS error between ground truth and simulation plotfiles.
static Real
ComputeRMS(const MultiFab& mf_truth,
           const MultiFab& mf_solution)
{
    if (mf_truth.nComp() != mf_solution.nComp()) {
        Abort("Ground truth and solution have different component counts.");
    }

    if (mf_truth.boxArray().size() == 0 ||
        mf_solution.boxArray().size() == 0) {
        Abort("Empty MultiFab.");
    }

    // Copy the solution onto the ground-truth BoxArray and
    // DistributionMapping. This permits different box decompositions,
    // provided both MultiFabs represent the same mesh resolution/domain.
    MultiFab mf_solution_aligned(
        mf_truth.boxArray(),
        mf_truth.DistributionMap(),
        mf_truth.nComp(),
        0);

    mf_solution_aligned.ParallelCopy(
        mf_solution, 0, 0, mf_truth.nComp(), 0, 0);

    // Difference = ground truth - solution.
    MultiFab dvec(
        mf_truth.boxArray(),
        mf_truth.DistributionMap(),
        mf_truth.nComp(),
        0);

    MultiFab::Copy(dvec, mf_truth, 0, 0, mf_truth.nComp(), 0);
    MultiFab::Subtract(dvec, mf_solution_aligned, 0, 0,
                       mf_truth.nComp(), 0);

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
        // Input lists.
        const std::string truth_list_file =
            "ground_truth_plt_files.txt";

        const std::string solution_list_file =
            "coarse_simulation_plt_files.txt";

        const std::string vars_file = "vars.txt";

        const std::string output_file = "rms_error.txt";

        // Read both lists.
        const std::vector<std::string> truth_files =
            ReadPlotFileList(truth_list_file);

        const std::vector<std::string> solution_files =
            ReadPlotFileList(solution_list_file);

        if (truth_files.size() != solution_files.size()) {
            Abort("Ground-truth and simulation lists have different "
                  "numbers of plotfiles.");
        }

        if (truth_files.empty()) {
            Abort("No plotfiles found in the input lists.");
        }

        if (ParallelDescriptor::IOProcessor()) {
            Print() << "Number of plotfile pairs: "
                    << truth_files.size() << '\n';
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
             iter < truth_files.size(); ++iter)
        {
            if (ParallelDescriptor::IOProcessor()) {
                Print() << "\nPair " << iter << '\n';
                Print() << "  Ground truth: "
                        << truth_files[iter] << '\n';
                Print() << "  Solution:     "
                        << solution_files[iter] << '\n';
            }

            // Read ground-truth plotfile.
            MultiFab mf_cc_truth;

            PlotFileData pf_truth(truth_files[iter]);

            ReadPlotFile(
                vars_file, pf_truth, mf_cc_truth);

            // Read simulation plotfile.
            MultiFab mf_cc_solution;

            PlotFileData pf_solution(solution_files[iter]);

            ReadPlotFile(
                vars_file, pf_solution, mf_cc_solution);

            // Compute RMS error for this time instant.
            const Real rms =
                ComputeRMS(mf_cc_truth, mf_cc_solution);

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
