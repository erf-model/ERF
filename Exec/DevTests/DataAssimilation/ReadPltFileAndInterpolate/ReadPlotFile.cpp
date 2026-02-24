#include <fstream>
#include <sstream>

#include <AMReX.H>
#include <AMReX_Vector.H>
#include <AMReX_MultiFab.H>
#include <AMReX_PlotFileUtil.H>
#include <AMReX_DistributionMapping.H>
#include <AMReX_BoxArray.H>

#include "ReadPlotFile.H"

#include <vector>
#include <algorithm>  // for std::remove_if, std::isspace

using namespace amrex;
// ------------------------------------------------------------
// Read variable names from a file
// ------------------------------------------------------------
Vector<std::string> 
ReadVarNames(const std::string& filename)
{
    Vector<std::string> varnames;
    std::ifstream infile(filename);
    if (!infile.is_open()) {
        amrex::Abort("Cannot open variable file: " + filename);
    }

    std::string line;
    while (std::getline(infile, line)) {
        // trim whitespace
        line.erase(line.begin(),
                   std::find_if(line.begin(), line.end(),
                                [](unsigned char ch){ return !std::isspace(ch); }));
        line.erase(std::find_if(line.rbegin(), line.rend(),
                                [](unsigned char ch){ return !std::isspace(ch); }).base(),
                   line.end());

        if (!line.empty()) varnames.push_back(line);
    }

    if (varnames.empty()) {
        amrex::Abort("Variable file is empty: " + filename);
    }

    return varnames;
}

void
ReadPlotFile(const std::string& var_filename,
             PlotFileData& pf,
             amrex::MultiFab& mf)
{
    // ------------------------------------------------------------
    // Read variable list from file
    // ------------------------------------------------------------
    Vector<std::string> varnames;
    {
        std::ifstream infile(var_filename);
        if (!infile.is_open()) {
            Abort("ReadPlotFile: Cannot open variable file: " + var_filename);
        }

        std::string line;
        while (std::getline(infile, line)) {
            if (!line.empty()) {
                // trim whitespace
                line.erase(0, line.find_first_not_of(" \t\r\n"));
                line.erase(line.find_last_not_of(" \t\r\n") + 1);

                if (!line.empty()) varnames.push_back(line);
            }
        }
    }

    if (varnames.empty()) {
        Abort("ReadPlotFile: No variables found in " + var_filename);
    }

    // ------------------------------------------------------------
    // Open plotfile
    // ------------------------------------------------------------
    const Vector<std::string>& var_names_pf = pf.varNames();

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
            Abort("ReadPlotFile: invalid variable name: " + v);
        }
    }

    // ------------------------------------------------------------
    // Define destination MultiFab (single level only)
    // ------------------------------------------------------------
    const int level = 0;

    BoxArray ba = pf.boxArray(level);
    DistributionMapping dm(ba);

    int ncomp = varnames.size();
    int ngrow = 0;

    mf.define(ba, dm, ncomp, ngrow);

    // ------------------------------------------------------------
    // Copy plotfile data → mf
    // ------------------------------------------------------------
    for (int comp = 0; comp < ncomp; ++comp)
    {
        const MultiFab& src = pf.get(level, varnames[comp]);
        MultiFab::Copy(mf, src, 0, comp, 1, 0);
    }
}
