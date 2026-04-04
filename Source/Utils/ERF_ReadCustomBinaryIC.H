#ifndef ERF_READCUSTOMBINARYIC_H_
#define ERF_READCUSTOMBINARYIC_H_

#include <array>
#include <fstream>

#include "ERF_DataStruct.H"

/*
 * Routine to read in the custom written binary file into vectors
 * of variables
 */

inline void
ReadCustomBinaryIC(const std::string filename,
                   amrex::Vector<amrex::Real>& latvec_h,
                   amrex::Vector<amrex::Real>& lonvec_h,
                   amrex::Vector<amrex::Real>& xvec_h,
                   amrex::Vector<amrex::Real>& yvec_h,
                   amrex::Vector<amrex::Real>& zvec_h,
                   amrex::Vector<amrex::Real>& rho_h,
                   amrex::Vector<amrex::Real>& uvel_h,
                   amrex::Vector<amrex::Real>& vvel_h,
                   amrex::Vector<amrex::Real>& wvel_h,
                   amrex::Vector<amrex::Real>& theta_h,
                   amrex::Vector<amrex::Real>& qv_h,
                   amrex::Vector<amrex::Real>& qc_h,
                   amrex::Vector<amrex::Real>& qr_h)
{
    int nx = 0, ny = 0, nz = 0, ndata = 0;
    float value;

    // Open the binary file in input mode
    std::ifstream infile(filename, std::ios::binary);
    if (!infile) {
        amrex::Abort("Error: Could not open file " + filename);
    }

    auto read_or_abort = [&](auto& dest, const std::string& field_name) {
        infile.read(reinterpret_cast<char*>(&dest), sizeof(dest));
        if (!infile) {
            amrex::Abort("Error reading " + field_name + " from " + filename);
        }
    };

    auto read_float_vector = [&](amrex::Vector<amrex::Real>& dest,
                                 int count,
                                 const std::string& field_name)
    {
        if (count < 0) {
            amrex::Abort("Error: negative element count for " + field_name + " in " + filename);
        }

        dest.reserve(dest.size() + count);
        for (int idx = 0; idx < count; ++idx) {
            read_or_abort(value, field_name);
            dest.emplace_back(value);
        }
    };

    // Read the four integers
    read_or_abort(nx, "nx");
    read_or_abort(ny, "ny");
    read_or_abort(nz, "nz");
    read_or_abort(ndata, "ndata");

    if (nx <= 0 || ny <= 0 || nz <= 0) {
        amrex::Abort("Error: invalid custom binary dimensions in " + filename);
    }

    if (ndata < 0 || ndata > 8) {
        amrex::Abort("Error: unsupported number of data fields in " + filename);
    }

    read_float_vector(latvec_h, nx*ny, "latitude");
    read_float_vector(lonvec_h, nx*ny, "longitude");
    read_float_vector(xvec_h, nx, "x coordinates");
    read_float_vector(yvec_h, ny, "y coordinates");
    read_float_vector(zvec_h, nz, "z coordinates");

    std::array<amrex::Vector<amrex::Real>*, 8> data_fields = {
        &rho_h, &uvel_h, &vvel_h, &wvel_h, &theta_h, &qv_h, &qc_h, &qr_h
    };

    // Read the file
    for (int idx = 0; idx < ndata; ++idx) {
        read_float_vector(*data_fields[idx], nx * ny * nz,
                          "data field " + std::to_string(idx));
    }
}
#endif
