#ifndef ERF_READCUSTOMBINARYIC_H_
#define ERF_READCUSTOMBINARYIC_H_

#include <array>
#include <fstream>
#include <limits>

#include "ERF_DataStruct.H"

/*
 * Routine to read in the custom written binary file into vectors
 * of variables
 */

inline void
ReadCustomBinaryIC (const std::string filename,
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
                                 amrex::Long count,
                                 const std::string& field_name)
    {
        if (count < 0) {
            amrex::Abort("Error: negative element count for " + field_name + " in " + filename);
        }

        dest.reserve(dest.size() + count);
        for (amrex::Long idx = 0; idx < count; ++idx) {
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

    // The dimensions come from the file: form the element counts in 64 bits, refusing any
    // whose byte count would overflow, and check the file holds that many values before
    // reserving room for them
    constexpr amrex::Long max_count = std::numeric_limits<amrex::Long>::max()
                                    / static_cast<amrex::Long>(sizeof(float));
    const amrex::Long nxy = static_cast<amrex::Long>(nx) * ny;
    if (nz > max_count / nxy) {
        amrex::Abort("Error: custom binary dimensions too large in " + filename);
    }
    const amrex::Long nxyz = nxy * nz;

    // 2 nxy + nx + ny + nz + ndata nxyz values, without overflowing
    amrex::Long nvalues = 2*nxy + nx + ny + nz;
    if (nvalues > max_count || (ndata > 0 && nxyz > (max_count - nvalues) / ndata)) {
        amrex::Abort("Error: custom binary dimensions too large in " + filename);
    }
    nvalues += ndata * nxyz;

    const auto data_start = infile.tellg();
    infile.seekg(0, std::ios::end);
    const auto bytes_left = static_cast<amrex::Long>(infile.tellg() - data_start);
    infile.seekg(data_start);
    if (!infile || bytes_left < nvalues * static_cast<amrex::Long>(sizeof(float))) {
        amrex::Abort("Error: " + filename + " is shorter than its dimensions say");
    }

    read_float_vector(latvec_h, nxy, "latitude");
    read_float_vector(lonvec_h, nxy, "longitude");
    read_float_vector(xvec_h, nx, "x coordinates");
    read_float_vector(yvec_h, ny, "y coordinates");
    read_float_vector(zvec_h, nz, "z coordinates");

    std::array<amrex::Vector<amrex::Real>*, 8> data_fields = {
        &rho_h, &uvel_h, &vvel_h, &wvel_h, &theta_h, &qv_h, &qc_h, &qr_h
    };

    // Read the file
    for (int idx = 0; idx < ndata; ++idx) {
        read_float_vector(*data_fields[idx], nxyz,
                          "data field " + std::to_string(idx));
    }
}
#endif
