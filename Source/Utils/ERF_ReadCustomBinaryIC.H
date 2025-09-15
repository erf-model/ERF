#ifndef ERF_READCUSTOMBINARYIC_H_
#define ERF_READCUSTOMBINARYIC_H_

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
    int nx, ny, nz, ndata;
    float value;

    // Open the binary file in input mode
    std::ifstream infile(filename, std::ios::binary);
    if (!infile) {
        std::cerr << "Error: Could not open file " << filename << std::endl;
    }

    // Read the four integers
    infile.read(reinterpret_cast<char*>(&nx), sizeof(int));
    infile.read(reinterpret_cast<char*>(&ny), sizeof(int));
    infile.read(reinterpret_cast<char*>(&nz), sizeof(int));
    infile.read(reinterpret_cast<char*>(&ndata), sizeof(int));

    for(int i=0; i<nx*ny; i++) {
        infile.read(reinterpret_cast<char*>(&value), sizeof(float));
        latvec_h.emplace_back(value);
    }

    for(int i=0; i<nx*ny; i++) {
        infile.read(reinterpret_cast<char*>(&value), sizeof(float));
        lonvec_h.emplace_back(value);
    }

    for(int i=0; i<nx; i++) {
        infile.read(reinterpret_cast<char*>(&value), sizeof(float));
        xvec_h.emplace_back(value);
    }

    for(int j=0; j<ny; j++) {
        infile.read(reinterpret_cast<char*>(&value), sizeof(float));
        yvec_h.emplace_back(value);
    }

    for(int k=0; k<nz; k++) {
        infile.read(reinterpret_cast<char*>(&value), sizeof(float));
        zvec_h.emplace_back(value);
    }

    amrex::Vector<amrex::Real>* data_h = nullptr; // Declare pointer outside the loop

    // Read the file
    for(int idx=0; idx<ndata; idx++){
        if(idx == 0){
            data_h = &rho_h;
        } else if (idx==1) {
            data_h = &uvel_h;
        } else if (idx==2) {
            data_h = &vvel_h;
        } else if (idx==3) {
            data_h = &wvel_h;
        } else if(idx==4) {
            data_h = &theta_h;
        } else if(idx==5) {
            data_h = &qv_h;
        } else if(idx==6) {
            data_h = &qc_h;
        } else if(idx==7) {
            data_h = &qr_h;
        }
        for(int k=0; k<nz; k++) {
            for(int j=0; j<ny; j++) {
                for(int i=0; i<nx; i++) {
                    infile.read(reinterpret_cast<char*>(&value), sizeof(float));
                    data_h->emplace_back(value);
                }
            }
        }
    }
    infile.close();
}
#endif
