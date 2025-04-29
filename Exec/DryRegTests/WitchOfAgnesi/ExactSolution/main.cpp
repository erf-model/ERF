#include <iostream>
#include <fstream>
#include <vector>
#include <complex>

std::complex<double> i_cplx(0.0, 1.0);

void
compute_hhat_of_k(const std::vector<double>& h_of_x,
                       const std::vector<double>& kvec,
                       const std::vector<double>& xvec,
                       const double dx,
                       std::vector<std::complex<double>>& hhat_of_k)
{


    for (int n = 0; n < kvec.size(); ++n) { // n = index for k
        std::complex<double> sum = 0.0;
        for (int j = 0; j < xvec.size(); ++j) { // j = index for x
            sum += h_of_x[j] * std::exp(-i_cplx * kvec[n] * xvec[j]) * dx;
        }
        hhat_of_k[n] = sum / (2.0 * M_PI);
    }
}

void write_vtk_structured_grid(const std::string& filename,
                               const std::vector<double>& x,
                               const std::vector<double>& z,
                               const std::vector<double>& h_of_x,
                               const std::vector<std::vector<double>>& var)
{
    const int nx = x.size();
    const int nz = z.size();

    std::ofstream ofs(filename);
    if (!ofs.is_open()) {
        std::cerr << "Error: Could not open file " << filename << "\n";
        return;
    }

    // Write VTK header
    ofs << "# vtk DataFile Version 3.0\n";
    ofs << "2D Structured Grid (x,z)\n";
    ofs << "ASCII\n";
    ofs << "DATASET STRUCTURED_GRID\n";
    ofs << "DIMENSIONS " << nx << " " << nz << " " << 1 <<"\n";
    ofs << "POINTS " << nx * nz << " float\n";

    // Write all point coordinates
    // Loop in VTK order: fastest x, then y, then z
    for (int k = 0; k < nz; ++k) {
        for (int i = 0; i < nx; ++i) {
            ofs << x[i] << " " << h_of_x[i] + z[k] << " " << 0.0 << "\n";
        }
    }

    // Write variable
    ofs << "\nPOINT_DATA " << nx * nz<< "\n";
    ofs << "SCALARS z_velocity float 1\n";
    ofs << "LOOKUP_TABLE default\n";

    for (int k = 0; k < nz; ++k) {
        for (int i = 0; i < nx; ++i) {
            ofs << var[i][k] << "\n";
        }
    }
    ofs.close();
}

double
compute_eta(const std::vector<std::complex<double>>& hhat_of_k,
                 const std::vector<double>& kvec,
                 const double& xval,
                 const double& zval,
                 const double dk,
                 const double N_BV,
                 const double U)
{

    std::complex<double> tmp(0.0,0.0);
    for(int i=0;i<kvec.size(); i++) {
        std::complex<double> m(0.0,0.0);
        double m_real=0.0;
        if(std::fabs(kvec[i]) > N_BV/U) {
            m = i_cplx*std::sqrt(kvec[i]*kvec[i] - N_BV*N_BV/(U*U));
        } else {
            m_real = std::sqrt(N_BV*N_BV/(U*U) - kvec[i]*kvec[i]);
            m_real = std::copysign(std::abs(m_real), kvec[i]);
            m = std::complex<double>(m_real, 0.0);
        }
        tmp = tmp + hhat_of_k[i]*std::exp(i_cplx*(m*zval+kvec[i]*xval));
    }
    tmp = tmp*dk;
    return std::real(tmp);
}


int main() {

    double h, xc, N_BV, U, L;
    h = 1.0;
    N_BV = 0.01; //1/s
    U = 10.0; //m/s
    L = 1000.0;

    double xmin, xmax;
    xmin = 0.0;
    xmax = 144e3;
    xc = (xmin+xmax)/2.0;
    int nx = 512;
    double dx = (xmax-xmin)/nx;

    double kmax = M_PI/dx;

    std::cout << "Value of kmax is " << kmax << std::endl;

    std::vector<double> xvec(nx, 0.0);
    std::vector<double> kvec(nx, 0.0);
    double dk = 2.0*kmax/nx;
    for(int i=0;i<nx; i++){
        xvec[i] = xmin + i*dx;
        kvec[i] = -kmax + i*dk;
    }

    std::ifstream infile("zvec.txt");
    std::vector<double> zvec;
    double value;

    while (infile >> value) {
        zvec.push_back(value);
    }

    int nz = zvec.size();

    std::vector<double> h_of_x;
    for(int i=0;i<xvec.size();i++){
        double x_by_L = (xvec[i] - xc)/L;
        double tmp = 1.0/(1.0 + x_by_L*x_by_L);
        h_of_x.push_back(tmp);
    }

    std::vector<std::complex<double>> hhat_of_k(nx, std::complex<double>(0.0, 0.0));

    // Compute hhat_of_k
    compute_hhat_of_k(h_of_x, kvec, xvec, dx, hhat_of_k);

    // Compute eta(x,y,z)
    std::vector<std::vector<double>> eta(nx, std::vector<double>(nz, 0.0));

    for(int k=0; k<nz; k++) {
        double zval = zvec[k];
        std::cout << "Doing k " << k << std::endl;
        for(int i=0; i<nx; i++) {
            double xval = xvec[i];
            eta[i][k] = compute_eta(hhat_of_k, kvec, xval, zval, dk, N_BV, U);
        }
    }

    std::vector<std::vector<double>> detadx(nx, std::vector<double>(nz, 0.0));
    for(int k=0; k<nz; k++) {
        for(int i=0; i<nx; i++) {
            if(i==0){
                detadx[i][k] = (eta[i+1][k]-eta[0][k])/dx;
            } else if(i==nx-1) {
                detadx[i][k] = (eta[i][k]-eta[i-1][k])/dx;
            } else {
                detadx[i][k] = (eta[i+1][k]-eta[i-1][k])/(2.0*dx);
            }
            detadx[i][k] = U*detadx[i][k];
        }
    }

    write_vtk_structured_grid("WoA_zvel.vtk", xvec, zvec, h_of_x, detadx);

    return 0;
}
