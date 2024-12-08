#ifndef ERF_FFT_UTILS_H_
#define ERF_FFT_UTILS_H_

#ifdef ERF_USE_FFT

using namespace amrex;

inline Array<std::pair<FFT::Boundary,FFT::Boundary>,AMREX_SPACEDIM>
get_fft_bc (Geometry const& lev_geom,
            Array<std::string,2*AMREX_SPACEDIM> l_domain_bc_type,
            Box const& bounding_box) noexcept
{
    Array<std::pair<FFT::Boundary,FFT::Boundary>,AMREX_SPACEDIM> r;

    for (int dir = 0; dir <= 1; dir++)
    {
        auto bc_type_lo = l_domain_bc_type[Orientation(dir,Orientation::low)];
        auto bc_type_hi = l_domain_bc_type[Orientation(dir,Orientation::high)];

        if ( lev_geom.isPeriodic(dir) &&
            (lev_geom.Domain().smallEnd(dir) == bounding_box.smallEnd(dir)) &&
            (lev_geom.Domain().bigEnd(dir)   == bounding_box.bigEnd(dir)) ) {
            r[dir] = std::make_pair(FFT::Boundary::periodic,FFT::Boundary::periodic);
             // amrex::Print() << "SETTING " << dir << " TO PERIODIC " << std::endl;
        }

        else if ( (lev_geom.Domain().smallEnd(dir) == bounding_box.smallEnd(dir)) &&
                  (lev_geom.Domain().bigEnd(dir)   == bounding_box.bigEnd(dir)) &&
                  (bc_type_lo != "Outflow" && bc_type_lo != "Open") &&
                  (bc_type_hi != "Outflow" && bc_type_hi != "Open") )
        {
            r[dir] = std::make_pair(FFT::Boundary::even,FFT::Boundary::even);
            // amrex::Print() << "SETTING " << dir << " TO EVEN EVEN " << std::endl;
        }

        else if ( (lev_geom.Domain().smallEnd(dir) == bounding_box.smallEnd(dir)) &&
                  (lev_geom.Domain().bigEnd(dir)   == bounding_box.bigEnd(dir)) &&
                  (bc_type_lo == "Outflow" || bc_type_lo == "Open") &&
                  (bc_type_hi == "Outflow" || bc_type_hi == "Open") )
        {
            r[dir] = std::make_pair(FFT::Boundary::odd,FFT::Boundary::odd);
            // amrex::Print() << "SETTING " << dir << " TO ODD ODD " << std::endl;
        }

        else if ( (lev_geom.Domain().smallEnd(dir)   == bounding_box.smallEnd(dir)) &&
                  (bc_type_lo == "Outflow" || bc_type_lo == "Open") )
        {
                r[dir] = std::make_pair(FFT::Boundary::odd,FFT::Boundary::even);
             // amrex::Print() << "SETTING " << dir << " TO ODD EVEN " << std::endl;
        }

        else if ( (lev_geom.Domain().bigEnd(dir)   == bounding_box.bigEnd(dir)) &&
                  (bc_type_hi == "Outflow" || bc_type_hi == "Open") )
        {
            r[dir] = std::make_pair(FFT::Boundary::even,FFT::Boundary::odd);
            // amrex::Print() << "SETTING " << dir << " TO EVEN ODD " << std::endl;
        }
        else
        {
            r[dir] = std::make_pair(FFT::Boundary::even,FFT::Boundary::even);
            // amrex::Print() << "SETTING " << dir << " TO EVEN EVEN " << std::endl;
        }
    } // dir

    //
    // Always Neumann in the vertical
    //
    r[2] = std::make_pair(FFT::Boundary::even,FFT::Boundary::even);

    return r;
}
#endif
#endif
