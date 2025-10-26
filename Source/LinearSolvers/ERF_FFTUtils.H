#ifndef ERF_FFT_UTILS_H_
#define ERF_FFT_UTILS_H_

#ifdef ERF_USE_FFT

using namespace amrex;

inline Array<std::pair<FFT::Boundary,FFT::Boundary>,AMREX_SPACEDIM>
get_fft_bc (Geometry const& lev_geom,
            Array<std::string,2*AMREX_SPACEDIM> l_domain_bc_type,
            Box const& bounding_box, bool use_real_bcs) noexcept
{
    Array<std::pair<FFT::Boundary,FFT::Boundary>,AMREX_SPACEDIM> r;

    FFT::Boundary first, second;

    for (int dir = 0; dir <= 2; dir++)
    {
        auto bc_type_lo = l_domain_bc_type[Orientation(dir,Orientation::low)];
        if ( lev_geom.isPeriodic(dir) &&
            (lev_geom.Domain().smallEnd(dir) == bounding_box.smallEnd(dir)) )
        {
            first = FFT::Boundary::periodic;
             // amrex::Print() << "SETTING FIRST OF " << dir << " TO PERIODIC " << std::endl;
        }

        else if ( (lev_geom.Domain().smallEnd(dir) == bounding_box.smallEnd(dir)) &&
                  (bc_type_lo == "Outflow" || bc_type_lo == "Open") && (dir == 2 || !use_real_bcs) )
        {
            first = FFT::Boundary::odd;
             // amrex::Print() << "SETTING FIRST OF " << dir << " TO ODD " << std::endl;
        }

        else
        {
            first = FFT::Boundary::even;
            // amrex::Print() << "SETTING FIRST OF " << dir << " TO EVEN " << std::endl;
        }

        auto bc_type_hi = l_domain_bc_type[Orientation(dir,Orientation::high)];
        if ( lev_geom.isPeriodic(dir) &&
            (lev_geom.Domain().bigEnd(dir) == bounding_box.bigEnd(dir)) )
        {
            second = FFT::Boundary::periodic;
             // amrex::Print() << "SETTING SECOND OF " << dir << " TO PERIODIC " << std::endl;
        }

        else if ( (lev_geom.Domain().bigEnd(dir) == bounding_box.bigEnd(dir)) &&
                  (bc_type_hi == "Outflow" || bc_type_hi == "Open") && (dir == 2 || !use_real_bcs) )
        {
            second = FFT::Boundary::odd;
             // amrex::Print() << "SETTING SECOND OF " << dir << " TO ODD " << std::endl;
        }

        else
        {
            second = FFT::Boundary::even;
            // amrex::Print() << "SETTING SECOND OF " << dir << " TO EVEN " << std::endl;
        }

        r[dir] = std::make_pair(first,second);

    } // dir

    return r;
}
#endif
#endif
