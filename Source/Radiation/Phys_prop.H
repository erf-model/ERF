
// Properties of aerosols that are used by radiation and other parameterizations.

#ifndef ERF_PHYS_PROP_H_
#define ERF_PHYS_PROP_H_

#include "YAKL_netcdf.h"
#include "rrtmgp_const.h"
#include "Rad_constants.H"
#include <complex>
#include <cstring>
#include <fstream>

using yakl::fortran::parallel_for;
using yakl::fortran::SimpleBounds;

class PhysProp {
  public:
    // Data from one input dataset is stored in a structure of type(physprop_type).
    struct physprop_t {
      std::string sourcefile; // Absolute pathname of data file.
      std::string opticsmethod; // one of {hygro, nonhygro}

      // for hygroscopic species of externally mixed aerosols
      real2d sw_hygro_ext;
      real2d sw_hygro_ssa;
      real2d sw_hygro_asm;
      real2d lw_hygro_abs;

      // for nonhygroscopic species of externally mixed aerosols
      real1d sw_nonhygro_ext;
      real1d sw_nonhygro_ssa;
      real1d sw_nonhygro_asm;
      real1d sw_nonhygro_scat;
      real1d sw_nonhygro_ascat;
      real1d lw_abs;

      // complex refractive index
      real1d refindex_real_aer_sw;
      real1d refindex_im_aer_sw;
      real1d refindex_real_aer_lw;
      real1d refindex_im_aer_lw;

      // for radius-dependent mass-specific quantities
      real2d r_sw_ext;
      real2d r_sw_scat;
      real2d r_sw_ascat;
      real2d r_lw_abs;
      real1d mu;

      // for modal optics
      real4d extpsw;     // specific extinction
      real4d abspsw;     // specific absorption
      real4d asmpsw;     // asymmetry factor
      real4d absplw;     // specific absorption
      real2d refrtabsw;  // table of real refractive indices for aerosols visible
      real2d refitabsw;  // table of imag refractive indices for aerosols visible
      real2d refrtablw;  // table of real refractive indices for aerosols infrared
      real2d refitablw;  // table of imag refractive indices for aerosols infrared

      // microphysics parameters.
      std::string aername;   // for output of number concentration
      real density_aer;      // density of aerosol (kg/m3)
      real hygro_aer;        // hygroscopicity of aerosol
      real dryrad_aer;       // number mode radius (m) of aerosol size distribution
      real dispersion_aer;   // geometric standard deviation of aerosol size distribution
      real num_to_mass_aer;  // ratio of number concentration to mass concentration (#/kg)
                                   // *** Is this actually (kg/#) ???
      //mode parameters
      int  ncoef;       // number of Chebyshev coefficients
      int  prefr;       // dimension in table of real refractive indices
      int  prefi;       // dimension in table of imag refractive indices
      real sigmag;      // geometric standard deviation of the number distribution for aerosol mode
      real dgnum;       // geometric dry mean diameter of the number distribution for aerosol mode
      real dgnumlo;     // lower limit of dgnum
      real dgnumhi;     // upper limit of dgnum
      real rhcrystal;   // crystalization relative humidity for mode
      real rhdeliques;  // deliquescence relative humidity for mode
    };

    // This module stores data in an array of physprop_type structures.  The way this data
    // is accessed outside the module is via a physprop ID, which is an index into the array.
    std::vector<physprop_t> physprop;

    //Temporary storage location for filenames in namelist, and construction of dynamic index
    // to properties.  The unique filenames specified in the namelist are the identifiers of
    // the properties.  Searching the uniquefilenames array provides the index into the physprop
    // array.
    std::vector<std::string> uniquefilenames;

  public:
    void physprop_accum_unique_files(const std::string& filename,
                                     const std::string& type)
    {
      // Count number of aerosols in input radname array.  Aerosols are identified
      // as strings with a ".nc" suffix.
      // Construct a cumulative list of unique filenames containing physical property data.

      // check if filename is either a bulk aerosol or a mode
      if (type == "A" || type == "M") {
        // check if this filename has been used by another aerosol.  If not
        // then add it to the list of unique names.
        if (physprop_get_id(filename) < 0)
           uniquefilenames.push_back(filename);
      }
    }

    // Read properties from the aerosol data files.
    // ***N.B.*** The calls to physprop_accum_unique_files must be made before calling
    //            this init routine.  physprop_accum_unique_files is responsible for building
    //            the list of files to be read here.
    void physprop_init()
    {
      int numphysprops = uniquefilenames.size();
      physprop.resize(numphysprops);

      for(auto fileindex = 0; fileindex < numphysprops; ++fileindex) {
          physprop[fileindex].sourcefile = uniquefilenames[fileindex];
          aerosol_optics_init(physprop[fileindex]);
      }
    }

    // Look for filename in the global list of unique filenames (module data uniquefilenames).
    // If found, return it's index in the list.  Otherwise return -1.
    int physprop_get_id(std::string filename)
    {
       auto physprop_id = -1;
       auto numphysprops = uniquefilenames.size();
       for(auto iphysprop = 0; iphysprop < numphysprops; ++iphysprop) {
         if(uniquefilenames[iphysprop] == filename) { 
           physprop_id = iphysprop;
           break;
         }
       }
      return physprop_id;
   }

   void get_sourcefile(int& id, std::string& sourcefile)
   {
     if (id <= 0 || id > physprop.size()) printf("illegal ID value %d\n", id);
     sourcefile = physprop[id].sourcefile;
   }

   void get_opticstype(int& id, std::string& opticstype)
   {
     if (id <= 0 || id > physprop.size()) printf("illegal ID value %d\n", id);
     opticstype  = physprop[id].opticsmethod;
   }

   void get_sw_hygro_ext(int& id, real2d& sw_hygro_ext)
   {
     if (id <= 0 || id > physprop.size()) printf("illegal ID value %d\n", id);
     sw_hygro_ext = physprop[id].sw_hygro_ext;
   }

   void get_sw_hygro_ssa(int& id, real2d& sw_hygro_ssa)
   {
     if (id <= 0 || id > physprop.size()) printf("illegal ID value %d\n", id);
     sw_hygro_ssa = physprop[id].sw_hygro_ssa;
   }

   void get_sw_hygro_asm(int& id, real2d& sw_hygro_asm)
   {
     if (id <= 0 || id > physprop.size()) printf("illegal ID value %d\n", id);
     sw_hygro_asm = physprop[id].sw_hygro_asm;
   }

   void get_lw_hygro_abs(int& id, real2d& lw_hygro_abs)
   {
     if (id <= 0 || id > physprop.size()) printf("illegal ID value %d\n", id);
     lw_hygro_abs = physprop[id].lw_hygro_abs;
   }

   void get_sw_nonhygro_ext(int& id, real1d& sw_nonhygro_ext)
   {
     if (id <= 0 || id > physprop.size()) printf("illegal ID value %d\n", id);
     sw_nonhygro_ext = physprop[id].sw_nonhygro_ext;
   }

   void get_sw_nonhygro_ssa(int& id, real1d& sw_nonhygro_ssa)
   {
     if (id <= 0 || id > physprop.size()) printf("illegal ID value %d\n", id);
     sw_nonhygro_ssa = physprop[id].sw_nonhygro_ssa;
   }

   void get_sw_nonhygro_asm(int& id, real1d& sw_nonhygro_asm)
   {
     if (id <= 0 || id > physprop.size()) printf("illegal ID value %d\n", id);
     sw_nonhygro_asm = physprop[id].sw_nonhygro_asm;
   }

   void get_sw_nonhygro_scat(int& id, real1d& sw_nonhygro_scat)
   {
     if (id <= 0 || id > physprop.size()) printf("illegal ID value %d\n", id);
     sw_nonhygro_scat = physprop[id].sw_nonhygro_scat;
   }

   void get_sw_nonhygro_ascat(int& id, real1d& sw_nonhygro_ascat)
   {
     if (id <= 0 || id > physprop.size()) printf("illegal ID value %d\n", id);
     sw_nonhygro_ascat = physprop[id].sw_nonhygro_ascat;
   }

   void get_lw_abs(int& id, real1d& lw_abs)
   {
     if (id <= 0 || id > physprop.size()) printf("illegal ID value %d\n", id);
     lw_abs = physprop[id].lw_abs;
   }

   void get_ref_real_aer_sw(int& id, real1d& ref_real_aer_sw)
   {
     if (id <= 0 || id > physprop.size()) printf("illegal ID value %d\n", id);
     ref_real_aer_sw = physprop[id].refindex_real_aer_sw;
   }

   void get_ref_real_aer_lw(int& id, real1d& ref_real_aer_lw)
   {
     if (id <= 0 || id > physprop.size()) printf("illegal ID value %d\n", id);
     ref_real_aer_lw = physprop[id].refindex_real_aer_lw;
   }

   void get_ref_im_aer_sw(int& id, real1d& ref_im_aer_sw)
   {
     if (id <= 0 || id > physprop.size()) printf("illegal ID value %d\n", id);
     ref_im_aer_sw = physprop[id].refindex_im_aer_sw;
   }

   void get_ref_im_aer_lw(int& id, real1d& ref_im_aer_lw)
   {
     if (id <= 0 || id > physprop.size()) printf("illegal ID value %d\n", id);
     ref_im_aer_lw = physprop[id].refindex_im_aer_lw;
   }

   void get_r_sw_ext(int& id, real2d& r_sw_ext)
   {
     if (id <= 0 || id > physprop.size()) printf("illegal ID value %d\n", id);
     r_sw_ext = physprop[id].r_sw_ext;
   }

   void get_r_sw_scat(int& id, real2d& r_sw_scat)
   {
     if (id <= 0 || id > physprop.size()) printf("illegal ID value %d\n", id);
     r_sw_scat = physprop[id].r_sw_scat;
   }

   void get_r_sw_ascat(int& id, real2d& r_sw_ascat)
   {
     if (id <= 0 || id > physprop.size()) printf("illegal ID value %d\n", id);
     r_sw_ascat = physprop[id].r_sw_ascat;
   }

   void get_r_lw_abs(int& id, real2d& r_lw_abs)
   {
     if (id <= 0 || id > physprop.size()) printf("illegal ID value %d\n", id);
     r_lw_abs = physprop[id].r_lw_abs;
   }

   void get_mu(int& id, real1d& mu)
   {
     if (id <= 0 || id > physprop.size()) printf("illegal ID value %d\n", id);
     mu = physprop[id].mu;
   }

   void get_extpsw(int& id, real4d& extpsw)
   {
     if (id <= 0 || id > physprop.size()) printf("illegal ID value %d\n", id);
     extpsw = physprop[id].extpsw;
   }

   void get_abspsw(int& id, real4d& abspsw)
   {
     if (id <= 0 || id > physprop.size()) printf("illegal ID value %d\n", id);
     abspsw = physprop[id].abspsw;
   }

   void get_asmpsw(int& id, real4d& asmpsw)
   {
     if (id <= 0 || id > physprop.size()) printf("illegal ID value %d\n", id);
     asmpsw = physprop[id].asmpsw;
   }

   void get_absplw(int& id, real4d& absplw)
   {
     if (id <= 0 || id > physprop.size()) printf("illegal ID value %d\n", id);
     absplw = physprop[id].absplw;
   }

   void get_refrtabsw(int& id, real2d& refrtabsw)
   {
     if (id <= 0 || id > physprop.size()) printf("illegal ID value %d\n", id);
     refrtabsw = physprop[id].refrtabsw;
   }

   void get_refitabsw(int& id, real2d& refitabsw)
   {
     if (id <= 0 || id > physprop.size()) printf("illegal ID value %d\n", id);
     refitabsw = physprop[id].refitabsw;
   }

   void get_refrtablw(int& id, real2d& refrtablw)
   {
     if (id <= 0 || id > physprop.size()) printf("illegal ID value %d\n", id);
     refrtablw = physprop[id].refrtablw;
   }

   void get_refitablw(int& id, real2d& refitablw)
   {
     if (id <= 0 || id > physprop.size()) printf("illegal ID value %d\n", id);
     refitablw = physprop[id].refitablw;
   }

   void get_aername(int& id, std::string& aername)
   {
     if (id <= 0 || id > physprop.size()) printf("illegal ID value %d\n", id);
     aername = physprop[id].aername;
   }

   void get_density_aer(int& id, real& density_aer)
   {
     if (id <= 0 || id > physprop.size()) printf("illegal ID value %d\n", id);
     density_aer = physprop[id].density_aer;
   }

   void get_hygro_aer(int& id, real& hygro_aer)
   {
     if (id <= 0 || id > physprop.size()) printf("illegal ID value %d\n", id);
     hygro_aer = physprop[id].hygro_aer;
   }

   void get_dryrad_aer(int& id, real& dryrad_aer)
   {
     if (id <= 0 || id > physprop.size()) printf("illegal ID value %d\n", id);
     dryrad_aer = physprop[id].dryrad_aer;
   }

   void get_dispersion_aer(int& id, real& dispersion_aer)
   {
     if (id <= 0 || id > physprop.size()) printf("illegal ID value %d\n", id);
     dispersion_aer = physprop[id].dispersion_aer;
   }

   void get_num_to_mass_aer(int& id, real& num_to_mass_aer)
   {
     if (id <= 0 || id > physprop.size()) printf("illegal ID value %d\n", id);
     num_to_mass_aer = physprop[id].num_to_mass_aer;
   }

   void get_ncoef(int& id, int& ncoef)
   {
     if (id <= 0 || id > physprop.size()) printf("illegal ID value %d\n", id);
     ncoef = physprop[id].ncoef;
   }

   void get_prefr(int& id, int& prefr)
   {
     if (id <= 0 || id > physprop.size()) printf("illegal ID value %d\n", id);
     prefr = physprop[id].prefr;
   }

   void get_prefi(int& id, int& prefi)
   {
     if (id <= 0 || id > physprop.size()) printf("illegal ID value %d\n", id);
     prefi = physprop[id].prefi;
   }

   void get_sigmag(int& id, real& sigmag)
   {
     if (id <= 0 || id > physprop.size()) printf("illegal ID value %d\n", id);
     sigmag = physprop[id].sigmag;
   }

   void get_dgnum(int& id, real& dgnum)
   {
     if (id <= 0 || id > physprop.size()) printf("illegal ID value %d\n", id);
     dgnum = physprop[id].dgnum;
   }

   void get_dgnumlo(int& id, real& dgnumlo)
   {
     if (id <= 0 || id > physprop.size()) printf("illegal ID value %d\n", id);
     dgnumlo = physprop[id].dgnumlo;
   }

   void get_dgnumhi(int& id, real& dgnumhi)
   {
     if (id <= 0 || id > physprop.size()) printf("illegal ID value %d\n", id);
     dgnumhi = physprop[id].dgnumhi;
   }

   void get_rhcrystal(int& id, real& rhcrystal)
   {
     if (id <= 0 || id > physprop.size()) printf("illegal ID value %d\n", id);
     rhcrystal = physprop[id].rhcrystal;
   }

   void get_rhdeliques(int& id, real& rhdeliques)
   {
     if (id <= 0 || id > physprop.size()) printf("illegal ID value %d\n", id);
     rhdeliques = physprop[id].rhdeliques;
   }

   // Determine the opticstype, then call the
   // appropriate routine to read the data.
   void aerosol_optics_init(physprop_t& phys_prop)
   {
      using charHost1d = FArray<char,1,yakl::memHost>;
      yakl::SimpleNetCDF prop;
      prop.open(phys_prop.sourcefile, yakl::NETCDF_MODE_READ);
      charHost1d temp;
      prop.read(temp, "opticsmethod");
      for (auto ichar = 1 ; ichar <= temp.extent(0); ichar++)
         if (!isspace(temp(ichar))) phys_prop.opticsmethod += temp(ichar);

      if(strcmp(phys_prop.opticsmethod.c_str(),"zero") == 0) {
           zero_optics_init(phys_prop, prop);
      } else if (strcmp(phys_prop.opticsmethod.c_str(),"hygro") == 0) {
           hygro_optics_init(phys_prop, prop);
      } else if (strcmp(phys_prop.opticsmethod.c_str(),"hygroscopic") == 0) {
           hygroscopic_optics_init(phys_prop, prop);
      } else if (strcmp(phys_prop.opticsmethod.c_str(),"nonhygro") == 0) {
          nonhygro_optics_init(phys_prop, prop);
      } else if (strcmp(phys_prop.opticsmethod.c_str(),"insoluble") == 0) {
          insoluble_optics_init(phys_prop, prop);
      } else if (strcmp(phys_prop.opticsmethod.c_str(),"volcanic_radius") == 0) {
          volcanic_radius_optics_init(phys_prop, prop);
      } else if (strcmp(phys_prop.opticsmethod.c_str(),"volcanic") == 0) {
          volcanic_optics_init(phys_prop, prop);
      } else if (strcmp(phys_prop.opticsmethod.c_str(),"modal") == 0) {
          modal_optics_init(phys_prop, prop);
      } else {
        amrex::Print() << "no options available\n";
      }
    }

    void hygro_optics_init(physprop_t& phys_prop, yakl::SimpleNetCDF& prop)
    {
      real1d frh;
      real2d fsw_ext;
      real2d fsw_ssa;
      real2d fsw_asm;
      real2d flw_abs;

      auto nrh       = RadConstants::nrh;
      auto nbnd      = prop.getDimSize( "lw_band" );
      auto nswbands  = prop.getDimSize( "sw_band" );

      prop.read( fsw_ext, "ext_sw");
      prop.read( fsw_ssa, "ssa_sw");
      prop.read( fsw_asm, "asm_sw");
      prop.read( flw_abs, "abs_lw");
      prop.read( frh,         "rh");

      real1d fswe("fswe",nrh),
             fsws("fsws",nrh),
             fswa("fswa",nrh);

     // interpolate onto cam's rh mesh
     for(auto kbnd = 0; kbnd < nswbands; ++kbnd) {
        for(auto krh = 0; krh < nrh; ++krh) {
          fswe(krh) = fsw_ext(krh, kbnd)/fsw_ext(1, kbnd);
          fsws(krh) = fsw_ssa(krh, kbnd)/fsw_ssa(1, kbnd);
          fswa(krh) = fsw_asm(krh, kbnd)/fsw_asm(1, kbnd);
        }

        // interpolation
        for(auto krh = 0; krh < nrh; ++krh) {
           auto rh = 1.0 / nrh * (krh - 1);
           phys_prop.sw_hygro_ext(krh, kbnd) =
              exp_interpol( frh, fswe, rh ) * fsw_ext(1, kbnd);
           phys_prop.sw_hygro_ssa(krh,kbnd) =
              lin_interpol( frh, fsws, rh ) * fsw_ssa(1, kbnd);
           phys_prop.sw_hygro_asm(krh,kbnd) =
              lin_interpol( frh, fswa, rh ) * fsw_asm(1, kbnd);
        }
     }

     // read refractive index data if available
     refindex_aer_init(phys_prop, prop);

     // read bulk aero props
     bulk_props_init(phys_prop, prop);
   }

   // Read optics data of type 'nonhygro'
   void zero_optics_init(physprop_t& phys_prop, yakl::SimpleNetCDF& prop)
   {
      auto nlwbands  = prop.getDimSize( "lw_band" );
      auto nswbands  = prop.getDimSize( "sw_band" );

      // perhaps this doesn't even need allocated.
      phys_prop.sw_nonhygro_ext = real1d("sw_nonhygro_ext", nswbands);
      phys_prop.sw_nonhygro_ssa = real1d("sw_nonhygro_ssa", nswbands);
      phys_prop.sw_nonhygro_asm = real1d("sw_nonhygro_asm", nswbands);
      phys_prop.lw_abs          = real1d("lwabs", nlwbands);

      yakl::memset(phys_prop.sw_nonhygro_ext, 0.);
      yakl::memset(phys_prop.sw_nonhygro_ssa, 0.);
      yakl::memset(phys_prop.sw_nonhygro_asm, 0.);
      yakl::memset(phys_prop.lw_abs, 0.);
   }

   void insoluble_optics_init(physprop_t& phys_prop, yakl::SimpleNetCDF& prop)
   {
     auto nbnd      = prop.getDimSize( "lw_band" );
     auto nswbands  = prop.getDimSize( "sw_band" );
     realHost2d ext_sw, ssa_sw, asm_sw, abs_lw;

     phys_prop.sw_nonhygro_ext = real1d("sw_nonhygro_ext", nswbands);
     phys_prop.sw_nonhygro_ssa = real1d("sw_nonhygro_ssa", nswbands);
     phys_prop.sw_nonhygro_asm = real1d("sw_nonhygro_asm", nswbands);
     phys_prop.lw_abs          = real1d("lw_abs", nbnd);

     prop.read( ext_sw, "ext_sw");
     prop.read( ssa_sw, "ssa_sw");
     prop.read( asm_sw, "asm_sw");
     prop.read( abs_lw, "abs_lw");

     parallel_for (SimpleBounds<1>(nswbands), YAKL_LAMBDA (int i) {
       phys_prop.sw_nonhygro_ext(i) = ext_sw(i,1);
       phys_prop.sw_nonhygro_ssa(i) = ssa_sw(i,1);
       phys_prop.sw_nonhygro_asm(i) = asm_sw(i,1);
     });

     parallel_for (SimpleBounds<1>(nbnd), YAKL_LAMBDA (int i) {
       phys_prop.lw_abs(i) = abs_lw(i,1);
     });

     // read refractive index data if available
     refindex_aer_init(phys_prop, prop);

     // read bulk aero props
     bulk_props_init(phys_prop, prop);
   }

   // Read optics data of type 'volcanic_radius'
   void volcanic_radius_optics_init(physprop_t& phys_prop, yakl::SimpleNetCDF& prop)
   {
     auto n_mu_samples = prop.getDimSize( "mu_samples" );
     auto nbnd         = prop.getDimSize( "lw_band" );
     auto nswbands     = prop.getDimSize( "sw_band" );

     prop.read( phys_prop.r_sw_ext,   "bext_sw");
     prop.read( phys_prop.r_sw_scat,  "bsca_sw");
     prop.read( phys_prop.r_sw_ascat, "basc_sw");
     prop.read( phys_prop.r_lw_abs,   "babs_lw");
     prop.read( phys_prop.mu,         "mu_samples");

     // read bulk aero props
     bulk_props_init(phys_prop, prop);
  }

  // Read optics data of type 'volcanic'
  void volcanic_optics_init(physprop_t& phys_prop, yakl::SimpleNetCDF& prop)
  {
    auto nbnd     = prop.getDimSize( "lw_band" );
    auto nswbands = prop.getDimSize( "sw_band" );

    prop.read( phys_prop.sw_nonhygro_ext,   "bext_sw");
    prop.read( phys_prop.sw_nonhygro_scat,  "bsca_sw");
    prop.read( phys_prop.sw_nonhygro_ascat, "basc_sw");
    prop.read( phys_prop.lw_abs,            "babs_lw");

    // read bulk aero props
    bulk_props_init(phys_prop, prop);
  }

  // Read optics data of type 'hygroscopic' and interpolate it to CAM's rh mesh.
  void hygroscopic_optics_init(physprop_t& phys_prop, yakl::SimpleNetCDF& prop)
  {
    // temp data from hygroscopic file before interpolation onto cam-rh-mesh
    int nfilerh; // number of rh values in file
    real1d frh;
    real2d fsw_ext, fsw_ssa, fsw_asm, flw_abs;

    auto nrh      = RadConstants::nrh;
    auto nlwbands = prop.getDimSize( "lw_band" );
    auto nswbands = prop.getDimSize( "sw_band" );

    phys_prop.sw_hygro_ext = real2d("sw_hygro_ext",nrh,nswbands);
    phys_prop.sw_hygro_ssa = real2d("sw_hygro_ssa",nrh,nswbands);
    phys_prop.sw_hygro_asm = real2d("sw_hygro_asm",nrh,nswbands);
    phys_prop.lw_hygro_abs = real2d("lw_hygro_abs",nrh,nlwbands);

    prop.read( fsw_ext, "ext_sw");
    prop.read( fsw_ssa, "ssa_sw");
    prop.read( fsw_asm, "asm_sw");
    prop.read( flw_abs, "abs_lw");
    prop.read( frh,     "rh");

    real1d fswe("",nrh), fsws("",nrh),
           fswa("",nrh), flwa("",nrh);

    // interpolate onto cam's rh mesh
    parallel_for (SimpleBounds<2> (nswbands, nrh), YAKL_LAMBDA (int kbnd, int krh) {
      fswe(krh) = fsw_ext(krh,kbnd) / fsw_ext(1,kbnd);
      fsws(krh) = fsw_ssa(krh,kbnd) / fsw_ssa(1,kbnd);
      fswa(krh) = fsw_asm(krh,kbnd) / fsw_asm(1,kbnd);

      auto rh = 1.0 / nrh * (krh - 1);
      phys_prop.sw_hygro_ext(krh,kbnd) = exp_interpol( frh, fswe, rh ) * fsw_ext(1, kbnd);
      phys_prop.sw_hygro_ssa(krh,kbnd) = lin_interpol( frh, fswe, rh ) * fsw_ssa(1, kbnd);
      phys_prop.sw_hygro_asm(krh,kbnd) = lin_interpol( frh, fswa, rh ) * fsw_asm(1, kbnd);
    });

    parallel_for (SimpleBounds<2> (nswbands, nrh), YAKL_LAMBDA (int kbnd, int krh) {
      flwa(krh) = flw_abs(krh,kbnd) / flw_abs(1,kbnd);
      auto rh = 1.0 / nrh * (krh - 1);
      phys_prop.lw_hygro_abs(krh,kbnd) = exp_interpol( frh, flwa, rh ) * flw_abs(1, kbnd);
    });

    // read refractive index data if available
    refindex_aer_init(phys_prop, prop);

    // bulk aero props
    bulk_props_init(phys_prop, prop);
  }

  // Read optics data of type 'nonhygro'
  void nonhygro_optics_init(physprop_t& phys_prop, yakl::SimpleNetCDF& prop)
  {
    auto nlwbands = prop.getDimSize( "lw_band" );
    auto nswbands = prop.getDimSize( "sw_band" );

    prop.read( phys_prop.sw_nonhygro_ext, "ext_sw");
    prop.read( phys_prop.sw_nonhygro_ssa, "ssa_sw");
    prop.read( phys_prop.sw_nonhygro_asm, "asm_sw");
    prop.read( phys_prop.lw_abs,          "abs_lw");

    // read refractive index data if available
    refindex_aer_init(phys_prop, prop);

    // read bulk aero props
    bulk_props_init(phys_prop, prop);
 }

 // Read refractive indices of aerosol
 void refindex_aer_init(physprop_t& phys_prop, yakl::SimpleNetCDF& prop)
 {
   prop.read(phys_prop.refindex_real_aer_sw, "refindex_real_aer_sw");
   prop.read(phys_prop.refindex_im_aer_sw, "refindex_im_aer_sw");

   prop.read(phys_prop.refindex_real_aer_lw, "refindex_real_aer_lw");
   prop.read(phys_prop.refindex_im_aer_lw, "refindex_im_aer_lw");
 }

 // Read optics data for modal aerosols
 void modal_optics_init(physprop_t& phys_prop, yakl::SimpleNetCDF& prop)
 {
     realHost5d extpsw, abspsw, asmpsw, absplw;
     auto nlwbnd = prop.getDimSize( "lw_band" );
     auto nswbnd = prop.getDimSize( "sw_band" );
     auto ncoef  = prop.getDimSize( "coef_number" );
     auto prefr  = prop.getDimSize( "refindex_real" );
     auto prefi  = prop.getDimSize( "refindex_im" );
     phys_prop.ncoef = ncoef;
     phys_prop.prefr = prefr;
     phys_prop.prefi = prefi;

     prop.read(extpsw, "extpsw" );
     prop.read(abspsw, "abspsw" );
     prop.read(asmpsw, "asmpsw" );
     prop.read(absplw, "absplw" );

     phys_prop.extpsw = real4d("extpsw", nswbnd, prefi, prefr, ncoef);
     phys_prop.abspsw = real4d("abspsw", nswbnd, prefi, prefr, ncoef);
     phys_prop.asmpsw = real4d("asmpsw", nswbnd, prefi, prefr, ncoef);
     phys_prop.absplw = real4d("absplw", nlwbnd, prefi, prefr, ncoef);

     parallel_for (SimpleBounds<4>(nswbnd, prefi, prefr, ncoef), 
                   YAKL_LAMBDA (int i, int j, int k, int l) {
        phys_prop.extpsw(i,j,k,l) = extpsw(i,1,j,k,l);
        phys_prop.abspsw(i,j,k,l) = abspsw(i,1,j,k,l);
        phys_prop.asmpsw(i,j,k,l) = asmpsw(i,1,j,k,l);
     });

     parallel_for (SimpleBounds<4>(nlwbnd, prefi, prefr, ncoef), 
                   YAKL_LAMBDA (int i, int j, int k, int l) {
        phys_prop.absplw(i,j,k,l) = absplw(i,1,j,k,l);
     });

     prop.read(phys_prop.refrtabsw,   "refindex_real_sw" );
     prop.read(phys_prop.refitabsw,   "refindex_im_sw" );
     prop.read(phys_prop.refrtablw,   "refindex_real_lw" );
     prop.read(phys_prop.refitablw,   "refindex_im_lw" );

     prop.read(phys_prop.sigmag,      "sigmag" );
     prop.read(phys_prop.dgnum,       "dgnum" );

     prop.read(phys_prop.dgnumlo,     "dgnumlo" );
     prop.read(phys_prop.dgnumhi,     "dgnumhi" );

     prop.read(phys_prop.rhcrystal,   "rhcrystal" );
     prop.read(phys_prop.rhdeliques,  "rhdeliques" );
 }

 void bulk_props_init(physprop_t& phys_prop, yakl::SimpleNetCDF& prop)
 {
     using charHost1d = FArray<char,1,yakl::memHost>;
     charHost1d temp;
     prop.read(temp, "name");
     phys_prop.aername = "";
     for (auto ichar = 1 ; ichar <= temp.extent(0); ichar++)
        if (!isspace(temp(ichar))) phys_prop.aername += temp(ichar);

     //  Read props for bulk aerosols
     prop.read( phys_prop.density_aer,       "density" );
     prop.read( phys_prop.dispersion_aer,    "sigma_logr" );
     prop.read( phys_prop.dryrad_aer,        "dryrad" );
     prop.read( phys_prop.hygro_aer,         "hygroscopicity" );
     prop.read( phys_prop.num_to_mass_aer,   "num_to_mass_ratio" );
 }

 // Purpose:
 //interpolates f(x) to point y
 //assuming f(x) = f(x0) exp a(x - x0)
 //where a = ( ln f(x1) - ln f(x0) ) / (x1 - x0)
 //x0 <= x <= x1
 //assumes x is monotonically increasing
 real exp_interpol(const real1d& x, const real1d& f, const real& y)
 {
   auto n = x.extent(0);
   int k;

   // find k such that x(k) < y =< x(k+1)
   // set k = 1 if y <= x(1)  and  k = n-1 if y > x(n)
   if (y <= x(1)) {
      k = 1;
   }
   else if (y >= x(n)) {
     k = n - 1;
   }
   else {
     k = 1;
     while (y > x(k+1) && k < n) {
       k = k + 1;
     }
   }

   // interpolate
   auto a = (  log( f(k+1) / f(k) )  ) / ( x(k+1) - x(k) );
   return f(k) * exp( a * (y - x(k)) );
 }

 //Purpose:
 // interpolates f(x) to point y
 // assuming f(x) = f(x0) + a * (x - x0)
 // where a = ( f(x1) - f(x0) ) / (x1 - x0)
 // x0 <= x <= x1
 // assumes x is monotonically increasing
 //
 real lin_interpol(const real1d& x, const real1d& f, const real& y)
 {
   auto n = x.extent(0);
   int k;
   // find k such that x(k) < y =< x(k+1)
   // set k = 1 if y <= x(1)  and  k = n-1 if y > x(n)
   if (y <= x(1)) {
     k = 1;
   }
   else if (y >= x(n)) {
     k = n - 1;
   }
   else {
     k = 1;
     while (y > x(k+1) && k < n) {
       k = k + 1;
     }
   }
   // interpolate
   auto a = (f(k+1)-f(k))/(x(k+1)-x(k));
   return f(k)+a*(y-x(k));
 }

 // Purpose:
 //   write out aerosol optical properties
 //   for a set of test rh values
 //   to test hygroscopic growth interpolation
 void aer_optics_log_rh(std::string name, const real1d& ext, const real1d& ssa, const real1d& asmin)
 {
   const int nrh_test = 36;
   int krh;
   real1d rh_test("rh_test", nrh_test);
   auto nrh = ext.extent(0);

   parallel_for (SimpleBounds<1> (nrh_test), YAKL_LAMBDA (int krh_test) {
     rh_test(krh_test) = sqrt(sqrt(sqrt(sqrt(((krh_test - 1.0) / (nrh_test - 1))))));
   });

   // loop through test rh values
   parallel_for (SimpleBounds<1> (nrh_test), YAKL_LAMBDA (int krh_test) {
     // find corresponding rh index
     auto rh = rh_test(krh_test);
     auto krh = std::min(floor( (rh) * nrh ) + 1, static_cast<real>(nrh - 1));
     auto wrh = (rh) *nrh - krh;
     auto exti = ext(krh + 1) * (wrh + 1) - ext(krh) * wrh;
     auto ssai = ssa(krh + 1) * (wrh + 1) - ssa(krh) * wrh;
     auto asmi = asmin(krh + 1) * (wrh + 1) - asmin(krh) * wrh;
   });
 }
};

#endif
