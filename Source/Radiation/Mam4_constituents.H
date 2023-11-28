//
// Provide constituent distributions and properties to the radiation and
// cloud microphysics routines.
//
// The logic to control which constituents are used in the climate calculations
// and which are used in diagnostic radiation calculations is contained in this module.
//
#include "YAKL_netcdf.h"
#include "rrtmgp_const.h"
#include "Phys_prop.H"
#include "Rad_constants.H"

#ifndef ERF_MAM4_CONSTITUENTS_H_
#define ERF_MAM4_CONSTITUENTS_H_

class MamConstituents {
  public:
    const int N_DIAG = 10;
    std::string iceopticsfile, liqopticsfile;
    std::string icecldoptics, liqcldoptics;
    const bool oldcldoptics = false;

    // max number of externally mixed entities in the climate/diag lists
    const int n_rad_cnst = 30; //N_RAD_CNST (hard coded here)

    // type to provide access to the components of a mode
    struct mode_component_t {
      int nspec;
      // For "source" variables below, value is:
      std::string source_num_a;               // source of interstitial number conc field
      std::string camname_num_a;              // name registered in pbuf or constituents for number mixing ratio of interstitial species
      std::string source_num_c;               // source of cloud borne number conc field
      std::string camname_num_c;              //  name registered in pbuf or constituents for number mixing ratio of cloud borne species
      std::vector<std::string> source_mmr_a;  // source of interstitial specie mmr fields
      std::vector<std::string> camname_mmr_a; // name registered in pbuf or constituents for mmr of interstitial components
      std::vector<std::string> source_mmr_c;  // source of cloud borne specie mmr fields
      std::vector<std::string> camname_mmr_c; // name registered in pbuf or constituents for mmr of cloud borne components
      std::vector<std::string> type;          // specie type (as used in MAM code)
      std::vector<std::string> props;         // file containing specie properties
      int idx_num_a;                          // index in pbuf or constituents for number mixing ratio of interstitial species
      int idx_num_c;                          // index in pbuf for number mixing ratio of interstitial species
      std::vector<int>  idx_mmr_a;            // index in pbuf or constituents for mmr of interstitial species
      std::vector<int>  idx_mmr_c;            // index in pbuf for mmr of interstitial species
      std::vector<int>  idx_props;            // ID used to access physical properties of mode species from phys_prop module
    };

    // type to provide access to all modes
    struct modes_t {
      int nmodes;
      std::vector<std::string> names;      // names used to identify a mode in the climate/diag lists
      std::vector<std::string> types;      // type of mode (as used in MAM code)
      std::vector<mode_component_t> comps; // components which define the mode
    };

    //Storage for gas components in the climate/diagnostic lists
    struct gas_t {
      std::string source;       // A for state (advected), N for pbuf (non-advected), Z for zero
      std::string camname;      // name of constituent in physics state or buffer
      std::string mass_name;    // name for mass per layer field in history output
      int idx;                  // index from constituents or from pbuf
    };

    struct  gaslist_t {
      int  ngas;
      std::string list_id;
      std::vector<gas_t> gas;  // dimension(ngas) where ngas = nradgas is from radconstants
    };

    //Storage for bulk aerosol components in the climate/diagnostic lists
    struct aerosol_t {
      std::string source;         // A for state (advected), N for pbuf (non-advected), Z for zero
      std::string camname;        // name of constituent in physics state or buffer
      std::string physprop_file;  // physprop filename
      std::string mass_name;      // name for mass per layer field in history output
      int idx;                   // index of constituent in physics state or buffer
      int physprop_id;           // ID used to access physical properties from phys_prop module
    };

   struct  aerlist_t {
     int  numaerosols;           // number of aerosols
     std::string list_id;
     std::vector<aerosol_t> aer; // dimension(numaerosols)
   };

   // storage for modal aerosol components in the climate/diagnostic lists
   struct modelist_t {
     int nmodes;                               // number of modes
     std::string list_id;
     std::vector<int> idx;                    // index of the mode in the mode definition object
     std::vector<std::string> physprop_files; // physprop filename
     std::vector<int> idx_props;             // index of the mode properties in the physprop object
   };

   // gasses used in climate/diagnostic calculations
   std::vector<gaslist_t> gaslist;

   // list of aerosols used in climate/diagnostic calcs
   std::vector<aerlist_t> aerosollist;

   // list of aerosol modes used in climate/diagnostic calcs
   std::vector<modelist_t> ma_list;

   // mode definitions
   modes_t modes;

   // physics prop
   PhysProp prop;

   const int num_mode_types = 8;
   const int num_spec_types = 8;
   const std::vector<std::string> mode_type_names {"accum", "aitken", "primary_carbon", "fine_seasalt",
                                                   "fine_dust", "coarse", "coarse_seasalt", "coarse_dust"};

   const std::vector<std::string> spec_type_names {"sulfate", "ammonium", "nitrate", "p-organic",
                                                   "s-organic", "black-c", "seasalt", "dust"};
  public:
   // constructor
   MamConstituents() {
     //
     // NOTE:
     // we hard-coded the use nmodes = 3, and nspeces = {6, 3, 3}
     //
     // allocate components that depend on nmodes
     int nmodes = 3;
     std::vector<int> num_spec{6, 3, 3};

     modes.nmodes = 3;
     modes.names.resize(nmodes);
     modes.types.resize(nmodes);
     modes.comps.resize(nmodes);

     modes.comps[0].nspec = num_spec[0];
     modes.names[0] = "mam3_mode_0";
     modes.types[0] = "accum";

     modes.comps[0].source_num_a = "N";
     modes.comps[0].camname_num_a = "num_a1";
     modes.comps[0].source_num_c = "N";
     modes.comps[0].camname_num_c = "num_c1";

     modes.comps[0].type.resize(num_spec[0]);
     modes.comps[0].type[0] = "sulfate";
     modes.comps[0].type[1] = "p-organic";
     modes.comps[0].type[2] = "s-organic";
     modes.comps[0].type[3] = "black-c";
     modes.comps[0].type[4] = "dust";
     modes.comps[0].type[5] = "seasalt";

     modes.comps[0].props.resize(num_spec[0]);
     modes.comps[0].props[0] = "data/sulfate_rrtmg_c080918.nc";
     modes.comps[0].props[1] = "data/ocpho_rrtmg_c101112.nc";
     modes.comps[0].props[2] = "data/ocphi_rrtmg_c100508.nc";
     modes.comps[0].props[3] = "data/bcpho_rrtmg_c100508.nc";
     modes.comps[0].props[4] = "data/dust_aeronet_rrtmg_c141106.nc";
     modes.comps[0].props[5] = "data/ssam_rrtmg_c100508.nc";

     modes.comps[1].nspec = num_spec[1];
     modes.names[1] = "mam3_mode_1";
     modes.types[1] = "aitken";

     modes.comps[1].source_num_a = "N";
     modes.comps[1].camname_num_a = "num_a2";
     modes.comps[1].source_num_c = "N";
     modes.comps[1].camname_num_c = "num_c2";

     modes.comps[1].type.resize(num_spec[1]);
     modes.comps[1].type[0] = "sulfate";
     modes.comps[1].type[1] = "s-organic";
     modes.comps[1].type[2] = "seasalt";

     modes.comps[1].props.resize(num_spec[1]);
     modes.comps[1].props[0] = "data/sulfate_rrtmg_c080918.nc";
     modes.comps[1].props[1] = "data/ocphi_rrtmg_c100508.nc";
     modes.comps[1].props[2] = "data/ssam_rrtmg_c100508.nc";

     modes.comps[2].nspec = num_spec[2];
     modes.names[2] = "mam3_mode_2";
     modes.types[2] = "coarse";

     modes.comps[2].source_num_a = "N";
     modes.comps[2].camname_num_a = "num_a3";
     modes.comps[2].source_num_c = "N";
     modes.comps[2].camname_num_c = "num_c3";

     modes.comps[2].type.resize(num_spec[2]);
     modes.comps[2].type[0] = "dust";
     modes.comps[2].type[1] = "seasalt";
     modes.comps[2].type[2] = "sulfate";

     modes.comps[2].props.resize(num_spec[2]);
     modes.comps[2].props[0] = "data/dust_aeronet_rrtmg_c141106.nc";
     modes.comps[2].props[1] = "data/ssam_rrtmg_c100508.nc";
     modes.comps[2].props[2] = "data/sulfate_rrtmg_c080918.nc";

     gaslist.resize(N_DIAG);
     aerosollist.resize(N_DIAG);
     ma_list.resize(N_DIAG);

     // Set the list_id fields which distinquish the climate and diagnostic lists
     for(auto i = 0; i < N_DIAG; ++i) {
        aerosollist[i].list_id = "";
        gaslist[i].list_id     = "";
        ma_list[i].list_id     = "";
     }

     gaslist[0].ngas = 8;
     gaslist[0].gas.resize(gaslist[0].ngas);
     gaslist[0].gas[0].source = "A";
     gaslist[0].gas[1].source = "N";
     gaslist[0].gas[2].source = "N";
     gaslist[0].gas[3].source = "N";
     gaslist[0].gas[4].source = "N";
     gaslist[0].gas[5].source = "N";
     gaslist[0].gas[6].source = "N";
     gaslist[0].gas[7].source = "N";

     gaslist[0].gas[0].camname = "Q";
     gaslist[0].gas[1].camname = "O3";
     gaslist[0].gas[2].camname = "O2";
     gaslist[0].gas[3].camname = "CO2";
     gaslist[0].gas[4].camname = "N2O";
     gaslist[0].gas[5].camname = "CH4";
     gaslist[0].gas[6].camname = "CFC11";
     gaslist[0].gas[7].camname = "CFC12";

     ma_list[0].nmodes = 3;
     ma_list[0].idx.resize(ma_list[0].nmodes);
     ma_list[0].idx[0] = 0;
     ma_list[0].idx[1] = 1;
     ma_list[0].idx[2] = 2;

     ma_list[0].physprop_files.resize(ma_list[0].nmodes);
     ma_list[0].idx_props.resize(ma_list[0].nmodes);

     ma_list[0].physprop_files.resize(ma_list[0].nmodes);
     ma_list[0].physprop_files[0] = "data/mam3_mode1_rrtmg_c110318.nc";
     ma_list[0].physprop_files[1] = "data/mam3_mode2_rrtmg_c110318.nc";
     ma_list[0].physprop_files[2] = "data/mam3_mode3_rrtmg_c110318.nc";

     ma_list[0].idx_props.resize(ma_list[0].nmodes);
     ma_list[0].idx_props[0] = 78384393;
     ma_list[0].idx_props[1] = 8192;
     ma_list[0].idx_props[2] = 78384392;

     // Start with the bulk aerosol species in the climate/diagnostic lists.
     // The physprop_accum_unique_files routine has the side effect of returning the number
     // of bulk aerosols in each list (they're identified by type='A').
     for(auto i = 0; i < ma_list[0].nmodes; ++i)
         prop.physprop_accum_unique_files(ma_list[0].physprop_files[i], "A");

     // Add physprop files for the species from the mode definitions.
     for(auto i = 0; i < modes.nmodes; ++i) {
         for(auto j = 0; j < modes.comps[i].nspec; ++j)
           prop.physprop_accum_unique_files(modes.comps[i].props[j], "M");
     }

     // initialize the gas and aerosol lists
     rad_cnst_init();
   }

   // The initialization of the gas and aerosol lists is finished by
   // 1) read the physprop files
   // 2) find the index of each constituent in the constituent or physics buffer arrays
   // 3) find the index of the aerosol constituents used to access its properties from the
   //    physprop module.
   inline
   void rad_cnst_init() {
     int num_aerosols;
     constexpr bool stricttest = true;

     // memory to point to if zero value requested
     //allocate(zero_cols(pcols,pver))
     //zero_cols = 0.;

     // Allocate storage for the physical properties of each aerosol; read properties from
     // the data files.
     prop.physprop_init();

     // Finish initializing the mode definitions.
     init_mode_comps();

     for(auto i = 0; i < ma_list.size(); ++i) {
       for(auto imode = 0; imode < ma_list[i].nmodes; ++imode) {
         // get the physprop_id from the phys_prop module
         ma_list[i].idx_props[imode] = prop.physprop_get_id(ma_list[i].physprop_files[imode]);
       }
     }
   }

 inline
 void init_mode_comps() {
   for(auto m = 0; m < modes.nmodes; ++m) {
      for(auto ispec = 0; ispec < modes.comps[m].nspec; ++ispec) {
         modes.comps[m].idx_props[ispec] = prop.physprop_get_id(modes.comps[m].props[ispec]);
         if (modes.comps[m].idx_props[ispec] == -1)
            amrex::Print() << "ERROR idx not found for " << modes.comps[m].props[ispec] << std::endl;
      }
    }
  }

  // Return pointer to mass mixing ratio for the gas from the specified
  // climate or diagnostic list.
  const inline
  void rad_cnst_get_gas(int list_idx, std::string gasname, real2d mmr) {
     gaslist_t list;
     if (list_idx >= 0 && list_idx <= N_DIAG) {
        list = gaslist[list_idx];
     } else {
        amrex::Print() << " list_idx= " << list_idx << "\n";
     }

     // Get index of gas in internal arrays.  rad_gas_index will abort if the
     // specified gasname is not recognized by the radiative transfer code.
     auto igas = RadConstants::rad_gas_index(gasname);

     // Get data source
     auto source = list.gas[igas].source;
     auto idx    = list.gas[igas].idx;

//   select case( source )
//   case ('A')
//     mmr => state%q(:,:,idx)
//   case ('N')
//     pbuf_get_field(pbuf, idx, mmr)
//   case ('Z')
//     mmr => zero_cols
//   end select
 }

 const
 void get_nmodes(int list_idx, int nmodes) {
    auto m_list = ma_list[list_idx];
    nmodes = m_list.nmodes;
 }

 const
 void get_ngas(int list_idx, int ngas) {
    auto g_list = gaslist[list_idx];
    ngas = g_list.ngas;
 }

 const
 void get_naero(int list_idx, int naero) {
    auto a_list = aerosollist[list_idx];
    naero = a_list.numaerosols;
 }

 const
 void get_gas_names(int list_idx, std::vector<std::string> gasnames, bool use_data_o3=false) {
   auto g_list = gaslist[list_idx];

   gasnames.reserve(g_list.ngas);

   for(auto i = 1; i < g_list.ngas; ++i) gasnames[i] = g_list.gas[i].camname;

   // Does the climate calculation use data ozone?
   // get index of O3 in gas list
   auto igas = RadConstants::rad_gas_index("O3");
   // Get data source
   auto source = g_list.gas[igas].source;
   use_data_o3 = false;
   if (source == "N") use_data_o3 = true;
 }

 const
 void get_aero_names(int list_idx, std::vector<std::string> aernames) {
   auto a_list = aerosollist[list_idx];

   aernames.reserve(a_list.numaerosols);

   for(auto i=0; i<a_list.numaerosols; ++i) aernames[i] = a_list.aer[i].camname;
 }

 const
 void get_mode_nspec(int list_idx, int m_idx, int nspec) {
   auto m_list = ma_list[list_idx];
   auto mm = m_list.idx[m_idx];

   // number of species in the mode
   nspec = modes.comps[mm].nspec;
 }

 // Return info about modal aerosol list
 const
 void rad_cnst_get_info_by_mode(const int& list_idx, const int& m_idx,
                                std::string mode_type = "", std::string num_name = "",
                                std::string num_name_cw = "", int nspec = 0) {
   auto m_list = ma_list[list_idx];

   // check for valid mode index
   auto nmodes = m_list.nmodes;
   if (m_idx < 1 || m_idx > nmodes)
      amrex::Print() <<  "ERROR - invalid mode index: " << m_idx << "\n";

   // get index into the mode definition object
   auto mm = m_list.idx[m_idx];

   // mode type
   if (mode_type.empty()) mode_type = modes.types[mm];

   // number of species in the mode
   if (nspec == 0) nspec = modes.comps[mm].nspec;

   // name of interstitial number mixing ratio
   if (num_name.empty()) num_name = modes.comps[mm].camname_num_a;

   // name of cloud borne number mixing ratio
   if (num_name_cw.empty()) num_name_cw = modes.comps[mm].camname_num_c;
 }

 // Return info about modal aerosol lists
 const
 void rad_cnst_get_info_by_mode_spec(const int& list_idx, const int& m_idx, const int& s_idx,
                                     std::string spec_type = "", std::string spec_name = "",
                                     std::string spec_name_cw = "") {
   auto m_list = ma_list[list_idx];

   // check for valid mode index
   auto nmodes = m_list.nmodes;
   if (m_idx < 1 || m_idx > nmodes)
      amrex::Print() <<  "ERROR - invalid mode index: " << m_idx << "\n";

   // get index into the mode definition object
   auto mm = m_list.idx[m_idx];

   // check for valid specie index
   auto nspec = modes.comps[mm].nspec;
   if (s_idx < 1 || s_idx > nspec)
       amrex::Print() << "ERROR - invalid specie index: " << s_idx << "\n";

   // specie type
   if (spec_type.empty()) spec_type = modes.comps[mm].type[s_idx];

   // interstitial specie name
   if (spec_name.empty()) spec_name = modes.comps[mm].camname_mmr_a[s_idx];

   // cloud borne specie name
   if (spec_name_cw.empty()) spec_name_cw = modes.comps[mm].camname_mmr_c[s_idx];
 }

 // Return info about modes in the specified climate/diagnostics list
 const
 void rad_cnst_get_info_by_spectype(const int& list_idx, const std::string& spectype,
                                    int mode_idx = 0, int spec_idx = 0) {
   auto m_list = ma_list[list_idx];

   // number of modes in specified list
   auto nmodes = m_list.nmodes;

   // loop through modes in specified climate/diagnostic list
   auto found_spectype = false;
   for(auto i = 0; i < nmodes; ++i) {

      // get index of the mode in the definition object
      auto m_idx = m_list.idx[i];

      // number of species in the mode
      auto nspec = modes.comps[m_idx].nspec;

      // loop through species looking for spectype
      for(auto ispec = 0; ispec < nspec; ++ispec) {
         if (modes.comps[m_idx].type[ispec] == spectype) {
            if (mode_idx) mode_idx = i;
            if (spec_idx) spec_idx = ispec;
            found_spectype = true;
            exit(0);
         }
      }
      if (found_spectype) exit(0);
   }

   if (!found_spectype) {
      mode_idx = -1;
      spec_idx = -1;
   }
 }

 const
 int rad_cnst_get_mode_idx(const int& list_idx, std::string& mode_type) {
   // if mode type not found return -1
   int mode_idx = -1;

   // specified mode list
   auto m_list = ma_list[list_idx];

   // number of modes in specified list
   auto nmodes = m_list.nmodes;

   // loop through modes in specified climate/diagnostic list
   for(auto i = 0; i < nmodes; ++i) {
      // get index of the mode in the definition object
      auto m_idx = m_list.idx[i];

      // look in mode definition object (modes) for the mode types
      if (modes.types[m_idx] == mode_type) {
         mode_idx = i;
         exit(0);
      }
   }
   return mode_idx;
 }

 const
 int rad_cnst_get_spec_idx(const int& list_idx, const int& mode_idx, const std::string& spec_type)
  {
    // if specie type not found return -1
    int spec_idx = -1;

    // modes in specified list
    auto m_list = ma_list[list_idx];

    // get index of the specified mode in the definition object
    auto m_idx = m_list.idx[mode_idx];

    // object containing the components of the mode
    auto mode_comps = modes.comps[m_idx];

    // number of species in specified mode
    auto nspec = mode_comps.nspec;

    // loop through species in specified mode
    for(auto i=0; i<nspec; ++i) {
      // look in mode definition object (modes) for the mode types
      if (mode_comps.type[i] == spec_type) {
         spec_idx = i;
         exit(0);
      }
    }
    return spec_idx;
  }

 // Output the mass per layer, and total column burdens for gas and aerosol
 // constituents in either the climate or diagnostic lists
 const
 void rad_cnst_out(int list_idx) {
   int ncol;
   int idx;
   std::string name, cbname, source;

   aerlist_t aerlist;
   // Associate pointer with requested aerosol list
   if (list_idx >= 0 && list_idx <= N_DIAG) {
      aerlist = aerosollist[list_idx];
   } else {
      amrex::Print() << "list_idx = " << list_idx << "\n";
   }

   auto naer = aerlist.numaerosols;
   for(auto i = 0; i < naer; ++i) {
      source = aerlist.aer[i].source;
      idx    = aerlist.aer[i].idx;
      name   = aerlist.aer[i].mass_name;
      // construct name for column burden field by replacing the 'm_' prefix by 'cb_'
      cbname = "cb_"; // name(3:len_trim(name))

//      select case( source )
//      case ('A')
//         mmr => state%q(:,:,idx)
//      case ('N')
//         call pbuf_get_field(pbuf, idx, mmr)
//      end select

//      mass(:ncol,:) = mmr(:ncol,:) * state%pdeldry(:ncol,:) * rga

//      cb(:ncol) = sum(mass(:ncol,:),2);

   }

   // Associate pointer with requested gas list
   auto g_list = gaslist[list_idx];

   auto ngas = g_list.ngas;
   for(auto i = 0; i < ngas; ++i) {
      source = g_list.gas[i].source;
      idx    = g_list.gas[i].idx;
      name   = g_list.gas[i].mass_name;
      cbname = "cb_"; // name(3:len_trim(name))

//      select case( source )
//      case ('A')
//         mmr => state%q(:,:,idx)
//      case ('N')
//         call pbuf_get_field(pbuf, idx, mmr)
//      end select

//      mass(:ncol,:) = mmr(:ncol,:) * state%pdeldry(:ncol,:) * rga

//      cb(:ncol) = sum(mass(:ncol,:),2)
   }
 }

 // Return pointer to mass mixing ratio for the aerosol from the specified
 // climate or diagnostic list.
 const
 void rad_cnst_get_aer_mmr_by_idx(const int& list_idx, const int& aer_idx, real2d& mmr)
 {
   aerlist_t aerlist;

   if (list_idx >= 0 && list_idx <= N_DIAG) {
      aerlist = aerosollist[list_idx];
   } else {
      amrex::Print() << " list_idx = " << list_idx << std::endl;
      exit(0);
   }

   // Check for valid input aerosol index
   if (aer_idx < 1 || aer_idx > aerlist.numaerosols) {
      amrex::Print() << " aer_idx= " <<  aer_idx << "  numaerosols= " << aerlist.numaerosols << std::endl;
      exit(0);
   }

   // Get data source
   auto source = aerlist.aer[aer_idx].source;
   auto idx    = aerlist.aer[aer_idx].idx;
 }

 // Return pointer to mass mixing ratio for the modal aerosol specie from the specified
 // climate or diagnostic list.
 const
 void rad_cnst_get_mam_mmr_by_idx(const int& list_idx, const int& mode_idx, const int& spec_idx,
                                  const std::string& phase, real2d& mmr)
 {
   int idx;
   std::string source;
   modelist_t mlist;

   if (list_idx >= 0 && list_idx <= N_DIAG) {
      mlist = ma_list[list_idx];
   } else {
      amrex::Print() << "list_idx =" << list_idx << std::endl;
      exit(0);
   }

   // Check for valid mode index
   if (mode_idx < 1  ||  mode_idx > mlist.nmodes) {
      amrex::Print() << " mode_idx= " <<  mode_idx << " nmodes= " << mlist.nmodes << std::endl;
      exit(0);
   }

   // Get the index for the corresponding mode in the mode definition object
   auto m_idx = mlist.idx[mode_idx];

   // Check for valid specie index
   if (spec_idx < 1  ||  spec_idx > modes.comps[m_idx].nspec) {
      amrex::Print() << " spec_idx= " << spec_idx << "  nspec= " << modes.comps[m_idx].nspec << std::endl;
      exit(0);
   }

   // Get data source
   if (phase == "a") {
      source = modes.comps[m_idx].source_mmr_a[spec_idx];
//      idx    = modes.comps[m_idx].idx_mmr_a[spec_idx];
   } else if (phase == "c") {
      source = modes.comps[m_idx].source_mmr_c[spec_idx];
//      idx    = modes.comps[m_idx].idx_mmr_c[spec_idx];
   } else {
      amrex::Print() << "phase= " << phase << std::endl;
      exit(0);
   }
 }

 const
 void rad_cnst_get_mam_mmr_idx(const int& mode_idx, const int& spec_idx, int idx)
 {
   modelist_t mlist;

   // assume climate list (i.e., species are in the constituent array)
   mlist = ma_list[0];

   // Check for valid mode index
   if (mode_idx < 1  ||  mode_idx > mlist.nmodes) {
      amrex::Print() << "mode_idx= " << mode_idx << " nmodes= " << mlist.nmodes << std::endl;
      exit(0);
   }

   // Get the index for the corresponding mode in the mode definition object
   auto m_idx = mlist.idx[mode_idx];

   // Check for valid specie index
   if (spec_idx < 1  ||  spec_idx > modes.comps[m_idx].nspec) {
      amrex::Print() << "spec_idx= " << spec_idx << " nspec= " << modes.comps[m_idx].nspec << std::endl;
      exit(0);
   }

   // Assume data source is interstitial since that's what's in the constituent array
 //  idx = modes.comps[m_idx].idx_mmr_a[spec_idx];
 }

 // Return pointer to number mixing ratio for the aerosol mode from the specified
 // climate or diagnostic list.
 const
 void rad_cnst_get_mode_num(int list_idx, int mode_idx, std::string phase, real2d& num)
 {
   modelist_t mlist;
   std::string source;
   int idx;

   if (list_idx >= 0 && list_idx <= N_DIAG) {
      mlist = ma_list[list_idx];
   } else {
      amrex::Print() << " list_idx =" << list_idx << std::endl;
      exit(0);
   }

   // Check for valid mode index
   if (mode_idx < 1  ||  mode_idx > mlist.nmodes) {
      amrex::Print() << " mode_idx= " << mode_idx << "  nmodes= " << mlist.nmodes << std::endl;
      exit(0);
   }

   // Get the index for the corresponding mode in the mode definition object
   auto m_idx = mlist.idx[mode_idx];

   // Get data source
   if (phase == "a") {
      source = modes.comps[m_idx].source_num_a;
//      idx    = modes.comps[m_idx].idx_num_a;
   } else if (phase == "c") {
      source = modes.comps[m_idx].source_num_c;
//      idx    = modes.comps[m_idx].idx_num_c;
   } else {
      amrex::Print() << " phase= " << phase << std::endl;
      exit(1);
   }
 }

 // Return constituent index of mode number mixing ratio for the aerosol mode in
 // the climate list.

 // This is a special routine to allow direct access to information in the
 // constituent array inside physics parameterizations that have been passed,
 // and are operating over the entire constituent array.  The interstitial phase
 // is assumed since that's what is contained in the constituent array.
 const
 void rad_cnst_get_mode_num_idx(const int& mode_idx, int& cnst_idx)
 {
   modelist_t mlist;
   // assume climate list
   mlist = ma_list[0];

   // Check for valid mode index
   if (mode_idx < 1  ||  mode_idx > mlist.nmodes) {
      amrex::Print() << " mode_idx= " << mode_idx << " nmodes= " << mlist.nmodes << std::endl;
      exit(0);
   }

   // Get the index for the corresponding mode in the mode definition object
   auto m_idx = mlist.idx[mode_idx];

   // Check that source is 'A' which means the index is for the constituent array
   auto source = modes.comps[m_idx].source_num_a;
   if (source != "A") {
      amrex::Print() << " source= " << source << std::endl;
      exit(0);
   }

   // Return index in constituent array
//   cnst_idx = modes.comps[m_idx].idx_num_a;
 }

 // Return the index of aerosol aer_name in the list specified by list_idx.
 const
 int rad_cnst_get_aer_idx(int list_idx, std::string aer_name)
 {
   aerlist_t aerlist;

   if (list_idx >= 0 && list_idx <= N_DIAG) {
      aerlist = aerosollist[list_idx];
   }
   else {
      amrex::Print() << " list_idx =" << list_idx << std::endl;
      exit(0);
   }

   // Get index in aerosol list for requested name
   auto aer_idx = -1;
   for(auto i = 1; i < aerlist.numaerosols; ++i) {
      if (aer_name == aerlist.aer[i].camname) {
         aer_idx = i;
         exit(0);
      }
   }

   if (aer_idx == -1) exit(0);

   return aer_idx;
 }

 // Return requested properties for the mode from the specified
 // climate or diagnostic list.
 const
 void get_mode_props(int list_idx, int mode_idx, real& sigmag,
                     real& rhcrystal, real& rhdeliques)
 {
    modelist_t mlist;
   if (list_idx >= 0 && list_idx <= N_DIAG)
      mlist = ma_list[list_idx];

   // Get the physprop index for the requested mode
   auto id = mlist.idx_props[mode_idx];

   prop.get_sigmag(id, sigmag);
   prop.get_rhcrystal(id, rhcrystal);
   prop.get_rhdeliques(id, rhdeliques);
 }

 const
 void get_mode_props(int list_idx, int mode_idx, real& sigmag, real2d& refrtablw,
                     real2d& refitablw, real4d& absplw)
 {
    modelist_t mlist;
   if (list_idx >= 0 && list_idx <= N_DIAG)
      mlist = ma_list[list_idx];

   // Get the physprop index for the requested mode
   auto id = mlist.idx_props[mode_idx];

   prop.get_sigmag(id, sigmag);
   prop.get_refrtablw(id, refrtablw);
   prop.get_refitablw(id, refitablw);
   prop.get_absplw(id, absplw);
 }

 const
 void get_mode_props(int list_idx, int mode_idx, real& sigmag, real2d& refrtabsw,
                     real2d& refitabsw, real4d& extpsw, real4d& abspsw, real4d& asmpsw)
 {
   modelist_t mlist;
   if (list_idx >= 0 && list_idx <= N_DIAG)
      mlist = ma_list[list_idx];

   // Get the physprop index for the requested mode
   auto id = mlist.idx_props[mode_idx];

   prop.get_sigmag(id, sigmag);
   prop.get_refrtabsw(id, refrtabsw);
   prop.get_refitabsw(id, refitabsw);
   prop.get_extpsw(id, extpsw);
   prop.get_abspsw(id, abspsw);
   prop.get_asmpsw(id, asmpsw);
 }

 const
 void get_mode_props(int list_idx, int mode_idx, int ncoef, int prefr, int prefi)
 {
    modelist_t mlist;
   if (list_idx >= 0 && list_idx <= N_DIAG)
      mlist = ma_list[list_idx];

   // Get the physprop index for the requested mode
   auto id = mlist.idx_props[mode_idx];

   prop.get_ncoef(id, ncoef);
   prop.get_prefr(id, prefr);
   prop.get_prefi(id, prefi);
 }

 const inline
 void get_mode_props(int list_idx, int mode_idx, real& dgnum, real& dgnumhi,
                     real& dgnumlo, real& sigmag)
 {
    modelist_t mlist;
   if (list_idx >= 0 && list_idx <= N_DIAG)
      mlist = ma_list[list_idx];

   // Get the physprop index for the requested mode
   auto id = mlist.idx_props[mode_idx];

   prop.get_dgnum(id, dgnum);
   prop.get_dgnumhi(id, dgnumhi);
   prop.get_dgnumlo(id, dgnumlo);
   prop.get_sigmag(id, sigmag);
 }

 const inline
 void get_mam_density_aer(int list_idx, int mode_idx, int spec_idx,
                          real& density_aer)
 {
   modelist_t mlist;
   if (list_idx >= 0 && list_idx <= N_DIAG)
      mlist = ma_list[list_idx];

   // Get the index for the corresponding mode in the mode definition object
   auto m_idx = mlist.idx[mode_idx];
   auto id = modes.comps[m_idx].idx_props[spec_idx];
   prop.get_density_aer(id, density_aer);
 }

 const inline
 void get_mam_hygro_aer(int list_idx, int mode_idx, int spec_idx,
                       real& hygro_aer)
 {
   modelist_t mlist;
   if (list_idx >= 0 && list_idx <= N_DIAG)
      mlist = ma_list[list_idx];

   // Get the index for the corresponding mode in the mode definition object
   auto m_idx = mlist.idx[mode_idx];
   auto id = modes.comps[m_idx].idx_props[spec_idx];
   prop.get_hygro_aer(id, hygro_aer);
 }

 // Return requested properties for the aerosol from the specified
 // climate or diagnostic list.
 const inline
 void get_mam_props(int list_idx, int mode_idx, int spec_idx, real& density_aer,
                    std::string& spectype, real& hygro_aer,
                    real1d& refindex_real_aer_sw, real1d& refindex_im_aer_sw)
 {
   modelist_t mlist;
   if (list_idx >= 0 && list_idx <= N_DIAG)
      mlist = ma_list[list_idx];

   // Get the index for the corresponding mode in the mode definition object
   auto m_idx = mlist.idx[mode_idx];
   auto id = modes.comps[m_idx].idx_props[spec_idx];

   prop.get_density_aer(id, density_aer);
   prop.get_hygro_aer(id, hygro_aer);
   prop.get_ref_real_aer_sw(id, refindex_real_aer_sw);
   prop.get_ref_im_aer_sw(id, refindex_im_aer_sw);
   spectype = modes.comps[m_idx].type[spec_idx];
 }

 const inline
 void get_mam_props_sw(int list_idx, int mode_idx, int spec_idx,
                    real& density_aer, real1d& refindex_real_aer_sw, real1d& refindex_im_aer_sw)
 {
   modelist_t mlist;
   if (list_idx >= 0 && list_idx <= N_DIAG)
      mlist = ma_list[list_idx];

   // Get the index for the corresponding mode in the mode definition object
   auto m_idx = mlist.idx[mode_idx];
   auto id = modes.comps[m_idx].idx_props[spec_idx];

   prop.get_density_aer(id, density_aer);
   prop.get_ref_real_aer_sw(id, refindex_real_aer_sw);
   prop.get_ref_im_aer_sw(id, refindex_im_aer_sw);
 }

 const inline
 void get_mam_props(int list_idx, int mode_idx, int spec_idx, real& density_aer)
 {
   modelist_t mlist;
   if (list_idx >= 0 && list_idx <= N_DIAG)
      mlist = ma_list[list_idx];

   auto m_idx = mlist.idx[mode_idx];
   auto id = modes.comps[m_idx].idx_props[spec_idx];

   prop.get_density_aer(id, density_aer);
 }

 const inline
 void get_mam_props_lw(int list_idx, int mode_idx, int spec_idx,
                       real& density_aer, real1d& refindex_real_aer_lw, real1d& refindex_im_aer_lw)
 {
   modelist_t mlist;
   if (list_idx >= 0 && list_idx <= N_DIAG)
      mlist = ma_list[list_idx];

   // Get the index for the corresponding mode in the mode definition object
   auto m_idx = mlist.idx[mode_idx];
   auto id = modes.comps[m_idx].idx_props[spec_idx];

   prop.get_density_aer(id, density_aer);
   prop.get_ref_real_aer_lw(id, refindex_real_aer_lw);
   prop.get_ref_im_aer_lw(id, refindex_im_aer_lw);
 }

 const inline
 void get_aer_opticstype(int list_idx, int aer_idx, std::string opticstype)
 {
   auto aerlist = aerosollist[list_idx];
   auto id = aerlist.aer[aer_idx].physprop_id;

   prop.get_opticstype(id, opticstype);
 }

 const inline
 void get_aer_sw_hygro_ext(int list_idx, int aer_idx, real2d& sw_hygro_ext)
 {
   auto aerlist = aerosollist[list_idx];
   auto id = aerlist.aer[aer_idx].physprop_id;

   prop.get_sw_hygro_ext(id, sw_hygro_ext);
 }

 const inline
 void get_aer_sw_hygro_ssa(int list_idx, int aer_idx, real2d& sw_hygro_ssa)
 {
   auto aerlist = aerosollist[list_idx];
   auto id = aerlist.aer[aer_idx].physprop_id;

   prop.get_sw_hygro_ssa(id, sw_hygro_ssa);
 }

 const inline
 void get_aer_sw_hygro_asm(int list_idx, int aer_idx, real2d& sw_hygro_asm)
 {
   auto aerlist = aerosollist[list_idx];
   auto id = aerlist.aer[aer_idx].physprop_id;

   prop.get_sw_hygro_asm(id, sw_hygro_asm);
 }

 const inline
 void get_aer_lw_hygro_abs(int list_idx, int aer_idx, real2d& lw_hygro_abs)
 {
   auto aerlist = aerosollist[list_idx];
   auto id = aerlist.aer[aer_idx].physprop_id;

   prop.get_lw_hygro_abs(id, lw_hygro_abs);
 }

 const inline
 void get_aer_sw_nonhygro_ext(int list_idx, int aer_idx, real1d& sw_nonhygro_ext)
 {
   auto aerlist = aerosollist[list_idx];
   auto id = aerlist.aer[aer_idx].physprop_id;

   prop.get_sw_nonhygro_ext(id, sw_nonhygro_ext);
 }

 const inline
 void get_aer_sw_nonhygro_ssa(int list_idx, int aer_idx, real1d& sw_nonhygro_ssa)
 {
   auto aerlist = aerosollist[list_idx];
   auto id = aerlist.aer[aer_idx].physprop_id;

   prop.get_sw_nonhygro_ssa(id, sw_nonhygro_ssa);
 }

 const inline
 void get_aer_sw_nonhygro_asm(int list_idx, int aer_idx, real1d& sw_nonhygro_asm)
 {
   auto aerlist = aerosollist[list_idx];
   auto id = aerlist.aer[aer_idx].physprop_id;

   prop.get_sw_nonhygro_asm(id, sw_nonhygro_asm);
 }

 const inline
 void get_aer_sw_nonhygro_scat(int list_idx, int aer_idx, real1d& sw_nonhygro_scat)
 {
   auto aerlist = aerosollist[list_idx];
   auto id = aerlist.aer[aer_idx].physprop_id;

   prop.get_sw_nonhygro_scat(id, sw_nonhygro_scat);
 }

 const inline
 void get_aer_sw_nonhygro_ascat(int list_idx, int aer_idx, real1d& sw_nonhygro_ascat)
 {
   auto aerlist = aerosollist[list_idx];
   auto id = aerlist.aer[aer_idx].physprop_id;

   prop.get_sw_nonhygro_ascat(id, sw_nonhygro_ascat);
 }

 const inline
 void get_aer_lw_abs(int list_idx, int aer_idx, real1d& lw_abs)
 {
   auto aerlist = aerosollist[list_idx];
   auto id = aerlist.aer[aer_idx].physprop_id;

   prop.get_lw_abs(id, lw_abs);
 }

 const inline
 void get_aer_refindex_aer_sw(int list_idx, int aer_idx, real1d& refindex_real_aer_sw, real1d& refindex_im_aer_sw)
 {
   auto aerlist = aerosollist[list_idx];
   auto id = aerlist.aer[aer_idx].physprop_id;

   prop.get_ref_real_aer_sw(id, refindex_real_aer_sw);
   prop.get_ref_im_aer_sw(id, refindex_im_aer_sw);
 }

 const inline
 void get_aer_refindex_aer_lw(int list_idx, int aer_idx, real1d& refindex_real_aer_lw, real1d& refindex_im_aer_lw)
 {
   auto aerlist = aerosollist[list_idx];
   auto id = aerlist.aer[aer_idx].physprop_id;

   prop.get_ref_real_aer_lw(id, refindex_real_aer_lw);
   prop.get_ref_im_aer_lw(id, refindex_im_aer_lw);
 }

 const inline
 void get_aer_r_sw_ext(int list_idx, int aer_idx, real2d& r_sw_ext)
 {
   auto aerlist = aerosollist[list_idx];
   auto id = aerlist.aer[aer_idx].physprop_id;

   prop.get_r_sw_ext(id, r_sw_ext);
 }

 const inline
 void get_aer_r_sw_scat(int list_idx, int aer_idx, real2d& r_sw_scat)
 {
   auto aerlist = aerosollist[list_idx];
   auto id = aerlist.aer[aer_idx].physprop_id;

   prop.get_r_sw_scat(id, r_sw_scat);
 }

 const inline
 void get_aer_r_sw_ascat(int list_idx, int aer_idx, real2d& r_sw_ascat)
 {
   auto aerlist = aerosollist[list_idx];
   auto id = aerlist.aer[aer_idx].physprop_id;

   prop.get_r_sw_ascat(id, r_sw_ascat);
 }

 const inline
 void get_aer_r_lw_abs(int list_idx, int aer_idx, real2d& r_lw_abs)
 {
   auto aerlist = aerosollist[list_idx];
   auto id = aerlist.aer[aer_idx].physprop_id;

   prop.get_r_lw_abs(id, r_lw_abs);
 }

 const inline
 void get_aer_mu(int list_idx, int aer_idx, real1d& mu)
 {
   auto aerlist = aerosollist[list_idx];
   auto id = aerlist.aer[aer_idx].physprop_id;

   prop.get_mu(id, mu);
 }

 const inline
 void get_aername(int list_idx, int aer_idx, std::string aername)
 {
   auto aerlist = aerosollist[list_idx];
   auto id = aerlist.aer[aer_idx].physprop_id;

   prop.get_aername(id, aername);
 }

 const inline
 void get_density_aer(int list_idx, int aer_idx, real density_aer)
 {
   auto aerlist = aerosollist[list_idx];
   auto id = aerlist.aer[aer_idx].physprop_id;

   prop.get_density_aer(id, density_aer);
 }

 const inline
 void get_hygro_aer(int list_idx, int aer_idx, real hygro_aer)
 {
   auto aerlist = aerosollist[list_idx];
   auto id = aerlist.aer[aer_idx].physprop_id;

   prop.get_hygro_aer(id, hygro_aer);
 }

 const inline
 void get_dryrad_aer(int list_idx, int aer_idx, real dryrad_aer)
 {
   auto aerlist = aerosollist[list_idx];
   auto id = aerlist.aer[aer_idx].physprop_id;

   prop.get_dryrad_aer(id, dryrad_aer);
 }

 const inline
 void get_dispersion_aer(int list_idx, int aer_idx, real dispersion_aer)
 {
   auto aerlist = aerosollist[list_idx];
   auto id = aerlist.aer[aer_idx].physprop_id;

   prop.get_dispersion_aer(id, dispersion_aer);
 }

 const inline
 void get_num_to_mass_aer(int list_idx, int aer_idx, real num_to_mass_aer)
 {
   auto aerlist = aerosollist[list_idx];
   auto id = aerlist.aer[aer_idx].physprop_id;

   prop.get_num_to_mass_aer(id, num_to_mass_aer);
 }
};
#endif
