/** \file nc_interface.H
 *
 *  Interface to NetCDF library
 *
 *  Defines convenience wrappers to interact with a NetCDF file in a more
 *  OOP-like manner.
 */

#ifndef NC_INTERFACE_H
#define NC_INTERFACE_H

#ifndef USE_NETCDF
#include <string>
#include <unordered_map>
#include <vector>

#include <netcdf.h>
#include <netcdf_par.h>

namespace ncutils {

//! Wrapper around NetCDF data types
struct NCDType
{
    static constexpr nc_type Int = NC_INT;
#ifdef AMREX_USE_FLOAT
    static constexpr nc_type Real = NC_FLOAT;
    using RType = float;
#else
    static constexpr nc_type Real = NC_DOUBLE;
    using RType = double;
#endif
};

//! Representation of NetCDF dimension
struct NCDim
{
    //! File/Group Identifier
    const int ncid;

    //! Dimension ID used with NetCDF API
    const int dimid;

    //! Name of this dimension
    [[nodiscard]] std::string name() const;

    //! Length of this dimension
    [[nodiscard]] size_t len() const;
};

//! Representation of a NetCDF variable
struct NCVar
{
    //! File/Group identifier
    const int ncid;

    //! Variable ID used with NetCDF API
    const int varid;

    //! Name of this variable
    [[nodiscard]] std::string name () const;

    //! Number of array dimensions for this variable
    [[nodiscard]] int ndim () const;

    //! Shape of the array (size in each array dimension)
    [[nodiscard]] std::vector<size_t> shape () const;

    //! Write out the entire variable
    void put (const double*) const;
    void put (const float*) const;
    void put (const int*) const;

    //! Write out a slice of data
    void
    put (const double*,
         const std::vector<size_t>&,
         const std::vector<size_t>&) const;

    //! Write out a slice of data with with strides (see hyperslab definition in
    //! NetCDF)
    void
    put (const double*,
         const std::vector<size_t>&,
         const std::vector<size_t>&,
         const std::vector<ptrdiff_t>&) const;
    //! Write out a slice of data
    void
    put (const float*,
         const std::vector<size_t>&,
         const std::vector<size_t>&) const;

    //! Write out a slice of data with with strides (see hyperslab definition in
    //! NetCDF)
    void
    put (const float*,
         const std::vector<size_t>&,
         const std::vector<size_t>&,
         const std::vector<ptrdiff_t>&) const;

    void put (const int*, const std::vector<size_t>&, const std::vector<size_t>&) const;

    void
    put (const int*,
         const std::vector<size_t>&,
         const std::vector<size_t>&,
         const std::vector<ptrdiff_t>&) const;

    void put (const char**, const std::vector<size_t>&, const std::vector<size_t>&) const;

    void
    put (const char** dptr,
         const std::vector<size_t>& start,
         const std::vector<size_t>& count,
         const std::vector<ptrdiff_t>& stride) const;

    //! Read the entire variable from file
    void get (double*) const;
    void get (float*) const;
    void get (int*) const;

    //! Read a chunk of data from the file
    void
    get (double*, const std::vector<size_t>&, const std::vector<size_t>&) const;

    //! Read a chunk of data with strides
    void
    get (double*,
         const std::vector<size_t>&,
         const std::vector<size_t>&,
         const std::vector<ptrdiff_t>&) const;

    //! Read a chunk of data from the file
    void
    get (float*, const std::vector<size_t>&, const std::vector<size_t>&) const;

    //! Read a chunk of data with strides
    void
    get (float*,
         const std::vector<size_t>&,
         const std::vector<size_t>&,
         const std::vector<ptrdiff_t>&) const;

    void
    get (int*, const std::vector<size_t>&, const std::vector<size_t>&) const;

    void
    get (int*,
         const std::vector<size_t>&,
         const std::vector<size_t>&,
         const std::vector<ptrdiff_t>&) const;

    //! Read a chunk of data from the file
    void
    get (char*, const std::vector<size_t>&, const std::vector<size_t>&) const;

    //! Read a chunk of data with strides
    void
    get (char*,
         const std::vector<size_t>&,
         const std::vector<size_t>&,
         const std::vector<ptrdiff_t>&) const;

    [[nodiscard]] bool has_attr (const std::string& name) const;

    void put_attr (const std::string& name, const std::string& value) const;

    void
    put_attr (const std::string& name, const std::vector<double>& value) const;

    void
    put_attr (const std::string& name, const std::vector<float>& value) const;

    void put_attr (const std::string& name, const std::vector<int>& value) const;

    [[nodiscard]] std::string get_attr (const std::string& name) const;

    void get_attr (const std::string& name, std::vector<double>& value) const;
    void get_attr (const std::string& name, std::vector<float>& value) const;
    void get_attr (const std::string& name, std::vector<int>& value) const;

    void par_access (int cmode) const; //Uncomment for parallel NetCDF
};

//! Representation of a NetCDF group
class NCGroup
{
public:
    //! Identifier used with NetCDF API calls
    const int ncid;

    //! Name of this group
    [[nodiscard]] std::string name () const;

    //! Full name for this group
    [[nodiscard]] std::string full_name () const;

    //! Number of sub-groups within this group
    [[nodiscard]] int num_groups () const;

    //! Number of dimensions in this group
    [[nodiscard]] int num_dimensions () const;

    //! Number of variables within this group
    [[nodiscard]] int num_variables () const;

    //! Number of attributes within this group
    [[nodiscard]] int num_attributes () const;

    //! Check if a group exists
    [[nodiscard]] bool has_group (const std::string&) const;

    //! Check if a dimension exists by name
    [[nodiscard]] bool has_dim (const std::string&) const;

    //! Check if a variable exists by name
    [[nodiscard]] bool has_var (const std::string&) const;

    //! Check if an attribute exists
    [[nodiscard]] bool has_attr (const std::string&) const;

    /** Get the group by name
     *
     *  Throws error if the group doesn't exist, use `has_group` to check
     */
    [[nodiscard]] NCGroup group (const std::string&) const;

    //! Get the dimension instance by name
    [[nodiscard]] NCDim dim (const std::string&) const;

    //! Get the variable instance by name
    [[nodiscard]] NCVar var (const std::string&) const;

    //! Define new dimension
    void def_dim (const std::string&, size_t len) const;

    //! Define an array
    void def_array (const std::string& name,
                    nc_type dtype,
                    const std::vector<std::string>&) const;

    //! Define a variable (wrapper for def_array)
    void def_var (const std::string& name,
                  const nc_type dtype,
                  const std::vector<std::string>& dnames) const
    {
        def_array(name, dtype, dnames);
    }

    void put_attr (const std::string& name, const std::string& value) const;
    void put_attr (const std::string& name, const std::vector<double>& value) const;
    void put_attr (const std::string& name, const std::vector<float>& value) const;
    void put_attr (const std::string& name, const std::vector<int>& value) const;

    [[nodiscard]] std::string get_attr (const std::string& name) const;
    void get_attr (const std::string& name, std::vector<double>& value) const;
    void get_attr (const std::string& name, std::vector<float>& value) const;
    void get_attr (const std::string& name, std::vector<int>& value) const;

    //! Return a list of all groups defined in this group
    [[nodiscard]] std::vector<NCGroup> all_groups () const;

    //! Return a list of all dimensions defined in this group
    [[nodiscard]] std::vector<NCDim> all_dims () const;

    //! Return a list of all variables defined in this group
    [[nodiscard]] std::vector<NCVar> all_vars () const;

    //! Enter definition mode (not needed for NetCDF4 format)
    void enter_def_mode () const;

    //! Exit definition mode
    void exit_def_mode () const;

protected:
    NCGroup (const int id) : ncid(id) {}
    NCGroup (const int id, const NCGroup* /*par*/) : ncid(id) {}
};

/** Representation of a NetCDF file
 *
 *  Provide wrappes to create and open file
 */
class NCFile : public NCGroup
{
public:
    static NCFile
    create (const std::string& name, const int cmode = NC_CLOBBER | NC_NETCDF4);

    static NCFile open (const std::string& name, const int cmode = NC_NOWRITE);

    static NCFile create_par (const std::string& name,
                              const int cmode = NC_CLOBBER | NC_NETCDF4 | NC_MPIIO,
                              MPI_Comm comm = MPI_COMM_WORLD,
                              MPI_Info info = MPI_INFO_NULL);

    static NCFile open_par (const std::string& name,
                            const int cmode = NC_NOWRITE,
                            MPI_Comm comm = MPI_COMM_WORLD,
                            MPI_Info info = MPI_INFO_NULL);

    ~NCFile ();

    void close ();

protected:
    NCFile (const int id) : NCGroup(id), is_open{true} {}

    bool is_open{false};
};

} // namespace ncutils

#else

namespace ncutils {

struct NCDim
{
    const int ncid{-1};
    const int dimid{-1};
};

struct NCVar
{
    const int ncid{-1};
    const int varid{-1};
};

class NCGroup
{
public:
    const int ncid{-1};
};

class NCFile : public NCGroup
{};

} // namespace ncutils

#endif

#endif /* NC_INTERFACE_H */
