#ifndef ERF_WRF_SURFACE_FIELDS_H_
#define ERF_WRF_SURFACE_FIELDS_H_

/**
 * \file ERF_WRFSurfaceFields.H
 *
 * Land and soil fields of a wrfinput file that no land model takes: leaf area, vegetation
 * fraction, deep soil temperature, roughness, albedo, background albedo, sea ice, snow,
 * and the soil temperature and moisture profiles with the depths of their layers.
 *
 * ERF::init_from_wrfinput fills this holder when the two-stream force-restore balance
 * (erf.radiation.seb_prognostic_enable) runs without a land model. Every field is
 * optional: a file without it leaves get() null. The holder keeps the values as the file
 * gives them, in the file's units; it applies no physics.
 *
 * A field that holds a value outside its catalog range, or one that is not finite, or
 * that the level's files do not cover, is reported and not held: get() is null, as for a
 * field the file does not have. The code that takes a field decides whether its absence
 * stops the run.
 *
 * The data are those of the initial time. They are written to the checkpoint and read
 * back on a restart. A regrid carries them onto the new grids: a column the old grids
 * of the level held keeps its values, and any other column takes those of the column of
 * the coarser level it lies in (piecewise constant, so a flag stays a flag, a soil
 * profile stays the profile of one column, and no value leaves the catalog range). A
 * level without a file of its own (made by tagging, at the start or by a regrid) takes
 * every field of the coarser level in the same way.
 *
 * The class has no NetCDF dependency: fill() takes the region a rank has read as an
 * FArrayBox, so the unit tests can drive it without a file.
 */

#include <array>
#include <iosfwd>
#include <limits>
#include <memory>
#include <string>

#include <AMReX_BoxArray.H>
#include <AMReX_DistributionMapping.H>
#include <AMReX_FArrayBox.H>
#include <AMReX_MultiFab.H>
#include <AMReX_Periodicity.H>
#include <AMReX_Vector.H>

class WRFSurfaceFields
{
public:
    /** The variables, in the order they are requested from the file. */
    enum class Field : int {
        LAI = 0,   ///< leaf area index [m2/m2]
        VEGFRA,    ///< vegetation fraction [percent]
        TMN,       ///< soil temperature at the lower boundary [K]
        ZNT,       ///< roughness length [m]; in wrfout files, not in a real.exe wrfinput
        ALBEDO,    ///< surface albedo [1]; in wrfout files, not in a real.exe wrfinput
        ALBBCK,    ///< background (snow-free) albedo [1]; what a real.exe wrfinput holds
        SEAICE,    ///< sea ice flag [1]
        SNOW,      ///< snow water equivalent [kg/m2]
        SNOWH,     ///< snow depth [m]
        TSLB,      ///< soil temperature, one component per soil layer, top layer first [K]
        SMOIS,     ///< soil moisture, one component per soil layer, top layer first [m3/m3]
        ZS,        ///< depth of the center of each soil layer, one column for the level [m]
        DZS,       ///< thickness of each soil layer, one column for the level [m]
        NumFields
    };

    static constexpr int num_fields = static_cast<int>(Field::NumFields);

    /** How a variable is laid out in the file. */
    enum class Kind : int {
        Surface,      ///< one value per column
        SoilField,    ///< one value per column and soil layer
        LayerColumn   ///< one value per soil layer, the same for every column (ZS, DZS)
    };

    /** A variable's name in the file, its layout, its units, and the range a file may hold. */
    struct Descriptor {
        Field       id;
        const char* nc_name;
        Kind        kind;
        const char* units;
        amrex::Real lo;    ///< smallest value a file may hold
        amrex::Real hi;    ///< largest value a file may hold
    };

    /**
     * The one list of the variables. The reader asks the file for every name in it, and
     * takes the layout of a name from its kind.
     *
     * The bounds leave out only values that cannot be the named quantity (a negative
     * temperature, a fraction above one). A field that a hand-built file declares and
     * leaves at zero passes; whether zero is usable is for the code that takes the field.
     * The upper bounds of LAI, TMN, ZNT, SNOW, SNOWH, TSLB, ZS and DZS are not physical
     * limits: they are set far above any value of the quantity so that a NetCDF fill
     * value (9.97e36), which is finite, is left out.
     */
    [[nodiscard]] static const std::array<Descriptor, num_fields>& catalog ();

    /** The descriptor of the variable with this name in the file, or null when it is not one of these. */
    [[nodiscard]] static const Descriptor* find (const std::string& nc_name);

    /** True when the third index of this variable in the file is the soil layer. */
    [[nodiscard]] static bool has_soil_layers (const std::string& nc_name)
    {
        const Descriptor* d = find(nc_name);
        return d != nullptr && d->kind != Kind::Surface;
    }

    void resize (int nlevs);

    [[nodiscard]] int num_levels () const { return static_cast<int>(m_fields.size()); }

    /** True when the level holds a field or a soil layer column. */
    [[nodiscard]] bool holds_any (int lev) const;

    /** Drop everything held at this level. Prints one line when something was held. */
    void clear (int lev);

    /** The field at this level, or null when it is not held. Always null for ZS and DZS. */
    [[nodiscard]] const amrex::MultiFab* get (int lev, Field id) const
    {
        if (lev < 0 || lev >= num_levels()) { return nullptr; }
        return m_fields[lev][static_cast<int>(id)].get();
    }

    [[nodiscard]] bool has (int lev, Field id) const { return get(lev, id) != nullptr; }

    /** Depth of the center of each soil layer [m], top layer first (ZS); empty when not held. */
    [[nodiscard]] const amrex::Vector<amrex::Real>& soil_depth (int lev) const
    {
        AMREX_ALWAYS_ASSERT(lev >= 0 && lev < num_levels());
        return m_soil_depth[lev];
    }

    /** Thickness of each soil layer [m], top layer first (DZS); empty when not held. */
    [[nodiscard]] const amrex::Vector<amrex::Real>& soil_thickness (int lev) const
    {
        AMREX_ALWAYS_ASSERT(lev >= 0 && lev < num_levels());
        return m_soil_thickness[lev];
    }

    /**
     * Take the soil layer depths (ZS) or thicknesses (DZS) from the column a rank has
     * read. Local: every rank reads the same column.
     *
     * A column that is zero in every layer is taken as not given. Otherwise the depths
     * must be zero or more and increase with the layer (a land model may put its first
     * level at the surface), and the thicknesses must be positive. A column that breaks
     * these rules, or holds a value that is not finite or above the catalog bound, is
     * reported and not held.
     *
     * @param lev level
     * @param id Field::ZS or Field::DZS
     * @param src the column as read; host-readable memory; its z range is the soil layers
     */
    void set_soil_layers (int lev, Field id, const amrex::FArrayBox& src);

    /**
     * Copy the part of one field that one file covers. Every rank calls it, once per file
     * of the level; no MPI call is made here.
     *
     * The valid cells of the level inside file_box take the file's values. A soil field
     * gets one component per layer, in the file's order. The values are counted against
     * the catalog range here; the report, the decision and the ghost cells are left to
     * finish_level().
     *
     * @param lev level
     * @param id field (not ZS or DZS)
     * @param ba2d the level's grids, one cell thick in z
     * @param dm the level's distribution
     * @param ng ghost cells of the stored field (none in z)
     * @param nlayers number of soil layers of the variable in the file (1 for a surface
     *                field); the same on every rank
     * @param src this rank's part of the field as read; memory a device kernel can read
     *            (the reader's fab); its low z index is the first layer
     * @param file_box the columns this file covers at this level (z is ignored)
     */
    void fill (int lev, Field id,
               const amrex::BoxArray& ba2d,
               const amrex::DistributionMapping& dm,
               const amrex::IntVect& ng,
               int nlayers,
               const amrex::FArrayBox& src,
               const amrex::Box& file_box);

    /**
     * Close the read of a level, after the last file. Collective: every rank calls it.
     *
     * Reduces the counts and ranges of fill() over the ranks (two reductions for the
     * level, whatever the number of fields), then for every field read: reports it and
     * fills its ghost cells; or, when a valid cell was not covered by a file, a value
     * lies outside the catalog range, a value is not finite, or two files disagree on
     * the number of layers, reports that with the counts and drops the field. Then drops
     * TSLB, SMOIS, ZS and DZS together when they disagree on the number of soil layers.
     *
     * @param lev level
     * @param period the level's periodicity
     */
    void finish_level (int lev, const amrex::Periodicity& period);

    /**
     * Put the level's fields on new grids (a new level, a regrid or a rebalance).
     * Collective: every rank calls it. Level lev-1 must already be on its new grids.
     *
     * A field the level holds keeps its values on the columns its old grids covered. The
     * other columns take the values of the column of level lev-1 they lie in. When level
     * lev-1 does not hold the field with the same number of components and the new grids
     * reach beyond the old ones, that field is dropped and one line says so.
     *
     * With take_coarser_fields, a field the level does not hold and level lev-1 does is
     * taken from level lev-1 on every column, and so are the soil layer columns when the
     * level has none. This is for a level without a file of its own.
     *
     * The soil layer columns do not depend on the grids and are kept.
     *
     * @param lev level
     * @param new_ba2d,new_dm the new grids, one cell thick in z, and their distribution
     * @param period the level's periodicity
     * @param ratio refinement ratio between level lev-1 and this level (z is ignored;
     *              not used at level 0)
     * @param take_coarser_fields take the fields the level lacks from level lev-1
     * @param verbose above zero: one line per field put on new grids, with the number of
     *                columns taken from level lev-1 and the range
     */
    void remake (int lev,
                 const amrex::BoxArray& new_ba2d,
                 const amrex::DistributionMapping& new_dm,
                 const amrex::Periodicity& period,
                 const amrex::IntVect& ratio,
                 bool take_coarser_fields,
                 int verbose = 0);

    /** Name of the checkpoint's text file that lists what each level holds. */
    static constexpr const char* checkpoint_header_name = "WRFSurfaceFieldsHeader";

    /**
     * Write the held fields into a checkpoint directory: one MultiFab file per field
     * under the level's directory, and one text file that lists what each level holds.
     * The text file is written by the I/O rank, also when nothing is held.
     */
    void write_checkpoint (const std::string& checkpoint_dir, int nlevels) const;

    /**
     * Read what write_checkpoint() wrote, onto the grids of this run (they need not be
     * those of the checkpoint). A checkpoint without the text file (written before this
     * class) leaves the holder empty, prints one line and returns false.
     *
     * Where this run's grids reach beyond those of the checkpoint, the rule of remake()
     * applies: the columns the checkpoint does not give come from level lev-1, or the
     * field is not restored and one line says so. No cell is left without a value.
     *
     * @param ba2d,dm the grids and distribution of each level of this run
     * @param ng ghost cells of the stored fields
     * @param period per level
     * @param ratio refinement ratio between level lev and lev+1, per level (z is ignored)
     */
    bool read_checkpoint (const std::string& checkpoint_dir, int nlevels,
                          const amrex::Vector<amrex::BoxArray>& ba2d,
                          const amrex::Vector<amrex::DistributionMapping>& dm,
                          const amrex::IntVect& ng,
                          const amrex::Vector<amrex::Periodicity>& period,
                          const amrex::Vector<amrex::IntVect>& ratio);

    /** Write and parse the text file; public for the unit tests. */
    void write_header (std::ostream& os, int nlevels) const;
    struct HeaderLevel {
        std::array<int, num_fields> ncomp{};   ///< 0: not held
        amrex::Vector<amrex::Real> soil_depth;
        amrex::Vector<amrex::Real> soil_thickness;
    };
    /** Returns false, leaving out empty, when the text is not a header of this class. */
    [[nodiscard]] static bool parse_header (std::istream& is, amrex::Vector<HeaderLevel>& out);

    /**
     * Fill the ghost cells of a field from its valid cells: a ghost cell that the
     * level's grids cover (also through a periodic boundary) takes the value there; one
     * beyond the bounding box of the level's grids takes the nearest column inside it, as
     * the reader's other surface fields do beyond the file; any other takes the nearest
     * column of its own box. Public for the unit tests.
     */
    static void fill_ghost_cells (amrex::MultiFab& mf, const amrex::Periodicity& period);

private:
    // What fill() has found on this rank for one field, until finish_level() reduces it
    struct Pending {
        amrex::Long n_filled     = 0;   // valid cell values taken from a file
        amrex::Long n_outside    = 0;   // of those, outside the catalog range
        amrex::Long n_not_finite = 0;   // of those, not finite
        amrex::Long n_uncovered  = 0;   // boxes whose part of a file this rank did not read
        bool        read         = false;   // a file of the level held the variable
        bool        conflict     = false;   // two files disagree (layer count, ZS or DZS)
        amrex::Real vmin =  std::numeric_limits<amrex::Real>::max();
        amrex::Real vmax = -std::numeric_limits<amrex::Real>::max();
    };

    void drop (int lev, Field id);

    // The field on the grids ba2d: the valid cells of src (the level's own values on
    // other grids; may be null) where they reach, and the column of level lev-1 it lies
    // in elsewhere. Null when the two together do not give every column.
    [[nodiscard]] std::unique_ptr<amrex::MultiFab>
    on_grids (int lev, Field id, const amrex::MultiFab* src,
              const amrex::BoxArray& ba2d, const amrex::DistributionMapping& dm,
              const amrex::IntVect& ng, const amrex::Periodicity& period,
              const amrex::IntVect& ratio) const;

    // (lev)(field); null when the field is not held
    amrex::Vector<std::array<std::unique_ptr<amrex::MultiFab>, num_fields>> m_fields;
    amrex::Vector<std::array<Pending, num_fields>> m_pending;
    // (lev)(layer)
    amrex::Vector<amrex::Vector<amrex::Real>> m_soil_depth;
    amrex::Vector<amrex::Vector<amrex::Real>> m_soil_thickness;
};

#endif
