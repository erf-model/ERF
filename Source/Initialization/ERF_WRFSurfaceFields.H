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
 * The data are those of the initial time. They are not checkpointed, they are dropped when
 * the level's grids change, and a level made by tagging (no file of its own) has none.
 * Code that needs them after the start must copy what it needs into its own state.
 *
 * The class has no NetCDF dependency: fill() takes the region a rank has read as an
 * FArrayBox, so the unit tests can drive it without a file.
 */

#include <array>
#include <cmath>
#include <limits>
#include <memory>
#include <string>

#include <AMReX_BoxArray.H>
#include <AMReX_DistributionMapping.H>
#include <AMReX_FArrayBox.H>
#include <AMReX_Gpu.H>
#include <AMReX_MultiFab.H>
#include <AMReX_ParallelDescriptor.H>
#include <AMReX_Periodicity.H>
#include <AMReX_Print.H>
#include <AMReX_Vector.H>

class WRFSurfaceFields
{
public:
    /** The fields, in the order they are requested from the file. */
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
        NumFields
    };

    static constexpr int num_fields = static_cast<int>(Field::NumFields);

    /** A field's name in the file, its units, and the range outside which the read stops. */
    struct Descriptor {
        Field       id;
        const char* nc_name;
        const char* units;
        amrex::Real lo;    ///< smallest value a file may hold
        amrex::Real hi;    ///< largest value a file may hold
        bool        soil;  ///< true: one value per soil layer
    };

    /**
     * The bounds stop only values that cannot be the named quantity (a negative
     * temperature, a fraction above one, a NaN). A field that a hand-built file declares
     * and leaves at zero passes; whether zero is usable is for the code that takes the field.
     * The upper bounds of LAI, TMN, ZNT, SNOW, SNOWH and TSLB are not physical limits: they
     * are set far above any value of the quantity so that a NetCDF fill value (9.97e36),
     * which is finite, is stopped.
     */
    [[nodiscard]] static const std::array<Descriptor, num_fields>& catalog ()
    {
        // Not named zero and one: those are constants of ERF_NumericalConstants.H
        constexpr amrex::Real lo0 = amrex::Real(0.0);
        constexpr amrex::Real hi1 = amrex::Real(1.0);
        static const std::array<Descriptor, num_fields> table{{
            {Field::LAI,    "LAI",    "m2/m2",  lo0,  amrex::Real(100.0),  false},
            {Field::VEGFRA, "VEGFRA", "%",      lo0,  amrex::Real(100.0),  false},
            {Field::TMN,    "TMN",    "K",      lo0,  amrex::Real(1000.0), false},
            {Field::ZNT,    "ZNT",    "m",      lo0,  amrex::Real(100.0),  false},
            {Field::ALBEDO, "ALBEDO", "1",      lo0,  hi1,                 false},
            {Field::ALBBCK, "ALBBCK", "1",      lo0,  hi1,                 false},
            {Field::SEAICE, "SEAICE", "1",      lo0,  hi1,                 false},
            {Field::SNOW,   "SNOW",   "kg/m2",  lo0,  amrex::Real(1.0e6),  false},
            {Field::SNOWH,  "SNOWH",  "m",      lo0,  amrex::Real(1.0e4),  false},
            {Field::TSLB,   "TSLB",   "K",      lo0,  amrex::Real(1000.0), true},
            {Field::SMOIS,  "SMOIS",  "m3/m3",  lo0,  hi1,                 true},
        }};
        return table;
    }

    /** The descriptor of the field with this name in the file, or null when it is not one of these. */
    [[nodiscard]] static const Descriptor* find (const std::string& nc_name)
    {
        for (const auto& d : catalog()) {
            if (nc_name == d.nc_name) { return &d; }
        }
        return nullptr;
    }

    void resize (int nlevs)
    {
        m_fields.resize(nlevs);
        m_soil_depth.resize(nlevs);
        m_soil_thickness.resize(nlevs);
    }

    [[nodiscard]] int num_levels () const { return static_cast<int>(m_fields.size()); }

    /** Drop everything held at this level (its grids are about to change). */
    void clear (int lev)
    {
        if (lev < 0 || lev >= num_levels()) { return; }
        for (auto& f : m_fields[lev]) { f.reset(); }
        m_soil_depth[lev].clear();
        m_soil_thickness[lev].clear();
    }

    /** The field at this level, or null when the level's file did not hold it. */
    [[nodiscard]] const amrex::MultiFab* get (int lev, Field id) const
    {
        if (lev < 0 || lev >= num_levels()) { return nullptr; }
        return m_fields[lev][static_cast<int>(id)].get();
    }

    [[nodiscard]] bool has (int lev, Field id) const { return get(lev, id) != nullptr; }

    /** Depth of the centre of each soil layer [m], top layer first (ZS); empty when not read. */
    [[nodiscard]] const amrex::Vector<amrex::Real>& soil_depth (int lev) const
    {
        AMREX_ALWAYS_ASSERT(lev >= 0 && lev < num_levels());
        return m_soil_depth[lev];
    }

    /** Thickness of each soil layer [m], top layer first (DZS); empty when not read. */
    [[nodiscard]] const amrex::Vector<amrex::Real>& soil_thickness (int lev) const
    {
        AMREX_ALWAYS_ASSERT(lev >= 0 && lev < num_levels());
        return m_soil_thickness[lev];
    }

    /**
     * Take the soil layer depths or thicknesses from the column a rank has read.
     *
     * @param lev level
     * @param thickness true for DZS, false for ZS
     * @param src the column as read; host memory; its z range is the soil layers
     *
     * A column that is zero in every layer is taken as not given and leaves the vector
     * empty. Otherwise every value must be positive, and the depths must increase.
     */
    void set_soil_layers (int lev, bool thickness, const amrex::FArrayBox& src)
    {
        AMREX_ALWAYS_ASSERT(lev >= 0 && lev < num_levels());
        const amrex::Box& sbx = src.box();
        AMREX_ALWAYS_ASSERT_WITH_MESSAGE(sbx.ok(),
            "WRFSurfaceFields::set_soil_layers: this rank holds no part of the soil layer column");
        auto& dst = thickness ? m_soil_thickness[lev] : m_soil_depth[lev];
        const char* name = thickness ? "DZS" : "ZS";
        const auto arr = src.const_array();
        const int ilo = sbx.smallEnd(0);
        const int jlo = sbx.smallEnd(1);
        const int nlay = sbx.length(2);
        dst.resize(nlay);
        bool all_zero = true;
        for (int n = 0; n < nlay; ++n) {
            dst[n] = arr(ilo, jlo, sbx.smallEnd(2) + n);
            // isfinite first: a comparison with a NaN raises the invalid-operation exception
            AMREX_ALWAYS_ASSERT_WITH_MESSAGE(std::isfinite(dst[n]),
                thickness ? "wrfinput DZS holds a value that is not finite"
                          : "wrfinput ZS holds a value that is not finite");
            all_zero = all_zero && (dst[n] == amrex::Real(0.0));
        }
        // A file written without a land model (ideal.exe, or a hand-built file) declares
        // the layers and leaves them at zero: that is a file without them.
        if (all_zero) {
            dst.clear();
            amrex::Print() << "   wrfinput " << name << " at level " << lev
                           << " is zero in every soil layer: taken as not given\n";
            return;
        }
        for (int n = 0; n < nlay; ++n) {
            AMREX_ALWAYS_ASSERT_WITH_MESSAGE(dst[n] > amrex::Real(0.0),
                thickness ? "wrfinput DZS: every soil layer thickness must be positive"
                          : "wrfinput ZS: every soil layer depth must be positive");
            // The top layer is first: the depths grow with the layer index
            AMREX_ALWAYS_ASSERT_WITH_MESSAGE(thickness || n == 0 || dst[n] > dst[n-1],
                "wrfinput ZS: the soil layer depths must increase, top layer first");
        }
        amrex::Print() << "   wrfinput soil layer " << (thickness ? "thicknesses DZS" : "depths ZS")
                       << " at level " << lev << ":";
        for (const amrex::Real v : dst) { amrex::Print() << " " << v; }
        amrex::Print() << " [m]\n";
    }

    /**
     * Stop the run when the soil data of a level do not agree on the number of layers:
     * TSLB, SMOIS, ZS and DZS, whichever of them were read. Local; every rank holds the
     * same counts.
     */
    void check_soil_layers (int lev) const
    {
        AMREX_ALWAYS_ASSERT(lev >= 0 && lev < num_levels());
        int nlay = 0;
        const auto agree = [&nlay] (int n) {
            if (n == 0) { return true; }   // not read
            if (nlay == 0) { nlay = n; }
            return n == nlay;
        };
        bool ok = agree(static_cast<int>(m_soil_depth[lev].size()));
        ok = agree(static_cast<int>(m_soil_thickness[lev].size())) && ok;
        for (const Field id : {Field::TSLB, Field::SMOIS}) {
            const amrex::MultiFab* mf = get(lev, id);
            ok = agree(mf ? mf->nComp() : 0) && ok;
        }
        AMREX_ALWAYS_ASSERT_WITH_MESSAGE(ok,
            "wrfinput TSLB, SMOIS, ZS and DZS must have the same number of soil layers");
    }

    /**
     * Copy one field from the region this rank has read. Collective: every rank calls it.
     *
     * A cell outside [i_lo,i_hi] x [j_lo,j_hi] (a ghost cell beyond the file) takes the
     * value of the nearest column inside, as the reader's other surface fields do. A soil
     * field gets one component per layer of the source's z range, in the file's order.
     * The values are checked against the catalog's bounds over the valid cells and the run
     * stops, naming the variable, when one lies outside or is not finite. The run also
     * stops when the source does not cover the rank's columns. Both checks are reduced
     * over the ranks first, so every rank stops together.
     *
     * Three reductions per call: the component count, the three counts of the checks, and
     * the range.
     *
     * @param lev level
     * @param id field
     * @param ba2d the level's grids, one cell thick in z
     * @param dm the level's distribution
     * @param ng ghost cells of the stored field (none in z)
     * @param period the level's periodicity
     * @param src this rank's part of the field; host memory; must cover the rank's
     *            columns in [i_lo,i_hi] x [j_lo,j_hi], ghost columns included
     * @param i_lo,i_hi,j_lo,j_hi the columns the file covers at this level
     */
    void fill (int lev, Field id,
               const amrex::BoxArray& ba2d,
               const amrex::DistributionMapping& dm,
               const amrex::IntVect& ng,
               const amrex::Periodicity& period,
               const amrex::FArrayBox& src,
               int i_lo, int i_hi, int j_lo, int j_hi)
    {
        AMREX_ALWAYS_ASSERT(lev >= 0 && lev < num_levels());
        const Descriptor& d = catalog()[static_cast<int>(id)];
        const amrex::Box& sbx = src.box();
        const int klo   = sbx.ok() ? sbx.smallEnd(2) : 0;
        int ncomp = (d.soil && sbx.ok()) ? sbx.length(2) : 1;
        // A rank without boxes at this level reads one column only; all ranks must agree
        // on the number of layers before the MultiFab is defined.
        amrex::ParallelDescriptor::ReduceIntMax(ncomp);

        auto& mf = m_fields[lev][static_cast<int>(id)];
        mf = std::make_unique<amrex::MultiFab>(ba2d, dm, ncomp, ng);

        amrex::Real vmin =  std::numeric_limits<amrex::Real>::max();
        amrex::Real vmax = -std::numeric_limits<amrex::Real>::max();
        // Valid cells outside the bounds, valid cells that are not finite, and cells
        // (ghost cells included) the source does not cover
        amrex::Long n_outside = 0;
        amrex::Long n_not_finite = 0;
        amrex::Long n_uncovered = 0;

        // The source is in host memory, so the copy is made on the host into a buffer
        // that is then sent to the field's own memory.
        for (amrex::MFIter mfi(*mf, false); mfi.isValid(); ++mfi) {
            const amrex::Box gbx = mfi.fabbox();
            const amrex::Box vbx = mfi.validbox();
            amrex::FArrayBox buf(gbx, ncomp, amrex::The_Pinned_Arena());
            const auto b = buf.array();
            const auto s = src.const_array();
            const amrex::Dim3 lo = amrex::lbound(gbx);
            const amrex::Dim3 hi = amrex::ubound(gbx);
            for (int n = 0; n < ncomp; ++n) {
                for (int j = lo.y; j <= hi.y; ++j) {
                    for (int i = lo.x; i <= hi.x; ++i) {
                        const int li = amrex::min(amrex::max(i, i_lo), i_hi);
                        const int lj = amrex::min(amrex::max(j, j_lo), j_hi);
                        if (!sbx.contains(amrex::IntVect(li, lj, klo + n))) {
                            // Counted and reported below; the source is not read here
                            ++n_uncovered;
                            b(i, j, lo.z, n) = amrex::Real(0.0);
                            continue;
                        }
                        const amrex::Real v = s(li, lj, klo + n);
                        b(i, j, lo.z, n) = v;
                        if (vbx.contains(amrex::IntVect(i, j, lo.z))) {
                            // isfinite first: a comparison with a NaN raises the
                            // invalid-operation exception, which a run may trap
                            if (!std::isfinite(v)) {
                                ++n_not_finite;
                            } else {
                                if (v < d.lo || v > d.hi) { ++n_outside; }
                                vmin = amrex::min(vmin, v);
                                vmax = amrex::max(vmax, v);
                            }
                        }
                    }
                }
            }
            amrex::Gpu::htod_memcpy_async((*mf)[mfi].dataPtr(), buf.dataPtr(), buf.nBytes());
            // The buffer leaves scope at the end of this pass
            amrex::Gpu::streamSynchronize();
        }

        amrex::Long counts[3] = {n_outside, n_not_finite, n_uncovered};
        amrex::ParallelDescriptor::ReduceLongSum(counts, 3);
        amrex::Real range[2] = {-vmin, vmax};
        amrex::ParallelDescriptor::ReduceRealMax(range, 2);
        vmin = -range[0];
        vmax =  range[1];

        if (counts[2] > 0) {
            amrex::Abort(std::string("wrfinput ") + d.nc_name +
                         ": the region read does not cover the level's columns");
        }
        if (counts[0] > 0 || counts[1] > 0) {
            amrex::Print() << "wrfinput " << d.nc_name << " at level " << lev << ": "
                           << counts[0] << " value(s) outside " << d.lo << " to " << d.hi
                           << " [" << d.units << "], " << counts[1] << " not finite";
            if (vmin <= vmax) {
                amrex::Print() << "; the finite values span " << vmin << " to " << vmax;
            }
            amrex::Print() << "\n";
            amrex::Abort(std::string("wrfinput ") + d.nc_name +
                         " holds a value outside its allowed range, or one that is not finite");
        }
        amrex::Print() << "   wrfinput surface field " << d.nc_name << " at level " << lev << ": "
                       << ncomp << (ncomp == 1 ? " component" : " components")
                       << ", min " << vmin << ", max " << vmax << " [" << d.units << "]\n";

        mf->FillBoundary(period);
    }

private:
    // (lev)(field); null when the level's file did not hold the field
    amrex::Vector<std::array<std::unique_ptr<amrex::MultiFab>, num_fields>> m_fields;
    // (lev)(layer)
    amrex::Vector<amrex::Vector<amrex::Real>> m_soil_depth;
    amrex::Vector<amrex::Vector<amrex::Real>> m_soil_thickness;
};

#endif
