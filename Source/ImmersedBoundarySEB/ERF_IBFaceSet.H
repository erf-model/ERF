#ifndef ERF_IB_FACE_SET_H
#define ERF_IB_FACE_SET_H

#include <AMReX_MultiFab.H>
#include <AMReX_Geometry.H>
#include <AMReX_GpuContainers.H>
#include <vector>
#include "ERF_IBSEBParams.H"

/**
 * @file ERF_IBFaceSet.H
 * @brief Compact list of the wall faces of resolved buildings on one level
 *
 * A wall face sits between a fluid cell (blanking < 0.5) and a solid cell
 * (blanking >= 0.5). Each face is stored once, on the rank that owns the
 * fluid cell, as a struct of device arrays: the fluid cell (i, j, k), the
 * face direction (0, 1, 2), the side of the fluid cell the solid is on
 * (-1: at i-1, +1: at i+1, and likewise for j and k), the building id, the
 * material id, the face area, the skin temperature, the slab temperatures,
 * the view fractions and the current fluxes. Faces are contiguous per local
 * fab in MFIter order, so a kernel over one fab's faces can write into that
 * fab's arrays.
 *
 * Every input of the balance is local to the fluid cell, so the per-step
 * update needs no communication; only the reports reduce across ranks.
 *
 * Building ids come from labelling the solid columns of the blanking
 * (4-connected in the horizontal), numbered in scan order; a face carries
 * the id of the column of its solid neighbour.
 *
 * Output goes through cell-centred fields on demand: scatter_diagnostics()
 * for the plotfile, save_state() / load_state() for the checkpoint, which
 * make restarts independent of the rank count because the list is rebuilt
 * from the blanking and refilled from the field.
 */
class IBFaceSet
{
public:
    IBFaceSet (const IBSEBParams& params, int lev) : m_params(params), m_lev(lev) {}

    /// Detect the faces of this rank from the cell-centred blanking (one ghost cell needed).
    void build (const amrex::MultiFab& blanking, const amrex::Geometry& geom);

    int n_faces () const { return m_nface; }
    int n_buildings () const { return m_nbld; }
    int n_layers () const { return m_params.n_slab_layers; }

    /// Components of the checkpoint state field: 6 face slots x (skin + layers).
    int state_ncomp () const { return 6 * (1 + n_layers()); }

    /// Number of faces touching each cell and their mean skin temperature (cell-centred, 1 comp each).
    void scatter_diagnostics (amrex::MultiFab& nfaces, amrex::MultiFab& tskin) const;

    /// Write / read the face state through a cell-centred field with state_ncomp() components.
    void save_state (amrex::MultiFab& state) const;
    void load_state (const amrex::MultiFab& state);

    /// Print the summary line and append one row per building to the CSV.
    void report (amrex::Real time, int step, bool write_csv) const;

    // Device arrays; the kernels of the later phases read and write them.
    amrex::Gpu::DeviceVector<int>         d_i, d_j, d_k, d_dir, d_side, d_bid, d_mat;
    amrex::Gpu::DeviceVector<amrex::Real> d_area, d_T_skin, d_T_slab;
    amrex::Gpu::DeviceVector<amrex::Real> d_f_sky, d_f_ground, d_f_bldg;
    amrex::Gpu::DeviceVector<amrex::Real> d_SW_abs, d_LW_net, d_H, d_LE, d_G, d_Q_ext;

    /// First face of each local fab (MFIter order); size nfabs + 1.
    const std::vector<int>& fab_start () const { return m_fab_start; }

private:
    IBSEBParams m_params;
    int m_lev   = 0;
    int m_nface = 0;
    int m_nbld  = 0;
    std::vector<int> m_fab_start;

    // Static per-building totals over all ranks, index 0 unused.
    std::vector<long>        m_bld_nface;
    std::vector<amrex::Real> m_bld_area;
    long m_nface_dir[3] = {0, 0, 0};
    amrex::Real m_area_total = 0.0;
};

#endif
