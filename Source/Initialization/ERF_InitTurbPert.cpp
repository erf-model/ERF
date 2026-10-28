/**
 * \file ERF_InitTurbPert.cpp
 */
#include <ERF.H>
#include <ERF_TerrainSurfaceSlab.H>
#include <AMReX_MultiFabUtil.H>
#include <ERF_TileNoZ.H>
#include <ERF_ProbCommon.H>

using namespace amrex;

namespace {

// A bottom nodal slab for every state box, including boxes stacked in z. This
// keeps terrain AGL calculations independent of the grid decomposition.
std::unique_ptr<MultiFab>
perturbation_ground (const MultiFab* z_nd, const MultiFab& state, const Geometry& geom)
{
    if (!z_nd) { return nullptr; }
    BoxArray ba = bottom_node_slab(state.boxArray(), geom);
    auto ground = std::make_unique<MultiFab>(ba, state.DistributionMap(), 1, 0);
    ground->setVal(bogus_large_value);
    ground->ParallelCopy(*z_nd, 0, 0, 1, 0, 0, geom.periodicity());
    return ground;
}

} // namespace

// Resolve each rectangular nest's lateral inflow edges using the latest
// available parent winds and the finest PBL ancestor's height. Run after the
// ancestors advance, so diagnostics exist on the first step of a new or restarted nest.
void
ERF::turbPert_update_automatic (int lev, double time)
{
    if (!turbPert.automatic_update_due(lev, time)) { return; }
    BL_PROFILE("ERF::turbPert_update_automatic()");
    AMREX_ALWAYS_ASSERT(lev > 0 && solverChoice.perturbation_auto);

    const int parent = lev-1;
    const auto& parent_state = vars_new[parent][Vars::cons];
    const auto& parent_u = vars_new[parent][Vars::xvel];
    const auto& parent_v = vars_new[parent][Vars::yvel];
    const Box& parent_domain = geom[parent].Domain();
    const Box& domain = geom[lev].Domain();
    const Real fraction = solverChoice.perturbation_auto_pbl_fraction;

    int pbl_level = parent;
    while (pbl_level >= 0 && solverChoice.turbChoice[pbl_level].pbl_type == PBLType::None) {
        --pbl_level;
    }
    if (pbl_level < 0) {
        Abort("erf.perturbation_auto requires a coarser level with a PBL scheme "
              "to supply the PBL height for level " + std::to_string(lev));
    }

    const MultiFab* pblh = nullptr;
    if (solverChoice.turbChoice[pbl_level].uses_native_shoc() && native_shoc_driver[pbl_level] &&
        native_shoc_driver[pbl_level]->has_native_diagnostics()) {
        pblh = &native_shoc_driver[pbl_level]->pblh_diagnostics();
    } else if (m_SurfaceLayer[Orientation::zlo()]) {
        // YSU/MRF supply their scheme's height via set_pblh, including when the
        // separate MOST estimator is disabled. MYNN uses most.pblh_calc.
        pblh = m_SurfaceLayer[Orientation::zlo()]->get_pblh(pbl_level);
    }
    if (!pblh) {
        Abort("erf.perturbation_auto requires a PBL-height diagnostic on PBL level " +
              std::to_string(pbl_level) + ": use a zlo surface_layer with "
              "erf.most.pblh_calc = MYNN25, YSU/MRF, or native SHOC");
    }

    // Gather the ancestor's column heights onto the parent layout. Each parent
    // column takes the diagnostic of the ancestor cell containing it. Coarsen
    // only in x/y: PBLH is a depth in metres AGL, independent of either level's
    // vertical spacing, domain height, and decomposition into stacked boxes.
    IntVect height_ratio(1);
    for (int level = pbl_level; level < parent; ++level) {
        height_ratio[0] *= refRatio(level)[0];
        height_ratio[1] *= refRatio(level)[1];
    }
    const int hklo = geom[pbl_level].Domain().smallEnd(2);
    BoxArray height_boxes = amrex::coarsen(parent_state.boxArray(), height_ratio);
    for (int b = 0; b < height_boxes.size(); ++b) {
        Box slab = height_boxes[b];
        slab.setRange(2, hklo);
        height_boxes.set(b, slab);
    }
    MultiFab sampled_pblh(height_boxes, parent_state.DistributionMap(), 1, 0);
    sampled_pblh.setVal(bogus_large_value);
    // Both the planar surface-layer diagnostic and native SHOC's 3D diagnostic
    // contain their column height at the source domain's bottom cell.
    sampled_pblh.ParallelCopy(*pblh, 0, 0, 1, 0, 0, geom[pbl_level].periodicity());
    const int hrx = height_ratio[0];
    const int hry = height_ratio[1];

    Vector<AutomaticPerturbationFace> faces;
    Vector<Box> sample_boxes;
    Vector<int> subdomain_ids;
    for (int isub = 0; isub < subdomains[lev].size(); ++isub) {
        const Box nest = subdomains[lev][isub].minimalBox();
        if (nest.smallEnd(2) != domain.smallEnd(2)) {
            Abort("Automatic perturbations require nests that reach the ground (level " +
                  std::to_string(lev) + ", subdomain " + std::to_string(isub) + ")");
        }
        for (int face = 0; face < 4; ++face) {
            const int dir = face % 2;
            const bool low = face < 2;
            const bool at_domain_edge = low ? nest.smallEnd(dir) == domain.smallEnd(dir)
                                            : nest.bigEnd(dir) == domain.bigEnd(dir);
            if (at_domain_edge) {
                const auto bc = phys_bc_type[Orientation(dir, low ? Orientation::low : Orientation::high)];
                if (geom[lev].isPeriodic(dir) || bc == ERF_BC::slip_wall ||
                    bc == ERF_BC::no_slip_wall || bc == ERF_BC::symmetry ||
                    bc == ERF_BC::surface_layer) { continue; }
            }
            AutomaticPerturbationFace control;
            control.region = turbPert.automatic_face_box(nest, face);
            control.face = face;
            if (!control.region.ok()) { continue; }
            faces.push_back(control);
            subdomain_ids.push_back(isub);

            // One parent column thick along this edge; use the full parent
            // depth to diagnose the winds in its boundary layer.
            Box sample = amrex::coarsen(nest, refRatio(parent));
            sample.setRange(dir, low ? sample.smallEnd(dir) : sample.bigEnd(dir));
            sample.setSmall(2, parent_domain.smallEnd(2));
            sample.setBig(2, parent_domain.bigEnd(2));
            sample_boxes.push_back(sample);
        }
    }

    auto parent_ground = perturbation_ground(z_phys_nd[parent].get(), parent_state, geom[parent]);
    const bool terrain = (z_phys_nd[parent] != nullptr);
    const auto pdx = geom[parent].CellSizeArray();
    const int pklo = parent_domain.smallEnd(2);
    // Per edge: integral u dz, integral v dz, depth, sum PBLH, column count,
    // and invalid diagnostic count. MPI sums keep the decision identical on all ranks.
    constexpr int nsum = 6;
    Vector<Real> sums(nsum*faces.size(), zero);
    for (MFIter mfi(parent_state); mfi.isValid(); ++mfi) {
        const auto ua = parent_u.const_array(mfi);
        const auto va = parent_v.const_array(mfi);
        const auto ha = sampled_pblh.const_array(mfi);
        const auto nd = terrain ? z_phys_nd[parent]->const_array(mfi) : Array4<const Real>{};
        const auto ground = terrain ? parent_ground->const_array(mfi) : Array4<const Real>{};
        for (int f = 0; f < faces.size(); ++f) {
            const Box bx = sample_boxes[f] & mfi.validbox();
            if (!bx.ok()) { continue; }
            ReduceOps<ReduceOpSum, ReduceOpSum, ReduceOpSum, ReduceOpSum, ReduceOpSum, ReduceOpSum> op;
            ReduceData<Real, Real, Real, Real, Real, Real> data(op);
            using Tuple = typename decltype(data)::Type;
            op.eval(bx, data, [=] AMREX_GPU_DEVICE (int i, int j, int k) noexcept -> Tuple {
                const Real h = ha(amrex::coarsen(i,hrx), amrex::coarsen(j,hry), hklo);
                const bool valid = h >= zero && h < amrex::Real(1.e10);
                const Real bottom = (k == pklo) ? one : zero;
                if (!valid) { return {zero, zero, zero, zero, bottom, one}; }
                Real zlo = (k-pklo)*pdx[2];
                Real zhi = zlo + pdx[2];
                if (terrain) {
                    const Real zs = amrex::Real(0.25)*(ground(i,j,pklo) + ground(i+1,j,pklo) +
                                                       ground(i,j+1,pklo) + ground(i+1,j+1,pklo));
                    zlo = amrex::Real(0.25)*(nd(i,j,k) + nd(i+1,j,k) + nd(i,j+1,k) + nd(i+1,j+1,k)) - zs;
                    zhi = amrex::Real(0.25)*(nd(i,j,k+1) + nd(i+1,j,k+1) + nd(i,j+1,k+1) + nd(i+1,j+1,k+1)) - zs;
                }
                const Real weight = amrex::max(zero, amrex::min(zhi, fraction*h) - amrex::max(zero, zlo));
                const Real u = amrex::Real(0.5)*(ua(i,j,k) + ua(i+1,j,k));
                const Real v = amrex::Real(0.5)*(va(i,j,k) + va(i,j+1,k));
                return {u*weight, v*weight, weight, bottom*h, bottom, zero};
            });
            const auto r = data.value(op);
            sums[nsum*f  ] += amrex::get<0>(r);
            sums[nsum*f+1] += amrex::get<1>(r);
            sums[nsum*f+2] += amrex::get<2>(r);
            sums[nsum*f+3] += amrex::get<3>(r);
            sums[nsum*f+4] += amrex::get<4>(r);
            sums[nsum*f+5] += amrex::get<5>(r);
        }
    }
    if (!sums.empty()) {
        ParallelAllReduce::Sum(sums.data(), static_cast<int>(sums.size()), ParallelContext::CommunicatorSub());
    }
    static constexpr const char* names[] = {"xlo", "ylo", "xhi", "yhi"};
    for (int f = 0; f < faces.size(); ++f) {
        Box footprint = sample_boxes[f];
        footprint.setRange(2, pklo);
        const Real count = sums[nsum*f+4];
        if (sums[nsum*f+5] > zero || count != Real(footprint.numPts())) {
            Abort("Automatic perturbations: missing/invalid PBL height from PBL level " +
                  std::to_string(pbl_level) + " on parent level " + std::to_string(parent) +
                  " along " + names[faces[f].face] + ". Check ancestor coverage and "
                  "erf.most.pblh_calc; the height must be diagnosed before the nest advances.");
        }
        const Real depth = sums[nsum*f+2];
        faces[f].u = depth > zero ? sums[nsum*f] / depth : zero;
        faces[f].v = depth > zero ? sums[nsum*f+1] / depth : zero;
        // One mean AGL cutoff for the entire edge, including all perturbation
        // layers. Column-local heights above only set the wind integration depth.
        faces[f].height = fraction * sums[nsum*f+3] / count;
        if (verbose > 0) {
            Print() << "Automatic perturbations: time=" << time << " level=" << lev
                    << " parent=" << parent << " pbl_level=" << pbl_level
                    << " subdomain=" << subdomain_ids[f]
                    << " edge=" << names[faces[f].face] << " inflow=" << faces[f].inflow()
                    << " u=" << faces[f].u << " v=" << faces[f].v
                    << " height_agl=" << faces[f].height << " m\n";
        }
    }

    auto fine_ground = perturbation_ground(z_phys_nd[lev].get(), vars_old[lev][Vars::cons], geom[lev]);
    turbPert.configure_automatic(lev, time, solverChoice.perturbation_auto_period, faces, geom[lev],
                                 z_phys_cc[lev].get(), z_phys_nd[lev].get(), fine_ground.get());
    // The mask kernels read the temporary ground slabs.
    Gpu::streamSynchronize();
}

/**
 * Update turbulent perturbation timing and amplitudes for one AMR level.
 *
 * @param lev Integer specifying the current level
 * @param local_dt Local timestep used to update the perturbation state
 */
void
ERF::turbPert_update (const int lev, const double local_dt)
{
    // Accessing data
    auto& lev_new = vars_new[lev];

    // Create aliases to state data to pass to calc_tpi_update
    int ncons = lev_new[Vars::cons].nComp();
    MultiFab cons_data(lev_new[Vars::cons], make_alias, 0, ncons);
    MultiFab xvel_data(lev_new[Vars::xvel], make_alias, 0, 1);
    MultiFab yvel_data(lev_new[Vars::yvel], make_alias, 0, 1);

    // Computing perturbation update time
    turbPert.calc_tpi_update(lev, local_dt, xvel_data, yvel_data, cons_data);

    Print() << "Successfully initialized turbulent perturbation update time and amplitude at level "
            << lev << " with type: " << turbPert.pt_type[lev] << "\n";
}

/**
 * Calculate the cell-centered turbulent perturbation amplitude.
 *
 * @param lev Integer specifying the current level
 */
void
ERF::turbPert_amplitude (int lev)
{
    // Accessing data
    auto& lev_new = vars_new[lev];

    // Creating local data
    int ncons = lev_new[Vars::cons].nComp();
    MultiFab cons_data(lev_new[Vars::cons], make_alias, 0, ncons);

    // Defining BoxArray type
    auto m_ixtype = cons_data.boxArray().ixType();

#ifdef _OPENMP
#pragma omp parallel if (amrex::Gpu::notInLaunchRegion())
#endif
    for (MFIter mfi(lev_new[Vars::cons], TileNoZ()); mfi.isValid(); ++mfi) {
        // Each tile adds the perturbation over its own cells only
        const Box &bx  = mfi.tilebox();
        const auto &cons_pert_arr = cons_data.array(mfi); // Address of perturbation array
        const amrex::Array4<const amrex::Real> &pert_cell = turbPert.pb_cell[lev].array(mfi); // per-cell perturbation stored in structure

        turbPert.apply_tpi(lev, bx, RhoTheta_comp, m_ixtype, cons_pert_arr, pert_cell);
    } // mfi
}

/**
 * Calculate the turbulent perturbation amplitude for vertical velocity.
 *
 * @param lev Integer specifying the current level
 */
void
ERF::turbPert_amplitude_w (int lev)
{
    // Accessing data
    auto& lev_new = vars_new[lev];

    // Creating local data
    int ncomps = lev_new[Vars::zvel].nComp();
    MultiFab wvel_data(lev_new[Vars::zvel], make_alias, 0, ncomps);

    // Defining BoxArray type
    auto m_ixtype = wvel_data.boxArray().ixType();

#ifdef _OPENMP
#pragma omp parallel if (amrex::Gpu::notInLaunchRegion())
#endif
    for (MFIter mfi(lev_new[Vars::zvel], TileNoZ()); mfi.isValid(); ++mfi) {
        // Each tile adds the perturbation over its own faces only
        const Box &bx  = mfi.tilebox();
        const auto &wvel_pert_arr = wvel_data.array(mfi); // Address of perturbation array
        const amrex::Array4<const amrex::Real> &pert_cell = turbPert.pb_cell[lev].array(mfi); // per-cell perturbation stored in structure

        turbPert.apply_tpi(lev, bx, -1, m_ixtype, wvel_pert_arr, pert_cell);
    } // mfi
}
