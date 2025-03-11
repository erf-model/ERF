#include <ERF_ForestDrag.H>

using namespace amrex;

/*
  Constructor to get the forest parameters:
  TreeType xc, yc, height, diameter, cd, lai, laimax
*/
ForestDrag::ForestDrag (std::string forestfile)
{
    std::ifstream file(forestfile, std::ios::in);
    if (!file.good()) {
        Abort("Cannot find forest file: " + forestfile);
    }
    // TreeType xc yc height diameter cd lai laimax
    Real value1, value2, value3, value4, value5, value6, value7, value8;
    while (file >> value1 >> value2 >> value3 >> value4 >> value5 >> value6 >>
           value7 >> value8) {
        m_type_forest.push_back(value1);
        m_x_forest.push_back(value2);
        m_y_forest.push_back(value3);
        m_height_forest.push_back(value4);
        m_diameter_forest.push_back(value5);
        m_cd_forest.push_back(value6);
        m_lai_forest.push_back(value7);
        m_laimax_forest.push_back(value8);
    }
    file.close();
}

void
ForestDrag::define_drag_field (const BoxArray& ba,
                               const DistributionMapping& dm,
                               Geometry& geom,
                               MultiFab* z_phys_cc,
                               MultiFab* z_phys_nd)
{
    // Geometry params
    const auto& dx = geom.CellSizeArray();
    const auto& prob_lo = geom.ProbLoArray();

    bool all_boxes_touch_bottom = true;
    for (int i = 0; i < ba.size(); i++) {
        if (ba[i].smallEnd(2) != geom.ProbLo(2)) {
            all_boxes_touch_bottom = false;
        }
    }
    AMREX_ALWAYS_ASSERT(all_boxes_touch_bottom);

    // Allocate the forest drag MF
    // NOTE: 1 ghost cell for averaging to faces
    m_forest_drag.reset();
    m_forest_drag = std::make_unique<MultiFab>(ba,dm,1,1);
    m_forest_drag->setVal(0.);

    // Loop over forest types and pre-compute factors
    for (unsigned ii = 0; ii < m_x_forest.size(); ++ii)
    {
        // Expose CPU data for GPU capture
        Real af; // Depends upon the type of forest (tf)
        Real treeZm = 0.0; // Only for forest type 2
        int  tf = int(m_type_forest[ii]);
        Real hf = m_height_forest[ii];
        Real xf = m_x_forest[ii];
        Real yf = m_y_forest[ii];
        Real df = m_diameter_forest[ii];
        Real cdf  = m_cd_forest[ii];
        Real laif = m_lai_forest[ii];
        Real laimaxf = m_laimax_forest[ii];
        if (tf == 1) {
            // Constant factor
            af = laif / hf;
        } else {
            // Discretize integral with 100 points and pre-compute
            int nk      = 100;
            Real ztree  = 0;
            Real expFun = 0;
            Real ratio  = 0;
            const Real dz = hf / Real(nk);
            treeZm = laimaxf * hf;
            for (int k(0); k<nk; ++k) {
                ratio = (hf - treeZm) / (hf - ztree);
                if (ztree < treeZm) {
                    expFun += std::pow(ratio, 6.0) *
                              std::exp(6 * (1 - ratio));
                } else {
                    expFun += std::pow(ratio, 0.5) *
                              std::exp(0.5 * (1 - ratio));
                }
                ztree += dz;
            }
            af = laif / (expFun * dz);
        }

        // Set the forest drag data
        for (MFIter mfi(*m_forest_drag); mfi.isValid(); ++mfi) {
            Box gtbx = mfi.growntilebox();
            const Array4<Real>& levelDrag  = m_forest_drag->array(mfi);
            const Array4<const Real>& z_cc = z_phys_cc->const_array(mfi);
            const Array4<const Real>& z_nd = z_phys_nd->const_array(mfi);

            ParallelFor(gtbx, [=] AMREX_GPU_DEVICE(int i, int j, int k) noexcept
            {
                // Physical positions of cell-centers
                const Real x = prob_lo[0] + (i + 0.5) * dx[0];
                const Real y = prob_lo[1] + (j + 0.5) * dx[1];

                // "z" is measured as distance from cell center to ground
                const Real z_sfc = 0.25 * ( z_nd(i,j  ,0) + z_nd(i+1,j  ,0)
                                           +z_nd(i,j+1,0) + z_nd(i+1,j+1,0));
                const Real z = std::max((z_cc(i,j,k)-z_sfc),0.0);

                // Proximity to the forest
                const Real radius = std::sqrt((x - xf) * (x - xf) +
                                              (y - yf) * (y - yf));

                // Hit for canopy region
                Real factor = 1;
                if ((z <= hf) && (radius <= (0.5 * df))) {
                    if (tf == 2) {
                        Real ratio = (hf - treeZm) / (hf - z);
                        if (z < treeZm) {
                            factor = std::pow(ratio, 6.0) *
                                     std::exp(6.0 * (1.0 - ratio));
                        } else if (z <= hf) {
                            factor = std::pow(ratio, 0.5) *
                                     std::exp(0.5 * (1.0 - ratio));
                        }
                    }
                    levelDrag(i, j, k) = cdf * af * factor;
                }
            });
        } // mfi
    } // ii (forest type)

    // Fillboundary for periodic ghost cell copy
    m_forest_drag->FillBoundary(geom.periodicity());

} // init_drag_field

