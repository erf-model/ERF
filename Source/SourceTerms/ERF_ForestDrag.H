#ifndef ERF_FORESTDRAG_H_
#define ERF_FORESTDRAG_H_

#include <memory>
#include <AMReX_MultiFab.H>

/*
  ForestDrag flow physics adapted from:
  Lalic & Mihailovic (2004)
  https://doi.org/10.1175/1520-0450(2004)043<0641:AERDLD>2.0.CO;2
 */
class ForestDrag
{
public:

    explicit ForestDrag (std::string forestfile);

    ~ForestDrag () = default;

    void
    define_drag_field (const amrex::BoxArray& ba,
                       const amrex::DistributionMapping& dm,
                       amrex::Geometry& geom,
                       amrex::MultiFab* z_phys_cc,
                       amrex::MultiFab* z_phys_nd);

    amrex::MultiFab*
    get_drag_field () { return m_forest_drag.get(); }

private:
    amrex::Vector<amrex::Real> m_type_forest;
    amrex::Vector<amrex::Real> m_x_forest;
    amrex::Vector<amrex::Real> m_y_forest;
    amrex::Vector<amrex::Real> m_height_forest;
    amrex::Vector<amrex::Real> m_diameter_forest;
    amrex::Vector<amrex::Real> m_cd_forest;
    amrex::Vector<amrex::Real> m_lai_forest;
    amrex::Vector<amrex::Real> m_laimax_forest;
    std::unique_ptr<amrex::MultiFab> m_forest_drag;
};
#endif
