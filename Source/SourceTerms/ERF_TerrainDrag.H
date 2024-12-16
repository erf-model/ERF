#ifndef ERF_TERRAINDRAG_H_
#define ERF_TERRAINDRAG_H_

#include <memory>
#include <AMReX_MultiFab.H>

/*
  Adapted from: https://agupubs.onlinelibrary.wiley.com/doi/full/10.1029/2020MS002141
 */
class TerrainDrag
{
public:

    explicit TerrainDrag (std::string terrainfile);

    ~TerrainDrag () = default;

    void
    define_terrain_blank_field (const amrex::BoxArray& ba,
                       const amrex::DistributionMapping& dm,
                       amrex::Geometry& geom,
                       amrex::MultiFab* z_phys_nd);

    amrex::MultiFab*
    get_terrain_blank_field () { return m_terrain_blank.get(); }

private:
    amrex::Vector<amrex::Real> m_x_terrain;
    amrex::Vector<amrex::Real> m_y_terrain;
    amrex::Vector<amrex::Real> m_height_terrain;
    std::unique_ptr<amrex::MultiFab> m_terrain_blank;
};
#endif
