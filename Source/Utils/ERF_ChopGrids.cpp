#include <ERF_Utils.H>

using namespace amrex;

/**
 * Decompose the base grids to avoid creating too many grids for the number of processors.
 *
 * @param domain Box specifying the domain to decompose.
 * @param decompose_in_z Whether to decompose in the z-direction.
 * @return BoxArray of the decomposed grids.
 */
BoxArray
ERFPostProcessBaseGrids (const Box& domain, bool decompose_in_z)
{
    //
    // This is used to avoid the case where the native amrex decomposition makes
    //     too many grids for the number of processors.
    //
    // The idea is to not override the user preference if expressed by max_grid_size
    //     but instead to ensure that the default behavior is what we want.
    //
    BoxArray ba0 = amrex::decompose(domain, ParallelDescriptor::NProcs(),
                                    {true,true,decompose_in_z});
    return ba0;
}

/**
 * Iteratively decompose grids in 2D until the target number of grids is reached.
 *
 * @param[in,out] ba BoxArray to be decomposed.
 * @param domain Box specifying the domain.
 * @param target_size Target number of grids.
 */
void
ChopGrids2D (BoxArray& ba, const Box& domain, int target_size)
{
    IntVect chunk = domain.length();

    while (ba.size() < target_size)
    {
        IntVect chunk_prev = chunk;

        // We only decompose in x and y, so this array holds only those two directions;
        //     sizing it with AMREX_SPACEDIM would leave a default {0,0} entry that sorts
        //     to the front and pushes the largest direction out of the loop below
        std::array<std::pair<int,int>,2>
            chunk_dir{std::make_pair(chunk[0],int(0)),
                      std::make_pair(chunk[1],int(1))};
        std::sort(chunk_dir.begin(), chunk_dir.end());

        // Try the largest direction first, then the smaller one
        for (int idx = 1; idx >= 0; idx--) {
            int idim = chunk_dir[idx].second;
            int new_chunk_size = chunk[idim] / 2;
            if (new_chunk_size != 0)
            {
                chunk[idim] = new_chunk_size;
                ba.maxSize(chunk);
                break;
            }
        }

        if (chunk == chunk_prev) {
            break;
        }
    }
}
