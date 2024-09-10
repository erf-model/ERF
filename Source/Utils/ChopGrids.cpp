#include <Utils.H>

using namespace amrex;

void
ChopGrids2D (BoxArray& ba, const Box& domain, int target_size)
{
    IntVect chunk = domain.length();

    while (ba.size() < target_size)
    {
        IntVect chunk_prev = chunk;

        std::array<std::pair<int,int>,AMREX_SPACEDIM>
            chunk_dir{std::make_pair(chunk[0],int(0)),
                      std::make_pair(chunk[1],int(1))};
        std::sort(chunk_dir.begin(), chunk_dir.end());

        // We only decompose in and y
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
