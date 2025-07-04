#include "ERF.H"

using namespace amrex;

void
ERF::make_subdomains(const BoxList& bl, Vector<BoxArray>& bins)
{
    Vector<BoxList> bins_bl;

    // Clear out any old bins
    bins.clear();

    // Iterate over boxes
    for (auto bx : bl)
    {
        bool added = false;

        // Try to add box to existing bin
        for (int j = 0; j < bins_bl.size(); ++j) {
            BoxList& bin = bins_bl[j];
            bool touches = false;

            for (auto& b : bin)
            {
                Box gbx(bx); gbx.grow(1);
                if (gbx.intersects(b)) {
                    touches = true;
                    break;
                }
            }

            if (touches) {
                bin.push_back(bx);
                added = true;
                break;
            }
        }

        // If box couldn't be added to existing bin, create new bin
        if (!added) {
            BoxList new_bin;
            new_bin.push_back(bx);
            bins_bl.push_back(new_bin);
        }
    }

    // Convert the BoxLists to BoxArrays
    for (int i = 0; i < bins_bl.size(); ++i) {
        bins.push_back(BoxArray(bins_bl[i]));
    }
}
