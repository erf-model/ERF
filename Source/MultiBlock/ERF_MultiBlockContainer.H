#ifndef ERF_MultiBlockContainer_H_
#define ERF_MultiBlockContainer_H_

#include <AMReX_NonLocalBC.H>
#include <ERF.H>

class MultiBlockContainer
{
public:

    // Constructor with vector of parameters
    MultiBlockContainer (const std::vector<amrex::RealBox>& rb_v,
                         std::vector<int> max_level_in_v,
                         const std::vector<amrex::Vector<int>>& n_cell_in_v,
                         std::vector<int> coord_v,
                         const std::vector<amrex::Vector<amrex::IntVect>>& ref_ratios_v,
                         const std::vector<amrex::Array<int,AMREX_SPACEDIM>>& is_per_v,
                         std::vector<std::string> prefix_v,
                         int max_step);
    // Destructor
    ~MultiBlockContainer ();

    // Initialize block data
    void InitializeBlocks ();

    // Set up BoxList vector for use with Communication Meta Data
    void SetBoxLists ();

    // Set up MB Communication Meta Data
    void SetBlockCommMetaData ();

    // Advance blocks
    void AdvanceBlocks ();

    // Fill MFs from other block(s)
    void FillPatchBlocks (int src_ind, int dst_ind);

private:
    int m_max_step;
    ERF erf1;
    ERF erf2;
    amrex::NonLocalBC::MultiBlockIndexMapping dtos;
    std::vector<std::vector<amrex::NonLocalBC::MultiBlockCommMetaData*>> cmd;
    std::vector<amrex::BoxList> blv;
};

#endif
