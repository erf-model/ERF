#include "MultiBlockContainer.H"
#include <AMReX_NonLocalBC.H>
#include <ERF.H>

using namespace amrex;

/**
 * Constructor for the MultiBlockContainer class capable of taking a vector of boxes as input.
 *
 * Inputs that are vectors are used to define the two ERF instances in the MultiBlock class.
 *
 * @param[in] rb_v Vector of RealBoxes to define this MultiBlock
 * @param[in] max_level_in_v Maximum level vector
 * @param[in] n_cell_in_v Number of cells vector
 * @param[in] coord_v Coordinate selection vector
 * @param[in] ref_ratios_v Refinement ratio vector
 * @param[in] is_per_v Defines whether the domain is periodic in each coordinate direction
 * @param[in] prefix_v Prefixes for ParmParse as a vector
 * @param[in] max_step Maximum number of timesteps to take
 */
MultiBlockContainer::MultiBlockContainer(const std::vector<RealBox>& rb_v,
                                         std::vector<int> max_level_in_v,
                                         const std::vector<Vector<int>>& n_cell_in_v,
                                         std::vector<int> coord_v,
                                         const std::vector<Vector<IntVect>>& ref_ratios_v,
                                         const std::vector<Array<int,AMREX_SPACEDIM>>& is_per_v,
                                         std::vector<std::string> prefix_v,
                                         int max_step)
: m_max_step(max_step),
  erf1(rb_v[0],max_level_in_v[0],n_cell_in_v[0],coord_v[0],ref_ratios_v[0],is_per_v[0],prefix_v[0]),
  erf2(rb_v[1],max_level_in_v[1],n_cell_in_v[1],coord_v[1],ref_ratios_v[1],is_per_v[1],prefix_v[1])
{
    // Store ptr to container to call member functions
    erf1.SetMultiBlockPointer(this);
    erf2.SetMultiBlockPointer(this);

    // Set the permutation/sign of dtos
    dtos.permutation = IntVect{AMREX_D_DECL(   0,   1,   2)};
    dtos.sign        = IntVect{AMREX_D_DECL(   1,   1,   1)};

    // Set offset of dtos (NOTE: i_dst =  i_src - i_off -> [0] - [1])
    Real dx = ( rb_v[0].hi(0) - rb_v[0].lo(0) ) / n_cell_in_v[0][0];
    Real dy = ( rb_v[0].hi(1) - rb_v[0].lo(1) ) / n_cell_in_v[0][1];
    Real dz = ( rb_v[0].hi(2) - rb_v[0].lo(2) ) / n_cell_in_v[0][2];
    int offx = amrex::Math::floor(( rb_v[0].lo(0) - rb_v[1].lo(0) ) / dx);
    int offy = amrex::Math::floor(( rb_v[0].lo(1) - rb_v[1].lo(1) ) / dy);
    int offz = amrex::Math::floor(( rb_v[0].lo(2) - rb_v[1].lo(2) ) / dz);
    // DEBUG
    //offx=0; offy=0; offz=0;
    dtos.offset = IntVect{AMREX_D_DECL(offx, offy, offz)};
}

/**
 * Destructor for the MultiBlockContainer
 */
MultiBlockContainer::~MultiBlockContainer()
{
}

/**
 * Initialize block data for the MultiBlockContainer
 */
void
MultiBlockContainer::InitializeBlocks()
{
    erf1.InitData();
    erf2.InitData();
}

/**
 * Set up BoxList vector for use with Communication Meta Data
 */
void
MultiBlockContainer::SetBoxLists()
{
    // Hard-coded bounds for now
    int nvars  = erf2.vars_new[0].size();
    int ndirs  = AMREX_SPACEDIM;

    for (int i(0); i<nvars; ++i) {
        // Get ghost cells, domain & grown box
        IntVect nghost = erf2.vars_new[0][i].nGrowVect();
        Box dom = erf2.domain_p[i];
        Box gbx = grow(dom,nghost);
        // Tmp BoxList
        BoxList bl;
        bl.clear();
        bl.set(gbx.ixType());
        for (int j(0); j<ndirs; ++j) {
            // Local box copies
            Box lgbx(gbx);
            Box ugbx(gbx);
            // Get lower & upper bound
            int se = dom.smallEnd(j) - 1;
            int be = dom.bigEnd(j)   + 1;
            // Modify bounds for nodal vars
            if (gbx.ixType().nodeCentered(j)) {
                se += 1;
                be -= 1;
            }
            // Populate BoxList
            bl.push_back( lgbx.setBig(j,se) );
            bl.push_back( ugbx.setSmall(j,be) );
        }
        blv.push_back(bl);
    }

    /*
    // DEBUG BOX LIST
    for (int i(0); i<nvars; ++i) {
        Print() << "DOM: " << erf2.domain_p[i] << "\n";
        Print() << "BA: " << erf2.vars_new[0][i].boxArray() << "\n";
        for (int j(0); j<6; ++j)
            Print() << (blv[i].data())[j] << "\n";

        Print() << "\n";
    }
    exit(0);
    */
}

/**
 * Set up MultiBlock Communication Meta Data
 */
void
MultiBlockContainer::SetBlockCommMetaData()
{
    // Hard-coded bounds for now
    int nvars  = erf2.vars_new[0].size(); // Destination MF
    int ndirs  = AMREX_SPACEDIM;

    // Loop over num_vars to set communicator
    for (int i(0); i<nvars; ++i) {
        // Make space
        cmd.push_back(std::vector<amrex::NonLocalBC::MultiBlockCommMetaData*>());
        // Get ghost cell vector for multifab growth
        IntVect nghost = erf2.vars_new[0][i].nGrowVect();
        for (int j(0); j<2*ndirs; ++j) {

            // Store temp ptr to communicator for i^th variable
            amrex::NonLocalBC::MultiBlockCommMetaData *cmd_tmp =
                new amrex::NonLocalBC::MultiBlockCommMetaData(erf2.vars_new[0][i], (blv[i].data())[j],
                                                              erf1.vars_new[0][i], nghost, dtos);
            // Populate cmd vector
            cmd[i].push_back(cmd_tmp);
        }
    }
}

/**
 * Advance blocks in the MultiBlockContainer by calling each timestep advance sequentially.
 */
void
MultiBlockContainer::AdvanceBlocks()
{
    Print() << "STARTING MAIN DRIVER FOR: " << m_max_step << " STEPS" << "\n";
    Print() << "\n";

    for (int step(1); step <= m_max_step; ++step) {
        Print() << "    STARTING ADVANCE DRIVER: " << step << "\n";
        Print() << "===================================="  << "\n";
        erf1.Evolve_MB(step,1);
        Print() << '\n';
        Print() << "        SECOND BLOCK STARTS         "  << "\n";
        Print() << "------------------------------------"  << "\n";
        erf2.Evolve_MB(step,1);
        Print() << "COMPLETE" << "\n";
        Print() << "\n";
    }
}

/**
 * Wrapper for ParallelCopy between classes
 */
void
MultiBlockContainer::FillPatchBlocks(int src_ind, int dst_ind)
{
    // Hard-coded bounds for now
    int ndirs  = AMREX_SPACEDIM;

    // Loop faces of box to perform ParallelCopy
    // NOTE - cmd built with ERF2 so uses dst_ind
    for (int j(0); j<2*ndirs; ++j)
        amrex::NonLocalBC::ParallelCopy(erf2.vars_new[0][dst_ind], erf1.vars_new[0][src_ind],
                                        *(cmd[dst_ind][j]), 0, 0, 1, dtos);
}


