#include "MultiBlockContainer.H"
#include <AMReX_NonLocalBC.H>
#include <ERF.H>

// Vector input constructor
MultiBlockContainer::MultiBlockContainer(const std::vector<amrex::RealBox>& rb_v,
                                         std::vector<int> max_level_in_v,
                                         const std::vector<amrex::Vector<int>>& n_cell_in_v,
                                         std::vector<int> coord_v,
                                         const std::vector<amrex::Vector<amrex::IntVect>>& ref_ratios_v,
                                         const std::vector<amrex::Array<int,AMREX_SPACEDIM>>& is_per_v,
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
    dtos.permutation = amrex::IntVect{AMREX_D_DECL(   0,   1,   2)};
    dtos.sign        = amrex::IntVect{AMREX_D_DECL(   1,   1,   1)};

    // Set offset of dtos (NOTE: i_dst =  i_src - i_off -> [0] - [1])
    amrex::Real dx = ( rb_v[0].hi(0) - rb_v[0].lo(0) ) / n_cell_in_v[0][0];
    amrex::Real dy = ( rb_v[0].hi(1) - rb_v[0].lo(1) ) / n_cell_in_v[0][1];
    amrex::Real dz = ( rb_v[0].hi(2) - rb_v[0].lo(2) ) / n_cell_in_v[0][2];
    int offx = amrex::Math::floor(( rb_v[0].lo(0) - rb_v[1].lo(0) ) / dx);
    int offy = amrex::Math::floor(( rb_v[0].lo(1) - rb_v[1].lo(1) ) / dy);
    int offz = amrex::Math::floor(( rb_v[0].lo(2) - rb_v[1].lo(2) ) / dz);
    // DEBUG
    //offx=0; offy=0; offz=0;
    dtos.offset = amrex::IntVect{AMREX_D_DECL(offx, offy, offz)};
}

// Destructor
MultiBlockContainer::~MultiBlockContainer()
{
}

// Initialize block data
void
MultiBlockContainer::InitializeBlocks()
{
    erf1.InitData();
    erf2.InitData();
}

// Set up BoxList vector for use with Communication Meta Data
void
MultiBlockContainer::SetBoxLists()
{
    // Hard-coded bounds for now
    int nvars  = erf2.vars_new[0].size();
    int ndirs  = AMREX_SPACEDIM;

    for (int i(0); i<nvars; ++i) {
        // Get ghost cells, domain & grown box
        amrex::IntVect nghost = erf2.vars_new[0][i].nGrowVect();
        amrex::Box dom = erf2.domain_p[i];
        amrex::Box gbx = grow(dom,nghost);
        // Tmp BoxList
        amrex::BoxList bl;
        bl.clear();
        bl.set(gbx.ixType());
        for (int j(0); j<ndirs; ++j) {
            // Local box copies
            amrex::Box lgbx(gbx);
            amrex::Box ugbx(gbx);
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
        amrex::Print() << "DOM: " << erf2.domain_p[i] << "\n";
        amrex::Print() << "BA: " << erf2.vars_new[0][i].boxArray() << "\n";
        for (int j(0); j<6; ++j)
            amrex::Print() << (blv[i].data())[j] << "\n";

        amrex::Print() << "\n";
    }
    exit(0);
    */
}

// Set up MB Communication Meta Data
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
        amrex::IntVect nghost = erf2.vars_new[0][i].nGrowVect();
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

// Advance blocks
void
MultiBlockContainer::AdvanceBlocks()
{
    amrex::Print() << "STARTING MAIN DRIVER FOR: " << m_max_step << " STEPS" << "\n";
    amrex::Print() << "\n";

    for (int step(1); step <= m_max_step; ++step) {
        amrex::Print() << "    STARTING ADVANCE DRIVER: " << step << "\n";
        amrex::Print() << "===================================="  << "\n";
        erf1.Evolve_MB(step,1);
        amrex::Print() << '\n';
        amrex::Print() << "        SECOND BLOCK STARTS         "  << "\n";
        amrex::Print() << "------------------------------------"  << "\n";
        erf2.Evolve_MB(step,1);
        amrex::Print() << "COMPLETE" << "\n";
        amrex::Print() << "\n";
    }
}

// Wrapper for ParallelCopy between classes
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


