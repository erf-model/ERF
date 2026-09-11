#include "ERF_SurfaceModel.H"

#include "AMReX_Box.H"
#include "AMReX_MFIter.H"
#include <AMReX_PlotFileUtil.H>
#include "ERF_TileNoZ.H"

using namespace amrex;

void SurfaceModel::apply_weight_average(int lev, amrex::MultiFab *lsm_data, amrex::MultiFab *urban_data)
{
    bool valid_land = (lsm_data != nullptr);
    bool valid_urban = (urban_data != nullptr);

    AMREX_ASSERT_WITH_MESSAGE(valid_land || valid_urban, "Need at least one pointer to apply weights");

    // make sure the weights have been calculated before applying them
    AMREX_ALWAYS_ASSERT(m_weights_updated);

    // Apply weights separately to support different underlying grids between lsm_data and urban_data
    if (valid_land) {
        for (MFIter mfi(*lsm_data, TileNoZ()); mfi.isValid(); ++mfi)
        {
            Box tbx = mfi.tilebox();

            auto lsm_arr = lsm_data->array(mfi);
            auto weights_arr = wavg[lev]->const_array(mfi);

            // TODO: can amrex::Mult be used here instead?
            ParallelFor(tbx, [=] AMREX_GPU_DEVICE (int i, int j, int k)
            {
                lsm_arr(i, j, k) *= weights_arr(i, j, 0, SurfaceModelType::LAND);
            });

            // TODO: do we always want to update the boundaries?
            lsm_data->FillBoundary(m_geom[lev].periodicity());
        }
    }

    if (valid_urban) {
        for (MFIter mfi(*urban_data, TileNoZ()); mfi.isValid(); ++mfi)
        {
            Box tbx = mfi.tilebox();

            auto urb_arr = urban_data->array(mfi);
            auto weights_arr = wavg[lev]->const_array(mfi);

            // TODO: can amrex::Mult be used here instead?
            ParallelFor(tbx, [=] AMREX_GPU_DEVICE (int i, int j, int k)
            {
                urb_arr(i, j, k) *= weights_arr(i, j, 0, SurfaceModelType::URBAN);
            });

            // TODO: do we always want to update the boundaries?
            urban_data->FillBoundary(m_geom[lev].periodicity());
        }
    }
}


void SurfaceModel::calculate_weight_average(int lev, amrex::MultiFab* const urban_frac)
{
    AMREX_ASSERT_WITH_MESSAGE(m_use_land || m_use_urban, "Must have at least one land or urban model enabled for weighted average");

    // Calculate weighted averages between Land and Urban
    calculate_simple_average(lev, urban_frac);

    const bool use_urban = m_use_urban;
    const bool use_land = m_use_land;
    const bool use_fluxes = m_export_fluxes;

    // Reset surface fluxes
    u_star[lev]->setVal(0.0);
    t_star[lev]->setVal(0.0);
    q_star[lev]->setVal(0.0);

    amrex::MultiFab* const outputs[] = {u_star[lev].get(), t_star[lev].get(), q_star[lev].get(), t_surf[lev].get()};

    const int nfields = (m_export_fluxes) ? 5 : 4;

    // Output weighted surface fluxes into ustar, tstar, qstar and surface temperature into tsurf
    // TODO: make sure grids of urban and LSM inputs match
    for (int field=0; field < nfields; field++)
    {
        int output_field = (use_fluxes && field > 0) ? field - 1 : field;
        int comp = (use_fluxes && field < 2) ? field : 0;

        // whether we have a valid LSM multifab
        bool valid_land = (use_land &&
                           lsm_fields[field] != -1 &&
                           lsm_data_lev[lev][lsm_fields[field]]);

        // whether we have a valid urban multifab
        bool valid_urban = (use_urban &&
                            urban_fields[field] != -1 &&
                            urban_data_lev[lev][urban_fields[field]]);

        for (MFIter mfi(*u_star[lev], TileNoZ()); mfi.isValid(); ++mfi)
        {
            Box tbx = mfi.tilebox();

            auto weights_arr = wavg[lev]->const_array(mfi);

            // Outputs for surface boundary condition
            auto output_arr = outputs[output_field]->array(mfi);
            auto lsm_data_arr = (valid_land) ? lsm_data_lev[lev][lsm_fields[field]]->const_array(mfi) : Array4<const Real>{};
            auto urban_data_arr = (valid_urban) ? urban_data_lev[lev][urban_fields[field]]->const_array(mfi) : Array4<const Real>{};

            ParallelFor(tbx, [=] AMREX_GPU_DEVICE (int i, int j, int k)
            {
                Real land = (lsm_data_arr) ? lsm_data_arr(i, j, -1) * weights_arr(i, j, 0, SurfaceModelType::LAND) : 0.0;
                Real urb  = (urban_data_arr) ? urban_data_arr(i, j, k) * weights_arr(i, j, 0, SurfaceModelType::URBAN) : 0.0;

                output_arr(i, j, k, comp) = land + urb;

            });
        }

        outputs[output_field]->FillBoundary(comp, 1, m_geom[lev].periodicity());
    }

    // Apply the weight fractions to the rest of the LSM and urban fields (if any)
    if (use_land) {
        for (int field = nfields; field < lsm_fields.size(); field++) {
            if (lsm_fields[field] == -1) continue;
            for (MFIter mfi(*lsm_data_lev[lev][lsm_fields[field]], TileNoZ()); mfi.isValid(); ++mfi)
            {
                Box tbx = mfi.tilebox();

                auto weights_arr = wavg[lev]->const_array(mfi);

                // whether we have a valid LSM multifab
                bool valid_land = (use_land &&
                                   lsm_fields[field] != -1 &&
                                   lsm_data_lev[lev][lsm_fields[field]]);
                auto lsm_data_arr = (valid_land) ? lsm_data_lev[lev][lsm_fields[field]]->array(mfi) : Array4<Real>{};
            
                ParallelFor(tbx, [=] AMREX_GPU_DEVICE (int i, int j, int k)
                {
                    lsm_data_arr(i, j, k) *= weights_arr(i, j, 0, SurfaceModelType::LAND);
                });
            }
        }
    }

    if (use_urban) {
        for (int field = nfields; field < urban_fields.size(); field++) {
            if (urban_fields[field] == -1) continue;
            for (MFIter mfi(*urban_data_lev[lev][urban_fields[field]], TileNoZ()); mfi.isValid(); ++mfi)
            {
                Box tbx = mfi.tilebox();

                auto weights_arr = wavg[lev]->const_array(mfi);

                // whether we have a valid urban multifab
                bool valid_urban = (use_urban &&
                                    urban_fields[field] != -1 &&
                                    urban_data_lev[lev][urban_fields[field]]);
                auto urban_data_arr = (valid_urban) ? urban_data_lev[lev][urban_fields[field]]->array(mfi) : Array4<Real>{};

                ParallelFor(tbx, [=] AMREX_GPU_DEVICE (int i, int j, int k)
                {
                    urban_data_arr(i, j, k) *= weights_arr(i, j, 0, SurfaceModelType::URBAN);
                });
            }

            if (urban_data_lev[lev][urban_fields[field]]) {
                urban_data_lev[lev][urban_fields[field]]->FillBoundary(m_geom[lev].periodicity());
            }
        }
    }

    weight_average_fields(lev, urban_frac);
}

void SurfaceModel::calculate_simple_average(int lev, amrex::MultiFab* const urban_frac)
{
    // If no urban fraction is available, use the only enabled model as the
    // complete surface.  A missing fraction is ambiguous when both models
    // are enabled, so fail explicitly rather than silently choosing a model.
    if (urban_frac == nullptr) {
        AMREX_ALWAYS_ASSERT_WITH_MESSAGE(
            !(m_use_land && m_use_urban),
            "Urban fraction is required when both land and urban models are enabled");
    }

    const bool use_land = m_use_land;
    const bool use_urban = m_use_urban;

    for (MFIter mfi(*wavg[lev], TileNoZ()); mfi.isValid(); ++mfi)
    {
        Box tbx = mfi.tilebox();
        // weights are defined only on k = 0
        tbx = tbx.makeSlab(2, 0);

        auto urban_frac_arr = (urban_frac != nullptr) ? urban_frac->const_array(mfi) : Array4<const Real>{};
        auto weights_arr = wavg[lev]->array(mfi);

        ParallelFor(tbx, [=] AMREX_GPU_DEVICE (int i, int j, int k)
        {
            if (urban_frac_arr) {
                // Weights are proportional to urban fraction coverage in the current cell.
                weights_arr(i, j, 0, SurfaceModelType::URBAN) = urban_frac_arr(i, j, 0);
                weights_arr(i, j, 0, SurfaceModelType::LAND) = 1.0 - urban_frac_arr(i, j, 0);
            } else {
                // With one model disabled, the enabled model covers the whole cell.
                weights_arr(i, j, 0, SurfaceModelType::URBAN) = use_urban ? 1.0 : 0.0;
                weights_arr(i, j, 0, SurfaceModelType::LAND) = use_land ? 1.0 : 0.0;
            }
        });
    }

    m_weights_updated = true;
}

void SurfaceModel::register_field_map(std::string name, const std::pair<int, int> &lsm_urb_map, bool fill_boundary)
{
    amrex::Print() << " adding mapping between LSM<->Urban fields <"<<lsm_urb_map.first << "," << lsm_urb_map.second << "> to common surface name " << name << std::endl;

    AMREX_ALWAYS_ASSERT(!(lsm_urb_map.first == -1 && lsm_urb_map.second == -1));
    if (fieldmap.find(name) == fieldmap.end()) {

        Field field;
        field.map = lsm_urb_map;
        field.mf_ind = -1;
        field.fill_bound = fill_boundary;

        // Create MF to hold output for this field

        amrex::Vector<std::unique_ptr<amrex::MultiFab>> mf_lev(m_nlevs);
        for (int lev = 0; lev < m_nlevs; lev++)
        {
            mf_lev[lev] = std::make_unique<amrex::MultiFab>(m_ba2d[lev], m_dmap[lev], 1, IntVect(1,1,0));
            mf_lev[lev]->setVal(0.0);
        }

        fields.push_back(std::move(mf_lev));
        field.mf_ind = fields.size() - 1;

        amrex::Print() << "    -- created at ind = " << field.mf_ind << std::endl;

        fieldmap.insert({name, field});
    } else {
        return;
    }
}

void SurfaceModel::register_field_map(std::string name, amrex::Vector<amrex::MultiFab*> &lsm_lev_mf, amrex::Vector<amrex::MultiFab*> &urb_lev_mf, bool fill_boundary)
{
    // Register a field using explicit data pointers, rather than indices into lsm_data and urban_data
    AMREX_ALWAYS_ASSERT(lsm_lev_mf.size() > 0);
    amrex::Print() << " adding mapping between LSM<->Urban mfs <"<<lsm_lev_mf[0] << "," << urb_lev_mf[0] << "> to common surface name " << name << std::endl;

    AMREX_ALWAYS_ASSERT(lsm_lev_mf.size() == urb_lev_mf.size());
    if (fieldmap.find(name) == fieldmap.end()) {

        Field field;
        field.map = std::pair<int,int>(-1,-1);
        field.mf_ind = -1;
        field.fill_bound = fill_boundary;
        field.lsm_ptr = lsm_lev_mf;
        field.urb_ptr = urb_lev_mf;

        // Create MF to hold output for this field

        amrex::Vector<std::unique_ptr<amrex::MultiFab>> mf_lev(m_nlevs);
        for (int lev = 0; lev < m_nlevs; lev++)
        {
            mf_lev[lev] = std::make_unique<amrex::MultiFab>(m_ba2d[lev], m_dmap[lev], 1, IntVect(1,1,0));
            mf_lev[lev]->setVal(0.0);
        }

        fields.push_back(std::move(mf_lev));
        field.mf_ind = fields.size() - 1;

        amrex::Print() << "    -- created at ind = " << field.mf_ind << std::endl;

        fieldmap.insert({name, field});
    } else {
        return;
    }
}

void SurfaceModel::weight_average_fields(int lev, amrex::MultiFab* const urban_frac)
{
    for (auto &field : fieldmap)
    {
        int mf_idx = field.second.mf_ind;
        AMREX_ASSERT(mf_idx != -1);

        int lsm_idx = field.second.map.first;
        int urb_idx = field.second.map.second;

        // whether we have a valid LSM multifab
        bool valid_land = (m_use_land &&
                           lsm_idx != -1 &&
                           lsm_data_lev[lev][lsm_idx]);

        // whether we have a valid urban multifab
        bool valid_urban = (m_use_urban &&
                            urb_idx != -1 &&
                            urban_data_lev[lev][urb_idx]);

        bool use_mf = false;
        if (lsm_idx == -1 && urb_idx == -1) {
            // use explicit MF ptrs rather than indices
            use_mf = true;
            valid_land = (m_use_land && field.second.lsm_ptr[lev]);
            valid_urban = (m_use_urban && field.second.urb_ptr[lev]);
        }

        for (MFIter mfi(*fields[mf_idx][lev], TileNoZ()); mfi.isValid(); ++mfi)
        {
            Box tbx = mfi.tilebox();
            Box b2d = makeSlab(tbx, 2, 0);

            auto weights_arr = wavg[lev]->const_array(mfi);

            // Calculate weight average into output
            auto output_arr = fields[mf_idx][lev]->array(mfi);
            const amrex::MultiFab *lsm_mf = (use_mf) ? field.second.lsm_ptr[lev] : (lsm_idx != -1 ? lsm_data_lev[lev][lsm_idx] : nullptr);
            const amrex::MultiFab *urb_mf = (use_mf) ? field.second.urb_ptr[lev] : (urb_idx != -1 ? urban_data_lev[lev][urb_idx] : nullptr);

            auto lsm_data_arr = (valid_land) ? lsm_mf->const_array(mfi) : Array4<const Real>{};
            auto urban_data_arr = (valid_urban) ? urb_mf->const_array(mfi) : Array4<const Real>{};

            if (valid_land && !valid_urban) {
                // use solely land value
                ParallelFor(tbx, [=] AMREX_GPU_DEVICE (int i, int j, int k)
                {
                    output_arr(i, j, k) = lsm_data_arr(i, j, k);
                });
            } else if (!valid_land && valid_urban) {
                // use solely urban value
                ParallelFor(tbx, [=] AMREX_GPU_DEVICE (int i, int j, int k)
                {
                    output_arr(i, j, k) = urban_data_arr(i, j, k);
                });
            } else {
                ParallelFor(tbx, [=] AMREX_GPU_DEVICE (int i, int j, int k)
                {
                    Real land = (lsm_data_arr) ? lsm_data_arr(i, j, k) * weights_arr(i, j, 0, SurfaceModelType::LAND) : 0.0;
                    Real urb  = (urban_data_arr) ? urban_data_arr(i, j, k) * weights_arr(i, j, 0, SurfaceModelType::URBAN) : 0.0;

                    output_arr(i, j, k) = land + urb;
                });
            }
        }
        if (field.second.fill_bound) {
            fields[mf_idx][lev]->FillBoundary(m_geom[lev].periodicity());
        }
        //outputs[output_field]->FillBoundary(comp, 1, m_geom[lev].periodicity());
    }
}


void SurfaceModel::write_output(const int &finest_lev, const amrex::Real &time, const std::string &plot_prefix, const amrex::Vector<int> &level_steps, const amrex::Vector<amrex::IntVect> &ref_ratio)
{
    std::string plotfilename = amrex::Concatenate(plot_prefix + "2D_", level_steps[0], 5);

    const int nfields = (m_export_fluxes) ? 5 : 4;
    const int nlsm_fields = (m_use_land) ? (lsm_fields.size() - nfields) : 0;
    const int nurb_fields = (m_use_urban) ? (urban_fields.size() - nfields) : 0;
    const int noutput = nfields + 2 + fieldmap.size(); // + nlsm_fields + nurb_fields; // MOST, lmask, urb frac, lsm fields, urban fields
    IntVect ng(0, 0, 0);

    amrex::Vector<std::string> varnames(nfields);

    amrex::Vector<amrex::MultiFab> fab(finest_lev+1);
    for (int lev = 0; lev <= finest_lev; lev++) {
        amrex::MultiFab* const outputs[] = {u_star[lev].get(), t_star[lev].get(), q_star[lev].get(), t_surf[lev].get()};
        fab[lev].define(outputs[0]->boxArray(), m_dmap[lev], noutput, ng);

        //fab[lev].setVal(0.0);

        // Output weighted surface fluxes into ustar, tstar, qstar and surface temperature into tsurf
        // TODO: make sure grids of urban and LSM inputs match
        for (int field=0; field < nfields; field++)
        {
            int output_field = (m_export_fluxes && field > 0) ? field - 1 : field;
            int comp = (m_export_fluxes && field < 2) ? field : 0;

            MultiFab::Copy(fab[lev], *(outputs[output_field]), comp, field, 1, ng);

            varnames[field] = field_names[field];
        }

        int nout = nfields;

        MultiFab lmask_tmp = amrex::ToMultiFab(*m_lmask[lev]); // iMultiFab -> MultiFab
        MultiFab::Copy(fab[lev], lmask_tmp, 0, nout, 1, ng);
        nout++;
        MultiFab::Copy(fab[lev], *(wavg[lev]), SurfaceModelType::URBAN, nout, 1, ng);
        nout++;

        if (lev == 0) {
            varnames.push_back("lmask");
            varnames.push_back("urb_frac");
        }

        // Add any mapped fields
        for (auto &field : fieldmap) {
            MultiFab::Copy(fab[lev], *(fields[field.second.mf_ind][lev]), 0, nout, 1, ng);
            if (lev == 0) {
                varnames.push_back(field.first);
            }
            nout++;
        }

        AMREX_ALWAYS_ASSERT(varnames.size() == noutput);
    }

    amrex::WriteMultiLevelPlotfile(plotfilename, finest_lev+1, GetVecOfConstPtrs(fab), varnames, m_geom2d, time, level_steps, ref_ratio);
}

// utility to skip to next line in Header
//  -- taken from ERF_Checkpoint
void
SurfaceModel::GotoNextLine (std::istream& is)
{
    constexpr std::streamsize bl_ignore_max { 100000 };
    is.ignore(bl_ignore_max, '\n');
}

void SurfaceModel::WriteCheckpoint(const std::string &checkpointname)
{
    auto check_start = amrex::second();

    // write header
    if (ParallelDescriptor::IOProcessor()) {

        amrex::Print() << " Writing SurfaceModel checkpoint " << std::endl;

        std::string HeaderFileName(checkpointname + "/SurfaceModel_Header");
        VisMF::IO_Buffer io_buffer(VisMF::IO_Buffer_Size);
        std::ofstream HeaderFile;
        HeaderFile.rdbuf()->pubsetbuf(io_buffer.dataPtr(), io_buffer.size());
        HeaderFile.open(HeaderFileName.c_str(), std::ofstream::out   |
                                                std::ofstream::trunc |
                                                std::ofstream::binary);
        if(! HeaderFile.good()) {
            FileOpenFailed(HeaderFileName);
        }

        HeaderFile.precision(17);

        // write out title line
        HeaderFile << "Checkpoint file for SurfaceModel\n";

        // write out number of levels
        HeaderFile << m_nlevs << "\n";

        // write out flags
        HeaderFile << m_use_urban << "\n";
        HeaderFile << m_use_land << "\n";
        HeaderFile << m_export_fluxes << "\n";
        HeaderFile << m_weights_updated << "\n";
        HeaderFile << "\n";

        // Write box arrays
        for (int lev = 0; lev < m_nlevs; lev++) {
            m_ba[lev].writeOn(HeaderFile);
            HeaderFile << '\n';

        }
        HeaderFile << '\n';
        for (int lev = 0; lev < m_nlevs; lev++) {
            m_ba2d[lev].writeOn(HeaderFile);
            HeaderFile << '\n';
        }
        HeaderFile << '\n';

        // Write out fields
        // LSM fields
        for (int i = 0; i < lsm_fields.size(); ++i) {
            HeaderFile << lsm_fields[i] << " ";
        }
        HeaderFile << '\n';

        // Urban fields
        for (int i = 0; i < urban_fields.size(); ++i) {
            HeaderFile << urban_fields[i] << " ";
        }
        HeaderFile << '\n';

        // Field mapping
        HeaderFile << '\n';
        HeaderFile << fieldmap.size() << "\n";
        for (auto &field : fieldmap) {
            // pointer fields are reconstructed at load
            HeaderFile << field.first << " " << field.second.mf_ind << " " << field.second.map.first << " " << field.second.map.second << " " << field.second.fill_bound << '\n';
        }
        HeaderFile << '\n';
    }

    amrex::ParallelDescriptor::Barrier();

    const std::string prefix = "SurfaceModel_";

    for (int lev = 0; lev < m_nlevs; lev++) {
        IntVect ng(1,1,0);

        {
            MultiFab mf(m_ba2d[lev],m_dmap[lev],2,ng);
            MultiFab::Copy(mf,*u_star[lev],0,0,2,ng);
            VisMF::Write(mf, MultiFabFileFullPrefix(lev, checkpointname, "Level_", prefix + "ustar"));

            MultiFab::Copy(mf,*wavg[lev],0,0,2,ng);
            VisMF::Write(mf, MultiFabFileFullPrefix(lev, checkpointname, "Level_", prefix + "wavg"));
        }

        MultiFab mf(m_ba2d[lev],m_dmap[lev],1,ng);

        MultiFab::Copy(mf,*t_star[lev],0,0,1,ng);
        VisMF::Write(mf, MultiFabFileFullPrefix(lev, checkpointname, "Level_", prefix + "tstar"));

        MultiFab::Copy(mf,*q_star[lev],0,0,1,ng);
        VisMF::Write(mf, MultiFabFileFullPrefix(lev, checkpointname, "Level_", prefix + "qstar"));

        MultiFab::Copy(mf,*t_surf[lev],0,0,1,ng);
        VisMF::Write(mf, MultiFabFileFullPrefix(lev, checkpointname, "Level_", prefix + "tsurf"));

        for (auto &field : fieldmap) {
            MultiFab::Copy(mf, *(fields[field.second.mf_ind][lev]), 0, 0, 1, ng);
            VisMF::Write(mf, MultiFabFileFullPrefix(lev, checkpointname, "Level_", prefix + "_f_" + field.first));
        }
    }

    auto check_end = amrex::second() - check_start;
    ParallelDescriptor::ReduceRealMax(check_end,ParallelDescriptor::IOProcessorNumber());
    amrex::Print() << "    SurfaceModel Checkpoint write time = " << check_end << " seconds." << '\n';
}

void SurfaceModel::ReadCheckpoint(const std::string &checkpointname)
{
    auto check_start = amrex::second();

    amrex::Print() << " Reading SurfaceModel checkpoint " << std::endl;

    // Header
    std::string File(checkpointname + "/SurfaceModel_Header");

    VisMF::IO_Buffer io_buffer(VisMF::GetIOBufferSize());

    Vector<char> fileCharPtr;
    ParallelDescriptor::ReadAndBcastFile(File, fileCharPtr);
    std::string fileCharPtrString(fileCharPtr.dataPtr());
    std::istringstream is(fileCharPtrString, std::istringstream::in);

    std::string line, word;

    int chk_ncomp_cons, chk_ncomp;

    // read in title line
    std::getline(is, line);

    // read in number of levels
    is >> m_nlevs;

    // read in flags
    is >> m_use_urban;
    is >> m_use_land;
    is >> m_export_fluxes;
    is >> m_weights_updated;
    GotoNextLine(is);

    // read in box arrays for each level
    m_ba.resize(m_nlevs);
    m_ba2d.resize(m_nlevs);
    for (int lev = 0; lev < m_nlevs; ++lev) {
        BoxArray ba;
        ba.readFrom(is);
        AMREX_ALWAYS_ASSERT(ba == m_ba[lev]);
        GotoNextLine(is);
    }
    GotoNextLine(is);
    for (int lev = 0; lev < m_nlevs; ++lev) {
        BoxArray ba;
        ba.readFrom(is);
        AMREX_ALWAYS_ASSERT(ba == m_ba2d[lev]);
        GotoNextLine(is);
    }
    GotoNextLine(is);

    // Read in LSM fields
    std::getline(is, line);
    {
        std::istringstream lis(line);
        lsm_fields.clear();
        while (lis >> word) {
            lsm_fields.push_back(std::stoi(word));
        }
    }

    // Read in Urban fields
    std::getline(is, line);
    {
        std::istringstream lis(line);
        urban_fields.clear();
        while (lis >> word) {
            urban_fields.push_back(std::stoi(word));
        }
    }
    GotoNextLine(is);

    // Read number of mapped fields
    int nfields = 0;
    is >> nfields;
    GotoNextLine(is);

    // Read any mapped fields
    AMREX_ALWAYS_ASSERT(nfields == fieldmap.size());

    // The field mapping should already be created before the restart, but verify consistency with the file
    for (int i = 0; i < nfields; i++) {
        std::getline(is, line);
        std::istringstream lis(line);

        std::string field_name;
        int lsm_ind;
        int urb_ind;
        bool fill_bound;

        lis >> field_name;
        lis >> word; // mf_ind (ignored for restart)
        lis >> word;
        lsm_ind = std::stoi(word);
        lis >> word;
        urb_ind = std::stoi(word);
        lis >> word;
        fill_bound = static_cast<bool>(std::stoi(word));

        AMREX_ALWAYS_ASSERT(fieldmap.find(field_name) != fieldmap.end());
        const auto &field = fieldmap.at(field_name);
        AMREX_ALWAYS_ASSERT(field.map.first == lsm_ind);
        AMREX_ALWAYS_ASSERT(field.map.second == urb_ind);
        AMREX_ALWAYS_ASSERT(field.fill_bound == fill_bound);
        // skip checking mf_ind because the ordering could have changed
    }

    const std::string prefix = "SurfaceModel_";

    IntVect ng = IntVect(1,1,0);
    for (int lev = 0; lev < m_nlevs; lev++) {
        {
            MultiFab mf(m_ba2d[lev],m_dmap[lev],2,ng);
            VisMF::Read(mf, MultiFabFileFullPrefix(lev, checkpointname, "Level_", prefix + "ustar"));
            MultiFab::Copy(*(u_star[lev]),mf,0,0,2,ng);

            VisMF::Read(mf, MultiFabFileFullPrefix(lev, checkpointname, "Level_", prefix + "wavg"));
            MultiFab::Copy(*wavg[lev],mf,0,0,2,ng);
        }

        MultiFab mf(m_ba2d[lev],m_dmap[lev],1,ng);

        VisMF::Read(mf, MultiFabFileFullPrefix(lev, checkpointname, "Level_", prefix + "tstar"));
        MultiFab::Copy(*t_star[lev],mf,0,0,1,ng);

        VisMF::Read(mf, MultiFabFileFullPrefix(lev, checkpointname, "Level_", prefix + "qstar"));
        MultiFab::Copy(*q_star[lev],mf,0,0,1,ng);

        VisMF::Read(mf, MultiFabFileFullPrefix(lev, checkpointname, "Level_", prefix + "tsurf"));
        MultiFab::Copy(*t_surf[lev],mf,0,0,1,ng);

        for (auto &field : fieldmap) {
            VisMF::Read(mf, MultiFabFileFullPrefix(lev, checkpointname, "Level_", prefix + "_f_" + field.first));
            MultiFab::Copy(*(fields[field.second.mf_ind][lev]), mf, 0, 0, 1, ng);
        }
    }

    auto check_end = amrex::second() - check_start;
    ParallelDescriptor::ReduceRealMax(check_end,ParallelDescriptor::IOProcessorNumber());
    amrex::Print() << "    SurfaceModel Checkpoint load time = " << check_end << " seconds." << '\n';
}
