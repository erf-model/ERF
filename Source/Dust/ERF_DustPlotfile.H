/**
 * @file ERF_DustPlotfile.H
 * @brief Dust plotfile writer for output on native dust grid.
 *
 * Writes dust state variables to AMReX plotfile format on the native
 * 2D dust grid (independent of the ERF AMR level hierarchy).
 *
 * Pattern: ERF_FirePlotfile.H/cpp in Source/Fire/:
 *   https://github.com/hgopalan/ERF/blob/ERF-Fire/Source/Fire/ERF_FirePlotfile.cpp
 *
 * The fire module uses VisMF::Write + WriteGenericPlotfileHeader rather than
 * WriteMultiLevelPlotfile because the fire/dust grid is not part of the
 * AMR hierarchy. The same approach is used here.
 */
#ifdef ERF_USE_DUST
class DustLayer;

void WriteDustPlotfile(const std::string& plotfile_prefix,
                       const DustLayer&   dust_layer,
                       amrex::Real        time,
                       int                step);
#endif
