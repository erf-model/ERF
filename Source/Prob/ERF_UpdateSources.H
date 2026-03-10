    void update_rhotheta_sources (
        const amrex::Real& /*time*/,
        amrex::MultiFab* src,
        const amrex::Geometry& geom,
        std::unique_ptr<amrex::MultiFab>& z_phys_cc) override;

    void update_rhoqt_sources (
        const amrex::Real& /*time*/,
        amrex::MultiFab* qsrc,
        const amrex::Geometry& geom,
        std::unique_ptr<amrex::MultiFab>& z_phys_cc) override;

    void update_w_subsidence (
        const amrex::Real& /*time*/,
        amrex::Vector<amrex::Real>& wbar,
        amrex::Gpu::DeviceVector<amrex::Real>& d_wbar,
        const amrex::MultiFab& state,
        const amrex::Geometry& geom,
        std::unique_ptr<amrex::MultiFab>& z_phys_nd) override;

    void update_geostrophic_profile (
        const amrex::Real& /*time*/,
        amrex::Vector<amrex::Real>& u_geos,
        amrex::Gpu::DeviceVector<amrex::Real>& d_u_geos,
        amrex::Vector<amrex::Real>& v_geos,
        amrex::Gpu::DeviceVector<amrex::Real>& d_v_geos,
        const amrex::Geometry& geom,
        std::unique_ptr<amrex::MultiFab>& z_phys_cc) override;
