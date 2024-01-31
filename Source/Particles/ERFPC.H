#ifndef ERF_PC_H_
#define ERF_PC_H_

#ifdef ERF_USE_PARTICLES

#include <string>
#include <AMReX_Particles.H>

struct ERFParticlesRealIdx
{
    enum {
        vx = 0,
        vy,
        vz,
        mass,
        ncomps
    };
};

struct ERFParticlesIntIdx
{
    enum {
        i = 0,
        j,
        k,
        ncomps
    };
};

namespace ERFParticleInitializations
{
    /* list of particle initializations */
    const std::string init_default = "default";
    const std::string init_uniform = "uniform";

}

namespace ERFParticleNames
{
    const std::string tracers = "tracer_particles";
    const std::string hydro = "hydro_particles";
}

struct ERFParticlesAssignor
{
    template <typename P>
    AMREX_GPU_HOST_DEVICE
    amrex::IntVect operator() (  P const& p,
                                 amrex::GpuArray<amrex::Real,AMREX_SPACEDIM> const& plo,
                                 amrex::GpuArray<amrex::Real,AMREX_SPACEDIM> const& dxi,
                                 const amrex::Box& domain ) const noexcept
    {
        amrex::IntVect iv(
            AMREX_D_DECL(   int(amrex::Math::floor((p.pos(0)-plo[0])*dxi[0])),
                            int(amrex::Math::floor((p.pos(1)-plo[1])*dxi[1])),
                            p.idata(ERFParticlesIntIdx::k) ) );
        iv[0] += domain.smallEnd()[0];
        iv[1] += domain.smallEnd()[1];
        return iv;
    }
};

template <typename P>
AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE
void update_location_idata (P& p,
                            amrex::GpuArray<amrex::Real,AMREX_SPACEDIM> const& plo,
                            amrex::GpuArray<amrex::Real,AMREX_SPACEDIM> const& dxi,
                            const amrex::Array4<amrex::Real const>&  height_arr)
{
    amrex::IntVect iv( int(amrex::Math::floor((p.pos(0)-plo[0])*dxi[0])),
                       int(amrex::Math::floor((p.pos(1)-plo[1])*dxi[1])),
                       p.idata(ERFParticlesIntIdx::k) );

    p.idata(ERFParticlesIntIdx::i) = iv[0];
    p.idata(ERFParticlesIntIdx::j) = iv[1];

    if (height_arr) {
        amrex::Real lx = (p.pos(0)-plo[0])*dxi[0] - static_cast<amrex::Real>(iv[0]);
        amrex::Real ly = (p.pos(1)-plo[1])*dxi[1] - static_cast<amrex::Real>(iv[1]);
        auto zlo = height_arr(iv[0]  ,iv[1]  ,iv[2]  ) * (1.0-lx) * (1.0-ly) +
                   height_arr(iv[0]+1,iv[1]  ,iv[2]  ) *      lx  * (1.0-ly) +
                   height_arr(iv[0]  ,iv[1]+1,iv[2]  ) * (1.0-lx) * ly +
                   height_arr(iv[0]+1,iv[1]+1,iv[2]  ) *      lx  * ly;
        auto zhi = height_arr(iv[0]  ,iv[1]  ,iv[2]+1) * (1.0-lx) * (1.0-ly) +
                   height_arr(iv[0]+1,iv[1]  ,iv[2]+1) *      lx  * (1.0-ly) +
                   height_arr(iv[0]  ,iv[1]+1,iv[2]+1) * (1.0-lx) * ly +
                   height_arr(iv[0]+1,iv[1]+1,iv[2]+1) *      lx  * ly;

        if (p.pos(2) > zhi) {
            p.idata(ERFParticlesIntIdx::k) += 1;
        } else if (p.pos(2) <= zlo) {
            p.idata(ERFParticlesIntIdx::k) -= 1;
        }
    }
}

class ERFPC : public amrex::ParticleContainer<  ERFParticlesRealIdx::ncomps,
                                                ERFParticlesIntIdx::ncomps,
                                                0,
                                                0,
                                                amrex::DefaultAllocator,
                                                ERFParticlesAssignor >
{
    public:

        /*! Constructor */
        ERFPC ( amrex::ParGDBBase* a_gdb,
                const std::string& a_name = "particles" )
            : amrex::ParticleContainer< ERFParticlesRealIdx::ncomps,
                                        ERFParticlesIntIdx::ncomps,
                                        0,
                                        0,
                                        amrex::DefaultAllocator,
                                        ERFParticlesAssignor> (a_gdb)
        {
            BL_PROFILE("ERFPCPC::ERFPC()");
            m_name = a_name;
            readInputs();
        }

        /*! Constructor */
        ERFPC ( const amrex::Geometry&            a_geom,
                const amrex::DistributionMapping& a_dmap,
                const amrex::BoxArray&            a_ba,
                const std::string&                a_name = "particles" )
            : amrex::ParticleContainer< ERFParticlesRealIdx::ncomps,
                                        ERFParticlesIntIdx::ncomps,
                                        0,
                                        0,
                                        amrex::DefaultAllocator,
                                        ERFParticlesAssignor> ( a_geom, a_dmap, a_ba )
        {
            BL_PROFILE("ERFPCPC::ERFPC()");
            m_name = a_name;
            readInputs();
        }

        /*! Initialize particles in domain */
        virtual void InitializeParticles (const std::unique_ptr<amrex::MultiFab>& a_ptr = nullptr);

        /*! Evolve particles for one time step */
        virtual void EvolveParticles (   int,
                                         amrex::Real,
                                         amrex::Vector<amrex::Vector<amrex::MultiFab>>&,
                                         const amrex::Vector<std::unique_ptr<amrex::MultiFab>>& );

        /*! Get real-type particle attribute names */
        virtual amrex::Vector<std::string> varNames () const
        {
            BL_PROFILE("ERFPCPC::varNames()");
            return {AMREX_D_DECL("xvel","yvel","zvel"),"mass"};
        }

        /*! Uses midpoint method to advance particles using flow velocity. */
        virtual void AdvectWithFlow (  amrex::MultiFab*,
                                       int,
                                       amrex::Real,
                                       const std::unique_ptr<amrex::MultiFab>& );

        /*! Uses midpoint method to advance particles falling under gravity. */
        virtual void AdvectWithGravity ( int,
                                         amrex::Real,
                                         const std::unique_ptr<amrex::MultiFab>& );

        /*! Specify if particles should advect with flow */
        inline void setAdvectWithFlow (bool a_flag)
        {
            BL_PROFILE("ERFPCPC::setAdvectWithFlow()");
            m_advect_w_flow = a_flag;
        }
        /*! Specify if particles fall under gravity */
        inline void setAdvectWithGravity (bool a_flag)
        {
            BL_PROFILE("ERFPCPC::setAdvectWithGravity()");
            m_advect_w_gravity = a_flag;
        }

    protected:

        bool m_advect_w_flow;               /*!< advect with flow velocity */
        bool m_advect_w_gravity;            /*!< advect under gravitational force */

        std::string m_name;                 /*!< name of this particle species */

        std::string m_initialization_type;  /*!< initial particle distribution type */
        int m_ppc_init;                     /*!< initial number of particles per cell */

        /*! read inputs from file */
        virtual void readInputs ();

        /*! Default particle initialization */
        void initializeParticlesUniformDistribution (const std::unique_ptr<amrex::MultiFab>& a_ptr = nullptr);

    private:

        /*! Default particle initialization */
        void initializeParticlesDefault (const std::unique_ptr<amrex::MultiFab>& a_ptr = nullptr);
        /*! Default initialization for tracer particles for WoA case (ref: AA) */
        void initializeParticlesDefaultTracersWoA (const std::unique_ptr<amrex::MultiFab>& a_ptr=nullptr);
        /*! Default initialization for hydro particles (ref: AA) */
        void initializeParticlesDefaultHydro (const std::unique_ptr<amrex::MultiFab>& a_ptr = nullptr);

};

#endif
#endif
