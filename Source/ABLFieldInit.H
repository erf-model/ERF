#ifndef ABLFIELDINIT_H
#define ABLFIELDINIT_H

#include "AMReX_Gpu.H"
#include "AMReX_MultiFab.H"
#include "ERF_Constants.H"
#include "IndexDefines.H"

/** Initialize ABL fields: The contents of this file were copied almost verbatim from amr-wind:
    https://github.com/Exawind/amr-wind
 */
struct ABLFieldInit {

  public:

    void init_params()
    {
        amrex::ParmParse pp_abl("ABL");

        // Temperature variation as a function of height
        pp_abl.getarr("temperature_heights", m_theta_heights);
        pp_abl.getarr("temperature_values", m_theta_values);

        AMREX_ALWAYS_ASSERT(m_theta_heights.size() == m_theta_values.size());
        int num_theta_values = m_theta_heights.size();

        pp_abl.query("perturb_velocity", m_perturb_vel);
        pp_abl.query("perturb_ref_height", m_ref_height);
        pp_abl.query("Uperiods", m_Uperiods);
        pp_abl.query("Vperiods", m_Vperiods);

        pp_abl.query("perturb_temperature", m_perturb_theta);
        pp_abl.query("random_gauss_mean", m_theta_gauss_mean);
        pp_abl.query("random_gauss_var", m_theta_gauss_var);
        pp_abl.query("cutoff_height", m_theta_cutoff_height);
        pp_abl.query("theta_amplitude", m_deltaT);

        pp_abl.query("init_tke", m_tke_init);
        pp_abl.query("init_qke", m_qke_init);

        m_thht_d.resize(num_theta_values);
        m_thvv_d.resize(num_theta_values);

        amrex::Gpu::copy(
            amrex::Gpu::hostToDevice, m_theta_heights.begin(),
            m_theta_heights.end(), m_thht_d.begin());
        amrex::Gpu::copy(
            amrex::Gpu::hostToDevice, m_theta_values.begin(), m_theta_values.end(),
            m_thvv_d.begin());

        amrex::ParmParse pp_prob("prob");

        pp_prob.query("rho_0", m_rho);
        pp_prob.query("U0"   , m_vel_x);
        pp_prob.query("V0"   , m_vel_y);
        pp_prob.query("W0"   , m_vel_z);
        pp_prob.query("U0_Pert_Mag", m_deltaU);
        pp_prob.query("V0_Pert_Mag", m_deltaV);
    }

    void operator()(
        const amrex::Box& bx,
        const amrex::Array4<amrex::Real>& state,
        const amrex::Array4<amrex::Real>& x_vel,
        const amrex::Array4<amrex::Real>& y_vel,
        const amrex::Array4<amrex::Real>& z_vel,
        const amrex::GeometryData& geom) const

    {
       const amrex::Real pi = PI;
       const auto& dx = geom.CellSize();
       const auto& problo = geom.ProbLo();
       const auto& probhi = geom.ProbHi();

       const bool perturb_vel = m_perturb_vel;
       const amrex::Real rho_init = m_rho;
       const amrex::Real umean = m_vel_x;
       const amrex::Real vmean = m_vel_y;
       const amrex::Real wmean = m_vel_z;
       const amrex::Real aval = m_Uperiods * 2.0 * pi / (probhi[1] - problo[1]);
       const amrex::Real bval = m_Vperiods * 2.0 * pi / (probhi[0] - problo[0]);
       const amrex::Real ufac = m_deltaU * std::exp(0.5) / m_ref_height;
       const amrex::Real vfac = m_deltaV * std::exp(0.5) / m_ref_height;
       const amrex::Real ref_height = m_ref_height;

       const amrex::Real tke_init = m_tke_init;
       const amrex::Real qke_init = m_qke_init;

       const int ntvals = m_theta_heights.size();
       const amrex::Real* th = m_thht_d.data();
       const amrex::Real* tv = m_thvv_d.data();

       amrex::ParallelFor(bx, [=] AMREX_GPU_DEVICE(int i, int j, int k) noexcept {
           const amrex::Real z = problo[2] + (k + 0.5) * dx[2];

           state(i, j, k, Rho_comp) = rho_init;

           amrex::Real theta = tv[0];
           for (int iz = 0; iz < ntvals - 1; ++iz) {
               if ((z > th[iz]) && (z <= th[iz + 1])) {
                   const amrex::Real slope =
                       (tv[iz + 1] - tv[iz]) / (th[iz + 1] - th[iz]);
                   theta = tv[iz] + (z - th[iz]) * slope;
               }
           }

           state(i, j, k, RhoTheta_comp) = theta;
           state(i, j, k, RhoTheta_comp) *= state(i, j, k, Rho_comp);
           state(i, j, k, RhoScalar_comp) = 0.;

           state(i, j, k, RhoKE_comp) = tke_init;
           state(i, j, k, RhoQKE_comp) = qke_init;
       });

     if (add_temperature_perturbations())
         perturb_temperature(bx,state,geom);

     const amrex::Box& xbx = amrex::surroundingNodes(bx,0);
     amrex::ParallelFor(xbx, [=] AMREX_GPU_DEVICE(int i, int j, int k) noexcept
         {
           x_vel(i, j, k) = umean;

           if (perturb_vel) {
               const amrex::Real y = problo[1] + (j + 0.5) * dx[1];
               const amrex::Real z = problo[2] + (k + 0.5) * dx[2];
               const amrex::Real yl = y - problo[1];
               const amrex::Real zl = z / ref_height;
               const amrex::Real damp = std::exp(-0.5 * zl * zl);

               x_vel(i, j, k) += ufac * damp * z * std::cos(aval * yl);
           }
         });

     const amrex::Box& ybx = amrex::surroundingNodes(bx,1);
     amrex::ParallelFor(ybx, [=] AMREX_GPU_DEVICE(int i, int j, int k) noexcept
         {
           y_vel(i, j, k) = vmean;

           if (perturb_vel) {
               const amrex::Real x = problo[0] + (i + 0.5) * dx[0];
               const amrex::Real z = problo[2] + (k + 0.5) * dx[2];
               const amrex::Real xl = x - problo[0];
               const amrex::Real zl = z / ref_height;
               const amrex::Real damp = std::exp(-0.5 * zl * zl);

               y_vel(i, j, k) += vfac * damp * z * std::cos(bval * xl);
        }
      });

      const amrex::Box& zbx = amrex::surroundingNodes(bx,2);
      amrex::ParallelFor(zbx, [=] AMREX_GPU_DEVICE(int i, int j, int k) noexcept
      {
          z_vel(i, j, k) = wmean;
      });
    }

    /** Add temperature perturbations
     *
     *  This uses amrex::Random and, therefore, executes the loop on CPU and
     *  pushes the field to device.
     */

    void perturb_temperature(const amrex::Box& bx,
                             const amrex::Array4<amrex::Real>& state,
                             const amrex::GeometryData& geom) const
    {
    /** Perturbations for the temperature field is adapted from the following
     * paper:
     *
     *  D. Munoz-Esparza, B. Kosovic, J. van Beeck, J. D. Mirocha, A stocastic
     *  perturbation method to generate inflow turbulence in large-eddy
     *  simulation models: Application to neutrally stratified atmospheric
     *  boundary layers. Physics of Fluids, Vol. 27, 2015.
     *
     */

    const auto& dx     = geom.CellSize();
    const auto& problo = geom.ProbLo();
    const auto theta_cutoff_height = m_theta_cutoff_height;
    const auto theta_gauss_mean = m_theta_gauss_mean;
    const auto theta_gauss_var = m_theta_gauss_var;
    const auto deltaT = m_deltaT;

    amrex::ParallelForRNG(
        bx, [=] AMREX_GPU_DEVICE(
            int i, int j, int k,
            const amrex::RandomEngine& engine) noexcept {
            const amrex::Real z = problo[2] + (k + 0.5) * dx[2];
            if (z < theta_cutoff_height)
                state(i, j, k, RhoTheta_comp) =
                    deltaT * amrex::RandomNormal(
                                 theta_gauss_mean, theta_gauss_var, engine);
        });
    }

    //! Flag indicating whether temperature field needs perturbations
    bool add_temperature_perturbations() const { return m_perturb_theta; }

    //! Initialize TKE field
    void init_tke(const amrex::Geometry& geom, amrex::MultiFab& tke) const;

  private:
    //! Initial velocity components
    amrex::Real m_vel_x;
    amrex::Real m_vel_y;
    amrex::Real m_vel_z;

    ///@{
    //! Lookup table containing temperature variation w.r.t. height above
    //! terrain
    amrex::Vector<amrex::Real> m_theta_heights;
    amrex::Vector<amrex::Real> m_theta_values;
    ///@}

    // Device copies of the above arrays
    amrex::Gpu::DeviceVector<amrex::Real> m_thht_d;
    amrex::Gpu::DeviceVector<amrex::Real> m_thvv_d;

    //! Initial density field
    amrex::Real m_rho;

    //! Amplitude of fluctuations in x-direction
    amrex::Real m_deltaU{1.0};

    //! Amplitude of fluctuations in y-direction
    amrex::Real m_deltaV{1.0};

    //! Number of sinusoidal waves in x-direction
    amrex::Real m_Uperiods{4.0};

    //! Number of sinusoidal waves in y-direction
    amrex::Real m_Vperiods{4.0};

    //! Reference height for velocity perturbations
    amrex::Real m_ref_height{50.0};

    //! Amplitude of temperature perturbations
    amrex::Real m_deltaT{0.8};

    //! Mean for Gaussian number generator
    amrex::Real m_theta_gauss_mean{0.0};

    //! Variance for Gaussian number generator
    amrex::Real m_theta_gauss_var{1.0};

    //! Cutoff height for temperature fluctuations
    amrex::Real m_theta_cutoff_height{1.0e16};

    //! Initial value for tke field
    amrex::Real m_tke_init{0.1};
    amrex::Real m_qke_init{0.2};

    //! Perturb initial velocity field with sinusoidal fluctuations
    bool m_perturb_vel{true};

    //! Perturb temperature field with random fluctuations
    bool m_perturb_theta{false};
};

#endif /* ABLFIELDINIT_H */
