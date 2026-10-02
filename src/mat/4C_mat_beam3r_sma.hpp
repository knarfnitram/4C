// This file is part of 4C multiphysics licensed under the
// GNU Lesser General Public License v3.0 or later.
//
// See the LICENSE.md file in the top-level for license information.
//
// SPDX-License-Identifier: LGPL-3.0-or-later

// Superelastic shape-memory alloy (Nitinol) material for Simo-Reissner beams. Reduced
// Auricchio-type model with multiple-variant (vM) and single-variant (vS) martensite fractions at
// constant temperature, phase-dependent stiffness, a return-mapping update of the phase fractions
// and the analytic consistent tangent.
//
// The axial force is computed from a single state per cross-section (exact for uniform axial
// strain). Bending and torsion are either elastic, use a single sectional state (KAPPALM), or are
// obtained by integrating the material law over fibres of a circular cross-section, each carrying
// an axial strain and a torsional shear strain.

#ifndef FOUR_C_MAT_BEAM3R_SMA_HPP
#define FOUR_C_MAT_BEAM3R_SMA_HPP

#include "4C_config.hpp"

#include "4C_comm_parobjectfactory.hpp"
#include "4C_linalg_fixedsizematrix.hpp"
#include "4C_mat_beam_elasthyper.hpp"
#include "4C_mat_beam_elasthyper_parameter.hpp"
#include "4C_mat_material_factory.hpp"
#include "4C_material_base.hpp"

#include <memory>
#include <vector>

FOUR_C_NAMESPACE_OPEN

namespace Discret
{
  class ParObject;
}

namespace Mat
{
  namespace PAR
  {
    class BeamReissnerSMAMaterialParams : public BeamElastHyperMaterialParameterGeneric
    {
     public:
      //! model for the bending and torsion response
      enum class MomentModel
      {
        elastic,    ///< no transformation from moments
        sectional,  ///< one SMA state per cross-section with transformation curvature KAPPALM
        fibers      ///< SMA law integrated over fibres of a circular cross-section
      };

      explicit BeamReissnerSMAMaterialParams(const Core::Mat::PAR::Parameter::Data& matdata);

      std::shared_ptr<Core::Mat::Material> create_material() override;

      // constant material temperature
      double get_temperature() const { return temperature_; }

      // Simplified Auricchio SMA parameters.
      double get_delta_s_ams() const { return delta_s_ams_; }
      double get_t0_sma() const { return t0_sma_; }
      double get_w_in() const { return w_in_; }
      double get_r_m() const { return r_m_; }

      double get_rs_f0() const { return rs_f0_; }
      double get_rs_r0() const { return rs_r0_; }
      double get_h_sf() const { return h_sf_; }
      double get_h_sr() const { return h_sr_; }
      double get_c_ts() const { return c_ts_; }
      double get_a_sf0() const { return a_sf0_; }
      double get_a_sf1() const { return a_sf1_; }
      double get_a_sr0() const { return a_sr0_; }
      double get_a_sr1() const { return a_sr1_; }
      double get_n_exp() const { return n_exp_; }

      double get_eps_l_n() const { return eps_l_n_; }
      double get_kappa_l_m() const { return kappa_l_m_; }
      bool get_torsion_sma() const { return torsion_sma_; }

      MomentModel get_moment_model() const { return moment_model_; }
      int get_num_fiber_rings() const { return num_fiber_rings_; }
      int get_num_fiber_sectors() const { return num_fiber_sectors_; }

      double get_rs_regularization() const { return rs_regularization_; }
      int get_local_staggered_maxiter() const { return local_staggered_maxiter_; }

      int get_local_substep_max() const { return local_substep_max_; }
      double get_local_substep_strain_limit() const { return local_substep_strain_limit_; }
      double get_local_substep_curvature_limit() const { return local_substep_curvature_limit_; }

      // Phase-dependent elastic parameters. YOUNG and either POISSONRATIO or SHEARMOD are
      // interpreted as austenite properties; YOUNGMART/SHEARMODMART are martensite properties.
      double get_youngs_modulus_austenite() const { return youngs_modulus_austenite_; }
      double get_youngs_modulus_martensite() const { return youngs_modulus_martensite_; }
      double get_shear_modulus_austenite() const { return shear_modulus_austenite_; }
      double get_shear_modulus_martensite() const { return shear_modulus_martensite_; }

      // Geometry/inertia data used directly by BeamSMAMaterial. The generic modal stiffness
      // getters intentionally throw because stiffness is phase-dependent and must be built in the
      // material object from the current Gauss-point phase fractions.
      double get_sma_density() const { return density_; }
      double get_sma_cross_section_area() const { return cross_section_area_; }
      double get_sma_shear_correction_factor() const { return shear_correction_factor_; }
      double get_sma_moment_inertia2() const { return area_moment_inertia_2_; }
      double get_sma_moment_inertia3() const { return area_moment_inertia_3_; }
      double get_sma_moment_inertia_polar() const { return area_moment_inertia_polar_; }

      bool use_consistent_tangent() const { return use_consistent_tangent_; }

      double get_axial_rigidity() const override
      {
        FOUR_C_THROW(
            "BeamReissnerSMAMaterialParams has no constant axial rigidity; use the phase-dependent "
            "SMA stiffness.");
        return -1.0;
      }

      double get_shear_rigidity2() const override
      {
        FOUR_C_THROW(
            "BeamReissnerSMAMaterialParams has no constant shear rigidity; use the phase-dependent "
            "SMA stiffness.");
        return -1.0;
      }

      double get_shear_rigidity3() const override
      {
        FOUR_C_THROW(
            "BeamReissnerSMAMaterialParams has no constant shear rigidity; use the phase-dependent "
            "SMA stiffness.");
        return -1.0;
      }

      double get_torsional_rigidity() const override
      {
        FOUR_C_THROW(
            "BeamReissnerSMAMaterialParams has no constant torsional rigidity; use the "
            "phase-dependent SMA stiffness.");
        return -1.0;
      }

      double get_bending_rigidity2() const override
      {
        FOUR_C_THROW(
            "BeamReissnerSMAMaterialParams has no constant bending rigidity; use the "
            "phase-dependent SMA stiffness.");
        return -1.0;
      }

      double get_bending_rigidity3() const override
      {
        FOUR_C_THROW(
            "BeamReissnerSMAMaterialParams has no constant bending rigidity; use the "
            "phase-dependent SMA stiffness.");
        return -1.0;
      }

      double get_translational_mass_inertia() const override
      {
        FOUR_C_THROW(
            "BeamReissnerSMAMaterialParams should not provide mass through the generic getter; "
            "BeamSMAMaterial handles it directly.");
        return -1.0;
      }

      double get_polar_mass_moment_of_inertia() const override
      {
        FOUR_C_THROW(
            "BeamReissnerSMAMaterialParams should not provide rotational inertia through the "
            "generic getter; BeamSMAMaterial handles it directly.");
        return -1.0;
      }

      double get_mass_moment_of_inertia2() const override
      {
        FOUR_C_THROW(
            "BeamReissnerSMAMaterialParams should not provide rotational inertia through the "
            "generic getter; BeamSMAMaterial handles it directly.");
        return -1.0;
      }

      double get_mass_moment_of_inertia3() const override
      {
        FOUR_C_THROW(
            "BeamReissnerSMAMaterialParams should not provide rotational inertia through the "
            "generic getter; BeamSMAMaterial handles it directly.");
        return -1.0;
      }

      double get_interaction_radius() const override
      {
        if (radius_interaction_ == -1.0)
          FOUR_C_THROW(
              "the radius of a beam which is to be used for interactions (contact, potentials, "
              "viscous drag in background fluid ...) has not been specified in the material "
              "definition!");
        return radius_interaction_;
      }

      // The following plasticity accessors must not be used by the SMA material.
      double get_yield_stress_n() const override
      {
        FOUR_C_THROW(
            "BeamReissnerSMAMaterialParams has no plastic yield force. "
            "Use SMA transformation parameters RSF0/RSR0/HSF/HSR instead.");
        return -1.0;
      }

      double get_yield_stress_m() const override
      {
        FOUR_C_THROW(
            "BeamReissnerSMAMaterialParams has no plastic yield moment. "
            "Use SMA transformation parameters and KAPPALM instead.");
        return -1.0;
      }

      double get_hardening_axial_rigidity() const override
      {
        FOUR_C_THROW("BeamReissnerSMAMaterialParams has no plastic axial hardening rigidity.");
        return -1.0;
      }

      double get_hardening_shear_rigidity2() const override
      {
        FOUR_C_THROW("BeamReissnerSMAMaterialParams has no plastic shear hardening rigidity.");
        return -1.0;
      }

      double get_hardening_shear_rigidity3() const override
      {
        FOUR_C_THROW("BeamReissnerSMAMaterialParams has no plastic shear hardening rigidity.");
        return -1.0;
      }

      double get_hardening_momental_rigidity() const override
      {
        FOUR_C_THROW("BeamReissnerSMAMaterialParams has no plastic moment hardening rigidity.");
        return -1.0;
      }

      bool get_torsion_plasticity() const override
      {
        FOUR_C_THROW(
            "BeamReissnerSMAMaterialParams has no torsion plasticity flag. "
            "Use TORSIONSMA instead.");
        return false;
      }

     private:
      double temperature_;
      double delta_s_ams_;
      double t0_sma_;
      double w_in_;
      double r_m_;

      double rs_f0_;
      double rs_r0_;
      double h_sf_;
      double h_sr_;
      double c_ts_;
      double a_sf0_;
      double a_sf1_;
      double a_sr0_;
      double a_sr1_;
      double n_exp_;

      double eps_l_n_;
      double kappa_l_m_;
      bool torsion_sma_;

      MomentModel moment_model_;
      int num_fiber_rings_;
      int num_fiber_sectors_;

      double rs_regularization_;
      int local_staggered_maxiter_;
      int local_substep_max_;
      double local_substep_strain_limit_;
      double local_substep_curvature_limit_;

      double youngs_modulus_austenite_;
      double youngs_modulus_martensite_;
      double shear_modulus_austenite_;
      double shear_modulus_martensite_;

      double density_;
      double cross_section_area_;
      double shear_correction_factor_;
      double area_moment_inertia_polar_;
      double area_moment_inertia_2_;
      double area_moment_inertia_3_;
      double radius_interaction_;

      bool use_consistent_tangent_;
    };
  }  // namespace PAR

  template <typename T>
  class BeamSMAMaterialType : public Core::Communication::ParObjectType
  {
   public:
    std::string name() const override { return typeid(this).name(); }

    static BeamSMAMaterialType& instance() { return instance_; }

    Core::Communication::ParObject* create(Core::Communication::UnpackBuffer& buffer) override;

   private:
    static BeamSMAMaterialType instance_;
  };

  template <typename T>
  class BeamSMAMaterial : public BeamElastHyperMaterial<T>
  {
   public:
    BeamSMAMaterial() = default;

    explicit BeamSMAMaterial(Mat::PAR::BeamReissnerSMAMaterialParams* params);

    void setup(int numgp_force, int numgp_moment) override;

    int unique_par_object_id() const override
    {
      return BeamSMAMaterialType<T>::instance().unique_par_object_id();
    }

    void pack(Core::Communication::PackBuffer& data) const override;
    void unpack(Core::Communication::UnpackBuffer& buffer) override;

    std::shared_ptr<Core::Mat::Material> clone() const override
    {
      return std::make_shared<BeamSMAMaterial>(*this);
    }

    void evaluate_force_contributions_to_stress(Core::LinAlg::Matrix<3, 1, T>& stressN,
        const Core::LinAlg::Matrix<3, 3, T>& C_N, const Core::LinAlg::Matrix<3, 1, T>& Gamma,
        const unsigned int gp) override;

    void evaluate_moment_contributions_to_stress(Core::LinAlg::Matrix<3, 1, T>& stressM,
        const Core::LinAlg::Matrix<3, 3, T>& C_M, const Core::LinAlg::Matrix<3, 1, T>& Cur,
        const unsigned int gp) override;

    void compute_constitutive_parameter(
        Core::LinAlg::Matrix<3, 3, T>& C_N, Core::LinAlg::Matrix<3, 3, T>& C_M) override;

    void get_constitutive_matrix_of_forces_material_frame(
        Core::LinAlg::Matrix<3, 3, T>& C_N) const override;

    void get_constitutive_matrix_of_moments_material_frame(
        Core::LinAlg::Matrix<3, 3, T>& C_M) const override;

    void get_stiffness_matrix_of_forces(Core::LinAlg::Matrix<3, 3, T>& stiffness_matrix,
        const Core::LinAlg::Matrix<3, 3, T>& C_N, const int gp) override;

    void get_stiffness_matrix_of_moments(Core::LinAlg::Matrix<3, 3, T>& stiffness_matrix,
        const Core::LinAlg::Matrix<3, 3, T>& C_M, const int gp) override;

    double get_translational_mass_inertia_factor() const override;

    void get_mass_moment_of_inertia_tensor_material_frame(
        Core::LinAlg::Matrix<3, 3>& J) const override;

    void get_mass_moment_of_inertia_tensor_material_frame(
        Core::LinAlg::Matrix<3, 3, Sacado::Fad::DFad<double>>& J) const override;

    double get_interaction_radius() const override;

    void update() override;
    void reset() override;

   private:
    using Vec3 = Core::LinAlg::Matrix<3, 1, T>;
    using Mat3 = Core::LinAlg::Matrix<3, 3, T>;

    //! the two independent sets of stress resultants with their own SMA state
    //! sets of strain measures with their own SMA state: q = Gamma (force), q = K (sectional
    //! moment) or q = (epsilon, gamma_12, gamma_13) of a single fibre
    enum class Mode
    {
      force,
      moment,
      fiber
    };

    //! position and area weight of a fibre in the cross-section
    struct Fiber
    {
      double y;
      double z;
      double weight;
    };

    //! outcome of the update of a single phase fraction within a step
    enum class Transformation
    {
      none,
      forward,
      reverse,
      lower_bound,
      upper_bound
    };

    struct LocalState
    {
      T vM{0.0};
      T vS{0.0};
      Core::LinAlg::Matrix<3, 1, T> dir{Core::LinAlg::Initialization::zero};
    };

    //! thermodynamic driving forces B_i = -dpsi/dv_i and their derivatives; index order [S, M]
    struct DrivingForces
    {
      double B[2] = {0.0, 0.0};
      double dB_dv[2][2] = {{0.0, 0.0}, {0.0, 0.0}};
      Vec3 dB_dq[2];
      Vec3 dB_dd[2];
      Vec3 stress;
      Vec3 dstress_dv[2];
      Mat3 dstress_dd;
    };

    static constexpr unsigned int idx_s = 0;
    static constexpr unsigned int idx_m = 1;

    // Below this single-variant martensite fraction, the direction of the transformation
    // strain/curvature may be re-chosen from the current trial stress resultant.
    static constexpr double direction_reset_tolerance_ = 1.0e-6;

    //! absolute tolerance of the phase fractions in the local update
    static constexpr double phase_fraction_tolerance_ = 1.0e-14;

    //! maximum number of safeguarded Newton iterations for one phase fraction
    static constexpr unsigned int max_root_iterations_ = 200;

    const Mat::PAR::BeamReissnerSMAMaterialParams& sma_params() const;

    static double sign(double x);
    static double vec_norm(const Vec3& v);
    static double dot(const Vec3& a, const Vec3& b);
    double vector_difference_norm(const Vec3& a, const Vec3& b) const;

    double current_temperature() const;

    //! transformation strain amplitude EPSLN (force) or curvature amplitude KAPPALM (moment)
    double transformation_amplitude(Mode mode) const;

    //! radius of the elastic domain of phase fraction i and its derivative w.r.t. v_i
    double transformation_radius(
        Mode mode, unsigned int i, bool forward, double v, double& dR_dv) const;

    //! factor converting the volume-specific thermodynamic parameters to the energy measure of
    //! the mode: cross-section area for the sectional modes, one for a fibre
    double energy_scaling(Mode mode) const;

    //! constitutive matrix and its first two derivatives w.r.t. the martensite fraction m
    void compute_phase_constitutive_matrices(
        double m, Mode mode, Mat3& C, Mat3& dC, Mat3& d2C) const;

    void evaluate_driving_forces(
        Mode mode, const Vec3& q, double vM, double vS, const Vec3& dir, DrivingForces& df) const;

    //! g = B_i - s R_i (s = +1 forward, -1 reverse) with v_i = v and its derivative w.r.t. v
    double consistency_function(Mode mode, const Vec3& q, const Vec3& dir, unsigned int i,
        double v_other, bool forward, double v, double& dg_dv) const;

    //! root of the consistency function on the bracket g(v_pos) > 0 > g(v_neg)
    double find_consistent_phase_fraction(Mode mode, const Vec3& q, const Vec3& dir, unsigned int i,
        double v_other, bool forward, double v_pos, double v_neg) const;

    //! return mapping of phase fraction i for fixed other phase fraction
    Transformation update_phase_fraction(Mode mode, const Vec3& q, const Vec3& dir, unsigned int i,
        double v_other, double v_old, double& v) const;

    void solve_local_state(Mode mode, const Vec3& q, const LocalState& old_state, const Vec3& dir,
        LocalState& new_state, Transformation status[2]) const;

    void choose_force_direction(const Vec3& Gamma, const LocalState& old_state, Vec3& dir) const;
    void choose_moment_direction(const Vec3& Cur, const LocalState& old_state, Vec3& dir) const;
    void choose_fiber_direction(const Vec3& e, const LocalState& old_state, Vec3& dir) const;
    void choose_direction(Mode mode, const Vec3& q, const LocalState& old_state, Vec3& dir) const;

    //! derivative of the transformation direction w.r.t. q if it is chosen from the trial stress
    //! in this (sub)step; returns false if the direction of the old state is kept
    bool direction_derivative(Mode mode, const Vec3& q, const LocalState& old_state,
        const Vec3& dir, Mat3& ddir_dq) const;

    //! (substepped) local update from q_start to q_end including the consistent tangent
    void integrate_local_state(Mode mode, const Vec3& q_start, const Vec3& q_end,
        const LocalState& state_start, LocalState& state_end, Vec3& stress, Mat3& C_alg) const;

    void local_response(Mode mode, const Vec3& q, const unsigned int gp, Vec3& stress, Mat3& C_alg);

    //! bending moments and torque by integration over the fibres of the cross-section
    void local_fiber_moment_response(
        const Vec3& Cur, const unsigned int gp, Vec3& stress, Mat3& C_alg);

    //! fibre layout of the circular cross-section (Gauss rings x equidistant sectors)
    const std::vector<Fiber>& fibers() const;

    //! fibre strains (epsilon, gamma_12, gamma_13) = B K for curvatures K = (K1, K2, K3)
    static void fiber_strains(const Fiber& fiber, const Vec3& Cur, Vec3& e);

    bool uses_fibers() const;

    std::vector<LocalState> force_state_conv_;
    std::vector<LocalState> force_state_new_;
    std::vector<LocalState> moment_state_conv_;
    std::vector<LocalState> moment_state_new_;

    //! fibre states of all moment Gauss points, index gp * num_fibers + fiber
    std::vector<LocalState> fiber_state_conv_;
    std::vector<LocalState> fiber_state_new_;

    mutable std::vector<Fiber> fibers_;

    std::vector<Vec3> gamma_conv_;
    std::vector<Vec3> gamma_new_;

    std::vector<Vec3> cur_conv_;
    std::vector<Vec3> cur_new_;


    std::vector<Mat3> c_n_alg_;
    std::vector<Mat3> c_m_alg_;

    unsigned int numgp_force_{0};
    unsigned int numgp_moment_{0};
  };
}  // namespace Mat

FOUR_C_NAMESPACE_CLOSE

#endif