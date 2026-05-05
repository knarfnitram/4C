// This file is part of 4C multiphysics licensed under the
// GNU Lesser General Public License v3.0 or later.
//
// See the LICENSE.md file in the top-level for license information.
//
// SPDX-License-Identifier: LGPL-3.0-or-later

// Reduced SMA beam material for Simo-Reissner beams.
// Includes prescribed-temperature support, phase-dependent stiffness,
// and optional numerical consistent tangent.

#ifndef FOUR_C_MAT_BEAM3R_SMA_HPP
#define FOUR_C_MAT_BEAM3R_SMA_HPP

#include "4C_config.hpp"

#include "4C_comm_parobjectfactory.hpp"
#include "4C_linalg_fixedsizematrix.hpp"
#include "4C_mat_beam_elasthyper.hpp"
#include "4C_mat_beam_elasthyper_parameter.hpp"
#include "4C_mat_material_factory.hpp"
#include "4C_material_base.hpp"
#include "4C_utils_local_newton.hpp"

#include <array>
#include <functional>
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
      explicit BeamReissnerSMAMaterialParams(const Core::Mat::PAR::Parameter::Data& matdata);

      std::shared_ptr<Core::Mat::Material> create_material() override;

      // Temperature.  TEMP is used unless the material object is provided with a
      // time dependent temperature function.
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

      double get_fb_regularization() const { return fb_regularization_; }
      double get_rs_regularization() const { return rs_regularization_; }
      double get_local_newton_tol() const { return local_newton_tol_; }
      int get_local_newton_maxiter() const { return local_newton_maxiter_; }

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

      double fb_regularization_;
      double rs_regularization_;
      double local_newton_tol_;
      int local_newton_maxiter_;
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

    // Optional hook.  Call this from the element before constitutive evaluation
    // if the prescribed temperature function depends on time.
    void set_current_time(double time) { current_time_ = time; }

    // Optional hook.  If not set, the constant material parameter TEMP is used.
    void set_temperature_function(std::function<double(double)> temperature_function)
    {
      temperature_function_ = temperature_function;
    }

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
    using Vec7 = Core::LinAlg::Matrix<7, 1, double>;
    using Mat7 = Core::LinAlg::Matrix<7, 7, double>;

    struct LocalState
    {
      T vM{0.0};
      T vS{0.0};
      Core::LinAlg::Matrix<3, 1, T> dir{Core::LinAlg::Initialization::zero};
    };

    const Mat::PAR::BeamReissnerSMAMaterialParams& sma_params() const;

    static double sign(double x);
    static double abs_reg(double a, double eps);
    static double fb(double a, double b, double eps);
    static double vec_norm(const Vec3& v);
    static double dot(const Vec3& a, const Vec3& b);
    double vector_difference_norm(const Vec3& a, const Vec3& b) const;
    template <typename ResidualFunc>
    void numerical_jacobian(ResidualFunc&& residual, const Vec7& x, const Vec7& R, Mat7& J) const;

    template <typename ResidualFunc>
    bool solve_nonlinear_system(ResidualFunc&& residual, Vec7& x) const;

    double current_temperature() const;
    double compute_rs(double bs, double vS) const;

    void compute_phase_constitutive_matrices(
        double vM, double vS, Mat3& C_N, Mat3& C_M, Mat3& dC_N, Mat3& dC_M) const;
    double quadratic_form(const Vec3& x, const Mat3& A) const;

    void choose_force_direction(const Vec3& Gamma, const LocalState& old_state, Vec3& dir) const;
    void choose_moment_direction(const Vec3& Cur, const LocalState& old_state, Vec3& dir) const;

    void compute_force_transformation_strain(double vS, const Vec3& dir, Vec3& gamma_tr) const;
    void compute_moment_transformation_curvature(double vS, const Vec3& dir, Vec3& kappa_tr) const;

    void assemble_force_residual(const Vec7& x, const LocalState& old_state, const Vec3& Gamma,
        const Vec3& dir, Vec7& R, Vec3* stress = nullptr) const;

    void assemble_moment_residual(const Vec7& x, const LocalState& old_state, const Vec3& Cur,
        const Vec3& dir, Vec7& R, Vec3* stress = nullptr) const;

    bool solve_force_state(
        const Vec3& Gamma, const LocalState& old_state, LocalState& new_state, Vec3& stress) const;

    bool solve_moment_state(
        const Vec3& Cur, const LocalState& old_state, LocalState& new_state, Vec3& stress) const;

    bool integrate_force_state_substepped(const Vec3& Gamma_start, const Vec3& Gamma_end,
        const LocalState& state_start, LocalState& state_end, Vec3& stress_end, int& nsub) const;

    bool integrate_moment_state_substepped(const Vec3& Cur_start, const Vec3& Cur_end,
        const LocalState& state_start, LocalState& state_end, Vec3& stress_end, int& nsub) const;

    void local_force_response(const Vec3& Gamma, const unsigned int gp, Vec3& stress, Mat3& C_alg);

    void local_moment_response(const Vec3& Cur, const unsigned int gp, Vec3& stress, Mat3& C_alg);

    std::vector<LocalState> force_state_conv_;
    std::vector<LocalState> force_state_new_;
    std::vector<LocalState> moment_state_conv_;
    std::vector<LocalState> moment_state_new_;

    std::vector<Vec3> gamma_conv_;
    std::vector<Vec3> gamma_new_;

    std::vector<Vec3> cur_conv_;
    std::vector<Vec3> cur_new_;


    std::vector<Mat3> c_n_alg_;
    std::vector<Mat3> c_m_alg_;

    unsigned int numgp_force_{0};
    unsigned int numgp_moment_{0};

    double current_time_{0.0};
    std::function<double(double)> temperature_function_;
  };
}  // namespace Mat

FOUR_C_NAMESPACE_CLOSE

#endif