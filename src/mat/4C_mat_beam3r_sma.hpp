// This file is part of 4C multiphysics licensed under the
// GNU Lesser General Public License v3.0 or later.
//
// See the LICENSE.md file in the top-level for license information.
//
// SPDX-License-Identifier: LGPL-3.0-or-later

// Reduced isothermal SMA beam material for Simo-Reissner beams.
// Prototype implementation based on the simplified Auricchio model
// with proportional-loading reduction, axial SMA + bending SMA,
// shear elastic, torsion elastic by default, no FAD support.

#ifndef FOUR_C_MAT_BEAM3R_SMA_HPP
#define FOUR_C_MAT_BEAM3R_SMA_HPP

#include "4C_config.hpp"

#include "4C_comm_parobjectfactory.hpp"
#include "4C_linalg_fixedsizematrix.hpp"
#include "4C_mat_beam_elasthyper.hpp"
#include "4C_mat_beam_elasthyper_parameter.hpp"
#include "4C_mat_material_factory.hpp"
#include "4C_material_base.hpp"

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

      double get_temperature() const { return temperature_; }
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


      double determine_shear_modulus(const Core::Mat::PAR::Parameter::Data& matdata) { return 0; }

      double get_axial_rigidity() const override
      {
        FOUR_C_THROW("BeamReissnerSMAMaterialParams has no get_axial_rigidity.");
        return -1.0;
      }
      double get_shear_rigidity2() const override
      {
        FOUR_C_THROW("BeamReissnerSMAMaterialParams has no get_shear_rigidity2.");
        return -1.0;
      }

      double get_shear_rigidity3() const override
      {
        FOUR_C_THROW("BeamReissnerSMAMaterialParams has no get_shear_rigidity3.");
        return -1.0;
      }

      double get_torsional_rigidity() const override
      {
        FOUR_C_THROW("BeamReissnerSMAMaterialParams has no get_torsional_rigidity.");
        return -1.0;
      }

      double get_bending_rigidity2() const override
      {
        FOUR_C_THROW("BeamReissnerSMAMaterialParams has no get_bending_rigidity2. ");
        return -1.0;
      }

      double get_bending_rigidity3() const override
      {
        FOUR_C_THROW("BeamReissnerSMAMaterialParams has no get_bending_rigidity3.");
        return -1.0;
      }

      double get_translational_mass_inertia() const override
      {
        FOUR_C_THROW("BeamReissnerSMAMaterialParams has no plastic yield moment. ");
        return -1.0;
      }

      double get_polar_mass_moment_of_inertia() const override
      {
        FOUR_C_THROW("BeamReissnerSMAMaterialParams has no get_polar_mass_moment_of_inertia  ");
        return -1.0;
      }

      double get_mass_moment_of_inertia2() const override
      {
        FOUR_C_THROW("BeamReissnerSMAMaterialParams has no get_mass_moment_of_inertia2 ");
        return -1.0;
      }

      double get_mass_moment_of_inertia3() const override
      {
        FOUR_C_THROW("BeamReissnerSMAMaterialParams has no get_mass_moment_of_inertia3");
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

      double radius_interaction_;
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

    /*Core::Materials::MaterialType material_type() const override
    {
      return Core::Materials::m_beam_reissner_sma;
    }*/

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

    void get_stiffness_matrix_of_forces(Core::LinAlg::Matrix<3, 3, T>& stiffness_matrix,
        const Core::LinAlg::Matrix<3, 3, T>& C_N, const int gp) override;

    void get_stiffness_matrix_of_moments(Core::LinAlg::Matrix<3, 3, T>& stiffness_matrix,
        const Core::LinAlg::Matrix<3, 3, T>& C_M, const int gp) override;

    void update() override;
    void reset() override;

   private:
    using Vec3 = Core::LinAlg::Matrix<3, 1, T>;
    using Mat3 = Core::LinAlg::Matrix<3, 3, T>;
    using Vec7 = std::array<double, 7>;
    using Mat7 = std::array<std::array<double, 7>, 7>;

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
    static double residual_norm(const Vec7& r);
    static bool solve_linear_7x7(Mat7& A, Vec7& b);

    template <typename ResidualFunc>
    void numerical_jacobian(ResidualFunc&& residual, const Vec7& x, const Vec7& R, Mat7& J) const;

    template <typename ResidualFunc>
    bool solve_nonlinear_system(ResidualFunc&& residual, Vec7& x) const;

    double compute_rs(double bs, double vS) const;

    void choose_force_direction(const Vec3& Gamma, const LocalState& old_state, Vec3& dir) const;
    void choose_moment_direction(
        const Vec3& Cur, const Mat3& C_M, const LocalState& old_state, Vec3& dir) const;

    void compute_force_transformation_strain(double vS, const Vec3& dir, Vec3& gamma_tr) const;
    void compute_moment_transformation_curvature(double vS, const Vec3& dir, Vec3& kappa_tr) const;

    void assemble_force_residual(const Vec7& x, const LocalState& old_state, const Vec3& Gamma,
        const Mat3& C_N, const Vec3& dir, Vec7& R, Vec3* stress = nullptr) const;

    void assemble_moment_residual(const Vec7& x, const LocalState& old_state, const Vec3& Cur,
        const Mat3& C_M, const Vec3& dir, Vec7& R, Vec3* stress = nullptr) const;

    bool solve_force_state(const Vec3& Gamma, const Mat3& C_N, const LocalState& old_state,
        LocalState& new_state, Vec3& stress) const;

    bool solve_moment_state(const Vec3& Cur, const Mat3& C_M, const LocalState& old_state,
        LocalState& new_state, Vec3& stress) const;

    void local_force_response(
        const Vec3& Gamma, const Mat3& C_N, const unsigned int gp, Vec3& stress, Mat3& C_alg);

    void local_moment_response(
        const Vec3& Cur, const Mat3& C_M, const unsigned int gp, Vec3& stress, Mat3& C_alg);

    std::vector<LocalState> force_state_conv_;
    std::vector<LocalState> force_state_new_;
    std::vector<LocalState> moment_state_conv_;
    std::vector<LocalState> moment_state_new_;

    std::vector<Mat3> c_n_alg_;
    std::vector<Mat3> c_m_alg_;

    unsigned int numgp_force_{0};
    unsigned int numgp_moment_{0};
  };
}  // namespace Mat

FOUR_C_NAMESPACE_CLOSE

#endif