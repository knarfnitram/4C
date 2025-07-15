// This file is part of 4C multiphysics licensed under the
// GNU Lesser General Public License v3.0 or later.
//
// See the LICENSE.md file in the top-level for license information.
//
// SPDX-License-Identifier: LGPL-3.0-or-later
#ifndef FOUR_C_MAT_BEAM3R_SUPERELASTIC_HPP
#define FOUR_C_MAT_BEAM3R_SUPERELASTIC_HPP

#include "4C_config.hpp"

#include "4C_linalg_fixedsizematrix.hpp"
#include "4C_mat_beam_elasthyper.hpp"
#include "4C_mat_beam_elasthyper_parameter.hpp"

#include <vector>

FOUR_C_NAMESPACE_OPEN

namespace Mat
{
  namespace PAR
  {
    class BeamReissnerNitinolMaterialParams : public BeamElastHyperMaterialParameterGeneric
    {
     public:
      BeamReissnerNitinolMaterialParams(const Core::Mat::PAR::Parameter::Data& matdata);

      double get_youngs_modulus_austenite() const { return E_A_; }
      double get_youngs_modulus_martensite() const { return E_M_; }
      double get_transformation_strain() const { return eps_L_; }
      double get_forward_start_stress() const { return sigma_s_; }
      double get_reverse_finish_stress() const { return sigma_f_; }
      double get_bending_transformation_curvature() const { return kappa_L_; }
      double get_moment_start() const { return M_s_; }
      double get_moment_finish() const { return M_f_; }

      double get_cross_section_area() const { return cross_section_area_; }
      double get_shear_modulus() const { return shear_modulus_; }
      double get_shear_correction_factor() const { return 1.0; }
      std::shared_ptr<Core::Mat::Material> create_material() override;
      /// @name derived: accessors to 'modal' constitutive parameters
      //@{
      double get_axial_rigidity() const override { return 0; }

      double get_shear_rigidity2() const override
      {
        return shear_modulus_ * cross_section_area_ * shear_correction_factor_;
      }

      double get_shear_rigidity3() const override
      {
        return shear_modulus_ * cross_section_area_ * shear_correction_factor_;
      }


      double get_torsional_rigidity() const override
      {
        return shear_modulus_ * area_moment_inertia_polar_;
      }

      double get_bending_rigidity2() const override { return 1 * E_A_ * area_moment_inertia_2_; }

      double get_bending_rigidity3() const override { return 1 * E_A_ * area_moment_inertia_3_; }


      double get_translational_mass_inertia() const override
      {
        return density_ * cross_section_area_;
      }


      double get_polar_mass_moment_of_inertia() const override
      {
        return density_ * (area_moment_inertia_2_ + area_moment_inertia_3_);
      }

      double get_mass_moment_of_inertia2() const override
      {
        return density_ * area_moment_inertia_2_;
      }

      double get_mass_moment_of_inertia3() const override
      {
        return density_ * area_moment_inertia_3_;
      }
      double mass_moment_of_inertia() const
      {
        return (area_moment_inertia_2_ + area_moment_inertia_3_);
      }

      double moment_of_inertia2() const { return area_moment_inertia_2_; }

      double moment_of_inertia3() const { return area_moment_inertia_3_; }

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
      double E_A_, E_M_, eps_L_, sigma_s_, sigma_f_;
      double kappa_L_, M_s_, M_f_;

      /// shear modulus
      const double shear_modulus_;
      /// mass density
      const double density_;

      /// cross-section area
      const double cross_section_area_;
      /// shear correction factor
      const double shear_correction_factor_;

      const double martensite_update_step_;

      /// polar/axial area moment of inertia
      const double area_moment_inertia_polar_;
      /// area moment of inertia w.r.t. first principal axis of inertia (i.e. second base vector)
      const double area_moment_inertia_2_;
      /// area moment of inertia w.r.t. second principal axis of inertia (i.e. third base vector)
      const double area_moment_inertia_3_;

      const double radius_interaction_;
    };


  }  // namespace PAR
  // Material class
  template <typename T>
  class BeamNitinolMaterial : public BeamElastHyperMaterial<T>
  {
   public:
    BeamNitinolMaterial() = default;
    explicit BeamNitinolMaterial(Mat::PAR::BeamReissnerNitinolMaterialParams* params);

    void setup(int numgp_force, int numgp_moment) override;
    void reset() override;
    void update() override;
    void get_stiffness_matrix_of_forces(Core::LinAlg::Matrix<3, 3, T>& stiffness_matrix,
        const Core::LinAlg::Matrix<3, 3, T>& C_N, const int gp);
    void evaluate_force_contributions_to_stress(Core::LinAlg::Matrix<3, 1, T>& stressN,
        const Core::LinAlg::Matrix<3, 3, T>& CN, const Core::LinAlg::Matrix<3, 1, T>& Gamma,
        const unsigned int gp) override;

    void evaluate_moment_contributions_to_stress(Core::LinAlg::Matrix<3, 1, T>& stressM,
        const Core::LinAlg::Matrix<3, 3, T>& CM, const Core::LinAlg::Matrix<3, 1, T>& Cur,
        const unsigned int gp) override;

    /** \brief get constitutive matrix relating stress force resultants and translational strain
     *         measures, expressed w.r.t. material frame
     */
    void get_constitutive_matrix_of_forces_material_frame(
        Core::LinAlg::Matrix<3, 3, T>& C_N) const override;

    void get_stiffness_matrix_of_moments(Core::LinAlg::Matrix<3, 3, T>& stiffM,
        const Core::LinAlg::Matrix<3, 3, T>& C_M, const Core::LinAlg::Matrix<3, 1, T>& Cur,
        const int gp) override;

    void compute_constitutive_parameter(
        Core::LinAlg::Matrix<3, 3, T>& C_N, Core::LinAlg::Matrix<3, 3, T>& C_M) override;

    void compute_constitutive_parameter(Core::LinAlg::Matrix<3, 3, T>& C_N,
        Core::LinAlg::Matrix<3, 3, T>& C_M, const Core::LinAlg::Matrix<3, 1, T>& Gamma,
        int gp) override;

    T compute_martensite_fraction(T M, T xi_prev) const;
    T compute_dxi_dsigma(T sigma) const;

   protected:
    double get_moment_start() const { return M_s_; }
    double get_moment_finish() const { return M_f_; }



   private:
    std::vector<T> xi_old_;                                 // Axial martensite fraction
    std::vector<T> xi_curr_;                                // current Axial martensite fraction
    std::vector<T> xi_m_;                                   // Bending martensite fraction
    std::vector<T> xi_m_curr_;                              // current Bending martensite prev
    std::vector<T> sigma_gp_;                               // Axial stress history
    std::vector<Core::LinAlg::Matrix<3, 1, T>> moment_gp_;  // Moment vector

    double E_A_, E_M_, eps_L_, sigma_s_, sigma_f_;
    double kappa_L_, M_s_, M_f_;
    double martensite_update_step_;
    double shear_modulus_, cross_section_area_, shear_correction_factor_;
    double torsional_rigidity_;
    double tol_ = 1e-6;

    /// Number of integration points for forces
    unsigned int numgp_force_;

    /// Number of integration points for moments
    unsigned int numgp_moment_;

    void set_parameter(Mat::PAR::BeamReissnerNitinolMaterialParams* parameter)
    {
      params_ = parameter;
    }
    /// my material parameters
    Mat::PAR::BeamReissnerNitinolMaterialParams* params_;
  };
}  // namespace Mat

FOUR_C_NAMESPACE_CLOSE

#endif