// This file is part of 4C multiphysics licensed under the
// GNU Lesser General Public License v3.0 or later.
//
// See the LICENSE.md file in the top-level for license information.
//
// SPDX-License-Identifier: LGPL-3.0-or-later

// beam_nitinol_material.cpp
#include "4C_mat_beam3r_superelastic.hpp"

#include "4C_comm_pack_helpers.hpp"
#include "4C_global_data.hpp"
#include "4C_mat_par_bundle.hpp"
#include "4C_utils_enum.hpp"
#include "4C_utils_fad.hpp"

#include <cmath>

FOUR_C_NAMESPACE_OPEN

namespace Mat
{
  namespace PAR
  {

    // Parameter constructor
    BeamReissnerNitinolMaterialParams::BeamReissnerNitinolMaterialParams(
        const Core::Mat::PAR::Parameter::Data& matdata)
        : BeamElastHyperMaterialParameterGeneric(matdata),
          E_A_(matdata.parameters.get<double>("YOUNGAUST")),
          E_M_(matdata.parameters.get<double>("YOUNGMART")),
          eps_L_(matdata.parameters.get<double>("EPS_L")),
          sigma_s_(matdata.parameters.get<double>("SIGMA_S")),
          sigma_f_(matdata.parameters.get<double>("SIGMA_F")),
          kappa_L_(matdata.parameters.get<double>("KAPPAL")),
          M_s_(matdata.parameters.get<double>("MSTART")),
          M_f_(matdata.parameters.get<double>("MFINISH")),
          shear_modulus_(determine_shear_modulus(matdata)),
          density_(matdata.parameters.get<double>("DENS")),
          cross_section_area_(matdata.parameters.get<double>("CROSSAREA")),
          shear_correction_factor_(matdata.parameters.get<double>("SHEAR_CORRECTION_FACTOR")),
          martensite_update_step_(matdata.parameters.get<double>("MARTENSITE_UPDATE_STEP")),
          area_moment_inertia_polar_(matdata.parameters.get<double>("MOMINPOL")),
          area_moment_inertia_2_(matdata.parameters.get<double>("MOMIN2")),
          area_moment_inertia_3_(matdata.parameters.get<double>("MOMIN3")),
          radius_interaction_(determine_default_interaction_radius(matdata))
    {
    }
  }  // namespace PAR

  std::shared_ptr<Core::Mat::Material> PAR::BeamReissnerNitinolMaterialParams::create_material()
  {
    return std::make_shared<BeamNitinolMaterial<double>>(this);
  }

  // Material constructor

  template <typename T>
  BeamNitinolMaterial<T>::BeamNitinolMaterial(PAR::BeamReissnerNitinolMaterialParams* params)
      : BeamElastHyperMaterial<T>(params),
        E_A_(params->get_youngs_modulus_austenite()),
        E_M_(params->get_youngs_modulus_martensite()),
        eps_L_(params->get_transformation_strain()),
        sigma_s_(params->get_forward_start_stress()),
        sigma_f_(params->get_reverse_finish_stress()),
        kappa_L_(params->get_bending_transformation_curvature()),
        M_s_(params->get_moment_start()),
        M_f_(params->get_moment_finish()),
        shear_modulus_(params->get_shear_modulus()),
        cross_section_area_(params->get_cross_section_area()),
        shear_correction_factor(params->get_cross_section_area())
  {
  }

  template <typename T>
  void BeamNitinolMaterial<T>::setup(int numgp_force, int numgp_moment)
  {
    xi_.resize(numgp_force, 0.0);
    sigma_gp_.resize(numgp_force, 0.0);
    xi_m_.resize(numgp_moment, 0.0);
    moment_gp_.resize(numgp_moment);
  }

  template <typename T>
  void BeamNitinolMaterial<T>::reset()
  {
    // nothing special for reset at this stage
  }

  template <typename T>
  void BeamNitinolMaterial<T>::update()
  {
  }

  template <typename T>
  void BeamNitinolMaterial<T>::evaluate_force_contributions_to_stress(
      Core::LinAlg::Matrix<3, 1, T>& stressN, const Core::LinAlg::Matrix<3, 3, T>& CN,
      const Core::LinAlg::Matrix<3, 1, T>& Gamma, const unsigned int gp)
  {
    T eps = Gamma(0);
    T xi = xi_[gp];
    T E_eff = (1 - xi) * E_A_ + xi * E_M_;
    T eps_tr = xi * eps_L_;
    T sigma = E_eff * (eps - eps_tr);
    if (std::abs(sigma) > sigma_s_ && xi < 1.0)
      xi = std::min(1.0, xi + 0.01);
    else if (std::abs(sigma) < sigma_f_ && xi > 0.0)
      xi = std::max(0.0, xi - 0.01);
    xi_[gp] = xi;
    E_eff = (1 - xi) * E_A_ + xi * E_M_;
    if (E_eff < 1e-8) std::cerr << "WARNING: E_eff approxx 0 at gp " << gp << std::endl;

    eps_tr = xi * eps_L_;
    sigma = E_eff * (eps - eps_tr);

    sigma_gp_[gp] = sigma;
    stressN(0) = sigma;
    stressN(1) = 0.0;
    stressN(2) = 0.0;
  }

  template <typename T>
  void BeamNitinolMaterial<T>::get_stiffness_matrix_of_moments(
      Core::LinAlg::Matrix<3, 3, T>& stiffM, const Core::LinAlg::Matrix<3, 3, T>& C_M, const int gp)
  {
    // Use current C_M as-is (e.g., computed per GP before)
    stiffM = C_M;
  }

  template <typename T>
  void BeamNitinolMaterial<T>::evaluate_moment_contributions_to_stress(
      Core::LinAlg::Matrix<3, 1, T>& stressM, const Core::LinAlg::Matrix<3, 3, T>& CM,
      const Core::LinAlg::Matrix<3, 1, T>& Cur, const unsigned int gp)
  {
    T kappa_mag = std::sqrt(Cur(1) * Cur(1) + Cur(2) * Cur(2));
    T xi = xi_m_[gp];

    T E_eff = (1.0 - xi) * E_A_ + xi * E_M_;
    T kappa_tr = xi * kappa_L_;
    T kappa_eff = kappa_mag - kappa_tr;
    T M = E_eff * kappa_eff;

    if (std::abs(M) > M_s_ && xi < 1.0)
      xi = std::min(1.0, xi + martensite_update_step_);
    else if (std::abs(M) < M_f_ && xi > 0.0)
      xi = std::max(0.0, xi - martensite_update_step_);

    xi_m_[gp] = xi;
    E_eff = (1.0 - xi) * E_A_ + xi * E_M_;
    kappa_tr = xi * kappa_L_;

    // Avoid division by zero: use normalized form robustly
    const T denom = std::max(kappa_mag, 1e-12);
    stressM(0) = Cur(0) * CM(0, 0);  // torsion unchanged
    stressM(1) = E_eff * (Cur(1) - kappa_tr * Cur(1) / denom);
    stressM(2) = E_eff * (Cur(2) - kappa_tr * Cur(2) / denom);

    moment_gp_[gp] = stressM;
  }

  template <typename T>
  void Mat::BeamNitinolMaterial<T>::get_constitutive_matrix_of_forces_material_frame(
      Core::LinAlg::Matrix<3, 3, T>& C_N) const
  {
    FOUR_C_THROW("Please do not call get_constitutive_matrix_of_forces_material_frame without gp");
  }
  template <typename T>
  void Mat::BeamNitinolMaterial<T>::compute_constitutive_parameter(
      Core::LinAlg::Matrix<3, 3, T>& C_N, Core::LinAlg::Matrix<3, 3, T>& C_M)
  {
    // FOUR_C_THROW("Please do not call compute_constitutive_parameter without gp");
  }
  template <typename T>
  void Mat::BeamNitinolMaterial<T>::compute_constitutive_parameter(
      Core::LinAlg::Matrix<3, 3, T>& C_N, Core::LinAlg::Matrix<3, 3, T>& C_M, int gp)
  {
    // Compute E_eff for the requested Gauss point
    T xi_force = xi_[gp];
    T E_eff_force = (1.0 - xi_force) * E_A_ + xi_force * E_M_;

    T A = this->cross_section_area_;
    T G = this->shear_modulus_;
    T kappa = this->shear_correction_factor;

    C_N.clear();
    C_N(0, 0) = E_eff_force * A;
    C_N(1, 1) = G * A * kappa;
    C_N(2, 2) = G * A * kappa;

    // Bending/torsion
    T xi_m = xi_m_[gp];
    T E_eff_m = (1.0 - xi_m) * E_A_ + xi_m * E_M_;

    C_M.clear();
    C_M(0, 0) = this->params().get_torsional_rigidity();
    C_M(1, 1) = E_eff_m * this->params().get_mass_moment_of_inertia2();
    C_M(2, 2) = E_eff_m * this->params().get_mass_moment_of_inertia3();
  }

}  // namespace Mat

template class Mat::BeamNitinolMaterial<double>;


FOUR_C_NAMESPACE_CLOSE
