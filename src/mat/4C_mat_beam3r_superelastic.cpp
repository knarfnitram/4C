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
        gamma_l_(params->get_transformation_strain()),
        sigma_fAS(300.0),
        sigma_sAS(500),
        sigma_fSA(params->get_reverse_finish_stress()),
        sigma_sSA(200),
        kappa_L_(params->get_bending_transformation_curvature()),
        M_s_(params->get_moment_start()),
        M_f_(params->get_moment_finish()),
        shear_modulus_(params->get_shear_modulus()),
        cross_section_area_(params->get_cross_section_area()),
        shear_correction_factor_(params->get_cross_section_area()),
        torsional_rigidity_(params->get_torsional_rigidity()),
        params_(params)
  {
  }

  template <typename T>
  void BeamNitinolMaterial<T>::setup(int numgp_force, int numgp_moment)
  {
    xi_old_.resize(numgp_force, 0.0);
    xi_curr_.resize(numgp_force, 0.0);
    sigma_gp_.resize(numgp_force, 0.0);
    delta_xi_s_AS.resize(numgp_force, 0.0);
    delta_xi_s_SA.resize(numgp_force, 0.0);
    xi_m_.resize(numgp_moment, 0.0);
    moment_gp_.resize(numgp_moment);
    xi_m_curr_.resize(numgp_moment, 0.0);
    numgp_moment_ = numgp_moment;
    numgp_force_ = numgp_force;
    std::cout << "i was called = " << std::endl;

    for (int gp = 0; gp < numgp_force; gp++)
    {
      xi_old_[gp] = 0.0;
      xi_curr_[gp] = 0.0;
      sigma_gp_[gp] = 0.0;
      delta_xi_s_AS[gp] = 0.0;
      delta_xi_s_SA[gp] = 0.0;
    }
    for (int gp = 0; gp < numgp_force; gp++)
    {
      xi_m_[gp] = 0.0;
      xi_m_curr_[gp] = 0.0;
      moment_gp_[gp] = 0.0;
    }
  }

  template <typename T>
  void BeamNitinolMaterial<T>::reset()
  {
    for (unsigned int gp = 0; gp < numgp_force_; gp++)
    {
      delta_xi_s_AS[gp] = 0;
      delta_xi_s_SA[gp] = 0;
      sigma_gp_[gp] = 0;
    }

    for (unsigned int gp = 0; gp < numgp_moment_; gp++)
    {
      xi_m_curr_[gp] = xi_m_[gp];
    }
  }

  template <typename T>
  void BeamNitinolMaterial<T>::update()
  {
    for (unsigned int gp = 0; gp < numgp_force_; gp++)
    {
      xi_old_[gp] = xi_old_[gp] + delta_xi_s_AS[gp] - delta_xi_s_SA[gp];
      xi_old_[gp] = std::clamp(xi_old_[gp], 0.0, 1.0);
      sigma_gp_[gp] = 0;
    }

    for (unsigned int gp = 0; gp < numgp_moment_; gp++)
    {
      xi_m_[gp] = xi_m_curr_[gp];
    }
  }
  template <typename T>
  T Mat::BeamNitinolMaterial<T>::compute_tangent_modulus(const T& gamma, const T& xi_new)
  {
    // Constants
    const T E = bulk_modulus_ + 2.0 * shear_modulus_;  // Effective modulus (1D version of Eq. 33)
    const T factor = std::sqrt(2.0 / 3.0) - alpha;

    // Compute F_AS
    T sigma_trial = E * (gamma - xi_new * gamma_l_);  // Log-strain formulation
    T F_AS = sigma_trial * factor - C_AS * temperature;

    // Compute F_f
    T R_f = sigma_sAS * (std::sqrt(2.0 / 3.0) + alpha) - C_AS * T_sAS;
    T F_f = F_AS - R_f;

    // Approximate dF_AS/dgamma = d(sigma)/dgamma * factor = E * factor
    T dxi_dgamma = 0.0;
    if (F_AS > 0.0 && F_f < 0.0)
    {
      dxi_dgamma = beta_AS * (1.0 - xi_new) * E * factor / (F_f * F_f);
    }

    // Final consistent tangent modulus
    return E * (1.0 - gamma_l_ * dxi_dgamma);
  }



  template <typename T>
  void BeamNitinolMaterial<T>::evaluate_force_contributions_to_stress(
      Core::LinAlg::Matrix<3, 1, T>& stressN, const Core::LinAlg::Matrix<3, 3, T>& CN,
      const Core::LinAlg::Matrix<3, 1, T>& Gamma, const unsigned int gp)
  {
    delta_xi_s_AS[gp] = compute_xiAS_s(Gamma, xi_old_[gp]);
    delta_xi_s_SA[gp] = compute_xiSA_s(Gamma, xi_old_[gp]);


    T xi_new = xi_old_[gp] + delta_xi_s_AS[gp] - delta_xi_s_SA[gp];
    xi_new = std::clamp(xi_new, 0.0, 1.0);

    // 2. Compute transformation strain
    T gamma_tr = xi_new * gamma_l_;
    T gamma_e = Gamma(0) - gamma_tr;

    T sigma = (bulk_modulus_ + 2.0 * shear_modulus_) * gamma_e;  // Eq. (*)
    sigma_gp_[gp] = sigma;

    sigma_gp_[gp] = compute_tangent_modulus(gamma_e, xi_new);
    stressN(0) = sigma_gp_[gp];
    stressN(1) = 0.0;
    stressN(2) = 0.0;
    // compute material stresses by multiplying strains with constitutive matrix
    // Axial strain from deformation gradient Gamma



    // compute material stresses by multiplying strains with constitutive matrix
    // stressN.multiply(CN, Gamma);
  }


  template <typename T>
  void BeamNitinolMaterial<T>::get_stiffness_matrix_of_moments(
      Core::LinAlg::Matrix<3, 3, T>& stiffM, const Core::LinAlg::Matrix<3, 3, T>& C_M,
      const Core::LinAlg::Matrix<3, 1, T>& Cur, const int gp)
  {
    // Use current C_M as-is (e.g., computed per GP before)
    stiffM = C_M;
  }


  template <typename T>
  void BeamNitinolMaterial<T>::get_stiffness_matrix_of_forces(
      Core::LinAlg::Matrix<3, 3, T>& stiffness_matrix, const Core::LinAlg::Matrix<3, 3, T>& C_N,
      const int gp)
  {
    T A = this->cross_section_area_;
    T G = this->shear_modulus_;
    T kappa = this->shear_correction_factor_;

    stiffness_matrix.clear();
    stiffness_matrix(0, 0) = sigma_gp_[gp];
    stiffness_matrix(1, 1) = G * A * kappa;
    stiffness_matrix(2, 2) = G * A * kappa;

    std::cout << "stiffness matrix:" << std::endl;
    stiffness_matrix.print(std::cout);
  }


  template <typename T>
  void BeamNitinolMaterial<T>::evaluate_moment_contributions_to_stress(
      Core::LinAlg::Matrix<3, 1, T>& stressM, const Core::LinAlg::Matrix<3, 3, T>& CM,
      const Core::LinAlg::Matrix<3, 1, T>& Cur, const unsigned int gp)
  {
    /*
    T kappa_mag = std::sqrt(Cur(1) * Cur(1) + Cur(2) * Cur(2));
    T xi = xi_m_[gp];

    T E_eff = (1.0 - xi) * E_A_ + xi * E_M_;
    T kappa_tr = xi * kappa_L_;
    T kappa_eff = kappa_mag - kappa_tr;
    T M = E_eff * kappa_eff;

    xi = compute_martensite_fraction(M, xi);

    xi_m_[gp] = xi;
    E_eff = (1.0 - xi) * E_A_ + xi * E_M_;
    kappa_tr = xi * kappa_L_;

    // Avoid division by zero: use normalized form robustly
    const T denom = std::max(kappa_mag, 1e-12);
    stressM(0) = Cur(0) * CM(0, 0);  // torsion unchanged
    stressM(1) = E_eff * (Cur(1) - kappa_tr * Cur(1) / denom);
    stressM(2) = E_eff * (Cur(2) - kappa_tr * Cur(2) / denom);

    moment_gp_[gp] = stressM;*/
    // compute material stresses by multiplying curvature with constitutive matrix
    stressM.multiply(CM, Cur);
  }

  template <typename T>
  void Mat::BeamNitinolMaterial<T>::get_constitutive_matrix_of_forces_material_frame(
      Core::LinAlg::Matrix<3, 3, T>& C_N) const
  {
    FOUR_C_THROW("Please do not call get_constitutive_matrix_of_forces_material_frame without gp");
  }

  template <typename T>
  T Mat::BeamNitinolMaterial<T>::compute_F_AS_from_gamma(const T& gamma)
  {
    T elasticModulus = bulk_modulus_ + 2.0 * shear_modulus_;

    T factor = std::sqrt(2.0 / 3.0) - alpha;

    T sigma = elasticModulus * gamma;

    T F_AS = factor * sigma - C_AS * temperature;

    return F_AS;
  }

  template <typename T>
  T Mat::BeamNitinolMaterial<T>::compute_F_AS_uniaxial(const Core::LinAlg::Matrix<3, 1, T>& stressN)
  {
    // Uniaxial stress assumed to be in stressN[0]
    T sigma = stressN(0);

    // Compute pressure: p = -1/3 * trace(stress)
    T p = -(1.0 / 3.0) * (stressN(0) + stressN(1) + stressN(2));

    // Compute ||t|| = sqrt(2/3) * sigma for uniaxial
    T norm_t = sigma * std::sqrt(2.0 / 3.0);

    // Final transformation function
    T F_AS = norm_t + 3.0 * alpha * p - C_AS * temperature;

    return F_AS;
  }

  template <typename T>
  T Mat::BeamNitinolMaterial<T>::compute_xiAS_s(const Core::LinAlg::Matrix<3, 1, T>& gamma, T xi)
  {
    auto F_base = compute_F_AS_from_gamma(gamma(0));

    // Equations (6) and (7): R_f^AS and R_i^AS
    T R_f = sigma_sAS * (std::sqrt(2.0 / 3.0) + alpha) - C_AS * T_sAS;
    T R_i = sigma_fAS * (std::sqrt(2.0 / 3.0) + alpha) - C_AS * T_fAS;

    // Equations (4) and (5): Initial and final transformation functions
    T F_i = F_base - R_i;
    T F_f = F_base - R_f;

    // Equation (8) + (11): Transformation condition check
    // Check transformation condition (Eq. 8)
    bool active = (F_i > 0.0 && F_f < 0.0 && F_base > 0.0);
    if (!active) return T(0.0);  // No contribution

    // Return-map solve (Newton for Δξ)
    T delta_xi = 0.0;
    int iter = 0;
    const T tol = 1e-10;
    const int max_iter = 25;

    while (iter++ < max_iter)
    {
      T res = F_f * F_f * delta_xi - beta_AS * (1.0 - xi) * F_base;
      if (std::abs(res) < tol) break;

      T dres = F_f * F_f;  // simple approx (∂res/∂delta_xi)
      delta_xi -= res / dres;
    }

    return delta_xi;
  }

  template <typename T>
  T Mat::BeamNitinolMaterial<T>::compute_xiSA_s(const Core::LinAlg::Matrix<3, 1, T>& gamma, T xi)
  {
    auto F_base = compute_F_AS_from_gamma(gamma(0));

    // R_f^SA and R_i^SA (Eq. 15 & 16)
    T R_f = sigma_sSA * (std::sqrt(2.0 / 3.0) + alpha) - C_SA * T_sSA;
    T R_i = sigma_fSA * (std::sqrt(2.0 / 3.0) + alpha) - C_SA * T_fSA;

    // Initial and final functions (Eq. 13–14)
    T F_i = F_base - R_i;
    T F_f = F_base - R_f;

    // Activation check (Eq. 17 + 20)
    bool active = (F_i < 0.0 && F_f > 0.0 && F_base < 0.0);
    T H_SA = active ? 1.0 : 0.0;

    // Return-map solve (Newton for Δξ)
    T delta_xi = 0.0;
    int iter = 0;
    const T tol = 1e-10;
    const int max_iter = 25;

    while (iter++ < max_iter)
    {
      T res = F_f * F_f * delta_xi - beta_SA * xi * F_base;
      if (std::abs(res) < tol) break;

      T dres = F_f * F_f;
      delta_xi -= res / dres;
    }

    return delta_xi;
  }

  /*
  template <typename T>
  T Mat::BeamNitinolMaterial<T>::compute_F_SS_uniaxial(const Core::LinAlg::Matrix<3, 1, T>& stressN)
  {
    // Uniaxial stress assumed to be in stressN[0]
    T sigma = stressN(0);

    // Pressure
    T p = -(1.0 / 3.0) * (stressN(0) + stressN(1) + stressN(2));

    // Deviatoric norm ||t|| = sqrt(2/3) * sigma
    T norm_t = sigma * std::sqrt(2.0 / 3.0);

    // f^SS = ||t|| + 3αp - C^SS T (Eq. 21)
    T F_SS = norm_t + 3.0 * alpha * p - C_SS * temperature;

    return F_SS;
  }

  template <typename T>
  bool Mat::BeamNitinolMaterial<T>::isReorientationActive(
      const Core::LinAlg::Matrix<3, 1, T>& stressN)
  {
    T F_SS = compute_F_SS_uniaxial(stressN);
    return (F_SS > 0.0);  // Eq. 24
  }
*/


  template <typename T>
  void Mat::BeamNitinolMaterial<T>::compute_constitutive_parameter(
      Core::LinAlg::Matrix<3, 3, T>& C_N, Core::LinAlg::Matrix<3, 3, T>& C_M)
  {
    for (unsigned int gp = 0; gp < numgp_force_; gp++)
    {
      /*
      // If plasticity for axial strains is enabled, get hardening constitutive parameters
      if (this->params().get_yield_stress_n() >= 0)
      {
        get_hardening_constitutive_matrix_of_forces_material_frame(c_n_eff_[gp]);
        get_effective_yield_stress_n(effyieldstress_n_[gp], this->params().get_yield_stress_n(),
            C_N(0, 0), c_n_eff_[gp](0, 0), gp);

      }*/
    }
    C_M.clear();
    C_M(0, 0) = this->params().get_torsional_rigidity();
    C_M(1, 1) = this->params().get_bending_rigidity2();
    C_M(2, 2) = this->params().get_bending_rigidity3();
    std::cout << "CM" << std::endl;
    C_M.print(std::cout);
  }

  template <typename T>
  void Mat::BeamNitinolMaterial<T>::compute_constitutive_parameter(
      Core::LinAlg::Matrix<3, 3, T>& C_N, Core::LinAlg::Matrix<3, 3, T>& C_M,
      const Core::LinAlg::Matrix<3, 1, T>& Gamma, const int gp)
  {
    // Final stiffness matrix assembly
    /*T A = this->cross_section_area_;
    T G = this->shear_modulus_;
    T kappa = this->shear_correction_factor_;
    C_N.clear();
    C_N(0, 0) = dσ_dε * A;
    C_N(1, 1) = G * A * kappa;
    C_N(2, 2) = G * A * kappa;

    // Moment matrix: keep it linear for now
    C_M.clear();
    C_M(0, 0) = this->params().get_torsional_rigidity();
    C_M(1, 1) = E_eff * this->params().get_mass_moment_of_inertia2();
    C_M(2, 2) = E_eff * this->params().get_mass_moment_of_inertia3();


    std::cout << std::scientific;

    std::cout << "[gp " << gp << "] START" << std::endl;
    C_M.print(std::cout);
    std::cout << "  eps      = " << eps << std::endl;
    std::cout << "  eps_L    = " << eps_L_ << std::endl;
    std::cout << "  E_A_     = " << E_A_ << ", E_M_ = " << E_M_ << std::endl;
    std::cout << "  sigma    = " << sigma << std::endl;*/
  }

}  // namespace Mat

template class Mat::BeamNitinolMaterial<double>;


FOUR_C_NAMESPACE_CLOSE
