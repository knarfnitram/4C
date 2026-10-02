// This file is part of 4C multiphysics licensed under the
// GNU Lesser General Public License v3.0 or later.
//
// See the LICENSE.md file in the top-level for license information.
//
// SPDX-License-Identifier: LGPL-3.0-or-later

#include "4C_mat_beam3r_sma.hpp"

#include "4C_comm_pack_helpers.hpp"
#include "4C_global_data.hpp"
#include "4C_mat_par_bundle.hpp"
#include "4C_utils_exceptions.hpp"
#include "4C_utils_fad.hpp"

#include <algorithm>
#include <cmath>

FOUR_C_NAMESPACE_OPEN

namespace
{
  /// Gauss-Legendre points and weights on [-1, 1]
  void gauss_legendre(unsigned int n, std::vector<double>& points, std::vector<double>& weights)
  {
    points.resize(n);
    weights.resize(n);
    for (unsigned int i = 0; i < n; ++i)
    {
      // Newton iteration for the i-th root of the Legendre polynomial P_n
      double x = std::cos(M_PI * (i + 0.75) / (n + 0.5));
      double dp = 0.0;
      for (unsigned int iter = 0; iter < 100; ++iter)
      {
        double p0 = 1.0;
        double p1 = x;
        for (unsigned int k = 2; k <= n; ++k)
        {
          const double p2 = ((2.0 * k - 1.0) * x * p1 - (k - 1.0) * p0) / k;
          p0 = p1;
          p1 = p2;
        }
        dp = n * (x * p1 - p0) / (x * x - 1.0);
        const double dx = p1 / dp;
        x -= dx;
        if (std::abs(dx) < 1.0e-15) break;
      }
      points[i] = x;
      weights[i] = 2.0 / ((1.0 - x * x) * dp * dp);
    }
  }
}  // namespace

/*-----------------------------------------------------------------------------------------------*/
/* Parameter class                                                                               */
/*-----------------------------------------------------------------------------------------------*/
Mat::PAR::BeamReissnerSMAMaterialParams::BeamReissnerSMAMaterialParams(
    const Core::Mat::PAR::Parameter::Data& matdata)
    : BeamElastHyperMaterialParameterGeneric(matdata),
      temperature_(matdata.parameters.get<double>("TEMP")),
      delta_s_ams_(matdata.parameters.get<double>("DSAMS")),
      t0_sma_(matdata.parameters.get<double>("T0SMA")),
      w_in_(matdata.parameters.get<double>("WIN")),
      r_m_(matdata.parameters.get<double>("RM")),
      rs_f0_(matdata.parameters.get<double>("RSF0")),
      rs_r0_(matdata.parameters.get<double>("RSR0")),
      h_sf_(matdata.parameters.get<double>("HSF")),
      h_sr_(matdata.parameters.get<double>("HSR")),
      c_ts_(matdata.parameters.get<double>("CTS")),
      a_sf0_(matdata.parameters.get<double>("ASF0")),
      a_sf1_(matdata.parameters.get<double>("ASF1")),
      a_sr0_(matdata.parameters.get<double>("ASR0")),
      a_sr1_(matdata.parameters.get<double>("ASR1")),
      n_exp_(matdata.parameters.get<double>("NEXP")),
      eps_l_n_(matdata.parameters.get<double>("EPSLN")),
      kappa_l_m_(matdata.parameters.get<double>("KAPPALM")),
      torsion_sma_(matdata.parameters.get<bool>("TORSIONSMA")),
      moment_model_(matdata.parameters.get<MomentModel>("MOMENTMODEL")),
      num_fiber_rings_(matdata.parameters.get<int>("NUMFIBERRINGS")),
      num_fiber_sectors_(matdata.parameters.get<int>("NUMFIBERSECTORS")),
      rs_regularization_(matdata.parameters.get<double>("RSREG")),
      local_staggered_maxiter_(matdata.parameters.get<int>("LOCALITER")),
      local_substep_max_(matdata.parameters.get<int>("LOCALSUBSTEP")),
      local_substep_strain_limit_(matdata.parameters.get<double>("LOCALSTRAINSTEP")),
      local_substep_curvature_limit_(matdata.parameters.get<double>("LOCALCURVSTEP")),
      youngs_modulus_austenite_(matdata.parameters.get<double>("YOUNG")),
      youngs_modulus_martensite_(matdata.parameters.get<double>("YOUNGMART")),
      shear_modulus_austenite_(Mat::PAR::determine_shear_modulus(matdata)),
      shear_modulus_martensite_(matdata.parameters.get<double>("SHEARMODMART")),
      density_(matdata.parameters.get<double>("DENS")),
      cross_section_area_(matdata.parameters.get<double>("CROSSAREA")),
      shear_correction_factor_(matdata.parameters.get<double>("SHEARCORR")),
      area_moment_inertia_polar_(matdata.parameters.get<double>("MOMINPOL")),
      area_moment_inertia_2_(matdata.parameters.get<double>("MOMIN2")),
      area_moment_inertia_3_(matdata.parameters.get<double>("MOMIN3")),
      radius_interaction_(matdata.parameters.get<double>("INTERACTIONRADIUS")),
      use_consistent_tangent_(matdata.parameters.get<bool>("CONSISTENTTANGENT"))
{
  if (eps_l_n_ < 0.0) FOUR_C_THROW("EPSLN must be non-negative.");
  if (kappa_l_m_ < 0.0) FOUR_C_THROW("KAPPALM must be non-negative.");

  switch (moment_model_)
  {
    case MomentModel::elastic:
      if (eps_l_n_ <= 0.0)
        FOUR_C_THROW(
            "Beam SMA material without transformation: set EPSLN > 0 or choose a MOMENTMODEL "
            "with transformation.");
      break;
    case MomentModel::sectional:
      if (kappa_l_m_ <= 0.0) FOUR_C_THROW("MOMENTMODEL sectional requires KAPPALM > 0.");
      break;
    case MomentModel::fibers:
    {
      if (eps_l_n_ <= 0.0) FOUR_C_THROW("MOMENTMODEL fibers requires EPSLN > 0.");
      if (num_fiber_rings_ < 2) FOUR_C_THROW("NUMFIBERRINGS must be at least 2.");
      if (num_fiber_sectors_ < 4) FOUR_C_THROW("NUMFIBERSECTORS must be at least 4.");

      // the fibre integration assumes a solid circular cross-section of area CROSSAREA
      const double r2 = cross_section_area_ / M_PI;
      const double I_circle = 0.25 * M_PI * r2 * r2;
      auto check = [&](double value, double expected, const char* name)
      {
        if (std::abs(value - expected) > 1.0e-3 * expected)
          FOUR_C_THROW(
              "MOMENTMODEL fibers assumes a solid circular cross-section: {} = {} does not match "
              "the value {} of a circle with CROSSAREA = {}.",
              name, value, expected, cross_section_area_);
      };
      check(area_moment_inertia_2_, I_circle, "MOMIN2");
      check(area_moment_inertia_3_, I_circle, "MOMIN3");
      check(area_moment_inertia_polar_, 2.0 * I_circle, "MOMINPOL");
      break;
    }
  }

  if (kappa_l_m_ > 0.0 && moment_model_ != MomentModel::sectional)
    FOUR_C_THROW("KAPPALM is only used by MOMENTMODEL sectional.");
  if (torsion_sma_ && moment_model_ != MomentModel::sectional)
    FOUR_C_THROW(
        "TORSIONSMA is only used by MOMENTMODEL sectional; the fibre model always includes "
        "torsion.");

  if (n_exp_ <= 0.0 || n_exp_ > 1.0) FOUR_C_THROW("NEXP must satisfy 0 < NEXP <= 1.");

  if (rs_regularization_ <= 0.0) FOUR_C_THROW("RSREG must be positive.");
  if (local_staggered_maxiter_ < 1) FOUR_C_THROW("LOCALITER must be positive.");
  if (local_substep_max_ < 1) FOUR_C_THROW("LOCALSUBSTEP must be positive.");
  if (local_substep_strain_limit_ < 0.0) FOUR_C_THROW("LOCALSTRAINSTEP must be non-negative.");
  if (local_substep_curvature_limit_ < 0.0) FOUR_C_THROW("LOCALCURVSTEP must be non-negative.");

  if (youngs_modulus_austenite_ <= 0.0) FOUR_C_THROW("YOUNG must be positive.");
  if (youngs_modulus_martensite_ <= 0.0) FOUR_C_THROW("YOUNGMART must be positive.");
  if (shear_modulus_austenite_ <= 0.0) FOUR_C_THROW("austenite shear modulus must be positive.");
  if (shear_modulus_martensite_ <= 0.0) FOUR_C_THROW("SHEARMODMART must be positive.");

  if (density_ < 0.0) FOUR_C_THROW("DENS must be non-negative.");
  if (cross_section_area_ <= 0.0) FOUR_C_THROW("CROSSAREA must be positive.");
  if (shear_correction_factor_ <= 0.0) FOUR_C_THROW("SHEARCORR must be positive.");
  if (area_moment_inertia_polar_ < 0.0) FOUR_C_THROW("MOMINPOL must be non-negative.");
  if (area_moment_inertia_2_ < 0.0) FOUR_C_THROW("MOMIN2 must be non-negative.");
  if (area_moment_inertia_3_ < 0.0) FOUR_C_THROW("MOMIN3 must be non-negative.");
}

std::shared_ptr<Core::Mat::Material> Mat::PAR::BeamReissnerSMAMaterialParams::create_material()
{
  std::shared_ptr<Core::Mat::Material> matobject;

  if (uses_fad())
  {
    FOUR_C_THROW("The SMA beam material is implemented for double only; FAD is not supported.");
  }

  matobject = std::make_shared<Mat::BeamSMAMaterial<double>>(this);
  return matobject;
}

/*-----------------------------------------------------------------------------------------------*/
/* ParObject type                                                                                */
/*-----------------------------------------------------------------------------------------------*/
template <typename T>
Core::Communication::ParObject* Mat::BeamSMAMaterialType<T>::create(
    Core::Communication::UnpackBuffer& buffer)
{
  auto* matobject = new Mat::BeamSMAMaterial<T>();
  matobject->unpack(buffer);
  return matobject;
}

template <typename T>
Mat::BeamSMAMaterialType<T> Mat::BeamSMAMaterialType<T>::instance_;

/*-----------------------------------------------------------------------------------------------*/
/* Utility functions                                                                             */
/*-----------------------------------------------------------------------------------------------*/
template <typename T>
double Mat::BeamSMAMaterial<T>::sign(double x)
{
  if (x > 0.0) return 1.0;
  if (x < 0.0) return -1.0;
  return 0.0;
}

template <typename T>
double Mat::BeamSMAMaterial<T>::vec_norm(const Vec3& v)
{
  return std::sqrt(dot(v, v));
}

template <typename T>
double Mat::BeamSMAMaterial<T>::vector_difference_norm(const Vec3& a, const Vec3& b) const
{
  double n2 = 0.0;
  for (unsigned int i = 0; i < 3; ++i)
  {
    const double d = Core::FADUtils::cast_to_double(a(i) - b(i));
    n2 += d * d;
  }
  return std::sqrt(n2);
}

template <typename T>
double Mat::BeamSMAMaterial<T>::dot(const Vec3& a, const Vec3& b)
{
  return Core::FADUtils::cast_to_double(a(0)) * Core::FADUtils::cast_to_double(b(0)) +
         Core::FADUtils::cast_to_double(a(1)) * Core::FADUtils::cast_to_double(b(1)) +
         Core::FADUtils::cast_to_double(a(2)) * Core::FADUtils::cast_to_double(b(2));
}

/*-----------------------------------------------------------------------------------------------*/
/* Material class                                                                                */
/*-----------------------------------------------------------------------------------------------*/
template <typename T>
Mat::BeamSMAMaterial<T>::BeamSMAMaterial(Mat::PAR::BeamReissnerSMAMaterialParams* params)
    : Mat::BeamElastHyperMaterial<T>(params)
{
}

template <typename T>
const Mat::PAR::BeamReissnerSMAMaterialParams& Mat::BeamSMAMaterial<T>::sma_params() const
{
  return static_cast<const Mat::PAR::BeamReissnerSMAMaterialParams&>(this->params());
}

template <typename T>
void Mat::BeamSMAMaterial<T>::setup(int numgp_force, int numgp_moment)
{
  numgp_force_ = static_cast<unsigned int>(numgp_force);
  numgp_moment_ = static_cast<unsigned int>(numgp_moment);

  force_state_conv_.assign(numgp_force_, LocalState{});
  force_state_new_.assign(numgp_force_, LocalState{});
  moment_state_conv_.assign(numgp_moment_, LocalState{});
  moment_state_new_.assign(numgp_moment_, LocalState{});

  c_n_alg_.assign(numgp_force_, Mat3(Core::LinAlg::Initialization::zero));
  c_m_alg_.assign(numgp_moment_, Mat3(Core::LinAlg::Initialization::zero));

  gamma_conv_.assign(numgp_force_, Vec3(Core::LinAlg::Initialization::zero));
  gamma_new_.assign(numgp_force_, Vec3(Core::LinAlg::Initialization::zero));

  cur_conv_.assign(numgp_moment_, Vec3(Core::LinAlg::Initialization::zero));
  cur_new_.assign(numgp_moment_, Vec3(Core::LinAlg::Initialization::zero));

  for (unsigned int gp = 0; gp < numgp_force_; ++gp)
  {
    force_state_conv_[gp].dir.put_scalar(0.0);
    force_state_conv_[gp].dir(0) = 1.0;
    force_state_new_[gp] = force_state_conv_[gp];
  }

  for (unsigned int gp = 0; gp < numgp_moment_; ++gp)
  {
    moment_state_conv_[gp].dir.put_scalar(0.0);
    moment_state_conv_[gp].dir(1) = 1.0;
    moment_state_new_[gp] = moment_state_conv_[gp];
  }

  // fibre states (only known once the material parameters are available; after unpack() they are
  // restored from the buffer)
  fiber_state_conv_.clear();
  if (this->parameter() != nullptr && uses_fibers())
  {
    LocalState initial;
    initial.dir(0) = 1.0;
    fiber_state_conv_.assign(numgp_moment_ * fibers().size(), initial);
  }
  fiber_state_new_ = fiber_state_conv_;
}

/*-----------------------------------------------------------------------------------------------*/
/* Packing and unpacking                                                                         */
/*-----------------------------------------------------------------------------------------------*/
template <typename T>
void Mat::BeamSMAMaterial<T>::pack(Core::Communication::PackBuffer& data) const
{
  int type = unique_par_object_id();
  add_to_pack(data, type);

  int matid = -1;
  if (this->parameter() != nullptr) matid = this->params().id();
  add_to_pack(data, matid);

  add_to_pack(data, numgp_force_);
  add_to_pack(data, numgp_moment_);

  add_to_pack(data, gamma_conv_);
  add_to_pack(data, cur_conv_);

  for (unsigned int gp = 0; gp < numgp_force_; ++gp)
  {
    add_to_pack(data, force_state_conv_[gp].vM);
    add_to_pack(data, force_state_conv_[gp].vS);
    add_to_pack(data, force_state_conv_[gp].dir);
  }

  for (unsigned int gp = 0; gp < numgp_moment_; ++gp)
  {
    add_to_pack(data, moment_state_conv_[gp].vM);
    add_to_pack(data, moment_state_conv_[gp].vS);
    add_to_pack(data, moment_state_conv_[gp].dir);
  }

  add_to_pack(data, fiber_state_conv_.size());
  for (const auto& state : fiber_state_conv_)
  {
    add_to_pack(data, state.vM);
    add_to_pack(data, state.vS);
    add_to_pack(data, state.dir);
  }
}

template <typename T>
void Mat::BeamSMAMaterial<T>::unpack(Core::Communication::UnpackBuffer& buffer)
{
  Core::Communication::extract_and_assert_id(buffer, unique_par_object_id());

  int matid = -1;
  extract_from_pack(buffer, matid);

  extract_from_pack(buffer, numgp_force_);
  extract_from_pack(buffer, numgp_moment_);

  this->setup(static_cast<int>(numgp_force_), static_cast<int>(numgp_moment_));

  // same order as in pack()
  extract_from_pack(buffer, gamma_conv_);
  extract_from_pack(buffer, cur_conv_);
  gamma_new_ = gamma_conv_;
  cur_new_ = cur_conv_;

  for (unsigned int gp = 0; gp < numgp_force_; ++gp)
  {
    extract_from_pack(buffer, force_state_conv_[gp].vM);
    extract_from_pack(buffer, force_state_conv_[gp].vS);
    extract_from_pack(buffer, force_state_conv_[gp].dir);
    force_state_new_[gp] = force_state_conv_[gp];
  }

  for (unsigned int gp = 0; gp < numgp_moment_; ++gp)
  {
    extract_from_pack(buffer, moment_state_conv_[gp].vM);
    extract_from_pack(buffer, moment_state_conv_[gp].vS);
    extract_from_pack(buffer, moment_state_conv_[gp].dir);
    moment_state_new_[gp] = moment_state_conv_[gp];
  }

  std::size_t num_fiber_states = 0;
  extract_from_pack(buffer, num_fiber_states);
  fiber_state_conv_.resize(num_fiber_states);
  for (auto& state : fiber_state_conv_)
  {
    extract_from_pack(buffer, state.vM);
    extract_from_pack(buffer, state.vS);
    extract_from_pack(buffer, state.dir);
  }
  fiber_state_new_ = fiber_state_conv_;

  this->set_parameter(nullptr);

  if (Global::Problem::instance()->materials() != nullptr)
    if (Global::Problem::instance()->materials()->num() != 0)
    {
      const int probinst = Global::Problem::instance()->materials()->get_read_from_problem();
      Core::Mat::PAR::Parameter* mat =
          Global::Problem::instance(probinst)->materials()->parameter_by_id(matid);
      this->set_parameter(static_cast<Mat::PAR::BeamReissnerSMAMaterialParams*>(mat));
    }
}

/*-----------------------------------------------------------------------------------------------*/
/* Constitutive helpers                                                                          */
/*-----------------------------------------------------------------------------------------------*/
template <typename T>
double Mat::BeamSMAMaterial<T>::current_temperature() const
{
  return sma_params().get_temperature();
}

template <typename T>
double Mat::BeamSMAMaterial<T>::transformation_amplitude(Mode mode) const
{
  return mode == Mode::moment ? sma_params().get_kappa_l_m() : sma_params().get_eps_l_n();
}

template <typename T>
double Mat::BeamSMAMaterial<T>::energy_scaling(Mode mode) const
{
  return mode == Mode::fiber ? 1.0 : sma_params().get_sma_cross_section_area();
}

template <typename T>
double Mat::BeamSMAMaterial<T>::transformation_radius(
    Mode mode, unsigned int i, bool forward, double v, double& dR_dv) const
{
  const double a = energy_scaling(mode);

  if (i == idx_m)
  {
    dR_dv = 0.0;
    return a * sma_params().get_r_m();
  }

  const double eps = sma_params().get_rs_regularization();
  const double n = sma_params().get_n_exp();
  const double vp = std::max(v + eps, eps);
  const double omp = std::max(1.0 - v + eps, eps);

  if (forward)
  {
    dR_dv = a * (sma_params().get_h_sf() + n * sma_params().get_a_sf0() * std::pow(vp, n - 1.0) -
                    n * sma_params().get_a_sf1() * std::pow(omp, n - 1.0));
    return a * (sma_params().get_rs_f0() + sma_params().get_h_sf() * v +
                   sma_params().get_a_sf0() * std::pow(vp, n) +
                   sma_params().get_a_sf1() * std::pow(omp, n));
  }

  const double dT = current_temperature() - sma_params().get_t0_sma();
  dR_dv = a * (-sma_params().get_h_sr() + n * sma_params().get_a_sr0() * std::pow(vp, n - 1.0) -
                  n * sma_params().get_a_sr1() * std::pow(omp, n - 1.0));
  return a * (sma_params().get_rs_r0() - sma_params().get_h_sr() * v -
                 sma_params().get_c_ts() * dT + sma_params().get_a_sr0() * std::pow(vp, n) +
                 sma_params().get_a_sr1() * std::pow(omp, n));
}

template <typename T>
void Mat::BeamSMAMaterial<T>::compute_phase_constitutive_matrices(
    double m, Mode mode, Mat3& C, Mat3& dC, Mat3& d2C) const
{
  C.clear();
  dC.clear();
  d2C.clear();

  // Reuss-type mixture of the compliances: 1/E(m) = 1/E_A + m (1/E_M - 1/E_A), m = vM + vS
  const double bE = 1.0 / sma_params().get_youngs_modulus_martensite() -
                    1.0 / sma_params().get_youngs_modulus_austenite();
  const double bG = 1.0 / sma_params().get_shear_modulus_martensite() -
                    1.0 / sma_params().get_shear_modulus_austenite();

  const double E =
      1.0 / std::max(1.0 / sma_params().get_youngs_modulus_austenite() + m * bE, 1.0e-30);
  const double G =
      1.0 / std::max(1.0 / sma_params().get_shear_modulus_austenite() + m * bG, 1.0e-30);

  const double dE = -E * E * bE;
  const double dG = -G * G * bG;
  const double d2E = 2.0 * E * E * E * bE * bE;
  const double d2G = 2.0 * G * G * G * bG * bG;

  // section factors of the three components
  double sE[3] = {0.0, 0.0, 0.0};
  double sG[3] = {0.0, 0.0, 0.0};
  if (mode == Mode::force)
  {
    const double A = sma_params().get_sma_cross_section_area();
    const double ks = sma_params().get_sma_shear_correction_factor();
    sE[0] = A;
    sG[1] = A * ks;
    sG[2] = A * ks;
  }
  else if (mode == Mode::fiber)
  {
    // stresses (sigma, tau_12, tau_13) of a fibre
    sE[0] = 1.0;
    sG[1] = 1.0;
    sG[2] = 1.0;
  }
  else
  {
    sG[0] = sma_params().get_sma_moment_inertia_polar();
    sE[1] = sma_params().get_sma_moment_inertia2();
    sE[2] = sma_params().get_sma_moment_inertia3();
  }

  for (unsigned int k = 0; k < 3; ++k)
  {
    C(k, k) = sE[k] * E + sG[k] * G;
    dC(k, k) = sE[k] * dE + sG[k] * dG;
    d2C(k, k) = sE[k] * d2E + sG[k] * d2G;
  }
}

template <typename T>
void Mat::BeamSMAMaterial<T>::evaluate_driving_forces(
    Mode mode, const Vec3& q, double vM, double vS, const Vec3& dir, DrivingForces& df) const
{
  Mat3 C, dC, d2C;
  compute_phase_constitutive_matrices(vM + vS, mode, C, dC, d2C);

  const double kappa = transformation_amplitude(mode);

  // elastic part of the strain resultant
  Vec3 q_el;
  for (unsigned int k = 0; k < 3; ++k) q_el(k) = q(k) - kappa * vS * dir(k);

  Vec3 Cd, P, d2C_q_el;
  df.stress.multiply(C, q_el);
  Cd.multiply(C, dir);
  P.multiply(dC, q_el);
  d2C_q_el.multiply(d2C, q_el);

  const double half_q_P = 0.5 * dot(q_el, P);
  const double Q = 0.5 * dot(q_el, d2C_q_el);
  const double dir_P = dot(dir, P);
  const double dir_Cd = dot(dir, Cd);

  // chemical and interaction part: theta = dS (T - T0) + W (1 - 2 (vM + vS)), scaled from energy
  // per volume to the energy measure of the mode
  const double a = energy_scaling(mode);
  const double w_in = a * sma_params().get_w_in();
  const double theta =
      a * sma_params().get_delta_s_ams() * (current_temperature() - sma_params().get_t0_sma()) +
      w_in * (1.0 - 2.0 * (vM + vS));

  // thermodynamic driving forces B = -dpsi/dv
  df.B[idx_s] = kappa * dot(dir, df.stress) - half_q_P - theta;
  df.B[idx_m] = -half_q_P - theta;

  df.dB_dv[idx_s][idx_s] = 2.0 * kappa * dir_P - kappa * kappa * dir_Cd - Q + 2.0 * w_in;
  df.dB_dv[idx_s][idx_m] = kappa * dir_P - Q + 2.0 * w_in;
  df.dB_dv[idx_m][idx_s] = df.dB_dv[idx_s][idx_m];
  df.dB_dv[idx_m][idx_m] = -Q + 2.0 * w_in;

  for (unsigned int k = 0; k < 3; ++k)
  {
    df.dB_dq[idx_s](k) = kappa * Cd(k) - P(k);
    df.dB_dq[idx_m](k) = -P(k);
    df.dstress_dv[idx_s](k) = P(k) - kappa * Cd(k);
    df.dstress_dv[idx_m](k) = P(k);
  }
}

template <typename T>
double Mat::BeamSMAMaterial<T>::consistency_function(Mode mode, const Vec3& q, const Vec3& dir,
    unsigned int i, double v_other, bool forward, double v, double& dg_dv) const
{
  const double vS = (i == idx_s) ? v : v_other;
  const double vM = (i == idx_s) ? v_other : v;

  DrivingForces df;
  evaluate_driving_forces(mode, q, vM, vS, dir, df);

  double dR_dv = 0.0;
  const double R = transformation_radius(mode, i, forward, v, dR_dv);
  const double s = forward ? 1.0 : -1.0;

  dg_dv = df.dB_dv[i][i] - s * dR_dv;
  return df.B[i] - s * R;
}

template <typename T>
double Mat::BeamSMAMaterial<T>::find_consistent_phase_fraction(Mode mode, const Vec3& q,
    const Vec3& dir, unsigned int i, double v_other, bool forward, double v_pos, double v_neg) const
{
  // The consistency function g is positive at v_pos and negative at v_neg. Newton's method is
  // safeguarded by bisection on this bracket, so it always converges.
  double x = forward ? v_pos : v_neg;

  for (unsigned int iter = 0; iter < max_root_iterations_; ++iter)
  {
    double dg = 0.0;
    const double g = consistency_function(mode, q, dir, i, v_other, forward, x, dg);

    if (g > 0.0)
      v_pos = x;
    else if (g < 0.0)
      v_neg = x;
    else
      return x;

    const double lower = std::min(v_pos, v_neg);
    const double upper = std::max(v_pos, v_neg);

    double x_new = (dg != 0.0) ? x - g / dg : 0.5 * (lower + upper);
    if (!(x_new > lower && x_new < upper)) x_new = 0.5 * (lower + upper);

    if (std::abs(x_new - x) <= phase_fraction_tolerance_ ||
        upper - lower <= phase_fraction_tolerance_)
      return x_new;

    x = x_new;
  }

  FOUR_C_THROW("Beam SMA: phase-fraction update did not converge within {} iterations.",
      max_root_iterations_);
}

template <typename T>
typename Mat::BeamSMAMaterial<T>::Transformation Mat::BeamSMAMaterial<T>::update_phase_fraction(
    Mode mode, const Vec3& q, const Vec3& dir, unsigned int i, double v_other, double v_old,
    double& v) const
{
  // admissible interval for this phase fraction: 0 <= v <= 1 - v_other
  const double upper = std::max(0.0, 1.0 - v_other);
  const double v_start = std::clamp(v_old, 0.0, upper);

  double dummy = 0.0;

  // forward transformation: B - R_f > 0 at the start of the step
  if (consistency_function(mode, q, dir, i, v_other, true, v_start, dummy) > 0.0)
  {
    if (v_start >= upper ||
        consistency_function(mode, q, dir, i, v_other, true, upper, dummy) >= 0.0)
    {
      v = upper;
      return Transformation::upper_bound;
    }
    v = find_consistent_phase_fraction(mode, q, dir, i, v_other, true, v_start, upper);
    return Transformation::forward;
  }

  // reverse transformation: B + R_r < 0 at the start of the step
  if (consistency_function(mode, q, dir, i, v_other, false, v_start, dummy) < 0.0)
  {
    if (v_start <= 0.0 || consistency_function(mode, q, dir, i, v_other, false, 0.0, dummy) <= 0.0)
    {
      v = 0.0;
      return Transformation::lower_bound;
    }
    v = find_consistent_phase_fraction(mode, q, dir, i, v_other, false, 0.0, v_start);
    return Transformation::reverse;
  }

  // no transformation; v_start < v_old only if the other phase fraction occupies the space
  v = v_start;
  return (v_start == v_old) ? Transformation::none : Transformation::upper_bound;
}

template <typename T>
void Mat::BeamSMAMaterial<T>::solve_local_state(Mode mode, const Vec3& q,
    const LocalState& old_state, const Vec3& dir, LocalState& new_state,
    Transformation status[2]) const
{
  const double vM_old = Core::FADUtils::cast_to_double(old_state.vM);
  const double vS_old = Core::FADUtils::cast_to_double(old_state.vS);

  double vM = vM_old;
  double vS = vS_old;

  // Staggered solution of the two (weakly coupled) phase-fraction updates. Each scalar update is a
  // return mapping that includes the bounds 0 <= vM, vS and vM + vS <= 1, so no Lagrange
  // multipliers are needed.
  const auto max_iter = static_cast<unsigned int>(sma_params().get_local_staggered_maxiter());
  bool converged = false;
  for (unsigned int iter = 0; iter < std::max(max_iter, 2u); ++iter)
  {
    const double vS_prev = vS;
    const double vM_prev = vM;

    status[idx_s] = update_phase_fraction(mode, q, dir, idx_s, vM, vS_old, vS);
    status[idx_m] = update_phase_fraction(mode, q, dir, idx_m, vS, vM_old, vM);

    if (std::abs(vS - vS_prev) + std::abs(vM - vM_prev) <= phase_fraction_tolerance_)
    {
      converged = true;
      break;
    }
  }

  if (!converged)
    FOUR_C_THROW("Beam SMA: staggered update of vM and vS did not converge within {} iterations.",
        std::max(max_iter, 2u));

  new_state.vM = vM;
  new_state.vS = vS;
  new_state.dir = dir;
}

template <typename T>
void Mat::BeamSMAMaterial<T>::choose_force_direction(
    const Vec3& Gamma, const LocalState& old_state, Vec3& dir) const
{
  dir.put_scalar(0.0);

  // The direction of the transformation strain is part of the state. As long as single-variant
  // martensite exists, keep it; otherwise the stored transformation strain would flip sign as soon
  // as the trial force changes sign (e.g. during reverse transformation).
  if (Core::FADUtils::cast_to_double(old_state.vS) > direction_reset_tolerance_ &&
      std::abs(Core::FADUtils::cast_to_double(old_state.dir(0))) > 0.5)
  {
    dir(0) = sign(Core::FADUtils::cast_to_double(old_state.dir(0)));
    return;
  }

  Mat3 C_N, dC_N, d2C_N;
  compute_phase_constitutive_matrices(
      Core::FADUtils::cast_to_double(old_state.vM + old_state.vS), Mode::force, C_N, dC_N, d2C_N);

  const double axial_trial =
      Core::FADUtils::cast_to_double(C_N(0, 0)) *
      (Core::FADUtils::cast_to_double(Gamma(0)) -
          sma_params().get_eps_l_n() * Core::FADUtils::cast_to_double(old_state.vS) *
              Core::FADUtils::cast_to_double(old_state.dir(0)));

  if (std::abs(axial_trial) > 1.0e-14)
    dir(0) = sign(axial_trial);
  else if (std::abs(Core::FADUtils::cast_to_double(old_state.dir(0))) > 0.5)
    dir(0) = old_state.dir(0);
  else
    dir(0) = 1.0;
}

template <typename T>
void Mat::BeamSMAMaterial<T>::choose_moment_direction(
    const Vec3& Cur, const LocalState& old_state, Vec3& dir) const
{
  dir.put_scalar(0.0);

  // Keep the direction of the transformation curvature as long as single-variant martensite
  // exists. Otherwise it would flip on unloading or rotate with the trial moment without any
  // dissipation. A new direction is only chosen from (practically) pure austenite.
  if (Core::FADUtils::cast_to_double(old_state.vS) > direction_reset_tolerance_)
  {
    dir = old_state.dir;
    if (!sma_params().get_torsion_sma()) dir(0) = 0.0;

    const double old_nrm = vec_norm(dir);
    if (old_nrm > 1.0e-14)
    {
      for (unsigned int i = 0; i < 3; ++i) dir(i) /= old_nrm;
      return;
    }
    dir.put_scalar(0.0);
  }

  Mat3 C_M, dC_M, d2C_M;
  compute_phase_constitutive_matrices(
      Core::FADUtils::cast_to_double(old_state.vM + old_state.vS), Mode::moment, C_M, dC_M, d2C_M);

  const double kappa_l = sma_params().get_kappa_l_m();
  Vec3 kappa_el;
  for (unsigned int i = 0; i < 3; ++i)
    kappa_el(i) = Cur(i) - kappa_l * old_state.vS * old_state.dir(i);

  Vec3 Mtrial;
  Mtrial.multiply(C_M, kappa_el);

  if (!sma_params().get_torsion_sma()) Mtrial(0) = 0.0;

  const double nrm = vec_norm(Mtrial);
  if (nrm > 1.0e-14)
  {
    for (unsigned int i = 0; i < 3; ++i) dir(i) = Mtrial(i) / nrm;
  }
  else
  {
    dir = old_state.dir;
    if (!sma_params().get_torsion_sma()) dir(0) = 0.0;

    const double old_nrm = vec_norm(dir);
    if (old_nrm > 1.0e-14)
    {
      for (unsigned int i = 0; i < 3; ++i) dir(i) /= old_nrm;
    }
    else
    {
      dir.put_scalar(0.0);
      dir(1) = 1.0;
    }
  }
}

template <typename T>
void Mat::BeamSMAMaterial<T>::choose_fiber_direction(
    const Vec3& e, const LocalState& old_state, Vec3& dir) const
{
  // The direction lives in the space of fibre strains (epsilon, gamma_12, gamma_13) with the norm
  // |d|^2 = d_0^2 + (d_1^2 + d_2^2) / 3. Then the uniaxial transformation strain EPSLN corresponds
  // to the transformation shear strain sqrt(3) EPSLN, and the work d . sigma of a direction aligned
  // with the stress is the von Mises stress sqrt(sigma^2 + 3 tau^2).
  if (Core::FADUtils::cast_to_double(old_state.vS) > direction_reset_tolerance_)
  {
    dir = old_state.dir;
    return;
  }

  Mat3 C, dC, d2C;
  compute_phase_constitutive_matrices(
      Core::FADUtils::cast_to_double(old_state.vM + old_state.vS), Mode::fiber, C, dC, d2C);

  Vec3 e_el;
  for (unsigned int k = 0; k < 3; ++k)
    e_el(k) = e(k) - sma_params().get_eps_l_n() * old_state.vS * old_state.dir(k);

  Vec3 stress_trial;
  stress_trial.multiply(C, e_el);

  const double sigma_eq =
      std::sqrt(stress_trial(0) * stress_trial(0) + 3.0 * stress_trial(1) * stress_trial(1) +
                3.0 * stress_trial(2) * stress_trial(2));
  if (sigma_eq > 1.0e-14)
  {
    dir(0) = stress_trial(0) / sigma_eq;
    dir(1) = 3.0 * stress_trial(1) / sigma_eq;
    dir(2) = 3.0 * stress_trial(2) / sigma_eq;
  }
  else
  {
    dir = old_state.dir;
  }
}

template <typename T>
void Mat::BeamSMAMaterial<T>::choose_direction(
    Mode mode, const Vec3& q, const LocalState& old_state, Vec3& dir) const
{
  switch (mode)
  {
    case Mode::force:
      choose_force_direction(q, old_state, dir);
      break;
    case Mode::moment:
      choose_moment_direction(q, old_state, dir);
      break;
    case Mode::fiber:
      choose_fiber_direction(q, old_state, dir);
      break;
  }
}

/*-----------------------------------------------------------------------------------------------*/
/* Local updates                                                                                 */
/*-----------------------------------------------------------------------------------------------*/
template <typename T>
void Mat::BeamSMAMaterial<T>::integrate_local_state(Mode mode, const Vec3& q_start,
    const Vec3& q_end, const LocalState& state_start, LocalState& state_end, Vec3& stress,
    Mat3& C_alg) const
{
  const double limit = (mode == Mode::moment) ? sma_params().get_local_substep_curvature_limit()
                                              : sma_params().get_local_substep_strain_limit();
  unsigned int nsub = 1;
  if (limit > 0.0)
  {
    const auto n = static_cast<int>(std::ceil(vector_difference_norm(q_end, q_start) / limit));
    nsub =
        static_cast<unsigned int>(std::max(1, std::min(n, sma_params().get_local_substep_max())));
  }

  // D[i] = d v_i / d q_end, carried through all substeps (exact chain rule)
  Vec3 D[2];

  LocalState state = state_start;
  DrivingForces df;

  for (unsigned int isub = 1; isub <= nsub; ++isub)
  {
    const double alpha = static_cast<double>(isub) / static_cast<double>(nsub);

    Vec3 q;
    for (unsigned int k = 0; k < 3; ++k) q(k) = q_start(k) + alpha * (q_end(k) - q_start(k));

    Vec3 dir;
    choose_direction(mode, q, state, dir);

    LocalState state_new;
    Transformation status[2];
    solve_local_state(mode, q, state, dir, state_new, status);

    const double v_new[2] = {
        Core::FADUtils::cast_to_double(state_new.vS), Core::FADUtils::cast_to_double(state_new.vM)};
    evaluate_driving_forces(mode, q, v_new[idx_m], v_new[idx_s], dir, df);

    // Linearization of the local update. Phase fractions that do not change keep their
    // sensitivity, those fixed at a bound have none, and the transforming ones follow from the
    // consistency conditions r_i = B_i(v, q) - s_i R_i(v_i) = 0.
    Vec3 D_new[2];
    bool active[2] = {false, false};
    double s[2] = {0.0, 0.0};
    for (unsigned int i = 0; i < 2; ++i)
    {
      if (status[i] == Transformation::none)
        D_new[i] = D[i];
      else if (status[i] == Transformation::forward || status[i] == Transformation::reverse)
      {
        active[i] = true;
        s[i] = (status[i] == Transformation::forward) ? 1.0 : -1.0;
      }
    }

    if (active[idx_s] || active[idx_m])
    {
      // K = dr_A/dv_A, rhs = dr_A/dq_end including the sensitivity of the inactive fraction
      double K[2][2] = {{0.0, 0.0}, {0.0, 0.0}};
      Vec3 rhs[2];
      for (unsigned int i = 0; i < 2; ++i)
      {
        if (!active[i]) continue;
        double dR_dv = 0.0;
        transformation_radius(mode, i, s[i] > 0.0, v_new[i], dR_dv);
        for (unsigned int j = 0; j < 2; ++j) K[i][j] = df.dB_dv[i][j];
        K[i][i] -= s[i] * dR_dv;

        for (unsigned int k = 0; k < 3; ++k)
        {
          rhs[i](k) = alpha * df.dB_dq[i](k);
          for (unsigned int j = 0; j < 2; ++j)
            if (!active[j]) rhs[i](k) += df.dB_dv[i][j] * D_new[j](k);
        }
      }

      if (active[idx_s] && active[idx_m])
      {
        const double det = K[0][0] * K[1][1] - K[0][1] * K[1][0];
        if (std::abs(det) > 1.0e-14 * (std::abs(K[0][0] * K[1][1]) + std::abs(K[0][1] * K[1][0])))
        {
          for (unsigned int k = 0; k < 3; ++k)
          {
            D_new[0](k) = -(K[1][1] * rhs[0](k) - K[0][1] * rhs[1](k)) / det;
            D_new[1](k) = -(-K[1][0] * rhs[0](k) + K[0][0] * rhs[1](k)) / det;
          }
        }
      }
      else
      {
        const unsigned int i = active[idx_s] ? idx_s : idx_m;
        if (K[i][i] != 0.0)
          for (unsigned int k = 0; k < 3; ++k) D_new[i](k) = -rhs[i](k) / K[i][i];
      }
    }

    D[0] = D_new[0];
    D[1] = D_new[1];
    state = state_new;
  }

  state_end = state;
  stress = df.stress;

  Mat3 C, dC, d2C;
  compute_phase_constitutive_matrices(
      Core::FADUtils::cast_to_double(state_end.vM + state_end.vS), mode, C, dC, d2C);

  C_alg = C;
  if (!sma_params().use_consistent_tangent()) return;

  // consistent tangent: dstress/dq = C + sum_i dstress/dv_i (x) dv_i/dq
  for (unsigned int i = 0; i < 2; ++i)
    for (unsigned int k = 0; k < 3; ++k)
      for (unsigned int l = 0; l < 3; ++l) C_alg(k, l) += df.dstress_dv[i](k) * D[i](l);
}

template <typename T>
void Mat::BeamSMAMaterial<T>::local_response(
    Mode mode, const Vec3& q, const unsigned int gp, Vec3& stress, Mat3& C_alg)
{
  const bool force = (mode == Mode::force);
  const LocalState& state_conv = force ? force_state_conv_[gp] : moment_state_conv_[gp];
  const Vec3& q_conv = force ? gamma_conv_[gp] : cur_conv_[gp];
  LocalState& state_new = force ? force_state_new_[gp] : moment_state_new_[gp];
  Vec3& q_new = force ? gamma_new_[gp] : cur_new_[gp];

  q_new = q;

  // without a transformation amplitude for this mode, the response is elastic with the
  // stiffness of the current phase composition
  if (transformation_amplitude(mode) <= 0.0)
  {
    state_new = state_conv;
    Mat3 dC, d2C;
    compute_phase_constitutive_matrices(
        Core::FADUtils::cast_to_double(state_conv.vM + state_conv.vS), mode, C_alg, dC, d2C);
    stress.multiply(C_alg, q);
    return;
  }

  integrate_local_state(mode, q_conv, q, state_conv, state_new, stress, C_alg);
}

template <typename T>
bool Mat::BeamSMAMaterial<T>::uses_fibers() const
{
  return sma_params().get_moment_model() ==
         Mat::PAR::BeamReissnerSMAMaterialParams::MomentModel::fibers;
}

template <typename T>
const std::vector<typename Mat::BeamSMAMaterial<T>::Fiber>& Mat::BeamSMAMaterial<T>::fibers() const
{
  if (!fibers_.empty()) return fibers_;

  // Gauss-Legendre rule in the radius (exact for the polynomial moments of area, inertia and polar
  // inertia) and equidistant sectors in the angle (exact for the trigonometric terms)
  const double radius = std::sqrt(sma_params().get_sma_cross_section_area() / M_PI);
  const auto n_rings = static_cast<unsigned int>(sma_params().get_num_fiber_rings());
  const auto n_sectors = static_cast<unsigned int>(sma_params().get_num_fiber_sectors());

  std::vector<double> xi, w;
  gauss_legendre(n_rings, xi, w);

  fibers_.reserve(n_rings * n_sectors);
  for (unsigned int i = 0; i < n_rings; ++i)
  {
    const double rho = 0.5 * radius * (1.0 + xi[i]);
    const double ring_weight = 0.5 * radius * w[i] * rho * 2.0 * M_PI / n_sectors;
    for (unsigned int j = 0; j < n_sectors; ++j)
    {
      const double phi = 2.0 * M_PI * (j + 0.5) / n_sectors;
      fibers_.push_back({rho * std::cos(phi), rho * std::sin(phi), ring_weight});
    }
  }
  return fibers_;
}

template <typename T>
void Mat::BeamSMAMaterial<T>::fiber_strains(const Fiber& fiber, const Vec3& Cur, Vec3& e)
{
  // axial strain from bending, engineering shear strains from (warping-free) torsion
  e(0) = fiber.z * Cur(1) - fiber.y * Cur(2);
  e(1) = -fiber.z * Cur(0);
  e(2) = fiber.y * Cur(0);
}

template <typename T>
void Mat::BeamSMAMaterial<T>::local_fiber_moment_response(
    const Vec3& Cur, const unsigned int gp, Vec3& stress, Mat3& C_alg)
{
  const std::vector<Fiber>& fiber_list = fibers();
  const std::size_t n_fibers = fiber_list.size();

  if (fiber_state_conv_.size() != numgp_moment_ * n_fibers)
    FOUR_C_THROW("Beam SMA: fibre states are not set up for {} fibres.", n_fibers);

  const Vec3& Cur_conv = cur_conv_[gp];
  cur_new_[gp] = Cur;

  stress.clear();
  C_alg.clear();

  for (std::size_t f = 0; f < n_fibers; ++f)
  {
    const Fiber& fiber = fiber_list[f];

    Vec3 e_conv, e;
    fiber_strains(fiber, Cur_conv, e_conv);
    fiber_strains(fiber, Cur, e);

    const std::size_t index = gp * n_fibers + f;
    Vec3 fiber_stress;
    Mat3 fiber_tangent;
    integrate_local_state(Mode::fiber, e_conv, e, fiber_state_conv_[index], fiber_state_new_[index],
        fiber_stress, fiber_tangent);

    // M = sum_f w B^T s and C_M = sum_f w B^T C_f B with e = B K
    const double B[3][3] = {{0.0, fiber.z, -fiber.y}, {-fiber.z, 0.0, 0.0}, {fiber.y, 0.0, 0.0}};
    double CB[3][3];
    for (unsigned int k = 0; k < 3; ++k)
      for (unsigned int l = 0; l < 3; ++l)
      {
        CB[k][l] = 0.0;
        for (unsigned int m = 0; m < 3; ++m) CB[k][l] += fiber_tangent(k, m) * B[m][l];
      }

    for (unsigned int k = 0; k < 3; ++k)
    {
      for (unsigned int m = 0; m < 3; ++m) stress(k) += fiber.weight * B[m][k] * fiber_stress(m);
      for (unsigned int l = 0; l < 3; ++l)
        for (unsigned int m = 0; m < 3; ++m) C_alg(k, l) += fiber.weight * B[m][k] * CB[m][l];
    }
  }
}

/*-----------------------------------------------------------------------------------------------*/
/* Interface methods                                                                             */
/*-----------------------------------------------------------------------------------------------*/
template <typename T>
void Mat::BeamSMAMaterial<T>::get_constitutive_matrix_of_forces_material_frame(
    Core::LinAlg::Matrix<3, 3, T>& C_N) const
{
  Mat3 dC, d2C;
  compute_phase_constitutive_matrices(0.0, Mode::force, C_N, dC, d2C);
}

template <typename T>
void Mat::BeamSMAMaterial<T>::get_constitutive_matrix_of_moments_material_frame(
    Core::LinAlg::Matrix<3, 3, T>& C_M) const
{
  Mat3 dC, d2C;
  compute_phase_constitutive_matrices(0.0, Mode::moment, C_M, dC, d2C);
}

template <typename T>
void Mat::BeamSMAMaterial<T>::compute_constitutive_parameter(
    Core::LinAlg::Matrix<3, 3, T>& C_N, Core::LinAlg::Matrix<3, 3, T>& C_M)
{
  get_constitutive_matrix_of_forces_material_frame(C_N);
  get_constitutive_matrix_of_moments_material_frame(C_M);

  for (unsigned int gp = 0; gp < numgp_force_; ++gp) c_n_alg_[gp] = C_N;
  for (unsigned int gp = 0; gp < numgp_moment_; ++gp) c_m_alg_[gp] = C_M;
}

template <typename T>
void Mat::BeamSMAMaterial<T>::evaluate_force_contributions_to_stress(
    Core::LinAlg::Matrix<3, 1, T>& stressN, const Core::LinAlg::Matrix<3, 3, T>&,
    const Core::LinAlg::Matrix<3, 1, T>& Gamma, const unsigned int gp)
{
  local_response(Mode::force, Gamma, gp, stressN, c_n_alg_[gp]);
}

template <typename T>
void Mat::BeamSMAMaterial<T>::evaluate_moment_contributions_to_stress(
    Core::LinAlg::Matrix<3, 1, T>& stressM, const Core::LinAlg::Matrix<3, 3, T>&,
    const Core::LinAlg::Matrix<3, 1, T>& Cur, const unsigned int gp)
{
  if (uses_fibers())
    local_fiber_moment_response(Cur, gp, stressM, c_m_alg_[gp]);
  else
    local_response(Mode::moment, Cur, gp, stressM, c_m_alg_[gp]);
}

template <typename T>
void Mat::BeamSMAMaterial<T>::get_stiffness_matrix_of_forces(
    Core::LinAlg::Matrix<3, 3, T>& stiffness_matrix, const Core::LinAlg::Matrix<3, 3, T>&,
    const int gp)
{
  stiffness_matrix = c_n_alg_[gp];
}

template <typename T>
void Mat::BeamSMAMaterial<T>::get_stiffness_matrix_of_moments(
    Core::LinAlg::Matrix<3, 3, T>& stiffness_matrix, const Core::LinAlg::Matrix<3, 3, T>&,
    const int gp)
{
  stiffness_matrix = c_m_alg_[gp];
}

template <typename T>
double Mat::BeamSMAMaterial<T>::get_translational_mass_inertia_factor() const
{
  return sma_params().get_sma_density() * sma_params().get_sma_cross_section_area();
}

template <typename T>
void Mat::BeamSMAMaterial<T>::get_mass_moment_of_inertia_tensor_material_frame(
    Core::LinAlg::Matrix<3, 3>& J) const
{
  J.clear();
  J(0, 0) = sma_params().get_sma_density() *
            (sma_params().get_sma_moment_inertia2() + sma_params().get_sma_moment_inertia3());
  J(1, 1) = sma_params().get_sma_density() * sma_params().get_sma_moment_inertia2();
  J(2, 2) = sma_params().get_sma_density() * sma_params().get_sma_moment_inertia3();
}

template <typename T>
void Mat::BeamSMAMaterial<T>::get_mass_moment_of_inertia_tensor_material_frame(
    Core::LinAlg::Matrix<3, 3, Sacado::Fad::DFad<double>>& J) const
{
  J.clear();
  J(0, 0) = sma_params().get_sma_density() *
            (sma_params().get_sma_moment_inertia2() + sma_params().get_sma_moment_inertia3());
  J(1, 1) = sma_params().get_sma_density() * sma_params().get_sma_moment_inertia2();
  J(2, 2) = sma_params().get_sma_density() * sma_params().get_sma_moment_inertia3();
}

template <typename T>
double Mat::BeamSMAMaterial<T>::get_interaction_radius() const
{
  return sma_params().get_interaction_radius();
}

template <typename T>
void Mat::BeamSMAMaterial<T>::update()
{
  force_state_conv_ = force_state_new_;
  moment_state_conv_ = moment_state_new_;
  fiber_state_conv_ = fiber_state_new_;

  gamma_conv_ = gamma_new_;
  cur_conv_ = cur_new_;
}

template <typename T>
void Mat::BeamSMAMaterial<T>::reset()
{
  force_state_new_ = force_state_conv_;
  moment_state_new_ = moment_state_conv_;
  fiber_state_new_ = fiber_state_conv_;

  gamma_new_ = gamma_conv_;
  cur_new_ = cur_conv_;
}

/*-----------------------------------------------------------------------------------------------*/
/* Explicit instantiations                                                                       */
/*-----------------------------------------------------------------------------------------------*/
template class Mat::BeamSMAMaterial<double>;
template class Mat::BeamSMAMaterialType<double>;

FOUR_C_NAMESPACE_CLOSE
