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
#include "4C_utils_local_newton.hpp"

#include <algorithm>
#include <cmath>
#include <limits>
#include <tuple>
#include <utility>

FOUR_C_NAMESPACE_OPEN

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
      fb_regularization_(matdata.parameters.get<double>("FBREG")),
      rs_regularization_(matdata.parameters.get<double>("RSREG")),
      local_newton_tol_(matdata.parameters.get<double>("LOCALTOL")),
      local_newton_maxiter_(matdata.parameters.get<int>("LOCALITER")),
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
  if (eps_l_n_ <= 0.0 && kappa_l_m_ <= 0.0)
  {
    FOUR_C_THROW(
        "Beam SMA material requires at least one positive transformation amplitude: EPSLN or "
        "KAPPALM.");
  }

  if (n_exp_ <= 0.0 || n_exp_ > 1.0) FOUR_C_THROW("NEXP must satisfy 0 < NEXP <= 1.");

  if (fb_regularization_ <= 0.0) FOUR_C_THROW("FBREG must be positive.");
  if (rs_regularization_ <= 0.0) FOUR_C_THROW("RSREG must be positive.");
  if (local_newton_tol_ <= 0.0) FOUR_C_THROW("LOCALTOL must be positive.");
  if (local_newton_maxiter_ < 1) FOUR_C_THROW("LOCALITER must be positive.");
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
    FOUR_C_THROW(
        "The reduced SMA beam material prototype is implemented for double only; FAD is not "
        "supported yet.");
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
double Mat::BeamSMAMaterial<T>::abs_reg(double a, double eps)
{
  return std::sqrt(a * a + eps * eps);
}

template <typename T>
double Mat::BeamSMAMaterial<T>::fb(double a, double b, double eps)
{
  return std::sqrt(a * a + b * b + 2.0 * eps * eps) + a - b;
}

template <typename T>
double Mat::BeamSMAMaterial<T>::vec_norm(const Vec3& v)
{
  return std::sqrt(Core::FADUtils::cast_to_double(v(0)) * Core::FADUtils::cast_to_double(v(0)) +
                   Core::FADUtils::cast_to_double(v(1)) * Core::FADUtils::cast_to_double(v(1)) +
                   Core::FADUtils::cast_to_double(v(2)) * Core::FADUtils::cast_to_double(v(2)));
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

template <typename T>
template <typename ResidualFunc>
void Mat::BeamSMAMaterial<T>::numerical_jacobian(
    ResidualFunc&& residual, const Vec7& x, const Vec7& R, Mat7& J) const
{
  J.clear();

  for (unsigned int j = 0; j < 7; ++j)
  {
    Vec7 xp = x;
    const double h = 1.0e-8 * std::max(1.0, std::abs(x(j)));
    xp(j) += h;

    Vec7 Rp(Core::LinAlg::Initialization::zero);
    residual(xp, Rp);

    for (unsigned int i = 0; i < 7; ++i) J(i, j) = (Rp(i) - R(i)) / h;
  }
}


template <typename T>
template <typename ResidualFunc>
bool Mat::BeamSMAMaterial<T>::solve_nonlinear_system(ResidualFunc&& residual, Vec7& x) const
{
  auto residual_and_jacobian = [&](Vec7 xin)
  {
    Vec7 R(Core::LinAlg::Initialization::zero);
    residual(xin, R);

    Mat7 J(Core::LinAlg::Initialization::zero);
    numerical_jacobian(residual, xin, R, J);

    return std::make_tuple(R, J);
  };

  x = Core::Utils::solve_local_newton_damped<double, Vec7>(residual_and_jacobian, x,
      sma_params().get_local_newton_tol(),
      static_cast<unsigned>(sma_params().get_local_newton_maxiter()),
      40,       // max line-search iterations
      1.0e-10,  // minimum step length
      1.0e-4);  // Armijo parameter

  return true;
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

  extract_from_pack(buffer, gamma_conv_);
  extract_from_pack(buffer, cur_conv_);

  gamma_new_ = gamma_conv_;
  cur_new_ = cur_conv_;

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
  if (temperature_function_) return temperature_function_(current_time_);
  return sma_params().get_temperature();
}

template <typename T>
double Mat::BeamSMAMaterial<T>::compute_rs(double bs, double vS) const
{
  const double eps = sma_params().get_rs_regularization();
  const double n = sma_params().get_n_exp();
  const double vp = std::max(vS + eps, eps);
  const double omp = std::max(1.0 - vS + eps, eps);
  const double dT = current_temperature() - sma_params().get_t0_sma();

  if (bs >= 0.0)
  {
    return sma_params().get_rs_f0() + sma_params().get_h_sf() * vS +
           sma_params().get_a_sf0() * std::pow(vp, n) + sma_params().get_a_sf1() * std::pow(omp, n);
  }

  return sma_params().get_rs_r0() - sma_params().get_h_sr() * vS - sma_params().get_c_ts() * dT +
         sma_params().get_a_sr0() * std::pow(vp, n) + sma_params().get_a_sr1() * std::pow(omp, n);
}

template <typename T>
void Mat::BeamSMAMaterial<T>::compute_phase_constitutive_matrices(
    double vM, double vS, Mat3& C_N, Mat3& C_M, Mat3& dC_N, Mat3& dC_M) const
{
  C_N.clear();
  C_M.clear();
  dC_N.clear();
  dC_M.clear();

  const double m = vM + vS;

  const double EA = sma_params().get_youngs_modulus_austenite();
  const double EM = sma_params().get_youngs_modulus_martensite();
  const double GA = sma_params().get_shear_modulus_austenite();
  const double GM = sma_params().get_shear_modulus_martensite();

  const double dinvE_dm = 1.0 / EM - 1.0 / EA;
  const double dinvG_dm = 1.0 / GM - 1.0 / GA;

  const double invE = 1.0 / EA + m * dinvE_dm;
  const double invG = 1.0 / GA + m * dinvG_dm;

  // During local Newton iterations the trial phase fractions can temporarily leave the
  // admissible interval. Do not throw inside residual evaluation; keep stiffness finite.
  const double invE_safe = std::max(invE, 1.0e-30);
  const double invG_safe = std::max(invG, 1.0e-30);

  const double E = 1.0 / invE_safe;
  const double G = 1.0 / invG_safe;

  // Since single- and multiple-variant martensite are represented by the same elastic phase here,
  // dC/dvM = dC/dvS = dC/dm.
  const double dE_dm = -E * E * dinvE_dm;
  const double dG_dm = -G * G * dinvG_dm;

  const double A = sma_params().get_sma_cross_section_area();
  const double ks = sma_params().get_sma_shear_correction_factor();
  const double J = sma_params().get_sma_moment_inertia_polar();
  const double I2 = sma_params().get_sma_moment_inertia2();
  const double I3 = sma_params().get_sma_moment_inertia3();

  C_N(0, 0) = E * A;
  C_N(1, 1) = G * A * ks;
  C_N(2, 2) = G * A * ks;

  C_M(0, 0) = G * J;
  C_M(1, 1) = E * I2;
  C_M(2, 2) = E * I3;

  dC_N(0, 0) = dE_dm * A;
  dC_N(1, 1) = dG_dm * A * ks;
  dC_N(2, 2) = dG_dm * A * ks;

  dC_M(0, 0) = dG_dm * J;
  dC_M(1, 1) = dE_dm * I2;
  dC_M(2, 2) = dE_dm * I3;
}

template <typename T>
double Mat::BeamSMAMaterial<T>::quadratic_form(const Vec3& x, const Mat3& A) const
{
  double value = 0.0;
  for (unsigned int i = 0; i < 3; ++i)
    for (unsigned int j = 0; j < 3; ++j)
      value += Core::FADUtils::cast_to_double(x(i)) * Core::FADUtils::cast_to_double(A(i, j)) *
               Core::FADUtils::cast_to_double(x(j));
  return value;
}

template <typename T>
void Mat::BeamSMAMaterial<T>::choose_force_direction(
    const Vec3& Gamma, const LocalState& old_state, Vec3& dir) const
{
  dir.put_scalar(0.0);

  Mat3 C_N(Core::LinAlg::Initialization::zero);
  Mat3 C_M_dummy(Core::LinAlg::Initialization::zero);
  Mat3 dC_N_dummy(Core::LinAlg::Initialization::zero);
  Mat3 dC_M_dummy(Core::LinAlg::Initialization::zero);
  compute_phase_constitutive_matrices(Core::FADUtils::cast_to_double(old_state.vM),
      Core::FADUtils::cast_to_double(old_state.vS), C_N, C_M_dummy, dC_N_dummy, dC_M_dummy);

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

  Mat3 C_N_dummy(Core::LinAlg::Initialization::zero);
  Mat3 C_M(Core::LinAlg::Initialization::zero);
  Mat3 dC_N_dummy(Core::LinAlg::Initialization::zero);
  Mat3 dC_M_dummy(Core::LinAlg::Initialization::zero);
  compute_phase_constitutive_matrices(Core::FADUtils::cast_to_double(old_state.vM),
      Core::FADUtils::cast_to_double(old_state.vS), C_N_dummy, C_M, dC_N_dummy, dC_M_dummy);

  Vec3 kappa_tr(Core::LinAlg::Initialization::zero);
  compute_moment_transformation_curvature(
      Core::FADUtils::cast_to_double(old_state.vS), old_state.dir, kappa_tr);

  Vec3 kappa_el(Core::LinAlg::Initialization::zero);
  for (unsigned int i = 0; i < 3; ++i) kappa_el(i) = Cur(i) - kappa_tr(i);

  Vec3 Mtrial(Core::LinAlg::Initialization::zero);
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
void Mat::BeamSMAMaterial<T>::compute_force_transformation_strain(
    double vS, const Vec3& dir, Vec3& gamma_tr) const
{
  gamma_tr.put_scalar(0.0);
  gamma_tr(0) = sma_params().get_eps_l_n() * vS * dir(0);
}

template <typename T>
void Mat::BeamSMAMaterial<T>::compute_moment_transformation_curvature(
    double vS, const Vec3& dir, Vec3& kappa_tr) const
{
  kappa_tr.put_scalar(0.0);

  if (sma_params().get_kappa_l_m() <= 0.0) return;

  if (sma_params().get_torsion_sma()) kappa_tr(0) = sma_params().get_kappa_l_m() * vS * dir(0);
  kappa_tr(1) = sma_params().get_kappa_l_m() * vS * dir(1);
  kappa_tr(2) = sma_params().get_kappa_l_m() * vS * dir(2);
}

template <typename T>
void Mat::BeamSMAMaterial<T>::assemble_force_residual(const Vec7& x, const LocalState& old_state,
    const Vec3& Gamma, const Vec3& dir, Vec7& R, Vec3* stress) const
{
  const double vM = x(0);
  const double dlamM = x(1);
  const double cM0 = x(2);
  const double vS = x(3);
  const double dlamS = x(4);
  const double cS0 = x(5);
  const double cMS = x(6);

  Mat3 C_N(Core::LinAlg::Initialization::zero);
  Mat3 C_M_dummy(Core::LinAlg::Initialization::zero);
  Mat3 dC_N(Core::LinAlg::Initialization::zero);
  Mat3 dC_M_dummy(Core::LinAlg::Initialization::zero);
  compute_phase_constitutive_matrices(vM, vS, C_N, C_M_dummy, dC_N, dC_M_dummy);

  Vec3 gamma_tr(Core::LinAlg::Initialization::zero);
  compute_force_transformation_strain(vS, dir, gamma_tr);

  Vec3 gamma_el(Core::LinAlg::Initialization::zero);
  for (unsigned int i = 0; i < 3; ++i) gamma_el(i) = Gamma(i) - gamma_tr(i);

  Vec3 N(Core::LinAlg::Initialization::zero);
  N.multiply(C_N, gamma_el);
  if (stress != nullptr) *stress = N;

  const double dpsi_el_dv = 0.5 * quadratic_form(gamma_el, dC_N);

  const double therm =
      sma_params().get_delta_s_ams() * (current_temperature() - sma_params().get_t0_sma()) +
      sma_params().get_w_in() * (1.0 - 2.0 * vM - 2.0 * vS);

  const double cM = cM0 + cMS;
  const double cS = cS0 + cMS;

  // B = -dW/dv.  The phase-dependent stiffness adds -0.5*q_el^T*(dC/dv)*q_el.
  const double BM = -dpsi_el_dv - therm - cM;
  const double BS = sma_params().get_eps_l_n() * dot(dir, N) - dpsi_el_dv - therm - cS;

  const double absBM = abs_reg(BM, sma_params().get_fb_regularization());
  const double absBS = abs_reg(BS, sma_params().get_fb_regularization());

  const double FM = absBM - sma_params().get_r_m();
  const double FS = absBS - compute_rs(BS, vS);

  R(0) = vM - Core::FADUtils::cast_to_double(old_state.vM) - dlamM * BM / absBM;
  R(1) = fb(FM, dlamM, sma_params().get_fb_regularization());
  R(2) = fb(cM0, vM, sma_params().get_fb_regularization());
  R(3) = vS - Core::FADUtils::cast_to_double(old_state.vS) - dlamS * BS / absBS;
  R(4) = fb(FS, dlamS, sma_params().get_fb_regularization());
  R(5) = fb(cS0, vS, sma_params().get_fb_regularization());
  R(6) = fb(vM + vS - 1.0, cMS, sma_params().get_fb_regularization());
}

template <typename T>
void Mat::BeamSMAMaterial<T>::assemble_moment_residual(const Vec7& x, const LocalState& old_state,
    const Vec3& Cur, const Vec3& dir, Vec7& R, Vec3* stress) const
{
  const double vM = x(0);
  const double dlamM = x(1);
  const double cM0 = x(2);
  const double vS = x(3);
  const double dlamS = x(4);
  const double cS0 = x(5);
  const double cMS = x(6);

  Mat3 C_N_dummy(Core::LinAlg::Initialization::zero);
  Mat3 C_M(Core::LinAlg::Initialization::zero);
  Mat3 dC_N_dummy(Core::LinAlg::Initialization::zero);
  Mat3 dC_M(Core::LinAlg::Initialization::zero);
  compute_phase_constitutive_matrices(vM, vS, C_N_dummy, C_M, dC_N_dummy, dC_M);

  Vec3 kappa_tr(Core::LinAlg::Initialization::zero);
  compute_moment_transformation_curvature(vS, dir, kappa_tr);

  Vec3 kappa_el(Core::LinAlg::Initialization::zero);
  for (unsigned int i = 0; i < 3; ++i) kappa_el(i) = Cur(i) - kappa_tr(i);

  Vec3 M(Core::LinAlg::Initialization::zero);
  M.multiply(C_M, kappa_el);
  if (stress != nullptr) *stress = M;

  const double dpsi_el_dv = 0.5 * quadratic_form(kappa_el, dC_M);

  const double therm =
      sma_params().get_delta_s_ams() * (current_temperature() - sma_params().get_t0_sma()) +
      sma_params().get_w_in() * (1.0 - 2.0 * vM - 2.0 * vS);

  const double cM = cM0 + cMS;
  const double cS = cS0 + cMS;

  const double BM = -dpsi_el_dv - therm - cM;
  const double BS = sma_params().get_kappa_l_m() * dot(dir, M) - dpsi_el_dv - therm - cS;

  const double absBM = abs_reg(BM, sma_params().get_fb_regularization());
  const double absBS = abs_reg(BS, sma_params().get_fb_regularization());

  const double FM = absBM - sma_params().get_r_m();
  const double FS = absBS - compute_rs(BS, vS);

  R(0) = vM - Core::FADUtils::cast_to_double(old_state.vM) - dlamM * BM / absBM;
  R(1) = fb(FM, dlamM, sma_params().get_fb_regularization());
  R(2) = fb(cM0, vM, sma_params().get_fb_regularization());
  R(3) = vS - Core::FADUtils::cast_to_double(old_state.vS) - dlamS * BS / absBS;
  R(4) = fb(FS, dlamS, sma_params().get_fb_regularization());
  R(5) = fb(cS0, vS, sma_params().get_fb_regularization());
  R(6) = fb(vM + vS - 1.0, cMS, sma_params().get_fb_regularization());
}

/*-----------------------------------------------------------------------------------------------*/
/* Local updates                                                                                 */
/*-----------------------------------------------------------------------------------------------*/
template <typename T>
bool Mat::BeamSMAMaterial<T>::solve_force_state(
    const Vec3& Gamma, const LocalState& old_state, LocalState& new_state, Vec3& stress) const
{
  if (sma_params().get_eps_l_n() <= 0.0)
  {
    new_state = old_state;
    Mat3 C_N(Core::LinAlg::Initialization::zero);
    Mat3 C_M_dummy(Core::LinAlg::Initialization::zero);
    Mat3 dC_N_dummy(Core::LinAlg::Initialization::zero);
    Mat3 dC_M_dummy(Core::LinAlg::Initialization::zero);
    compute_phase_constitutive_matrices(Core::FADUtils::cast_to_double(old_state.vM),
        Core::FADUtils::cast_to_double(old_state.vS), C_N, C_M_dummy, dC_N_dummy, dC_M_dummy);
    stress.multiply(C_N, Gamma);
    return true;
  }

  Vec3 dir(Core::LinAlg::Initialization::zero);
  choose_force_direction(Gamma, old_state, dir);

  Vec7 Rpred(Core::LinAlg::Initialization::zero);
  Vec7 xpred(Core::LinAlg::Initialization::zero);
  xpred(0) = Core::FADUtils::cast_to_double(old_state.vM);
  xpred(1) = 0.0;
  xpred(2) = 0.0;
  xpred(3) = Core::FADUtils::cast_to_double(old_state.vS);
  xpred(4) = 0.0;
  xpred(5) = 0.0;
  xpred(6) = 0.0;

  Vec3 stress_pred(Core::LinAlg::Initialization::zero);
  assemble_force_residual(xpred, old_state, Gamma, dir, Rpred, &stress_pred);

  // Use the same residual quantities for the elastic predictor as in the final residual.
  const double FM_pred = Rpred(1);
  const double FS_pred = Rpred(4);

  if (FM_pred <= sma_params().get_local_newton_tol() &&
      FS_pred <= sma_params().get_local_newton_tol() && xpred(0) >= -1.0e-12 &&
      xpred(3) >= -1.0e-12 && xpred(0) + xpred(3) <= 1.0 + 1.0e-12)
  {
    new_state = old_state;
    new_state.dir = dir;
    stress = stress_pred;
    return true;
  }

  Vec7 x(Core::LinAlg::Initialization::zero);
  x(0) = Core::FADUtils::cast_to_double(old_state.vM);
  x(1) = 0.0;
  x(2) = (x(0) <= 1.0e-12 ? -1.0e-8 : 0.0);
  x(3) = Core::FADUtils::cast_to_double(old_state.vS);
  x(4) = 0.0;
  x(5) = (x(3) <= 1.0e-12 ? -1.0e-8 : 0.0);
  x(6) = (x(0) + x(3) >= 1.0 - 1.0e-12 ? 1.0e-8 : 0.0);

  auto residual = [&](const Vec7& xin, Vec7& Rout)
  { assemble_force_residual(xin, old_state, Gamma, dir, Rout, nullptr); };

  if (!solve_nonlinear_system(residual, x)) return false;

  new_state.vM = std::max(0.0, std::min(1.0, x(0)));
  new_state.vS = std::max(0.0, std::min(1.0, x(3)));
  new_state.dir = dir;

  Vec7 Rfinal(Core::LinAlg::Initialization::zero);
  assemble_force_residual(x, old_state, Gamma, dir, Rfinal, &stress);
  return true;
}

template <typename T>
bool Mat::BeamSMAMaterial<T>::solve_moment_state(
    const Vec3& Cur, const LocalState& old_state, LocalState& new_state, Vec3& stress) const
{
  if (sma_params().get_kappa_l_m() <= 0.0)
  {
    new_state = old_state;
    Mat3 C_N_dummy(Core::LinAlg::Initialization::zero);
    Mat3 C_M(Core::LinAlg::Initialization::zero);
    Mat3 dC_N_dummy(Core::LinAlg::Initialization::zero);
    Mat3 dC_M_dummy(Core::LinAlg::Initialization::zero);
    compute_phase_constitutive_matrices(Core::FADUtils::cast_to_double(old_state.vM),
        Core::FADUtils::cast_to_double(old_state.vS), C_N_dummy, C_M, dC_N_dummy, dC_M_dummy);
    stress.multiply(C_M, Cur);
    return true;
  }

  Vec3 dir(Core::LinAlg::Initialization::zero);
  choose_moment_direction(Cur, old_state, dir);

  Vec7 xpred(Core::LinAlg::Initialization::zero);
  xpred(0) = Core::FADUtils::cast_to_double(old_state.vM);
  xpred(1) = 0.0;
  xpred(2) = 0.0;
  xpred(3) = Core::FADUtils::cast_to_double(old_state.vS);
  xpred(4) = 0.0;
  xpred(5) = 0.0;
  xpred(6) = 0.0;

  Vec7 Rpred(Core::LinAlg::Initialization::zero);
  Vec3 stress_pred(Core::LinAlg::Initialization::zero);
  assemble_moment_residual(xpred, old_state, Cur, dir, Rpred, &stress_pred);

  const double FM_pred = Rpred(1);
  const double FS_pred = Rpred(4);

  if (FM_pred <= sma_params().get_local_newton_tol() &&
      FS_pred <= sma_params().get_local_newton_tol() && xpred(0) >= -1.0e-12 &&
      xpred(3) >= -1.0e-12 && xpred(0) + xpred(3) <= 1.0 + 1.0e-12)
  {
    new_state = old_state;
    new_state.dir = dir;
    stress = stress_pred;
    return true;
  }

  Vec7 x(Core::LinAlg::Initialization::zero);
  x(0) = Core::FADUtils::cast_to_double(old_state.vM);
  x(1) = 0.0;
  x(2) = (x(0) <= 1.0e-12 ? -1.0e-8 : 0.0);
  x(3) = Core::FADUtils::cast_to_double(old_state.vS);
  x(4) = 0.0;
  x(5) = (x(3) <= 1.0e-12 ? -1.0e-8 : 0.0);
  x(6) = (x(0) + x(3) >= 1.0 - 1.0e-12 ? 1.0e-8 : 0.0);

  auto residual = [&](const Vec7& xin, Vec7& Rout)
  { assemble_moment_residual(xin, old_state, Cur, dir, Rout, nullptr); };

  if (!solve_nonlinear_system(residual, x)) return false;

  new_state.vM = std::max(0.0, std::min(1.0, x(0)));
  new_state.vS = std::max(0.0, std::min(1.0, x(3)));
  new_state.dir = dir;

  Vec7 Rfinal(Core::LinAlg::Initialization::zero);
  assemble_moment_residual(x, old_state, Cur, dir, Rfinal, &stress);
  return true;
}


template <typename T>
void Mat::BeamSMAMaterial<T>::local_force_response(
    const Vec3& Gamma, const unsigned int gp, Vec3& stress, Mat3& C_alg)
{
  const LocalState state_start = force_state_conv_[gp];
  const Vec3 gamma_start = gamma_conv_[gp];

  LocalState state_end = state_start;
  int nsub = 1;

  if (!integrate_force_state_substepped(gamma_start, Gamma, state_start, state_end, stress, nsub))
    FOUR_C_THROW("Beam SMA force update failed.");

  force_state_new_[gp] = state_end;
  gamma_new_[gp] = Gamma;

  Mat3 C_N(Core::LinAlg::Initialization::zero);
  Mat3 C_M_dummy(Core::LinAlg::Initialization::zero);
  Mat3 dC_N_dummy(Core::LinAlg::Initialization::zero);
  Mat3 dC_M_dummy(Core::LinAlg::Initialization::zero);
  compute_phase_constitutive_matrices(Core::FADUtils::cast_to_double(state_end.vM),
      Core::FADUtils::cast_to_double(state_end.vS), C_N, C_M_dummy, dC_N_dummy, dC_M_dummy);

  // If substepping was necessary, return the frozen-phase tangent for robustness.
  // A fully consistent tangent would need to differentiate the entire substepped integration path.
  if (!sma_params().use_consistent_tangent() || nsub > 1)
  {
    C_alg = C_N;
    return;
  }

  C_alg.clear();
  for (unsigned int j = 0; j < 3; ++j)
  {
    const double h = 1.0e-8 * std::max(1.0, std::abs(Core::FADUtils::cast_to_double(Gamma(j))));

    Vec3 Gp = Gamma;
    Gp(j) += h;

    LocalState pert_state = state_start;
    Vec3 stress_pert(Core::LinAlg::Initialization::zero);
    int nsub_pert = 1;

    if (!integrate_force_state_substepped(
            gamma_start, Gp, state_start, pert_state, stress_pert, nsub_pert))
      FOUR_C_THROW("Beam SMA force tangent update failed.");

    for (unsigned int i = 0; i < 3; ++i) C_alg(i, j) = (stress_pert(i) - stress(i)) / h;
  }
}

template <typename T>
bool Mat::BeamSMAMaterial<T>::integrate_force_state_substepped(const Vec3& Gamma_start,
    const Vec3& Gamma_end, const LocalState& state_start, LocalState& state_end, Vec3& stress_end,
    int& nsub) const
{
  const double dgamma_norm = vector_difference_norm(Gamma_end, Gamma_start);

  nsub = 1;
  const double limit = sma_params().get_local_substep_strain_limit();
  if (limit > 0.0)
  {
    nsub = static_cast<int>(std::ceil(dgamma_norm / limit));
    nsub = std::max(1, std::min(nsub, sma_params().get_local_substep_max()));
  }

  LocalState state_old = state_start;
  LocalState state_new = state_start;
  Vec3 stress_sub(Core::LinAlg::Initialization::zero);

  for (int isub = 1; isub <= nsub; ++isub)
  {
    const double alpha = static_cast<double>(isub) / static_cast<double>(nsub);

    Vec3 Gamma_sub(Core::LinAlg::Initialization::zero);
    for (unsigned int i = 0; i < 3; ++i)
      Gamma_sub(i) = Gamma_start(i) + alpha * (Gamma_end(i) - Gamma_start(i));

    if (!solve_force_state(Gamma_sub, state_old, state_new, stress_sub))
    {
      FOUR_C_THROW("Beam SMA force substep update failed at substep {} of {}.", isub, nsub);
    }

    state_old = state_new;
  }

  state_end = state_new;
  stress_end = stress_sub;
  return true;
}

template <typename T>
void Mat::BeamSMAMaterial<T>::local_moment_response(
    const Vec3& Cur, const unsigned int gp, Vec3& stress, Mat3& C_alg)
{
  const LocalState state_start = moment_state_conv_[gp];
  const Vec3 cur_start = cur_conv_[gp];

  LocalState state_end = state_start;
  int nsub = 1;

  if (!integrate_moment_state_substepped(cur_start, Cur, state_start, state_end, stress, nsub))
    FOUR_C_THROW("Beam SMA moment update failed.");

  moment_state_new_[gp] = state_end;
  cur_new_[gp] = Cur;

  Mat3 C_N_dummy(Core::LinAlg::Initialization::zero);
  Mat3 C_M(Core::LinAlg::Initialization::zero);
  Mat3 dC_N_dummy(Core::LinAlg::Initialization::zero);
  Mat3 dC_M_dummy(Core::LinAlg::Initialization::zero);
  compute_phase_constitutive_matrices(Core::FADUtils::cast_to_double(state_end.vM),
      Core::FADUtils::cast_to_double(state_end.vS), C_N_dummy, C_M, dC_N_dummy, dC_M_dummy);

  // If substepping was necessary, return the frozen-phase tangent for robustness.
  if (!sma_params().use_consistent_tangent() || nsub > 1)
  {
    C_alg = C_M;
    return;
  }

  C_alg.clear();
  for (unsigned int j = 0; j < 3; ++j)
  {
    const double h = 1.0e-8 * std::max(1.0, std::abs(Core::FADUtils::cast_to_double(Cur(j))));

    Vec3 Cp = Cur;
    Cp(j) += h;

    LocalState pert_state = state_start;
    Vec3 stress_pert(Core::LinAlg::Initialization::zero);
    int nsub_pert = 1;

    if (!integrate_moment_state_substepped(
            cur_start, Cp, state_start, pert_state, stress_pert, nsub_pert))
      FOUR_C_THROW("Beam SMA moment tangent update failed.");

    for (unsigned int i = 0; i < 3; ++i) C_alg(i, j) = (stress_pert(i) - stress(i)) / h;
  }
}

template <typename T>
bool Mat::BeamSMAMaterial<T>::integrate_moment_state_substepped(const Vec3& Cur_start,
    const Vec3& Cur_end, const LocalState& state_start, LocalState& state_end, Vec3& stress_end,
    int& nsub) const
{
  const double dcur_norm = vector_difference_norm(Cur_end, Cur_start);

  nsub = 1;
  const double limit = sma_params().get_local_substep_curvature_limit();
  if (limit > 0.0)
  {
    nsub = static_cast<int>(std::ceil(dcur_norm / limit));
    nsub = std::max(1, std::min(nsub, sma_params().get_local_substep_max()));
  }

  LocalState state_old = state_start;
  LocalState state_new = state_start;
  Vec3 stress_sub(Core::LinAlg::Initialization::zero);

  for (int isub = 1; isub <= nsub; ++isub)
  {
    const double alpha = static_cast<double>(isub) / static_cast<double>(nsub);

    Vec3 Cur_sub(Core::LinAlg::Initialization::zero);
    for (unsigned int i = 0; i < 3; ++i)
      Cur_sub(i) = Cur_start(i) + alpha * (Cur_end(i) - Cur_start(i));

    if (!solve_moment_state(Cur_sub, state_old, state_new, stress_sub))
    {
      FOUR_C_THROW("Beam SMA moment substep update failed at substep {} of {}.", isub, nsub);
    }

    state_old = state_new;
  }

  state_end = state_new;
  stress_end = stress_sub;
  return true;
}


/*-----------------------------------------------------------------------------------------------*/
/* Interface methods                                                                             */
/*-----------------------------------------------------------------------------------------------*/
template <typename T>
void Mat::BeamSMAMaterial<T>::get_constitutive_matrix_of_forces_material_frame(
    Core::LinAlg::Matrix<3, 3, T>& C_N) const
{
  Mat3 C_M_dummy(Core::LinAlg::Initialization::zero);
  Mat3 dC_N_dummy(Core::LinAlg::Initialization::zero);
  Mat3 dC_M_dummy(Core::LinAlg::Initialization::zero);
  compute_phase_constitutive_matrices(0.0, 0.0, C_N, C_M_dummy, dC_N_dummy, dC_M_dummy);
}

template <typename T>
void Mat::BeamSMAMaterial<T>::get_constitutive_matrix_of_moments_material_frame(
    Core::LinAlg::Matrix<3, 3, T>& C_M) const
{
  Mat3 C_N_dummy(Core::LinAlg::Initialization::zero);
  Mat3 dC_N_dummy(Core::LinAlg::Initialization::zero);
  Mat3 dC_M_dummy(Core::LinAlg::Initialization::zero);
  compute_phase_constitutive_matrices(0.0, 0.0, C_N_dummy, C_M, dC_N_dummy, dC_M_dummy);
}

template <typename T>
void Mat::BeamSMAMaterial<T>::compute_constitutive_parameter(
    Core::LinAlg::Matrix<3, 3, T>& C_N, Core::LinAlg::Matrix<3, 3, T>& C_M)
{
  Mat3 dC_N_dummy(Core::LinAlg::Initialization::zero);
  Mat3 dC_M_dummy(Core::LinAlg::Initialization::zero);
  compute_phase_constitutive_matrices(0.0, 0.0, C_N, C_M, dC_N_dummy, dC_M_dummy);

  for (unsigned int gp = 0; gp < numgp_force_; ++gp) c_n_alg_[gp] = C_N;
  for (unsigned int gp = 0; gp < numgp_moment_; ++gp) c_m_alg_[gp] = C_M;
}

template <typename T>
void Mat::BeamSMAMaterial<T>::evaluate_force_contributions_to_stress(
    Core::LinAlg::Matrix<3, 1, T>& stressN, const Core::LinAlg::Matrix<3, 3, T>&,
    const Core::LinAlg::Matrix<3, 1, T>& Gamma, const unsigned int gp)
{
  Mat3 C_alg(Core::LinAlg::Initialization::zero);
  local_force_response(Gamma, gp, stressN, C_alg);
  c_n_alg_[gp] = C_alg;
}

template <typename T>
void Mat::BeamSMAMaterial<T>::evaluate_moment_contributions_to_stress(
    Core::LinAlg::Matrix<3, 1, T>& stressM, const Core::LinAlg::Matrix<3, 3, T>&,
    const Core::LinAlg::Matrix<3, 1, T>& Cur, const unsigned int gp)
{
  Mat3 C_alg(Core::LinAlg::Initialization::zero);
  local_moment_response(Cur, gp, stressM, C_alg);
  c_m_alg_[gp] = C_alg;
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

  gamma_conv_ = gamma_new_;
  cur_conv_ = cur_new_;
}

template <typename T>
void Mat::BeamSMAMaterial<T>::reset()
{
  force_state_new_ = force_state_conv_;
  moment_state_new_ = moment_state_conv_;

  gamma_new_ = gamma_conv_;
  cur_new_ = cur_conv_;
}

/*-----------------------------------------------------------------------------------------------*/
/* Explicit instantiations                                                                       */
/*-----------------------------------------------------------------------------------------------*/
template class Mat::BeamSMAMaterial<double>;
template class Mat::BeamSMAMaterialType<double>;

FOUR_C_NAMESPACE_CLOSE