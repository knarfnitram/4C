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
#include <limits>

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
      local_newton_maxiter_(matdata.parameters.get<int>("LOCALITER"))
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
double Mat::BeamSMAMaterial<T>::dot(const Vec3& a, const Vec3& b)
{
  return Core::FADUtils::cast_to_double(a(0)) * Core::FADUtils::cast_to_double(b(0)) +
         Core::FADUtils::cast_to_double(a(1)) * Core::FADUtils::cast_to_double(b(1)) +
         Core::FADUtils::cast_to_double(a(2)) * Core::FADUtils::cast_to_double(b(2));
}

template <typename T>
double Mat::BeamSMAMaterial<T>::residual_norm(const Vec7& r)
{
  double n2 = 0.0;
  for (double ri : r) n2 += ri * ri;
  return std::sqrt(n2);
}

template <typename T>
bool Mat::BeamSMAMaterial<T>::solve_linear_7x7(Mat7& A, Vec7& b)
{
  constexpr int n = 7;
  for (int k = 0; k < n; ++k)
  {
    int pivot = k;
    double maxabs = std::abs(A[k][k]);
    for (int i = k + 1; i < n; ++i)
    {
      const double cand = std::abs(A[i][k]);
      if (cand > maxabs)
      {
        maxabs = cand;
        pivot = i;
      }
    }

    if (maxabs < 1.0e-14) return false;

    if (pivot != k)
    {
      std::swap(A[pivot], A[k]);
      std::swap(b[pivot], b[k]);
    }

    for (int i = k + 1; i < n; ++i)
    {
      const double factor = A[i][k] / A[k][k];
      A[i][k] = 0.0;
      for (int j = k + 1; j < n; ++j) A[i][j] -= factor * A[k][j];
      b[i] -= factor * b[k];
    }
  }

  for (int i = n - 1; i >= 0; --i)
  {
    double sum = b[i];
    for (int j = i + 1; j < n; ++j) sum -= A[i][j] * b[j];
    b[i] = sum / A[i][i];
  }

  return true;
}

template <typename T>
template <typename ResidualFunc>
void Mat::BeamSMAMaterial<T>::numerical_jacobian(
    ResidualFunc&& residual, const Vec7& x, const Vec7& R, Mat7& J) const
{
  for (auto& row : J) row.fill(0.0);

  for (unsigned int j = 0; j < 7; ++j)
  {
    Vec7 xp = x;
    const double h = 1.0e-8 * std::max(1.0, std::abs(x[j]));
    xp[j] += h;

    Vec7 Rp{};
    residual(xp, Rp);

    for (unsigned int i = 0; i < 7; ++i) J[i][j] = (Rp[i] - R[i]) / h;
  }
}

template <typename T>
template <typename ResidualFunc>
bool Mat::BeamSMAMaterial<T>::solve_nonlinear_system(ResidualFunc&& residual, Vec7& x) const
{
  const double tol = sma_params().get_local_newton_tol();
  const int maxiter = sma_params().get_local_newton_maxiter();

  Vec7 R{};
  residual(x, R);
  double nR = residual_norm(R);
  if (nR < tol) return true;

  for (int iter = 0; iter < maxiter; ++iter)
  {
    Mat7 J{};
    numerical_jacobian(residual, x, R, J);

    Vec7 dx{};
    for (unsigned int i = 0; i < 7; ++i) dx[i] = -R[i];

    if (!solve_linear_7x7(J, dx)) return false;

    bool accepted = false;
    double alpha = 1.0;
    Vec7 xtrial{};
    Vec7 Rtrial{};

    for (int ls = 0; ls < 20; ++ls)
    {
      for (unsigned int i = 0; i < 7; ++i) xtrial[i] = x[i] + alpha * dx[i];
      residual(xtrial, Rtrial);

      const double nTrial = residual_norm(Rtrial);
      if (nTrial <= (1.0 - 1.0e-4 * alpha) * nR)
      {
        x = xtrial;
        R = Rtrial;
        nR = nTrial;
        accepted = true;
        break;
      }

      alpha *= 0.5;
    }

    if (!accepted) return false;
    if (nR < tol) return true;
  }

  return false;
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
double Mat::BeamSMAMaterial<T>::compute_rs(double bs, double vS) const
{
  const double eps = sma_params().get_rs_regularization();
  const double n = sma_params().get_n_exp();
  const double vp = std::max(vS + eps, eps);
  const double omp = std::max(1.0 - vS + eps, eps);
  const double dT = sma_params().get_temperature() - sma_params().get_t0_sma();

  if (bs >= 0.0)
  {
    return sma_params().get_rs_f0() + sma_params().get_h_sf() * vS +
           sma_params().get_a_sf0() * std::pow(vp, n) + sma_params().get_a_sf1() * std::pow(omp, n);
  }

  return sma_params().get_rs_r0() - sma_params().get_h_sr() * vS - sma_params().get_c_ts() * dT +
         sma_params().get_a_sr0() * std::pow(vp, n) + sma_params().get_a_sr1() * std::pow(omp, n);
}

template <typename T>
void Mat::BeamSMAMaterial<T>::choose_force_direction(
    const Vec3& Gamma, const LocalState& old_state, Vec3& dir) const
{
  dir.put_scalar(0.0);

  const double axial_trial =
      this->params().get_axial_rigidity() *
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
    const Vec3& Cur, const Mat3& C_M, const LocalState& old_state, Vec3& dir) const
{
  dir.put_scalar(0.0);

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
    const Vec3& Gamma, const Mat3& C_N, const Vec3& dir, Vec7& R, Vec3* stress) const
{
  const double vM = x[0];
  const double dlamM = x[1];
  const double cM0 = x[2];
  const double vS = x[3];
  const double dlamS = x[4];
  const double cS0 = x[5];
  const double cMS = x[6];

  Vec3 gamma_tr(Core::LinAlg::Initialization::zero);
  compute_force_transformation_strain(vS, dir, gamma_tr);

  Vec3 gamma_el(Core::LinAlg::Initialization::zero);
  for (unsigned int i = 0; i < 3; ++i) gamma_el(i) = Gamma(i) - gamma_tr(i);

  Vec3 N(Core::LinAlg::Initialization::zero);
  N.multiply(C_N, gamma_el);
  if (stress != nullptr) *stress = N;

  const double therm = sma_params().get_delta_s_ams() *
                           (sma_params().get_temperature() - sma_params().get_t0_sma()) +
                       sma_params().get_w_in() * (1.0 - 2.0 * vM - 2.0 * vS);

  const double cM = cM0 + cMS;
  const double cS = cS0 + cMS;

  const double BM = -therm - cM;
  const double BS = sma_params().get_eps_l_n() * dot(dir, N) - therm - cS;

  const double absBM = abs_reg(BM, sma_params().get_fb_regularization());
  const double absBS = abs_reg(BS, sma_params().get_fb_regularization());

  const double FM = absBM - sma_params().get_r_m();
  const double FS = absBS - compute_rs(BS, vS);

  R[0] = vM - Core::FADUtils::cast_to_double(old_state.vM) - dlamM * BM / absBM;
  R[1] = fb(FM, dlamM, sma_params().get_fb_regularization());
  R[2] = fb(cM0, vM, sma_params().get_fb_regularization());
  R[3] = vS - Core::FADUtils::cast_to_double(old_state.vS) - dlamS * BS / absBS;
  R[4] = fb(FS, dlamS, sma_params().get_fb_regularization());
  R[5] = fb(cS0, vS, sma_params().get_fb_regularization());
  R[6] = fb(vM + vS - 1.0, cMS, sma_params().get_fb_regularization());
}

template <typename T>
void Mat::BeamSMAMaterial<T>::assemble_moment_residual(const Vec7& x, const LocalState& old_state,
    const Vec3& Cur, const Mat3& C_M, const Vec3& dir, Vec7& R, Vec3* stress) const
{
  const double vM = x[0];
  const double dlamM = x[1];
  const double cM0 = x[2];
  const double vS = x[3];
  const double dlamS = x[4];
  const double cS0 = x[5];
  const double cMS = x[6];

  Vec3 kappa_tr(Core::LinAlg::Initialization::zero);
  compute_moment_transformation_curvature(vS, dir, kappa_tr);

  Vec3 kappa_el(Core::LinAlg::Initialization::zero);
  for (unsigned int i = 0; i < 3; ++i) kappa_el(i) = Cur(i) - kappa_tr(i);

  Vec3 M(Core::LinAlg::Initialization::zero);
  M.multiply(C_M, kappa_el);
  if (stress != nullptr) *stress = M;

  const double therm = sma_params().get_delta_s_ams() *
                           (sma_params().get_temperature() - sma_params().get_t0_sma()) +
                       sma_params().get_w_in() * (1.0 - 2.0 * vM - 2.0 * vS);

  const double cM = cM0 + cMS;
  const double cS = cS0 + cMS;

  const double BM = -therm - cM;
  const double BS = sma_params().get_kappa_l_m() * dot(dir, M) - therm - cS;

  const double absBM = abs_reg(BM, sma_params().get_fb_regularization());
  const double absBS = abs_reg(BS, sma_params().get_fb_regularization());

  const double FM = absBM - sma_params().get_r_m();
  const double FS = absBS - compute_rs(BS, vS);

  R[0] = vM - Core::FADUtils::cast_to_double(old_state.vM) - dlamM * BM / absBM;
  R[1] = fb(FM, dlamM, sma_params().get_fb_regularization());
  R[2] = fb(cM0, vM, sma_params().get_fb_regularization());
  R[3] = vS - Core::FADUtils::cast_to_double(old_state.vS) - dlamS * BS / absBS;
  R[4] = fb(FS, dlamS, sma_params().get_fb_regularization());
  R[5] = fb(cS0, vS, sma_params().get_fb_regularization());
  R[6] = fb(vM + vS - 1.0, cMS, sma_params().get_fb_regularization());
}

/*-----------------------------------------------------------------------------------------------*/
/* Local updates                                                                                 */
/*-----------------------------------------------------------------------------------------------*/
template <typename T>
bool Mat::BeamSMAMaterial<T>::solve_force_state(const Vec3& Gamma, const Mat3& C_N,
    const LocalState& old_state, LocalState& new_state, Vec3& stress) const
{
  if (sma_params().get_eps_l_n() <= 0.0)
  {
    new_state = old_state;
    Vec3 gamma_el(Core::LinAlg::Initialization::zero);
    for (unsigned int i = 0; i < 3; ++i) gamma_el(i) = Gamma(i);
    stress.multiply(C_N, gamma_el);
    return true;
  }

  Vec3 dir(Core::LinAlg::Initialization::zero);
  choose_force_direction(Gamma, old_state, dir);

  Vec7 Rpred{};
  Vec7 xpred{};
  xpred[0] = Core::FADUtils::cast_to_double(old_state.vM);
  xpred[1] = 0.0;
  xpred[2] = 0.0;
  xpred[3] = Core::FADUtils::cast_to_double(old_state.vS);
  xpred[4] = 0.0;
  xpred[5] = 0.0;
  xpred[6] = 0.0;

  Vec3 stress_pred(Core::LinAlg::Initialization::zero);
  assemble_force_residual(xpred, old_state, Gamma, C_N, dir, Rpred, &stress_pred);

  const double therm_pred = sma_params().get_delta_s_ams() *
                                (sma_params().get_temperature() - sma_params().get_t0_sma()) +
                            sma_params().get_w_in() * (1.0 - 2.0 * xpred[0] - 2.0 * xpred[3]);
  const double BM_pred = -therm_pred;
  const double BS_pred = sma_params().get_eps_l_n() * dot(dir, stress_pred) - therm_pred;
  const double FM_pred =
      abs_reg(BM_pred, sma_params().get_fb_regularization()) - sma_params().get_r_m();
  const double FS_pred =
      abs_reg(BS_pred, sma_params().get_fb_regularization()) - compute_rs(BS_pred, xpred[3]);

  if (FM_pred <= sma_params().get_local_newton_tol() &&
      FS_pred <= sma_params().get_local_newton_tol() && xpred[0] >= -1.0e-12 &&
      xpred[3] >= -1.0e-12 && xpred[0] + xpred[3] <= 1.0 + 1.0e-12)
  {
    new_state = old_state;
    new_state.dir = dir;
    stress = stress_pred;
    return true;
  }

  Vec7 x{};
  x[0] = Core::FADUtils::cast_to_double(old_state.vM);
  x[1] = 0.0;
  x[2] = (x[0] <= 1.0e-12 ? -1.0e-8 : 0.0);
  x[3] = Core::FADUtils::cast_to_double(old_state.vS);
  x[4] = 0.0;
  x[5] = (x[3] <= 1.0e-12 ? -1.0e-8 : 0.0);
  x[6] = (x[0] + x[3] >= 1.0 - 1.0e-12 ? 1.0e-8 : 0.0);

  auto residual = [&](const Vec7& xin, Vec7& Rout)
  { assemble_force_residual(xin, old_state, Gamma, C_N, dir, Rout, nullptr); };

  if (!solve_nonlinear_system(residual, x)) return false;

  new_state.vM = std::max(0.0, std::min(1.0, x[0]));
  new_state.vS = std::max(0.0, std::min(1.0, x[3]));
  new_state.dir = dir;

  Vec7 Rfinal{};
  assemble_force_residual(x, old_state, Gamma, C_N, dir, Rfinal, &stress);
  return true;
}

template <typename T>
bool Mat::BeamSMAMaterial<T>::solve_moment_state(const Vec3& Cur, const Mat3& C_M,
    const LocalState& old_state, LocalState& new_state, Vec3& stress) const
{
  if (sma_params().get_kappa_l_m() <= 0.0)
  {
    new_state = old_state;
    Vec3 kappa_el(Core::LinAlg::Initialization::zero);
    for (unsigned int i = 0; i < 3; ++i) kappa_el(i) = Cur(i);
    stress.multiply(C_M, kappa_el);
    return true;
  }

  Vec3 dir(Core::LinAlg::Initialization::zero);
  choose_moment_direction(Cur, C_M, old_state, dir);

  Vec7 xpred{};
  xpred[0] = Core::FADUtils::cast_to_double(old_state.vM);
  xpred[1] = 0.0;
  xpred[2] = 0.0;
  xpred[3] = Core::FADUtils::cast_to_double(old_state.vS);
  xpred[4] = 0.0;
  xpred[5] = 0.0;
  xpred[6] = 0.0;

  Vec7 Rpred{};
  Vec3 stress_pred(Core::LinAlg::Initialization::zero);
  assemble_moment_residual(xpred, old_state, Cur, C_M, dir, Rpred, &stress_pred);

  const double therm_pred = sma_params().get_delta_s_ams() *
                                (sma_params().get_temperature() - sma_params().get_t0_sma()) +
                            sma_params().get_w_in() * (1.0 - 2.0 * xpred[0] - 2.0 * xpred[3]);
  const double BM_pred = -therm_pred;
  const double BS_pred = sma_params().get_kappa_l_m() * dot(dir, stress_pred) - therm_pred;
  const double FM_pred =
      abs_reg(BM_pred, sma_params().get_fb_regularization()) - sma_params().get_r_m();
  const double FS_pred =
      abs_reg(BS_pred, sma_params().get_fb_regularization()) - compute_rs(BS_pred, xpred[3]);

  if (FM_pred <= sma_params().get_local_newton_tol() &&
      FS_pred <= sma_params().get_local_newton_tol() && xpred[0] >= -1.0e-12 &&
      xpred[3] >= -1.0e-12 && xpred[0] + xpred[3] <= 1.0 + 1.0e-12)
  {
    new_state = old_state;
    new_state.dir = dir;
    stress = stress_pred;
    return true;
  }

  Vec7 x{};
  x[0] = Core::FADUtils::cast_to_double(old_state.vM);
  x[1] = 0.0;
  x[2] = (x[0] <= 1.0e-12 ? -1.0e-8 : 0.0);
  x[3] = Core::FADUtils::cast_to_double(old_state.vS);
  x[4] = 0.0;
  x[5] = (x[3] <= 1.0e-12 ? -1.0e-8 : 0.0);
  x[6] = (x[0] + x[3] >= 1.0 - 1.0e-12 ? 1.0e-8 : 0.0);

  auto residual = [&](const Vec7& xin, Vec7& Rout)
  { assemble_moment_residual(xin, old_state, Cur, C_M, dir, Rout, nullptr); };

  if (!solve_nonlinear_system(residual, x)) return false;

  new_state.vM = std::max(0.0, std::min(1.0, x[0]));
  new_state.vS = std::max(0.0, std::min(1.0, x[3]));
  new_state.dir = dir;

  Vec7 Rfinal{};
  assemble_moment_residual(x, old_state, Cur, C_M, dir, Rfinal, &stress);
  return true;
}

template <typename T>
void Mat::BeamSMAMaterial<T>::local_force_response(
    const Vec3& Gamma, const Mat3& C_N, const unsigned int gp, Vec3& stress, Mat3& C_alg)
{
  const LocalState& old_state = force_state_conv_[gp];
  LocalState new_state = old_state;

  if (!solve_force_state(Gamma, C_N, old_state, new_state, stress))
    FOUR_C_THROW("Beam SMA force update failed.");

  force_state_new_[gp] = new_state;

  C_alg.clear();
  C_alg(1, 1) = C_N(1, 1);
  C_alg(2, 2) = C_N(2, 2);

  if (sma_params().get_eps_l_n() <= 0.0)
  {
    C_alg(0, 0) = C_N(0, 0);
    return;
  }

  const double h = 1.0e-8 * std::max(1.0, std::abs(Core::FADUtils::cast_to_double(Gamma(0))));
  Vec3 Gp = Gamma;
  Gp(0) += h;

  LocalState pert_state = old_state;
  Vec3 stress_pert(Core::LinAlg::Initialization::zero);
  if (!solve_force_state(Gp, C_N, old_state, pert_state, stress_pert))
    FOUR_C_THROW("Beam SMA force tangent update failed.");

  C_alg(0, 0) = (stress_pert(0) - stress(0)) / h;
}

template <typename T>
void Mat::BeamSMAMaterial<T>::local_moment_response(
    const Vec3& Cur, const Mat3& C_M, const unsigned int gp, Vec3& stress, Mat3& C_alg)
{
  const LocalState& old_state = moment_state_conv_[gp];
  LocalState new_state = old_state;

  if (!solve_moment_state(Cur, C_M, old_state, new_state, stress))
    FOUR_C_THROW("Beam SMA moment update failed.");

  moment_state_new_[gp] = new_state;

  C_alg.clear();

  if (!sma_params().get_torsion_sma()) C_alg(0, 0) = C_M(0, 0);

  if (sma_params().get_kappa_l_m() <= 0.0)
  {
    C_alg = C_M;
    return;
  }

  const int jbeg = sma_params().get_torsion_sma() ? 0 : 1;
  for (int j = jbeg; j < 3; ++j)
  {
    const double h = 1.0e-8 * std::max(1.0, std::abs(Core::FADUtils::cast_to_double(Cur(j))));
    Vec3 Cp = Cur;
    Cp(j) += h;

    LocalState pert_state = old_state;
    Vec3 stress_pert(Core::LinAlg::Initialization::zero);
    if (!solve_moment_state(Cp, C_M, old_state, pert_state, stress_pert))
      FOUR_C_THROW("Beam SMA moment tangent update failed.");

    for (int i = 0; i < 3; ++i) C_alg(i, j) = (stress_pert(i) - stress(i)) / h;
  }
}

/*-----------------------------------------------------------------------------------------------*/
/* Interface methods                                                                             */
/*-----------------------------------------------------------------------------------------------*/
template <typename T>
void Mat::BeamSMAMaterial<T>::compute_constitutive_parameter(
    Core::LinAlg::Matrix<3, 3, T>& C_N, Core::LinAlg::Matrix<3, 3, T>& C_M)
{
  this->BeamElastHyperMaterial<T>::get_constitutive_matrix_of_forces_material_frame(C_N);
  this->BeamElastHyperMaterial<T>::get_constitutive_matrix_of_moments_material_frame(C_M);

  for (unsigned int gp = 0; gp < numgp_force_; ++gp) c_n_alg_[gp] = C_N;
  for (unsigned int gp = 0; gp < numgp_moment_; ++gp) c_m_alg_[gp] = C_M;
}

template <typename T>
void Mat::BeamSMAMaterial<T>::evaluate_force_contributions_to_stress(
    Core::LinAlg::Matrix<3, 1, T>& stressN, const Core::LinAlg::Matrix<3, 3, T>& C_N,
    const Core::LinAlg::Matrix<3, 1, T>& Gamma, const unsigned int gp)
{
  Mat3 C_alg(Core::LinAlg::Initialization::zero);
  local_force_response(Gamma, C_N, gp, stressN, C_alg);
  c_n_alg_[gp] = C_alg;
}

template <typename T>
void Mat::BeamSMAMaterial<T>::evaluate_moment_contributions_to_stress(
    Core::LinAlg::Matrix<3, 1, T>& stressM, const Core::LinAlg::Matrix<3, 3, T>& C_M,
    const Core::LinAlg::Matrix<3, 1, T>& Cur, const unsigned int gp)
{
  Mat3 C_alg(Core::LinAlg::Initialization::zero);
  local_moment_response(Cur, C_M, gp, stressM, C_alg);
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
void Mat::BeamSMAMaterial<T>::update()
{
  force_state_conv_ = force_state_new_;
  moment_state_conv_ = moment_state_new_;
}

template <typename T>
void Mat::BeamSMAMaterial<T>::reset()
{
  force_state_new_ = force_state_conv_;
  moment_state_new_ = moment_state_conv_;
}

/*-----------------------------------------------------------------------------------------------*/
/* Explicit instantiations                                                                       */
/*-----------------------------------------------------------------------------------------------*/
template class Mat::BeamSMAMaterial<double>;
template class Mat::BeamSMAMaterialType<double>;

FOUR_C_NAMESPACE_CLOSE