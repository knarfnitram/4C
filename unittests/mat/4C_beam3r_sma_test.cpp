// This file is part of 4C multiphysics licensed under the
// GNU Lesser General Public License v3.0 or later.
//
// See the LICENSE.md file in the top-level for license information.
//
// SPDX-License-Identifier: LGPL-3.0-or-later

#include <gtest/gtest.h>

#include "4C_io_input_parameter_container.hpp"
#include "4C_linalg_fixedsizematrix.hpp"
#include "4C_mat_beam3r_sma.hpp"

#include <cmath>
#include <memory>

namespace
{
  using namespace FourC;
  using Params = Mat::PAR::BeamReissnerSMAMaterialParams;
  using Vec3 = Core::LinAlg::Matrix<3, 1>;
  using Mat3 = Core::LinAlg::Matrix<3, 3>;

  //! parameters of a superelastic wire with solid circular cross-section of area A
  std::shared_ptr<Params> make_params(double A, Params::MomentModel moment_model)
  {
    const double I = 0.25 * A * A / M_PI;

    Core::IO::InputParameterContainer container;
    container.add("YOUNG", 53198.6);
    container.add("SHEARMOD", -1.0);
    container.add("POISSONRATIO", 0.215);
    container.add("YOUNGMART", 24666.3);
    container.add("SHEARMODMART", 9016.0);
    container.add("DENS", 6450.0);
    container.add("CROSSAREA", A);
    container.add("SHEARCORR", 1.0);
    container.add("MOMINPOL", 2.0 * I);
    container.add("MOMIN2", I);
    container.add("MOMIN3", I);
    container.add("FAD", false);
    container.add("INTERACTIONRADIUS", std::sqrt(A / M_PI));
    container.add("TEMP", 10.0);
    container.add("DSAMS", 0.31);
    container.add("T0SMA", -48.15);
    container.add("WIN", -0.1);
    container.add("RM", 1.0e20);
    container.add("RSF0", 1.0);
    container.add("RSR0", 17.0);
    container.add("HSF", 1.0);
    container.add("HSR", 1.0);
    container.add("CTS", 0.2);
    container.add("ASF0", 2.0);
    container.add("ASF1", -0.5);
    container.add("ASR0", -0.01);
    container.add("ASR1", 9.0);
    container.add("NEXP", 0.1);
    container.add("EPSLN", 0.0615);
    container.add("MOMENTMODEL", moment_model);
    container.add("NUMFIBERRINGS", 6);
    container.add("NUMFIBERSECTORS", 16);
    container.add("KAPPALM", 0.0);
    container.add("TORSIONSMA", false);
    container.add("RSREG", 1.0e-4);
    container.add("LOCALITER", 50);
    container.add("CONSISTENTTANGENT", true);
    container.add("LOCALSUBSTEP", 10);
    container.add("LOCALSTRAINSTEP", 0.0);
    container.add("LOCALCURVSTEP", 0.0);

    return std::make_shared<Params>(Core::Mat::PAR::Parameter::Data{.parameters = container});
  }

  class BeamSMAMaterialTest : public ::testing::Test
  {
   protected:
    void SetUp() override
    {
      params_ = make_params(1.0, Params::MomentModel::fibers);
      material_ = std::make_shared<Mat::BeamSMAMaterial<double>>(params_.get());
      material_->setup(1, 1);
      material_->compute_constitutive_parameter(C_N_, C_M_);
    }

    Vec3 force(const Vec3& Gamma, Mat3* tangent = nullptr)
    {
      Vec3 N;
      material_->evaluate_force_contributions_to_stress(N, C_N_, Gamma, 0);
      if (tangent) material_->get_stiffness_matrix_of_forces(*tangent, C_N_, 0);
      return N;
    }

    Vec3 moment(const Vec3& Cur, Mat3* tangent = nullptr)
    {
      Vec3 M;
      material_->evaluate_moment_contributions_to_stress(M, C_M_, Cur, 0);
      if (tangent) material_->get_stiffness_matrix_of_moments(*tangent, C_M_, 0);
      return M;
    }

    //! central finite-difference derivative of a stress resultant w.r.t. the strain measure
    template <typename Response>
    static Mat3 finite_difference(Response response, const Vec3& q, double h)
    {
      Mat3 D;
      for (unsigned int j = 0; j < 3; ++j)
      {
        Vec3 qp = q, qm = q;
        qp(j) += h;
        qm(j) -= h;
        const Vec3 sp = response(qp);
        const Vec3 sm = response(qm);
        for (unsigned int i = 0; i < 3; ++i) D(i, j) = (sp(i) - sm(i)) / (2.0 * h);
      }
      return D;
    }

    static void expect_matrix_near(const Mat3& A, const Mat3& B, double relative_tolerance)
    {
      double scale = 0.0;
      for (unsigned int i = 0; i < 3; ++i)
        for (unsigned int j = 0; j < 3; ++j) scale = std::max(scale, std::abs(B(i, j)));
      for (unsigned int i = 0; i < 3; ++i)
        for (unsigned int j = 0; j < 3; ++j)
          EXPECT_NEAR(A(i, j), B(i, j), relative_tolerance * scale)
              << "entry (" << i << "," << j << ")";
    }

    std::shared_ptr<Params> params_;
    std::shared_ptr<Mat::BeamSMAMaterial<double>> material_;
    Mat3 C_N_;
    Mat3 C_M_;
  };

  //! the fibre integration reproduces the elastic section stiffness G J, E I2, E I3 exactly
  TEST_F(BeamSMAMaterialTest, FiberSectionElasticStiffness)
  {
    Vec3 Cur;
    Cur(0) = 1.0e-4;
    Cur(1) = 2.0e-4;
    Cur(2) = -3.0e-4;

    Mat3 tangent;
    const Vec3 M = moment(Cur, &tangent);

    Vec3 M_expected;
    M_expected.multiply(C_M_, Cur);
    for (unsigned int i = 0; i < 3; ++i) EXPECT_NEAR(M(i), M_expected(i), 1.0e-12 * C_M_(0, 0));
    expect_matrix_near(tangent, C_M_, 1.0e-12);
  }

  //! bending is antisymmetric in the curvature
  TEST_F(BeamSMAMaterialTest, FiberSectionAntisymmetry)
  {
    Vec3 Cur;
    Cur(2) = 0.12;
    const Vec3 M_pos = moment(Cur);
    Cur.scale(-1.0);
    const Vec3 M_neg = moment(Cur);

    EXPECT_GT(M_pos(2), 0.0);
    for (unsigned int i = 0; i < 3; ++i) EXPECT_NEAR(M_pos(i), -M_neg(i), 1.0e-10 * M_pos(2));
  }

  //! consistent tangent of the axial state on the forward and on the reverse branch
  TEST_F(BeamSMAMaterialTest, ConsistentTangentForce)
  {
    auto response = [&](const Vec3& q) { return force(q); };

    // forward transformation
    Vec3 Gamma;
    Gamma(0) = 0.03;
    Gamma(1) = 1.0e-3;
    Gamma(2) = -5.0e-4;
    Mat3 tangent;
    force(Gamma, &tangent);
    EXPECT_LT(tangent(0, 0), 0.05 * C_N_(0, 0));  // on the plateau
    expect_matrix_near(tangent, finite_difference(response, Gamma, 1.0e-7), 1.0e-6);

    // reverse transformation after unloading from a converged transformed state
    Gamma(0) = 0.07;
    force(Gamma);
    material_->update();
    Gamma(0) = 0.045;
    force(Gamma, &tangent);
    EXPECT_LT(tangent(0, 0), 0.05 * C_N_(0, 0));
    expect_matrix_near(tangent, finite_difference(response, Gamma, 1.0e-7), 1.0e-6);
  }

  //! consistent tangent of the fibre section under combined bending and torsion
  TEST_F(BeamSMAMaterialTest, ConsistentTangentFiberSection)
  {
    auto response = [&](const Vec3& q) { return moment(q); };

    Vec3 Cur;
    Cur(0) = 0.06;
    Cur(1) = 0.08;
    Cur(2) = -0.05;
    Mat3 tangent;
    moment(Cur, &tangent);
    expect_matrix_near(tangent, finite_difference(response, Cur, 1.0e-8), 1.0e-5);

    // partial unloading: the outer fibres transform back
    material_->update();
    Cur.scale(0.6);
    moment(Cur, &tangent);
    expect_matrix_near(tangent, finite_difference(response, Cur, 1.0e-8), 1.0e-5);
  }

  //! the volume-specific parameters are scaled with the cross-section area for the axial force
  TEST(BeamSMAMaterialScaling, AxialForceScalesWithArea)
  {
    auto params_1 = make_params(1.0, Params::MomentModel::elastic);
    auto params_2 = make_params(0.25, Params::MomentModel::elastic);
    Mat::BeamSMAMaterial<double> material_1(params_1.get());
    Mat::BeamSMAMaterial<double> material_2(params_2.get());
    material_1.setup(1, 1);
    material_2.setup(1, 1);

    Mat3 C_N, C_M;
    material_1.compute_constitutive_parameter(C_N, C_M);
    material_2.compute_constitutive_parameter(C_N, C_M);

    for (const double strain : {0.003, 0.03, 0.075})
    {
      Vec3 Gamma, N_1, N_2;
      Gamma(0) = strain;
      material_1.evaluate_force_contributions_to_stress(N_1, C_N, Gamma, 0);
      material_2.evaluate_force_contributions_to_stress(N_2, C_N, Gamma, 0);
      EXPECT_NEAR(N_2(0), 0.25 * N_1(0), 1.0e-10 * N_1(0)) << "strain " << strain;
    }
  }
}  // namespace
