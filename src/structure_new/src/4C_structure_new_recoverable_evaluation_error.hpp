// This file is part of 4C multiphysics licensed under the
// GNU Lesser General Public License v3.0 or later.
//
// See the LICENSE.md file in the top-level for license information.
//
// SPDX-License-Identifier: LGPL-3.0-or-later

#ifndef FOUR_C_STRUCTURE_NEW_RECOVERABLE_EVALUATION_ERROR_HPP
#define FOUR_C_STRUCTURE_NEW_RECOVERABLE_EVALUATION_ERROR_HPP

#include "4C_config.hpp"

#include "4C_utils_exceptions.hpp"

FOUR_C_NAMESPACE_OPEN

namespace Solid
{
  /**
   * @brief Error signaling that the current (trial) state could not be evaluated, but that the
   * problem might be resolved by a smaller time step.
   *
   * Model evaluators can throw this error instead of a hard error (FOUR_C_THROW) if the current
   * state can not be evaluated, e.g., because a Newton iterate is strongly distorted. The implicit
   * structural time integration catches this error and treats the current time step as a failed
   * nonlinear solve, i.e., the action defined by DIVERCONT (e.g., a time step size reduction) is
   * performed.
   *
   * @note When this error leaves a model evaluator, it has to be thrown on all MPI ranks, since the
   * catching code continues with collective operations. Code that might throw this error only on
   * some ranks has to catch it locally and rethrow it on all ranks.
   */
  class RecoverableEvaluationError : public Core::Exception
  {
   public:
    using Core::Exception::Exception;
  };
}  // namespace Solid

FOUR_C_NAMESPACE_CLOSE

#endif
