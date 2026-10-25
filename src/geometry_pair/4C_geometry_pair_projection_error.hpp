// This file is part of 4C multiphysics licensed under the
// GNU Lesser General Public License v3.0 or later.
//
// See the LICENSE.md file in the top-level for license information.
//
// SPDX-License-Identifier: LGPL-3.0-or-later

#ifndef FOUR_C_GEOMETRY_PAIR_PROJECTION_ERROR_HPP
#define FOUR_C_GEOMETRY_PAIR_PROJECTION_ERROR_HPP

#include "4C_config.hpp"

#include "4C_utils_exceptions.hpp"

FOUR_C_NAMESPACE_OPEN

namespace GeometryPair
{
  /**
   * @brief Error signaling that the Gauss points of a segment could not be projected for the
   * current configuration, or that the segments of a line overlap.
   *
   * It is thrown instead of a hard error (FOUR_C_THROW), if
   * GEOMETRY_PAIR_SEGMENTATION_PROJECTION_FAILURE_ACTION is set to reduce_time_step. The code
   * evaluating the pairs can catch it and treat the current state as not evaluable, e.g., to reduce
   * the time step size. Such failures typically occur for strongly distorted Newton iterates.
   */
  class ProjectionError : public Core::Exception
  {
   public:
    using Core::Exception::Exception;
  };
}  // namespace GeometryPair

FOUR_C_NAMESPACE_CLOSE

#endif
