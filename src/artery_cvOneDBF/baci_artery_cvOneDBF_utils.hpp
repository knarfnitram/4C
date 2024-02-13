/*----------------------------------------------------------------------*/
/*! \file

\brief  Utils functions partitioned coupled 1D artery and 3d fluid simulations

\level 3

*/
/*----------------------------------------------------------------------*/
#include "baci_config.hpp"

#include <Epetra_Comm.h>

#ifndef BACI_ARTERY_CVONEDBF_UTILS_HPP
#define BACI_ARTERY_CVONEDBF_UTILS_HPP


BACI_NAMESPACE_OPEN

namespace ARTCV::UTILS
{

  // function which ensures that a code gets only executed on 1 proc.
  void executeSerial(const Epetra_Comm& comm, const std::function<void()>& code);

}  // namespace ARTCV::UTILS

BACI_NAMESPACE_CLOSE

#endif