/*----------------------------------------------------------------------*/
/*! \file

\brief  Utils functions partitioned coupled 1D artery and 3d fluid simulations

\level 3

*/
/*----------------------------------------------------------------------*/

#include "baci_config.hpp"

#include "baci_artery_cvOneDBF_utils.hpp"

#include <Epetra_Comm.h>

BACI_NAMESPACE_OPEN

namespace ARTCV::UTILS
{
  void executeSerial(const Epetra_Comm& comm, const std::function<void()>& code)
  {
    if (comm.MyPID() == 0)
    {
      code();
    }
  }
}  // namespace ARTCV::UTILS

BACI_NAMESPACE_CLOSE