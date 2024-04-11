/*-----------------------------------------------------------*/
/*! \file

\brief input parameters for cvOneD arteries coupled with Fluid-Problem

\level 3

*/
/*-----------------------------------------------------------*/

#ifndef FOUR_C_INPAR_ARTERY_CVONED_HPP
#define FOUR_C_INPAR_ARTERY_CVONED_HPP

#include "baci_config.hpp"

#include "baci_utils_parameter_list.hpp"
BACI_NAMESPACE_OPEN
// Forward declaration.
namespace INPUT
{
  class ConditionDefinition;
}

namespace INPAR
{
  namespace ARTCV
  {

    /// set the ARTCV parameters
    void SetValidParameters(Teuchos::RCP<Teuchos::ParameterList> list);

    /// set specific ARTCV conditions
    void SetValidConditions(std::vector<Teuchos::RCP<INPUT::ConditionDefinition>>& condlist);

  }  // namespace ARTCV

}  // namespace INPAR

/*----------------------------------------------------------------------*/
BACI_NAMESPACE_CLOSE

#endif
