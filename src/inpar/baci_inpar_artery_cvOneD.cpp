/*----------------------------------------------------------------------*/
/*! \file
\brief Input parameters for cvOneD arteries coupled with Fluid-Problem


\level 3
*/


#include "baci_inpar_artery_cvOneD.hpp"

#include "baci_inpar.hpp"
#include "baci_inpar_validparameters.hpp"

BACI_NAMESPACE_OPEN


void INPAR::ARTCV::SetValidParameters(Teuchos::RCP<Teuchos::ParameterList> list)
{
  using namespace INPUT;
  using Teuchos::setStringToIntegralParameter;
  using Teuchos::tuple;

  Teuchos::ParameterList& arteryCV = list->sublist("FLUID ARTERY CVONED COUPLING", false,
      "Parameter list for the FLUID coupled with CVONED ARTERY Problem");

  StringParameter("cvOneD_Inputfile", "", "Absolute Path to cvOneD input file", &arteryCV);

  DoubleParameter("TOL_COUPLE_Q", 1e-4, "Tolerance for convergence check of flow rate", &arteryCV);
  DoubleParameter("TOL_COUPLE_P", 1e-4, "Tolerance for convergence check of pressure", &arteryCV);

  IntParameter("Couple_Iter", 5, "Numer of Partitioned Iterations", &arteryCV);

  // Teuchos::ParameterList &beaminteraction = list->sublist("", false, "");
}

BACI_NAMESPACE_CLOSE