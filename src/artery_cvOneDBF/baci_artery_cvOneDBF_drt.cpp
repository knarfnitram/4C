/*----------------------------------------------------------------------*/
/*! \file

\brief  Main control routine for partitioned coupled 1D artery and 3d fluid simulations

\level 3

*/
/*----------------------------------------------------------------------*/

#include "baci_config.hpp"

#include "baci_artery_cvOneDBF_drt.hpp"

#include "baci_artery_cvOneDBF_Partitioned.hpp"
#include "baci_global_data.hpp"

BACI_NAMESPACE_OPEN

void artery_cvOneDBF_drt()
{
  ARTCV::PartitionAlg partition_solver = ARTCV::PartitionAlg();

  partition_solver.Initialize_Fluid();
  partition_solver.Initialize_Artery();

  partition_solver.Check_Input();

  // run the fluid simulation
  partition_solver.Start_Solve();

  //
  partition_solver.Perform_Baci_tests();
}

BACI_NAMESPACE_CLOSE
