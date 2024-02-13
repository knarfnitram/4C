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
  partition_solver.Initialize_Coupling();



  // run the fluid simulation


  //  fluidalgo->FluidField()->PrepareTimeStep();

  //  fluidalgo->FluidField()->Solve();

  //  fluidalgo->FluidField()->StatisticsAndOutput();



  // Setup the fluid simulation
  // Setup the cvOneDBF simulation

  // initialize a Newton Method for partioned coupling

  // fluidalgo->FluidField()->Integrate();
  //    fluidalgo->FluidField()->TimeLoop();

  // perform result tests if required
  // DRT::Problem::Instance()->AddFieldTest(fluidalgo->FluidField()->CreateFieldTest());
  // DRT::Problem::Instance()->TestAll(comm);

  // start the Newton iteration
}

BACI_NAMESPACE_CLOSE
