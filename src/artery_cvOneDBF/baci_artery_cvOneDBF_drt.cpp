/*----------------------------------------------------------------------*/
/*! \file

\brief  Main control routine for partitioned coupled 1D artery and 3d fluid simulations

\level 3

*/
/*----------------------------------------------------------------------*/

#include "baci_artery_cvOneDBF_drt.h"

#include "baci_adapter_fld_base_algorithm.H"
#include "baci_fluid_implicit_integration.H"
#include "baci_fluid_utils.H"
#include "baci_lib_discret.H"
#include "baci_lib_globalproblem.H"

BACI_NAMESPACE_OPEN

void artery_cvOneDBF_drt()
{
  // access to some parameter lists
  // const Teuchos::ParameterList& probtype = DRT::Problem::Instance()->ProblemTypeParams();
  const Teuchos::ParameterList& fdyn = DRT::Problem::Instance()->FluidDynamicParams();

  // create instance of fluid basis algorithm
  Teuchos::RCP<ADAPTER::FluidBaseAlgorithm> fluidalgo =
      Teuchos::rcp(new ADAPTER::FluidBaseAlgorithm(fdyn, fdyn, "fluid", false));

  // read the restart information, set vectors and variables
  if (DRT::Problem::Instance()->Restart()) dserror("Currently we do not have a propper restart.");

  auto comm = fluidalgo->FluidField()->Discretization()->Comm().MyPID();
  // const Epetra_Comm& comm = DRT::Problem::Instance()->GetDis("fluid")->Comm();

  if (comm == 0)
  {
    std::cout << "---------------------------------------------------------------------------------"
              << std::endl;
    std::cout << "-------------------- Welcome to the Partioned artery coupling -------------------"
              << std::endl;
    std::cout << "---------------------------------------------------------------------------------"
              << std::endl;
  }

  if (comm == 0)
  {
  }



  std::vector<DRT::Condition*> flowratecond;
  std::string condstring;
  condstring = "SurfFlowRate";
  auto discret = fluidalgo->FluidField()->Discretization();
  discret->GetCondition("SurfFlowRate", flowratecond);

  const Teuchos::RCP<const Epetra_Vector> velnp = fluidalgo->FluidField()->Velnp();
  if (not velnp.is_valid_ptr())
  {
    dserror("velnp of your fluid problem is not initialized.");
  }
  const Teuchos::RCP<Epetra_Vector> test = Teuchos::RCP(new Epetra_Vector(*velnp));
  // test->Update(1.0,*velnp,0.0);
  //  if no flowrate condition is present we do not compute anything
  if ((int)flowratecond.size() == 0) return;
  auto physicalType = fluidalgo->FluidField()->PhysicalType();

  const std::map<int, double> flowrates = FLD::UTILS::ComputeFlowRates(
      *discret, test, Teuchos::null, Teuchos::null, condstring, physicalType);

  const std::map<int, double> meanPressure = FLD::UTILS::ComputeMeanPressure(
      *discret, test, condstring, INPAR::FLUID::PhysicalType::incompressible);

  std::cout << "flowrates";
  for (const auto& [key, value] : flowrates) std::cout << '[' << key << "] = " << value << "; ";
  std::cout << '\n';

  std::cout << "meanPressure";
  for (const auto& [key, value] : meanPressure) std::cout << '[' << key << "] = " << value << "; ";
  std::cout << '\n';

  // run the fluid simulation


  fluidalgo->FluidField()->PrepareTimeStep();

  fluidalgo->FluidField()->Solve();

  fluidalgo->FluidField()->StatisticsAndOutput();



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
