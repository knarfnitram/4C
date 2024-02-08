/*----------------------------------------------------------------------*/
/*! \file

\brief  Main Algorithm for partitioned coupled 1D artery and 3d fluid simulations

\level 3

*/
/*----------------------------------------------------------------------*/
#include "baci_artery_cvOneDBF_Partitioned.hpp"

#include "baci_adapter_fld_base_algorithm.H"
#include "baci_fluid_implicit_integration.H"
#include "baci_fluid_utils.H"
#include "baci_lib_discret.H"
#include "baci_lib_globalproblem.H"

BACI_NAMESPACE_OPEN

namespace ARTCV
{


  PartitionAlg::PartitionAlg() {}

  void PartitionAlg::Initialize_Fluid(void)
  {
    // access to some parameter lists
    // const Teuchos::ParameterList& probtype = DRT::Problem::Instance()->ProblemTypeParams();
    const Teuchos::ParameterList& fdyn = DRT::Problem::Instance()->FluidDynamicParams();

    // create instance of fluid basis algorithm
    Teuchos::RCP<ADAPTER::FluidBaseAlgorithm> fluidalgo_ =
        Teuchos::rcp(new ADAPTER::FluidBaseAlgorithm(fdyn, fdyn, "fluid", false));

    // read the restart information, set vectors and variables
    if (DRT::Problem::Instance()->Restart()) dserror("Currently we do not have a propper restart.");
  }

  void PartitionAlg::Initialize_Artery() {}

  void PartitionAlg::Print_Logo()
  {
    const std::string fluid_disname = "fluid";
    const Epetra_Comm& comm = DRT::Problem::Instance()->GetDis(fluid_disname)->Comm();
    if (comm.MyPID() == 0)
    {
      std::cout
          << "---------------------------------------------------------------------------------"
          << std::endl;
      std::cout
          << "-------------------- Welcome to the Partioned artery coupling -------------------"
          << std::endl;
      std::cout
          << "---------------------------------------------------------------------------------"
          << std::endl;
    }
  }

  void PartitionAlg::Post_Process_Fluid()
  {
    std::vector<DRT::Condition*> flowratecond;
    std::string condstring;

    condstring = "SurfFlowRate";
    auto discret_fluid = fluidalgo_->FluidField()->Discretization();
    discret_fluid->GetCondition("SurfFlowRate", flowratecond);

    const Teuchos::RCP<const Epetra_Vector> velnp = fluidalgo_->FluidField()->Velnp();
    if (not velnp.is_valid_ptr())
    {
      dserror("velnp of your fluid problem is not initialized.");
    }
    const Teuchos::RCP<Epetra_Vector> test = Teuchos::RCP(new Epetra_Vector(*velnp));

    //  if no flowrate condition is present we do not compute anything
    if ((int)flowratecond.size() == 0) return;
    auto physicalType = fluidalgo_->FluidField()->PhysicalType();

    const std::map<int, double> flowrates = FLD::UTILS::ComputeFlowRates(
        *discret_fluid, test, Teuchos::null, Teuchos::null, condstring, physicalType);

    const std::map<int, double> meanPressure = FLD::UTILS::ComputeMeanPressure(
        *discret_fluid, test, condstring, INPAR::FLUID::PhysicalType::incompressible);

    std::cout << "flowrates";
    for (const auto& [key, value] : flowrates) std::cout << '[' << key << "] = " << value << "; ";
    std::cout << '\n';

    std::cout << "meanPressure";
    for (const auto& [key, value] : meanPressure)
      std::cout << '[' << key << "] = " << value << "; ";
    std::cout << '\n';
  }


}  // namespace ARTCV

BACI_NAMESPACE_CLOSE