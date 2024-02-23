/*----------------------------------------------------------------------*/
/*! \file

\brief  Main Algorithm for partitioned coupled 1D artery and 3d fluid simulations

\level 3

*/
/*----------------------------------------------------------------------*/
#include "baci_config.hpp"

#include "baci_adapter_algorithmbase.hpp"

#include <NOX.H>
#include <NOX_Epetra_Interface_Required.H>
#include <Teuchos_RCPDecl.hpp>

#ifndef BACI_ARTERY_CVONEDBF_PARTITIONED_HPP
#define BACI_ARTERY_CVONEDBF_PARTITIONED_HPP

class Epetra_Vector;
class FillType;
class OneDSolverInterface;
class cvOneDSynchronizer;
class cvOneDOptions;

BACI_NAMESPACE_OPEN

namespace ADAPTER
{
  class FluidBaseAlgorithm;
}

namespace ARTCV
{

  class PartitionAlg : public ADAPTER::AlgorithmBase, public ::NOX::Epetra::Interface::Required
  {
   public:
    //! standard Constructor();
    PartitionAlg();

    //! Restart
    void ReadRestart(int step) override;

    //! Initialize the fluid field
    void Initialize_Fluid(void);

    //! Initialize the artery field
    void Initialize_Artery(void);

    //! Initialize the artery field
    void Initialize_Coupling(void);

    //! Synchronize 1D solver at step
    void Synch_Step(int step);

    //! Set the Pressure on the 3D Domain
    void Set_Neumann_Pressure(const double& pressure);

    //! Output the Logo of the problem
    void Print_Logo(void);

    //! Perform checks to ensure that the problems are set up correctly
    void Check_Input(void);

    //! Solve artery problem;
    void Artery_Solve(void);

    //! evaluate all quantities needed from 3D fluid after Newton
    void Post_Process_Fluid(double& flowrate, double& pressure);

    //! evaluate all quantities needed from 1D artery after Newton
    void Post_Process_Artery(void);

    void Check_Convergence(void);

    void Start_Solve(void);

    void Timeloop(const Teuchos::RCP<NOX::Epetra::Interface::Required>& interface);

    /// compute FSI interface residual S^{-1}(F(d)) - d
    bool computeF(const Epetra_Vector& x, Epetra_Vector& F, const FillType fillFlag) override;

    void Perform_Baci_tests(void);

   private:
    //! 3D fluid solver
    Teuchos::RCP<ADAPTER::FluidBaseAlgorithm> fluidalgo_;

    //! 1D artery solver
    Teuchos::RCP<OneDSolverInterface> myOneDSolver_;

    //! 1D artery Synchronizer to exchange Data
    Teuchos::RCP<cvOneDSynchronizer> cvOneDSynchronizer_;

    //! Options of the 1D artery solver
    Teuchos::RCP<cvOneDOptions> opts_;

    //! Number of Newton iterations of artery - needed for Postprocessing the Solution
    int new_iter_artery;

    //! Number of Newton iterations of artery - needed for Postprocessing the Solution
    int time_step_artery;

    int stepmax_;

    // communicator
    const Epetra_Comm& comm_;
  };

}  // namespace ARTCV

BACI_NAMESPACE_CLOSE

#endif
