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

#ifndef FOUR_C_ARTERY_CVONEDBF_PARTITIONED_HPP
#define FOUR_C_ARTERY_CVONEDBF_PARTITIONED_HPP

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

    //! Update the pseudo Neumann pressure of condition with ID
    void Set_Neumann_Pressure(const double& pressure, const int ID);

    //! Update the dirichlet flowrate of condition with ID
    void Set_Coupling_Flowrate(const double& flowrate, const int ID);


    //! Output the Logo of the problem
    void Print_Logo(void);

    //! Perform checks to ensure that the problems are set up correctly
    void Check_Input(void);

    //! Solve artery problem;
    void Artery_Solve(void);

    //! evaluate all quantities needed from 3D fluid after Newton
    void Post_Process_Fluid();

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

    // list containing all the coupling ids
    std::vector<int> id_list;

    // list of the maximum coupling ids
    int coupling_id_max;

    // maximum coupling iterations for partitioned solve
    const int coupled_iter_max;

    // different values for pressure
    std::vector<double> p_3d;
    std::vector<double> p_1d;

    // different values for flowrate
    std::vector<double> q_3d;
    std::vector<double> q_1d;

    // stopping tolerances
    const double tol_q;
    const double tol_p;

    // tolerances compared to last iteration
    const double tol_q_c;
    const double tol_p_c;
  };

}  // namespace ARTCV

BACI_NAMESPACE_CLOSE

#endif
