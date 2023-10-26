/*----------------------------------------------------------------------*/
/*! \file

\brief Implementation of NOX::Group for poroelast

\level 1

*/

/*----------------------------------------------------------------------*/
#include "baci_poroelast_nox_group.H"

#include "baci_poroelast_monolithic.H"
#include "baci_utils_exceptions.H"

/*----------------------------------------------------------------------*
 *----------------------------------------------------------------------*/
NOX::POROELAST::Group::Group(::POROELAST::Monolithic& mporo, Teuchos::ParameterList& printParams,
    const Teuchos::RCP<NOX::Epetra::Interface::Required>& i, const NOX::Epetra::Vector& x,
    const Teuchos::RCP<NOX::Epetra::LinearSystem>& linSys)
    : NOX::Epetra::Group(printParams, i, x, linSys), mporo_(mporo)
{
}


/*----------------------------------------------------------------------*
 *----------------------------------------------------------------------*/
void NOX::POROELAST::Group::CaptureSystemState()
{
  // we know we already have the first linear system calculated
  // mporo_.SetupSystemMatrix();
  // mporo_.SetupRHS(RHSVector.getEpetraVector(), true);
  // mporo_.Extractor()->ExtractVector(*(mporo_.RHS()), 1)->Scale(-1.0);
  // RHSVector.getEpetraVector().Update(1.0, *(mporo_.RHS()), 0.0);
  // mporo_->Extractor()->ExtractVector(F, 1)->Scale(-1.0);
  // std::cout << *(mporo_.RHS()) << std::endl;
  RHSVector.getEpetraVector().Update(-1.0, *(mporo_.RHS()), 0.0);
  sharedLinearSystem.getObject(this);
  isValidJacobian = true;
  isValidRHS = true;
}


/*----------------------------------------------------------------------*
 *----------------------------------------------------------------------*/
NOX::Abstract::Group::ReturnType NOX::POROELAST::Group::computeF()
{
  NOX::Abstract::Group::ReturnType ret = NOX::Epetra::Group::computeF();
  if (ret == NOX::Abstract::Group::Ok)
  {
    if (not isValidJacobian)
    {
      mporo_.SetupSystemMatrix();
      sharedLinearSystem.getObject(this);
      isValidJacobian = true;
    }
  }
  return ret;
}


/*----------------------------------------------------------------------*
 *----------------------------------------------------------------------*/
/*NOX::Abstract::Group::ReturnType NOX::POROELAST::Group::computeJacobian()
{
  NOX::Abstract::Group::ReturnType ret = NOX::Epetra::Group::computeJacobian();
  if (ret == NOX::Abstract::Group::Ok)
  {
    if (not isValidRHS)
    {
      mporo_.SetupRHS(RHSVector.getEpetraVector(), false);
      isValidRHS = true;
    }
  }
  return ret;
}*/


/*----------------------------------------------------------------------*
 *----------------------------------------------------------------------*/
NOX::Abstract::Group::ReturnType NOX::POROELAST::Group::computeNewton(Teuchos::ParameterList& p)
{
  // mporo_.ScaleSystem(RHSVector.getEpetraVector());
  NOX::Abstract::Group::ReturnType status = NOX::Epetra::Group::computeNewton(p);
  // mporo_.UnscaleSolution(NewtonVector.getEpetraVector(), RHSVector.getEpetraVector());

  // check return value of computeNewton call
  if (status == NOX::Abstract::Group::NotConverged || status == NOX::Abstract::Group::Failed)
    dserror("NOX::FSI::Group::computeNewton: linear solver not converged...");

  return status;
}
