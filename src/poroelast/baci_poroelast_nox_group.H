/*----------------------------------------------------------------------*/
/*! \file

\brief Implementation of NOX::Group for Poroelast

\level 1

*/

/*----------------------------------------------------------------------*/

#ifndef BACI_POROELAST_NOX_GROUP_H
#define BACI_POROELAST_NOX_GROUP_H

#include <NOX_Epetra_Group.H>

// forward declarations
namespace POROELAST
{
  class Monolithic;
}

namespace NOX
{
  namespace POROELAST
  {
    /// Special NOX group that always sets Jacobian and RHS at the same time.
    class Group : public NOX::Epetra::Group
    {
     public:
      Group(::POROELAST::Monolithic& mporo,     ///< monolithic POROELAST interface
          Teuchos::ParameterList& printParams,  ///< printing parameters
          const Teuchos::RCP<NOX::Epetra::Interface::Required>& i,  ///< NOX interface
          const NOX::Epetra::Vector& x,                             ///< initial guess
          const Teuchos::RCP<NOX::Epetra::LinearSystem>& linSys     ///< linear system
      );

      /// fetch the known Jacobian and RHS from the field solvers
      void CaptureSystemState();

      /// calculate the RHS vector
      virtual NOX::Abstract::Group::ReturnType computeF();

      /// calculate the Jacobian matrix
      virtual NOX::Abstract::Group::ReturnType computeJacobian();

      NOX::Abstract::Group::ReturnType computeNewton(Teuchos::ParameterList& p);

     private:
      ::POROELAST::Monolithic& mporo_;
    };
  }  // namespace POROELAST
}  // namespace NOX

#endif
