/*----------------------------------------------------------------------*/
/*! \file

\brief  Main Algorithm for partitioned coupled 1D artery and 3d fluid simulations

\level 3

*/
/*----------------------------------------------------------------------*/
#include <Teuchos_RCPDecl.hpp>

#ifndef BACI_PARTITIONALG_H
#define BACI_PARTITIONALG_H

namespace BACI
{

  namespace ADAPTER
  {
    class FluidBaseAlgorithm;
  }

  namespace ARTCV
  {

    class PartitionAlg
    {
     public:
      //! Constructor();
      PartitionAlg();

      //! Initialize the fluid field
      void Initialize_Fluid(void);

      //! Initialize the artery field
      void Initialize_Artery(void);

      //! Output the Logo of the problem
      void Print_Logo(void);

      void Check_Input(void);

      //! evaluate all quantities needed form 3D fluid
      void Post_Process_Fluid(void);

      //! evaluate all quantities needed from 1D artery
      void Post_Process_Artery(void);

     private:
      Teuchos::RCP<ADAPTER::FluidBaseAlgorithm> fluidalgo_;
    };

  }  // namespace ARTCV

#endif  // BACI_PARTITION_H
}