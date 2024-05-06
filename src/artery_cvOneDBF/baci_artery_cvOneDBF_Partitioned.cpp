/*----------------------------------------------------------------------*/
/*! \file

\brief  Main Algorithm for partitioned coupled 1D artery and 3d fluid simulations

\level 3

*/
/*----------------------------------------------------------------------*/
#include "baci_config.hpp"

#include "baci_artery_cvOneDBF_Partitioned.hpp"

#include "baci_adapter_fld_base_algorithm.hpp"
#include "baci_artery_cvOneDBF_utils.hpp"
#include "baci_fluid_implicit_integration.hpp"
#include "baci_fluid_utils.hpp"
#include "baci_global_data.hpp"
#include "baci_lib_discret.hpp"

// Trilinos
#include <NOX_Epetra_Interface_Required.H>

// svOneD
#include <cvOneDSynchronizer.h>
#include <OneDSolverInterface.h>

BACI_NAMESPACE_OPEN

namespace ARTCV
{


  PartitionAlg::PartitionAlg()
      : AlgorithmBase(GLOBAL::Problem::Instance()->GetDis("fluid")->Comm(),
            GLOBAL::Problem::Instance()->FSIDynamicParams()),
        comm_(GLOBAL::Problem::Instance()->GetDis("fluid")->Comm())
  {
  }

  void PartitionAlg::Initialize_Fluid()
  {
    // set up parameter list for the fluid
    const Teuchos::ParameterList& fdyn = GLOBAL::Problem::Instance()->FluidDynamicParams();

    // create instance of fluid basis algorithm
    fluidalgo_ = Teuchos::rcp(new ADAPTER::FluidBaseAlgorithm(fdyn, fdyn, "fluid", false));

    auto fluid_dis = fluidalgo_->FluidField()->Discretization();

    stepmax_ = fdyn.get<int>("NUMSTEP");

    // read the restart information, set vectors and variables
    if (GLOBAL::Problem::Instance()->Restart())
      dserror("Currently we do not have a propper restart.");
  }

  void PartitionAlg::Set_Neumann_Pressure(const double& pressure)
  {
    // Get discretization
    const Teuchos::RCP<DRT::Discretization>& fluid_dis = fluidalgo_->FluidField()->Discretization();
    for (auto& [name, cond] : fluid_dis->GetAllConditions())
    {
      if (name == (std::string) "LineNeumann" || name == (std::string) "SurfaceNeumann" ||
          name == (std::string) "VolumeNeumann")
      {
        const std::string* type = cond->Get<std::string>("type");
        if (type->compare("neum_pseudo_orthopressure") == 0)
        {
          const std::vector<double> val{pressure, 0.0, 0.0};
          cond->Add("val", val);
        }
      }
    }
  }

  void PartitionAlg::Set_Coupling_Flowrate(const double& flowrate)
  {
    // Get discretization and search for Volumetric Surface Flow conditions
    const Teuchos::RCP<DRT::Discretization>& fluid_dis = fluidalgo_->FluidField()->Discretization();
    for (auto& [name, cond] : fluid_dis->GetAllConditions())
    {
      if (name == (std::string) "VolumetricSurfaceFlowCond")
      {
          // Replace the function value
        const double val = -1.0 * flowrate;
        cond->Add("Val", val);
      }
    }
  }



  void PartitionAlg::Check_Input(void)
  {
    UTILS::executeSerial(comm_,
        [&]
        {
          const Teuchos::ParameterList& fdyn = GLOBAL::Problem::Instance()->FluidDynamicParams();
          fdyn.get<int>("NUMSTEP");
          fdyn.get<double>("TIMESTEP");

          if ((fdyn.get<double>("TIMESTEP") - opts_->timeStep) >
              std::numeric_limits<double>::epsilon())
          {
            dserror("You are using wrong Time step sizes");
          }
          if ((fdyn.get<int>("NUMSTEP") - opts_->maxStep) > std::numeric_limits<int>::epsilon())
          {
            dserror("maxStep and NUMSTEP must be same in the Inputfiles");
          }

          if (fluidalgo_.is_null())
          {
            dserror("fluid algo is empty");
          }
          auto fluid_dis = fluidalgo_->FluidField()->Discretization();

          int number_of_neum_pseudo_orthopressure_conditions = 0;

          for (auto& [name, cond] : fluid_dis->GetAllConditions())
          {
            if (name == (std::string) "LineNeumann" || name == (std::string) "SurfaceNeumann" ||
                name == (std::string) "VolumeNeumann")
            {
              if (cond->Get<std::string>("type")->compare("neum_pseudo_orthopressure"))
              {
                number_of_neum_pseudo_orthopressure_conditions++;
              }
            }
          }

          if (number_of_neum_pseudo_orthopressure_conditions > 1)
          {
            dserror(
                "There can only be 1 neum_pseudo_orthopressure."
                "Sorry for that.");
          }

          std::cout << "Performed all checks" << std::endl;
        });
  }

  void PartitionAlg::Initialize_Artery()
  {
    new_iter_artery = 0;
    const Teuchos::ParameterList& art_params = GLOBAL::Problem::Instance()->Artery_cvOneDParams();


    string inputfile = art_params.get<string>("cvOneD_Inputfile");
    std::cout << "inputfile:" << inputfile << std::endl;


    UTILS::executeSerial(comm_,
        [&]()
        {
          // create model manager
          myOneDSolver_ = Teuchos::rcp(new OneDSolverInterface());

          // Create Solver Options of cvOneD
          opts_ = Teuchos::rcp(new cvOneDOptions());

          // Read Model From File
          myOneDSolver_->readModel(inputfile, opts_.get());

          // perform model check
          opts_->check();

          // 2=char * coupling_type="3D_1D";

          cvOneDSynchronizer_ =
              Teuchos::rcp(new cvOneDSynchronizer(opts_->maxStep + 1, opts_->timeStep, 2, 3));

          myOneDSolver_->setupModeluntilNewton(opts_.get(), cvOneDSynchronizer_.get());
        });
  }

  void PartitionAlg::Artery_Solve()
  {
    UTILS::executeSerial(comm_,
        [&]() {

        });
  }
  void PartitionAlg::Synch_Step(int step)
  {
    UTILS::executeSerial(comm_, [&]() { myOneDSolver_->SynchronizeDataofStep(step); });
  }

  void PartitionAlg::Print_Logo()
  {
    const std::string fluid_disname = "fluid";
    const Epetra_Comm& comm = GLOBAL::Problem::Instance()->GetDis(fluid_disname)->Comm();
    if (comm.MyPID() == 0)
    {
      std::cout
          << "---------------------------------------------------------------------------------"
          << std::endl;
      std::cout
          << "-------------------- Welcome to the Partioned Artery coupling -------------------"
          << std::endl;
      std::cout
          << "---------------------------------------------------------------------------------"
          << std::endl;
    }
  }

  void PartitionAlg::Post_Process_Fluid(double& flowrate, double& pressure, const int cond_id)
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

    std::map<int, double> flowrates = FLD::UTILS::ComputeFlowRates(
        *discret_fluid, test, Teuchos::null, Teuchos::null, condstring, physicalType);

    std::map<int, double> meanPressure = FLD::UTILS::ComputeMeanPressure(
        *discret_fluid, test, condstring, INPAR::FLUID::PhysicalType::incompressible);

    pressure = meanPressure[cond_id];
    flowrate = flowrates[cond_id];

    std::cout << "flowrates";
    for (const auto& [key, value] : flowrates) std::cout << '[' << key << "] = " << value << "; ";
    std::cout << '\n';

    std::cout << "meanPressure";
    for (const auto& [key, value] : meanPressure)
      std::cout << '[' << key << "] = " << value << "; ";
    std::cout << '\n';
  }

  void PartitionAlg::Initialize_Coupling() {}

  void PartitionAlg::Start_Solve(void)
  {
    const int coupling_id_max = 3;
    double p_3d[coupling_id_max] = {0.0, 0.0, 0.0};
    double p_1d[coupling_id_max] = {0.0, 0.0, 0.0};
    double q_3d[coupling_id_max] = {0.0, 0.0, 0.0};
    double q_1d[coupling_id_max] = {0.0, 0.0, 0.0};

    int coupled_iter_max = 10;
    std::vector<int> partition_iteration_count;
    std::vector<double> vec_q_norm;
    std::vector<double> vec_p_norm;
    std::vector<double> vec_p_norm_rel;

    // Assupmtion: condition with 1 equals 1d artery flow into the domain,
    // Conditiond with id 2 equals outflow of domain
    int id_3d_1d = 2;
    int id_1d_3d = 1;
    // const int cond_id = 2;

    // Both time loops start from time_step 1
    for (int time_step = 1; time_step <= stepmax_; time_step++)
    {
      int iter = 0;

      // we take here the next step, because this we do not call initally prep time step.
      const double t_next = fluidalgo_->FluidField()->Time() + fluidalgo_->FluidField()->Dt();
      const double t_prev = fluidalgo_->FluidField()->Time();

      std::cout << "TIMES: " << t_next << " " << t_prev << std::endl;

      UTILS::executeSerial(comm_, [&]() { myOneDSolver_->UpdateTimeStep(); });

      double p_norm_prev[coupling_id_max] = {1e7, 1e7, 1e7};
      double p_norm[coupling_id_max] = {1e7, 1e7, 1e7};
      double q_norm[coupling_id_max] = {1e7, 1e7, 1e7};
      double p_norm_rel[coupling_id_max] = {1e7, 1e7, 1e7};


      while (iter < coupled_iter_max)
      {
        if (comm_.MyPID() == 0)
        {
          new_iter_artery = 0;
          while (true)
          {
            // Newton-Raphson Iterations...
            if (not myOneDSolver_->Do_Newton_Step(&new_iter_artery))
            {
              break;
            }
          }

          std::cout << "done solvin" << std::endl;
          myOneDSolver_->SynchronizeDataofStep(time_step);
          std::cout << " Synchronization completed. get values at t_next:" << t_next << std::endl;
          p_1d[id_3d_1d] = cvOneDSynchronizer_->Get_1d_p_at_t(t_next, id_3d_1d);
          q_1d[id_3d_1d] = cvOneDSynchronizer_->Get_1d_q_at_t(t_next, id_3d_1d);
          p_1d[id_1d_3d] = cvOneDSynchronizer_->Get_1d_p_at_t(t_next, id_1d_3d);
          q_1d[id_1d_3d] = cvOneDSynchronizer_->Get_1d_q_at_t(t_next, id_1d_3d);
          // std::cout << "q_1d:" << q_1d << std::endl;
        }

        // synch pressure to all procs
        comm_.Broadcast(p_1d, coupling_id_max, 0);
        comm_.Broadcast(q_1d, coupling_id_max, 0);

        fluidalgo_->FluidField()->SetTimeStep(t_prev, time_step);

        // call coupling condition for outflow
        // if (id_3d_1d == 2)
        //{
        if (p_1d[id_3d_1d] > 0) {
            std::cout<<"Set_Neumann_Pressure"<<std::endl;
            Set_Neumann_Pressure(p_1d[id_3d_1d]);
        }
        //}


        // call coupling condition for changing inflow
        // if (id_3d_1d == 1)
        //{
        Set_Coupling_Flowrate(q_1d[id_1d_3d]);
        //}

        fluidalgo_->FluidField()->PrepareTimeStep();
        fluidalgo_->FluidField()->Solve();

        Post_Process_Fluid(q_3d[id_3d_1d], p_3d[id_3d_1d], id_3d_1d);
        Post_Process_Fluid(q_3d[id_1d_3d], p_3d[id_1d_3d], id_1d_3d);

        UTILS::executeSerial(comm_,
            [&]()
            {
              cvOneDSynchronizer_->Set_3d_q_at_t(t_next, q_3d[id_3d_1d], id_3d_1d);
              cvOneDSynchronizer_->Set_3d_p_at_t(t_next, p_3d[id_3d_1d], id_3d_1d);
              cvOneDSynchronizer_->Set_3d_q_at_t(t_next, q_3d[id_1d_3d], id_1d_3d);
              cvOneDSynchronizer_->Set_3d_p_at_t(t_next, p_3d[id_1d_3d], id_1d_3d);
            });
        Synch_Step(time_step);
          for (int i = 0; i < coupling_id_max; ++i) {
              p_norm[i] = (p_3d[i] - p_1d[i]) * (p_3d[i] - p_1d[i]);
              q_norm[i] = (abs(q_3d[i]) - abs(q_1d[i])) *
                                 (abs(q_3d[i]) - abs(q_1d[i]));
          }


        /*std::cout << "iter: " << iter << " p_1d: " << p_1d << " p_3d " << p_3d << std::endl;
        std::cout << "iter: " << iter << " q_1d: " << q_1d << " q_3d " << q_3d << std::endl;
        std::cout << "p_norm: " << p_norm << " q_norm: " << q_norm << std::endl;*/

        UTILS::executeSerial(comm_, [&]() { cvOneDSynchronizer_->Print(); });

        // TODO this needs to be done for every condition
        if (p_norm[0] < 1e-4 and q_norm[0] < 1e-4 and p_norm[1] < 1e-4 and q_norm[1] < 1e-4 and iter)
        {
          std::cout << "p_norm and q_norm converged" << std::endl;
          break;
        }

        if (abs(p_norm_prev[0] - p_norm[0]) < 1e-6 and abs(p_norm_prev[1] - p_norm[1]) < 1e-6 and iter > 0 and q_norm[0] < 1e-10 and q_norm[1] < 1e-10)
        {
          std::cout << "p_rel_prev stayed same... breaking" << std::endl;
          break;
        }
          for (int i = 0; i < coupling_id_max; ++i) {
              p_norm_prev[i] = p_norm[i];
          }


        iter++;
      }

      partition_iteration_count.push_back(iter);
      vec_q_norm.push_back(q_norm[1]);
      vec_p_norm.push_back(p_norm[1]);
      vec_p_norm_rel.push_back(p_norm[1] / p_3d[1]);

      fluidalgo_->FluidField()->Update();
      fluidalgo_->FluidField()->StatisticsAndOutput();

      comm_.Barrier();
      UTILS::executeSerial(
          comm_, [&]() { myOneDSolver_->UpdateSolution(new_iter_artery, time_step_artery); });
    }

    if (comm_.MyPID() == 0)
    {
      std::cout << "Partitioned iteration count: " << std::endl;
      for (int n : partition_iteration_count) std::cout << n << ' ';
      std::cout << '\n';
      std::cout << "vec_q_norm: " << std::endl;
      for (double n : vec_q_norm) std::cout << n << ' ';
      std::cout << '\n';
      std::cout << "vec_p_norm: " << std::endl;
      for (double n : vec_p_norm) std::cout << n << ' ';
      std::cout << '\n';
      std::cout << "vec_p_norm_rel " << std::endl;
      for (double n : vec_p_norm_rel) std::cout << n << ' ';
      std::cout << '\n';
    }
    // Post Process Solution of Artery
    UTILS::executeSerial(comm_, [&]() { myOneDSolver_->DoPostProcessing(); });
  }

  void PartitionAlg::ReadRestart(int step) { dserror("Restart needs to be implemented!."); }

  bool PartitionAlg::computeF(const Epetra_Vector& x, Epetra_Vector& F, const FillType fillFlag)
  {
    return false;
  }

  void PartitionAlg::Perform_Baci_tests(void)
  {
    GLOBAL::Problem::Instance()->AddFieldTest(fluidalgo_->FluidField()->CreateFieldTest());
    GLOBAL::Problem::Instance()->TestAll(comm_);
    return;
  }


}  // namespace ARTCV

BACI_NAMESPACE_CLOSE