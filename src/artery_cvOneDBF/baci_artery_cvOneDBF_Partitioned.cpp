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
        comm_(GLOBAL::Problem::Instance()->GetDis("fluid")->Comm()),
        coupled_iter_max(
            GLOBAL::Problem::Instance()->Artery_cvOneDParams().get<int>("COUPLE_ITER")),
        tol_q(GLOBAL::Problem::Instance()->Artery_cvOneDParams().get<double>("TOL_COUPLE_Q")),
        tol_p(GLOBAL::Problem::Instance()->Artery_cvOneDParams().get<double>("TOL_COUPLE_P")),
        tol_q_c(
            GLOBAL::Problem::Instance()->Artery_cvOneDParams().get<double>("TOL_COUPLE_CHANGE_Q")),
        tol_p_c(
            GLOBAL::Problem::Instance()->Artery_cvOneDParams().get<double>("TOL_COUPLE_CHANGE_P")),
        tol_q_rel(
            GLOBAL::Problem::Instance()->Artery_cvOneDParams().get<double>("TOL_COUPLE_Q_REL")),
        tol_p_rel(
            GLOBAL::Problem::Instance()->Artery_cvOneDParams().get<double>("TOL_COUPLE_P_REL"))
  {
  }

  void PartitionAlg::Initialize_Fluid()
  {
    // set up parameter list for the fluid
    const Teuchos::ParameterList& fdyn = GLOBAL::Problem::Instance()->FluidDynamicParams();

    // Create Coupling Pressure and Inflow Conditions according
    Teuchos::RCP<DRT::Discretization> actdis = GLOBAL::Problem::Instance()->GetDis("fluid");
    std::vector<DRT::Condition*> neumann_coupling_conditions;
    actdis->GetCondition("SurfaceCoupling1DArteryNeumannPressure", neumann_coupling_conditions);

    // The conditions we create here have an additional parameter called "ID"
    // with this they can be distinguished from "normal" Neumann or Inflowrate conditions
    for (DRT::Condition* ncc : neumann_coupling_conditions)
    {
      // Add for every neumann pressure artery coupling a new surface condition
      id_list.push_back(*ncc->Get<int>("ID"));
      auto new_cond = Teuchos::rcp(new DRT::Condition(*ncc, DRT::Condition::SurfaceNeumann));
      new_cond->Add("type", string("neum_pseudo_orthopressure"));

      actdis->SetCondition("SurfaceNeumann", new_cond);

      // Add for every neumann pressure the surface for evaluation
      auto new_flowrate_cond =
          Teuchos::rcp(new DRT::Condition(*ncc, DRT::Condition::FlowRateThroughSurface_3D));
      new_flowrate_cond->Add("ConditionID", int(*ncc->Get<int>("ID")));
      actdis->SetCondition("SurfFlowRate", new_flowrate_cond);
    }

    std::vector<DRT::Condition*> inflow_coupling_conditions;
    actdis->GetCondition("SurfaceCoupling1DArteryDirichletFlow", inflow_coupling_conditions);
    if (inflow_coupling_conditions.empty() and neumann_coupling_conditions.empty())
    {
      dserror(
          "I have not found a SurfaceCoupling1DArteryNeumannPressure or "
          "SurfaceCoupling1DArteryDirichletFlow Condition. Please provide one of those coupling "
          "conditions.");
    }

    // Add for every inflow artery coupling a new surface condition
    for (DRT::Condition* icc : inflow_coupling_conditions)
    {
      id_list.push_back(*icc->Get<int>("ID"));
      auto new_cond =
          Teuchos::rcp(new DRT::Condition(*icc, DRT::Condition::VolumetricSurfaceFlowCond));

      // Add Additional ConditionID, which is the same as the normal coupling ID
      new_cond->Add("ConditionID", *icc->Get<int>("ID"));
      actdis->SetCondition("VolumetricSurfaceFlowCond", new_cond);

      // Add flowrate and mean pressure conditions for evaluation
      auto new_flowrate_cond =
          Teuchos::rcp(new DRT::Condition(*icc, DRT::Condition::FlowRateThroughSurface_3D));
      new_flowrate_cond->Add("ConditionID", *icc->Get<int>("ID"));
      actdis->SetCondition("SurfFlowRate", new_flowrate_cond);
    }

    // create instance of fluid basis algorithm
    fluidalgo_ = Teuchos::rcp(new ADAPTER::FluidBaseAlgorithm(fdyn, fdyn, "fluid", false));

    auto fluid_dis = fluidalgo_->FluidField()->Discretization();

    stepmax_ = fdyn.get<int>("NUMSTEP");

    // read the restart information, set vectors and variables
    if (GLOBAL::Problem::Instance()->Restart()) dserror("Currently we do not support restart.");

    // check for unique ID by sorting
    if (!std::is_sorted(id_list.begin(), id_list.end())) std::sort(id_list.begin(), id_list.end());

    auto adjacent_element =
        std::adjacent_find(id_list.begin(), id_list.end(), std::not_equal_to<int>());

    if (adjacent_element == id_list.end() and id_list.size() > 1)
    {
      for (int id : id_list)
      {
        std::cout << id << std::endl;
      }
      dserror(
          "It seems that some Coupling IDs have not a unique numbering. Please provide unique IDs "
          "on Coupling Conditions.");
    }

    // we allocate the size of the highest id + 1 (place at 1 remains unused)
    coupling_id_max = id_list.back() + 1;

    // resize the arrays accordingly and intitalize

    p_3d.resize(coupling_id_max, 0.0);
    p_1d.resize(coupling_id_max, 0.0);
    q_3d.resize(coupling_id_max, 0.0);
    q_1d.resize(coupling_id_max, 0.0);
  }

  void PartitionAlg::Set_Neumann_Pressure(const double& pressure, const int ID)
  {
    if (pressure > 0)
    {
      // Get discretization
      const Teuchos::RCP<DRT::Discretization>& fluid_dis =
          fluidalgo_->FluidField()->Discretization();
      for (auto& [name, cond] : fluid_dis->GetAllConditions())
      {
        if (name == (std::string) "SurfaceNeumann" and ((cond->GetIf<int>("ID")) != nullptr))
        {
          if (*(cond->Get<int>("ID")) == ID)
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
    }
  }

  void PartitionAlg::Set_Coupling_Flowrate(const double& flowrate, const int ID)
  {
    // Get discretization and search for Volumetric Surface Flow conditions
    const Teuchos::RCP<DRT::Discretization>& fluid_dis = fluidalgo_->FluidField()->Discretization();
    for (auto& [name, cond] : fluid_dis->GetAllConditions())
    {
      // check if we have a flow condition which was created due to the coupled problem
      if (name == (std::string) "VolumetricSurfaceFlowCond" and (cond->GetIf<int>("ID") != nullptr))
      {
        // we for certainly have found a Condition, which needs to be udpated
        // check if we have the ID generated by originally Coupling Problem
        if (*cond->Get<int>("ID") == ID)
        {
          // Replace the function value
          const double val = -1.0 * flowrate;
          cond->Add("Val", val);
        }
      }
    }
  }



  void PartitionAlg::Check_Input(void)
  {
    UTILS::executeSerial(comm_,
        [&]
        {
          const Teuchos::ParameterList& fdyn = GLOBAL::Problem::Instance()->FluidDynamicParams();

          // TODO this check seems not to work...
          if (std::abs((fdyn.get<double>("TIMESTEP") - opts_->timeStep)) >
              std::numeric_limits<double>::epsilon())
          {
            dserror("Currently we are assuming, that the time step size must be same.");
          }

          if (std::abs((fdyn.get<int>("NUMSTEP") - opts_->maxStep)) <
              std::numeric_limits<int>::epsilon())
          {
            dserror("maxStep and NUMSTEP must be same in the Inputfiles");
          }

          if (fluidalgo_.is_null())
          {
            dserror("fluid algo is empty");
          }
          auto fluid_dis = fluidalgo_->FluidField()->Discretization();

          int number_of_neum_pseudo_orthopressure_conditions = 0;

          std::vector<int> flowrate_ids;

          for (auto& [name, cond] : fluid_dis->GetAllConditions())
          {
            cond->Print(cout);
            if (name == (std::string) "LineNeumann" || name == (std::string) "SurfaceNeumann" ||
                name == (std::string) "VolumeNeumann")
            {
              if (cond->Get<std::string>("type")->compare("neum_pseudo_orthopressure"))
              {
                number_of_neum_pseudo_orthopressure_conditions++;
              }
            }
            if (name == (std::string) "SurfFlowRate" and (cond->GetIf<int>("ID") != nullptr))
            {
              flowrate_ids.push_back(*cond->Get<int>("ID"));
            }
          }

          // sort flow rate ids
          std::sort(flowrate_ids.begin(), flowrate_ids.end());

          if (std::adjacent_find(flowrate_ids.begin(), flowrate_ids.end()) != flowrate_ids.end())
          {
            dserror(
                " I found multiple Surf Flow Rate Conditions. Please make sure that they are "
                "distinct. The First n- Coupling Conditions IDs are reserved. ");
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

          cvOneDSynchronizer_ = Teuchos::rcp(
              new cvOneDSynchronizer(opts_->maxStep + 1, opts_->timeStep, coupling_id_max - 1, 3));

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

    std::map<int, double> flowrates = FLD::UTILS::ComputeFlowRates(
        *discret_fluid, test, Teuchos::null, Teuchos::null, condstring, physicalType);

    std::map<int, double> meanPressure = FLD::UTILS::ComputeMeanPressure(
        *discret_fluid, test, condstring, INPAR::FLUID::PhysicalType::incompressible);

    // store data back, since flowrate condition is set up based
    // the ids
    for (int id : id_list)
    {
      p_3d[id] = meanPressure[id];
      q_3d[id] = flowrates[id];
    }


    /*std::cout << "flowrates";
    for (const auto& [key, value] : flowrates) std::cout << '[' << key << "] = " << value << "; ";
    std::cout << '\n';

    std::cout << "meanPressure";
    for (const auto& [key, value] : meanPressure)
      std::cout << '[' << key << "] = " << value << "; ";
    std::cout << '\n';*/
  }

  void PartitionAlg::Initialize_Coupling() {}

  void PartitionAlg::Start_Solve(void)
  {
    // iteration count
    std::vector<int> partition_iteration_count;

    // final norms
    std::vector<double> vec_q_norm;
    std::vector<double> vec_p_norm;
    std::vector<double> vec_p_norm_rel;

    // norms for conditions
    std::vector<double> p_norm;
    std::vector<double> q_norm;

    // norms of last iteration
    std::vector<double> p_norm_prev;
    std::vector<double> q_norm_prev;

    // indicate if condition converged
    std::vector<bool> converged_condition;
    std::vector<bool> norms_stayed_same;

    // Both time loops start from time_step 1
    for (int time_step = 1; time_step <= stepmax_; time_step++)
    {
      int iter = 0;

      // initialize norms with values
      p_norm.resize(coupling_id_max, 1e7);
      q_norm.resize(coupling_id_max, 1e7);

      p_norm_prev.resize(coupling_id_max, 1e6);
      q_norm_prev.resize(coupling_id_max, 1e6);



      // we take here the next step, because this we do not call initally prep time step.
      const double t_next = fluidalgo_->FluidField()->Time() + fluidalgo_->FluidField()->Dt();
      const double t_prev = fluidalgo_->FluidField()->Time();

      UTILS::executeSerial(comm_, [&]() { myOneDSolver_->UpdateTimeStep(); });

      while (iter < coupled_iter_max)
      {
        // reset the convergence arrays
        converged_condition.resize(coupling_id_max, false);
        norms_stayed_same.resize(coupling_id_max, false);

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

          myOneDSolver_->SynchronizeDataofStep(time_step);
          for (int id = 1; id < coupling_id_max; ++id)
          {
            p_1d[id] = cvOneDSynchronizer_->Get_1d_p_at_t(t_next, id);
            q_1d[id] = cvOneDSynchronizer_->Get_1d_q_at_t(t_next, id);
          }
        }


        //  synch pressure to all other procs
        comm_.Broadcast(p_1d.data(), coupling_id_max, 0);
        comm_.Broadcast(q_1d.data(), coupling_id_max, 0);

        fluidalgo_->FluidField()->SetTimeStep(t_prev, time_step);

        // Loop through conditions and either set the neuman or flowrate according to condtion
        for (int id = 1; id < coupling_id_max; ++id)
        {
          Set_Neumann_Pressure(p_1d[id], id);
          Set_Coupling_Flowrate(q_1d[id], id);
        }

        const Teuchos::RCP<DRT::Discretization>& fluid_dis =
            fluidalgo_->FluidField()->Discretization();
        fluid_dis->FillComplete();

        fluidalgo_->FluidField()->PrepareTimeStep();
        fluidalgo_->FluidField()->Solve();

        Post_Process_Fluid();

        UTILS::executeSerial(comm_,
            [&]()
            {
              for (int id = 1; id < coupling_id_max; ++id)
              {
                cvOneDSynchronizer_->Set_3d_q_at_t(t_next, q_3d[id], id);
                cvOneDSynchronizer_->Set_3d_p_at_t(t_next, p_3d[id], id);
              }
            });
        Synch_Step(time_step);
        for (int i = 1; i < coupling_id_max; ++i)
        {
          p_norm[i] = (p_3d[i] - p_1d[i]) * (p_3d[i] - p_1d[i]);
          q_norm[i] = (abs(q_3d[i]) - abs(q_1d[i])) * (abs(q_3d[i]) - abs(q_1d[i]));
        }


        UTILS::executeSerial(comm_, [&]() { cvOneDSynchronizer_->Print(); });

        // create norms for every condition
        for (int i = 1; i < coupling_id_max; ++i)
        {
          if (p_norm[i] < tol_p and q_norm[i] < tol_q)
          {
            converged_condition[i] = true;
          }
          if (sqrt(p_norm[i]) / (std::pow(std::max(p_3d[i], p_3d[i]), 1)) < tol_p_rel and
              sqrt(q_norm[i]) / (std::pow(std::max(q_3d[i], q_3d[i]), 1)) < tol_q_rel)
          {
            converged_condition[i] = true;
          }

          if ((abs(p_norm[i] - p_norm_prev[i])) < tol_p_c and
              (abs(q_norm[i] - q_norm_prev[i])) < tol_q_c)
          {
            norms_stayed_same[i] = true;
          }
        }

        if ((long unsigned int)std::accumulate(converged_condition.begin() + 1,
                converged_condition.end(), 0) == converged_condition.size() - 1)
        {
          if (comm_.MyPID() == 0)
          {
            std::cout << "p_norm and q_norm converged" << std::endl;
          }
          break;
        }
        // TODO this makes no sense, always breaks at iter=5!
        if (((long unsigned int)std::accumulate(norms_stayed_same.begin() + 1,
                 norms_stayed_same.end(), 0) == norms_stayed_same.size() - 1) and
            iter > 4)
        {
          if (comm_.MyPID() == 0)
          {
            std::cout << "p_rel_prev stayed same... breaking" << std::endl;
          }
          break;
        }
        p_norm_prev = p_norm;
        q_norm_prev = q_norm;

        iter++;
      }

      partition_iteration_count.push_back(iter);

      // TODO save this for per condition?
      vec_q_norm.push_back(q_norm[id_list[0]]);
      vec_p_norm.push_back(p_norm[id_list[0]]);
      vec_p_norm_rel.push_back(p_norm[id_list[0]] / p_3d[id_list[0]]);

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