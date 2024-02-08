/*----------------------------------------------------------------------*/
/*! \file

\brief Unit tests for integration of cvOneDBFSolver with baci

\level 1
*/
// End doxygen header.

#include <gtest/gtest.h>

#include "baci_so3_hex8.H"

#include <cvOneDSynchronizer.h>
#include <OneDSolverInterface.h>

#include <iostream>

namespace
{
  using namespace BACI;


  class OneDSolverTest : public ::testing::Test
  {
   protected:
    void SetUp() override {}

    //! testing parameters
    static constexpr double TOL = 1.0e-4;

    //! read input file
    // TODO fix absolute file path
    string inputFile =
        "/home/a11bmafr/software/baci/baci/unittests/cvOneDBFSolver/cv_velocity_Coupling.in";
  };

  TEST_F(OneDSolverTest, Setup)
  {
    // create model manager
    OneDSolverInterface myOneDSolver = OneDSolverInterface();

    // Create Solver Options of cvOneD
    cvOneDOptions* opts = new cvOneDOptions();

    // Read Model From File
    myOneDSolver.readModel(inputFile, opts);

    // perform model check
    opts->check();

    // Create Model and Run Simulation
    // myOneDSolver.setupModelManager(opts);

    // Delete Options
    delete opts;

    // currently this test does not anything to test?
    // If the setup of the Model Manger fails, an error is generated and detected.
  }

  TEST_F(OneDSolverTest, SetupuntilNewton)
  {
    // create model manager
    OneDSolverInterface myOneDSolver = OneDSolverInterface();

    // Create Solver Options of cvOneD
    cvOneDOptions* opts = new cvOneDOptions();

    // Read Model From File
    myOneDSolver.readModel(inputFile, opts);

    // perform model check
    opts->check();

    // set up dummy synch
    cvOneDSynchronizer* pSynchronizer = new cvOneDSynchronizer();
    // Create Model and Run Simulation
    myOneDSolver.setupModeluntilNewton(opts, pSynchronizer);

    myOneDSolver.UpdateTimeStep();
    int iter = 0;
    int step = 1;
    myOneDSolver.Do_Newton_Step(&iter);
    myOneDSolver.UpdateSolution(iter, step);
  }


  TEST_F(OneDSolverTest, cvOneDSynchronizer)
  {
    // This is a small integration test to check if the velocity is right set and evaluated in the
    // one d artery library

    // create model manager
    OneDSolverInterface myOneDSolver = OneDSolverInterface();

    // Create Solver Options of cvOneD
    cvOneDOptions* opts = new cvOneDOptions();

    // Read Model From File
    myOneDSolver.readModel(inputFile, opts);

    // perform model check
    opts->check();

    const double dt = 0.01;
    const int steps = 10;
    // create the Synchronizer for data coupling
    cvOneDSynchronizer* pSynchronizer = new cvOneDSynchronizer(steps + 1, dt);

    // Set the 3d coupling values at the time steps
    for (int i = 1; i < steps; ++i)
    {
      pSynchronizer->Set_3d_q_at_t(dt * i, 0.002 * i);
    }

    // Create Model and Run Simulation
    myOneDSolver.createAndRunModel(opts, pSynchronizer);

    // Delete Options
    delete opts;

    // test if we have written back the solution
    // values can be found in velocity_CouplingARTERY_pressure.dat
    // they should correspond to the first row exept the first entry since the initial state is not
    // extracted?
    std::array<double, steps> solution{{0.0, 112773.79784974219, 109676.68791493628,
        100722.24699280947, 83542.696068952733, 57956.489380590065, 26243.264610449332,
        -7495.9862413965748, -38436.425457495294, -77023.547305151558}};

    // compare the results
    for (int i = 0; i < steps; ++i)
    {
      EXPECT_NEAR(solution[i], pSynchronizer->Get_1d_p_at_t(i * dt), OneDSolverTest::TOL);
    }
  }


  TEST_F(OneDSolverTest, partitioned)
  {
    /*double q_1d = 0;
    double q_3d = 0;
    double p_3d = 0;
    double p_1d = 0;
    double t = 0;
    double t_end = 1;
    double tol = 0.03;*/

    // Initialize the OneDSolver
    // create model manager
    OneDSolverInterface myOneDSolver = OneDSolverInterface();

    // Create Solver Options of cvOneD
    cvOneDOptions* opts = new cvOneDOptions();

    // Read Model From File
    myOneDSolver.readModel(inputFile, opts);

    // perform model check
    opts->check();

    const double dt = 0.01;
    const int steps = 10;
    // create the Synchronizer for data coupling
    cvOneDSynchronizer* pSynchronizer = new cvOneDSynchronizer(steps + 1, dt);

    // Create Model and Run Simulation
    myOneDSolver.setupModeluntilNewton(opts, pSynchronizer);

    myOneDSolver.UpdateTimeStep();
    // for every time step
    //  while(t < t_end) {

    //
    //      while (abs(q_3d - q_1d) > tol and abs(p_3d - p_1d) > tol) {



    // set pressure in simulation

    // perform baci solve

    // get outflow q3d and p3d


    // set q_1d simulation from q_3d

    /*int iter = 0;
    int time_step = 0;
    while (true)
    {
      // Newton-Raphson Iterations...
    /  if (!myOneDSolver.Do_Newton_Step(&iter))
      {
        break;
      }


    }  // End while
    myOneDSolver.SynchronizeDataofStep(time_step);

    // get p_1d



    //     itermaxx++;
  }*/

    //}
  }


}  // namespace
