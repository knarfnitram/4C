/*----------------------------------------------------------------------*/
/*! \file

 \brief  monolithic fluid split poroelasticity algorithms

\level 2

 *------------------------------------------------------------------------------------------------*/

#include "baci_poroelast_monolithicfluidsplit.H"

#include "baci_adapter_fld_poro.H"
#include "baci_adapter_str_fpsiwrapper.H"
#include "baci_coupling_adapter_converter.H"
#include "baci_fluid_utils_mapextractor.H"
#include "baci_fsi_overlapprec_fsiamg.H"
#include "baci_fsi_statustest.H"
#include "baci_linalg_matrixtransform.H"
#include "baci_linalg_utils_sparse_algebra_assemble.H"
#include "baci_linalg_utils_sparse_algebra_create.H"
#include "baci_poroelast_nox_group.H"
#include "baci_structure_aux.H"

#include <Teuchos_TimeMonitor.hpp>

#define FLUIDSPLITAMG


POROELAST::MonolithicFluidSplit::MonolithicFluidSplit(const Epetra_Comm& comm,
    const Teuchos::ParameterList& timeparams,
    Teuchos::RCP<CORE::LINALG::MapExtractor> porosity_splitter)
    : MonolithicSplit(comm, timeparams, porosity_splitter)
{
  fggtransform_ = Teuchos::rcp(new CORE::LINALG::MatrixRowColTransform);
  fgitransform_ = Teuchos::rcp(new CORE::LINALG::MatrixRowTransform);
  figtransform_ = Teuchos::rcp(new CORE::LINALG::MatrixColTransform);
  cfggtransform_ = Teuchos::rcp(new CORE::LINALG::MatrixRowTransform);
  csggtransform_ = Teuchos::rcp(new CORE::LINALG::MatrixColTransform);
  cfgitransform_ = Teuchos::rcp(new CORE::LINALG::MatrixRowTransform);
  csigtransform_ = Teuchos::rcp(new CORE::LINALG::MatrixColTransform);

  // Recovering of Lagrange multiplier happens on structure field
  lambda_ = Teuchos::rcp(new Epetra_Vector(*FluidField()->Interface()->FSICondMap()));
}

void POROELAST::MonolithicFluidSplit::SetupSystem()
{
  {
    // create combined map
    std::vector<Teuchos::RCP<const Epetra_Map>> vecSpaces;

    vecSpaces.push_back(StructureField()->DofRowMap());
#ifdef FLUIDSPLITAMG
    vecSpaces.push_back(FluidField()->DofRowMap());
#else
    vecSpaces.push_back(FluidField()->Interface()->OtherMap());
#endif

    if (vecSpaces[0]->NumGlobalElements() == 0) dserror("No structure equation. Panic.");
    if (vecSpaces[1]->NumGlobalElements() == 0) dserror("No fluid equation. Panic.");

    // full Poroelasticity-map
    fullmap_ = CORE::LINALG::MultiMapExtractor::MergeMaps(vecSpaces);
    // full Poroelasticity-blockmap
    blockrowdofmap_->Setup(*fullmap_, vecSpaces);
  }

  // Switch fluid to interface split block matrix
  FluidField()->UseBlockMatrix(true);

  SetupCouplingAndMatrices();

  BuildCombinedDBCMap();

  SetupEquilibration();
}
/*
void POROELAST::MonolithicFluidSplit::SetupRHS(Epetra_Vector& f, bool firstcall)
{
  TEUCHOS_FUNC_TIME_MONITOR("POROELAST::MonolithicFluidSplit::SetupRHS");

  // create full monolithic rhs vector
  rhs_ = Teuchos::rcp(new Epetra_Vector(*DofRowMap(), true));

  SetupVector(*rhs_, StructureField()->RHS(), FluidField()->RHS(), FluidField()->ResidualScaling());
  // Extractor()->AddVector(*StructureField()->RHS(), 0, *rhs_,FluidField()->ResidualScaling());
  // Extractor()->AddVector(*FluidField()->RHS(), 1, *rhs_,FluidField()->ResidualScaling());
  if (firstcall and evaluateinterface_)
  {
    // add additional rhs-terms depending on the solution on the interface
    // of the previous time step

    // get time integration parameters of structure an fluid time integrators
    // to enable consistent time integration among the fields
    double stiparam = StructureField()->TimIntParam();
    double ftiparam = FluidField()->TimIntParam();

    Teuchos::RCP<CORE::LINALG::BlockSparseMatrixBase> blockf = FluidField()->BlockSystemMatrix();
    Teuchos::RCP<CORE::LINALG::BlockSparseMatrixBase> k_sf = StructFluidCouplingBlockMatrix();
    if (k_sf == Teuchos::null) dserror("expect coupling block matrix");

    CORE::LINALG::SparseMatrix& fig = blockf->Matrix(0, 1);
    CORE::LINALG::SparseMatrix& fgg = blockf->Matrix(1, 1);
    CORE::LINALG::SparseMatrix& kig = k_sf->Matrix(0, 1);
    CORE::LINALG::SparseMatrix& kgg = k_sf->Matrix(1, 1);

    Teuchos::RCP<Epetra_Vector> fveln = FluidField()->ExtractInterfaceVeln();
    double timescale = FluidField()->TimeScaling();
    double scale = FluidField()->ResidualScaling();

    Teuchos::RCP<Epetra_Vector> rhs = Teuchos::rcp(new Epetra_Vector(fig.RowMap()));

    fig.Apply(*fveln, *rhs);
    rhs->Scale(timescale * Dt());

#ifdef FLUIDSPLITAMG
    rhs = FluidField()->Interface()->InsertOtherVector(rhs);
#endif

    Extractor()->AddVector(*rhs, 1, *rhs_);  // add fluid contributions to 'f'

    rhs = Teuchos::rcp(new Epetra_Vector(fgg.RowMap()));

    fgg.Apply(*fveln, *rhs);
    rhs->Scale(scale * timescale * Dt());
    rhs->Scale(
        (1.0 - stiparam) / (1.0 - ftiparam));  // scale 'rhs' due to consistent time integration

    rhs = FluidToStructureAtInterface(rhs);
    rhs = StructureField()->Interface()->InsertFSICondVector(rhs);

    Extractor()->AddVector(*rhs, 0, *rhs_);  // add structure contributions to 'f'

    rhs = Teuchos::rcp(new Epetra_Vector(kig.RowMap()));

    kig.Apply(*fveln, *rhs);
    rhs->Scale(timescale * Dt());

    rhs = StructureField()->Interface()->InsertOtherVector(rhs);

    Extractor()->AddVector(*rhs, 0, *rhs_);  // add structure contributions to 'f'

    rhs = Teuchos::rcp(new Epetra_Vector(kgg.RowMap()));

    kgg.Apply(*fveln, *rhs);
    rhs->Scale(timescale * Dt());

    rhs = StructureField()->Interface()->InsertFSICondVector(rhs);

    Extractor()->AddVector(*rhs, 0, *rhs_);  // add structure contributions to 'f'
  }

  // store interface force onto the structure to know it in the next time step as previous force
  // in order to recover the Lagrange multiplier
  fgcur_ = FluidField()->Interface()->ExtractFSICondVector(FluidField()->RHS());
  f.Update(-1.0, *rhs_, 0.0);
}
*/

void POROELAST::MonolithicFluidSplit::SetupRHS(bool firstcall)
{
  TEUCHOS_FUNC_TIME_MONITOR("POROELAST::MonolithicFluidSplit::SetupRHS");

  // create full monolithic rhs vector
  rhs_ = Teuchos::rcp(new Epetra_Vector(*DofRowMap(), true));

  SetupVector(*rhs_, StructureField()->RHS(), FluidField()->RHS(), FluidField()->ResidualScaling());

  if (firstcall and evaluateinterface_)
  {
    // add additional rhs-terms depending on the solution on the interface
    // of the previous time step

    // get time integration parameters of structure an fluid time integrators
    // to enable consistent time integration among the fields
    double stiparam = StructureField()->TimIntParam();
    double ftiparam = FluidField()->TimIntParam();

    Teuchos::RCP<CORE::LINALG::BlockSparseMatrixBase> blockf = FluidField()->BlockSystemMatrix();
    Teuchos::RCP<CORE::LINALG::BlockSparseMatrixBase> k_sf = StructFluidCouplingBlockMatrix();
    if (k_sf == Teuchos::null) dserror("expect coupling block matrix");

    CORE::LINALG::SparseMatrix& fig = blockf->Matrix(0, 1);
    CORE::LINALG::SparseMatrix& fgg = blockf->Matrix(1, 1);
    CORE::LINALG::SparseMatrix& kig = k_sf->Matrix(0, 1);
    CORE::LINALG::SparseMatrix& kgg = k_sf->Matrix(1, 1);

    Teuchos::RCP<Epetra_Vector> fveln = FluidField()->ExtractInterfaceVeln();
    double timescale = FluidField()->TimeScaling();
    double scale = FluidField()->ResidualScaling();

    Teuchos::RCP<Epetra_Vector> rhs = Teuchos::rcp(new Epetra_Vector(fig.RowMap()));

    fig.Apply(*fveln, *rhs);
    rhs->Scale(timescale * Dt());

#ifdef FLUIDSPLITAMG
    rhs = FluidField()->Interface()->InsertOtherVector(rhs);
#endif

    Extractor()->AddVector(*rhs, 1, *rhs_);  // add fluid contributions to 'f'

    rhs = Teuchos::rcp(new Epetra_Vector(fgg.RowMap()));

    fgg.Apply(*fveln, *rhs);
    rhs->Scale(scale * timescale * Dt());
    rhs->Scale(
        (1.0 - stiparam) / (1.0 - ftiparam));  // scale 'rhs' due to consistent time integration

    rhs = FluidToStructureAtInterface(rhs);
    rhs = StructureField()->Interface()->InsertFSICondVector(rhs);

    Extractor()->AddVector(*rhs, 0, *rhs_);  // add structure contributions to 'f'

    rhs = Teuchos::rcp(new Epetra_Vector(kig.RowMap()));

    kig.Apply(*fveln, *rhs);
    rhs->Scale(timescale * Dt());

    rhs = StructureField()->Interface()->InsertOtherVector(rhs);

    Extractor()->AddVector(*rhs, 0, *rhs_);  // add structure contributions to 'f'

    rhs = Teuchos::rcp(new Epetra_Vector(kgg.RowMap()));

    kgg.Apply(*fveln, *rhs);
    rhs->Scale(timescale * Dt());

    rhs = StructureField()->Interface()->InsertFSICondVector(rhs);

    Extractor()->AddVector(*rhs, 0, *rhs_);  // add structure contributions to 'f'
  }

  // store interface force onto the structure to know it in the next time step as previous force
  // in order to recover the Lagrange multiplier
  fgcur_ = FluidField()->Interface()->ExtractFSICondVector(FluidField()->RHS());
}
void POROELAST::MonolithicFluidSplit::DoTimeStepNew(
    const Teuchos::RCP<NOX::Epetra::Interface::Required>& interface)
{
  // counter and print header
  // predict solution of both field (call the adapter)
  MonolithicSplit::PrepareTimeStep();

  TEUCHOS_FUNC_TIME_MONITOR("POROELAST::MonolithicFluidSplit::TimeStep");

  // Get the top level parameter list
  Teuchos::ParameterList& nlParams = noxparameterlist_;

  // sublistsc

  Teuchos::ParameterList& dirParams = nlParams.sublist("Direction");
  // Teuchos::ParameterList& solverOptions = nlParams.sublist("Solver Options");
  Teuchos::ParameterList& newtonParams = dirParams.sublist("Newton");
  // newtonParams.set("Forcing Term Method", "Constant");
  Teuchos::ParameterList& lsParams = newtonParams.sublist("Linear Solver");
  // Teuchos::ParameterList& newtonParams= dirParams.sublist("Forcing Term Initial Tolerance");

  Teuchos::ParameterList& searchParams = nlParams.sublist("Line Search");
  // searchParams.set("Method", "Backtrack");
  searchParams.set<int>("Forcing Term Initial Tolerance", 0.01);
  Teuchos::ParameterList& printParams = nlParams.sublist("Printing");
  printParams.set("MyPID", Comm().MyPID());
  zeros_ = CORE::LINALG::CreateVector(*DofRowMap(), true);
  // Evaluate(zeros_, false);
  Evaluate(Teuchos::null, true);
  if (iterinc_.is_null()) iterinc_ = Teuchos::rcp(new Epetra_Vector(*DofRowMap(), true));
  // EvaluateNOX(iterinc_);
  // SetupRHS(true);
  // std::cout<<*iterinc_<< std::endl;
  // SetupRHS(true);
  std::cout << *systemmatrix_ << std::endl;
  std::cout << *rhs_ << std::endl;

  // SetupSystemMatrix(*systemmatrix_);
  std::cout << "After setupsystemmatrix" << std::endl;
  std::cout << *systemmatrix_ << std::endl;
  std::cout << *rhs_ << std::endl;

  // FluidField()->Evaluate(Teuchos::null);
  // StructureField()->Evaluate(Teuchos::null);
  // EvaluateNOX(Teuchos::null);
  //  Get initial guess.
  //  The initial system is there, so we can happily extract the
  //  initial guess. (The Dirichlet conditions are already build in!)
  Teuchos::RCP<Epetra_Vector> initial_guess = Teuchos::rcp(new Epetra_Vector(*DofRowMap(), true));
  // Evaluate(initial_guess,true);

  InitialGuess(initial_guess);

  NOX::Epetra::Vector noxSoln(initial_guess, NOX::Epetra::Vector::CreateView);

  // Create the linear system
  Teuchos::RCP<NOX::Epetra::LinearSystem> linSys = CreateLinearSystem(nlParams, noxSoln);

  // const Teuchos::RCP<NOX::Epetra::Interface::Required>& interface =
  // static_cast<Teuchos::RCP<NOX::Epetra::Interface::Required>>(this);
  //  Create the Group
  Teuchos::RCP<NOX::POROELAST::Group> grp =
      Teuchos::rcp(new NOX::POROELAST::Group(*this, printParams, interface, noxSoln, linSys));

  // Convergence Tests
  Teuchos::RCP<NOX::StatusTest::Combo> combo = this->CreateStatusTest(nlParams, grp);

  // Create the solver
  Teuchos::RCP<NOX::Solver::Generic> solver =
      NOX::Solver::buildSolver(grp, combo, Teuchos::RCP<Teuchos::ParameterList>(&nlParams, false));

  nox_prev_ = Teuchos::rcp(new Epetra_Vector(*DofRowMap(), true));

  if (zeros_.is_null()) zeros_ = Teuchos::rcp(new Epetra_Vector(*DofRowMap(), true));
  auto zeros_2 = Teuchos::rcp(new Epetra_Vector(*DofRowMap(), true));

  // we know we already have the first linear system calculated
  // SetupSolver();

  std::cout << "systemmatrix_ in time step after dirichlet" << std::endl;
  std::cout << *systemmatrix_ << std::endl;

  CORE::LINALG::ApplyDirichletToSystem(*systemmatrix_, *zeros_2, *rhs_, *zeros_, *CombinedDBCMap());
  grp->CaptureSystemState();

  std::ostringstream oss;
  if (grp->getF().length() == 0) dserror("well thats bad");
  // solve the whole thing


  std::cout << "systemmatrix_ in time step (after diricht) " << std::endl;
  std::cout << *systemmatrix_ << std::endl;
  noxstatus_ = solver->solve();
  noxiter_ = solver->getNumIterations();
  grp->getX().print(oss);
  // iterinc_->Update(1.0,grp->getX().createMultiVector(2,NOX::DeepCopy),0);
  std::cout << oss.str() << std::endl;
  // BuildConvergenceNorms();
  // PrintNewtonIter();

  // Newton-Raphson iteration
  // Solve();

  // calculate stresses, strains, energies
  constexpr bool force_prepare = false;
  PrepareOutput(force_prepare);

  // update all single field solvers
  Update();

  // write output to screen and files
  Output();
}


void POROELAST::MonolithicFluidSplit::TimeLoopNew(
    const Teuchos::RCP<NOX::Epetra::Interface::Required>& interface)
{
  while (NotFinished())
  {
    // solve one time step
    DoTimeStepNew(interface);
  }
}

void POROELAST::MonolithicFluidSplit::SetupSystemMatrix(CORE::LINALG::BlockSparseMatrixBase& mat)
{
  TEUCHOS_FUNC_TIME_MONITOR("POROELAST::MonolithicFluidSplit::SetupSystemMatrix");

  Teuchos::RCP<CORE::LINALG::SparseMatrix> s = StructureField()->SystemMatrix();
  if (s == Teuchos::null) dserror("expect structure matrix");
  Teuchos::RCP<CORE::LINALG::BlockSparseMatrixBase> f = FluidField()->BlockSystemMatrix();
  if (f == Teuchos::null) dserror("expect fluid block matrix");

  mat.Matrix(0, 1).Zero();
  mat.Matrix(1, 0).Zero();
#ifdef FLUIDSPLITAMG
  mat.Matrix(1, 1).Zero();
#endif

  /*----------------------------------------------------------------------*/

  // build block matrix
  // The maps of the block matrix have to match the maps of the blocks we
  // insert here.

  /*----------------------------------------------------------------------*/
  // structural part k_sf (3nxn)
  // build mechanical-fluid block

  // create empty matrix
  Teuchos::RCP<CORE::LINALG::BlockSparseMatrixBase> k_sf = StructFluidCouplingBlockMatrix();
  if (k_sf == Teuchos::null) dserror("expect coupling block matrix");

  // call the element and calculate the matrix block
  ApplyStrCouplMatrix(k_sf);

  /*----------------------------------------------------------------------*/
  // fluid part k_fs ( (3n+1)x3n )
  // build fluid-mechanical block

  // create empty matrix
  Teuchos::RCP<CORE::LINALG::BlockSparseMatrixBase> k_fs = FluidStructCouplingBlockMatrix();
  if (k_fs == Teuchos::null) dserror("expect coupling block matrix");

  // call the element and calculate the matrix block
  ApplyFluidCouplMatrix(k_fs);

  /*----------------------------------------------------------------------*/

  k_fs->Complete();
  k_sf->Complete();

  s->UnComplete();

  /*----------------------------------------------------------------------*/

  if (evaluateinterface_)
  {
    double timescale = FluidField()->TimeScaling();

    (*figtransform_)(f->FullRowMap(), f->FullColMap(), f->Matrix(0, 1), timescale,
        CORE::ADAPTER::CouplingSlaveConverter(*icoupfs_), k_fs->Matrix(0, 1), true, true);

    (*csggtransform_)(f->FullRowMap(), f->FullColMap(), k_sf->Matrix(1, 1), timescale,
        CORE::ADAPTER::CouplingSlaveConverter(*icoupfs_), *s, true, true);

    (*csigtransform_)(f->FullRowMap(), f->FullColMap(), k_sf->Matrix(0, 1), timescale,
        CORE::ADAPTER::CouplingSlaveConverter(*icoupfs_), *s, true, true);
  }

  /*----------------------------------------------------------------------*/
  // pure fluid part
  // uncomplete because the fluid interface can have more connections than the
  // structural one. (Tet elements in fluid can cause this.) We should do
  // this just once...
#ifdef FLUIDSPLITAMG
  mat.Matrix(1, 1).Add(f->Matrix(0, 0), false, 1., 0.0);
  Teuchos::RCP<CORE::LINALG::SparseMatrix> eye =
      CORE::LINALG::Eye(*FluidField()->Interface()->FSICondMap());
  mat.Matrix(1, 1).Add(*eye, false, 1., 1.0);
#else
  f->Matrix(0, 0).UnComplete();
  mat.Assign(1, 1, View, f->Matrix(0, 0));
#endif

  // fluid coupling part
  mat.Matrix(1, 0).Add(k_fs->Matrix(0, 0), false, 1.0, 0.0);
  mat.Matrix(1, 0).Add(k_fs->Matrix(0, 1), false, 1.0, 1.0);

  // pure structure part
  mat.Assign(0, 0, CORE::LINALG::View, *s);

  // structure coupling part
  mat.Matrix(0, 1).Add(k_sf->Matrix(0, 0), false, 1.0, 0.0);
  mat.Matrix(0, 1).Add(k_sf->Matrix(1, 0), false, 1.0, 1.0);
  /*----------------------------------------------------------------------*/
  // done. make sure all blocks are filled.
  mat.Complete();


  fgicur_ = Teuchos::rcp(new CORE::LINALG::SparseMatrix(f->Matrix(1, 0)));
  fggcur_ = Teuchos::rcp(new CORE::LINALG::SparseMatrix(f->Matrix(1, 1)));
  cgicur_ = Teuchos::rcp(new CORE::LINALG::SparseMatrix(k_fs->Matrix(1, 0)));
  cggcur_ = Teuchos::rcp(new CORE::LINALG::SparseMatrix(k_fs->Matrix(1, 1)));
  mat.Merge(true);
  mat.ApplyDirichlet(*(combinedDBCMap_), true);
}

void POROELAST::MonolithicFluidSplit::SetupVector(Epetra_Vector& f,
    Teuchos::RCP<const Epetra_Vector> sv, Teuchos::RCP<const Epetra_Vector> fv, double fluidscale)
{
  // extract the inner and boundary dofs of all fields

  Teuchos::RCP<Epetra_Vector> fov = FluidField()->Interface()->ExtractOtherVector(fv);
#ifdef FLUIDSPLITAMG
  fov = FluidField()->Interface()->InsertOtherVector(fov);
#endif

  Extractor()->InsertVector(*sv, 0, f);

  Extractor()->InsertVector(*fov, 1, f);  // add fluid contributions to 'f'
}

void POROELAST::MonolithicFluidSplit::ExtractFieldVectors(Teuchos::RCP<const Epetra_Vector> x,
    Teuchos::RCP<const Epetra_Vector>& sx, Teuchos::RCP<const Epetra_Vector>& fx, bool firstcall)
{
  TEUCHOS_FUNC_TIME_MONITOR("POROELAST::MonolithicFluidSplit::ExtractFieldVectors");

  // process structure unknowns
  sx = Extractor()->ExtractVector(x, 0);

  // process fluid unknowns
  if (evaluateinterface_)
  {
    Teuchos::RCP<const Epetra_Vector> scx = StructureField()->Interface()->ExtractFSICondVector(sx);

    Teuchos::RCP<Epetra_Vector> fcx = StructureToFluidAtInterface(scx);
    Teuchos::RCP<const Epetra_Vector> fox = Extractor()->ExtractVector(x, 1);
#ifdef FLUIDSPLITAMG
    fox = FluidField()->Interface()->ExtractOtherVector(fox);
#endif

    {
      double timescale = FluidField()->TimeScaling();
      fcx->Scale(timescale);
    }

    Teuchos::RCP<Epetra_Vector> f = FluidField()->Interface()->InsertOtherVector(fox);
    FluidField()->Interface()->InsertFSICondVector(fcx, f);

    auto zeros = Teuchos::rcp(new const Epetra_Vector(f->Map(), true));
    CORE::LINALG::ApplyDirichletToSystem(
        *f, *zeros, *(FluidField()->GetDBCMapExtractor()->CondMap()));

    fx = f;

    // Store field vectors to know them later on as previous quantities
    Teuchos::RCP<Epetra_Vector> sox = StructureField()->Interface()->ExtractOtherVector(sx);
    if (solipre_ != Teuchos::null)
      ddiinc_->Update(1.0, *sox, -1.0, *solipre_, 0.0);  // compute current iteration increment
    else
      ddiinc_ = Teuchos::rcp(new Epetra_Vector(*sox));  // first iteration increment

    solipre_ = sox;  // store current step increment

    if (solgvelpre_ != Teuchos::null)
      duginc_->Update(1.0, *fcx, -1.0, *solgvelpre_, 0.0);  // compute current iteration increment
    else
      duginc_ = Teuchos::rcp(new Epetra_Vector(*fcx));  // first iteration increment

    solgvelpre_ = fcx;  // store current step increment

    if (solivelpre_ != Teuchos::null)
      duiinc_->Update(1.0, *fox, -1.0, *solivelpre_, 0.0);  // compute current iteration increment
    else
      duiinc_ = Teuchos::rcp(new Epetra_Vector(*fox));  // first iteration increment

    solivelpre_ = fox;  // store current step increment
  }
  else
    fx = Extractor()->ExtractVector(x, 1);
}

void POROELAST::MonolithicFluidSplit::RecoverLagrangeMultiplierAfterTimeStep()
{
  TEUCHOS_FUNC_TIME_MONITOR("POROELAST::MonolithicFluidSplit::RecoverLagrangeMultiplier");

  if (evaluateinterface_)
  {
    // get time integration parameter of structural time integrator
    // to enable consistent time integration among the fields
    double ftiparam = FluidField()->TimIntParam();
    double timescale = FluidField()->TimeScaling();

    // store the product F_{\GammaI} \Delta u_I^{n+1} in here
    Teuchos::RCP<Epetra_Vector> fgiddi =
        CORE::LINALG::CreateVector(*FluidField()->Interface()->FSICondMap(), true);
    // compute the above mentioned product
    fgicur_->Multiply(false, *duiinc_, *fgiddi);

    // store the product C_{\GammaI} \Delta d_I^{n+1} in here
    Teuchos::RCP<Epetra_Vector> sgiddi =
        CORE::LINALG::CreateVector(*FluidField()->Interface()->FSICondMap(), true);
    // compute the above mentioned product
    cgicur_->Multiply(false, *ddiinc_, *sgiddi);

    // store the product F_{\Gamma\Gamma} \Delta u_\Gamma^{n+1} in here
    Teuchos::RCP<Epetra_Vector> sggddg =
        CORE::LINALG::CreateVector(*FluidField()->Interface()->FSICondMap(), true);
    // compute the above mentioned product
    fggcur_->Multiply(false, *duginc_, *sggddg);

    // store the prodcut C_{\Gamma\Gamma} \Delta d_\Gamma^{n+1} in here
    Teuchos::RCP<Epetra_Vector> cggddg =
        CORE::LINALG::CreateVector(*FluidField()->Interface()->FSICondMap(), true);
    // compute the above mentioned product
    cggcur_->Multiply(false, *duginc_, *cggddg);
    cggddg->Scale(1.0 / timescale);

    // Update the Lagrange multiplier:
    /* \lambda^{n+1} =  1/b * [ - a*\lambda^n - f_\Gamma^S
     *                          - F_{\Gamma I} \Delta u_I
     *                          - C_{\Gamma I} \Delta d_I
     *                          - F_{\Gamma\Gamma} \Delta u_\Gamma]
     *                          - C_{\Gamma\Gamma} * Delta t / 2 * \Delta u_\Gamma]
     */
    lambda_->Update(1.0, *fgcur_, -ftiparam);
    lambda_->Update(-1.0, *fgiddi, -1.0, *sggddg, 1.0);
    lambda_->Update(-1.0, *fgiddi, -1.0, *cggddg, 1.0);
    lambda_->Scale(1 / (1.0 - ftiparam));  // entire Lagrange multiplier is divided by (1.-ftiparam)
  }
}

Teuchos::RCP<NOX::StatusTest::Combo> POROELAST::MonolithicFluidSplit::CreateStatusTest(
    Teuchos::ParameterList& nlParams, Teuchos::RCP<NOX::Epetra::Group> grp)
{
  // --------------------------------------------------------------------
  // Setup the test framework
  // --------------------------------------------------------------------
  // Create the top-level test combo
  Teuchos::RCP<NOX::StatusTest::Combo> combo =
      Teuchos::rcp(new NOX::StatusTest::Combo(NOX::StatusTest::Combo::OR));

  // Create test combo for convergence of residuals and iterative increments
  Teuchos::RCP<NOX::StatusTest::Combo> converged =
      Teuchos::rcp(new NOX::StatusTest::Combo(NOX::StatusTest::Combo::AND));

  // Create some other plausibility tests
  Teuchos::RCP<NOX::StatusTest::MaxIters> maxiters =
      Teuchos::rcp(new NOX::StatusTest::MaxIters(nlParams.get<int>("Max Iterations")));
  Teuchos::RCP<NOX::StatusTest::FiniteValue> fv = Teuchos::rcp(new NOX::StatusTest::FiniteValue);

  // Add single tests to the top-level test combo
  combo->addStatusTest(fv);
  combo->addStatusTest(converged);
  combo->addStatusTest(maxiters);

  // Start filling the 'converged' combo here
  // require one solve
  converged->addStatusTest(
      Teuchos::rcp(new NOX::FSI::MinIters(nlParams.get<int>("Min Iterations"))));

  // TODO implement logic of the norms here

  // --------------------------------------------------------------------
  // setup tests for structural displacement field
  // --------------------------------------------------------------------
  // create NOX::StatusTest::Combo for structural displacement field
  Teuchos::RCP<NOX::StatusTest::Combo> structcombo =
      Teuchos::rcp(new NOX::StatusTest::Combo(NOX::StatusTest::Combo::AND));

  // create Norm-objects for each norm that has to be tested
  Teuchos::RCP<NOX::FSI::PartialNormF> structureDisp_L2 = Teuchos::rcp(new NOX::FSI::PartialNormF(
      "DISPL residual", *Extractor(), 0, nlParams.get<double>("Tol res dis"),
      NOX::Abstract::Vector::TwoNorm, NOX::FSI::PartialNormF::Scaled));
  Teuchos::RCP<NOX::FSI::PartialNormF> structureDisp_inf = Teuchos::rcp(new NOX::FSI::PartialNormF(
      "DISPL residual", *Extractor(), 0, nlParams.get<double>("Tol res dis"),
      NOX::Abstract::Vector::MaxNorm, NOX::FSI::PartialNormF::Unscaled));
  Teuchos::RCP<NOX::FSI::PartialNormUpdate> structureDispUpdate_L2 =
      Teuchos::rcp(new NOX::FSI::PartialNormUpdate("DISPL update", *Extractor(), 0,
          nlParams.get<double>("Tol inc dis"), NOX::Abstract::Vector::TwoNorm,
          NOX::FSI::PartialNormUpdate::Scaled));
  Teuchos::RCP<NOX::FSI::PartialNormUpdate> structureDispUpdate_inf =
      Teuchos::rcp(new NOX::FSI::PartialNormUpdate("DISPL update", *Extractor(), 0,
          nlParams.get<double>("Tol inc dis"), NOX::Abstract::Vector::MaxNorm,
          NOX::FSI::PartialNormUpdate::Unscaled));

  // tests needed to adapt relative tolerance of the linear solver
  AddStatusTest(structureDisp_L2);

  // add norm-tests to structural displacement NOX::StatusTest::Combo
  structcombo->addStatusTest(structureDisp_L2);
  structcombo->addStatusTest(structureDisp_inf);
  structcombo->addStatusTest(structureDispUpdate_L2);
  structcombo->addStatusTest(structureDispUpdate_inf);

  // add structural displacement test combo to top-level test combo
  converged->addStatusTest(structcombo);

  // ---------- end of structural displacement field tests

  // --------------------------------------------------------------------
  // setup tests for interface
  // --------------------------------------------------------------------
  // build mapextractor
  // TODO find out how to calculate the norms of Jac and dx without using Do
  std::vector<Teuchos::RCP<const Epetra_Map>> interface;
  interface.push_back(DofRowMap());
  interface.push_back(Teuchos::null);
  CORE::LINALG::MultiMapExtractor interfaceextract(*DofRowMap(), interface);

  // create NOX::StatusTest::Combo for interface
  Teuchos::RCP<NOX::StatusTest::Combo> interfacecombo =
      Teuchos::rcp(new NOX::StatusTest::Combo(NOX::StatusTest::Combo::AND));

  // create Norm-objects for each norm that has to be tested
  Teuchos::RCP<NOX::FSI::PartialNormF> interfaceTest_L2 = Teuchos::rcp(new NOX::FSI::PartialNormF(
      "GAMMA residual", interfaceextract, 0, nlParams.get<double>("Tol res poro"),
      NOX::Abstract::Vector::TwoNorm, NOX::FSI::PartialNormF::Scaled));
  Teuchos::RCP<NOX::FSI::PartialNormF> interfaceTest_inf = Teuchos::rcp(new NOX::FSI::PartialNormF(
      "GAMMA residual", interfaceextract, 0, nlParams.get<double>("Tol res poro"),
      NOX::Abstract::Vector::MaxNorm, NOX::FSI::PartialNormF::Unscaled));
  Teuchos::RCP<NOX::FSI::PartialNormUpdate> interfaceTestUpdate_L2 =
      Teuchos::rcp(new NOX::FSI::PartialNormUpdate("poro update", interfaceextract, 0,
          nlParams.get<double>("Tol inc poro"), NOX::Abstract::Vector::TwoNorm,
          NOX::FSI::PartialNormUpdate::Scaled));
  Teuchos::RCP<NOX::FSI::PartialNormUpdate> interfaceTestUpdate_inf =
      Teuchos::rcp(new NOX::FSI::PartialNormUpdate("poro update", interfaceextract, 0,
          nlParams.get<double>("Tol inc poro"), NOX::Abstract::Vector::MaxNorm,
          NOX::FSI::PartialNormUpdate::Unscaled));

  // tests needed to adapt relative tolerance of the linear solver
  // AddStatusTest(interfaceTest_L2);

  // add norm-tests to interface NOX::StatusTest::Combo
  interfacecombo->addStatusTest(interfaceTest_L2);
  interfacecombo->addStatusTest(interfaceTest_inf);
  interfacecombo->addStatusTest(interfaceTestUpdate_L2);
  interfacecombo->addStatusTest(interfaceTestUpdate_inf);

  // add interface test combo to top-level test combo
  // converged->addStatusTest(interfacecombo);
  Teuchos::RCP<NOX::StatusTest::NormF> absresid = Teuchos::rcp(new NOX::StatusTest::NormF(
      1.0e-8, NOX::Abstract::Vector::TwoNorm, NOX::StatusTest::NormF::Unscaled));
  converged->addStatusTest(absresid);

  Teuchos::RCP<NOX::StatusTest::NormUpdate> update =
      Teuchos::rcp(new NOX::StatusTest::NormUpdate(1.0e-6, NOX::StatusTest::NormUpdate::Unscaled));
  converged->addStatusTest(update);
  // ---------- end of interface tests

  // --------------------------------------------------------------------
  // setup tests for fluid velocity field
  // --------------------------------------------------------------------
  // build mapextractor
  std::vector<Teuchos::RCP<const Epetra_Map>> fluid_vel;
  Teuchos::RCP<const Epetra_Vector> rhs_f = Extractor()->ExtractVector(rhs_, 1);
  auto fluid_vel_block_map = FluidField()->ExtractVelocityPart(rhs_f)->Map();
  Teuchos::RCP<Epetra_Map> vecmap = Teuchos::rcp(
      new Epetra_Map(fluid_vel_block_map.NumGlobalElements(), fluid_vel_block_map.NumMyElements(),
          fluid_vel_block_map.MyGlobalElements(), 0, fluid_vel_block_map.Comm()));
  fluid_vel.push_back(vecmap);
  fluid_vel.push_back(Teuchos::null);
  CORE::LINALG::MultiMapExtractor fluid_vel_extract(*DofRowMap(), fluid_vel);

  // create NOX::StatusTest::Combo for fluid velocity field
  Teuchos::RCP<NOX::StatusTest::Combo> fluidvelcombo =
      Teuchos::rcp(new NOX::StatusTest::Combo(NOX::StatusTest::Combo::AND));

  // create Norm-objects for each norm that has to be tested
  Teuchos::RCP<NOX::FSI::PartialNormF> innerFluidVel_L2 = Teuchos::rcp(new NOX::FSI::PartialNormF(
      "vel res l2", fluid_vel_extract, 0, nlParams.get<double>("Tol res vel"),
      NOX::Abstract::Vector::TwoNorm, NOX::FSI::PartialNormF::Scaled));
  Teuchos::RCP<NOX::FSI::PartialNormF> innerFluidVel_inf = Teuchos::rcp(new NOX::FSI::PartialNormF(
      "VELOC residual", fluid_vel_extract, 0, nlParams.get<double>("Tol res vel"),
      NOX::Abstract::Vector::MaxNorm, NOX::FSI::PartialNormF::Unscaled));
  Teuchos::RCP<NOX::FSI::PartialNormUpdate> innerFluidVelUpdate_L2 =
      Teuchos::rcp(new NOX::FSI::PartialNormUpdate("VELOC update", fluid_vel_extract, 0,
          nlParams.get<double>("Tol inc vel"), NOX::Abstract::Vector::TwoNorm,
          NOX::FSI::PartialNormUpdate::Scaled));
  Teuchos::RCP<NOX::FSI::PartialNormUpdate> innerFluidVelUpdate_inf =
      Teuchos::rcp(new NOX::FSI::PartialNormUpdate("VELOC update", fluid_vel_extract, 0,
          nlParams.get<double>("Tol inc vel"), NOX::Abstract::Vector::MaxNorm,
          NOX::FSI::PartialNormUpdate::Unscaled));

  // tests needed to adapt relative tolerance of the linear solver
  AddStatusTest(innerFluidVel_L2);

  // add norm-tests to fluid velocity NOX::StatusTest::Combo
  fluidvelcombo->addStatusTest(innerFluidVel_L2);
  fluidvelcombo->addStatusTest(innerFluidVel_inf);
  fluidvelcombo->addStatusTest(innerFluidVelUpdate_L2);
  fluidvelcombo->addStatusTest(innerFluidVelUpdate_inf);

  // add fluid velocity test combo to top-level test combo
  converged->addStatusTest(fluidvelcombo);
  // ---------- end of fluid velocity field tests

  // --------------------------------------------------------------------
  // setup tests for fluid pressure field
  // --------------------------------------------------------------------
  // build mapextractor
  std::vector<Teuchos::RCP<const Epetra_Map>> fluidpress;

  Teuchos::RCP<const Epetra_Vector> rhs_pres = Extractor()->ExtractVector(rhs_, 1);
  auto fluid_pres_block_map = FluidField()->ExtractPressurePart(rhs_pres)->Map();
  Teuchos::RCP<Epetra_Map> pres_map = Teuchos::rcp(
      new Epetra_Map(fluid_vel_block_map.NumGlobalElements(), fluid_vel_block_map.NumMyElements(),
          fluid_vel_block_map.MyGlobalElements(), 0, fluid_vel_block_map.Comm()));

  fluidpress.push_back(pres_map);
  fluidpress.push_back(Teuchos::null);
  CORE::LINALG::MultiMapExtractor fluidpressextract(*DofRowMap(), fluidpress);

  // create NOX::StatusTest::Combo for fluid pressure field
  Teuchos::RCP<NOX::StatusTest::Combo> fluidpresscombo =
      Teuchos::rcp(new NOX::StatusTest::Combo(NOX::StatusTest::Combo::AND));

  // create Norm-objects for each norm that has to be tested
  Teuchos::RCP<NOX::FSI::PartialNormF> fluidPress_L2 = Teuchos::rcp(new NOX::FSI::PartialNormF(
      "PRESS residual", fluidpressextract, 0, nlParams.get<double>("Tol res pres"),
      NOX::Abstract::Vector::TwoNorm, NOX::FSI::PartialNormF::Scaled));
  Teuchos::RCP<NOX::FSI::PartialNormF> fluidPress_inf = Teuchos::rcp(new NOX::FSI::PartialNormF(
      "PRESS residual", fluidpressextract, 0, nlParams.get<double>("Tol res pres"),
      NOX::Abstract::Vector::MaxNorm, NOX::FSI::PartialNormF::Unscaled));
  Teuchos::RCP<NOX::FSI::PartialNormUpdate> fluidPressUpdate_L2 =
      Teuchos::rcp(new NOX::FSI::PartialNormUpdate("PRESS update", fluidpressextract, 0,
          nlParams.get<double>("Tol inc pres"), NOX::Abstract::Vector::TwoNorm,
          NOX::FSI::PartialNormUpdate::Scaled));
  Teuchos::RCP<NOX::FSI::PartialNormUpdate> fluidPressUpdate_inf =
      Teuchos::rcp(new NOX::FSI::PartialNormUpdate("PRESS update", fluidpressextract, 0,
          nlParams.get<double>("Tol inc pres"), NOX::Abstract::Vector::MaxNorm,
          NOX::FSI::PartialNormUpdate::Unscaled));

  // tests needed to adapt relative tolerance of the linear solver
  AddStatusTest(fluidPress_L2);

  // add norm-tests to fluid pressure NOX::StatusTest::Combo
  fluidpresscombo->addStatusTest(fluidPress_L2);
  fluidpresscombo->addStatusTest(fluidPress_inf);
  fluidpresscombo->addStatusTest(fluidPressUpdate_L2);
  fluidpresscombo->addStatusTest(fluidPressUpdate_inf);

  // add fluid pressure test combo to top-level test combo
  converged->addStatusTest(fluidpresscombo);
  // ---------- end of fluid pressure field tests

  // Finally, return the test combo
  return combo;
}
