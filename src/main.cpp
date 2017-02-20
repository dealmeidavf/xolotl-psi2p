/**
 * Main.c, currently only able to load clusters
 */
#include <cstdlib>
#include <iostream>
#include <fstream>
#include <cassert>
#include <Reactant.h>
#include <PetscSolver.h>
#include <mpi.h>
#include <MPIUtils.h>
#include <Options.h>
#include <xolotlPerf.h>
#include <IMaterialFactory.h>
#include <TemperatureHandlerFactory.h>
#include <VizHandlerRegistryFactory.h>
#include <INetworkLoader.h>
#include <IReactionNetwork.h>
#include <SolverHandlerFactory.h>
#include <ISolverHandler.h>
#include <IReactionHandlerFactory.h>
#include <ctime>

using namespace std;
using std::shared_ptr;

namespace mainHelper
{

//*------ Global helper functions -----------------------------------------------*

//................................................................................
 void 
 printStartMessage() 
 {
  std::cout << "Starting Xolotl Plasma-Surface Interactions Simulator" << std::endl;
  // Print date and time
  std::time_t currentTime = std::time(NULL);
  std::cout << std::asctime( std::localtime( &currentTime ) ); // << std::endl;
 }

//................................................................................
 std::shared_ptr<xolotlFactory::IMaterialFactory> 
 initMaterial( Options& options ) 
 {
// Create the material factory
  auto materialOpt = options.getMaterial();
  auto dimOpt      = options.getDimensionNumber();
  auto materialFactory = 
       xolotlFactory::IMaterialFactory::createMaterialFactory( materialOpt, dimOpt );

// Initialize it with the options
  materialFactory->initializeMaterial( options );

  return materialFactory;
 }

//................................................................................
 bool 
 initTemp(Options &options) 
 {
  bool tempInitOK = xolotlFactory::initializeTempHandler(options);
  if (!tempInitOK) 
  {
   std::cerr << "Unable to initialize requested temperature. Aborting" 
             << std::endl;
   return EXIT_FAILURE;
  } 
  else
  {
   return tempInitOK;
  }
 }

//................................................................................
 bool 
 initViz( bool opts ) 
 {
  bool vizInitOK = xolotlFactory::initializeVizHandler(opts);
  if (!vizInitOK) 
  {
   std::cerr << "Unable to initialize requested visualization infrastructure. Aborting"
             << std::endl;
   return EXIT_FAILURE;
  } 
  else
  {
   return vizInitOK;
  }
 }

//................................................................................
 std::shared_ptr<xolotlSolver::PetscSolver> 
 setUpSolver( std::shared_ptr<xolotlPerf::IHandlerRegistry>    handlerRegistry,
              std::shared_ptr<xolotlFactory::IMaterialFactory> material,
              std::shared_ptr<xolotlCore::ITemperatureHandler> tempHandler,
              std::shared_ptr<xolotlCore::IReactionNetwork>    networkHandler,
              std::shared_ptr<xolotlSolver::ISolverHandler>    solvHandler,
              Options&                                         options )  
 {
 // Initialize the solver handler
  solvHandler->initializeHandlers( material, tempHandler, networkHandler, options );
 
 // Setup the solver
  auto solverInitTimer = handlerRegistry->getTimer("initSolver");
  solverInitTimer->start();
 
  std::shared_ptr<xolotlSolver::PetscSolver> solver;
  solver = std::make_shared<xolotlSolver::PetscSolver>( handlerRegistry );

  solver->setCommandLineOptions( options.getPetscArgc(), options.getPetscArgv() );

  solver->initialize( solvHandler );
 
  solverInitTimer->stop();
 
  return solver;
 }

//................................................................................
 void 
 launchPetscSolver( std::shared_ptr<xolotlSolver::PetscSolver>    solver,
                    std::shared_ptr<xolotlPerf::IHandlerRegistry> handlerRegistry )
 {
  xolotlPerf::IHardwareCounter::SpecType hwctrSpec;
  hwctrSpec.push_back(xolotlPerf::IHardwareCounter::FPOps);
  hwctrSpec.push_back(xolotlPerf::IHardwareCounter::Cycles);
  hwctrSpec.push_back(xolotlPerf::IHardwareCounter::L3CacheMisses);
 
// Launch the PetscSolver
  auto solverTimer = handlerRegistry->getTimer("solve");
  auto solverHwctr = handlerRegistry->getHardwareCounter("solve", hwctrSpec);
 
  solverTimer->start();
  solverHwctr->start();
 
// The action is here
  solver->solve();
 
  solverHwctr->stop();
  solverTimer->stop();
 }

//*------ End global helper functions -------------------------------------------*
}


//******** Main Program **********************************************************
int 
main(int argc, char** argv) 
{

// Local Declarations
 int rank;

// Check the command line arguments.
// Skip the executable name before parsing
 argc -= 1; // one for the executable name
 argv += 1; // one for the executable name

 Options opts;
 opts.readParams(argv);
 if (!opts.shouldRun()) return opts.getExitCode(); 

// Skip the name of the parameter file that was just used.
// The arguments should be empty now.
 argc -= 1;
 argv += 1;

// Extract the argument for the file name
 std::string networkFilename = opts.getNetworkFilename();
 assert(!networkFilename.empty());

//--------------------------------------------------------------------------------
// Let the games begin
 try 
 {
// Set up our performance data infrastructure.
  xolotlPerf::initialize(opts.getPerfHandlerType());

// Initialize MPI. We do this instead of leaving it to some
// other package (e.g., PETSc), because we want to avoid problems
// with overlapping Timer scopes.
  MPI_Init(&argc, &argv);

// Get the MPI rank
  MPI_Comm_rank(MPI_COMM_WORLD, &rank);

// Print the start message on root process
  if (rank == 0) mainHelper::printStartMessage();

// Set up the material infrastructure that is used to calculate flux
  auto material = mainHelper::initMaterial(opts);
// Set up the temperature infrastructure
  bool tempInitOK = mainHelper::initTemp(opts);
  if ( !tempInitOK ) return EXIT_FAILURE;
// Set up the visualization infrastructure.
  bool vizInitOK = mainHelper::initViz( opts.useVizStandardHandlers() );
  if ( !vizInitOK ) return EXIT_FAILURE;

// Access the temperature handler registry to get the temperature
  auto tempHandler = xolotlFactory::getTemperatureHandler();

// Access our performance handler registry to obtain a Timer
// measuring the runtime of the entire program.
  auto handlerRegistry = xolotlPerf::getHandlerRegistry();
  auto totalTimer      = handlerRegistry->getTimer("total");
  totalTimer->start();

// Initialize and get the solver handler
  bool dimOK = xolotlFactory::initializeDimension( opts );
  if ( !dimOK ) return EXIT_FAILURE;
  auto solvHandler = xolotlFactory::getSolverHandler();

// Create the network handler factory
  auto networkFactory =
  xolotlFactory::IReactionHandlerFactory::createNetworkFactory( opts.getMaterial() );

// Setup and load the reaction network
  auto networkLoadTimer = handlerRegistry->getTimer("loadNetwork");
  networkLoadTimer->start();

  networkFactory->initializeReactionNetwork( opts, handlerRegistry );

  networkLoadTimer->stop();

// Get the network handler
  auto networkHandler = networkFactory->getNetworkHandler();

// Setup the solver
  auto solver = mainHelper::setUpSolver( handlerRegistry, material, tempHandler, 
                                         networkHandler, solvHandler, opts      );

// Here is the big action:
// Launch the PetscSolver
  mainHelper::launchPetscSolver( solver, handlerRegistry );

// Finalize our use of the solver.
  auto solverFinalizeTimer = handlerRegistry->getTimer("solverFinalize");
  solverFinalizeTimer->start();

  solver->finalize();

  solverFinalizeTimer->stop();

  totalTimer->stop();

// Report statistics of the performance data collected during the run just completed.
  xolotlPerf::PerfObjStatsMap<xolotlPerf::ITimer::ValType>               timerStats;
  xolotlPerf::PerfObjStatsMap<xolotlPerf::IEventCounter::ValType>        counterStats;
  xolotlPerf::PerfObjStatsMap<xolotlPerf::IHardwareCounter::CounterType> hwCtrStats;

  handlerRegistry->collectStatistics( timerStats, counterStats, hwCtrStats );

  if (rank == 0) handlerRegistry->reportStatistics( std::cout, timerStats, 
                                                    counterStats, hwCtrStats);

  } // end of try block

  catch (const std::exception& e) 
  {
   std::cerr << e.what() << std::endl;
   std::cerr << "Aborting." << std::endl;
   return EXIT_FAILURE;
  } 

  catch (const std::string& error) 
  {
   std::cout << error << std::endl;
   std::cout << "Aborting." << std::endl;
   return EXIT_FAILURE;
  }

// finalize our use of MPI
 MPI_Finalize();

 return EXIT_SUCCESS;
}
//******** End Main Program ******************************************************
