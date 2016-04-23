//---------------------------------------------------------------------------//
//!
//! \file   tstRandomNumberGenerator.cpp
//! \author Alex Robinson
//! \brief  Random number generator class unit tests
//!
//---------------------------------------------------------------------------//

// Std Lib Includes
#include <iostream>
#include <set>

// Trilinos Includes
#include <Teuchos_UnitTestHarness.hpp>
#include <Teuchos_ScalarTraits.hpp>
#include <Teuchos_VerboseObject.hpp>

// FRENSIE Includes
#include "Utility_UnitTestHarnessExtensions.hpp"
#include "Utility_RandomNumberGenerator.hpp"
#include "Utility_GlobalOpenMPSession.hpp"

//---------------------------------------------------------------------------//
// Instantiation Macros.
//---------------------------------------------------------------------------//
#define UNIT_TEST_INSTANTIATION( type, name ) \
  TEUCHOS_UNIT_TEST_TEMPLATE_1_INSTANT( type, name, float ) \
  TEUCHOS_UNIT_TEST_TEMPLATE_1_INSTANT( type, name, double ) 

//---------------------------------------------------------------------------//
// Tests.
//---------------------------------------------------------------------------//
// Check that the random number generator streams can be initialized
TEUCHOS_UNIT_TEST( RandomNumberGenerator, createStreams )
{
  TEST_ASSERT( !Utility::RandomNumberGenerator::hasStreams() );

  // Initialize the random number generator
  Utility::RandomNumberGenerator::createStreams();

  TEST_ASSERT( Utility::RandomNumberGenerator::hasStreams() );
}

//---------------------------------------------------------------------------//
// Check that the random number generator can be initialized
TEUCHOS_UNIT_TEST_TEMPLATE_1_DECL( RandomNumberGenerator,
				   initialize,
				   ScalarType )
{
  Utility::RandomNumberGenerator::initialize();

  // An exception will be thrown if the initialization failed
  double random_number = 
    Utility::RandomNumberGenerator::getRandomNumber<ScalarType>();
}

UNIT_TEST_INSTANTIATION( RandomNumberGenerator, initialize );

//---------------------------------------------------------------------------//
// Check that a fake stream can be set
TEUCHOS_UNIT_TEST_TEMPLATE_1_DECL( RandomNumberGenerator,
				   setFakeStream,
				   ScalarType )
{
  // Create the fake stream
  std::vector<double> fake_stream( 3 );
  fake_stream[0] = 0.2;
  fake_stream[1] = 0.4;
  fake_stream[2] = 0.6;

  // Set the fake stream
  Utility::RandomNumberGenerator::setFakeStream( fake_stream );

  // Test the the fake stream returns the correct values
  ScalarType random_number = 
    Utility::RandomNumberGenerator::getRandomNumber<ScalarType>();
  TEST_EQUALITY_CONST( random_number, (ScalarType)0.2 );
  
  random_number = 
    Utility::RandomNumberGenerator::getRandomNumber<ScalarType>();
  
  TEST_EQUALITY_CONST( random_number, (ScalarType)0.4 );

  random_number = 
    Utility::RandomNumberGenerator::getRandomNumber<ScalarType>();
  
  TEST_EQUALITY_CONST( random_number, (ScalarType)0.6 );

  // Unset the fake stream
  Utility::RandomNumberGenerator::unsetFakeStream();
}

UNIT_TEST_INSTANTIATION( RandomNumberGenerator, setFakeStream );

//---------------------------------------------------------------------------//
// Check that the random number generator can be initialized to a new history
TEUCHOS_UNIT_TEST( RandomNumberGenerator, initialize_history )
{
  Teuchos::Array<unsigned long long> local_random_numbers( 
		 Utility::GlobalOpenMPSession::getRequestedNumberOfThreads() );
  
  // Initialize the generator to a particular history depending on the process
#pragma omp parallel num_threads( Utility::GlobalOpenMPSession::getRequestedNumberOfThreads() )
  {
    unsigned history_number = Teuchos::GlobalMPISession::getRank()*
      Utility::GlobalOpenMPSession::getRequestedNumberOfThreads() +
      Utility::GlobalOpenMPSession::getThreadId();
    
    Utility::RandomNumberGenerator::initialize( history_number );
    
    // Generate a random number
    local_random_numbers[Utility::GlobalOpenMPSession::getThreadId()] = 
      Utility::RandomNumberGenerator::getRandomNumber<unsigned long long>();
  }

  // Retrieve the random numbers generated by the other processes and store
  // them in an array
  Teuchos::Array<int> all_random_numbers( 
		 Teuchos::GlobalMPISession::getNProc()*
		 Utility::GlobalOpenMPSession::getRequestedNumberOfThreads() );
  
  for( unsigned i = 0; i < local_random_numbers.size(); ++i )
  {
    unsigned lower_bound = 
      i*Teuchos::GlobalMPISession::getNProc();

    unsigned length = Teuchos::GlobalMPISession::getNProc();
    
    Teuchos::GlobalMPISession::allGather( 
				      local_random_numbers[i],
			              all_random_numbers(lower_bound,length) );
  }
  
  Teuchos::GlobalMPISession::barrier();
  
  std::cout << local_random_numbers << std::endl;

  if( Teuchos::GlobalMPISession::getRank() == 0 )
    std::cout << all_random_numbers << std::endl;
  
  // Store all of the array elements in a set
  std::set<int> random_set;

  for( int i = 0; i < all_random_numbers.size(); ++i )
    random_set.insert( all_random_numbers[i] );

  TEST_EQUALITY( all_random_numbers.size(), random_set.size() );
}

//---------------------------------------------------------------------------//
// Custom main function
//---------------------------------------------------------------------------//
int main( int argc, char** argv )
{
  Teuchos::CommandLineProcessor& clp = Teuchos::UnitTestRepository::getCLP();
 
  int threads = 1;
 
  clp.setOption( "threads",
		 &threads,
		 "Number of threads to use" );
  
  Teuchos::CommandLineProcessor::EParseCommandLineReturn parse_return = 
    clp.parse(argc,argv);

  const Teuchos::RCP<Teuchos::FancyOStream> out = Teuchos::VerboseObjectBase::getDefaultOStream();

  if ( parse_return != Teuchos::CommandLineProcessor::PARSE_SUCCESSFUL ) {
    *out << "\nEnd Result: TEST FAILED" << std::endl;
    return parse_return;
  }

  // Set up the global OpenMP session
  if( Utility::GlobalOpenMPSession::isOpenMPUsed() )
    Utility::GlobalOpenMPSession::setNumberOfThreads( threads );

  // Run the unit tests
  Teuchos::GlobalMPISession mpiSession( &argc, &argv );

  out->setProcRankAndSize( mpiSession.getRank(), mpiSession.getNProc() );
  out->setOutputToRootOnly( 0 );

  const bool success = Teuchos::UnitTestRepository::runUnitTests(*out);

  if (success)
    *out << "\nEnd Result: TEST PASSED" << std::endl;
  else
    *out << "\nEnd Result: TEST FAILED" << std::endl;

  clp.printFinalTimerSummary(out.ptr());

  return (success ? 0 : 1);
}

//---------------------------------------------------------------------------//
// end tstRandomNumberGenerator.cpp
//---------------------------------------------------------------------------//

