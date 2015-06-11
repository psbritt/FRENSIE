//---------------------------------------------------------------------------//
//!
//! \file   tstHardElasticElectronScatteringDistribution.cpp
//! \author Luke Kersting
//! \brief  Hard elastic electron scattering distribution unit tests
//!
//---------------------------------------------------------------------------//

// Std Lib Includes
#include <iostream>

// Trilinos Includes
#include <Teuchos_UnitTestHarness.hpp>
#include <Teuchos_VerboseObject.hpp>
#include <Teuchos_RCP.hpp>

// FRENSIE Includes
#include "MonteCarlo_HardElasticElectronScatteringDistribution.hpp"
#include "Data_ACEFileHandler.hpp"
#include "Data_XSSEPRDataExtractor.hpp"
#include "Utility_RandomNumberGenerator.hpp"
#include "Utility_HistogramDistribution.hpp"
#include "Utility_UnitTestHarnessExtensions.hpp"

//---------------------------------------------------------------------------//
// Testing Variables.
//---------------------------------------------------------------------------//

Teuchos::RCP<MonteCarlo::HardElasticElectronScatteringDistribution> 
  ace_basic_elastic_distribution;

//---------------------------------------------------------------------------//
// Tests
//---------------------------------------------------------------------------//
// Check that the screening angle can be evaluated
TEUCHOS_UNIT_TEST( HardElasticElectronScatteringDistribution, 
                   evaluateScreeningFactor )
{
  // Set energy in MeV
  double energy = 1.0;

  // Calculate scrrening angle
  double screening_factor = 
    ace_basic_elastic_distribution->evaluateScreeningFactor( energy );

  // Test
  TEST_FLOATING_EQUALITY( screening_factor, 2.195957718728E-04, 1e-12 );

}

//---------------------------------------------------------------------------//
// Check that the screened analytical function angle can be evaluated
TEUCHOS_UNIT_TEST( HardElasticElectronScatteringDistribution, 
                   evaluateScreenedScatteringAngle )
{
  // Set fake random number stream
  std::vector<double> fake_stream( 1 );
  fake_stream[0] = 0.5;

  Utility::RandomNumberGenerator::setFakeStream( fake_stream );

  // Set energy in MeV
  double energy = 1.0;

  // Calculate screening angle
  double scattering_angle_cosine = 
    ace_basic_elastic_distribution->evaluateScreenedScatteringAngle( energy );

  // Test
  TEST_FLOATING_EQUALITY( scattering_angle_cosine, 0.999999501136, 1e-12 );
}

//---------------------------------------------------------------------------//
// Check that sampleAndRecordTrialsImpl can be evaluated from the distribution
TEUCHOS_UNIT_TEST( HardElasticElectronScatteringDistribution, 
                   sampleAndRecordTrialsImpl_distribution )
{
  // Set fake random number stream
  std::vector<double> fake_stream( 2 );
  fake_stream[0] = 9.9990E-01; // sample angle from distribution
  fake_stream[1] = 0.5; // sample mu = 0.9874366113907

  Utility::RandomNumberGenerator::setFakeStream( fake_stream );

  MonteCarlo::ElectronState electron( 0 );
  electron.setEnergy( 1.0e-3 );
  electron.setDirection( 0.0, 0.0, 1.0 );

  double scattering_angle_cosine;
  unsigned trials = 10;

  // sampleAndRecordTrialsImpl from distribution
  ace_basic_elastic_distribution->sampleAndRecordTrialsImpl( 
                                                electron.getEnergy(),
                                                scattering_angle_cosine,
                                                trials );

  // Test
  TEST_FLOATING_EQUALITY( scattering_angle_cosine, 9.874339332031E-01, 1e-12 );
  TEST_EQUALITY_CONST( trials, 11 );
}
//---------------------------------------------------------------------------//
// Check that sampleAndRecordTrialsImpl can be evaluated from the screened analytical function
TEUCHOS_UNIT_TEST( HardElasticElectronScatteringDistribution, 
                   sampleAndRecordTrialsImpl_analytical )
{
  // Set fake random number stream
  std::vector<double> fake_stream( 2 );
  fake_stream[0] = 9.9991E-01;
  fake_stream[1] = 0.5;

  Utility::RandomNumberGenerator::setFakeStream( fake_stream );
  
  MonteCarlo::ElectronState electron( 0 );
  electron.setEnergy( 1.1e-3 );
  electron.setDirection( 0.0, 0.0, 1.0 );

  double scattering_angle_cosine;
  unsigned trials = 10;

  // sampleAndRecordTrialsImpl from analytical function
  ace_basic_elastic_distribution->sampleAndRecordTrialsImpl( 
                                                electron.getEnergy(),
                                                scattering_angle_cosine,
                                                trials );;

  // Test
  TEST_FLOATING_EQUALITY( scattering_angle_cosine, 0.999999500000, 1e-12 );
  TEST_EQUALITY_CONST( trials, 11 );

}

//---------------------------------------------------------------------------//
// Check sample can be evaluated from the distribution
TEUCHOS_UNIT_TEST( HardElasticElectronScatteringDistribution, 
                   sample_distribution )
{
  // Set fake random number stream
  std::vector<double> fake_stream( 2 );
  fake_stream[0] = 9.9990E-01; // sample angle from distribution
  fake_stream[1] = 0.5; // sample mu = 0.9874366113907

  Utility::RandomNumberGenerator::setFakeStream( fake_stream );

  MonteCarlo::ElectronState electron( 0 );
  electron.setEnergy( 1.0e-3 );
  electron.setDirection( 0.0, 0.0, 1.0 );

  double scattering_angle_cosine, outgoing_energy;

  // sampleAndRecordTrialsImpl from distribution
  ace_basic_elastic_distribution->sample( electron.getEnergy(),
                                          outgoing_energy,
                                          scattering_angle_cosine );

  // Test
  TEST_FLOATING_EQUALITY( scattering_angle_cosine, 9.874339332031E-01, 1e-12 );
  TEST_FLOATING_EQUALITY( outgoing_energy, 1.0e-3, 1e-12 );
}
//---------------------------------------------------------------------------//
// Check that sample can be evaluated from the screened analytical function
TEUCHOS_UNIT_TEST( HardElasticElectronScatteringDistribution, 
                   sample_analytical )
{
  // Set fake random number stream
  std::vector<double> fake_stream( 2 );
  fake_stream[0] = 9.9991E-01;
  fake_stream[1] = 0.5;

  Utility::RandomNumberGenerator::setFakeStream( fake_stream );
  
  MonteCarlo::ElectronState electron( 0 );
  electron.setEnergy( 1.1e-3 );
  electron.setDirection( 0.0, 0.0, 1.0 );

  double scattering_angle_cosine, outgoing_energy;

  // sampleAndRecordTrialsImpl from analytical function
  ace_basic_elastic_distribution->sample( electron.getEnergy(),
                                          outgoing_energy,
                                          scattering_angle_cosine );

  // Test
  TEST_FLOATING_EQUALITY( scattering_angle_cosine, 0.999999500000, 1e-12 );
  TEST_FLOATING_EQUALITY( outgoing_energy, 1.1e-3, 1e-12 );

}
//---------------------------------------------------------------------------//
// Check sampleAndRecordTrials can be evaluated from the distribution
TEUCHOS_UNIT_TEST( HardElasticElectronScatteringDistribution, 
                   sampleAndRecordTrials_distribution )
{
  // Set fake random number stream
  std::vector<double> fake_stream( 2 );
  fake_stream[0] = 9.9990E-01; // sample angle from distribution
  fake_stream[1] = 0.5; // sample mu = 0.9874366113907

  Utility::RandomNumberGenerator::setFakeStream( fake_stream );

  MonteCarlo::ElectronState electron( 0 );
  electron.setEnergy( 1.0e-3 );
  electron.setDirection( 0.0, 0.0, 1.0 );

  double scattering_angle_cosine, outgoing_energy;
  unsigned trials = 10;

  // sampleAndRecordTrialsImpl from distribution
  ace_basic_elastic_distribution->sampleAndRecordTrials( 
                                          electron.getEnergy(),
                                          outgoing_energy,
                                          scattering_angle_cosine,
                                          trials );

  // Test
  TEST_FLOATING_EQUALITY( scattering_angle_cosine, 9.874339332031E-01, 1e-12 );
  TEST_FLOATING_EQUALITY( outgoing_energy, 1.0e-3, 1e-12 );
}
//---------------------------------------------------------------------------//
// Check that sampleAndRecordTrials can be evaluated from the screened analytical function
TEUCHOS_UNIT_TEST( HardElasticElectronScatteringDistribution, 
                   sampleAndRecordTrials_analytical )
{
  // Set fake random number stream
  std::vector<double> fake_stream( 2 );
  fake_stream[0] = 9.9991E-01;
  fake_stream[1] = 0.5;

  Utility::RandomNumberGenerator::setFakeStream( fake_stream );
  
  MonteCarlo::ElectronState electron( 0 );
  electron.setEnergy( 1.1e-3 );
  electron.setDirection( 0.0, 0.0, 1.0 );

  double scattering_angle_cosine, outgoing_energy;
  unsigned trials = 10;

  // sampleAndRecordTrialsImpl from analytical function
  ace_basic_elastic_distribution->sampleAndRecordTrials( 
                                          electron.getEnergy(),
                                          outgoing_energy,
                                          scattering_angle_cosine,
                                          trials );

  // Test
  TEST_FLOATING_EQUALITY( scattering_angle_cosine, 0.999999500000, 1e-12 );
  TEST_FLOATING_EQUALITY( outgoing_energy, 1.1e-3, 1e-12 );

}
//---------------------------------------------------------------------------//
// Check that the angle can be evaluated from the distribution
TEUCHOS_UNIT_TEST( HardElasticElectronScatteringDistribution, 
                   ScatterElectron_distribution )
{
  // Set fake random number stream
  std::vector<double> fake_stream( 2 );
  fake_stream[0] = 9.9990E-01; // sample angle from distribution
  fake_stream[1] = 0.5; // sample mu = 0.9874366113907

  Utility::RandomNumberGenerator::setFakeStream( fake_stream );

  MonteCarlo::ParticleBank bank;
  MonteCarlo::SubshellType shell_of_interaction;
  
  MonteCarlo::ElectronState electron( 0 );
  electron.setEnergy( 1.0e-3 );
  electron.setDirection( 0.0, 0.0, 1.0 );

  // Analytically scatter electron
  ace_basic_elastic_distribution->scatterElectron( electron, 
                                                   bank, 
                                                   shell_of_interaction );

  // Test
  TEST_FLOATING_EQUALITY( electron.getZDirection(), 9.874339332031E-01, 1e-12 );
  TEST_FLOATING_EQUALITY( electron.getEnergy(), 1.0e-3, 1e-12 );

}
//---------------------------------------------------------------------------//
// Check that the angle can be evaluated from the screened analytical function
TEUCHOS_UNIT_TEST( HardElasticElectronScatteringDistribution, 
                   ScatterElectron_analytical )
{
  // Set fake random number stream
  std::vector<double> fake_stream( 2 );
  fake_stream[0] = 9.9991E-01;
  fake_stream[1] = 0.5;

  Utility::RandomNumberGenerator::setFakeStream( fake_stream );

  MonteCarlo::ParticleBank bank;
  MonteCarlo::SubshellType shell_of_interaction;
  
  MonteCarlo::ElectronState electron( 0 );
  electron.setEnergy( 1.1e-3 );
  electron.setDirection( 0.0, 0.0, 1.0 );

  // Analytically scatter electron
  ace_basic_elastic_distribution->scatterElectron( electron, 
                                                   bank, 
                                                   shell_of_interaction );

  // Test
  TEST_FLOATING_EQUALITY( electron.getZDirection(), 0.999999500000, 1e-12 );
  TEST_FLOATING_EQUALITY( electron.getEnergy(), 1.1e-3, 1e-12 );

}
//---------------------------------------------------------------------------//
// Check that the angle can be evaluated from the distribution
TEUCHOS_UNIT_TEST( HardElasticElectronScatteringDistribution, 
                   ScatterAdjointElectron_distribution )
{
  // Set fake random number stream
  std::vector<double> fake_stream( 2 );
  fake_stream[0] = 9.9990E-01; // sample angle from distribution
  fake_stream[1] = 0.5; // sample mu = 0.9874366113907

  Utility::RandomNumberGenerator::setFakeStream( fake_stream );

  MonteCarlo::ParticleBank bank;
  MonteCarlo::SubshellType shell_of_interaction;
  
  MonteCarlo::AdjointElectronState adjoint_electron( 0 );
  adjoint_electron.setEnergy( 1.0e-3 );
  adjoint_electron.setDirection( 0.0, 0.0, 1.0 );

  // Analytically scatter electron
  ace_basic_elastic_distribution->scatterAdjointElectron( 
                                                   adjoint_electron, 
                                                   bank, 
                                                   shell_of_interaction );

  // Test
  TEST_FLOATING_EQUALITY( adjoint_electron.getZDirection(), 
                          9.874339332031E-01, 
                          1e-12 );
  TEST_FLOATING_EQUALITY( adjoint_electron.getEnergy(), 1.0e-3, 1e-12 );

}
//---------------------------------------------------------------------------//
// Check that the angle can be evaluated from the screened analytical function
TEUCHOS_UNIT_TEST( HardElasticElectronScatteringDistribution, 
                   ScatterAdjointElectron_analytical )
{
  // Set fake random number stream
  std::vector<double> fake_stream( 2 );
  fake_stream[0] = 9.9991E-01;
  fake_stream[1] = 0.5;

  Utility::RandomNumberGenerator::setFakeStream( fake_stream );

  MonteCarlo::ParticleBank bank;
  MonteCarlo::SubshellType shell_of_interaction;
  
  MonteCarlo::AdjointElectronState adjoint_electron( 0 );
  adjoint_electron.setEnergy( 1.1e-3 );
  adjoint_electron.setDirection( 0.0, 0.0, 1.0 );

  // Analytically scatter adjoint electron
  ace_basic_elastic_distribution->scatterAdjointElectron( 
                                                   adjoint_electron, 
                                                   bank, 
                                                   shell_of_interaction );

  // Test
  TEST_FLOATING_EQUALITY( adjoint_electron.getZDirection(), 
                          0.999999500000, 
                          1e-12 );
  TEST_FLOATING_EQUALITY( adjoint_electron.getEnergy(), 1.1e-3, 1e-12 );

}

//---------------------------------------------------------------------------//
// Custom main function
//---------------------------------------------------------------------------//
int main( int argc, char** argv )
{
  std::string test_ace_file_name, test_ace_table_name;
  
  Teuchos::CommandLineProcessor& clp = Teuchos::UnitTestRepository::getCLP();

  clp.setOption( "test_ace_file",
		 &test_ace_file_name,
		 "Test ACE file name" );
  clp.setOption( "test_ace_table",
		 &test_ace_table_name,
		 "Test ACE table name" );

  const Teuchos::RCP<Teuchos::FancyOStream> out = 
    Teuchos::VerboseObjectBase::getDefaultOStream();

  Teuchos::CommandLineProcessor::EParseCommandLineReturn parse_return = 
    clp.parse(argc,argv);

  if ( parse_return != Teuchos::CommandLineProcessor::PARSE_SUCCESSFUL ) {
    *out << "\nEnd Result: TEST FAILED" << std::endl;
    return parse_return;
  }
  
  // Create a file handler and data extractor
  Teuchos::RCP<Data::ACEFileHandler> ace_file_handler( 
				 new Data::ACEFileHandler( test_ace_file_name,
							   test_ace_table_name,
							   1u ) );
  Teuchos::RCP<Data::XSSEPRDataExtractor> xss_data_extractor(
                            new Data::XSSEPRDataExtractor( 
				      ace_file_handler->getTableNXSArray(),
				      ace_file_handler->getTableJXSArray(),
				      ace_file_handler->getTableXSSArray() ) );

  // Extract the elastic scattering information data block (ELASI)
  Teuchos::ArrayView<const double> elasi_block(
				      xss_data_extractor->extractELASIBlock() );
  
  // Extract the number of tabulated distributions
  int size = elasi_block.size()/3;

  // Extract the energy grid for elastic scattering angular distributions
  Teuchos::Array<double> energy_grid(elasi_block(0,size));

  // Extract the table lengths for elastic scattering angular distributions
  Teuchos::Array<double> table_length(elasi_block(size,size));

  // Extract the offsets for elastic scattering angular distributions
  Teuchos::Array<double> offset(elasi_block(2*size,size));

  // Extract the elastic scattering angular distributions block (elas)
  Teuchos::ArrayView<const double> elas_block = 
    xss_data_extractor->extractELASBlock();

  // Create the elastic scattering distributions
  Teuchos::Array<Utility::Pair<double,Teuchos::RCP<const Utility::TabularOneDDistribution> > >
    elastic_scattering_distribution( size );
  
  for( unsigned n = 0; n < size; ++n )
  {
    elastic_scattering_distribution[n].first = energy_grid[n];

    elastic_scattering_distribution[n].second.reset( 
	  new Utility::HistogramDistribution(
		 elas_block( offset[n], table_length[n] ),
		 elas_block( offset[n] + 1 + table_length[n], table_length[n]-1 ),
         true ) );
/*	  new Utility::TabularDistribution<Utility::LinLin>(
		 elas_block( offset[n], table_length[n] ),
		 elas_block( offset[n] + 1 + table_length[n], table_length[n] ) ) );
*/
  }  

  // Get the atomic number 
  const int atomic_number = xss_data_extractor->extractAtomicNumber();

  // Create the distributions
  ace_basic_elastic_distribution.reset(
		new MonteCarlo::HardElasticElectronScatteringDistribution(
						    atomic_number,
						    elastic_scattering_distribution ) );

  // Clear setup data
  ace_file_handler.reset();
  xss_data_extractor.reset();

  // Initialize the random number generator
  Utility::RandomNumberGenerator::createStreams();

  // Run the unit tests
  Teuchos::GlobalMPISession mpiSession( &argc, &argv );

  const bool success = Teuchos::UnitTestRepository::runUnitTests( *out );

  if (success)
    *out << "\nEnd Result: TEST PASSED" << std::endl;
  else
    *out << "\nEnd Result: TEST FAILED" << std::endl;

  clp.printFinalTimerSummary(out.ptr());

  return (success ? 0 : 1);  					    
}

//---------------------------------------------------------------------------//
// end tstHardElasticElectronScatteringDistribution.cpp
//---------------------------------------------------------------------------//
