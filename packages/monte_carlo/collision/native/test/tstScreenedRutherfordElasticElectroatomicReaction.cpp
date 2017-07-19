//---------------------------------------------------------------------------//
//!
//! \file   tstScreenedRutherfordElasticElectroatomicReaction.cpp
//! \author Luke Kersting
//! \brief  Moment preserving Elastic electroatomic reaction unit tests
//!
//---------------------------------------------------------------------------//

// Std Lib Includes
#include <iostream>

// Trilinos Includes
#include <Teuchos_UnitTestHarness.hpp>
#include <Teuchos_VerboseObject.hpp>
#include <Teuchos_RCP.hpp>
#include <Teuchos_ArrayRCP.hpp>

// FRENSIE Includes
#include "Data_ElectronPhotonRelaxationDataContainer.hpp"
#include "MonteCarlo_ScreenedRutherfordElasticElectroatomicReaction.hpp"
#include "MonteCarlo_ElasticElectronScatteringDistributionNativeFactory.hpp"
#include "Utility_RandomNumberGenerator.hpp"
#include "Utility_HistogramDistribution.hpp"
#include "Utility_UnitTestHarnessExtensions.hpp"
#include "Utility_TabularOneDDistribution.hpp"

typedef MonteCarlo::ElasticElectronScatteringDistributionNativeFactory 
    NativeFactory;

//---------------------------------------------------------------------------//
// Testing Variables.
//---------------------------------------------------------------------------//

std::shared_ptr<MonteCarlo::ElectroatomicReaction>
    rutherford_elastic_reaction;

//---------------------------------------------------------------------------//
// Tests
//---------------------------------------------------------------------------//
// Check that the reaction type can be returned
TEUCHOS_UNIT_TEST( ScreenedRutherfordElasticElectroatomicReaction, getReactionType )
{
  TEST_EQUALITY_CONST( rutherford_elastic_reaction->getReactionType(),
		       MonteCarlo::SCREENED_RUTHERFORD_ELASTIC_ELECTROATOMIC_REACTION );
}

//---------------------------------------------------------------------------//
// Check that the threshold energy can be returned
TEUCHOS_UNIT_TEST( ScreenedRutherfordElasticElectroatomicReaction, getThresholdEnergy )
{
  TEST_EQUALITY_CONST( rutherford_elastic_reaction->getThresholdEnergy(),
                       1.75 );
}

//---------------------------------------------------------------------------//
// Check that the number of electrons emitted from the rxn can be returned
TEUCHOS_UNIT_TEST( ScreenedRutherfordElasticElectroatomicReaction, getNumberOfEmittedElectrons )
{
  TEST_EQUALITY_CONST( rutherford_elastic_reaction->getNumberOfEmittedElectrons(1e-3),
		       0u );

  TEST_EQUALITY_CONST( rutherford_elastic_reaction->getNumberOfEmittedElectrons(20.0),
		       0u );
}

//---------------------------------------------------------------------------//
// Check that the number of photons emitted from the rxn can be returned
TEUCHOS_UNIT_TEST( ScreenedRutherfordElasticElectroatomicReaction, getNumberOfEmittedPhotons )
{
  TEST_EQUALITY_CONST( rutherford_elastic_reaction->getNumberOfEmittedPhotons(1e-3),
		       0u );

  TEST_EQUALITY_CONST( rutherford_elastic_reaction->getNumberOfEmittedPhotons(20.0),
		       0u );
}

//---------------------------------------------------------------------------//
// Check that the moment preserving cross section can be returned
TEUCHOS_UNIT_TEST( ScreenedRutherfordElasticElectroatomicReaction,
                   getCrossSection )
{

  double cross_section =
    rutherford_elastic_reaction->getCrossSection( 1.0E-05 );

  TEST_FLOATING_EQUALITY( cross_section, 0.0, 1e-12 );

  cross_section =
    rutherford_elastic_reaction->getCrossSection( 1.0E-03 );

  TEST_FLOATING_EQUALITY( cross_section, 0.0, 1e-12 );

  cross_section =
    rutherford_elastic_reaction->getCrossSection( 1.75 );

  TEST_FLOATING_EQUALITY( cross_section, 5.822907080362080E+03, 1e-12 );

  cross_section =
    rutherford_elastic_reaction->getCrossSection( 1.0E+05 );

  TEST_FLOATING_EQUALITY( cross_section, 3.855069977822990E+05, 1e-12 );
}


//---------------------------------------------------------------------------//
// Check that the elastic reaction can be simulated
TEUCHOS_UNIT_TEST( ScreenedRutherfordElasticElectroatomicReaction, react )
{
  MonteCarlo::ElectronState electron( 0 );
  electron.setEnergy( 20.0 );
  electron.setDirection( 0.0, 0.0, 1.0 );

  MonteCarlo::ParticleBank bank;

  Data::SubshellType shell_of_interaction;

  rutherford_elastic_reaction->react( electron, bank, shell_of_interaction );

  TEST_EQUALITY_CONST( electron.getEnergy(), 20.0 );
  TEST_ASSERT( electron.getZDirection() < 1.0 );
  TEST_ASSERT( electron.getZDirection() > 0.0 );
  TEST_ASSERT( bank.isEmpty() );
  TEST_EQUALITY_CONST( shell_of_interaction, Data::UNKNOWN_SUBSHELL );
}

//---------------------------------------------------------------------------//
// Custom setup
//---------------------------------------------------------------------------//
UTILITY_CUSTOM_TEUCHOS_UNIT_TEST_SETUP_BEGIN();

std::string test_native_file_name;

UTILITY_CUSTOM_TEUCHOS_UNIT_TEST_COMMAND_LINE_OPTIONS()
{
  clp().setOption( "test_native_file",
                   &test_native_file_name,
                   "Test Native file name" );
}

UTILITY_CUSTOM_TEUCHOS_UNIT_TEST_DATA_INITIALIZATION()
{
  // Create reaction
  {
    // Get native data container
    Data::ElectronPhotonRelaxationDataContainer data_container =
        Data::ElectronPhotonRelaxationDataContainer( test_native_file_name );

    // Create cutoff distribution
    std::shared_ptr<const MonteCarlo::CutoffElasticElectronScatteringDistribution>
        cutoff_elastic_distribution;

    double cutoff_angle_cosine = 0.9;
    double evaluation_tol = 1e-7;
    bool correlated_sampling_mode_on = true;

    NativeFactory::createCutoffElasticDistribution(
        cutoff_elastic_distribution,
        data_container,
        cutoff_angle_cosine,
        correlated_sampling_mode_on,
        evaluation_tol );

    // Create the screened rutherford distribution
    std::shared_ptr<const MonteCarlo::ScreenedRutherfordElasticElectronScatteringDistribution>
        rutherford_elastic_distribution;

    NativeFactory::createScreenedRutherfordElasticDistribution(
        rutherford_elastic_distribution,
        cutoff_elastic_distribution,
        data_container.getAtomicNumber() );

    Teuchos::ArrayRCP<double> energy_grid;
    energy_grid.assign(
        data_container.getElectronEnergyGrid().begin(),
        data_container.getElectronEnergyGrid().end() );

    Teuchos::ArrayRCP<double> cross_section;
    cross_section.assign(
        data_container.getScreenedRutherfordElasticCrossSection().begin(),
        data_container.getScreenedRutherfordElasticCrossSection().end() );

    unsigned threshold_index(
        data_container.getScreenedRutherfordElasticCrossSectionThresholdEnergyIndex() );

    // Create the reaction
    rutherford_elastic_reaction.reset(
      new MonteCarlo::ScreenedRutherfordElasticElectroatomicReaction<Utility::LinLin>(
                energy_grid,
                cross_section,
                threshold_index,
                rutherford_elastic_distribution ) );
  }

  // Initialize the random number generator
  Utility::RandomNumberGenerator::createStreams();
}

UTILITY_CUSTOM_TEUCHOS_UNIT_TEST_SETUP_END();

//---------------------------------------------------------------------------//
// end tstCutoffElasticElectroatomicReaction.cpp
//---------------------------------------------------------------------------//