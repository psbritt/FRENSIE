//---------------------------------------------------------------------------//
//!
//! \file   tstCoupledElasticAdjointElectroatomicReaction.cpp
//! \author Luke Kersting
//! \brief  Coupled Elastic adjoint electroatomic reaction unit tests
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
#include "Data_AdjointElectronPhotonRelaxationDataContainer.hpp"
#include "MonteCarlo_CoupledElasticAdjointElectroatomicReaction.hpp"
#include "MonteCarlo_CoupledElasticElectronScatteringDistribution.hpp"
#include "MonteCarlo_ElasticElectronScatteringDistributionNativeFactory.hpp"
#include "Utility_RandomNumberGenerator.hpp"
#include "Utility_HistogramDistribution.hpp"
#include "Utility_UnitTestHarnessExtensions.hpp"
#include "Utility_TabularDistribution.hpp"

//---------------------------------------------------------------------------//
// Testing Variables.
//---------------------------------------------------------------------------//

Teuchos::RCP<MonteCarlo::CoupledElasticAdjointElectroatomicReaction<Utility::LinLin> > coupled_elastic_reaction;

//---------------------------------------------------------------------------//
// Tests
//---------------------------------------------------------------------------//
// Check that the reaction type can be returned
TEUCHOS_UNIT_TEST( CoupledElasticAdjointElectroatomicReaction, getReactionType )
{
  TEST_EQUALITY_CONST( coupled_elastic_reaction->getReactionType(),
		       MonteCarlo::COUPLED_ELASTIC_ADJOINT_ELECTROATOMIC_REACTION );
}

//---------------------------------------------------------------------------//
// Check that the threshold energy can be returned
TEUCHOS_UNIT_TEST( CoupledElasticAdjointElectroatomicReaction, getThresholdEnergy )
{
  TEST_EQUALITY_CONST( coupled_elastic_reaction->getThresholdEnergy(),
                       1.0e-5 );
}

//---------------------------------------------------------------------------//
// Check that the number of electrons emitted from the rxn can be returned
TEUCHOS_UNIT_TEST( CoupledElasticAdjointElectroatomicReaction, getNumberOfEmittedElectrons )
{
  TEST_EQUALITY_CONST( coupled_elastic_reaction->getNumberOfEmittedElectrons(1e-3),
                       0u );

  TEST_EQUALITY_CONST( coupled_elastic_reaction->getNumberOfEmittedElectrons(20.0),
                       0u );
}

//---------------------------------------------------------------------------//
// Check that the number of photons emitted from the rxn can be returned
TEUCHOS_UNIT_TEST( CoupledElasticAdjointElectroatomicReaction, getNumberOfEmittedPhotons )
{
  TEST_EQUALITY_CONST( coupled_elastic_reaction->getNumberOfEmittedPhotons(1e-3),
                       0u );

  TEST_EQUALITY_CONST( coupled_elastic_reaction->getNumberOfEmittedPhotons(20.0),
                       0u );
}

//---------------------------------------------------------------------------//
// Check that the coupled cross section can be returned
TEUCHOS_UNIT_TEST( CoupledElasticAdjointElectroatomicReaction,
                   getCrossSection )
{

  double cross_section =
    coupled_elastic_reaction->getCrossSection( 1.0E-05 );

  TEST_FLOATING_EQUALITY( cross_section, 2.74896E+08, 1e-12 );

  cross_section =
    coupled_elastic_reaction->getCrossSection( 1.0E-03 );

  TEST_FLOATING_EQUALITY( cross_section, 2.80490481543817E+06, 1e-12 );

  cross_section =
    coupled_elastic_reaction->getCrossSection( 20.0 );

  TEST_FLOATING_EQUALITY( cross_section, 1.3022122514987E+04, 1e-12 );
}

//---------------------------------------------------------------------------//
// Check that the elastic reaction can be simulated
TEUCHOS_UNIT_TEST( CoupledElasticAdjointElectroatomicReaction, react )
{
  MonteCarlo::AdjointElectronState electron( 0 );
  electron.setEnergy( 20.0 );
  electron.setDirection( 0.0, 0.0, 1.0 );

  MonteCarlo::ParticleBank bank;

  Data::SubshellType shell_of_interaction;

  coupled_elastic_reaction->react( electron, bank, shell_of_interaction );

  TEST_EQUALITY_CONST( electron.getEnergy(), 20.0 );
  TEST_ASSERT( electron.getZDirection() < 2.0 );
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
    Data::AdjointElectronPhotonRelaxationDataContainer data_container =
        Data::AdjointElectronPhotonRelaxationDataContainer( test_native_file_name );

    bool correlated_sampling_mode_on = true;
    double evaluation_tol = 1e-7;

    Teuchos::ArrayRCP<double> energy_grid;
    energy_grid.assign(
        data_container.getAdjointElectronEnergyGrid().begin(),
        data_container.getAdjointElectronEnergyGrid().end() );

    Teuchos::ArrayRCP<double> cutoff_cross_section;
    cutoff_cross_section.assign(
        data_container.getAdjointCutoffElasticCrossSection().begin(),
        data_container.getAdjointCutoffElasticCrossSection().end() );

    Teuchos::ArrayRCP<double> total_cross_section;
    total_cross_section.assign(
        data_container.getAdjointTotalElasticCrossSection().begin(),
        data_container.getAdjointTotalElasticCrossSection().end() );

    // Create coupled distribution
    std::shared_ptr<const MonteCarlo::CoupledElasticElectronScatteringDistribution>
        coupled_elastic_distribution;
    MonteCarlo::ElasticElectronScatteringDistributionNativeFactory::createCoupledElasticDistribution<Utility::LinLinLog>(
        coupled_elastic_distribution,
        energy_grid,
        cutoff_cross_section,
        total_cross_section,
        data_container,
        correlated_sampling_mode_on,
        evaluation_tol );

    // Create the reaction
    coupled_elastic_reaction.reset(
      new MonteCarlo::CoupledElasticAdjointElectroatomicReaction<Utility::LinLin>(
        energy_grid,
        total_cross_section,
        data_container.getAdjointTotalElasticCrossSectionThresholdEnergyIndex(),
        coupled_elastic_distribution ) );
  }

  // Initialize the random number generator
  Utility::RandomNumberGenerator::createStreams();
}

UTILITY_CUSTOM_TEUCHOS_UNIT_TEST_SETUP_END(); 

//---------------------------------------------------------------------------//
// end tstCoupledElasticAdjointElectroatomicReaction.cpp
//---------------------------------------------------------------------------//
