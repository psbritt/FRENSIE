//---------------------------------------------------------------------------//
//!
//! \file   tstAbsorptionAdjointElectroatomicReaction.cpp
//! \author Luke Kersting
//! \brief  Absorption adjoint electroatomic reaction unit tests
//!
//---------------------------------------------------------------------------//

// Std Lib Includes
#include <iostream>

// Trilinos Includes
#include <Teuchos_UnitTestHarness.hpp>
#include <Teuchos_VerboseObject.hpp>
#include <Teuchos_RCP.hpp>

// FRENSIE Includes
#include "MonteCarlo_AbsorptionAdjointElectroatomicReaction.hpp"
#include "Data_AdjointElectronPhotonRelaxationDataContainer.hpp"
#include "Utility_UnitTestHarnessExtensions.hpp"

//---------------------------------------------------------------------------//
// Testing Variables
//---------------------------------------------------------------------------//

std::shared_ptr<MonteCarlo::AdjointElectroatomicReaction> absorption_reaction;

//---------------------------------------------------------------------------//
// Tests
//---------------------------------------------------------------------------//
// Check that the reaction type can be returned
TEUCHOS_UNIT_TEST( AbsorptionAdjointElectroatomicReaction, getReactionType )
{
  TEST_EQUALITY_CONST( absorption_reaction->getReactionType(),
                       MonteCarlo::TOTAL_ADJOINT_ELECTROATOMIC_REACTION );
}

//---------------------------------------------------------------------------//
// Check that the threshold energy can be returned
TEUCHOS_UNIT_TEST( AbsorptionAdjointElectroatomicReaction, getThresholdEnergy )
{
  TEST_EQUALITY_CONST( absorption_reaction->getThresholdEnergy(),
                       1e-5 );
}

//---------------------------------------------------------------------------//
// Check that the number of electrons emitted from the rxn can be returned
TEUCHOS_UNIT_TEST( AbsorptionAdjointElectroatomicReaction,
                   getNumberOfEmittedElectrons )
{
  TEST_EQUALITY_CONST( absorption_reaction->getNumberOfEmittedElectrons( 1e-3 ),
                       0u );

  TEST_EQUALITY_CONST( absorption_reaction->getNumberOfEmittedElectrons( 20.0 ),
                       0u );
}

//---------------------------------------------------------------------------//
// Check that the number of photons emitted from the rxn can be returned
TEUCHOS_UNIT_TEST( AbsorptionAdjointElectroatomicReaction,
                   getNumberOfEmittedPhotons )
{
  TEST_EQUALITY_CONST( absorption_reaction->getNumberOfEmittedPhotons( 1e-3 ),
                       0u );

  TEST_EQUALITY_CONST( absorption_reaction->getNumberOfEmittedPhotons( 20.0 ),
                       0u );
}

//---------------------------------------------------------------------------//
// Check that the cross section can be returned
TEUCHOS_UNIT_TEST( AbsorptionAdjointElectroatomicReaction, getCrossSection )
{
  double cross_section = absorption_reaction->getCrossSection( 1e-5 );
  TEST_FLOATING_EQUALITY( cross_section, 4.6179443997604473e+01, 1e-12 );

  cross_section = absorption_reaction->getCrossSection( 2e-2 );
  TEST_FLOATING_EQUALITY( cross_section, 1.8837112085990035, 1e-12 );

  cross_section = absorption_reaction->getCrossSection( 20.0 );
  TEST_FLOATING_EQUALITY( cross_section, 7.7113235533702451e-01, 1e-12 );
}

//---------------------------------------------------------------------------//
// Check that the absorption reaction can be simulated
TEUCHOS_UNIT_TEST( AbsorptionAdjointElectroatomicReaction, react )
{
  MonteCarlo::AdjointElectronState electron( 0 );
  electron.setEnergy( 20.0 );
  electron.setDirection( 0.0, 0.0, 1.0 );

  MonteCarlo::ParticleBank bank;

  Data::SubshellType shell_of_interaction;

  absorption_reaction->react( electron, bank, shell_of_interaction );

  TEST_ASSERT( electron.isGone() );
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
  // Create Native Reaction
  {
    // Get native data container
    Data::AdjointElectronPhotonRelaxationDataContainer data_container =
      Data::AdjointElectronPhotonRelaxationDataContainer( test_native_file_name );

    // Get the energy grid
    Teuchos::ArrayRCP<double> energy_grid;
    energy_grid.assign(
        data_container.getAdjointElectronEnergyGrid().begin(),
        data_container.getAdjointElectronEnergyGrid().end() );

    // Get the cross section (use the brem cross sections as a filler)
    Teuchos::ArrayRCP<double> cross_section;
    cross_section.assign(
        data_container.getAdjointBremsstrahlungElectronCrossSection().begin(),
        data_container.getAdjointBremsstrahlungElectronCrossSection().end() );

    // Create the reaction
    absorption_reaction.reset(
      new MonteCarlo::AbsorptionAdjointElectroatomicReaction<Utility::LinLin>(
        energy_grid,
        cross_section,
        data_container.getAdjointBremsstrahlungElectronCrossSectionThresholdEnergyIndex(),
        MonteCarlo::TOTAL_ADJOINT_ELECTROATOMIC_REACTION ) );
  }
}

UTILITY_CUSTOM_TEUCHOS_UNIT_TEST_SETUP_END();

//---------------------------------------------------------------------------//
// end tstAbsorptionAdjointElectroatomicReaction.cpp
//---------------------------------------------------------------------------//
