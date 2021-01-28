//---------------------------------------------------------------------------//
//!
//! \file   tstGenericHistogramImportanceParticleDistribution.cpp
//! \author Philip Britt
//! \brief  Generic histogram importance particle distribution unit tests
//!
//---------------------------------------------------------------------------//

// Std Lib Includes
#include <memory>

// FRENSIE Includes
#include "MonteCarlo_GenericHistogramImportanceParticleDistribution.hpp"
#include "MonteCarlo_IndependentPhaseSpaceDimensionDistribution.hpp"
#include "MonteCarlo_ImportanceSampledIndependentPhaseSpaceDimensionDistribution.hpp"
#include "MonteCarlo_PhotonState.hpp"
#include "Utility_BasicCartesianCoordinateConversionPolicy.hpp"
#include "Utility_BasicSphericalCoordinateConversionPolicy.hpp"
#include "Utility_HistogramDistribution.hpp"
#include "Utility_DeltaDistribution.hpp"
#include "Utility_RandomNumberGenerator.hpp"
#include "Utility_PhysicalConstants.hpp"
#include "Utility_UnitTestHarnessWithMain.hpp"
#include "ArchiveTestHelpers.hpp"
#include "MonteCarlo_PhaseSpaceDimension.hpp"

//---------------------------------------------------------------------------//
// Testing Types
//---------------------------------------------------------------------------//

typedef TestArchiveHelper::TestArchives TestArchives;

//---------------------------------------------------------------------------//
// Testing Variables
//---------------------------------------------------------------------------//

std::shared_ptr< MonteCarlo::GenericHistogramImportanceParticleDistribution > distribution;
double number_of_mesh_elements = 2;
double number_of_direction_elements = 8*4;
std::shared_ptr< Utility::DeltaDistribution> time_distribution;
std::shared_ptr< Utility::HistogramDistribution> mesh_distribution;
std::shared_ptr< Utility::HistogramDistribution> direction_distribution;
std::shared_ptr< Utility::HistogramDistribution> energy_distribution;
std::map< MonteCarlo::PhaseSpaceDimension, std::vector< double > > importance_distribution_boundaries;
std::vector< MonteCarlo::PhaseSpaceDimension > phase_space_dimension_importance_order;

std::shared_ptr< MonteCarlo::PhaseSpaceDimensionDistribution > independent_time_distribution;
std::map<  MonteCarlo::PhaseSpaceDimension, std::vector< std::shared_ptr<  MonteCarlo::PhaseSpaceDimensionDistribution> > > importance_dimension_map;

//---------------------------------------------------------------------------//
// Tests
//---------------------------------------------------------------------------//

// Test that it constructs
FRENSIE_UNIT_TEST( GenericHistogramImportanceParticleDistribution, constructor )
{
  distribution = std::make_shared<MonteCarlo::GenericHistogramImportanceParticleDistribution>("Generic importance sampled test distribution");
}

// Test that it sets independent non-importance sampled distributions properly
FRENSIE_UNIT_TEST( GenericHistogramImportanceParticleDistribution, setIndependentDimensionDistribution)
{
  distribution->setIndependentDimensionDistribution(independent_time_distribution, true);

  MonteCarlo::PhotonState photon( 0 );
  distribution->sample(photon);

  FRENSIE_CHECK_EQUAL(photon.getWeight(), 1.0);
  FRENSIE_CHECK_EQUAL(photon.getTime(), 3.0);
}

FRENSIE_UNIT_TEST( GenericHistogramImportanceParticleDistribution, setImportanceDimensionDistributions)
{
  distribution->setImportanceDimensionDistributions(importance_dimension_map,
                                                    importance_distribution_boundaries,
                                                    phase_space_dimension_importance_order);
  

}

//---------------------------------------------------------------------------//
// Custom setup
//---------------------------------------------------------------------------//
FRENSIE_CUSTOM_UNIT_TEST_SETUP_BEGIN();

FRENSIE_CUSTOM_UNIT_TEST_INIT()
{
  distribution = std::make_shared< MonteCarlo::GenericHistogramImportanceParticleDistribution >("non-archive test distribution");

  std::vector<double> x_planes = {0.0, 1.0, 2.0};
  std::vector<double> y_planes = {0.0, 1.0};
  std::vector<double> z_planes = {0.0, 1.0};
  std::shared_ptr< Utility::StructuredHexMesh > mesh = std::make_shared< Utility::StructuredHexMesh >(x_planes, y_planes, z_planes); 

  std::shared_ptr< Utility::PQLAQuadrature > direction_discretization = std::make_shared< Utility::PQLAQuadrature >(2);

  distribution->setMeshIndexDimensionDistributionObject( mesh );
  distribution->setDirectionIndexDimensionDistributionObject( direction_discretization );

  // Real time distribution
  time_distribution = std::make_shared< Utility::DeltaDistribution >(3.0);
  independent_time_distribution = std::make_shared<MonteCarlo::IndependentPhaseSpaceDimensionDistribution<MonteCarlo::TIME_DIMENSION>>(time_distribution);
  // 2*32*3 = 192 energy importance distributions that need to be formed

  // Real mesh distribution;
  std::vector<double> mesh_index_boundaries = {0.0, 1.0, 2.0};
  std::vector<double> mesh_distribution_values = {1.0, 3.0};
  mesh_distribution = std::make_shared< Utility::HistogramDistribution >( mesh_index_boundaries, mesh_distribution_values);

  // Real direction distribution
  std::vector<double> direction_index_boundaries;
  std::vector<double> direction_distribution_values;
  for(size_t i = 0; i < number_of_direction_elements + 1; ++i)
  {
    direction_index_boundaries.push_back(i);
  }
  for(size_t i = 0; i < number_of_direction_elements; ++i)
  {
    direction_distribution_values.push_back(i+1);
  }
  direction_distribution = std::make_shared< Utility::HistogramDistribution >( direction_index_boundaries, direction_distribution_values);

  // Real energy distribution;
  std::vector<double> energy_boundaries = {10.0, 11.0, 20.0, 22.0};
  std::vector<double> energy_values = {5.0, 2.0, 6.0};
  energy_distribution = std::make_shared< Utility::HistogramDistribution >(energy_boundaries, energy_values);

  
  importance_distribution_boundaries[MonteCarlo::SPATIAL_INDEX_DIMENSION] = mesh_index_boundaries;
  importance_distribution_boundaries[MonteCarlo::DIRECTION_INDEX_DIMENSION] = direction_index_boundaries;
  importance_distribution_boundaries[MonteCarlo::ENERGY_DIMENSION] = energy_boundaries;

  phase_space_dimension_importance_order.push_back(MonteCarlo::SPATIAL_INDEX_DIMENSION);
  phase_space_dimension_importance_order.push_back(MonteCarlo::DIRECTION_INDEX_DIMENSION);
  phase_space_dimension_importance_order.push_back(MonteCarlo::ENERGY_DIMENSION);

  std::vector< std::shared_ptr<MonteCarlo::PhaseSpaceDimensionDistribution> > mesh_importance_distributions;
  std::vector< std::shared_ptr<MonteCarlo::PhaseSpaceDimensionDistribution> > direction_importance_distributions;
  std::vector< std::shared_ptr<MonteCarlo::PhaseSpaceDimensionDistribution> > energy_importance_distributions;

  std::vector<double> mesh_importance_distribution_values;
  for( size_t mesh_index = 0; mesh_index < 2; ++mesh_index)
  {
    std::vector<double> spatial_direction_elements;
    double spatial_element_integrated_direction = 0;
    for( size_t direction_index = 0;  direction_index < 32; ++direction_index)
    {
      std::vector<double> direction_energy_elements;
      double direction_element_integrated_energy = 0;
      for( size_t energy_index = 0; energy_index < 3; ++energy_index)
      {
        direction_energy_elements.push_back( energy_index + 3*direction_index + 32*3*mesh_index + 1);
        direction_element_integrated_energy += energy_values[energy_index]*(energy_boundaries[energy_index+1] - energy_boundaries[energy_index]);
      }
      std::shared_ptr< Utility::HistogramDistribution > energy_importance_distribution = std::make_shared< Utility::HistogramDistribution >(energy_boundaries, direction_energy_elements);
      std::shared_ptr< MonteCarlo::PhaseSpaceDimensionDistribution > local_energy_importance_distribution = 
        std::make_shared< MonteCarlo::ImportanceSampledIndependentPhaseSpaceDimensionDistribution< MonteCarlo::ENERGY_DIMENSION > >( energy_distribution, energy_importance_distribution );
      energy_importance_distributions.push_back(local_energy_importance_distribution);

      spatial_direction_elements.push_back( direction_element_integrated_energy );
      spatial_element_integrated_direction += direction_element_integrated_energy;
    }
    importance_dimension_map[MonteCarlo::ENERGY_DIMENSION] = energy_importance_distributions;

    std::shared_ptr< Utility::HistogramDistribution > direction_importance_distribution = std::make_shared< Utility::HistogramDistribution >( direction_index_boundaries, spatial_direction_elements);
    std::shared_ptr< MonteCarlo::PhaseSpaceDimensionDistribution > local_direction_importance_distribution = 
      std::make_shared< MonteCarlo::ImportanceSampledIndependentPhaseSpaceDimensionDistribution< MonteCarlo::DIRECTION_INDEX_DIMENSION > >( direction_distribution, direction_importance_distribution);
    direction_importance_distributions.push_back(local_direction_importance_distribution);
    
    mesh_importance_distribution_values.push_back(spatial_element_integrated_direction);
  }
  importance_dimension_map[MonteCarlo::DIRECTION_INDEX_DIMENSION] = direction_importance_distributions;

  std::shared_ptr< Utility::HistogramDistribution > mesh_importance_distribution = std::make_shared< Utility::HistogramDistribution >( mesh_index_boundaries, mesh_importance_distribution_values);
  std::shared_ptr< MonteCarlo::PhaseSpaceDimensionDistribution > local_mesh_importance_distribution =
    std::make_shared< MonteCarlo::ImportanceSampledIndependentPhaseSpaceDimensionDistribution< MonteCarlo::SPATIAL_INDEX_DIMENSION > >( mesh_distribution, mesh_importance_distribution);
  mesh_importance_distributions.push_back(local_mesh_importance_distribution);

  importance_dimension_map[MonteCarlo::SPATIAL_INDEX_DIMENSION] = mesh_importance_distributions;

  // Initialize the random number generator
  Utility::RandomNumberGenerator::createStreams();
}

FRENSIE_CUSTOM_UNIT_TEST_SETUP_END();

//---------------------------------------------------------------------------//
// end tstGenericHistogramImportanceParticleDistribution.cpp
//---------------------------------------------------------------------------//
