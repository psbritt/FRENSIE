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

std::shared_ptr< Utility::StructuredHexMesh > mesh_object;
std::shared_ptr< Utility::PQLAQuadrature > direction_discretization_object;
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

  // Should be in mesh cell 1, direction cell 18, energy cell 3
  std::vector<double> fake_random_stream = {0.5, 0.5, 0.5, 0.5,  // random numbers for mesh sampling
                                            0.5, 0.5, 0.5,       // random numbers for direction sampling
                                            0.7};                // random number for energy sampling
  Utility::RandomNumberGenerator::setFakeStream(fake_random_stream);
  MonteCarlo::PhotonState photon( 0 );
  distribution->sample(photon);

  FRENSIE_CHECK_EQUAL(photon.getXPosition(), 1.5);
  FRENSIE_CHECK_EQUAL(photon.getYPosition(), 0.5);
  FRENSIE_CHECK_EQUAL(photon.getZPosition(), 0.5);

  FRENSIE_CHECK_FLOATING_EQUALITY(photon.getWeight(), 7.525974025974025983e-01, 1e-15);

  Utility::RandomNumberGenerator::unsetFakeStream();
}

FRENSIE_UNIT_TEST(GenericHistogramImportanceParticleDistribution, sampleAndRecordTrials)
{

  MonteCarlo::ParticleDistribution::DimensionCounterMap counter_map;

  distribution->initializeDimensionCounters(counter_map);

  std::vector<double> fake_random_stream = {0.5, 0.5, 0.5, 0.5,  // random numbers for mesh sampling
                                            0.5, 0.5, 0.5,       // random numbers for direction sampling
                                            0.7};                // random number for energy sampling
  Utility::RandomNumberGenerator::setFakeStream(fake_random_stream);
  MonteCarlo::PhotonState photon( 0 );
  distribution->sampleAndRecordTrials(photon,
                                      counter_map);

  FRENSIE_CHECK_EQUAL(photon.getXPosition(), 1.5);
  FRENSIE_CHECK_EQUAL(photon.getYPosition(), 0.5);
  FRENSIE_CHECK_EQUAL(photon.getZPosition(), 0.5);

  FRENSIE_CHECK_FLOATING_EQUALITY(photon.getWeight(),  7.525974025974025983e-01, 1e-15);
  FRENSIE_CHECK_FLOATING_EQUALITY(photon.getEnergy(), 18.411184210526315, 1e-15);

  Utility::RandomNumberGenerator::unsetFakeStream();
}

//---------------------------------------------------------------------------//
// Check that the distribution can be archived
FRENSIE_UNIT_TEST_TEMPLATE_EXPAND( StandardParticleDistribution,
                                   archive,
                                   TestArchives )
{
  FETCH_TEMPLATE_PARAM( 0, RawOArchive );
  FETCH_TEMPLATE_PARAM( 1, RawIArchive );

  typedef typename std::remove_pointer<RawOArchive>::type OArchive;
  typedef typename std::remove_pointer<RawIArchive>::type IArchive;

  std::string archive_base_name( "test_generic_histogram_importance_particle_distribution" );
  std::ostringstream archive_ostream;

  {
    std::unique_ptr<OArchive> oarchive;

    createOArchive( archive_base_name, archive_ostream, oarchive );

    std::shared_ptr<const MonteCarlo::ParticleDistribution> copy_distribution = distribution;

    FRENSIE_REQUIRE_NO_THROW( (*oarchive) << BOOST_SERIALIZATION_NVP(copy_distribution) );
  }

  // Copy the archive ostream to an istream
  std::istringstream archive_istream( archive_ostream.str() );

  // Load the archived distributions
  std::unique_ptr<IArchive> iarchive;

  createIArchive( archive_istream, iarchive );

  std::shared_ptr<const MonteCarlo::ParticleDistribution> new_distribution;

  FRENSIE_REQUIRE_NO_THROW( (*iarchive) >> BOOST_SERIALIZATION_NVP(new_distribution) );

  iarchive.reset();

  // Set the random number generator stream
  std::vector<double> fake_random_stream = {0.5, 0.5, 0.5, 0.5,  // random numbers for mesh sampling
                                            0.5, 0.5, 0.5,       // random numbers for direction sampling
                                            0.7};                // random number for energy sampling

  Utility::RandomNumberGenerator::setFakeStream( fake_random_stream );

  MonteCarlo::PhotonState photon( 0 );

  new_distribution->sample( photon );

  FRENSIE_CHECK_EQUAL(photon.getXPosition(), 1.5);
  FRENSIE_CHECK_EQUAL(photon.getYPosition(), 0.5);
  FRENSIE_CHECK_EQUAL(photon.getZPosition(), 0.5);

  FRENSIE_CHECK_FLOATING_EQUALITY(photon.getWeight(), 7.525974025974025983e-01, 1e-15);
  Utility::RandomNumberGenerator::unsetFakeStream();
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
  mesh_object = std::make_shared< Utility::StructuredHexMesh >(x_planes, y_planes, z_planes); 

  direction_discretization_object = std::make_shared< Utility::PQLAQuadrature >(2);

  distribution->setMeshIndexDimensionDistributionObject( mesh_object );
  distribution->setDirectionIndexDimensionDistributionObject( direction_discretization_object );

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
        direction_element_integrated_energy += direction_energy_elements[energy_index]*(energy_boundaries[energy_index+1] - energy_boundaries[energy_index]);
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
