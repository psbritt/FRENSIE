//---------------------------------------------------------------------------//
//!
//! \file   tstIndependentPhaseSpaceDimensionDistribution_indexDimensions.cpp
//! \author Philip Britt
//! \brief  Independent phase space dimension distribution unit tests
//!
//---------------------------------------------------------------------------//

// Std Lib Includes
#include <iostream>
#include <memory>

// FRENSIE Includes
#include "MonteCarlo_IndependentPhaseSpaceDimensionDistribution.hpp"
#include "MonteCarlo_PhaseSpaceDimensionTraits.hpp"
#include "Utility_BasicCartesianCoordinateConversionPolicy.hpp"
#include "Utility_UniformDistribution.hpp"
#include "Utility_DeltaDistribution.hpp"
#include "Utility_DiscreteDistribution.hpp"
#include "Utility_HistogramDistribution.hpp"
#include "Utility_ExponentialDistribution.hpp"
#include "Utility_RandomNumberGenerator.hpp"
#include "Utility_UnitTestHarnessWithMain.hpp"
#include "ArchiveTestHelpers.hpp"

//---------------------------------------------------------------------------//
// Testing Types
//---------------------------------------------------------------------------//

using namespace MonteCarlo;

typedef std::tuple<std::integral_constant<PhaseSpaceDimension,SPATIAL_INDEX_DIMENSION>,
                   std::integral_constant<PhaseSpaceDimension,DIRECTION_INDEX_DIMENSION>
                  > TestPhaseSpaceDimensions;

typedef TestArchiveHelper::TestArchives TestArchives;

//---------------------------------------------------------------------------//
// Testing Variables.
//---------------------------------------------------------------------------//
std::shared_ptr<const Utility::SpatialCoordinateConversionPolicy>
spatial_coord_conversion_policy( new Utility::BasicCartesianCoordinateConversionPolicy );

std::shared_ptr<const Utility::DirectionalCoordinateConversionPolicy>
directional_coord_conversion_policy( new Utility::BasicCartesianCoordinateConversionPolicy );

//---------------------------------------------------------------------------//
// Tests.
//---------------------------------------------------------------------------//
// Test that the dimension can be returned
FRENSIE_UNIT_TEST_TEMPLATE( IndependentPhaseSpaceDimensionDistribution_indexDimension,
                            getDimension,
                            TestPhaseSpaceDimensions )
{
  FETCH_TEMPLATE_PARAM( 0, WrappedDimension );
  constexpr PhaseSpaceDimension Dimension = WrappedDimension::value;
  
  std::shared_ptr<const Utility::UnivariateDistribution> basic_distribution(
                           new Utility::UniformDistribution( 0.5, 1.5, 0.5 ) );

  std::shared_ptr<const MonteCarlo::PhaseSpaceDimensionDistribution>
    dimension_distribution( new MonteCarlo::IndependentPhaseSpaceDimensionDistribution<Dimension>( basic_distribution ) );

  FRENSIE_CHECK_EQUAL( dimension_distribution->getDimension(), Dimension );
}
                         
//---------------------------------------------------------------------------//
// Test that the dimension class can be returned
FRENSIE_UNIT_TEST_TEMPLATE( IndependentPhaseSpaceDimensionDistribution_indexDimension,
                            getDimensionClass,
                            TestPhaseSpaceDimensions )
{
  FETCH_TEMPLATE_PARAM( 0, WrappedDimension );
  constexpr PhaseSpaceDimension Dimension = WrappedDimension::value;
  
  std::shared_ptr<const Utility::UnivariateDistribution> basic_distribution(
                           new Utility::UniformDistribution( 0.5, 1.5, 0.5 ) );

  std::shared_ptr<const MonteCarlo::PhaseSpaceDimensionDistribution>
    dimension_distribution( new MonteCarlo::IndependentPhaseSpaceDimensionDistribution<Dimension>( basic_distribution ) );

  FRENSIE_CHECK_EQUAL( dimension_distribution->getDimensionClass(),
                       MonteCarlo::PhaseSpaceDimensionTraits<Dimension>::getClass() );
}

//---------------------------------------------------------------------------//
// Test that the distribution type name can be returned
FRENSIE_UNIT_TEST_TEMPLATE( IndependentPhaseSpaceDimensionDistribution_indexDimension,
                            getDistributionTypeName,
                            TestPhaseSpaceDimensions )
{
  FETCH_TEMPLATE_PARAM( 0, WrappedDimension );
  constexpr PhaseSpaceDimension Dimension = WrappedDimension::value;
  
  std::shared_ptr<const Utility::UnivariateDistribution> basic_distribution(
                           new Utility::UniformDistribution( 0.5, 1.5, 0.5 ) );

  std::shared_ptr<const MonteCarlo::PhaseSpaceDimensionDistribution>
    dimension_distribution( new MonteCarlo::IndependentPhaseSpaceDimensionDistribution<Dimension>( basic_distribution ) );

  FRENSIE_CHECK_EQUAL( dimension_distribution->getDistributionTypeName(),
                       "Uniform Distribution" );
}

//---------------------------------------------------------------------------//
// Test if the distribution is independent
FRENSIE_UNIT_TEST_TEMPLATE( IndependentPhaseSpaceDimensionDistribution_indexDimension,
                            isIndependent,
                            TestPhaseSpaceDimensions )
{
  FETCH_TEMPLATE_PARAM( 0, WrappedDimension );
  constexpr PhaseSpaceDimension Dimension = WrappedDimension::value;
  
  std::shared_ptr<const Utility::UnivariateDistribution> basic_distribution(
                           new Utility::UniformDistribution( 0.5, 1.5, 0.5 ) );

  std::shared_ptr<const MonteCarlo::PhaseSpaceDimensionDistribution>
    dimension_distribution( new MonteCarlo::IndependentPhaseSpaceDimensionDistribution<Dimension>( basic_distribution ) );

  FRENSIE_CHECK( dimension_distribution->isIndependent() );
}

//---------------------------------------------------------------------------//
// Test if the distribution is dependent on another dimension
FRENSIE_UNIT_TEST_TEMPLATE( IndependentPhaseSpaceDimensionDistribution_indexDimension,
                            isDependentOnDimension,
                            TestPhaseSpaceDimensions )
{
  FETCH_TEMPLATE_PARAM( 0, WrappedDimension );
  constexpr PhaseSpaceDimension Dimension = WrappedDimension::value;
  
  std::shared_ptr<const Utility::UnivariateDistribution> basic_distribution(
                           new Utility::UniformDistribution( 0.5, 1.5, 0.5 ) );

  std::shared_ptr<const MonteCarlo::PhaseSpaceDimensionDistribution>
    dimension_distribution( new MonteCarlo::IndependentPhaseSpaceDimensionDistribution<Dimension>( basic_distribution ) );

  FRENSIE_CHECK( !dimension_distribution->isDependentOnDimension( MonteCarlo::SPATIAL_INDEX_DIMENSION ) );
  FRENSIE_CHECK( !dimension_distribution->isDependentOnDimension( MonteCarlo::DIRECTION_INDEX_DIMENSION ) );
  FRENSIE_CHECK( !dimension_distribution->isDependentOnDimension( MonteCarlo::TERTIARY_SPATIAL_DIMENSION ) );
  FRENSIE_CHECK( !dimension_distribution->isDependentOnDimension( MonteCarlo::PRIMARY_DIRECTIONAL_DIMENSION ) );
  FRENSIE_CHECK( !dimension_distribution->isDependentOnDimension( MonteCarlo::SECONDARY_DIRECTIONAL_DIMENSION ) );
  FRENSIE_CHECK( !dimension_distribution->isDependentOnDimension( MonteCarlo::TERTIARY_DIRECTIONAL_DIMENSION ) );
  FRENSIE_CHECK( !dimension_distribution->isDependentOnDimension( MonteCarlo::ENERGY_DIMENSION ) );
  FRENSIE_CHECK( !dimension_distribution->isDependentOnDimension( MonteCarlo::TIME_DIMENSION ) );
  FRENSIE_CHECK( !dimension_distribution->isDependentOnDimension( MonteCarlo::WEIGHT_DIMENSION ) );
  FRENSIE_CHECK( !dimension_distribution->isDependentOnDimension( MonteCarlo::SPATIAL_INDEX_DIMENSION ) );
  FRENSIE_CHECK( !dimension_distribution->isDependentOnDimension( MonteCarlo::DIRECTION_INDEX_DIMENSION ) );
}

//---------------------------------------------------------------------------//
// Test if the distribution is continuous
FRENSIE_UNIT_TEST_TEMPLATE( IndependentPhaseSpaceDimensionDistribution_indexDimension,
                            isContinuous,
                            TestPhaseSpaceDimensions )
{
  FETCH_TEMPLATE_PARAM( 0, WrappedDimension );
  constexpr PhaseSpaceDimension Dimension = WrappedDimension::value;
  
  std::shared_ptr<const Utility::UnivariateDistribution> basic_distribution(
                           new Utility::UniformDistribution( 0.5, 1.5, 0.5 ) );

  std::shared_ptr<const MonteCarlo::PhaseSpaceDimensionDistribution>
    dimension_distribution( new MonteCarlo::IndependentPhaseSpaceDimensionDistribution<Dimension>( basic_distribution ) );

  FRENSIE_CHECK( dimension_distribution->isContinuous() );

  basic_distribution.reset( new Utility::DeltaDistribution( 1.0 ) );

  dimension_distribution.reset( new MonteCarlo::IndependentPhaseSpaceDimensionDistribution<Dimension>( basic_distribution ) );

  FRENSIE_CHECK( !dimension_distribution->isContinuous() );
}

//---------------------------------------------------------------------------//
// Test if the distribution is tabular
FRENSIE_UNIT_TEST_TEMPLATE( IndependentPhaseSpaceDimensionDistribution_indexDimension,
                            isTabular,
                            TestPhaseSpaceDimensions )
{
  FETCH_TEMPLATE_PARAM( 0, WrappedDimension );
  constexpr PhaseSpaceDimension Dimension = WrappedDimension::value;
  
  std::shared_ptr<const Utility::UnivariateDistribution> basic_distribution(
                           new Utility::UniformDistribution( 0.5, 1.5, 0.5 ) );

  std::shared_ptr<const MonteCarlo::PhaseSpaceDimensionDistribution>
    dimension_distribution( new MonteCarlo::IndependentPhaseSpaceDimensionDistribution<Dimension>( basic_distribution ) );

  FRENSIE_CHECK( dimension_distribution->isTabular() );

  basic_distribution.reset( new Utility::ExponentialDistribution( 1.0, 1.0 ) );

  dimension_distribution.reset( new MonteCarlo::IndependentPhaseSpaceDimensionDistribution<Dimension>( basic_distribution ) );

  FRENSIE_CHECK( !dimension_distribution->isTabular() );
}

//---------------------------------------------------------------------------//
// Test if the distribution is uniform
FRENSIE_UNIT_TEST_TEMPLATE( IndependentPhaseSpaceDimensionDistribution_indexDimension,
                            isUniform,
                            TestPhaseSpaceDimensions )
{
  FETCH_TEMPLATE_PARAM( 0, WrappedDimension );
  constexpr PhaseSpaceDimension Dimension = WrappedDimension::value;
  
  std::shared_ptr<const Utility::UnivariateDistribution> basic_distribution(
                           new Utility::UniformDistribution( 0.5, 1.5, 0.5 ) );

  std::shared_ptr<const MonteCarlo::PhaseSpaceDimensionDistribution>
    dimension_distribution( new MonteCarlo::IndependentPhaseSpaceDimensionDistribution<Dimension>( basic_distribution ) );

  FRENSIE_CHECK( dimension_distribution->isUniform() );

  basic_distribution.reset( new Utility::ExponentialDistribution( 1.0, 1.0 ) );

  dimension_distribution.reset( new MonteCarlo::IndependentPhaseSpaceDimensionDistribution<Dimension>( basic_distribution ) );

  FRENSIE_CHECK( !dimension_distribution->isUniform() );
}

//---------------------------------------------------------------------------//
// Test if the distribution has the specified form
FRENSIE_UNIT_TEST_TEMPLATE( IndependentPhaseSpaceDimensionDistribution_indexDimension,
                            hasForm,
                            TestPhaseSpaceDimensions )
{
  FETCH_TEMPLATE_PARAM( 0, WrappedDimension );
  constexpr PhaseSpaceDimension Dimension = WrappedDimension::value;
  
  std::shared_ptr<const Utility::UnivariateDistribution> basic_distribution(
                           new Utility::UniformDistribution( 0.5, 1.5, 0.5 ) );

  std::shared_ptr<const MonteCarlo::PhaseSpaceDimensionDistribution>
    dimension_distribution( new MonteCarlo::IndependentPhaseSpaceDimensionDistribution<Dimension>( basic_distribution ) );

  FRENSIE_CHECK( dimension_distribution->hasForm( Utility::UNIFORM_DISTRIBUTION ) );
  FRENSIE_CHECK( !dimension_distribution->hasForm( Utility::EXPONENTIAL_DISTRIBUTION ) );

  basic_distribution.reset( new Utility::ExponentialDistribution( 1.0, 1.0 ) );

  dimension_distribution.reset( new MonteCarlo::IndependentPhaseSpaceDimensionDistribution<Dimension>( basic_distribution ) );

  FRENSIE_CHECK( !dimension_distribution->hasForm( Utility::UNIFORM_DISTRIBUTION ) );
  FRENSIE_CHECK( dimension_distribution->hasForm( Utility::EXPONENTIAL_DISTRIBUTION ) );
}

//---------------------------------------------------------------------------//
// Test if the distribution can be evaluated without a cascade
FRENSIE_UNIT_TEST_TEMPLATE( IndependentPhaseSpaceDimensionDistribution_indexDimension,
                            evaluateWithoutCascade,
                            TestPhaseSpaceDimensions )
{
  FETCH_TEMPLATE_PARAM( 0, WrappedDimension );
  constexpr PhaseSpaceDimension Dimension = WrappedDimension::value;
  
  std::vector<double> bin_boundaries {0.0, 1.0, 2.0};
  std::vector<double> bin_values {2.0, 3.0};

  std::shared_ptr<const Utility::UnivariateDistribution> basic_distribution(
                           new Utility::HistogramDistribution( bin_boundaries, bin_values ) );

  std::shared_ptr<const MonteCarlo::PhaseSpaceDimensionDistribution>
    dimension_distribution( new MonteCarlo::IndependentPhaseSpaceDimensionDistribution<Dimension>( basic_distribution ) );

  MonteCarlo::PhaseSpacePoint point( spatial_coord_conversion_policy,
                                     directional_coord_conversion_policy );
  setCoordinate<Dimension>( point, 0.01 );
  
  FRENSIE_CHECK_EQUAL( dimension_distribution->evaluateWithoutCascade( point ),
                       2.0 );

  setCoordinate<Dimension>( point, 0.5 );

  FRENSIE_CHECK_EQUAL( dimension_distribution->evaluateWithoutCascade( point ),
                       2.0 );

  setCoordinate<Dimension>( point, 1.0 );

  FRENSIE_CHECK_EQUAL( dimension_distribution->evaluateWithoutCascade( point ),
                       3.0 );


  FRENSIE_CHECK_THROW(setCoordinate<Dimension>( point, 72 ), std::runtime_error);
  
}

//---------------------------------------------------------------------------//
// Test if the distribution can be sampled without a cascade
FRENSIE_UNIT_TEST_TEMPLATE( IndependentPhaseSpaceDimensionDistribution_indexDimension,
                            sampleWithoutCascade,
                            TestPhaseSpaceDimensions )
{
  FETCH_TEMPLATE_PARAM( 0, WrappedDimension );
  constexpr PhaseSpaceDimension Dimension = WrappedDimension::value;
  
  std::vector<double> bin_boundaries {0.0, 1.0, 2.0};
  std::vector<double> bin_values {2.0, 3.0};

  std::shared_ptr<const Utility::UnivariateDistribution> basic_distribution(
                           new Utility::HistogramDistribution( bin_boundaries, bin_values ) );

  std::shared_ptr<const MonteCarlo::PhaseSpaceDimensionDistribution>
    dimension_distribution( new MonteCarlo::IndependentPhaseSpaceDimensionDistribution<Dimension>( basic_distribution ) );

  MonteCarlo::PhaseSpacePoint point( spatial_coord_conversion_policy,
                                     directional_coord_conversion_policy );

  // First random number is used to sample the index. The ones after that are for determining where in that element

  std::vector<size_t> expected_samples = {0, 1, 1};
  std::vector<double> index_fake_stream = {0.0, 0.5, 1-1e-15};
  std::vector<double> other_fake_stream = {0.5, 0.5, 0.5};
  for(size_t sample = 0; sample < expected_samples.size(); ++sample)
  {
    std::vector<double> fake_stream(4);
    fake_stream[0] = index_fake_stream[sample];
    fake_stream.insert( fake_stream.end(), other_fake_stream.begin(), other_fake_stream.end() );
    Utility::RandomNumberGenerator::setFakeStream( fake_stream );
    dimension_distribution->sampleWithoutCascade( point );
    FRENSIE_CHECK_EQUAL( getCoordinate<Dimension>( point ), expected_samples[sample] );
    FRENSIE_CHECK_EQUAL( getCoordinateWeight<Dimension>( point ), 1.0 );
    Utility::RandomNumberGenerator::unsetFakeStream();
  }
}

//---------------------------------------------------------------------------//
// Test if the distribution can be sampled without a cascade and the
// trials can be counted
FRENSIE_UNIT_TEST_TEMPLATE( IndependentPhaseSpaceDimensionDistribution_indexDimension,
                            sampleAndRecordTrialsWithoutCascade,
                            TestPhaseSpaceDimensions )
{
  FETCH_TEMPLATE_PARAM( 0, WrappedDimension );
  constexpr PhaseSpaceDimension Dimension = WrappedDimension::value;

  std::vector<double> bin_boundaries {0.0, 1.0, 2.0};
  std::vector<double> bin_values {2.0, 3.0};

  std::shared_ptr<const Utility::UnivariateDistribution> basic_distribution(
                           new Utility::HistogramDistribution( bin_boundaries, bin_values ) );

  std::shared_ptr<const MonteCarlo::PhaseSpaceDimensionDistribution>
    dimension_distribution( new MonteCarlo::IndependentPhaseSpaceDimensionDistribution<Dimension>( basic_distribution ) );

  MonteCarlo::PhaseSpacePoint point( spatial_coord_conversion_policy,
                                     directional_coord_conversion_policy );

  // First random number is used to sample the index. The ones after that are for determining where in that element

  typename MonteCarlo::IndependentPhaseSpaceDimensionDistribution<Dimension>::Counter trials = 0;  
  std::vector<size_t> expected_samples = {0, 1, 1};
  std::vector<double> index_fake_stream = {0.0, 0.5, 1-1e-15};
  std::vector<double> other_fake_stream = {0.5, 0.5, 0.5};
  for(size_t sample = 0; sample < expected_samples.size(); ++sample)
  {
    std::vector<double> fake_stream(4);
    fake_stream[0] = index_fake_stream[sample];
    fake_stream.insert( fake_stream.end(), other_fake_stream.begin(), other_fake_stream.end() );
    Utility::RandomNumberGenerator::setFakeStream( fake_stream );
    dimension_distribution->sampleAndRecordTrialsWithoutCascade( point, trials );
    FRENSIE_CHECK_EQUAL( getCoordinate<Dimension>( point ), expected_samples[sample] );
    FRENSIE_CHECK_EQUAL( getCoordinateWeight<Dimension>( point ), 1.0 );
    // Just counts the number of times a distribution has been sampled
    FRENSIE_CHECK_EQUAL( trials, sample+1 );
    Utility::RandomNumberGenerator::unsetFakeStream();
  }
}

//---------------------------------------------------------------------------//
// Test that the dimension value can be set and weighted appropriately
FRENSIE_UNIT_TEST_TEMPLATE( IndependentPhaseSpaceDimensionDistribution_indexDimension,
                            setDimensionValueAndApplyWeight,
                            TestPhaseSpaceDimensions )
{
  FETCH_TEMPLATE_PARAM( 0, WrappedDimension );
  constexpr PhaseSpaceDimension Dimension = WrappedDimension::value;
  
  std::vector<double> bin_boundaries {0.0, 1.0, 2.0};
  std::vector<double> bin_values {2.0, 3.0};

  std::shared_ptr<const Utility::UnivariateDistribution> basic_distribution(
                           new Utility::HistogramDistribution( bin_boundaries, bin_values ) );

  std::shared_ptr<const MonteCarlo::PhaseSpaceDimensionDistribution>
    dimension_distribution( new MonteCarlo::IndependentPhaseSpaceDimensionDistribution<Dimension>( basic_distribution ) );

  MonteCarlo::PhaseSpacePoint point( spatial_coord_conversion_policy,
                                     directional_coord_conversion_policy );

  // weight values are just values of PDF from above distribution

  dimension_distribution->setDimensionValueAndApplyWeight( point, 0 );

  FRENSIE_CHECK_EQUAL( getCoordinate<Dimension>( point ), 0 );
  FRENSIE_CHECK_FLOATING_EQUALITY( getCoordinateWeight<Dimension>( point ), 2.0/5.0, 1e-14 );

  dimension_distribution->setDimensionValueAndApplyWeight( point, 1 );

  FRENSIE_CHECK_EQUAL( getCoordinate<Dimension>( point ), 1 );
  FRENSIE_CHECK_FLOATING_EQUALITY( getCoordinateWeight<Dimension>( point ), 3.0/5.0, 1e-14 );
}

//---------------------------------------------------------------------------//
// Check that the distribution can be archived
FRENSIE_UNIT_TEST_TEMPLATE_EXPAND( PhaseSpaceDimension,
                                   archive,
                                   TestArchives )
{
  FETCH_TEMPLATE_PARAM( 0, RawOArchive );
  FETCH_TEMPLATE_PARAM( 1, RawIArchive );

  typedef typename std::remove_pointer<RawOArchive>::type OArchive;
  typedef typename std::remove_pointer<RawIArchive>::type IArchive;

  std::string archive_base_name( "test_independent_phase_dimension_distribution" );
  std::ostringstream archive_ostream;

  {
    std::unique_ptr<OArchive> oarchive;

    createOArchive( archive_base_name, archive_ostream, oarchive );

    std::vector<double> bin_boundaries {0.0, 1.0, 2.0};
    std::vector<double> bin_values {2.0, 3.0};

    std::shared_ptr<const Utility::UnivariateDistribution> basic_distribution(
                            new Utility::HistogramDistribution( bin_boundaries, bin_values ) );

    std::shared_ptr<const MonteCarlo::PhaseSpaceDimensionDistribution>
      spatial_index_dimension_distribution( new MonteCarlo::IndependentSpatialIndexDimensionDistribution( basic_distribution ) );

    std::shared_ptr<const MonteCarlo::PhaseSpaceDimensionDistribution>
      direction_index_dimension_distribution( new MonteCarlo::IndependentDirectionIndexDimensionDistribution( basic_distribution ) );

    FRENSIE_REQUIRE_NO_THROW( (*oarchive) << BOOST_SERIALIZATION_NVP(spatial_index_dimension_distribution) );
    FRENSIE_REQUIRE_NO_THROW( (*oarchive) << BOOST_SERIALIZATION_NVP(direction_index_dimension_distribution) );
  }

  // Copy the archive ostream to an istream
  std::istringstream archive_istream( archive_ostream.str() );

  // Load the archived distributions
  std::unique_ptr<IArchive> iarchive;

  createIArchive( archive_istream, iarchive );

  std::shared_ptr<const MonteCarlo::PhaseSpaceDimensionDistribution>
    spatial_index_dimension_distribution,
    direction_index_dimension_distribution;

  FRENSIE_REQUIRE_NO_THROW( (*iarchive) >> BOOST_SERIALIZATION_NVP(spatial_index_dimension_distribution) );
  FRENSIE_REQUIRE_NO_THROW( (*iarchive) >> BOOST_SERIALIZATION_NVP(direction_index_dimension_distribution) );

  iarchive.reset();

  {
    FRENSIE_CHECK_EQUAL( spatial_index_dimension_distribution->getDimension(),
                         SPATIAL_INDEX_DIMENSION );

    MonteCarlo::PhaseSpacePoint point( spatial_coord_conversion_policy,
                                     directional_coord_conversion_policy );

    setCoordinate<SPATIAL_INDEX_DIMENSION>( point, 0 );
  
    FRENSIE_CHECK_EQUAL( spatial_index_dimension_distribution->evaluateWithoutCascade( point ),
                         2.0 );

    setCoordinate<SPATIAL_INDEX_DIMENSION>( point, 1 );
    
    FRENSIE_CHECK_EQUAL( spatial_index_dimension_distribution->evaluateWithoutCascade( point ),
                         3.0 );
  }

  {
    FRENSIE_CHECK_EQUAL( direction_index_dimension_distribution->getDimension(),
                         DIRECTION_INDEX_DIMENSION );

    MonteCarlo::PhaseSpacePoint point( spatial_coord_conversion_policy,
                                     directional_coord_conversion_policy );
    setCoordinate<DIRECTION_INDEX_DIMENSION>( point, 0 );
  
    FRENSIE_CHECK_EQUAL( direction_index_dimension_distribution->evaluateWithoutCascade( point ),
                         2.0 );

    setCoordinate<DIRECTION_INDEX_DIMENSION>( point, 1 );
    
    FRENSIE_CHECK_EQUAL( direction_index_dimension_distribution->evaluateWithoutCascade( point ),
                         3.0 );
  }
}

//---------------------------------------------------------------------------//
// Custom setup
//---------------------------------------------------------------------------//
FRENSIE_CUSTOM_UNIT_TEST_SETUP_BEGIN();

FRENSIE_CUSTOM_UNIT_TEST_INIT()
{
  // Initialize dimension discretizations
  std::vector<double> x_planes( {0.0, 0.5, 1.0} ),
                      y_planes( {0.0, 0.5} ),
                      z_planes( {0.0, 0.5} );

  std::shared_ptr<Utility::StructuredHexMesh> source_mesh = std::make_shared<Utility::StructuredHexMesh>(x_planes, y_planes, z_planes);

  MonteCarlo::PhaseSpaceDimensionTraits<MonteCarlo::SPATIAL_INDEX_DIMENSION>::setMesh(source_mesh);

  std::shared_ptr<Utility::PQLAQuadrature> source_direction_discretization = std::make_shared<Utility::PQLAQuadrature>(2);

  MonteCarlo::PhaseSpaceDimensionTraits<MonteCarlo::DIRECTION_INDEX_DIMENSION>::setDirectionDiscretization(source_direction_discretization);

  // Initialize the random number generator
  Utility::RandomNumberGenerator::createStreams();
}

FRENSIE_CUSTOM_UNIT_TEST_SETUP_END();

//---------------------------------------------------------------------------//
// end tstIndependentPhaseSpaceDimensionDistribution_indexDimension.cpp
//---------------------------------------------------------------------------//
