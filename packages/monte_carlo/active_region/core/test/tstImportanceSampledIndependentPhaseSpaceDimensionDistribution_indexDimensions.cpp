//---------------------------------------------------------------------------//
//!
//! \file   tstImportanceSampledIndependentPhaseSpaceDimensionDistribution_indexDimensions.cpp
//! \author Philip Britt
//! \brief  Importance sampled indep. phase space dimension dist. unit tests
//!
//---------------------------------------------------------------------------//

// Std Lib Includes
#include <iostream>
#include <memory>

// FRENSIE Includes
#include "MonteCarlo_ImportanceSampledIndependentPhaseSpaceDimensionDistribution.hpp"
#include "MonteCarlo_PhaseSpaceDimensionTraits.hpp"
#include "Utility_BasicCartesianCoordinateConversionPolicy.hpp"
#include "Utility_HistogramDistribution.hpp"
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
// Testing Variables
//---------------------------------------------------------------------------//
std::shared_ptr<const Utility::SpatialCoordinateConversionPolicy>
spatial_coord_conversion_policy( new Utility::BasicCartesianCoordinateConversionPolicy );

std::shared_ptr<const Utility::DirectionalCoordinateConversionPolicy>
directional_coord_conversion_policy( new Utility::BasicCartesianCoordinateConversionPolicy );

//---------------------------------------------------------------------------//
// Tests.
//---------------------------------------------------------------------------//
// Test that the dimension can be returned
FRENSIE_UNIT_TEST_TEMPLATE(
                   ImportanceSampledIndependentPhaseSpaceDimensionDistribution_indexDimensions,
                   getDimension,
                   TestPhaseSpaceDimensions )
{
  FETCH_TEMPLATE_PARAM( 0, WrappedDimension );
  constexpr PhaseSpaceDimension Dimension = WrappedDimension::value;
  
  std::vector<double> bin_boundaries {0.0, 1.0, 2.0};
  std::vector<double> bin_values {2.0, 3.0};
  std::shared_ptr<const Utility::HistogramDistribution> true_distribution = std::make_shared< Utility::HistogramDistribution >( bin_boundaries, bin_values);
  // Only bin_values needs to change
  bin_values[0] = 1.0;
  bin_values[1] = 5.0;
  std::shared_ptr<const Utility::HistogramDistribution> importance_distribution = std::make_shared< Utility::HistogramDistribution >( bin_boundaries, bin_values);

  std::shared_ptr<const MonteCarlo::PhaseSpaceDimensionDistribution>
    dimension_distribution( new MonteCarlo::ImportanceSampledIndependentPhaseSpaceDimensionDistribution<Dimension>( true_distribution,
                                    importance_distribution ) );

  FRENSIE_CHECK_EQUAL( dimension_distribution->getDimension(), Dimension );
}

//---------------------------------------------------------------------------//
// Test that the dimension can be returned
FRENSIE_UNIT_TEST_TEMPLATE(
                   ImportanceSampledIndependentPhaseSpaceDimensionDistribution_indexDimensions,
                   getDimensionClass,
                   TestPhaseSpaceDimensions )
{
  FETCH_TEMPLATE_PARAM( 0, WrappedDimension );
  constexpr PhaseSpaceDimension Dimension = WrappedDimension::value;

  std::vector<double> bin_boundaries {0.0, 1.0, 2.0};
  std::vector<double> bin_values {2.0, 3.0};
  std::shared_ptr<const Utility::HistogramDistribution> true_distribution = std::make_shared< Utility::HistogramDistribution >( bin_boundaries, bin_values);
  // Only bin_values needs to change
  bin_values[0] = 1.0;
  bin_values[1] = 5.0;
  std::shared_ptr<const Utility::HistogramDistribution> importance_distribution = std::make_shared< Utility::HistogramDistribution >( bin_boundaries, bin_values);

  std::shared_ptr<const MonteCarlo::PhaseSpaceDimensionDistribution>
    dimension_distribution( new MonteCarlo::ImportanceSampledIndependentPhaseSpaceDimensionDistribution<Dimension>( true_distribution,
                                    importance_distribution ) );

  FRENSIE_CHECK_EQUAL( dimension_distribution->getDimensionClass(),
                       MonteCarlo::PhaseSpaceDimensionTraits<Dimension>::getClass() );
  
}

//---------------------------------------------------------------------------//
// Test that the distribution type name can be returned
FRENSIE_UNIT_TEST_TEMPLATE(
                   ImportanceSampledIndependentPhaseSpaceDimensionDistribution_indexDimensions,
                   getDistributionTypeName,
                   TestPhaseSpaceDimensions )
{
  FETCH_TEMPLATE_PARAM( 0, WrappedDimension );
  constexpr PhaseSpaceDimension Dimension = WrappedDimension::value;

  std::vector<double> bin_boundaries {0.0, 1.0, 2.0};
  std::vector<double> bin_values {2.0, 3.0};
  std::shared_ptr<const Utility::HistogramDistribution> true_distribution = std::make_shared< Utility::HistogramDistribution >( bin_boundaries, bin_values);
  // Only bin_values needs to change
  bin_values[0] = 1.0;
  bin_values[1] = 5.0;
  std::shared_ptr<const Utility::HistogramDistribution> importance_distribution = std::make_shared< Utility::HistogramDistribution >( bin_boundaries, bin_values);

  std::shared_ptr<const MonteCarlo::PhaseSpaceDimensionDistribution>
    dimension_distribution( new MonteCarlo::ImportanceSampledIndependentPhaseSpaceDimensionDistribution<Dimension>( true_distribution,
                                    importance_distribution ) );

  FRENSIE_CHECK_EQUAL( dimension_distribution->getDistributionTypeName(),
                       "Histogram Distribution" );
}

//---------------------------------------------------------------------------//
// Test if the distribution is independent
FRENSIE_UNIT_TEST_TEMPLATE(
                   ImportanceSampledIndependentPhaseSpaceDimensionDistribution_indexDimensions,
                   isIndependent,
                   TestPhaseSpaceDimensions )
{
  FETCH_TEMPLATE_PARAM( 0, WrappedDimension );
  constexpr PhaseSpaceDimension Dimension = WrappedDimension::value;

  std::vector<double> bin_boundaries {0.0, 1.0, 2.0};
  std::vector<double> bin_values {2.0, 3.0};
  std::shared_ptr<const Utility::HistogramDistribution> true_distribution = std::make_shared< Utility::HistogramDistribution >( bin_boundaries, bin_values);
  // Only bin_values needs to change
  bin_values[0] = 1.0;
  bin_values[1] = 5.0;
  std::shared_ptr<const Utility::HistogramDistribution> importance_distribution = std::make_shared< Utility::HistogramDistribution >( bin_boundaries, bin_values);

  std::shared_ptr<const MonteCarlo::PhaseSpaceDimensionDistribution>
    dimension_distribution( new MonteCarlo::ImportanceSampledIndependentPhaseSpaceDimensionDistribution<Dimension>( true_distribution,
                                    importance_distribution ) );

  FRENSIE_CHECK( dimension_distribution->isIndependent() );
  
}

//---------------------------------------------------------------------------//
// Test if the distribution is dependent on another dimension
FRENSIE_UNIT_TEST_TEMPLATE(
                   ImportanceSampledIndependentPhaseSpaceDimensionDistribution_indexDimensions,
                   isDependentOnDimension,
                   TestPhaseSpaceDimensions )
{
  FETCH_TEMPLATE_PARAM( 0, WrappedDimension );
  constexpr PhaseSpaceDimension Dimension = WrappedDimension::value;

  std::vector<double> bin_boundaries {0.0, 1.0, 2.0};
  std::vector<double> bin_values {2.0, 3.0};
  std::shared_ptr<const Utility::HistogramDistribution> true_distribution = std::make_shared< Utility::HistogramDistribution >( bin_boundaries, bin_values);
  // Only bin_values needs to change
  bin_values[0] = 1.0;
  bin_values[1] = 5.0;
  std::shared_ptr<const Utility::HistogramDistribution> importance_distribution = std::make_shared< Utility::HistogramDistribution >( bin_boundaries, bin_values);

  std::shared_ptr<const MonteCarlo::PhaseSpaceDimensionDistribution>
    dimension_distribution( new MonteCarlo::ImportanceSampledIndependentPhaseSpaceDimensionDistribution<Dimension>( true_distribution,
                                    importance_distribution ) );

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
FRENSIE_UNIT_TEST_TEMPLATE(
                   ImportanceSampledIndependentPhaseSpaceDimensionDistribution_indexDimensions,
                   isContinuous,
                   TestPhaseSpaceDimensions )
{
  FETCH_TEMPLATE_PARAM( 0, WrappedDimension );
  constexpr PhaseSpaceDimension Dimension = WrappedDimension::value;

  std::vector<double> bin_boundaries {0.0, 1.0, 2.0};
  std::vector<double> bin_values {2.0, 3.0};
  std::shared_ptr<const Utility::HistogramDistribution> true_distribution = std::make_shared< Utility::HistogramDistribution >( bin_boundaries, bin_values);
  // Only bin_values needs to change
  bin_values[0] = 1.0;
  bin_values[1] = 5.0;
  std::shared_ptr<const Utility::HistogramDistribution> importance_distribution = std::make_shared< Utility::HistogramDistribution >( bin_boundaries, bin_values);

  std::shared_ptr<const MonteCarlo::PhaseSpaceDimensionDistribution>
    dimension_distribution( new MonteCarlo::ImportanceSampledIndependentPhaseSpaceDimensionDistribution<Dimension>( true_distribution,
                                    importance_distribution ) );

  FRENSIE_CHECK( dimension_distribution->isContinuous() );
}

//---------------------------------------------------------------------------//
// Test if the distribution is tabular
FRENSIE_UNIT_TEST_TEMPLATE(
                   ImportanceSampledIndependentPhaseSpaceDimensionDistribution_indexDimensions,
                   isTabular,
                   TestPhaseSpaceDimensions )
{
  FETCH_TEMPLATE_PARAM( 0, WrappedDimension );
  constexpr PhaseSpaceDimension Dimension = WrappedDimension::value;

  std::vector<double> bin_boundaries {0.0, 1.0, 2.0};
  std::vector<double> bin_values {2.0, 3.0};
  std::shared_ptr<const Utility::HistogramDistribution> true_distribution = std::make_shared< Utility::HistogramDistribution >( bin_boundaries, bin_values);
  // Only bin_values needs to change
  bin_values[0] = 1.0;
  bin_values[1] = 5.0;
  std::shared_ptr<const Utility::HistogramDistribution> importance_distribution = std::make_shared< Utility::HistogramDistribution >( bin_boundaries, bin_values);

  std::shared_ptr<const MonteCarlo::PhaseSpaceDimensionDistribution>
    dimension_distribution( new MonteCarlo::ImportanceSampledIndependentPhaseSpaceDimensionDistribution<Dimension>( true_distribution,
                                    importance_distribution ) );

  FRENSIE_CHECK( dimension_distribution->isTabular() );
}

//---------------------------------------------------------------------------//
// Test if the distribution is uniform
FRENSIE_UNIT_TEST_TEMPLATE(
                   ImportanceSampledIndependentPhaseSpaceDimensionDistribution_indexDimensions,
                   isUniform,
                   TestPhaseSpaceDimensions )
{
  FETCH_TEMPLATE_PARAM( 0, WrappedDimension );
  constexpr PhaseSpaceDimension Dimension = WrappedDimension::value;

  std::vector<double> bin_boundaries {0.0, 1.0, 2.0};
  std::vector<double> bin_values {2.0, 3.0};
  std::shared_ptr<const Utility::HistogramDistribution> true_distribution = std::make_shared< Utility::HistogramDistribution >( bin_boundaries, bin_values);
  // Only bin_values needs to change
  bin_values[0] = 1.0;
  bin_values[1] = 5.0;
  std::shared_ptr<const Utility::HistogramDistribution> importance_distribution = std::make_shared< Utility::HistogramDistribution >( bin_boundaries, bin_values);

  std::shared_ptr<const MonteCarlo::PhaseSpaceDimensionDistribution>
    dimension_distribution( new MonteCarlo::ImportanceSampledIndependentPhaseSpaceDimensionDistribution<Dimension>( true_distribution,
                                    importance_distribution ) );

  FRENSIE_CHECK( !dimension_distribution->isUniform() );
}

//---------------------------------------------------------------------------//
// Test if the distribution has the specified form
FRENSIE_UNIT_TEST_TEMPLATE(
                   ImportanceSampledIndependentPhaseSpaceDimensionDistribution_indexDimensions,
                   hasForm,
                   TestPhaseSpaceDimensions )
{
  FETCH_TEMPLATE_PARAM( 0, WrappedDimension );
  constexpr PhaseSpaceDimension Dimension = WrappedDimension::value;

  std::vector<double> bin_boundaries {0.0, 1.0, 2.0};
  std::vector<double> bin_values {2.0, 3.0};
  std::shared_ptr<const Utility::HistogramDistribution> true_distribution = std::make_shared< Utility::HistogramDistribution >( bin_boundaries, bin_values);
  // Only bin_values needs to change
  bin_values[0] = 1.0;
  bin_values[1] = 5.0;
  std::shared_ptr<const Utility::HistogramDistribution> importance_distribution = std::make_shared< Utility::HistogramDistribution >( bin_boundaries, bin_values);

  std::shared_ptr<const MonteCarlo::PhaseSpaceDimensionDistribution>
    dimension_distribution( new MonteCarlo::ImportanceSampledIndependentPhaseSpaceDimensionDistribution<Dimension>( true_distribution,
                                    importance_distribution ) );

  FRENSIE_CHECK( dimension_distribution->hasForm( Utility::HISTOGRAM_DISTRIBUTION ) );
}

//---------------------------------------------------------------------------//
// Test if the distribution can be evaluated without a cascade
FRENSIE_UNIT_TEST_TEMPLATE(
                   ImportanceSampledIndependentPhaseSpaceDimensionDistribution_indexDimensions,
                   evaluateWithoutCascade,
                   TestPhaseSpaceDimensions )
{
  FETCH_TEMPLATE_PARAM( 0, WrappedDimension );
  constexpr PhaseSpaceDimension Dimension = WrappedDimension::value;

  std::vector<double> bin_boundaries {0.0, 1.0, 2.0};
  std::vector<double> bin_values {2.0, 3.0};
  std::shared_ptr<const Utility::HistogramDistribution> true_distribution = std::make_shared< Utility::HistogramDistribution >( bin_boundaries, bin_values);
  // Only bin_values needs to change
  bin_values[0] = 1.0;
  bin_values[1] = 5.0;
  std::shared_ptr<const Utility::HistogramDistribution> importance_distribution = std::make_shared< Utility::HistogramDistribution >( bin_boundaries, bin_values);

  std::shared_ptr<const MonteCarlo::PhaseSpaceDimensionDistribution>
    dimension_distribution( new MonteCarlo::ImportanceSampledIndependentPhaseSpaceDimensionDistribution<Dimension>( true_distribution,
                                    importance_distribution ) );

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
FRENSIE_UNIT_TEST_TEMPLATE(
                   ImportanceSampledIndependentPhaseSpaceDimensionDistribution_indexDimensions,
                   sampleWithoutCascade,
                   TestPhaseSpaceDimensions )
{
  FETCH_TEMPLATE_PARAM( 0, WrappedDimension );
  constexpr PhaseSpaceDimension Dimension = WrappedDimension::value;

  std::vector<double> bin_boundaries {0.0, 1.0, 2.0};
  std::vector<double> bin_values {2.0, 3.0};
  std::shared_ptr<const Utility::HistogramDistribution> true_distribution = std::make_shared< Utility::HistogramDistribution >( bin_boundaries, bin_values);
  // Only bin_values needs to change
  bin_values[0] = 1.0;
  bin_values[1] = 5.0;
  std::shared_ptr<const Utility::HistogramDistribution> importance_distribution = std::make_shared< Utility::HistogramDistribution >( bin_boundaries, bin_values);

  std::shared_ptr<const MonteCarlo::PhaseSpaceDimensionDistribution>
    dimension_distribution( new MonteCarlo::ImportanceSampledIndependentPhaseSpaceDimensionDistribution<Dimension>( true_distribution,
                                    importance_distribution ) );

  MonteCarlo::PhaseSpacePoint point( spatial_coord_conversion_policy,
                                     directional_coord_conversion_policy );

  // First random number is used to sample the index. The ones after that are for determining where in that element

  std::vector<size_t> expected_samples = {0, 1, 1};
  std::vector<double> expected_weights = {(2.0/5.0)/(1.0/6.0), (3.0/5.0)/(5.0/6.0), (3.0/5.0)/(5.0/6.0)};
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
    FRENSIE_CHECK_FLOATING_EQUALITY( getCoordinateWeight<Dimension>( point ), expected_weights[sample], 1e-14);
    Utility::RandomNumberGenerator::unsetFakeStream();
  }
}

//---------------------------------------------------------------------------//
// Test if the distribution can be sampled without a cascade
FRENSIE_UNIT_TEST_TEMPLATE(
                   ImportanceSampledIndependentPhaseSpaceDimensionDistribution_indexDimensions,
                   sampleAndRecordTrialsWithoutCascade,
                   TestPhaseSpaceDimensions )
{
  FETCH_TEMPLATE_PARAM( 0, WrappedDimension );
  constexpr PhaseSpaceDimension Dimension = WrappedDimension::value;

  std::vector<double> bin_boundaries {0.0, 1.0, 2.0};
  std::vector<double> bin_values {2.0, 3.0};
  std::shared_ptr<const Utility::HistogramDistribution> true_distribution = std::make_shared< Utility::HistogramDistribution >( bin_boundaries, bin_values);
  // Only bin_values needs to change
  bin_values[0] = 1.0;
  bin_values[1] = 5.0;
  std::shared_ptr<const Utility::HistogramDistribution> importance_distribution = std::make_shared< Utility::HistogramDistribution >( bin_boundaries, bin_values);

  std::shared_ptr<const MonteCarlo::PhaseSpaceDimensionDistribution>
    dimension_distribution( new MonteCarlo::ImportanceSampledIndependentPhaseSpaceDimensionDistribution<Dimension>( true_distribution,
                                    importance_distribution ) );


  MonteCarlo::PhaseSpacePoint point( spatial_coord_conversion_policy,
                                     directional_coord_conversion_policy );

  typename MonteCarlo::IndependentPhaseSpaceDimensionDistribution<Dimension>::Counter trials = 0;  
  std::vector<size_t> expected_samples = {0, 1, 1};
  std::vector<double> expected_weights = {(2.0/5.0)/(1.0/6.0), (3.0/5.0)/(5.0/6.0), (3.0/5.0)/(5.0/6.0)};
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
    FRENSIE_CHECK_FLOATING_EQUALITY( getCoordinateWeight<Dimension>( point ), expected_weights[sample], 1e-14);
    // Just counts the number of times a distribution has been sampled
    FRENSIE_CHECK_EQUAL( trials, sample+1 );
    Utility::RandomNumberGenerator::unsetFakeStream();
  }
}

//---------------------------------------------------------------------------//
// Test that the dimension value can be set and weighted appropriately
FRENSIE_UNIT_TEST_TEMPLATE(
                   ImportanceSampledIndependentPhaseSpaceDimensionDistribution_indexDimensions,
                   setDimensionValueAndApplyWeight,
                   TestPhaseSpaceDimensions )
{
  FETCH_TEMPLATE_PARAM( 0, WrappedDimension );
  constexpr PhaseSpaceDimension Dimension = WrappedDimension::value;

  std::vector<double> bin_boundaries {0.0, 1.0, 2.0};
  std::vector<double> bin_values {2.0, 3.0};
  std::shared_ptr<const Utility::HistogramDistribution> true_distribution = std::make_shared< Utility::HistogramDistribution >( bin_boundaries, bin_values);
  // Only bin_values needs to change
  bin_values[0] = 1.0;
  bin_values[1] = 5.0;
  std::shared_ptr<const Utility::HistogramDistribution> importance_distribution = std::make_shared< Utility::HistogramDistribution >( bin_boundaries, bin_values);

  std::shared_ptr<const MonteCarlo::PhaseSpaceDimensionDistribution>
    dimension_distribution( new MonteCarlo::ImportanceSampledIndependentPhaseSpaceDimensionDistribution<Dimension>( true_distribution,
                                    importance_distribution ) );

  MonteCarlo::PhaseSpacePoint point( spatial_coord_conversion_policy,
                                     directional_coord_conversion_policy );

  // weight values are just values of PDF from above distribution

  dimension_distribution->setDimensionValueAndApplyWeight( point, 0 );

  FRENSIE_CHECK_EQUAL( getCoordinate<Dimension>( point ), 0 );
  FRENSIE_CHECK_FLOATING_EQUALITY( getCoordinateWeight<Dimension>( point ), (2.0/5.0), 1e-14 );

  dimension_distribution->setDimensionValueAndApplyWeight( point, 1 );

  FRENSIE_CHECK_EQUAL( getCoordinate<Dimension>( point ), 1 );
  FRENSIE_CHECK_FLOATING_EQUALITY( getCoordinateWeight<Dimension>( point ), (3.0/5.0), 1e-14 );
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
    std::shared_ptr<const Utility::HistogramDistribution> true_distribution = std::make_shared< Utility::HistogramDistribution >( bin_boundaries, bin_values);
    // Only bin_values needs to change
    bin_values[0] = 1.0;
    bin_values[1] = 5.0;
    std::shared_ptr<const Utility::HistogramDistribution> importance_distribution = std::make_shared< Utility::HistogramDistribution >( bin_boundaries, bin_values);

    std::shared_ptr<const MonteCarlo::PhaseSpaceDimensionDistribution>
      spatial_index_dimension_distribution( new MonteCarlo::ImportanceSampledIndependentSpatialIndexDimensionDistribution( true_distribution, importance_distribution ) );

    std::shared_ptr<const MonteCarlo::PhaseSpaceDimensionDistribution>
      direction_index_dimension_distribution( new MonteCarlo::ImportanceSampledIndependentDirectionIndexDimensionDistribution( true_distribution, importance_distribution ) );

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

    // First random number is used to sample the index. The ones after that are for determining where in that element

    std::vector<size_t> expected_samples = {0, 1, 1};
    std::vector<double> expected_weights = {(2.0/5.0)/(1.0/6.0), (3.0/5.0)/(5.0/6.0), (3.0/5.0)/(5.0/6.0)};
    std::vector<double> index_fake_stream = {0.0, 0.5, 1-1e-15};
    std::vector<double> other_fake_stream = {0.5, 0.5, 0.5};
    for(size_t sample = 0; sample < expected_samples.size(); ++sample)
    {
      std::vector<double> fake_stream(4);
      fake_stream[0] = index_fake_stream[sample];
      fake_stream.insert( fake_stream.end(), other_fake_stream.begin(), other_fake_stream.end() );
      Utility::RandomNumberGenerator::setFakeStream( fake_stream );
      spatial_index_dimension_distribution->sampleWithoutCascade( point );
      FRENSIE_CHECK_EQUAL( getCoordinate<SPATIAL_INDEX_DIMENSION>( point ), expected_samples[sample] );
      FRENSIE_CHECK_FLOATING_EQUALITY( getCoordinateWeight<SPATIAL_INDEX_DIMENSION>( point ), expected_weights[sample], 1e-14);
      Utility::RandomNumberGenerator::unsetFakeStream();
    }
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
    // First random number is used to sample the index. The ones after that are for determining where in that element

    std::vector<size_t> expected_samples = {0, 1, 1};
    std::vector<double> expected_weights = {(2.0/5.0)/(1.0/6.0), (3.0/5.0)/(5.0/6.0), (3.0/5.0)/(5.0/6.0)};
    std::vector<double> index_fake_stream = {0.0, 0.5, 1-1e-15};
    std::vector<double> other_fake_stream = {0.5, 0.5, 0.5};
    for(size_t sample = 0; sample < expected_samples.size(); ++sample)
    {
      std::vector<double> fake_stream(4);
      fake_stream[0] = index_fake_stream[sample];
      fake_stream.insert( fake_stream.end(), other_fake_stream.begin(), other_fake_stream.end() );
      Utility::RandomNumberGenerator::setFakeStream( fake_stream );
      direction_index_dimension_distribution->sampleWithoutCascade( point );
      FRENSIE_CHECK_EQUAL( getCoordinate<DIRECTION_INDEX_DIMENSION>( point ), expected_samples[sample] );
      FRENSIE_CHECK_FLOATING_EQUALITY( getCoordinateWeight<DIRECTION_INDEX_DIMENSION>( point ), expected_weights[sample], 1e-14);
      Utility::RandomNumberGenerator::unsetFakeStream();
    }
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
// end tstImportanceSampledIndependentPhaseSpaceDimensionDistribution_indexDimensions.cpp
//---------------------------------------------------------------------------//
