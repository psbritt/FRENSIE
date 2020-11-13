//---------------------------------------------------------------------------//
//!
//! \file   tstImportanceSampledDependentPhaseSpaceDimensionDistribution_indexDimensions.cpp
//! \author Alex Robinson
//! \brief  Importance sampled dependent phase space dimension dist. unit tests
//!
//---------------------------------------------------------------------------//

// Std Lib Includes
#include <iostream>
#include <memory>

// FRENSIE Includes
#include "MonteCarlo_ImportanceSampledDependentPhaseSpaceDimensionDistribution.hpp"
#include "MonteCarlo_PhaseSpaceDimensionTraits.hpp"
#include "Utility_BasicCartesianCoordinateConversionPolicy.hpp"
#include "Utility_HistogramDistribution.hpp"
#include "Utility_HistogramFullyTabularBasicBivariateDistribution.hpp"
#include "Utility_RandomNumberGenerator.hpp"
#include "Utility_UnitTestHarnessWithMain.hpp"
#include "ArchiveTestHelpers.hpp"

//---------------------------------------------------------------------------//
// Testing Types
//---------------------------------------------------------------------------//

using namespace MonteCarlo;

typedef std::tuple<
  /* primary spatial dimension */
  std::tuple<std::integral_constant<PhaseSpaceDimension,SPATIAL_INDEX_DIMENSION>,
             std::integral_constant<PhaseSpaceDimension,DIRECTION_INDEX_DIMENSION> >,
  std::tuple<std::integral_constant<PhaseSpaceDimension,DIRECTION_INDEX_DIMENSION>,
             std::integral_constant<PhaseSpaceDimension,SPATIAL_INDEX_DIMENSION> >
 > TestPhaseSpaceDimensions;

typedef TestArchiveHelper::TestArchives TestArchives;

//---------------------------------------------------------------------------//
// Testing Variables
//---------------------------------------------------------------------------//
std::shared_ptr<const Utility::SpatialCoordinateConversionPolicy>
spatial_coord_conversion_policy( new Utility::BasicCartesianCoordinateConversionPolicy );

std::shared_ptr<const Utility::DirectionalCoordinateConversionPolicy>
directional_coord_conversion_policy( new Utility::BasicCartesianCoordinateConversionPolicy );

std::shared_ptr<const Utility::BasicBivariateDistribution> raw_distribution;

std::shared_ptr<const Utility::BasicBivariateDistribution>
raw_importance_distribution;

//---------------------------------------------------------------------------//
// Tests.
//---------------------------------------------------------------------------//
// Test that the dimension can be returned
FRENSIE_UNIT_TEST_TEMPLATE_EXPAND(
                     ImportanceSampledDependentPhaseSpaceDimensionDistribution_indexDimensions,
                     getDimension,
                     TestPhaseSpaceDimensions )
{
  FETCH_TEMPLATE_PARAM( 0, WrappedParentDimension );
  FETCH_TEMPLATE_PARAM( 1, WrappedDimension );

  constexpr PhaseSpaceDimension ParentDimension = WrappedParentDimension::value;
  constexpr PhaseSpaceDimension Dimension = WrappedDimension::value;

  std::unique_ptr<const MonteCarlo::PhaseSpaceDimensionDistribution>
    dimension_distribution( new MonteCarlo::ImportanceSampledDependentPhaseSpaceDimensionDistribution<ParentDimension,Dimension>(
                                               raw_distribution,
                                               raw_importance_distribution ) );

  FRENSIE_CHECK_EQUAL( dimension_distribution->getDimension(), Dimension );
}

//---------------------------------------------------------------------------//
// Test that the dimension class can be returned
FRENSIE_UNIT_TEST_TEMPLATE_EXPAND(
                     ImportanceSampledDependentPhaseSpaceDimensionDistribution_indexDimensions,
                     getDimensionClass,
                     TestPhaseSpaceDimensions )
{
  FETCH_TEMPLATE_PARAM( 0, WrappedParentDimension );
  FETCH_TEMPLATE_PARAM( 1, WrappedDimension );

  constexpr PhaseSpaceDimension ParentDimension = WrappedParentDimension::value;
  constexpr PhaseSpaceDimension Dimension = WrappedDimension::value;

  std::unique_ptr<const MonteCarlo::PhaseSpaceDimensionDistribution>
    dimension_distribution( new MonteCarlo::ImportanceSampledDependentPhaseSpaceDimensionDistribution<ParentDimension,Dimension>(
                             raw_distribution, raw_importance_distribution ) );

  FRENSIE_CHECK_EQUAL( dimension_distribution->getDimensionClass(),
                       MonteCarlo::PhaseSpaceDimensionTraits<Dimension>::getClass() );
}

//---------------------------------------------------------------------------//
// Test that the indep. dimension can be returned
FRENSIE_UNIT_TEST_TEMPLATE_EXPAND(
                     ImportanceSampledDependentPhaseSpaceDimensionDistribution_indexDimensions,
                     getParentDimension,
                     TestPhaseSpaceDimensions )
{
  FETCH_TEMPLATE_PARAM( 0, WrappedParentDimension );
  FETCH_TEMPLATE_PARAM( 1, WrappedDimension );

  constexpr PhaseSpaceDimension ParentDimension = WrappedParentDimension::value;
  constexpr PhaseSpaceDimension Dimension = WrappedDimension::value;

  std::unique_ptr<const MonteCarlo::DependentPhaseSpaceDimensionDistribution<ParentDimension,Dimension> >
    dimension_distribution( new MonteCarlo::ImportanceSampledDependentPhaseSpaceDimensionDistribution<ParentDimension,Dimension>(
                             raw_distribution, raw_importance_distribution ) );

  FRENSIE_CHECK_EQUAL( dimension_distribution->getParentDimension(),
                       ParentDimension );
}

//---------------------------------------------------------------------------//
// Test that the indep. dimension class can be returned
FRENSIE_UNIT_TEST_TEMPLATE_EXPAND(
                     ImportanceSampledDependentPhaseSpaceDimensionDistribution_indexDimensions,
                     getParentDimensionClass,
                     TestPhaseSpaceDimensions )
{
  FETCH_TEMPLATE_PARAM( 0, WrappedParentDimension );
  FETCH_TEMPLATE_PARAM( 1, WrappedDimension );

  constexpr PhaseSpaceDimension ParentDimension = WrappedParentDimension::value;
  constexpr PhaseSpaceDimension Dimension = WrappedDimension::value;

  std::unique_ptr<const MonteCarlo::DependentPhaseSpaceDimensionDistribution<ParentDimension,Dimension> >
    dimension_distribution( new MonteCarlo::ImportanceSampledDependentPhaseSpaceDimensionDistribution<ParentDimension,Dimension>(
                             raw_distribution, raw_importance_distribution ) );

  FRENSIE_CHECK_EQUAL( dimension_distribution->getParentDimensionClass(),
                       MonteCarlo::PhaseSpaceDimensionTraits<ParentDimension>::getClass() );
}

//---------------------------------------------------------------------------//
// Test that the distribution type name can be returned
FRENSIE_UNIT_TEST_TEMPLATE_EXPAND(
                     ImportanceSampledDependentPhaseSpaceDimensionDistribution_indexDimensions,
                     getDistributionTypeName,
                     TestPhaseSpaceDimensions )
{
  FETCH_TEMPLATE_PARAM( 0, WrappedParentDimension );
  FETCH_TEMPLATE_PARAM( 1, WrappedDimension );

  constexpr PhaseSpaceDimension ParentDimension = WrappedParentDimension::value;
  constexpr PhaseSpaceDimension Dimension = WrappedDimension::value;

  std::unique_ptr<const MonteCarlo::PhaseSpaceDimensionDistribution>
    dimension_distribution( new MonteCarlo::ImportanceSampledDependentPhaseSpaceDimensionDistribution<ParentDimension,Dimension>(
                             raw_distribution, raw_importance_distribution ) );

  FRENSIE_CHECK_EQUAL( dimension_distribution->getDistributionTypeName(),
                       "BasicBivariateDistribution" );
}

//---------------------------------------------------------------------------//
// Test that the distribution type name can be returned
FRENSIE_UNIT_TEST_TEMPLATE_EXPAND(
                     ImportanceSampledDependentPhaseSpaceDimensionDistribution_indexDimensions,
                     isIndependent,
                     TestPhaseSpaceDimensions )
{
  FETCH_TEMPLATE_PARAM( 0, WrappedParentDimension );
  FETCH_TEMPLATE_PARAM( 1, WrappedDimension );

  constexpr PhaseSpaceDimension ParentDimension = WrappedParentDimension::value;
  constexpr PhaseSpaceDimension Dimension = WrappedDimension::value;

  std::unique_ptr<const MonteCarlo::PhaseSpaceDimensionDistribution>
    dimension_distribution( new MonteCarlo::ImportanceSampledDependentPhaseSpaceDimensionDistribution<ParentDimension,Dimension>(
                             raw_distribution, raw_importance_distribution ) );

  FRENSIE_CHECK( !dimension_distribution->isIndependent() );
}

//---------------------------------------------------------------------------//
// Test if the distribution is dependent on another dimension
FRENSIE_UNIT_TEST_TEMPLATE_EXPAND(
                     ImportanceSampledDependentPhaseSpaceDimensionDistribution_indexDimensions,
                     isDependentOnDimension,
                     TestPhaseSpaceDimensions )
{
  FETCH_TEMPLATE_PARAM( 0, WrappedParentDimension );
  FETCH_TEMPLATE_PARAM( 1, WrappedDimension );

  constexpr PhaseSpaceDimension ParentDimension = WrappedParentDimension::value;
  constexpr PhaseSpaceDimension Dimension = WrappedDimension::value;

  std::unique_ptr<const MonteCarlo::PhaseSpaceDimensionDistribution>
    dimension_distribution( new MonteCarlo::ImportanceSampledDependentPhaseSpaceDimensionDistribution<ParentDimension,Dimension>(
                             raw_distribution, raw_importance_distribution ) );

  FRENSIE_CHECK( dimension_distribution->isDependentOnDimension( ParentDimension ) );
  FRENSIE_CHECK( !dimension_distribution->isDependentOnDimension( Dimension ) );
}

//---------------------------------------------------------------------------//
// Test if the distribution is continuous
FRENSIE_UNIT_TEST_TEMPLATE_EXPAND(
                     ImportanceSampledDependentPhaseSpaceDimensionDistribution_indexDimensions,
                     isContinuous,
                     TestPhaseSpaceDimensions )
{
  FETCH_TEMPLATE_PARAM( 0, WrappedParentDimension );
  FETCH_TEMPLATE_PARAM( 1, WrappedDimension );

  constexpr PhaseSpaceDimension ParentDimension = WrappedParentDimension::value;
  constexpr PhaseSpaceDimension Dimension = WrappedDimension::value;

  std::unique_ptr<const MonteCarlo::PhaseSpaceDimensionDistribution>
    dimension_distribution( new MonteCarlo::ImportanceSampledDependentPhaseSpaceDimensionDistribution<ParentDimension,Dimension>(
                             raw_distribution, raw_importance_distribution ) );

  FRENSIE_CHECK( dimension_distribution->isContinuous() );
}

//---------------------------------------------------------------------------//
// Test if the distribution is tabular
FRENSIE_UNIT_TEST_TEMPLATE_EXPAND(
                     ImportanceSampledDependentPhaseSpaceDimensionDistribution_indexDimensions,
                     isTabular,
                     TestPhaseSpaceDimensions )
{
  FETCH_TEMPLATE_PARAM( 0, WrappedParentDimension );
  FETCH_TEMPLATE_PARAM( 1, WrappedDimension );

  constexpr PhaseSpaceDimension ParentDimension = WrappedParentDimension::value;
  constexpr PhaseSpaceDimension Dimension = WrappedDimension::value;

  std::unique_ptr<const MonteCarlo::PhaseSpaceDimensionDistribution>
    dimension_distribution( new MonteCarlo::ImportanceSampledDependentPhaseSpaceDimensionDistribution<ParentDimension,Dimension>(
                             raw_distribution, raw_importance_distribution ) );

  FRENSIE_CHECK( dimension_distribution->isTabular() );

}

//---------------------------------------------------------------------------//
// Test if the distribution is uniform
FRENSIE_UNIT_TEST_TEMPLATE_EXPAND(
                     ImportanceSampledDependentPhaseSpaceDimensionDistribution_indexDimensions,
                     isUniform,
                     TestPhaseSpaceDimensions )
{
  FETCH_TEMPLATE_PARAM( 0, WrappedParentDimension );
  FETCH_TEMPLATE_PARAM( 1, WrappedDimension );

  constexpr PhaseSpaceDimension ParentDimension = WrappedParentDimension::value;
  constexpr PhaseSpaceDimension Dimension = WrappedDimension::value;

  std::unique_ptr<const MonteCarlo::PhaseSpaceDimensionDistribution>
    dimension_distribution( new MonteCarlo::ImportanceSampledDependentPhaseSpaceDimensionDistribution<ParentDimension,Dimension>(
                             raw_distribution, raw_importance_distribution ) );

  FRENSIE_CHECK( !dimension_distribution->isUniform() );

}

//---------------------------------------------------------------------------//
// Test if the distribution is uniform
FRENSIE_UNIT_TEST_TEMPLATE_EXPAND(
                     ImportanceSampledDependentPhaseSpaceDimensionDistribution_indexDimensions,
                     hasForm,
                     TestPhaseSpaceDimensions )
{
  FETCH_TEMPLATE_PARAM( 0, WrappedParentDimension );
  FETCH_TEMPLATE_PARAM( 1, WrappedDimension );

  constexpr PhaseSpaceDimension ParentDimension = WrappedParentDimension::value;
  constexpr PhaseSpaceDimension Dimension = WrappedDimension::value;

  std::unique_ptr<const MonteCarlo::PhaseSpaceDimensionDistribution>
    dimension_distribution( new MonteCarlo::ImportanceSampledDependentPhaseSpaceDimensionDistribution<ParentDimension,Dimension>(
                             raw_distribution, raw_importance_distribution ) );

  FRENSIE_CHECK( !dimension_distribution->hasForm( Utility::DELTA_DISTRIBUTION) );
  FRENSIE_CHECK( !dimension_distribution->hasForm( Utility::UNIFORM_DISTRIBUTION ) );

}

//---------------------------------------------------------------------------//
// Test if the distribution can be evaluated without a cascade
FRENSIE_UNIT_TEST_TEMPLATE_EXPAND(
                     ImportanceSampledDependentPhaseSpaceDimensionDistribution_indexDimensions,
                     evaluateWithoutCascade,
                     TestPhaseSpaceDimensions )
{
  FETCH_TEMPLATE_PARAM( 0, WrappedParentDimension );
  FETCH_TEMPLATE_PARAM( 1, WrappedDimension );

  constexpr PhaseSpaceDimension ParentDimension = WrappedParentDimension::value;
  constexpr PhaseSpaceDimension Dimension = WrappedDimension::value;

  std::unique_ptr<const MonteCarlo::PhaseSpaceDimensionDistribution>
    index_index_dimension_distribution( new MonteCarlo::ImportanceSampledDependentPhaseSpaceDimensionDistribution<ParentDimension,Dimension>( raw_distribution, raw_importance_distribution ) );

  MonteCarlo::PhaseSpacePoint point( spatial_coord_conversion_policy,
                                     directional_coord_conversion_policy );

  std::vector<double> independent_values {0, 0.5, 1, 2-1e-14};
  std::vector<double> second_independent_values {0, 0.5, 1, 2-1e-14};
  std::vector<std::vector<double>> expected_results {{2.0, 2.0, 3.0, 3.0}, {2.0, 2.0, 3.0, 3.0}, {4.0, 4.0, 5.0, 5.0}, {4.0, 4.0, 5.0, 5.0}};

  for(size_t ind_value = 0; ind_value < independent_values.size(); ++ind_value)
  {
    setCoordinate<ParentDimension>(point, independent_values[ind_value]);
    for(size_t sec_ind_value = 0; sec_ind_value < second_independent_values.size(); ++sec_ind_value)
    {
      setCoordinate<Dimension>(point, second_independent_values[sec_ind_value]);
      FRENSIE_CHECK_EQUAL(index_index_dimension_distribution->evaluateWithoutCascade( point ), expected_results[ind_value][sec_ind_value]);
    }
    
    
  }
  // Should throw exception when index is out of bounds to counter roundoff errors resulting in undefined behavior
  FRENSIE_CHECK_THROW(setCoordinate<Dimension>( point, 72 ), std::runtime_error);
}

//---------------------------------------------------------------------------//
// Test if the distribution can be sampled without a cascade
FRENSIE_UNIT_TEST_TEMPLATE_EXPAND(
                     ImportanceSampledDependentPhaseSpaceDimensionDistribution_indexDimensions,
                     sampleWithoutCascade,
                     TestPhaseSpaceDimensions )
{
  FETCH_TEMPLATE_PARAM( 0, WrappedParentDimension );
  FETCH_TEMPLATE_PARAM( 1, WrappedDimension );

  constexpr PhaseSpaceDimension ParentDimension = WrappedParentDimension::value;
  constexpr PhaseSpaceDimension Dimension = WrappedDimension::value;

  std::unique_ptr<const MonteCarlo::PhaseSpaceDimensionDistribution>
    dimension_distribution( new MonteCarlo::ImportanceSampledDependentPhaseSpaceDimensionDistribution<ParentDimension,Dimension>( raw_distribution, raw_importance_distribution ) );

  MonteCarlo::PhaseSpacePoint point( spatial_coord_conversion_policy,
                                     directional_coord_conversion_policy );

  std::vector<double> index_values = {0, 1};
  std::vector<double> initial_fake_stream = {0.0, 0.5, 1-1e-15};
  std::vector<double> index_additional_fake_stream = {0.5, 0.5, 0.5};
  std::vector<std::vector<size_t>> results = {{0, 1, 1}, {0, 1, 1}};
  std::vector<std::vector<double>> expected_weights = {{(2.0/5.0)/(6.0/13.0), (3.0/5.0)/(7.0/13.0), (3.0/5.0)/(7.0/13.0)},
                                                       {(4.0/9.0)/(8.0/17.0), (5.0/9.0)/(9.0/17.0), (5.0/9.0)/(9.0/17.0)}};

  for(size_t index = 0; index < index_values.size(); ++index)
  {
    setCoordinate<ParentDimension>( point, index_values[index] );
    for(size_t sample = 0; sample < initial_fake_stream.size(); ++sample)
    {
      std::vector<double> fake_stream(4);
      fake_stream[0] = initial_fake_stream[sample];
      fake_stream.insert( fake_stream.end(), index_additional_fake_stream.begin(), index_additional_fake_stream.end() );
      Utility::RandomNumberGenerator::setFakeStream( fake_stream );
      dimension_distribution->sampleWithoutCascade( point );
      FRENSIE_CHECK_FLOATING_EQUALITY(getCoordinate<Dimension>( point ), results[index][sample], 1e-15);
      FRENSIE_CHECK_FLOATING_EQUALITY(getCoordinateWeight<Dimension>( point ), expected_weights[index][sample], 1e-15);
    }
    Utility::RandomNumberGenerator::unsetFakeStream();
  }
}

//---------------------------------------------------------------------------//
// Test if the distribution can be sampled without a cascade
FRENSIE_UNIT_TEST_TEMPLATE_EXPAND(
                     ImportanceSampledDependentPhaseSpaceDimensionDistribution_indexDimensions,
                     sampleAndRecordTrialsWithoutCascade,
                     TestPhaseSpaceDimensions )
{
  FETCH_TEMPLATE_PARAM( 0, WrappedParentDimension );
  FETCH_TEMPLATE_PARAM( 1, WrappedDimension );

  constexpr PhaseSpaceDimension ParentDimension = WrappedParentDimension::value;
  constexpr PhaseSpaceDimension Dimension = WrappedDimension::value;

  std::unique_ptr<const MonteCarlo::PhaseSpaceDimensionDistribution>
    dimension_distribution( new MonteCarlo::ImportanceSampledDependentPhaseSpaceDimensionDistribution<ParentDimension,Dimension>( raw_distribution, raw_importance_distribution ) );

  MonteCarlo::PhaseSpacePoint point( spatial_coord_conversion_policy,
                                     directional_coord_conversion_policy );

  // First random number is used to sample the index. The ones after that are for determining where in that element

  MonteCarlo::PhaseSpaceDimensionDistribution::Counter trials = 0;

  std::vector<double> index_values = {0.5, 1.0};
  std::vector<double> initial_fake_stream = {0.0, 0.5, 1-1e-15};
  std::vector<double> index_additional_fake_stream = {0.5, 0.5, 0.5};
  std::vector<std::vector<size_t>> results = {{0, 1, 1}, {0, 1, 1}};
  std::vector<std::vector<double>> expected_weights = {{(2.0/5.0)/(6.0/13.0), (3.0/5.0)/(7.0/13.0), (3.0/5.0)/(7.0/13.0)},
                                                       {(4.0/9.0)/(8.0/17.0), (5.0/9.0)/(9.0/17.0), (5.0/9.0)/(9.0/17.0)}};

  size_t counter = 0;
  for(size_t index = 0; index < index_values.size(); ++index)
  {
    setCoordinate<ParentDimension>( point, index_values[index] );
    for(size_t sample = 0; sample < initial_fake_stream.size(); ++sample)
    {
      std::vector<double> fake_stream(4);
      fake_stream[0] = initial_fake_stream[sample];
      fake_stream.insert( fake_stream.end(), index_additional_fake_stream.begin(), index_additional_fake_stream.end() );
      Utility::RandomNumberGenerator::setFakeStream( fake_stream );
      dimension_distribution->sampleAndRecordTrialsWithoutCascade( point, trials );
      FRENSIE_CHECK_EQUAL(getCoordinate<Dimension>( point ), results[index][sample]);
      FRENSIE_CHECK_FLOATING_EQUALITY(getCoordinateWeight<Dimension>( point ), expected_weights[index][sample], 1e-15);
      counter += 1;
      FRENSIE_CHECK_EQUAL(trials, counter);
    }
    Utility::RandomNumberGenerator::unsetFakeStream();
  }

}

//---------------------------------------------------------------------------//
// Test that the dimension value can be set and weighed approperiately
FRENSIE_UNIT_TEST_TEMPLATE_EXPAND(
                     ImportanceSampledDependentPhaseSpaceDimensionDistribution_indexDimensions,
                     setDimensionValueAndApplyWeight,
                     TestPhaseSpaceDimensions )
{
  FETCH_TEMPLATE_PARAM( 0, WrappedParentDimension );
  FETCH_TEMPLATE_PARAM( 1, WrappedDimension );

  constexpr PhaseSpaceDimension ParentDimension = WrappedParentDimension::value;
  constexpr PhaseSpaceDimension Dimension = WrappedDimension::value;

  std::unique_ptr<const MonteCarlo::PhaseSpaceDimensionDistribution>
    dimension_distribution( new MonteCarlo::ImportanceSampledDependentPhaseSpaceDimensionDistribution<ParentDimension,Dimension>( raw_distribution, raw_importance_distribution ) );

  MonteCarlo::PhaseSpacePoint point( spatial_coord_conversion_policy,
                                     directional_coord_conversion_policy );

  setCoordinate<ParentDimension>( point, 0 );

  dimension_distribution->setDimensionValueAndApplyWeight( point, 0 );

  FRENSIE_CHECK_EQUAL( getCoordinate<Dimension>( point ), 0);
  FRENSIE_CHECK_FLOATING_EQUALITY(getCoordinateWeight<Dimension>( point ), 2.0/5.0, 1e-14);

  dimension_distribution->setDimensionValueAndApplyWeight( point, 1 );

  FRENSIE_CHECK_EQUAL( getCoordinate<Dimension>( point ), 1);
  FRENSIE_CHECK_FLOATING_EQUALITY(getCoordinateWeight<Dimension>( point ), 3.0/5.0, 1e-14);

  setCoordinate<ParentDimension>( point, 1 );

  dimension_distribution->setDimensionValueAndApplyWeight( point, 0 );

  FRENSIE_CHECK_EQUAL( getCoordinate<Dimension>( point ), 0);
  FRENSIE_CHECK_FLOATING_EQUALITY(getCoordinateWeight<Dimension>( point ), 4.0/9.0, 1e-14);

  dimension_distribution->setDimensionValueAndApplyWeight( point, 1 );

  FRENSIE_CHECK_EQUAL( getCoordinate<Dimension>( point ), 1);
  FRENSIE_CHECK_FLOATING_EQUALITY(getCoordinateWeight<Dimension>( point ), 5.0/9.0, 1e-14);

}

//---------------------------------------------------------------------------//
// Check that the distribution can be archived
FRENSIE_UNIT_TEST_TEMPLATE_EXPAND( ImportanceSampledDependentPhaseSpaceDimensionDistribution_indexDimensions,
                                   archive,
                                   TestArchives )
{
  FETCH_TEMPLATE_PARAM( 0, RawOArchive );
  FETCH_TEMPLATE_PARAM( 1, RawIArchive );

  typedef typename std::remove_pointer<RawOArchive>::type OArchive;
  typedef typename std::remove_pointer<RawIArchive>::type IArchive;

  std::string archive_base_name( "test_importance_sampled_dependent_phase_dimension_distribution_index_dimensions" );
  std::ostringstream archive_ostream;

  {
    std::unique_ptr<OArchive> oarchive;

    createOArchive( archive_base_name, archive_ostream, oarchive );

    // Create the raw distribution
    std::shared_ptr<const Utility::BasicBivariateDistribution>
      raw_distribution;

    {
      std::vector<double> grid ({0, 1, 2});

      std::vector<std::shared_ptr<const Utility::TabularUnivariateDistribution> >
        secondary_dists( 3 );
      
      std::vector<double> secondary_dep_values( {2.0, 3.0} );
      // Create the secondary distribution in the first bin
      secondary_dists[0].reset( new Utility::HistogramDistribution( grid, secondary_dep_values ) );
      
      secondary_dep_values[0] = 4.0;
      secondary_dep_values[1] = 5.0;
      // Create the secondary distribution in the second bin
      secondary_dists[1].reset( new Utility::HistogramDistribution( grid, secondary_dep_values ) );

      secondary_dists[2] = secondary_dists[0];

      std::shared_ptr<Utility::HistogramFullyTabularBasicBivariateDistribution> local_raw_distribution = 
        std::make_shared<Utility::HistogramFullyTabularBasicBivariateDistribution>(grid, secondary_dists);

      local_raw_distribution->limitToPrimaryIndepLimits();
    
      raw_distribution = local_raw_distribution ;
    }

    // Create the importance distribution
    std::shared_ptr<const Utility::BasicBivariateDistribution>
      importance_distribution;
    
    {
      std::vector<double> grid ({0, 1, 2});

      std::vector<std::shared_ptr<const Utility::TabularUnivariateDistribution> >
        secondary_dists( 3 );
      
      std::vector<double> secondary_dep_values( {6.0, 7.0} );
      // Create the secondary distribution in the first bin
      secondary_dists[0].reset( new Utility::HistogramDistribution( grid, secondary_dep_values ) );
      
      secondary_dep_values[0] = 8.0;
      secondary_dep_values[1] = 9.0;
      // Create the secondary distribution in the second bin
      secondary_dists[1].reset( new Utility::HistogramDistribution( grid, secondary_dep_values ) );

      secondary_dists[2] = secondary_dists[0];

      std::shared_ptr<Utility::HistogramFullyTabularBasicBivariateDistribution> local_raw_distribution = 
        std::make_shared<Utility::HistogramFullyTabularBasicBivariateDistribution>(grid, secondary_dists);

      local_raw_distribution->limitToPrimaryIndepLimits();
    
      importance_distribution = local_raw_distribution ;
    }

    // Create the dependent primary spatial distributions
    std::shared_ptr<const MonteCarlo::PhaseSpaceDimensionDistribution>
      spatial_index_dependent_direction_index_dimension_distribution( new MonteCarlo::ImportanceSampledSpatialIndexDependentDirectionIndexDimensionDistribution( raw_distribution, importance_distribution ) );

    std::shared_ptr<const MonteCarlo::PhaseSpaceDimensionDistribution>
      direction_index_dependent_spatial_index_dimension_distribution( new MonteCarlo::ImportanceSampledDirectionIndexDependentSpatialIndexDimensionDistribution( raw_distribution, importance_distribution ) );

    FRENSIE_REQUIRE_NO_THROW( (*oarchive) << BOOST_SERIALIZATION_NVP(spatial_index_dependent_direction_index_dimension_distribution) );
    FRENSIE_REQUIRE_NO_THROW( (*oarchive) << BOOST_SERIALIZATION_NVP(direction_index_dependent_spatial_index_dimension_distribution) );
  }

  // Copy the archive ostream to an istream
  std::istringstream archive_istream( archive_ostream.str() );

  // Load the archived distributions
  std::unique_ptr<IArchive> iarchive;

  createIArchive( archive_istream, iarchive );

  std::shared_ptr<const MonteCarlo::PhaseSpaceDimensionDistribution>
    spatial_index_dependent_direction_index_dimension_distribution,
    direction_index_dependent_spatial_index_dimension_distribution;

  FRENSIE_REQUIRE_NO_THROW( (*iarchive) >> BOOST_SERIALIZATION_NVP(spatial_index_dependent_direction_index_dimension_distribution) );
  FRENSIE_REQUIRE_NO_THROW( (*iarchive) >> BOOST_SERIALIZATION_NVP(direction_index_dependent_spatial_index_dimension_distribution) );

  iarchive.reset();

  {
    FRENSIE_CHECK_EQUAL( spatial_index_dependent_direction_index_dimension_distribution->getDimension(),
                          DIRECTION_INDEX_DIMENSION );

    auto concrete_spatial_index_dependent_direction_index_dimension_distribution = dynamic_cast<const MonteCarlo::ImportanceSampledSpatialIndexDependentDirectionIndexDimensionDistribution*>( spatial_index_dependent_direction_index_dimension_distribution.get() );

    FRENSIE_CHECK_EQUAL( concrete_spatial_index_dependent_direction_index_dimension_distribution->getParentDimension(),
                          SPATIAL_INDEX_DIMENSION );

    MonteCarlo::PhaseSpacePoint point( spatial_coord_conversion_policy,
                                        directional_coord_conversion_policy );

    std::vector<double> index_values = {0, 1};
    std::vector<double> initial_fake_stream = {0.0, 0.5, 1-1e-15};
    std::vector<double> index_additional_fake_stream = {0.5, 0.5, 0.5};
    std::vector<std::vector<size_t>> results = {{0, 1, 1}, {0, 1, 1}};
    std::vector<std::vector<double>> expected_weights = {{(2.0/5.0)/(6.0/13.0), (3.0/5.0)/(7.0/13.0), (3.0/5.0)/(7.0/13.0)},
                                                        {(4.0/9.0)/(8.0/17.0), (5.0/9.0)/(9.0/17.0), (5.0/9.0)/(9.0/17.0)}};

    for(size_t index = 0; index < index_values.size(); ++index)
    {
      setCoordinate<SPATIAL_INDEX_DIMENSION>( point, index_values[index] );
      for(size_t sample = 0; sample < initial_fake_stream.size(); ++sample)
      {
        std::vector<double> fake_stream(4);
        fake_stream[0] = initial_fake_stream[sample];
        fake_stream.insert( fake_stream.end(), index_additional_fake_stream.begin(), index_additional_fake_stream.end() );
        Utility::RandomNumberGenerator::setFakeStream( fake_stream );
        spatial_index_dependent_direction_index_dimension_distribution->sampleWithoutCascade( point );
        FRENSIE_CHECK_FLOATING_EQUALITY(getCoordinate<DIRECTION_INDEX_DIMENSION>( point ), results[index][sample], 1e-15);
        FRENSIE_CHECK_FLOATING_EQUALITY(getCoordinateWeight<DIRECTION_INDEX_DIMENSION>( point ), expected_weights[index][sample], 1e-15);
      }
      Utility::RandomNumberGenerator::unsetFakeStream();
    }
  }

  {
    FRENSIE_CHECK_EQUAL( direction_index_dependent_spatial_index_dimension_distribution->getDimension(),
                          SPATIAL_INDEX_DIMENSION );

    auto concrete_direction_index_dependent_spatial_index_dimension_distribution = dynamic_cast<const MonteCarlo::ImportanceSampledDirectionIndexDependentSpatialIndexDimensionDistribution*>( direction_index_dependent_spatial_index_dimension_distribution.get() );

    FRENSIE_CHECK_EQUAL( concrete_direction_index_dependent_spatial_index_dimension_distribution->getParentDimension(),
                          DIRECTION_INDEX_DIMENSION );

    MonteCarlo::PhaseSpacePoint point( spatial_coord_conversion_policy,
                                        directional_coord_conversion_policy );

    std::vector<double> index_values = {0, 1};
    std::vector<double> initial_fake_stream = {0.0, 0.5, 1-1e-15};
    std::vector<double> index_additional_fake_stream = {0.5, 0.5, 0.5};
    std::vector<std::vector<size_t>> results = {{0, 1, 1}, {0, 1, 1}};
    std::vector<std::vector<double>> expected_weights = {{(2.0/5.0)/(6.0/13.0), (3.0/5.0)/(7.0/13.0), (3.0/5.0)/(7.0/13.0)},
                                                        {(4.0/9.0)/(8.0/17.0), (5.0/9.0)/(9.0/17.0), (5.0/9.0)/(9.0/17.0)}};

    for(size_t index = 0; index < index_values.size(); ++index)
    {
      setCoordinate<DIRECTION_INDEX_DIMENSION>( point, index_values[index] );
      for(size_t sample = 0; sample < initial_fake_stream.size(); ++sample)
      {
        std::vector<double> fake_stream(4);
        fake_stream[0] = initial_fake_stream[sample];
        fake_stream.insert( fake_stream.end(), index_additional_fake_stream.begin(), index_additional_fake_stream.end() );
        Utility::RandomNumberGenerator::setFakeStream( fake_stream );
        direction_index_dependent_spatial_index_dimension_distribution->sampleWithoutCascade( point );
        FRENSIE_CHECK_FLOATING_EQUALITY(getCoordinate<SPATIAL_INDEX_DIMENSION>( point ), results[index][sample], 1e-15);
        FRENSIE_CHECK_FLOATING_EQUALITY(getCoordinateWeight<SPATIAL_INDEX_DIMENSION>( point ), expected_weights[index][sample], 1e-15);
      }
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
  std::vector<double> primary_grid( {0, 1, 2} );
  // Create the distribution
  {
    std::vector<std::shared_ptr<const Utility::TabularUnivariateDistribution> >
      secondary_dists( 3 );

    std::vector<double> secondary_ind_grid( {0, 1, 2} );
    
    std::vector<double> secondary_dep_values( {2.0, 3.0} );
    // Create the secondary distribution in the first bin
    secondary_dists[0].reset( new Utility::HistogramDistribution( secondary_ind_grid, secondary_dep_values ) );
    
    secondary_dep_values[0] = 4.0;
    secondary_dep_values[1] = 5.0;
    // Create the secondary distribution in the second bin
    secondary_dists[1].reset( new Utility::HistogramDistribution( secondary_ind_grid, secondary_dep_values ) );

    // Create the secondary distribution in the third bin
    secondary_dists[2] = secondary_dists[1];

    Utility::HistogramFullyTabularBasicBivariateDistribution* local_raw_distribution =
      new Utility::HistogramFullyTabularBasicBivariateDistribution( primary_grid, secondary_dists );

    local_raw_distribution->limitToPrimaryIndepLimits();
    
    raw_distribution.reset( local_raw_distribution );
  }

  // Create the fully tabular importance distribution
  {
    std::vector<std::shared_ptr<const Utility::TabularUnivariateDistribution> >
      secondary_dists( 3 );

    std::vector<double> secondary_ind_grid( {0, 1, 2} );
    
    std::vector<double> secondary_dep_values( {6.0, 7.0} );

    // Create the secondary distribution in the first bin
    secondary_dists[0].reset( new Utility::HistogramDistribution( secondary_ind_grid, secondary_dep_values ) );

    secondary_dep_values[0] = 8.0;
    secondary_dep_values[1] = 9.0;

    // Create the secondary distribution in the second bin
    secondary_dists[1].reset( new Utility::HistogramDistribution( secondary_ind_grid, secondary_dep_values ) );
    
    // Create the secondary distribution in the third bin
    secondary_dists[2] = secondary_dists[1];
    
    Utility::HistogramFullyTabularBasicBivariateDistribution*
      local_raw_importance_distribution =
      new Utility::HistogramFullyTabularBasicBivariateDistribution( primary_grid, secondary_dists );
    
    local_raw_importance_distribution->limitToPrimaryIndepLimits();
    
    raw_importance_distribution.reset( local_raw_importance_distribution );
  }
    
  // Initialize the random number generator
  Utility::RandomNumberGenerator::createStreams();
}

FRENSIE_CUSTOM_UNIT_TEST_SETUP_END();

//---------------------------------------------------------------------------//
// end tstImportanceSampledDependentPhaseSpaceDimensionDistribution_indexDimensions.cpp
//---------------------------------------------------------------------------//
