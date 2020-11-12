//---------------------------------------------------------------------------//
//!
//! \file   tstDependentPhaseSpaceDimensionDistribution.cpp
//! \author Alex Robinson
//! \brief  Dependent phase space dimension dist. unit tests
//!
//---------------------------------------------------------------------------//

// Std Lib Includes
#include <iostream>
#include <memory>

// FRENSIE Includes
#include "MonteCarlo_DependentPhaseSpaceDimensionDistribution.hpp"
#include "MonteCarlo_PhaseSpaceDimensionTraits.hpp"
#include "Utility_BasicCartesianCoordinateConversionPolicy.hpp"
#include "Utility_DeltaDistribution.hpp"
#include "Utility_HistogramDistribution.hpp"
#include "Utility_UniformDistribution.hpp"
#include "Utility_HistogramFullyTabularBasicBivariateDistribution.hpp"
#include "Utility_RandomNumberGenerator.hpp"
#include "Utility_UnitTestHarnessWithMain.hpp"
#include "ArchiveTestHelpers.hpp"

//---------------------------------------------------------------------------//
// Testing Types
//---------------------------------------------------------------------------//

using namespace MonteCarlo;



  /* spatial index dimension */
typedef std::tuple<
  std::tuple<std::integral_constant<PhaseSpaceDimension,SPATIAL_INDEX_DIMENSION>,
             std::integral_constant<PhaseSpaceDimension,PRIMARY_DIRECTIONAL_DIMENSION> >,
  std::tuple<std::integral_constant<PhaseSpaceDimension,SPATIAL_INDEX_DIMENSION>,
             std::integral_constant<PhaseSpaceDimension,SECONDARY_DIRECTIONAL_DIMENSION> >,
  std::tuple<std::integral_constant<PhaseSpaceDimension,SPATIAL_INDEX_DIMENSION>,
             std::integral_constant<PhaseSpaceDimension,TERTIARY_DIRECTIONAL_DIMENSION> >,
  std::tuple<std::integral_constant<PhaseSpaceDimension,SPATIAL_INDEX_DIMENSION>,
             std::integral_constant<PhaseSpaceDimension,ENERGY_DIMENSION> >,
  std::tuple<std::integral_constant<PhaseSpaceDimension,SPATIAL_INDEX_DIMENSION>,
             std::integral_constant<PhaseSpaceDimension,TIME_DIMENSION> >,

  /* direction index dimension */
  std::tuple<std::integral_constant<PhaseSpaceDimension,DIRECTION_INDEX_DIMENSION>,
             std::integral_constant<PhaseSpaceDimension,PRIMARY_SPATIAL_DIMENSION> >,
  std::tuple<std::integral_constant<PhaseSpaceDimension,DIRECTION_INDEX_DIMENSION>,
             std::integral_constant<PhaseSpaceDimension,SECONDARY_SPATIAL_DIMENSION> >,
  std::tuple<std::integral_constant<PhaseSpaceDimension,DIRECTION_INDEX_DIMENSION>,
             std::integral_constant<PhaseSpaceDimension,TERTIARY_SPATIAL_DIMENSION> >,
  std::tuple<std::integral_constant<PhaseSpaceDimension,DIRECTION_INDEX_DIMENSION>,
             std::integral_constant<PhaseSpaceDimension,ENERGY_DIMENSION> >,
  std::tuple<std::integral_constant<PhaseSpaceDimension,DIRECTION_INDEX_DIMENSION>,
             std::integral_constant<PhaseSpaceDimension,TIME_DIMENSION> >
 > TestIndexRealPhaseSpaceDimensions;

typedef std::tuple<
  /* spatial index dimension */
  std::tuple<std::integral_constant<PhaseSpaceDimension,PRIMARY_DIRECTIONAL_DIMENSION>,
             std::integral_constant<PhaseSpaceDimension,SPATIAL_INDEX_DIMENSION> >,
  std::tuple<std::integral_constant<PhaseSpaceDimension,SECONDARY_DIRECTIONAL_DIMENSION>,
             std::integral_constant<PhaseSpaceDimension,SPATIAL_INDEX_DIMENSION> >,
  std::tuple<std::integral_constant<PhaseSpaceDimension,TERTIARY_DIRECTIONAL_DIMENSION>,
             std::integral_constant<PhaseSpaceDimension,SPATIAL_INDEX_DIMENSION> >,
  std::tuple<std::integral_constant<PhaseSpaceDimension,ENERGY_DIMENSION>,
             std::integral_constant<PhaseSpaceDimension,SPATIAL_INDEX_DIMENSION> >,
  std::tuple<std::integral_constant<PhaseSpaceDimension,TIME_DIMENSION>,
             std::integral_constant<PhaseSpaceDimension,SPATIAL_INDEX_DIMENSION> >,
  /* direction index dimension */
  std::tuple<std::integral_constant<PhaseSpaceDimension,PRIMARY_SPATIAL_DIMENSION>,
             std::integral_constant<PhaseSpaceDimension,DIRECTION_INDEX_DIMENSION> >,
  std::tuple<std::integral_constant<PhaseSpaceDimension,SECONDARY_SPATIAL_DIMENSION>,
             std::integral_constant<PhaseSpaceDimension,DIRECTION_INDEX_DIMENSION> >,
  std::tuple<std::integral_constant<PhaseSpaceDimension,TERTIARY_SPATIAL_DIMENSION>,
             std::integral_constant<PhaseSpaceDimension,DIRECTION_INDEX_DIMENSION> >,
  std::tuple<std::integral_constant<PhaseSpaceDimension,ENERGY_DIMENSION>,
             std::integral_constant<PhaseSpaceDimension,DIRECTION_INDEX_DIMENSION> >,
  std::tuple<std::integral_constant<PhaseSpaceDimension,TIME_DIMENSION>,
             std::integral_constant<PhaseSpaceDimension,DIRECTION_INDEX_DIMENSION> >
  > TestRealIndexPhaseSpaceDimensions;

typedef std::tuple<
  std::tuple<std::integral_constant<PhaseSpaceDimension,DIRECTION_INDEX_DIMENSION>,
             std::integral_constant<PhaseSpaceDimension,SPATIAL_INDEX_DIMENSION> >,
  std::tuple<std::integral_constant<PhaseSpaceDimension,SPATIAL_INDEX_DIMENSION>,
             std::integral_constant<PhaseSpaceDimension,DIRECTION_INDEX_DIMENSION> >
  > TestIndexIndexPhaseSpaceDimensions;
typedef TestArchiveHelper::TestArchives TestArchives;


typedef decltype(std::tuple_cat(TestRealIndexPhaseSpaceDimensions(), 
                                TestIndexRealPhaseSpaceDimensions(),
                                TestIndexIndexPhaseSpaceDimensions())) TestPhaseSpaceDimensions;
//---------------------------------------------------------------------------//
// Testing Variables
//---------------------------------------------------------------------------//
std::shared_ptr<const Utility::SpatialCoordinateConversionPolicy>
spatial_coord_conversion_policy( new Utility::BasicCartesianCoordinateConversionPolicy );

std::shared_ptr<const Utility::DirectionalCoordinateConversionPolicy>
directional_coord_conversion_policy( new Utility::BasicCartesianCoordinateConversionPolicy );

std::shared_ptr<const Utility::BasicBivariateDistribution> index_real_raw_distribution;
std::shared_ptr<const Utility::BasicBivariateDistribution> real_index_raw_distribution;
std::shared_ptr<const Utility::BasicBivariateDistribution> index_index_raw_distribution;

//---------------------------------------------------------------------------//
// Tests.
//---------------------------------------------------------------------------//
// Test that the dimension can be returned (behavior of underlying distribution doesn't matter)
FRENSIE_UNIT_TEST_TEMPLATE_EXPAND( DependentPhaseSpaceDimensionDistribution_indexDimension,
                                   getDimension,
                                   TestPhaseSpaceDimensions )
{
  FETCH_TEMPLATE_PARAM( 0, WrappedParentDimension );
  FETCH_TEMPLATE_PARAM( 1, WrappedDimension );

  constexpr PhaseSpaceDimension ParentDimension = WrappedParentDimension::value;
  constexpr PhaseSpaceDimension Dimension = WrappedDimension::value;

  std::unique_ptr<const MonteCarlo::PhaseSpaceDimensionDistribution>
    index_real_dimension_distribution( new MonteCarlo::DependentPhaseSpaceDimensionDistribution<ParentDimension,Dimension>( index_real_raw_distribution ) );

  FRENSIE_CHECK_EQUAL(index_real_dimension_distribution->getDimension(), Dimension );

}

//---------------------------------------------------------------------------//
// Test that the dimension class can be returned (behavior of underlying distribution doesn't matter)
FRENSIE_UNIT_TEST_TEMPLATE_EXPAND( DependentPhaseSpaceDimensionDistribution_indexDimension,
                                   getDimensionClass,
                                   TestPhaseSpaceDimensions )
{
  FETCH_TEMPLATE_PARAM( 0, WrappedParentDimension );
  FETCH_TEMPLATE_PARAM( 1, WrappedDimension );

  constexpr PhaseSpaceDimension ParentDimension = WrappedParentDimension::value;
  constexpr PhaseSpaceDimension Dimension = WrappedDimension::value;
  
std::unique_ptr<const MonteCarlo::PhaseSpaceDimensionDistribution>
    index_real_dimension_distribution( new MonteCarlo::DependentPhaseSpaceDimensionDistribution<ParentDimension,Dimension>( index_real_raw_distribution ) );

  FRENSIE_CHECK_EQUAL( index_real_dimension_distribution->getDimensionClass(),
                       MonteCarlo::PhaseSpaceDimensionTraits<Dimension>::getClass() );
}

//---------------------------------------------------------------------------//
// Test that the indep. dimension can be returned (behavior of underlying distribution doesn't matter)
FRENSIE_UNIT_TEST_TEMPLATE_EXPAND( DependentPhaseSpaceDimensionDistribution_indexDimension,
                                   getParentDimension,
                                   TestPhaseSpaceDimensions )
{
  FETCH_TEMPLATE_PARAM( 0, WrappedParentDimension );
  FETCH_TEMPLATE_PARAM( 1, WrappedDimension );

  constexpr PhaseSpaceDimension ParentDimension = WrappedParentDimension::value;
  constexpr PhaseSpaceDimension Dimension = WrappedDimension::value;
  
  std::unique_ptr<const MonteCarlo::DependentPhaseSpaceDimensionDistribution<ParentDimension,Dimension> >
    index_real_dimension_distribution( new MonteCarlo::DependentPhaseSpaceDimensionDistribution<ParentDimension,Dimension>( index_real_raw_distribution ) );

  FRENSIE_CHECK_EQUAL( index_real_dimension_distribution->getParentDimension(),
                       ParentDimension );
}

//---------------------------------------------------------------------------//
// Test that the indep. dimension class can be returned (behavior of underlying distribution doesn't matter)
FRENSIE_UNIT_TEST_TEMPLATE_EXPAND( DependentPhaseSpaceDimensionDistribution_indexDimension,
                                   getParentDimensionClass,
                                   TestPhaseSpaceDimensions )
{
  FETCH_TEMPLATE_PARAM( 0, WrappedParentDimension );
  FETCH_TEMPLATE_PARAM( 1, WrappedDimension );

  constexpr PhaseSpaceDimension ParentDimension = WrappedParentDimension::value;
  constexpr PhaseSpaceDimension Dimension = WrappedDimension::value;
  
  std::unique_ptr<const MonteCarlo::DependentPhaseSpaceDimensionDistribution<ParentDimension,Dimension> >
    index_real_dimension_distribution( new MonteCarlo::DependentPhaseSpaceDimensionDistribution<ParentDimension,Dimension>( index_real_raw_distribution ) );

  FRENSIE_CHECK_EQUAL( index_real_dimension_distribution->getParentDimensionClass(),
                       MonteCarlo::PhaseSpaceDimensionTraits<ParentDimension>::getClass() );
}

//---------------------------------------------------------------------------//
// Test that the distribution type name can be returned (behavior of underlying distribution doesn't matter)
FRENSIE_UNIT_TEST_TEMPLATE_EXPAND( DependentPhaseSpaceDimensionDistribution_indexDimension,
                                   getDistributionTypeName,
                                   TestPhaseSpaceDimensions )
{
  FETCH_TEMPLATE_PARAM( 0, WrappedParentDimension );
  FETCH_TEMPLATE_PARAM( 1, WrappedDimension );

  constexpr PhaseSpaceDimension ParentDimension = WrappedParentDimension::value;
  constexpr PhaseSpaceDimension Dimension = WrappedDimension::value;
  
  std::unique_ptr<const MonteCarlo::PhaseSpaceDimensionDistribution>
    index_real_dimension_distribution( new MonteCarlo::DependentPhaseSpaceDimensionDistribution<ParentDimension,Dimension>( index_real_raw_distribution ) );

  FRENSIE_CHECK_EQUAL( index_real_dimension_distribution->getDistributionTypeName(),
                       "BasicBivariateDistribution" );
}

//---------------------------------------------------------------------------//
// Test that the distribution is not independent (behavior of underlying distribution doesn't matter)
FRENSIE_UNIT_TEST_TEMPLATE_EXPAND( DependentPhaseSpaceDimensionDistribution_indexDimension,
                                   isIndependent,
                                   TestPhaseSpaceDimensions )
{
  FETCH_TEMPLATE_PARAM( 0, WrappedParentDimension );
  FETCH_TEMPLATE_PARAM( 1, WrappedDimension );

  constexpr PhaseSpaceDimension ParentDimension = WrappedParentDimension::value;
  constexpr PhaseSpaceDimension Dimension = WrappedDimension::value;
  
  std::unique_ptr<const MonteCarlo::PhaseSpaceDimensionDistribution>
    index_real_dimension_distribution( new MonteCarlo::DependentPhaseSpaceDimensionDistribution<ParentDimension,Dimension>( index_real_raw_distribution ) );

  FRENSIE_CHECK( !index_real_dimension_distribution->isIndependent() );
}

//---------------------------------------------------------------------------//
// Test if the distribution is dependent on another dimension (behavior of underlying distribution doesn't matter)
FRENSIE_UNIT_TEST_TEMPLATE_EXPAND( DependentPhaseSpaceDimensionDistribution_indexDimension,
                                   isDependentOnDimension,
                                   TestPhaseSpaceDimensions )
{
  FETCH_TEMPLATE_PARAM( 0, WrappedParentDimension );
  FETCH_TEMPLATE_PARAM( 1, WrappedDimension );

  constexpr PhaseSpaceDimension ParentDimension = WrappedParentDimension::value;
  constexpr PhaseSpaceDimension Dimension = WrappedDimension::value;
  
  std::unique_ptr<const MonteCarlo::PhaseSpaceDimensionDistribution>
    index_real_dimension_distribution( new MonteCarlo::DependentPhaseSpaceDimensionDistribution<ParentDimension,Dimension>( index_real_raw_distribution ) );

  FRENSIE_CHECK( index_real_dimension_distribution->isDependentOnDimension( ParentDimension ) );
  FRENSIE_CHECK( !index_real_dimension_distribution->isDependentOnDimension( Dimension ) );
}

//---------------------------------------------------------------------------//
// Test if the distribution is continuous (behavior of underlying distribution doesn't matter)
FRENSIE_UNIT_TEST_TEMPLATE_EXPAND( DependentPhaseSpaceDimensionDistribution_indexDimension,
                                   isContinuous,
                                   TestPhaseSpaceDimensions )
{
  FETCH_TEMPLATE_PARAM( 0, WrappedParentDimension );
  FETCH_TEMPLATE_PARAM( 1, WrappedDimension );

  constexpr PhaseSpaceDimension ParentDimension = WrappedParentDimension::value;
  constexpr PhaseSpaceDimension Dimension = WrappedDimension::value;
  
  std::unique_ptr<const MonteCarlo::PhaseSpaceDimensionDistribution>
    index_real_dimension_distribution( new MonteCarlo::DependentPhaseSpaceDimensionDistribution<ParentDimension,Dimension>( index_real_raw_distribution ) );

  FRENSIE_CHECK( index_real_dimension_distribution->isContinuous() );
}

//---------------------------------------------------------------------------//
// Test if the distribution is tabular
FRENSIE_UNIT_TEST_TEMPLATE_EXPAND( DependentPhaseSpaceDimensionDistribution_indexDimension,
                                   isTabular,
                                   TestPhaseSpaceDimensions )
{
  FETCH_TEMPLATE_PARAM( 0, WrappedParentDimension );
  FETCH_TEMPLATE_PARAM( 1, WrappedDimension );

  constexpr PhaseSpaceDimension ParentDimension = WrappedParentDimension::value;
  constexpr PhaseSpaceDimension Dimension = WrappedDimension::value;
  
  std::unique_ptr<const MonteCarlo::PhaseSpaceDimensionDistribution>
    index_real_dimension_distribution( new MonteCarlo::DependentPhaseSpaceDimensionDistribution<ParentDimension,Dimension>( index_real_raw_distribution ) );

  FRENSIE_CHECK( index_real_dimension_distribution->isTabular() );
}

//---------------------------------------------------------------------------//
// Test if the distribution can be evaluated without a cascade
FRENSIE_UNIT_TEST_TEMPLATE_EXPAND( DependentPhaseSpaceDimensionDistribution_indexDimension,
                                   evaluateWithoutCascade_IndexReal,
                                   TestIndexRealPhaseSpaceDimensions )
{
  FETCH_TEMPLATE_PARAM( 0, WrappedParentDimension );
  FETCH_TEMPLATE_PARAM( 1, WrappedDimension );

  constexpr PhaseSpaceDimension ParentDimension = WrappedParentDimension::value;
  constexpr PhaseSpaceDimension Dimension = WrappedDimension::value;
  
  std::unique_ptr<const MonteCarlo::PhaseSpaceDimensionDistribution>
    index_real_dimension_distribution( new MonteCarlo::DependentPhaseSpaceDimensionDistribution<ParentDimension,Dimension>( index_real_raw_distribution ) );

  MonteCarlo::PhaseSpacePoint point( spatial_coord_conversion_policy,
                                     directional_coord_conversion_policy );

  std::vector<double> independent_values {0, 1, 1.5, 2-1e-13};
  double second_independent_value = 0.6;
  std::vector<double> expected_results {0.5, 0.4, 0.4, 0.4};

  for(size_t ind_value = 0; ind_value < independent_values.size(); ++ind_value)
  {
    setCoordinate<ParentDimension>(point, independent_values[ind_value]);
    setCoordinate<Dimension>(point, second_independent_value);
    FRENSIE_CHECK_EQUAL(index_real_dimension_distribution->evaluateWithoutCascade( point ), expected_results[ind_value]);
  }
  
  // Should throw exception when index is out of bounds to counter roundoff errors resulting in undefined behavior
  FRENSIE_CHECK_THROW(setCoordinate<ParentDimension>( point, 72 ), std::runtime_error);

}

//---------------------------------------------------------------------------//
// Test if the distribution can be evaluated without a cascade
FRENSIE_UNIT_TEST_TEMPLATE_EXPAND( DependentPhaseSpaceDimensionDistribution_indexDimension,
                                   evaluateWithoutCascade_RealIndex,
                                   TestRealIndexPhaseSpaceDimensions )
{
  FETCH_TEMPLATE_PARAM( 0, WrappedParentDimension );
  FETCH_TEMPLATE_PARAM( 1, WrappedDimension );

  constexpr PhaseSpaceDimension ParentDimension = WrappedParentDimension::value;
  constexpr PhaseSpaceDimension Dimension = WrappedDimension::value;
  
  std::unique_ptr<const MonteCarlo::PhaseSpaceDimensionDistribution>
    real_index_dimension_distribution( new MonteCarlo::DependentPhaseSpaceDimensionDistribution<ParentDimension,Dimension>( real_index_raw_distribution ) );

  MonteCarlo::PhaseSpacePoint point( spatial_coord_conversion_policy,
                                     directional_coord_conversion_policy );

  std::vector<double> independent_values {0.4, 0.6};
  std::vector<double> second_independent_values {0, 0.5, 1, 2-1e-14};
  std::vector<std::vector<double>> expected_results {{2.0, 2.0, 3.0, 3.0}, {1.0, 1.0, 5.0, 5.0}};

  for(size_t ind_value = 0; ind_value < independent_values.size(); ++ind_value)
  {
    setCoordinate<ParentDimension>(point, independent_values[ind_value]);
    for(size_t sec_ind_value = 0; sec_ind_value < second_independent_values.size(); ++sec_ind_value)
    {
      setCoordinate<Dimension>(point, second_independent_values[sec_ind_value]);
      FRENSIE_CHECK_EQUAL(real_index_dimension_distribution->evaluateWithoutCascade( point ), expected_results[ind_value][sec_ind_value]);
    }
    
    
  }
  // Should throw exception when index is out of bounds to counter roundoff errors resulting in undefined behavior
  FRENSIE_CHECK_THROW(setCoordinate<Dimension>( point, 72 ), std::runtime_error);

}

//---------------------------------------------------------------------------//
// Test if the distribution can be evaluated without a cascade
FRENSIE_UNIT_TEST_TEMPLATE_EXPAND( DependentPhaseSpaceDimensionDistribution_indexDimension,
                                   evaluateWithoutCascade_IndexIndex,
                                   TestIndexIndexPhaseSpaceDimensions )
{
  FETCH_TEMPLATE_PARAM( 0, WrappedParentDimension );
  FETCH_TEMPLATE_PARAM( 1, WrappedDimension );

  constexpr PhaseSpaceDimension ParentDimension = WrappedParentDimension::value;
  constexpr PhaseSpaceDimension Dimension = WrappedDimension::value;
  
  std::unique_ptr<const MonteCarlo::PhaseSpaceDimensionDistribution>
    index_index_dimension_distribution( new MonteCarlo::DependentPhaseSpaceDimensionDistribution<ParentDimension,Dimension>( index_index_raw_distribution ) );

  MonteCarlo::PhaseSpacePoint point( spatial_coord_conversion_policy,
                                     directional_coord_conversion_policy );

  std::vector<double> independent_values {0, 0.5, 1, 2-1e-14};
  std::vector<double> second_independent_values {0, 0.5, 1, 2-1e-14};
  std::vector<std::vector<double>> expected_results {{2.0, 2.0, 3.0, 3.0}, {2.0, 2.0, 3.0, 3.0}, {1.0, 1.0, 5.0, 5.0}, {1.0, 1.0, 5.0, 5.0}};

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

FRENSIE_UNIT_TEST_TEMPLATE_EXPAND( DependentPhaseSpaceDimensionDistribution_indexDimension,
                                   sampleWithoutCascade_IndexReal,
                                   TestIndexRealPhaseSpaceDimensions )
{
  FETCH_TEMPLATE_PARAM( 0, WrappedParentDimension );
  FETCH_TEMPLATE_PARAM( 1, WrappedDimension );

  constexpr PhaseSpaceDimension ParentDimension = WrappedParentDimension::value;
  constexpr PhaseSpaceDimension Dimension = WrappedDimension::value;

  std::unique_ptr<const MonteCarlo::PhaseSpaceDimensionDistribution>
    index_real_dimension_distribution( new MonteCarlo::DependentPhaseSpaceDimensionDistribution<ParentDimension,Dimension>( index_real_raw_distribution ) );

  MonteCarlo::PhaseSpacePoint point( spatial_coord_conversion_policy,
                                     directional_coord_conversion_policy );

  std::vector<double> index_values = {0, 0.5, 1, 1.5};
  std::vector<double> initial_fake_stream = {0.0, 0.5, 1-1e-15};
  std::vector<double> index_additional_fake_stream = {0.5, 0.5, 0.5};
  std::vector<std::vector<double>> results = {{0.5, 0.7, 0.8999999999999996}, {0.5, 0.7, 0.8999999999999996}, {0.6, 0.7, 0.7999999999999998}, {0.6, 0.7, 0.7999999999999998}};

  for(size_t index = 0; index < index_values.size(); ++index)
  {
    setCoordinate<ParentDimension>( point, index_values[index] );
    for(size_t sample = 0; sample < initial_fake_stream.size(); ++sample)
    {
      std::vector<double> fake_stream(4);
      fake_stream[0] = initial_fake_stream[sample];
      fake_stream.insert( fake_stream.end(), index_additional_fake_stream.begin(), index_additional_fake_stream.end() );
      Utility::RandomNumberGenerator::setFakeStream( fake_stream );
      index_real_dimension_distribution->sampleWithoutCascade( point );
      FRENSIE_CHECK_FLOATING_EQUALITY(getCoordinate<Dimension>( point ), results[index][sample], 1e-15);
    }
    Utility::RandomNumberGenerator::unsetFakeStream();
  }
}

FRENSIE_UNIT_TEST_TEMPLATE_EXPAND( DependentPhaseSpaceDimensionDistribution_indexDimension,
                                   sampleWithoutCascade_RealIndex,
                                   TestRealIndexPhaseSpaceDimensions )
{
  FETCH_TEMPLATE_PARAM( 0, WrappedParentDimension );
  FETCH_TEMPLATE_PARAM( 1, WrappedDimension );

  constexpr PhaseSpaceDimension ParentDimension = WrappedParentDimension::value;
  constexpr PhaseSpaceDimension Dimension = WrappedDimension::value;

  std::unique_ptr<const MonteCarlo::PhaseSpaceDimensionDistribution>
    real_index_dimension_distribution( new MonteCarlo::DependentPhaseSpaceDimensionDistribution<ParentDimension,Dimension>( real_index_raw_distribution ) );

  MonteCarlo::PhaseSpacePoint point( spatial_coord_conversion_policy,
                                     directional_coord_conversion_policy );

  std::vector<double> real_values = {0.4, 0.6};
  std::vector<double> initial_fake_stream = {0.0, 0.5, 1-1e-15};
  std::vector<double> index_additional_fake_stream = {0.5, 0.5, 0.5};
  std::vector<std::vector<size_t>> results = {{0, 1, 1}, {0, 1, 1}};

  for(size_t index = 0; index < real_values.size(); ++index)
  {
    setCoordinate<ParentDimension>( point, real_values[index] );
    for(size_t sample = 0; sample < initial_fake_stream.size(); ++sample)
    {
      std::vector<double> fake_stream(4);
      fake_stream[0] = initial_fake_stream[sample];
      fake_stream.insert( fake_stream.end(), index_additional_fake_stream.begin(), index_additional_fake_stream.end() );
      Utility::RandomNumberGenerator::setFakeStream( fake_stream );
      real_index_dimension_distribution->sampleWithoutCascade( point );
      FRENSIE_CHECK_EQUAL(getCoordinate<Dimension>( point ), results[index][sample]);
    }
    Utility::RandomNumberGenerator::unsetFakeStream();
  }
}

FRENSIE_UNIT_TEST_TEMPLATE_EXPAND( DependentPhaseSpaceDimensionDistribution_indexDimension,
                                   sampleWithoutCascade_IndexIndex,
                                   TestIndexIndexPhaseSpaceDimensions )
{
  FETCH_TEMPLATE_PARAM( 0, WrappedParentDimension );
  FETCH_TEMPLATE_PARAM( 1, WrappedDimension );

  constexpr PhaseSpaceDimension ParentDimension = WrappedParentDimension::value;
  constexpr PhaseSpaceDimension Dimension = WrappedDimension::value;

  std::unique_ptr<const MonteCarlo::PhaseSpaceDimensionDistribution>
    index_index_dimension_distribution( new MonteCarlo::DependentPhaseSpaceDimensionDistribution<ParentDimension,Dimension>( index_index_raw_distribution ) );

  MonteCarlo::PhaseSpacePoint point( spatial_coord_conversion_policy,
                                     directional_coord_conversion_policy );

  std::vector<double> index_values = {0.5, 1.0};
  std::vector<double> initial_fake_stream = {0.0, 0.5, 1-1e-15};
  std::vector<double> index_additional_fake_stream = {0.5, 0.5, 0.5};
  std::vector<std::vector<size_t>> results = {{0, 1, 1}, {0, 1, 1}};

  for(size_t index = 0; index < index_values.size(); ++index)
  {
    setCoordinate<ParentDimension>( point, index_values[index] );
    for(size_t sample = 0; sample < initial_fake_stream.size(); ++sample)
    {
      std::vector<double> fake_stream(4);
      fake_stream[0] = initial_fake_stream[sample];
      fake_stream.insert( fake_stream.end(), index_additional_fake_stream.begin(), index_additional_fake_stream.end() );
      Utility::RandomNumberGenerator::setFakeStream( fake_stream );
      index_index_dimension_distribution->sampleWithoutCascade( point );
      FRENSIE_CHECK_EQUAL(getCoordinate<Dimension>( point ), results[index][sample]);
    }
    Utility::RandomNumberGenerator::unsetFakeStream();
  }
}

//---------------------------------------------------------------------------//
// Test if the distribution can be sampled without a cascade and the
// trials can be counted
FRENSIE_UNIT_TEST_TEMPLATE_EXPAND( DependentPhaseSpaceDimensionDistribution_indexDimension,
                                   sampleAndRecordTrialsWithoutCascade_IndexReal,
                                   TestIndexRealPhaseSpaceDimensions )
{

  FETCH_TEMPLATE_PARAM( 0, WrappedParentDimension );
  FETCH_TEMPLATE_PARAM( 1, WrappedDimension );

  constexpr PhaseSpaceDimension ParentDimension = WrappedParentDimension::value;
  constexpr PhaseSpaceDimension Dimension = WrappedDimension::value;

  std::unique_ptr<const MonteCarlo::PhaseSpaceDimensionDistribution>
    index_real_dimension_distribution( new MonteCarlo::DependentPhaseSpaceDimensionDistribution<ParentDimension,Dimension>( index_real_raw_distribution ) );

  MonteCarlo::PhaseSpacePoint point( spatial_coord_conversion_policy,
                                     directional_coord_conversion_policy );

  // First random number is used to sample the index. The ones after that are for determining where in that element

  MonteCarlo::PhaseSpaceDimensionDistribution::Counter trials = 0;

  std::vector<double> index_values = {0, 0.5, 1, 1.5};
  std::vector<double> initial_fake_stream = {0.0, 0.5, 1-1e-15};
  std::vector<double> index_additional_fake_stream = {0.5, 0.5, 0.5};
  std::vector<std::vector<double>> results = {{0.5, 0.7, 0.8999999999999996}, {0.5, 0.7, 0.8999999999999996}, {0.6, 0.7, 0.7999999999999998}, {0.6, 0.7, 0.7999999999999998}};

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
      index_real_dimension_distribution->sampleAndRecordTrialsWithoutCascade( point, trials );
      FRENSIE_CHECK_FLOATING_EQUALITY(getCoordinate<Dimension>( point ), results[index][sample], 1e-15);
      counter += 1;
      FRENSIE_CHECK_EQUAL(trials, counter);
    }
    Utility::RandomNumberGenerator::unsetFakeStream();
  }
}

//---------------------------------------------------------------------------//
// Test if the distribution can be sampled without a cascade and the
// trials can be counted
FRENSIE_UNIT_TEST_TEMPLATE_EXPAND( DependentPhaseSpaceDimensionDistribution_indexDimension,
                                   sampleAndRecordTrialsWithoutCascade_RealIndex,
                                   TestRealIndexPhaseSpaceDimensions )
{

  FETCH_TEMPLATE_PARAM( 0, WrappedParentDimension );
  FETCH_TEMPLATE_PARAM( 1, WrappedDimension );

  constexpr PhaseSpaceDimension ParentDimension = WrappedParentDimension::value;
  constexpr PhaseSpaceDimension Dimension = WrappedDimension::value;

  std::unique_ptr<const MonteCarlo::PhaseSpaceDimensionDistribution>
    real_index_dimension_distribution( new MonteCarlo::DependentPhaseSpaceDimensionDistribution<ParentDimension,Dimension>( real_index_raw_distribution ) );

  MonteCarlo::PhaseSpacePoint point( spatial_coord_conversion_policy,
                                     directional_coord_conversion_policy );

  // First random number is used to sample the index. The ones after that are for determining where in that element

  MonteCarlo::PhaseSpaceDimensionDistribution::Counter trials = 0;

  std::vector<double> real_values = {0.4, 0.6};
  std::vector<double> initial_fake_stream = {0.0, 0.5, 1-1e-15};
  std::vector<double> index_additional_fake_stream = {0.5, 0.5, 0.5};
  std::vector<std::vector<size_t>> results = {{0, 1, 1}, {0, 1, 1}};

  size_t counter = 0;
  for(size_t index = 0; index < real_values.size(); ++index)
  {
    setCoordinate<ParentDimension>( point, real_values[index] );
    for(size_t sample = 0; sample < initial_fake_stream.size(); ++sample)
    {
      std::vector<double> fake_stream(4);
      fake_stream[0] = initial_fake_stream[sample];
      fake_stream.insert( fake_stream.end(), index_additional_fake_stream.begin(), index_additional_fake_stream.end() );
      Utility::RandomNumberGenerator::setFakeStream( fake_stream );
      real_index_dimension_distribution->sampleAndRecordTrialsWithoutCascade( point, trials );
      FRENSIE_CHECK_EQUAL(getCoordinate<Dimension>( point ), results[index][sample]);
      counter += 1;
      FRENSIE_CHECK_EQUAL(trials, counter);
    }
    Utility::RandomNumberGenerator::unsetFakeStream();
  }
}

//---------------------------------------------------------------------------//
// Test if the distribution can be sampled without a cascade and the
// trials can be counted
FRENSIE_UNIT_TEST_TEMPLATE_EXPAND( DependentPhaseSpaceDimensionDistribution_indexDimension,
                                   sampleAndRecordTrialsWithoutCascade_IndexIndex,
                                   TestIndexIndexPhaseSpaceDimensions )
{

  FETCH_TEMPLATE_PARAM( 0, WrappedParentDimension );
  FETCH_TEMPLATE_PARAM( 1, WrappedDimension );

  constexpr PhaseSpaceDimension ParentDimension = WrappedParentDimension::value;
  constexpr PhaseSpaceDimension Dimension = WrappedDimension::value;

  std::unique_ptr<const MonteCarlo::PhaseSpaceDimensionDistribution>
    index_index_dimension_distribution( new MonteCarlo::DependentPhaseSpaceDimensionDistribution<ParentDimension,Dimension>( index_index_raw_distribution ) );

  MonteCarlo::PhaseSpacePoint point( spatial_coord_conversion_policy,
                                     directional_coord_conversion_policy );

  // First random number is used to sample the index. The ones after that are for determining where in that element

  MonteCarlo::PhaseSpaceDimensionDistribution::Counter trials = 0;

  std::vector<double> index_values = {0.5, 1.0};
  std::vector<double> initial_fake_stream = {0.0, 0.5, 1-1e-15};
  std::vector<double> index_additional_fake_stream = {0.5, 0.5, 0.5};
  std::vector<std::vector<size_t>> results = {{0, 1, 1}, {0, 1, 1}};

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
      index_index_dimension_distribution->sampleAndRecordTrialsWithoutCascade( point, trials );
      FRENSIE_CHECK_EQUAL(getCoordinate<Dimension>( point ), results[index][sample]);
      counter += 1;
      FRENSIE_CHECK_EQUAL(trials, counter);
    }
    Utility::RandomNumberGenerator::unsetFakeStream();
  }
}

FRENSIE_UNIT_TEST_TEMPLATE_EXPAND( DependentPhaseSpaceDimensionDistribution_indexDimension,
                                   setDimensionValueAndApplyWeight_IndexIndex,
                                   TestIndexIndexPhaseSpaceDimensions )
{
  FETCH_TEMPLATE_PARAM( 0, WrappedParentDimension );
  FETCH_TEMPLATE_PARAM( 1, WrappedDimension );

  constexpr PhaseSpaceDimension ParentDimension = WrappedParentDimension::value;
  constexpr PhaseSpaceDimension Dimension = WrappedDimension::value;
  std::unique_ptr<const MonteCarlo::PhaseSpaceDimensionDistribution>
    index_index_dimension_distribution( new MonteCarlo::DependentPhaseSpaceDimensionDistribution<ParentDimension,Dimension>( index_index_raw_distribution ) );

  MonteCarlo::PhaseSpacePoint point( spatial_coord_conversion_policy,
                                     directional_coord_conversion_policy );

  setCoordinate<ParentDimension>( point, 0 );

  index_index_dimension_distribution->setDimensionValueAndApplyWeight( point, 0 );

  FRENSIE_CHECK_EQUAL( getCoordinate<Dimension>( point ), 0);
  FRENSIE_CHECK_FLOATING_EQUALITY(getCoordinateWeight<Dimension>( point ), 2.0/5.0, 1e-14);

  index_index_dimension_distribution->setDimensionValueAndApplyWeight( point, 1 );

  FRENSIE_CHECK_EQUAL( getCoordinate<Dimension>( point ), 1);
  FRENSIE_CHECK_FLOATING_EQUALITY(getCoordinateWeight<Dimension>( point ), 3.0/5.0, 1e-14);

  setCoordinate<ParentDimension>( point, 1 );

  index_index_dimension_distribution->setDimensionValueAndApplyWeight( point, 0 );

  FRENSIE_CHECK_EQUAL( getCoordinate<Dimension>( point ), 0);
  FRENSIE_CHECK_FLOATING_EQUALITY(getCoordinateWeight<Dimension>( point ), 1.0/6.0, 1e-14);

  index_index_dimension_distribution->setDimensionValueAndApplyWeight( point, 1 );

  FRENSIE_CHECK_EQUAL( getCoordinate<Dimension>( point ), 1);
  FRENSIE_CHECK_FLOATING_EQUALITY(getCoordinateWeight<Dimension>( point ), 5.0/6.0, 1e-14);
}

//---------------------------------------------------------------------------//
// Check that the distribution can be archived
FRENSIE_UNIT_TEST_TEMPLATE_EXPAND( DependentPhaseSpaceDimensionDistribution_indexDimension,
                                   archive_IndexIndex,
                                   TestArchives )
{
  FETCH_TEMPLATE_PARAM( 0, RawOArchive );
  FETCH_TEMPLATE_PARAM( 1, RawIArchive );

  typedef typename std::remove_pointer<RawOArchive>::type OArchive;
  typedef typename std::remove_pointer<RawIArchive>::type IArchive;

  std::string archive_base_name( "test_dependent_phase_dimension_distribution" );
  std::ostringstream archive_ostream;

  {
    std::unique_ptr<OArchive> oarchive;

    createOArchive( archive_base_name, archive_ostream, oarchive );

    // Create the fully tabular distribution
    std::shared_ptr<const Utility::BasicBivariateDistribution>
      raw_distribution;
    
    {
      std::vector<double> primary_grid( {0, 1, 2} );
  
      // Index secondary distributions
      std::vector<std::shared_ptr<const Utility::TabularUnivariateDistribution> >
        index_secondary_dists( 3 );

      std::vector<double> bin_boundaries {0.0, 1.0, 2.0};
      std::vector<double> bin_values {2.0, 3.0};

      index_secondary_dists[0] = std::make_shared< Utility::HistogramDistribution >( bin_boundaries, bin_values);

      // Only bin_values needs to change
      bin_values[0] = 1.0;
      bin_values[1] = 5.0;
        
      index_secondary_dists[1] = std::make_shared< Utility::HistogramDistribution >( bin_boundaries, bin_values);

      index_secondary_dists[2] = index_secondary_dists[0];

      Utility::HistogramFullyTabularBasicBivariateDistribution*
        local_raw_distribution = new Utility::HistogramFullyTabularBasicBivariateDistribution( primary_grid, index_secondary_dists );

      local_raw_distribution->limitToPrimaryIndepLimits();

      raw_distribution.reset( local_raw_distribution );
    }
    // Create the dependent time distributions
    std::shared_ptr<const MonteCarlo::PhaseSpaceDimensionDistribution>
      direction_spatial_index_dimension_distribution( new MonteCarlo::DirectionIndexDependentSpatialIndexDimensionDistribution( raw_distribution ) );
    
    std::shared_ptr<const MonteCarlo::PhaseSpaceDimensionDistribution>
      spatial_direction_index_dimension_distribution( new MonteCarlo::SpatialIndexDependentDirectionIndexDimensionDistribution( raw_distribution ) );

    FRENSIE_REQUIRE_NO_THROW( (*oarchive) << BOOST_SERIALIZATION_NVP(direction_spatial_index_dimension_distribution) );
    FRENSIE_REQUIRE_NO_THROW( (*oarchive) << BOOST_SERIALIZATION_NVP(spatial_direction_index_dimension_distribution) );
  }

  // Copy the archive ostream to an istream
  std::istringstream archive_istream( archive_ostream.str() );

  // Load the archived distributions
  std::unique_ptr<IArchive> iarchive;

  createIArchive( archive_istream, iarchive );

  std::shared_ptr<const MonteCarlo::PhaseSpaceDimensionDistribution>
    direction_spatial_index_dimension_distribution,
    spatial_direction_index_dimension_distribution;

  FRENSIE_REQUIRE_NO_THROW( (*iarchive) >> BOOST_SERIALIZATION_NVP(direction_spatial_index_dimension_distribution) );
  FRENSIE_REQUIRE_NO_THROW( (*iarchive) >> BOOST_SERIALIZATION_NVP(spatial_direction_index_dimension_distribution) );

    iarchive.reset();

    {
      FRENSIE_CHECK_EQUAL( direction_spatial_index_dimension_distribution->getDimension(),
                           SPATIAL_INDEX_DIMENSION );

      auto concrete_direction_spatial_index_dimension_distribution = 
        dynamic_cast<const MonteCarlo::DirectionIndexDependentSpatialIndexDimensionDistribution*>( direction_spatial_index_dimension_distribution.get() );

      FRENSIE_CHECK_EQUAL( concrete_direction_spatial_index_dimension_distribution->getParentDimension(),
                           DIRECTION_INDEX_DIMENSION );

      MonteCarlo::PhaseSpacePoint point( spatial_coord_conversion_policy,
                                         directional_coord_conversion_policy );

      // Parent dimension value outside of distribution bounds
      setCoordinate<DIRECTION_INDEX_DIMENSION>( point, 0 );
      setCoordinate<SPATIAL_INDEX_DIMENSION>( point, 0 );

      FRENSIE_CHECK_EQUAL( direction_spatial_index_dimension_distribution->evaluateWithoutCascade( point ),
                           2 );

      setCoordinate<SPATIAL_INDEX_DIMENSION>( point, 1 );
      
      FRENSIE_CHECK_EQUAL( direction_spatial_index_dimension_distribution->evaluateWithoutCascade( point ),
                           3 );
    }

    {
      FRENSIE_CHECK_EQUAL( spatial_direction_index_dimension_distribution->getDimension(),
                           DIRECTION_INDEX_DIMENSION );

      auto concrete_spatial_direction_index_dimension_distribution = 
        dynamic_cast<const MonteCarlo::SpatialIndexDependentDirectionIndexDimensionDistribution*>( spatial_direction_index_dimension_distribution.get() );

      FRENSIE_CHECK_EQUAL( concrete_spatial_direction_index_dimension_distribution->getParentDimension(),
                           SPATIAL_INDEX_DIMENSION );

      MonteCarlo::PhaseSpacePoint point( spatial_coord_conversion_policy,
                                         directional_coord_conversion_policy );

      // Parent dimension value outside of distribution bounds
      setCoordinate<SPATIAL_INDEX_DIMENSION>( point, 0 );
      setCoordinate<DIRECTION_INDEX_DIMENSION>( point, 0 );

      FRENSIE_CHECK_EQUAL( spatial_direction_index_dimension_distribution->evaluateWithoutCascade( point ),
                           2 );

      setCoordinate<DIRECTION_INDEX_DIMENSION>( point, 1 );
      
      FRENSIE_CHECK_EQUAL( spatial_direction_index_dimension_distribution->evaluateWithoutCascade( point ),
                           3 );
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
                      y_planes( {0.0, 0.5, 1.0} ),
                      z_planes( {0.0, 0.5, 1.0} );

  std::shared_ptr<Utility::StructuredHexMesh> source_mesh = std::make_shared<Utility::StructuredHexMesh>(x_planes, y_planes, z_planes);

  MonteCarlo::PhaseSpaceDimensionTraits<MonteCarlo::SPATIAL_INDEX_DIMENSION>::setMesh(source_mesh);

  std::shared_ptr<Utility::PQLAQuadrature> source_direction_discretization = std::make_shared<Utility::PQLAQuadrature>(2);

  MonteCarlo::PhaseSpaceDimensionTraits<MonteCarlo::DIRECTION_INDEX_DIMENSION>::setDirectionDiscretization(source_direction_discretization);
  
  // 3 different distributions - 1 for index-real dependency, one for real-index, and one for index-index
  // index-real

  // Indexed primary grid
  std::vector<double> index_primary_grid( {0, 1, 2} );
  
  std::vector<std::shared_ptr<const Utility::TabularUnivariateDistribution> >
    real_secondary_dists( 3 );

  // Real secondary distributions
  // Create the secondary distribution in the first bin
  real_secondary_dists[0].reset( new Utility::UniformDistribution( 0.5, 0.9, 0.5 ) );
  
  // Create the secondary distribution in the second bin
  real_secondary_dists[1].reset( new Utility::UniformDistribution( 0.6, 0.8, 0.4 ) );
  
  // Create the secondary distribution in the third bin
  real_secondary_dists[2] = real_secondary_dists[0];

  std::shared_ptr<Utility::HistogramFullyTabularBasicBivariateDistribution> index_real_local_raw_distribution = 
    std::make_shared<Utility::HistogramFullyTabularBasicBivariateDistribution> (index_primary_grid, real_secondary_dists);

  index_real_local_raw_distribution->limitToPrimaryIndepLimits();

  index_real_raw_distribution = index_real_local_raw_distribution;

  // real-index
  
  // Real primary grid
  std::vector<double> real_primary_grid( {0.1, 0.5, 0.9} );
  
  // Index secondary distributions
  std::vector<std::shared_ptr<const Utility::TabularUnivariateDistribution> >
    index_secondary_dists( 3 );

  std::vector<double> bin_boundaries {0.0, 1.0, 2.0};
  std::vector<double> bin_values {2.0, 3.0};

  index_secondary_dists[0] = std::make_shared< Utility::HistogramDistribution >( bin_boundaries, bin_values);

  // Only bin_values needs to change
  bin_values[0] = 1.0;
  bin_values[1] = 5.0;
    
  index_secondary_dists[1] = std::make_shared< Utility::HistogramDistribution >( bin_boundaries, bin_values);

  index_secondary_dists[2] = index_secondary_dists[0];

  std::shared_ptr<Utility::HistogramFullyTabularBasicBivariateDistribution> real_index_local_raw_distribution = 
    std::make_shared<Utility::HistogramFullyTabularBasicBivariateDistribution> (real_primary_grid, index_secondary_dists);

  real_index_local_raw_distribution->limitToPrimaryIndepLimits();

  real_index_raw_distribution = real_index_local_raw_distribution;

  // index-index
  std::shared_ptr<Utility::HistogramFullyTabularBasicBivariateDistribution> index_index_local_raw_distribution = 
    std::make_shared<Utility::HistogramFullyTabularBasicBivariateDistribution> (index_primary_grid, index_secondary_dists);

  index_index_local_raw_distribution->limitToPrimaryIndepLimits();

  index_index_raw_distribution = index_index_local_raw_distribution;

  // Initialize the random number generator
  Utility::RandomNumberGenerator::createStreams();
}

FRENSIE_CUSTOM_UNIT_TEST_SETUP_END();

//---------------------------------------------------------------------------//
// end tstDependentPhaseSpaceDimensionDistribution.cpp
//---------------------------------------------------------------------------//
