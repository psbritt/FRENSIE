//---------------------------------------------------------------------------//
//!
//! \file   tstDependentPhaseSpaceDimensionDistribution_indexDimensions.cpp
//! \author Philip Britt
//! \brief  Dependent phase space dimension dist. unit tests (dimensions w/ integer values)
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
#include "Utility_UniformDistribution.hpp"
#include "Utility_HistogramFullyTabularBasicBivariateDistribution.hpp"
#include "Utility_RandomNumberGenerator.hpp"
#include "Utility_UnitTestHarnessWithMain.hpp"
#include "ArchiveTestHelpers.hpp"
#include "Utility_DiscreteDistribution.hpp"
//---------------------------------------------------------------------------//
// Testing Types
//---------------------------------------------------------------------------//

using namespace MonteCarlo;

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
  std::tuple<std::integral_constant<PhaseSpaceDimension,DIRECTION_INDEX_DIMENSION>,
             std::integral_constant<PhaseSpaceDimension,SPATIAL_INDEX_DIMENSION> >,
  /* spatial index parent dimension */
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
  std::tuple<std::integral_constant<PhaseSpaceDimension,SPATIAL_INDEX_DIMENSION>,
             std::integral_constant<PhaseSpaceDimension,DIRECTION_INDEX_DIMENSION> >,
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
             std::integral_constant<PhaseSpaceDimension,DIRECTION_INDEX_DIMENSION> >,
  /* direction index parent dimension */
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
  /* secondary directional dimension */
 > TestPhaseSpaceDimensions;

typedef TestArchiveHelper::TestArchives TestArchives;

//---------------------------------------------------------------------------//
// Testing Variables
//---------------------------------------------------------------------------//
std::shared_ptr<const Utility::SpatialCoordinateConversionPolicy>
spatial_coord_conversion_policy( new Utility::BasicCartesianCoordinateConversionPolicy );

std::shared_ptr<const Utility::DirectionalCoordinateConversionPolicy>
directional_coord_conversion_policy( new Utility::BasicCartesianCoordinateConversionPolicy );

std::shared_ptr<const Utility::BasicBivariateDistribution> spatial_parent_raw_distribution;
std::shared_ptr<const Utility::BasicBivariateDistribution> direction_parent_raw_distribution;
std::shared_ptr<const Utility::BasicBivariateDistribution> spatial_direction_raw_distribution;
std::shared_ptr<const Utility::BasicBivariateDistribution> direction_spatial_raw_distribution;

//---------------------------------------------------------------------------//
// Tests.
//---------------------------------------------------------------------------//
// Test that the dimension can be returned
FRENSIE_UNIT_TEST_TEMPLATE_EXPAND( DependentPhaseSpaceDimensionDistribution,
                                   getDimension,
                                   TestPhaseSpaceDimensions )
{
  FETCH_TEMPLATE_PARAM( 0, WrappedParentDimension );
  FETCH_TEMPLATE_PARAM( 1, WrappedDimension );

  constexpr PhaseSpaceDimension ParentDimension = WrappedParentDimension::value;
  constexpr PhaseSpaceDimension Dimension = WrappedDimension::value;

}

//---------------------------------------------------------------------------//
// Test that the dimension class can be returned
FRENSIE_UNIT_TEST_TEMPLATE_EXPAND( DependentPhaseSpaceDimensionDistribution,
                                   getDimensionClass,
                                   TestPhaseSpaceDimensions )
{
  FETCH_TEMPLATE_PARAM( 0, WrappedParentDimension );
  FETCH_TEMPLATE_PARAM( 1, WrappedDimension );

  constexpr PhaseSpaceDimension ParentDimension = WrappedParentDimension::value;
  constexpr PhaseSpaceDimension Dimension = WrappedDimension::value;

}

//---------------------------------------------------------------------------//
// Test that the indep. dimension can be returned
FRENSIE_UNIT_TEST_TEMPLATE_EXPAND( DependentPhaseSpaceDimensionDistribution,
                                   getParentDimension,
                                   TestPhaseSpaceDimensions )
{
  FETCH_TEMPLATE_PARAM( 0, WrappedParentDimension );
  FETCH_TEMPLATE_PARAM( 1, WrappedDimension );

  constexpr PhaseSpaceDimension ParentDimension = WrappedParentDimension::value;
  constexpr PhaseSpaceDimension Dimension = WrappedDimension::value;

}

//---------------------------------------------------------------------------//
// Test that the indep. dimension class can be returned
FRENSIE_UNIT_TEST_TEMPLATE_EXPAND( DependentPhaseSpaceDimensionDistribution,
                                   getParentDimensionClass,
                                   TestPhaseSpaceDimensions )
{
  FETCH_TEMPLATE_PARAM( 0, WrappedParentDimension );
  FETCH_TEMPLATE_PARAM( 1, WrappedDimension );

  constexpr PhaseSpaceDimension ParentDimension = WrappedParentDimension::value;
  constexpr PhaseSpaceDimension Dimension = WrappedDimension::value;

}

//---------------------------------------------------------------------------//
// Test that the distribution type name can be returned
FRENSIE_UNIT_TEST_TEMPLATE_EXPAND( DependentPhaseSpaceDimensionDistribution,
                                   getDistributionTypeName,
                                   TestPhaseSpaceDimensions )
{
  FETCH_TEMPLATE_PARAM( 0, WrappedParentDimension );
  FETCH_TEMPLATE_PARAM( 1, WrappedDimension );

  constexpr PhaseSpaceDimension ParentDimension = WrappedParentDimension::value;
  constexpr PhaseSpaceDimension Dimension = WrappedDimension::value;

}

//---------------------------------------------------------------------------//
// Test that the distribution type name can be returned
FRENSIE_UNIT_TEST_TEMPLATE_EXPAND( DependentPhaseSpaceDimensionDistribution,
                                   isIndependent,
                                   TestPhaseSpaceDimensions )
{
  FETCH_TEMPLATE_PARAM( 0, WrappedParentDimension );
  FETCH_TEMPLATE_PARAM( 1, WrappedDimension );

  constexpr PhaseSpaceDimension ParentDimension = WrappedParentDimension::value;
  constexpr PhaseSpaceDimension Dimension = WrappedDimension::value;

}

//---------------------------------------------------------------------------//
// Test if the distribution is dependent on another dimension
FRENSIE_UNIT_TEST_TEMPLATE_EXPAND( DependentPhaseSpaceDimensionDistribution,
                                   isDependentOnDimension,
                                   TestPhaseSpaceDimensions )
{
  FETCH_TEMPLATE_PARAM( 0, WrappedParentDimension );
  FETCH_TEMPLATE_PARAM( 1, WrappedDimension );

  constexpr PhaseSpaceDimension ParentDimension = WrappedParentDimension::value;
  constexpr PhaseSpaceDimension Dimension = WrappedDimension::value;

}

//---------------------------------------------------------------------------//
// Test if the distribution is continuous
FRENSIE_UNIT_TEST_TEMPLATE_EXPAND( DependentPhaseSpaceDimensionDistribution,
                                   isContinuous,
                                   TestPhaseSpaceDimensions )
{
  FETCH_TEMPLATE_PARAM( 0, WrappedParentDimension );
  FETCH_TEMPLATE_PARAM( 1, WrappedDimension );

  constexpr PhaseSpaceDimension ParentDimension = WrappedParentDimension::value;
  constexpr PhaseSpaceDimension Dimension = WrappedDimension::value;

}

//---------------------------------------------------------------------------//
// Test if the distribution is tabular
FRENSIE_UNIT_TEST_TEMPLATE_EXPAND( DependentPhaseSpaceDimensionDistribution,
                                   isTabular,
                                   TestPhaseSpaceDimensions )
{
  FETCH_TEMPLATE_PARAM( 0, WrappedParentDimension );
  FETCH_TEMPLATE_PARAM( 1, WrappedDimension );

  constexpr PhaseSpaceDimension ParentDimension = WrappedParentDimension::value;
  constexpr PhaseSpaceDimension Dimension = WrappedDimension::value;

}

//---------------------------------------------------------------------------//
// Test if the distribution is uniform
FRENSIE_UNIT_TEST_TEMPLATE_EXPAND( DependentPhaseSpaceDimensionDistribution,
                                   isUniform,
                                   TestPhaseSpaceDimensions )
{
  FETCH_TEMPLATE_PARAM( 0, WrappedParentDimension );
  FETCH_TEMPLATE_PARAM( 1, WrappedDimension );

  constexpr PhaseSpaceDimension ParentDimension = WrappedParentDimension::value;
  constexpr PhaseSpaceDimension Dimension = WrappedDimension::value;
}

//---------------------------------------------------------------------------//
// Test if the distribution is uniform
FRENSIE_UNIT_TEST_TEMPLATE_EXPAND( DependentPhaseSpaceDimensionDistribution,
                                   hasForm,
                                   TestPhaseSpaceDimensions )
{
  FETCH_TEMPLATE_PARAM( 0, WrappedParentDimension );
  FETCH_TEMPLATE_PARAM( 1, WrappedDimension );

  constexpr PhaseSpaceDimension ParentDimension = WrappedParentDimension::value;
  constexpr PhaseSpaceDimension Dimension = WrappedDimension::value;
  
}

//---------------------------------------------------------------------------//
// Test if the distribution can be evaluated without a cascade
FRENSIE_UNIT_TEST_TEMPLATE_EXPAND( DependentPhaseSpaceDimensionDistribution,
                                   evaluateWithoutCascade,
                                   TestPhaseSpaceDimensions )
{
  FETCH_TEMPLATE_PARAM( 0, WrappedParentDimension );
  FETCH_TEMPLATE_PARAM( 1, WrappedDimension );

  constexpr PhaseSpaceDimension ParentDimension = WrappedParentDimension::value;
  constexpr PhaseSpaceDimension Dimension = WrappedDimension::value;

}

//---------------------------------------------------------------------------//
// Test if the distribution can be sampled without a cascade
FRENSIE_UNIT_TEST_TEMPLATE_EXPAND( DependentPhaseSpaceDimensionDistribution,
                                   sampleWithoutCascade,
                                   TestPhaseSpaceDimensions )
{
  FETCH_TEMPLATE_PARAM( 0, WrappedParentDimension );
  FETCH_TEMPLATE_PARAM( 1, WrappedDimension );

  constexpr PhaseSpaceDimension ParentDimension = WrappedParentDimension::value;
  constexpr PhaseSpaceDimension Dimension = WrappedDimension::value;


}

//---------------------------------------------------------------------------//
// Test if the distribution can be sampled without a cascade and the trials
// can be counted
FRENSIE_UNIT_TEST_TEMPLATE_EXPAND( DependentPhaseSpaceDimensionDistribution,
                                   sampleAndRecordTrialsWithoutCascade,
                                   TestPhaseSpaceDimensions )
{
  FETCH_TEMPLATE_PARAM( 0, WrappedParentDimension );
  FETCH_TEMPLATE_PARAM( 1, WrappedDimension );

  constexpr PhaseSpaceDimension ParentDimension = WrappedParentDimension::value;
  constexpr PhaseSpaceDimension Dimension = WrappedDimension::value;


}

//---------------------------------------------------------------------------//
// Test that the dimension value can be set and weighed approperiately
FRENSIE_UNIT_TEST_TEMPLATE_EXPAND( DependentPhaseSpaceDimensionDistribution,
                                   setDimensionValueAndApplyWeight,
                                   TestPhaseSpaceDimensions )
{

}

//---------------------------------------------------------------------------//
// Check that the distribution can be archived
FRENSIE_UNIT_TEST_TEMPLATE_EXPAND( DependentPhaseSpaceDimensionDistribution_indexDimension,
                                   archive,
                                   TestArchives )
{

}

//---------------------------------------------------------------------------//
// Custom setup
//---------------------------------------------------------------------------//
FRENSIE_CUSTOM_UNIT_TEST_SETUP_BEGIN();

FRENSIE_CUSTOM_UNIT_TEST_INIT()
{
  // Set up relevant index dimension objects
    std::vector<double> x_planes( {0.0, 0.5, 1.0} ),
                        y_planes( {0.0, 0.5} ),
                        z_planes( {0.0, 0.5} );

  std::shared_ptr<Utility::StructuredHexMesh> hex_mesh = std::make_shared<Utility::StructuredHexMesh>(x_planes, y_planes, z_planes);

  std::shared_ptr<Utility::PQLAQuadrature> direction_discretization = std::make_shared<Utility::PQLAQuadrature>(2);

  MonteCarlo::PhaseSpaceDimensionTraits<MonteCarlo::SPATIAL_INDEX_DIMENSION>::setMesh(hex_mesh);

  MonteCarlo::PhaseSpaceDimensionTraits<MonteCarlo::DIRECTION_INDEX_DIMENSION>::setDirectionDiscretization(direction_discretization);

  // Index dimensions require integer grids - each case must be considered
  // Create distribution where parent distribution is index (only 2 mesh voxels)
  std::vector<size_t> spatial_index_primary_grid {0, 1 };

  std::vector<std::shared_ptr<const Utility::TabularUnivariateDistribution> >
    secondary_dists( 2 );

  // Create the secondary distribution in the first bin
  secondary_dists[0] = std::make_shared<Utility::UniformDistribution>( 0.5, 0.9, 0.5 );
  // Create the secondary distribution in the second bin
  secondary_dists[1] = std::make_shared<Utility::UniformDistribution>( 0.6, 0.8, 0.4 );


  std::shared_ptr<Utility::HistogramFullyTabularBasicBivariateDistribution> local_spatial_index_parent_raw_distribution = 
    std::make_shared<Utility::HistogramFullyTabularBasicBivariateDistribution>(spatial_index_primary_grid, secondary_dists);

  std::vector<size_t> independent_values{1, 2, 3};
  std::vector<double> dependent_values{ 3.0, 5.0, 6.0};
  std::shared_ptr<Utility::DiscreteDistribution> distribution_test = std::make_shared<Utility::DiscreteDistribution>(independent_values, dependent_values);
  
  size_t sample = distribution_test->sample();
  // Create distribution where parent distribution is direction index
  std::vector<size_t> direction_index_primary_grid(direction_discretization->getNumberOfTriangles());
  for(size_t direction_index = 0; direction_index < direction_index_primary_grid.size(); direction_index++)
  {
    direction_index_primary_grid[direction_index] = direction_index;
  }

  std::shared_ptr<const Utility::BasicBivariateDistribution> parent_index_raw_distribution;

  // Create distribution where child distribution is index
  std::shared_ptr<const Utility::BasicBivariateDistribution> child_index_raw_distribution;

  // Create distribution where both are indexed
  std::shared_ptr<const Utility::BasicBivariateDistribution> both_index_raw_distribution;
  
  // Initialize the random number generator
  Utility::RandomNumberGenerator::createStreams();
}

FRENSIE_CUSTOM_UNIT_TEST_SETUP_END();

//---------------------------------------------------------------------------//
// end tstDependentPhaseSpaceDimensionDistribution.cpp
//---------------------------------------------------------------------------//
