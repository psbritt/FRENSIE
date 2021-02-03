//---------------------------------------------------------------------------//
//!
//! \file   MonteCarlo_GenericHistogramImportanceParticleDistribution.cpp
//! \author Philip Britt
//! \brief  Generic thesis related importance sampled source
//!
//---------------------------------------------------------------------------//

// std includes
#include<iterator>

// FRENSIE Includes
#include "FRENSIE_Archives.hpp"
#include "MonteCarlo_GenericHistogramImportanceParticleDistribution.hpp"
#include "MonteCarlo_PhaseSpacePoint.hpp"
#include "Utility_BasicCartesianCoordinateConversionPolicy.hpp"
#include "Utility_BasicSphericalCoordinateConversionPolicy.hpp"
#include "Utility_DeltaDistribution.hpp"
#include "Utility_UniformDistribution.hpp"
#include "Utility_ToStringTraits.hpp"
#include "Utility_LoggingMacros.hpp"
#include "Utility_DesignByContract.hpp"
#include "Utility_HistogramDistribution.hpp"
#include "MonteCarlo_ImportanceSampledIndependentPhaseSpaceDimensionDistribution.hpp"

namespace MonteCarlo{
// Default Constructor
GenericHistogramImportanceParticleDistribution::GenericHistogramImportanceParticleDistribution( )
  : GenericHistogramImportanceParticleDistribution( "Generic Histogram Importance Sampled Particle Distribution" )
{ /* ... */ }

// Constructor
GenericHistogramImportanceParticleDistribution::GenericHistogramImportanceParticleDistribution( const std::string& name )
  : ParticleDistribution( name ),
    d_spatial_coord_conversion_policy( new Utility::BasicCartesianCoordinateConversionPolicy ),
    d_directional_coord_conversion_policy( new Utility::BasicCartesianCoordinateConversionPolicy ),
    d_spatial_dimension_mesh(),
    d_direction_dimension_discretization()
{ /* ... */ }

// Set a dimension distribution
void GenericHistogramImportanceParticleDistribution::setIndependentDimensionDistribution(
                        const std::shared_ptr< PhaseSpaceDimensionDistribution >& dimension_distribution,
                        const bool finished_with_independent  )
{
// Check that pointer is valid
testPrecondition(dimension_distribution);

d_dimension_order.push_back( dimension_distribution->getDimension() );

std::vector< std::shared_ptr< PhaseSpaceDimensionDistribution > > temp_vector;

temp_vector.push_back(dimension_distribution);

d_dimension_distributions[dimension_distribution->getDimension()] = temp_vector;

d_independent_finished = finished_with_independent;

}

void GenericHistogramImportanceParticleDistribution::setImportanceDimensionDistributions(
                        const std::map< PhaseSpaceDimension, std::vector< std::shared_ptr< PhaseSpaceDimensionDistribution> > > importance_sampled_distributions,
                        const std::map< PhaseSpaceDimension, std::vector< double > > importance_distribution_boundaries,
                        const std::vector< PhaseSpaceDimension > importance_distribution_order)
{
  for(auto vector_it = importance_distribution_order.begin(); vector_it != importance_distribution_order.end(); ++vector_it)
  {
    auto current_distribution_vector = importance_sampled_distributions.find(*vector_it)->second;
    auto current_distribution_boundary_vector = importance_distribution_boundaries.find(*vector_it)->second;
    if( vector_it != importance_distribution_order.end() - 1)
    {
      auto next_distribution_vector = importance_sampled_distributions.find(*(std::next(vector_it)))->second;
      // Defensive pre-processing to make sure distributions are right size
      if( next_distribution_vector.size() != current_distribution_vector.size()*(current_distribution_boundary_vector.size()-1))
        THROW_EXCEPTION(std::runtime_error, "Mismatched boundary/distribution vector size on" << *(std::next(vector_it)));
    }


    d_dimension_order.push_back(*vector_it);
    d_dimension_distributions[*vector_it] = current_distribution_vector;
    d_dimension_bounds[*vector_it] = current_distribution_boundary_vector;
  }
}

void GenericHistogramImportanceParticleDistribution::setMeshIndexDimensionDistributionObject( std::shared_ptr< Utility::StructuredHexMesh> mesh_object )
{
  d_spatial_dimension_mesh = std::make_shared<Utility::StructuredHexMesh>( mesh_object->getXPlanesCopy(),
                                                                           mesh_object->getYPlanesCopy(),
                                                                           mesh_object->getZPlanesCopy() );
  PhaseSpaceDimensionTraits<SPATIAL_INDEX_DIMENSION>::setMesh(mesh_object);
}

void GenericHistogramImportanceParticleDistribution::setDirectionIndexDimensionDistributionObject( std::shared_ptr< Utility::PQLAQuadrature> quadrature_pointer )
{
  d_direction_dimension_discretization = quadrature_pointer;
  PhaseSpaceDimensionTraits<DIRECTION_INDEX_DIMENSION>::setDirectionDiscretization(quadrature_pointer);
}

// Return the dimension distribution type name
std::string GenericHistogramImportanceParticleDistribution::getDimensionDistributionTypeName(
                                    const PhaseSpaceDimension dimension ) const
{
  return "Does not make sense to implement with this signature - ignoring";
}

// Check if the distribution is spatially uniform (somewhere)
/*! \details This method will simply check if the region where the
 * spatial distributions are defined correspond to uniform spatial
 * distribution of particles.
 */
bool GenericHistogramImportanceParticleDistribution::isSpatiallyUniform() const
{
  return true;
}

// Check if the distribution is directionally uniform (isotropic)
/*! \details This method will simply check if the solid angle where
 * the directional distributions are defined correspond to an isotropic
 * distribution of particle directions.
 */
bool GenericHistogramImportanceParticleDistribution::isDirectionallyUniform() const
{
  return true;
}

// Initialize dimension counter map
/*! \details A counter for each dimension will be created and set to 0.
 */
void GenericHistogramImportanceParticleDistribution::initializeDimensionCounters(
                                            DimensionCounterMap& trials ) const
{
  trials[PRIMARY_SPATIAL_DIMENSION] = 0;
  trials[SECONDARY_SPATIAL_DIMENSION] = 0;
  trials[TERTIARY_SPATIAL_DIMENSION] = 0;

  trials[PRIMARY_DIRECTIONAL_DIMENSION] = 0;
  trials[SECONDARY_DIRECTIONAL_DIMENSION] = 0;
  trials[TERTIARY_DIRECTIONAL_DIMENSION] = 0;

  trials[ENERGY_DIMENSION] = 0;
  trials[TIME_DIMENSION] = 0;
  trials[WEIGHT_DIMENSION] = 0;

  trials[SPATIAL_INDEX_DIMENSION] = 0;
  trials[DIRECTION_INDEX_DIMENSION] = 0;
}

// Evaluate the distribution at the desired phase space point
/*! \details The dimension distribution dependency tree must be constructed
 * before the particle distribution can be evaluated. All dimensions except
 * for the weight dimension will be used.
 */
double GenericHistogramImportanceParticleDistribution::evaluate(
                                          const ParticleState& particle ) const
{


  // Convert the particle to a phase space point
  PhaseSpacePoint phase_space_point( particle,
                                     d_spatial_coord_conversion_policy,
                                     d_directional_coord_conversion_policy );

  double distribution_value = 1.0;

  // Evaluate the distribution at the phase space point
  std::map<PhaseSpaceDimension,std::vector< std::shared_ptr<PhaseSpaceDimensionDistribution> > >::const_iterator
    dimension_dist_it = d_dimension_distributions.begin();

  // Each dimension distribution will be evaluated individually (without
  // cascading to dependent distributions)
  while( dimension_dist_it != d_dimension_distributions.end() )
  {
    if( dimension_dist_it->second[0]->getDimension() != WEIGHT_DIMENSION )
    {
      distribution_value *=
        dimension_dist_it->second[0]->evaluateWithoutCascade( phase_space_point );
    }

    ++dimension_dist_it;
  }

  return distribution_value;

}

// Sample a particle state from the distribution
/*! \details The dimension distribution dependency tree must be constructed
 * before the particle distribution can be evaluated.
 */
void GenericHistogramImportanceParticleDistribution::sample( ParticleState& particle ) const
{

  // Create the dimension sampling functor
  DimensionSamplingFunction dimension_sample_functor =
    std::bind<void>( &PhaseSpaceDimensionDistribution::sampleWithCascade,
                     std::placeholders::_1,
                     std::placeholders::_2 );

  this->sampleImpl( dimension_sample_functor, particle );
}

// Sample a particle state from the dist. and record the number of trials
/*! \details The dimension distribution dependency tree must be constructed
 * before the particle distribution can be evaluated.
 */
void GenericHistogramImportanceParticleDistribution::sampleAndRecordTrials(
                                       ParticleState& particle,
                                       DimensionCounterMap& trials ) const
{

  // Create the dimension sampling functor
  DimensionSamplingFunction dimension_sample_functor =
    std::bind<void>(
            &PhaseSpaceDimensionDistribution::sampleAndRecordTrialsWithCascade,
            std::placeholders::_1,
            std::placeholders::_2,
            std::ref( trials ) );

  this->sampleImpl( dimension_sample_functor, particle );
}

// Sample a particle state with the desired dimension value
/*! \details The dimension distribution dependency tree must be constructed
 * before the particle distribution can be evaluated.
 */
void GenericHistogramImportanceParticleDistribution::sampleWithDimensionValue(
                                           ParticleState& particle,
                                           const PhaseSpaceDimension dimension,
                                           const double dimension_value ) const
{

  // Create the dimension sampling functor
  DimensionSamplingFunction dimension_sample_functor =
      std::bind<void>(
       &PhaseSpaceDimensionDistribution::sampleWithCascadeUsingDimensionValue,
       std::placeholders::_1,
       std::placeholders::_2,
       dimension,
       dimension_value );

  this->sampleImpl( dimension_sample_functor, particle );
}

// Sample a particle state with the desired dim. value and record trials
/*! \details The dimension distribution dependency tree must be constructed
 * before the particle distribution can be evaluated.
 */
void GenericHistogramImportanceParticleDistribution::sampleWithDimensionValueAndRecordTrials(
                                           ParticleState& particle,
                                           DimensionCounterMap& trials,
                                           const PhaseSpaceDimension dimension,
                                           const double dimension_value ) const
{

  // Create the dimension sampling functor
  DimensionSamplingFunction dimension_sample_functor =
      std::bind<void>(
       &PhaseSpaceDimensionDistribution::sampleAndRecordTrialsWithCascadeUsingDimensionValue,
       std::placeholders::_1,
       std::placeholders::_2,
       std::ref( trials ),
       dimension,
       dimension_value );

  this->sampleImpl( dimension_sample_functor, particle );
}

EXPLICIT_CLASS_SAVE_LOAD_INST( GenericHistogramImportanceParticleDistribution );

} // end MonteCarlo namespace

BOOST_CLASS_EXPORT_IMPLEMENT( MonteCarlo::GenericHistogramImportanceParticleDistribution );

//---------------------------------------------------------------------------//
// end MonteCarlo_GenericHistogramImportanceParticleDistribution.cpp
//---------------------------------------------------------------------------//
