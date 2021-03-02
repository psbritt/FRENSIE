//---------------------------------------------------------------------------//
//!
//! \file   MonteCarlo_GenericHistogramImportanceParticleDistribution_def.hpp
//! \author Philip Britt
//! \brief  Generic histogram importance particle distribution template definitions
//!
//---------------------------------------------------------------------------//

#ifndef MONTE_CARLO_GENERIC_HISTOGRAM_IMPORTANCE_PARTICLE_DISTRIBUTION_DEF_HPP
#define MONTE_CARLO_GENERIC_HISTOGRAM_IMPORTANCE_PARTICLE_DISTRIBUTION_DEF_HPP

// FRENSIE includes
#include "MonteCarlo_PhaseSpaceDimensionTraits.hpp"

// std includes
#include <iterator>

namespace MonteCarlo{

// Sample the particle state using the desired dimension sampling functor
template<typename DimensionSamplingFunctor>
inline void GenericHistogramImportanceParticleDistribution::sampleImpl(
                         DimensionSamplingFunctor& dimension_sampling_function,
                         ParticleState& particle ) const
{
  
  // Initialize a phase space point
  PhaseSpacePoint phase_space_sample( d_spatial_coord_conversion_policy,
                                      d_directional_coord_conversion_policy );

  size_t current_distribution_index = 0;
  for( auto order_it = d_dimension_order.begin(); order_it != d_dimension_order.end(); ++order_it)
  {
    std::cout << *order_it << " index: " << current_distribution_index << std::endl;
    auto distribution_vector = d_dimension_distributions.find(*order_it)->second;
    dimension_sampling_function(
                *(distribution_vector[current_distribution_index]),
                phase_space_sample );

    std::cout << "Size of distribution vector for " << *order_it << ": " << distribution_vector.size() << std::endl;
    // check if next dimension has more than 1 distribution
    if( distribution_vector.size() > 1 && order_it != d_dimension_order.end()-1 )
    {
      const std::vector<double>& dimension_boundaries_vector = d_dimension_bounds.find(*order_it)->second;
      // Get the value of the sample so a search can be done
      double sample_value = phase_space_sample.getCoordinate(*order_it);
      size_t temp_index = 0;
      for(auto vec_it = dimension_boundaries_vector.begin(); vec_it != dimension_boundaries_vector.end()-1; ++vec_it)
      {
        if( *vec_it <= sample_value && sample_value < *(std::next(vec_it)))
        {
         current_distribution_index = temp_index + (dimension_boundaries_vector.size()-1)*current_distribution_index;
         break;
        }
        temp_index += 1;
      }
    }
    
  }
  
  // Convert the sampled phase space point to a particle state. This will
  // use the spatial and directional conversion policies
  phase_space_sample.setParticleState( particle );

}

// Save the state to an archive
template<typename Archive>
void GenericHistogramImportanceParticleDistribution::save( Archive& ar, const unsigned version ) const
{
  // Save the base class member data
  ar & BOOST_SERIALIZATION_BASE_OBJECT_NVP( ParticleDistribution );

  // Save the local member data
  ar & BOOST_SERIALIZATION_NVP( d_spatial_coord_conversion_policy );
  ar & BOOST_SERIALIZATION_NVP( d_directional_coord_conversion_policy );
  ar & BOOST_SERIALIZATION_NVP( d_dimension_order );
  ar & BOOST_SERIALIZATION_NVP( d_dimension_distributions );
  ar & BOOST_SERIALIZATION_NVP( d_dimension_bounds );
  ar & BOOST_SERIALIZATION_NVP( d_spatial_dimension_mesh );
  ar & BOOST_SERIALIZATION_NVP( d_direction_dimension_discretization );  
}

// Load the data from an archive
template<typename Archive>
void GenericHistogramImportanceParticleDistribution::load( Archive& ar, const unsigned version )
{
  // Load the base class member data
  ar & BOOST_SERIALIZATION_BASE_OBJECT_NVP( ParticleDistribution );

  // Load the local member data
  ar & BOOST_SERIALIZATION_NVP( d_spatial_coord_conversion_policy );
  ar & BOOST_SERIALIZATION_NVP( d_directional_coord_conversion_policy );
  ar & BOOST_SERIALIZATION_NVP( d_dimension_order );
  ar & BOOST_SERIALIZATION_NVP( d_dimension_distributions );
  ar & BOOST_SERIALIZATION_NVP( d_dimension_bounds );
  ar & BOOST_SERIALIZATION_NVP( d_spatial_dimension_mesh );
  ar & BOOST_SERIALIZATION_NVP( d_direction_dimension_discretization ); 

  if (d_spatial_dimension_mesh)
  {
    PhaseSpaceDimensionTraits<SPATIAL_INDEX_DIMENSION>::setMesh(d_spatial_dimension_mesh);
  }
  if (d_direction_dimension_discretization)
  {
    PhaseSpaceDimensionTraits<DIRECTION_INDEX_DIMENSION>::setDirectionDiscretization(d_direction_dimension_discretization);
  }
}
  
} // end MonteCarlo namespace

#endif // end MONTE_CARLO_STANDARD_PARTICLE_DISTRIBUTION_DEF_HPP

//---------------------------------------------------------------------------//
// end MonteCarlo_GenericHistogramImportanceParticleDistribution_def.hpp
//---------------------------------------------------------------------------//
