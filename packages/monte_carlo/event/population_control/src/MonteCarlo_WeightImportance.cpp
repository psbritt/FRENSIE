//---------------------------------------------------------------------------//
//!
//! \file   MonteCarlo_Importance.cpp
//! \author Philip Britt
//! \brief  Importance mesh class definition
//!
//---------------------------------------------------------------------------//

// std includes
#include <cmath>
#include <functional>

// FRENSIE Includes
#include "FRENSIE_Archives.hpp"
#include "MonteCarlo_WeightImportance.hpp"
#include "Utility_ExceptionTestMacros.hpp"
#include "Utility_RandomNumberGenerator.hpp"

namespace MonteCarlo{

WeightImportance::WeightImportance()
: d_use_non_importance_weight_transforms(false)
{ 
  d_max_split = 5;
}

// Update the particle state and bank
void WeightImportance::checkParticleWithPopulationController( ParticleState& particle,
                                                              ParticleBank& bank ) const
{
  if( this->isParticleInWeightImportanceDiscretization( particle ) )
  {
    double weight_importance = this->getWeightImportance(particle);
    if( d_use_non_importance_weight_transforms )
    {
      weight_importance *= particle.getWeight()/particle.getImportanceWeightTransform();
    }
    double particle_weight = particle.getWeight();

    // Don't do anything if the weight is within a tolerance (Mostly to keep from splitting/terminating on birth)
    if( fabs( weight_importance-particle_weight) > d_weight_importance_tolerance)
    {
      double weight_importance_fraction = particle_weight / weight_importance;

      // Split probabilistically
      if( weight_importance_fraction > 1 )
      {
        // Split particle into lower possible number of emergent particles
        if( weight_importance_fraction > d_max_split )
        {
            this->splitParticle( particle,
                                  bank,
                                  d_max_split );
        }
        else
        {
          this->splitParticle(particle,
                              bank,
                              weight_importance_fraction);
        }
      }
      // Terminate
      else
      {
        this->terminateParticle(particle, 1-weight_importance_fraction); 
      }
    }
  }
}

void WeightImportance::setMaxSplit( const unsigned max_split_integer )
{
  testPrecondition( max_split_integer > 0 );
  d_max_split = max_split_integer;
}

void WeightImportance::setNonImportanceWeightTransform( const bool use_non_importance_weight_transforms)
{
  d_use_non_importance_weight_transforms = use_non_importance_weight_transforms;
}

} // end MonteCarlo namespace

EXPLICIT_CLASS_SERIALIZE_INST( MonteCarlo::WeightImportance );

//---------------------------------------------------------------------------//
// end MonteCarlo_Importance.cpp
//---------------------------------------------------------------------------//