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
{
  // Set default value for max split parameter (able to be changed)
  d_max_split = 5;
  d_use_IC_weight_transforms = false;
}

// Update the particle state and bank
void WeightImportance::checkParticleWithPopulationController( ParticleState& particle,
                                                              ParticleBank& bank ) const
{
  if( this->isParticleInWeightImportanceDiscretization( particle ) )
  {
    double weight_importance = this->getWeightImportance(particle);
    if( d_use_IC_weight_transforms )
    {
      weight_importance *= particle.getICWeightTransform();
    }
    double particle_weight = particle.getWeight();

    // Don't do anything if the weight is within a tolerance (Mostly to keep from splitting/terminating on birth)
    if( fabs( weight_importance-particle_weight) > d_weight_importance_tolerance)
    {
      double weight_importance_fraction = particle_weight / weight_importance;

      // Split probabilistically
      if( weight_importance_fraction > 1 )
      {
        double rounded_weight_importance_fraction;
        // Split particle into lower possible number of emergent particles
        if( weight_importance_fraction > d_max_split )
        {
            this->splitParticle( particle,
                                  bank,
                                  d_max_split );
        }
        else
        {
          if( Utility::RandomNumberGenerator::getRandomNumber<double>() < 1-std::fmod(weight_importance_fraction, 1))
          {
            rounded_weight_importance_fraction = std::floor(weight_importance_fraction);
          }
          // Split particle into greater possible number of emergent particles
          else
          {
            rounded_weight_importance_fraction = std::ceil(weight_importance_fraction);
          }
          this->splitParticle(particle,
                              bank,
                              static_cast<unsigned>(rounded_weight_importance_fraction),
                              1/weight_importance_fraction);
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

void WeightImportance::setICWeightTransform( const bool use_IC_weight_transform)
{
  d_use_IC_weight_transforms = use_IC_weight_transform;
}

} // end MonteCarlo namespace

EXPLICIT_CLASS_SERIALIZE_INST( MonteCarlo::WeightImportance );

//---------------------------------------------------------------------------//
// end MonteCarlo_Importance.cpp
//---------------------------------------------------------------------------//