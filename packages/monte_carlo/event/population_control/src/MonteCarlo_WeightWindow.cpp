//---------------------------------------------------------------------------//
//!
//! \file   MonteCarlo_WeightWindowMesh.cpp
//! \author Philip Britt
//! \brief  Weight window mesh class definition
//!
//---------------------------------------------------------------------------//

// std includes
#include <math.h>

// FRENSIE Includes
#include "FRENSIE_Archives.hpp"
#include "MonteCarlo_WeightWindow.hpp"
#include "Utility_ExceptionTestMacros.hpp"

namespace MonteCarlo{

WeightWindowBase::WeightWindowBase()
{ 
  d_max_split = 5;
}

// Update the particle state and bank
void WeightWindowBase::checkParticleWithPopulationController( ParticleState& particle,
                                                              ParticleBank& bank ) const
{

  //! Make sure there is a weight window where this particle is.
  if(this->isParticleInWeightWindowDiscretization( particle ))
  {

    const WeightWindow& window = this->getWeightWindow(particle);
    double weight = particle.getWeight();

    if(weight > window.upper_weight)
    {
      // return number after decimal
      unsigned number_of_particles = static_cast<unsigned>(floor(weight/window.upper_weight) + 1);
      if( number_of_particles > d_max_split ) number_of_particles = d_max_split;
      this->splitParticle(particle, bank, number_of_particles);
    }
    else if(weight < window.lower_weight)
    {
      this->terminateParticle(particle, 
                              1 - (weight/window.survival_weight));
    }

  }
}

void WeightWindowBase::setMaxSplit( const unsigned max_split_integer )
{
  testPrecondition( max_split_integer > 0 );
  d_max_split = max_split_integer;
}

} // end MonteCarlo namespace

EXPLICIT_CLASS_SERIALIZE_INST( MonteCarlo::WeightWindowBase );

//---------------------------------------------------------------------------//
// end MonteCarlo_WeightWindowMesh.cpp
//---------------------------------------------------------------------------//