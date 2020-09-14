//---------------------------------------------------------------------------//
//!
//! \file   MonteCarlo_WeightWindow.hpp
//! \author Philip Britt
//! \brief  Importance class declaration
//!
//---------------------------------------------------------------------------//

#ifndef MONTE_CARLO_IMPORTANCE_HPP
#define MONTE_CARLO_IMPORTANCE_HPP

// FRENSIE Includes
#include "MonteCarlo_PopulationControl.hpp"


namespace MonteCarlo{

//! The weight window base class
class Importance: public PopulationControl
{

public:

  //! Constructor
  Importance()
  { /* ... */ }

  //! Destructor
  virtual ~Importance()
  { /* ... */ }

  void checkParticleWithPopulationController( ParticleState& particle, 
                                              ParticleBank& bank) const;

protected:

  virtual double getImportance( ParticleState& particle ) const = 0;

  virtual bool isParticleInImportanceDiscretization( ParticleState& particle ) const = 0;

private:

  // Declare the boost serialization access object as a friend
  friend class boost::serialization::access;

  // Serialize the data
  template<typename Archive>
  void serialize( Archive& ar, const unsigned version )
  { /* ... */ }

};

} // end MonteCarlo namespace

BOOST_CLASS_VERSION( MonteCarlo::Importance, 0 );
BOOST_SERIALIZATION_ASSUME_ABSTRACT( MonteCarlo::Importance );
EXTERN_EXPLICIT_CLASS_SERIALIZE_INST( MonteCarlo, Importance );

#endif // end MONTE_CARLO_WEIGHT_WINDOW_HPP

//---------------------------------------------------------------------------//
// end MonteCarlo_WeightWindow.hpp
//---------------------------------------------------------------------------//