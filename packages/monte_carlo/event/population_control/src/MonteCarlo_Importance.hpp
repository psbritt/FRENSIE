//---------------------------------------------------------------------------//
//!
//! \file   MonteCarlo_Importance.hpp
//! \author Philip Britt
//! \brief  Importance class declaration
//!
//---------------------------------------------------------------------------//

#ifndef MONTE_CARLO_IMPORTANCE_HPP
#define MONTE_CARLO_IMPORTANCE_HPP

// FRENSIE Includes
#include "MonteCarlo_PopulationControl.hpp"


namespace MonteCarlo{

//! The importance base class
class Importance: public PopulationControl
{

public:

  //! Constructor
  Importance()
  : d_is_max_split_set(false)
  { /* ... */ }

  //! Destructor
  virtual ~Importance()
  { /* ... */ }

  void checkParticleWithPopulationController( ParticleState& particle, 
                                              ParticleBank& bank) const;

  void setMaxSplit( const unsigned max_split_integer );

  virtual double getImportance( ParticleState& particle ) const = 0;

  virtual bool isParticleInImportanceDiscretization( ParticleState& particle ) const = 0;

private:

  // Declare the boost serialization access object as a friend
  friend class boost::serialization::access;

  // Serialize the data
  template<typename Archive>
  void serialize( Archive& ar, const unsigned version )
  { 
    // Serialize the base class data
    ar & BOOST_SERIALIZATION_BASE_OBJECT_NVP( PopulationControl );

    ar & BOOST_SERIALIZATION_NVP( d_is_max_split_set );
  }

  //! max splitting value may not be recommended for importances in all cases.
  bool d_is_max_split_set;

};

} // end MonteCarlo namespace

BOOST_SERIALIZATION_CLASS_VERSION( Importance, MonteCarlo, 0 );
BOOST_SERIALIZATION_ASSUME_ABSTRACT_CLASS( Importance, MonteCarlo );
EXTERN_EXPLICIT_CLASS_SERIALIZE_INST( MonteCarlo, Importance );

#endif // end MONTE_CARLO_WEIGHT_WINDOW_HPP

//---------------------------------------------------------------------------//
// end MonteCarlo_Importance.hpp
//---------------------------------------------------------------------------//
