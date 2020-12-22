//---------------------------------------------------------------------------//
//!
//! \file   MonteCarlo_WeightImportance.hpp
//! \author Philip Britt
//! \brief  Weight importance class declaration
//!
//---------------------------------------------------------------------------//

#ifndef MONTE_CARLO_WEIGHT_IMPORTANCE_HPP
#define MONTE_CARLO_WEIGHT_IMPORTANCE_HPP

// FRENSIE Includes
#include "MonteCarlo_PopulationControl.hpp"


namespace MonteCarlo{

//! The weight importance base class
class WeightImportance: public PopulationControl
{

public:

  //! Constructor
  WeightImportance();

  //! Destructor
  virtual ~WeightImportance()
  { /* ... */ }

  void checkParticleWithPopulationController( ParticleState& particle, 
                                              ParticleBank& bank) const;

  void setMaxSplit( const unsigned max_split_integer );
  
  // Change IC because not capturing only IC
  void setNonImportanceWeightTransform( const bool use_non_importance_weight_transforms );

  virtual double getWeightImportance( ParticleState& particle ) const = 0;

  virtual bool isParticleInWeightImportanceDiscretization( ParticleState& particle ) const = 0;

private:

  // Declare the boost serialization access object as a friend
  friend class boost::serialization::access;

  // Serialize the data
  template<typename Archive>
  void serialize( Archive& ar, const unsigned version )
  { 
    // Serialize the base class data
    ar & BOOST_SERIALIZATION_BASE_OBJECT_NVP( PopulationControl );
  }

  double d_weight_importance_tolerance = 1e-13;

  bool d_use_non_importance_weight_transforms;
};

} // end MonteCarlo namespace

BOOST_SERIALIZATION_CLASS_VERSION( WeightImportance, MonteCarlo, 0 );
BOOST_SERIALIZATION_ASSUME_ABSTRACT_CLASS( WeightImportance, MonteCarlo );
EXTERN_EXPLICIT_CLASS_SERIALIZE_INST( MonteCarlo, WeightImportance );

#endif // end MONTE_CARLO_WEIGHT_IMPORTANCE_HPP

//---------------------------------------------------------------------------//
// end MonteCarlo_WeightImportance.hpp
//---------------------------------------------------------------------------//
