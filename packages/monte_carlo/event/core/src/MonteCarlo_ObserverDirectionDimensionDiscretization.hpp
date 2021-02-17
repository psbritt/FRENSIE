//---------------------------------------------------------------------------//
//!
//! \file   MonteCarlo_ObserverDirectionDimensionDiscretization.hpp
//! \author Philip Britt
//! \brief  Direction dimension discretization specialization 
//!
//---------------------------------------------------------------------------//

#ifndef MONTE_CARLO_OBSERVER_DIRECTION_DIMENSION_DISCRETIZATION_HPP
#define MONTE_CARLO_OBSERVER_DIRECTION_DIMENSION_DISCRETIZATION_HPP

// FRENSIE includes
#include "MonteCarlo_ObserverPhaseSpaceDimensionDiscretization.hpp"
#include "Utility_ExplicitSerializationTemplateInstantiationMacros.hpp"
#include "Utility_SerializationHelpers.hpp"

// Boost Includes
#include <boost/any.hpp>
#include <boost/serialization/split_member.hpp>
#include <boost/serialization/version.hpp>
#include <boost/serialization/assume_abstract.hpp>
#include <boost/serialization/export.hpp>
#include <boost/serialization/shared_ptr.hpp>

namespace MonteCarlo{

//! The direction dimension discretization type enumeration
enum ObserverDirectionDiscretizationType{
  PQLA = 0,
};

class ObserverDirectionDimensionDiscretization : public ObserverPhaseSpaceDimensionDiscretization
{

public: 



  //! Constructor
  ObserverDirectionDimensionDiscretization()
  { /* ... */ }

  //! Destructor
  virtual ~ObserverDirectionDimensionDiscretization()
  { /* ... */}

  //! Class that always returns direction dimension
  ObserverPhaseSpaceDimension getDimension() const final override;

  //! Always returns Direction as name
  std::string getDimensionName() const final override;

  static std::string name()
  {
    return "Direction";
  };

  //! Return number of bins - depends on discretization type
  virtual size_t getNumberOfBins() const = 0;

  //! Always returns true. Direction discretizations currently cover entire direction unit sphere
  bool isValueInDiscretization( const ObserverParticleStateWrapper& particle_state_wrapper ) const final override;

  bool isValueInDiscretization( const boost::any& any_value ) const final override;

  //! Check if range intersects discretization - always returns true
  bool doesRangeIntersectDiscretization( const ObserverParticleStateWrapper& particle_state_wrapper ) const final override;

  //! calculate bin index for direction
  virtual void calculateBinIndicesOfValue( const ObserverParticleStateWrapper& particle_state_wrapper,
                                          BinIndexArray& bin_indices) const = 0;
                                      
  //! Calculate the index of bins that the value falls in
  virtual void calculateBinIndicesOfValue( const ObserverParticleStateWrapper& particle_state_wrapper,
                                          BinIndexWeightPairArray& bin_indices_and_weights ) const = 0;

    //! Calculate the index of bins that the value falls in
  virtual void calculateBinIndicesOfValue( const boost::any& any_value,
                                           BinIndexArray& bin_indices ) const = 0;

    //! Calculate the index of bins that the value range falls in
  virtual void calculateBinIndicesOfRange( const ObserverParticleStateWrapper& particle_state_wrapper,
                                           BinIndexWeightPairArray& bin_indices_and_weights ) const = 0;

  //! Print the boundaries of a bin
  virtual void printBoundariesOfBin( std::ostream& os,
                                     const size_t bin_index ) const = 0;

  //! Print the dimension discretization
  virtual void print( std::ostream& os ) const = 0;

  private:

  // Serialize the data
  template<typename Archive>
  void serialize( Archive& ar, const unsigned version )
  { 
    ar & BOOST_SERIALIZATION_BASE_OBJECT_NVP( ObserverPhaseSpaceDimensionDiscretization );
  }

  // Declare the boost serialization access object as a friend
  friend class boost::serialization::access;
};

} // end MonteCarlo namespace
namespace boost{
namespace serialization{

//! Serialize the direction discretization type enum
template<typename Archive>
void serialize( Archive& archive,
                MonteCarlo::ObserverDirectionDiscretizationType& discretization_type,
                const unsigned version )
{
  if( Archive::is_saving::value )
    archive & (int)discretization_type;
  else
  {
    int raw_discretization_type;

    archive & raw_discretization_type;

    switch( raw_discretization_type )
    {
      BOOST_SERIALIZATION_ENUM_CASE( MonteCarlo::PQLA, int, discretization_type );
      default:
      {
        THROW_EXCEPTION( std::logic_error,
                         "Unable to deserialize direction discretization type "
                         << raw_discretization_type << "!" );
      }
    }
  }
}

} // end serialization namespace

} // end boost namespace

BOOST_SERIALIZATION_CLASS_VERSION( ObserverDirectionDimensionDiscretization , MonteCarlo, 0 );
BOOST_SERIALIZATION_ASSUME_ABSTRACT_CLASS( ObserverDirectionDimensionDiscretization, MonteCarlo );
EXTERN_EXPLICIT_CLASS_SERIALIZE_INST( MonteCarlo, ObserverDirectionDimensionDiscretization);

#endif // end MONTE_CARLO_OBSERVER_DIRECTION_DIMENSION_DISCRETIZATION_HPP

//---------------------------------------------------------------------------//
// end MonteCarlo_ObserverDirectionDimensionDiscretization.hpp
//---------------------------------------------------------------------------//