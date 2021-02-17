//---------------------------------------------------------------------------//
//!
//! \file   MonteCarlo_PQLATypeObserverDirectionDimensionDiscretization.hpp
//! \author Philip Britt
//! \brief  PQLA Direction Discretization declaration
//!
//---------------------------------------------------------------------------//

#ifndef MONTE_CARLO_PQLA_TYPE_OBSERVER_DIRECTION_DIMENSION_DISCRETIZATION
#define MONTE_CARLO_PQLA_TYPE_OBSERVER_DIRECTION_DIMENSION_DISCRETIZATION

// STD includes
#include <memory>

// FRENSIE includes
#include "MonteCarlo_ObserverDirectionDimensionDiscretization.hpp"
#include "Utility_ExplicitSerializationTemplateInstantiationMacros.hpp"
#include "Utility_SerializationHelpers.hpp"
#include "Utility_PQLAQuadrature.hpp"

// Boost Includes
#include <boost/any.hpp>
#include <boost/serialization/split_member.hpp>
#include <boost/serialization/version.hpp>
#include <boost/serialization/assume_abstract.hpp>
#include <boost/serialization/export.hpp>
#include <boost/serialization/shared_ptr.hpp>


namespace MonteCarlo{


class PQLATypeObserverDirectionDimensionDiscretization: public ObserverDirectionDimensionDiscretization
{

  public:

  //! Default constructor (should throw DBC error)
  PQLATypeObserverDirectionDimensionDiscretization()
  :d_pqla_quadrature_handler(0)
  { /* ... */ }

  //! Constructor
  PQLATypeObserverDirectionDimensionDiscretization(unsigned quadrature_order,
                                                   bool forward_bin = true);

  //! Destructor
  ~PQLATypeObserverDirectionDimensionDiscretization()
  { /* ... */ }

  //! Return number of bins
  size_t getNumberOfBins() const final override;

  //! calculate bin index for direction
  void calculateBinIndicesOfValue( const ObserverParticleStateWrapper& particle_state_wrapper,
                                          BinIndexArray& bin_indices) const final override;
                                      
  //! Calculate the index of bins that the value falls in
  void calculateBinIndicesOfValue( const ObserverParticleStateWrapper& particle_state_wrapper,
                                          BinIndexWeightPairArray& bin_indices_and_weights ) const final override;

  //! Calculate the index of bins that the value falls in
  void calculateBinIndicesOfValue(  const boost::any& any_value,
                                            BinIndexArray& bin_indices ) const final override;

    //! Calculate the index of bins that the value range falls in
  void calculateBinIndicesOfRange( const ObserverParticleStateWrapper& particle_state_wrapper,
                                  BinIndexWeightPairArray& bin_indices_and_weights ) const final override;

  //! Print the boundaries of a bin. Will be implemented with source biasing changes due to needing the same information for both.
  void printBoundariesOfBin( std::ostream& os,
                             const size_t bin_index ) const final override;

  //! Print the dimension discretization. Will be implemented with source biasing changes due to needing the same information for both.
  void print( std::ostream& os ) const final override;
  
  private: 

  //! Return the triangle bin index
  unsigned returnTriangleBin( const ObserverParticleStateWrapper& particle_state_wrapper ) const;

  //! PQLA object that handles PQLA math
  std::shared_ptr<Utility::PQLAQuadrature> d_pqla_quadrature_handler;

  //! Stores whether or not we're reverse binning (for VR purposes)
  bool d_forward_bin;
  
  // Declare the boost serialization access object as a friend
  friend class boost::serialization::access;

  // Serialize the observer discretization
  template<typename Archive>
  void serialize( Archive& ar, const unsigned version )
  { 
    ar & BOOST_SERIALIZATION_BASE_OBJECT_NVP( ObserverDirectionDimensionDiscretization );
    ar & BOOST_SERIALIZATION_NVP( d_pqla_quadrature_handler );
    ar & BOOST_SERIALIZATION_NVP( d_forward_bin );
  }

};

} // end MonteCarlo namespace

BOOST_SERIALIZATION_CLASS_VERSION( PQLATypeObserverDirectionDimensionDiscretization , MonteCarlo, 0 );
BOOST_SERIALIZATION_CLASS_EXPORT_STANDARD_KEY(PQLATypeObserverDirectionDimensionDiscretization, MonteCarlo);
EXTERN_EXPLICIT_CLASS_SERIALIZE_INST( MonteCarlo, PQLATypeObserverDirectionDimensionDiscretization);

#endif // end MONTE_CARLO_PQLA_OBSERVER_DIRECTION_DIMENSION_DISCRETIZATION

//---------------------------------------------------------------------------//
// end MonteCarlo_PQLATypeObserverDirectionDimensionDiscretization.hpp
//---------------------------------------------------------------------------//