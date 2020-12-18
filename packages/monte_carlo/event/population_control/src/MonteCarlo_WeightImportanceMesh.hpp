//---------------------------------------------------------------------------//
//!
//! \file   MonteCarlo_WeightImportanceMesh.hpp
//! \author Philip Britt
//! \brief  Weight importance mesh class declaration
//!
//---------------------------------------------------------------------------//

#ifndef MONTE_CARLO_WEIGHT_IMPORTANCE_MESH_HPP
#define MONTE_CARLO_WEIGHT_IMPORTANCE_MESH_HPP

// Std Lib Includes
#include <memory>

// FRENSIE Includes
#include "MonteCarlo_WeightImportance.hpp"
#include "Utility_Mesh.hpp"
#include "Utility_Map.hpp"

namespace MonteCarlo{

//! The importance mesh class
class WeightImportanceMesh : public WeightImportance
{

public:

  //! Constructor
  WeightImportanceMesh();

  //! Destructor
  ~WeightImportanceMesh()
  { /* ... */ }

  //! Set the mesh for a particle
  void setMesh( const std::shared_ptr<const Utility::Mesh> mesh );

  //! Set the discretization map for the importance mesh (vector index is discretization index)
  void setWeightImportanceMap( std::unordered_map<Utility::Mesh::ElementHandle, std::vector<double>>& weight_importance_map );

  double getWeightImportance( ParticleState& particle) const final override;

  bool isParticleInWeightImportanceDiscretization( ParticleState& particle ) const final override;

  std::shared_ptr<const Utility::Mesh> getMesh() const;

  const std::unordered_map<Utility::Mesh::ElementHandle, std::vector<double>>& getWeightImportanceMap() const;

private:

  // Declare the boost serialization access object as a friend
  friend class boost::serialization::access;

  std::shared_ptr<const Utility::Mesh > d_mesh;

  //! Map that contains importance windows. First key is the index of the mesh element, second key is the index of the discretization
  std::unordered_map<Utility::Mesh::ElementHandle, std::vector<double>> d_weight_importance_map;

  // Serialize the data
  template<typename Archive>
  void serialize( Archive& ar, const unsigned version )
  {   
    // Serialize the base class data
    ar & BOOST_SERIALIZATION_BASE_OBJECT_NVP( WeightImportance );
    // Serialize the member data
    ar & BOOST_SERIALIZATION_NVP( d_mesh );
    ar & BOOST_SERIALIZATION_NVP( d_weight_importance_map );
  }

};

} // end MonteCarlo namespace

  
BOOST_SERIALIZATION_CLASS_VERSION( WeightImportanceMesh, MonteCarlo, 0 );
EXTERN_EXPLICIT_CLASS_SERIALIZE_INST( MonteCarlo, WeightImportanceMesh );

#endif // end MONTE_CARLO_WEIGHT_IMPORTANCE_MESH_HPP

//---------------------------------------------------------------------------//
// end MonteCarlo_ImportanceMesh.hpp
//---------------------------------------------------------------------------//
