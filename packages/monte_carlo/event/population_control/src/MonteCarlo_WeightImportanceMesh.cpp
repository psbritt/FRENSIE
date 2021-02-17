//---------------------------------------------------------------------------//
//!
//! \file   MonteCarlo_WeightImportanceMesh.cpp
//! \author Philip Britt
//! \brief  Weight importance mesh class definition
//!
//---------------------------------------------------------------------------//

// FRENSIE Includes
#include "FRENSIE_Archives.hpp"
#include "MonteCarlo_WeightImportanceMesh.hpp"
#include "Utility_ExceptionTestMacros.hpp"
#include "Utility_DesignByContract.hpp"

namespace MonteCarlo{

// Constructor
WeightImportanceMesh::WeightImportanceMesh()
{ /* ... */ }

// Set the mesh for a particle
void WeightImportanceMesh::setMesh(const std::shared_ptr<const Utility::Mesh> mesh)
{
  d_mesh = mesh;
}

void WeightImportanceMesh::setWeightImportanceMap( std::unordered_map<Utility::Mesh::ElementHandle, std::vector<double>>& weight_importance_map )
{
  testPrecondition(d_mesh);
  d_weight_importance_map = weight_importance_map;
}

double WeightImportanceMesh::getWeightImportance( ParticleState& particle) const
{
  ObserverParticleStateWrapper observer_particle(particle);
  ObserverPhaseSpaceDimensionDiscretization::BinIndexArray discretization_index;
  this->calculateBinIndicesOfPoint(observer_particle, discretization_index);
  return d_weight_importance_map.at(d_mesh->whichElementIsPointIn(particle.getPosition()))[discretization_index[0]];
}

bool WeightImportanceMesh::isParticleInWeightImportanceDiscretization( ParticleState& particle ) const
{
  ObserverParticleStateWrapper observer_particle(particle);

  return(d_mesh->isPointInMesh(particle.getPosition()) && this->isPointInObserverPhaseSpace(observer_particle));
}

std::shared_ptr<const Utility::Mesh> WeightImportanceMesh::getMesh() const
{
  return d_mesh;
}

const std::unordered_map<Utility::Mesh::ElementHandle, std::vector<double>>& WeightImportanceMesh::getWeightImportanceMap() const
{
  return d_weight_importance_map;
}

} // end MonteCarlo namespace

EXPLICIT_CLASS_SERIALIZE_INST( MonteCarlo::WeightImportanceMesh );
BOOST_SERIALIZATION_CLASS_EXPORT_IMPLEMENT( WeightImportanceMesh , MonteCarlo );

//---------------------------------------------------------------------------//
// end MonteCarlo_WeightImportanceMesh.cpp
//---------------------------------------------------------------------------//
