
//---------------------------------------------------------------------------//
//!
//! \file   MonteCarlo_PhaseSpaceDimensionTraits.cpp
//! \author Philip Britt
//! \brief  Particle source dimension traits class specializations
//!
//---------------------------------------------------------------------------//

// FRENSIE Includes
#include "MonteCarlo_PhaseSpaceDimensionTraits.hpp"

namespace MonteCarlo{
  
  std::shared_ptr<Utility::StructuredHexMesh> PhaseSpaceDimensionTraits<SPATIAL_INDEX_DIMENSION>::s_mesh;
  std::shared_ptr<Utility::PQLAQuadrature> PhaseSpaceDimensionTraits<DIRECTION_INDEX_DIMENSION>::s_direction_discretization;

} // end MonteCarlo namespace

//---------------------------------------------------------------------------//
// end MonteCarlo_PhaseSpaceDimensionTraits.cpp
//---------------------------------------------------------------------------//