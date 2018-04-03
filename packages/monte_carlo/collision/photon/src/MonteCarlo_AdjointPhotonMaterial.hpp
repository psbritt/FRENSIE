//---------------------------------------------------------------------------//
//!
//! \file   MonteCarlo_AdjointPhotonMaterial.hpp
//! \author Alex Robinson
//! \brief  Adjoint photon material class declaration
//!
//---------------------------------------------------------------------------//

#ifndef MONTE_CARLO_ADJOINT_PHOTON_MATERIAL_HPP
#define MONTE_CARLO_ADJOINT_PHOTON_MATERIAL_HPP

// FRENSIE Includes
#include "MonteCarlo_AdjointMaterial.hpp"
#include "MonteCarlo_AdjointPhotoatom.hpp"

namespace MonteCarlo{

//! The adjoint photon material class
class AdjointPhotonMaterial : public AdjointMaterial<AdjointPhotoatom>
{
  // Typedef for QuantityTraits
  typedef Utility::QuantityTraits<double> QT;

  // Typedef for the base type
  typedef AdjointMaterial<AdjointPhotoatom> BaseType;

public:

  //! The material handle type
  typedef typename BaseType::InternalMaterialHandle InternalMaterialHandle;

  //! The adjoint photoatom name map type
  typedef typename BaseType::ScatteringCenterNameMap AdjointPhotoatomNameMap;

  //! Constructor
  AdjointPhotonMaterial(
                  const InternalMaterialHandle id,
                  const double density,
                  const AdjointPhotoatomNameMap& adjoint_photoatom_name_map,
                  const std::vector<double>& adjoint_photoatom_fractions,
                  const std::vector<std::string>& adjoint_photoatom_names );

  //! Destructor
  ~AdjointPhotonMaterial()
  { /* ... */ }
};
  
} // end MonteCarlo namespace

#endif // end MONTE_CARLO_ADJOINT_PHOTON_MATERIAL_HPP

//---------------------------------------------------------------------------//
// end MonteCarlo_AdjointPhotonMaterial.hpp
//---------------------------------------------------------------------------//
