//---------------------------------------------------------------------------//
//!
//! \file   Facemc_NeutronScatteringAngularDistribution.hpp
//! \author Alex Robinson, Alex Bennett
//! \brief  The neutron scattering angular distribution base class declaration
//!
//---------------------------------------------------------------------------//

#ifndef FACEMC_NEUTRON_SCATTERING_ANGULAR_DISTRIBUTION
#define FACEMC_NEUTRON_SCATTERING_ANGULAR_DISTRIBUTION

// Trilinos Includse
#include <Teuchos_Array.hpp>

// FRENSIE Includes
#include "Facemc_NeutronState.hpp"
#include "Utility_PhysicalConstants.hpp"
#include "Utility_RandomNumberGenerator.hpp"
#include "Utility_OneDDistribution.hpp"

namespace Facemc{

//! The angular scattering distribution base class
class NeutronScatteringAngularDistribution
{

public:

  //! Typedef for energy dependent angular distribution
  typedef Teuchos::Array<Utility::Pair<double,
			       Teuchos::RCP<Utility::OneDDistribution> > > 
  AngularDistribution;
  
  //! Constructor
  NeutronScatteringAngularDistribution( const AngularDistribution& dist );

  //! Destructor
  virtual ~NeutronScatteringAngularDistribution()
  { /* ... */ }

  //! Sample a scattering angle cosine
  virtual double sampleAngleCosine( const double energy ) const;

private:

  // The angular distribution
  AngularDistribution d_angular_distribution;
};

} // end Facemc namespace

#endif // end FACEMC_NEUTRON_SCATTERING_ANGULAR_DISTRIBUTION

//---------------------------------------------------------------------------//
// end Facemc_NeutronScatteringAngularDistribution.hpp
//---------------------------------------------------------------------------//
