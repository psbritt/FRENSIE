//---------------------------------------------------------------------------//
//!
//! \file   MonteCarlo_IncoherentPhotonScatteringDistributionFactory.hpp
//! \author Alex Robinson
//! \brief  The incoherent photon scattering distribution factory base class
//!
//---------------------------------------------------------------------------//

#ifndef MONTE_CARLO_INCOHERENT_PHOTON_SCATTERING_DISTRIBUTION_FACTORY_HPP
#define MONTE_CARLO_INCOHERENT_PHOTON_SCATTERING_DISTRIBUTION_FACTORY_HPP

// Trilinos Includes
#include <Teuchos_RCP.hpp>

// FRENSIE Includes
#include "MonteCarlo_IncoherentPhotonScatteringDistribution.hpp"

namespace MonteCarlo{

//! The incoherent photon scattering distribution factory base class.
class IncoherentPhotonScatteringDistributionFactory
{

public:

  //! Create a Klein-Nishina distribution
  static void createKleinNishinaDistribution( 
		    Teuchos::RCP<const IncoherentPhotonScatteringDistribution>&
		    incoherent_distribution,
		    const double kahn_sampling_cutoff_energy );
};

} // end MonteCarlo namespace

#endif // end MONTE_CARLO_INCOHERENT_PHOTON_SCATTERING_DISTRIBUTION_FACTORY_HPP

//---------------------------------------------------------------------------//
// end MonteCarlo_IncoherentPhotonScatteringDistributionFactory.hpp
//---------------------------------------------------------------------------//
