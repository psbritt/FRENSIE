//---------------------------------------------------------------------------//
//!
//! \file   DataGen_SubshellIncoherentCrossSectionEvaluator.hpp
//! \author Alex Robinson
//! \brief  Subshell incoherent cross section evaluator declaration
//!
//---------------------------------------------------------------------------//

#ifndef DATA_GEN_SUBSHELL_INCOHERENT_CROSS_SECTION_EVALUATOR_HPP
#define DATA_GEN_SUBSHELL_INCOHERENT_CROSS_SECTION_EVALUATOR_HPP

// Boost Includes
#include <boost/scoped_ptr.hpp>

// Trilinos Includes
#include <Teuchos_RCP.hpp>

// FRENSIE Includes
#include "Utility_OneDDistribution.hpp"
#include "Utility_KleinNishinaDistribution.hpp"
#include "Utility_GaussKronrodQuadratureKernel.hpp"
#include "MonteCarlo_SubshellType.hpp"

namespace DataGen{

//! The subshell incoherent cross section class
class SubshellIncoherentCrossSectionEvaluator
{

public:

  //! Default constructor (free electron)
  SubshellIncoherentCrossSectionEvaluator();

  //! Constructor (bound electron)
  SubshellIncoherentCrossSectionEvaluator(
                           const MonteCarlo::SubshellType subshell,
			   const double num_electrons_in_subshell,
			   const double subshell_binding_energy,
                           const Teuchos::RCP<const Utility::OneDDistribution>&
			   subshell_occupation_number );

  //! Return the maximum electron momentum projection
  double getMaxElectronMomentumProjection( const double energy,
					   const double angle_cosine ) const;

  //! Evaluate the differential adjoint incoherent cross section (dc/dx)
  double evaluateDifferentialCrossSection(
		 const double inverse_energy_gain_ratio,
		 const Utility::KleinNishinaDistribution& distribution ) const;

private:

  // The subshell
  MonteCarlo::SubshellType d_subshell;

  // The number of electrons that occupy the subshell
  double d_num_electrons_in_subshell;

  // The binding energy of the subshell (MeV)
  double d_subshell_binding_energy

  // The occupation number of the subshell
  Teuchos::RCP<const Utility::OneDDistribution> d_subshell_occupation_number;
};

} // end DataGen namespace

#endif // end DATA_GEN_SUBSHELL_INCOHERENT_CROSS_SECTION_EVALUATOR_HPP

//---------------------------------------------------------------------------//
// end DataGen_SubshellIncoherentCrossSectionEvaluator.hpp
//---------------------------------------------------------------------------//
