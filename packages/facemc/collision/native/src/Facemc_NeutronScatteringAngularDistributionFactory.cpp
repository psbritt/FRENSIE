//---------------------------------------------------------------------------//
//! 
//! \file   Facemc_NeutronScatteringAngularDistributionFactory.cpp
//! \author Alex Robinson, Alex Bennett
//! \brief  Neutron scattering angular distribution factory class declaration
//!
//---------------------------------------------------------------------------//

// Trilinos Includes
#include <Teuchos_ArrayView.hpp>

// FRENSIE Includes
#include "Facemc_NeutronScatteringAngularDistributionFactory.hpp"
#include "Utility_UniformDistribution.hpp"
#include "Utility_HistogramDistribution.hpp"
#include "Utility_TabularDistribution.hpp"
#include "Utility_ThirtyTwoEquiprobableBinDistribution.hpp"
#include "Utility_ExceptionTestMacros.hpp"
#include "Utility_ContractException.hpp"

namespace Facemc{

// Initialize the static member data
Teuchos::RCP<Utility::OneDDistribution> 
NeutronScatteringAngularDistributionFactory::isotropic_angle_cosine_dist(
			  new Utility::UniformDistribution( -1.0, 1.0, 1.0 ) );

// Create the angular distribution
void NeutronScatteringAngularDistributionFactory::createDistribution(
	     const Teuchos::ArrayView<const double>& and_block_array,
	     const unsigned and_block_array_start_index,
	     const std::string& table_name,
	     const NuclearReactionType reaction,
	     Teuchos::RCP<NeutronScatteringAngularDistribution>& distribution )
{
  // Make sure the and block array is "valid"
  testPrecondition( and_block_array.size() > 0 );

  // Get the number of energies at which the angular distribution is tabulated
  unsigned num_tabulated_energies = static_cast<unsigned>(and_block_array[0]);

  // Get the energy grid
  Teuchos::ArrayView<const double> energy_grid = 
    and_block_array( 1, num_tabulated_energies);
  
  // Get the location of the angular distribution for each energy
  Teuchos::ArrayView<const double> distribution_indices = 
    and_block_array( num_tabulated_energies + 1,
		     num_tabulated_energies );
  
  // Initialize the angular distribution array
  NeutronScatteringAngularDistribution::AngularDistribution
    angular_distribution( num_tabulated_energies );

  for( unsigned i = 0u; i < energy_grid.size(); ++i )
  {
    angular_distribution[i].first = energy_grid[i];

    int distribution_index = 
      static_cast<int>( distribution_indices[i] );

    // Thirty two equiprobable bin distribution
    if( distribution_index > 0 )
    {
      Teuchos::ArrayView<const double> bin_boundaries = 
	and_block_array( distribution_index, 33 );

      angular_distribution[i].second.reset( 
	 new Utility::ThirtyTwoEquiprobableBinDistribution( bin_boundaries ) );
    }

    // Tabular distribution
    else if( distribution_index < 0 )
    {
      // Distribution index is relative to beginning of AND block - subtract
      // off start index of portion of and block for given MT #.
      distribution_index = abs(distribution_index) - 1 -
	and_block_array_start_index;
 
      unsigned interpolation_flag = 
	and_block_array[distribution_index];
      
      unsigned number_of_points_in_dist = 
	and_block_array[distribution_index + 1];
      
      Teuchos::ArrayView<const double> scattering_angle_cosine_grid = 
	and_block_array( distribution_index + 2,
			 number_of_points_in_dist );

      // Ignore the last evaluated point in the PDF
      Teuchos::ArrayView<const double> pdf;

      switch( interpolation_flag )
      {
      case 1u: // histogram interpolation
	pdf = and_block_array( 
			     distribution_index + 2 + number_of_points_in_dist,
			     number_of_points_in_dist - 1u );
	
	angular_distribution[i].second.reset( 
	      new Utility::HistogramDistribution( scattering_angle_cosine_grid,
						  pdf ) );
	break;
	
      case 2u: // Linear-Linear interpolation
	pdf = and_block_array( 
			     distribution_index + 2 + number_of_points_in_dist,
			     number_of_points_in_dist );
	
	angular_distribution[i].second.reset( 
	                 new Utility::TabularDistribution<Utility::LinLin>(
						  scattering_angle_cosine_grid,
						  pdf ) );
						 
	break;
	
      default:
	TEST_FOR_EXCEPTION( true,
			    std::runtime_error,
			    "Unknown interpolation flag found in table " 
			    << table_name << 
			    " for angular distribution of MT = "
			    << reaction << ": "
			    << interpolation_flag << "\n" );
      }
    }
    
    // Isotropic distribution
    else
    {
      angular_distribution[i].second = isotropic_angle_cosine_dist;
    }
  }

  // Create the angular distribution
  distribution.reset( 
	    new NeutronScatteringAngularDistribution( angular_distribution ) );
}

// Create an isotropic angular distribution
void NeutronScatteringAngularDistributionFactory::createIsotropicDistribution(
	     Teuchos::RCP<NeutronScatteringAngularDistribution>& distribution )
{
  NeutronScatteringAngularDistribution::AngularDistribution
    angular_distribution( 2 );

  angular_distribution[0].first = 1e-11;
  angular_distribution[0].second = isotropic_angle_cosine_dist;

  angular_distribution[1].first = 20.0;
  angular_distribution[1].second = isotropic_angle_cosine_dist;

  distribution.reset( 
	    new NeutronScatteringAngularDistribution( angular_distribution ) );
}

} // end Facemc namespace

//---------------------------------------------------------------------------//
// end Facemc_NeutronScatteringAngularDistributionFactory.cpp
//---------------------------------------------------------------------------//
