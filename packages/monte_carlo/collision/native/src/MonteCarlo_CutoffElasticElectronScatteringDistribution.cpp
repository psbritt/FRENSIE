//---------------------------------------------------------------------------//
//!
//! \file   MonteCarlo_CutoffElasticElectronScatteringDistribution.cpp
//! \author Luke Kersting
//! \brief  The electron cutoff elastic scattering distribution definition
//!
//---------------------------------------------------------------------------//

// FRENSIE Includes
#include "MonteCarlo_CutoffElasticElectronScatteringDistribution.hpp"
#include "Utility_RandomNumberGenerator.hpp"
#include "Utility_SearchAlgorithms.hpp"
#include "Utility_DirectionHelpers.hpp"
#include "Utility_KinematicHelpers.hpp"
#include "Utility_ElasticElectronTraits.hpp"


namespace MonteCarlo{

// Basic Constructor
/*! \details The basic constructor should only be used when sampling from the
 *  full cutoff elastic distribution. The other constructor should be used
 *  when sampling from a subrange of the cutoff distribution.
 */
CutoffElasticElectronScatteringDistribution::CutoffElasticElectronScatteringDistribution(
    const std::shared_ptr<TwoDDist>& scattering_distribution,
    const bool correlated_sampling_mode_on )
  : d_full_cutoff_distribution( scattering_distribution ),
    d_partial_cutoff_distribution( scattering_distribution ),
    d_cutoff_angle_cosine( scattering_distribution->getUpperBoundOfConditionalIndepVar(
        scattering_distribution->getLowerBoundOfPrimaryIndepVar() ) )
{
  // Make sure the array is valid
  testPrecondition( d_full_cutoff_distribution.use_count() > 0 );
  // Make sure the bool is valid
  testPrecondition( correlated_sampling_mode_on == 0 ||
                    correlated_sampling_mode_on == 1 );

  if( correlated_sampling_mode_on )
  {
    // Set the correlated unit based sample routine
    d_sample_func = std::bind<double>(
         &TwoDDist::sampleSecondaryConditionalExact,
         std::cref( *d_full_cutoff_distribution ),
         std::placeholders::_1 );
  }
  else
  {
    // Set the stochastic unit based sample routine
    d_sample_func = std::bind<double>(
         &TwoDDist::sampleSecondaryConditional,
         std::cref( *d_full_cutoff_distribution ),
         std::placeholders::_1 );
  }
}

// Constructor
/*! \details The constructor should be used when sampling from a subrange of
 *  the cutoff distribution. The full cutoff distribution is still needed for
 *  calculating the partial cutoff cross section.
 */
CutoffElasticElectronScatteringDistribution::CutoffElasticElectronScatteringDistribution(
    const std::shared_ptr<TwoDDist>& full_scattering_distribution,
    const std::shared_ptr<TwoDDist>& partial_scattering_distribution,
    const double cutoff_angle_cosine,
    const bool correlated_sampling_mode_on )
  : d_full_cutoff_distribution( full_scattering_distribution ),
    d_partial_cutoff_distribution( partial_scattering_distribution ),
    d_cutoff_angle_cosine( cutoff_angle_cosine )
{
  // Make sure the array is valid
  testPrecondition( d_full_cutoff_distribution.use_count() > 0 );
  // Make sure the array is valid
  testPrecondition( d_partial_cutoff_distribution.use_count() > 0 );
  // Make sure the cutoff_angle_cosine is valid
  testPrecondition( cutoff_angle_cosine >= -1.0 );
  testPrecondition( cutoff_angle_cosine <= 1.0 );
  // Make sure the bool is valid
  testPrecondition( correlated_sampling_mode_on == 0 ||
                    correlated_sampling_mode_on == 1 );

  if( correlated_sampling_mode_on )
  {
    // Set the correlated unit based sample routine
    d_sample_func = std::bind<double>(
         &TwoDDist::sampleSecondaryConditionalExact,
         std::cref( *d_partial_cutoff_distribution ),
         std::placeholders::_1 );
  }
  else
  {
    // Set the stochastic unit based sample routine
    d_sample_func = std::bind<double>(
         &TwoDDist::sampleSecondaryConditional,
         std::cref( *d_partial_cutoff_distribution ),
         std::placeholders::_1 );
  }
}

// Evaluate the cutoff cross section ratio
double CutoffElasticElectronScatteringDistribution::evaluateCutoffCrossSectionRatio(
    const double incoming_energy ) const
{
  // Make sure the energy is valid
  testPrecondition( incoming_energy > 0.0 );

  return d_full_cutoff_distribution->evaluateSecondaryConditionalCDFExact(
                                                        incoming_energy,
                                                        d_cutoff_angle_cosine );
}

// Evaluate the distribution
double CutoffElasticElectronScatteringDistribution::evaluate(
                            const double incoming_energy,
                            const double scattering_angle_cosine ) const
{
  // Make sure the energy and angle are valid
  testPrecondition( incoming_energy > 0.0 );
  testPrecondition( scattering_angle_cosine >= -1.0 );
  testPrecondition( scattering_angle_cosine <= 1.0 );

  if ( scattering_angle_cosine > d_cutoff_angle_cosine )
    return 0.0;
  else
    return d_partial_cutoff_distribution->evaluateExact( incoming_energy,
                                                         scattering_angle_cosine );
}

// Evaluate the PDF
double CutoffElasticElectronScatteringDistribution::evaluatePDF(
                            const double incoming_energy,
                            const double scattering_angle_cosine ) const
{
  // Make sure the energy and angle are valid
  testPrecondition( incoming_energy > 0.0 );
  testPrecondition( scattering_angle_cosine >= -1.0 );
  testPrecondition( scattering_angle_cosine <= 1.0 );

  if ( scattering_angle_cosine > d_cutoff_angle_cosine )
    return 0.0;
  else
    return d_partial_cutoff_distribution->evaluateSecondaryConditionalPDFExact(
                        incoming_energy,
                        scattering_angle_cosine );
}

// Evaluate the CDF
double CutoffElasticElectronScatteringDistribution::evaluateCDF(
                            const double incoming_energy,
                            const double scattering_angle_cosine ) const
{
  // Make sure the energy and angle are valid
  testPrecondition( incoming_energy > 0.0 );
  testPrecondition( scattering_angle_cosine >= -1.0 );
  testPrecondition( scattering_angle_cosine <= 1.0 );

  if ( scattering_angle_cosine >= d_cutoff_angle_cosine )
    return 1.0;
  else
    return d_partial_cutoff_distribution->evaluateSecondaryConditionalCDFExact(
                        incoming_energy,
                        scattering_angle_cosine );
}

// Sample an outgoing energy and direction from the distribution
void CutoffElasticElectronScatteringDistribution::sample(
    const double incoming_energy,
    double& outgoing_energy,
    double& scattering_angle_cosine ) const
{
  // The outgoing energy is always equal to the incoming energy
  outgoing_energy = incoming_energy;

  unsigned trial_dummy;

  // Sample an outgoing direction
  this->sampleAndRecordTrialsImpl( incoming_energy,
                                   scattering_angle_cosine,
                                   trial_dummy );
}

// Sample an outgoing energy and direction and record the number of trials
void CutoffElasticElectronScatteringDistribution::sampleAndRecordTrials(
                        const double incoming_energy,
                        double& outgoing_energy,
                        double& scattering_angle_cosine,
                        unsigned& trials ) const
{
  // The outgoing energy is always equal to the incoming energy
  outgoing_energy = incoming_energy;

  // Sample an outgoing direction
  this->sampleAndRecordTrialsImpl( incoming_energy,
                                   scattering_angle_cosine,
                                   trials );
}

// Randomly scatter the electron
void CutoffElasticElectronScatteringDistribution::scatterElectron(
                     ElectronState& electron,
                     ParticleBank& bank,
                     Data::SubshellType& shell_of_interaction ) const
{
  double scattering_angle_cosine;

  unsigned trial_dummy;

  // Sample an outgoing direction
  this->sampleAndRecordTrialsImpl( electron.getEnergy(),
                                   scattering_angle_cosine,
                                   trial_dummy );

  shell_of_interaction =Data::UNKNOWN_SUBSHELL;

  // Set the new direction
  electron.rotateDirection( scattering_angle_cosine,
                            this->sampleAzimuthalAngle() );
}

// Randomly scatter the adjoint electron
void CutoffElasticElectronScatteringDistribution::scatterAdjointElectron(
                     AdjointElectronState& adjoint_electron,
                     ParticleBank& bank,
                     Data::SubshellType& shell_of_interaction ) const
{
  double scattering_angle_cosine;

  unsigned trial_dummy;

  // Sample an outgoing direction
  this->sampleAndRecordTrialsImpl( adjoint_electron.getEnergy(),
                                   scattering_angle_cosine,
                                   trial_dummy );

  shell_of_interaction = Data::UNKNOWN_SUBSHELL;

  // Set the new direction
  adjoint_electron.rotateDirection( scattering_angle_cosine,
                                    this->sampleAzimuthalAngle() );
}

// Sample an outgoing direction from the distribution
void CutoffElasticElectronScatteringDistribution::sampleAndRecordTrialsImpl(
        const double incoming_energy,
        double& scattering_angle_cosine,
        unsigned& trials ) const
{
  // Make sure the incoming energy is valid
  testPrecondition( incoming_energy > 0.0 );

  // Increment the number of trials
  ++trials;

  // sample the scattering angle cosine
  scattering_angle_cosine = d_sample_func( incoming_energy );

  // Make sure the scattering angle cosine is valid
  testPostcondition( scattering_angle_cosine >= -1.0 );
  testPostcondition( scattering_angle_cosine <= d_cutoff_angle_cosine );
}

} // end MonteCarlo namespace

//---------------------------------------------------------------------------//
// end MonteCarlo_CutoffElasticElectronScatteringDistribution.cpp
//---------------------------------------------------------------------------//
