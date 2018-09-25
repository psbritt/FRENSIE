//---------------------------------------------------------------------------//
//!
//! \file   DataGen_ENDLElectronPhotonRelaxationDataGenerator.hpp
//! \author Alex Robinson, Luke Kersting
//! \brief  The ENDL electron-photon-relaxation data generator class
//!         definition
//!
//---------------------------------------------------------------------------//

// Boost Includes
#include <boost/function.hpp>
#include <boost/bind.hpp>

// FRENSIE Includes
#include "DataGen_ENDLElectronPhotonRelaxationDataGenerator.hpp"
#include "DataGen_MomentPreservingElectronDataGenerator.hpp"
#include "DataGen_FormFactorEvaluator.hpp"
#include "DataGen_ScatteringFunctionEvaluator.hpp"
#include "DataGen_OccupationNumberEvaluator.hpp"
#include "DataGen_ElasticElectronMomentsEvaluator.hpp"
#include "DataGen_ComptonProfileGenerator.hpp"
#include "MonteCarlo_ComptonProfileHelpers.hpp"
#include "MonteCarlo_ComptonProfileSubshellConverterFactory.hpp"
#include "MonteCarlo_ElasticElectronScatteringDistributionNativeFactory.hpp"
#include "Data_SubshellType.hpp"
#include "Utility_UniformDistribution.hpp"
#include "Utility_SloanRadauQuadrature.hpp"
#include "Utility_StandardHashBasedGridSearcher.hpp"
#include "Utility_AtomicMomentumUnit.hpp"
#include "Utility_InverseAngstromUnit.hpp"
#include "Utility_QuantityTraits.hpp"
#include "Utility_StaticOutputFormatter.hpp"
#include "Utility_ExceptionTestMacros.hpp"
#include "Utility_ExceptionCatchMacros.hpp"
#include "Utility_LoggingMacros.hpp"
#include "Utility_DesignByContract.hpp"

namespace DataGen{

// Constructor
ENDLElectronPhotonRelaxationDataGenerator::ENDLElectronPhotonRelaxationDataGenerator(
     const std::shared_ptr<const Data::ENDLDataContainer>& endl_data_container,
     const double min_photon_energy,
     const double max_photon_energy,
     const double min_electron_energy,
     const double max_electron_energy )
  : ENDLElectronPhotonRelaxationDataGenerator( endl_data_container,
                                               min_photon_energy,
                                               max_photon_energy,
                                               min_electron_energy,
                                               max_electron_energy,
                                               true )
{ /* ... */ }

// Constructor (check for valid min/max particle energies)
ENDLElectronPhotonRelaxationDataGenerator::ENDLElectronPhotonRelaxationDataGenerator(
     const std::shared_ptr<const Data::ENDLDataContainer>& endl_data_container,
     const double min_photon_energy,
     const double max_photon_energy,
     const double min_electron_energy,
     const double max_electron_energy,
     const bool check_photon_energies )
  : ElectronPhotonRelaxationDataGenerator( endl_data_container->getAtomicNumber(),
                                           endl_data_container->getAtomicWeight(),
                                           min_photon_energy,
                                           max_photon_energy,
                                           min_electron_energy,
                                           max_electron_energy ),
    d_endl_data_container( endl_data_container )
{
  // Make sure the endl data is valid
  testPrecondition( endl_data_container.get() );

  if( check_photon_energies )
  {
    // Find the min table energy
    double table_min_photon_energy = this->findTableMinPhotonEnergy();

    if( min_photon_energy < table_min_photon_energy )
    {
      this->setMinPhotonEnergy( table_min_photon_energy );

      FRENSIE_LOG_WARNING( "the min photon energy requested is below the "
                           "ENDL table min photon energy! The ENDL table's "
                           "min photon energy (" << table_min_photon_energy <<
                           ") will be used instead." );
    }

    // Find the max table energy
    double table_max_photon_energy = this->findTableMaxPhotonEnergy();

    if( max_photon_energy > table_max_photon_energy )
    {
      this->setMaxPhotonEnergy( table_max_photon_energy );

      FRENSIE_LOG_WARNING( "the max photon energy requested is above the "
                           "ENDL table max photon energy! The ENDL table's "
                           "max photon energy (" << table_max_photon_energy <<
                           ") will be used instead." );
    }
  }

  // Check if the min electron energy is below the endl table min energy
  const double table_min_electron_energy =
    endl_data_container->getElasticEnergyGrid().front();

  if( min_electron_energy < table_min_electron_energy )
  {
    this->setMinElectronEnergy( table_min_electron_energy );

    FRENSIE_LOG_WARNING( "the min electron energy requested is below the "
                         << "ENDL table min electron energy! The ENDL "
                         << "table's min " << "electron energy ("
                         << table_min_electron_energy
                         << ") will be used instead." );
  }

  // Check if the max electron energy is above the endl table max energy
  const double table_max_electron_energy =
    endl_data_container->getElasticEnergyGrid().back();

  if( max_electron_energy > table_max_electron_energy )
  {
    this->setMaxElectronEnergy( table_max_electron_energy );

    FRENSIE_LOG_WARNING( "the max electron energy requested is above the "
                         << "ENDL table max electron energy! The ENDL table's "
                         << "max electron energy ("
                         << table_max_electron_energy
                         << ") will be used instead." );
  }
}

// Basic Constructor
ENDLElectronPhotonRelaxationDataGenerator::ENDLElectronPhotonRelaxationDataGenerator(
                          const std::shared_ptr<const Data::ENDLDataContainer>&
                          endl_data_container )
  : ENDLElectronPhotonRelaxationDataGenerator( endl_data_container,
                                               1e-3,
                                               1.0,
                                               endl_data_container->getElasticEnergyGrid().front(),
                                               endl_data_container->getElasticEnergyGrid().back(),
                                               false )
{
  this->setMinPhotonEnergy( this->findTableMinPhotonEnergy() );
  this->setMaxPhotonEnergy( this->findTableMaxPhotonEnergy() );
}

// Constructor (existing data container)
ENDLElectronPhotonRelaxationDataGenerator::ENDLElectronPhotonRelaxationDataGenerator(
                          const std::shared_ptr<const Data::ENDLDataContainer>&
                          endl_data_container,
                          const boost::filesystem::path& file_name_with_path )
  : ElectronPhotonRelaxationDataGenerator( file_name_with_path ),
    d_endl_data_container( endl_data_container )
{
  TEST_FOR_EXCEPTION( endl_data_container->getAtomicNumber() !=
                      this->getDataContainer().getAtomicNumber(),
                      std::runtime_error,
                      "The atomic number of the endl data container ("
                      << endl_data_container->getAtomicNumber() <<
                      ") and the native electron-photon-relaxation data "
                      "container ("
                      << endl_data_container->getAtomicNumber() <<
                      ") do not match!" );
}

// Find the table min photon energy
double ENDLElectronPhotonRelaxationDataGenerator::findTableMinPhotonEnergy() const
{
  // Find the min table energy
  double table_min_photon_energy = Utility::QuantityTraits<double>::inf();

  if( d_endl_data_container->getCoherentCrossSectionEnergyGrid().front() <
      table_min_photon_energy )
  {
    table_min_photon_energy =
      d_endl_data_container->getCoherentCrossSectionEnergyGrid().front();
  }

  if( d_endl_data_container->getIncoherentCrossSectionEnergyGrid().front() <
      table_min_photon_energy )
  {
    table_min_photon_energy =
      d_endl_data_container->getIncoherentCrossSectionEnergyGrid().front();
  }

  return table_min_photon_energy;
}

// Fine the table max photon energy
double ENDLElectronPhotonRelaxationDataGenerator::findTableMaxPhotonEnergy() const
{
  // Find the max table energy
  double table_max_photon_energy = -Utility::QuantityTraits<double>::inf();

  if( d_endl_data_container->getCoherentCrossSectionEnergyGrid().back() >
      table_max_photon_energy )
  {
    table_max_photon_energy =
      d_endl_data_container->getCoherentCrossSectionEnergyGrid().back();
  }

  if( d_endl_data_container->getIncoherentCrossSectionEnergyGrid().back() >
      table_max_photon_energy )
  {
    table_max_photon_energy =
      d_endl_data_container->getIncoherentCrossSectionEnergyGrid().back();
  }

  return table_max_photon_energy;
}

// Populate the electron-photon-relaxation data container
void ENDLElectronPhotonRelaxationDataGenerator::populateEPRDataContainer()
{
  // Set the relaxation data
  FRENSIE_LOG_PARTIAL_NOTIFICATION( Utility::Bold( "Setting the relaxation data" ) << "..." );
  FRENSIE_FLUSH_ALL_LOGS();

  this->setRelaxationData();

  FRENSIE_LOG_NOTIFICATION( Utility::BoldGreen( "done." ) );

  // Set the Compton profile data
  FRENSIE_LOG_PARTIAL_NOTIFICATION( Utility::Bold( "Setting the Compton profile data " ) << "..." );
  FRENSIE_FLUSH_ALL_LOGS();

  this->setComptonProfileData();

  FRENSIE_LOG_NOTIFICATION( Utility::BoldGreen( "done." ) );

  // Set the occupation number data
  FRENSIE_LOG_PARTIAL_NOTIFICATION( Utility::Bold( "Setting the occupation number data " ) << "..." );
  FRENSIE_FLUSH_ALL_LOGS();

  this->setOccupationNumberData();

  FRENSIE_LOG_NOTIFICATION( Utility::BoldGreen( "done." ) );

  // Set the Waller-Hartree scattering function data
  FRENSIE_LOG_PARTIAL_NOTIFICATION( Utility::Bold( "Setting the Waller-Hartree scattering function data " ) << "..." );
  FRENSIE_FLUSH_ALL_LOGS();

  this->setWallerHartreeScatteringFunctionData();

  FRENSIE_LOG_NOTIFICATION( Utility::BoldGreen( "done." ) );

  // Set the Waller-Hartree atomic form factor data
  FRENSIE_LOG_PARTIAL_NOTIFICATION( Utility::Bold( "Setting the Waller-Hartree atomic form factor data " ) << "..." );
  FRENSIE_FLUSH_ALL_LOGS();

  this->setWallerHartreeAtomicFormFactorData();

  FRENSIE_LOG_NOTIFICATION( Utility::BoldGreen( "done." ) );

  // Set the photon data
  FRENSIE_LOG_NOTIFICATION( Utility::Bold( "Setting the photon data: " ) );

  this->setPhotonData();

  FRENSIE_LOG_NOTIFICATION( Utility::BoldGreen( "done." ) );

  // Set the electron data
  FRENSIE_LOG_NOTIFICATION( Utility::Bold("Setting the electron data: " ) );

  this->setElectronData();

  FRENSIE_LOG_NOTIFICATION( Utility::BoldGreen( "done." ) );
}

// Repopulate the electron elastic data
void ENDLElectronPhotonRelaxationDataGenerator::repopulateElectronElasticDataImpl(
                         const double max_electron_energy,
                         const double cutoff_angle_cosine,
                         const double tabular_evaluation_tol,
                         const unsigned number_of_moment_preserving_angles,
                         const MonteCarlo::TwoDGridType two_d_grid,
                         const MonteCarlo::TwoDInterpolationType two_d_interp )
{
  TEST_FOR_EXCEPTION( max_electron_energy <= 0.0,
                      std::runtime_error,
                      "The max electron energy must be greater than "
                      << this->getMinElectronEnergy() << "!" );

  TEST_FOR_EXCEPTION( cutoff_angle_cosine > 1.0,
                      std::runtime_error,
                      "The cutoff angle cosine must be between -1.0 and "
                      "1.0!" );
  TEST_FOR_EXCEPTION( cutoff_angle_cosine <= -1.0,
                      std::runtime_error,
                      "The cutoff angle cosine must be between -1.0 and "
                      "1.0!" );

  TEST_FOR_EXCEPTION( tabular_evaluation_tol <= 0.0,
                      std::runtime_error,
                      "The tabular evaluation tolerance must be between 0.0 "
                      "and 1.0!" );
  TEST_FOR_EXCEPTION( tabular_evaluation_tol >= 1.0,
                      std::runtime_error,
                      "The tabular evaluation tolerance must be between 0.0 "
                      "and 1.0!" );

  Data::ElectronPhotonRelaxationVolatileDataContainer& data_container =
    this->getVolatileDataContainer();

  // Get the elastic angular energy grid
  std::vector<double> angular_energy_grid(
    data_container.getElasticAngularEnergyGrid() );

  if( max_electron_energy != angular_energy_grid.back() )
  {
    // Get the elastic cutoff data
    std::map<double,std::vector<double> > elastic_pdf(
      data_container.getCutoffElasticPDF() );
    std::map<double,std::vector<double> > elastic_angle(
      data_container.getCutoffElasticAngles() );

    // Get the upper boundary of the max energy
    std::vector<double>::iterator energy_bin =
      Utility::Search::binaryUpperBound( angular_energy_grid.begin(),
                                         angular_energy_grid.end(),
                                         max_electron_energy );

    if( *energy_bin != max_electron_energy )
    {
      std::vector<double> angles, pdf;
      double max_cutoff_angle_cosine = 1.0;

      // Get the angular grid and pdf at the max energy
      if( two_d_grid == MonteCarlo::DIRECT_GRID ||
          two_d_grid == MonteCarlo::UNIT_BASE_GRID )
      {
        if ( two_d_interp == MonteCarlo::LOGLOGLOG_INTERPOLATION )
        {
          MonteCarlo::ElasticElectronScatteringDistributionNativeFactory::getAngularGridAndPDF<Utility::LogNudgedLogCosLog,Utility::Direct>(
            angles,
            pdf,
            elastic_angle,
            elastic_pdf,
            max_electron_energy,
            max_cutoff_angle_cosine,
            tabular_evaluation_tol );
        }
        else if ( two_d_interp == MonteCarlo::LINLINLIN_INTERPOLATION )
        {
          MonteCarlo::ElasticElectronScatteringDistributionNativeFactory::getAngularGridAndPDF<Utility::LinLinLin,Utility::Direct>(
            angles,
            pdf,
            elastic_angle,
            elastic_pdf,
            max_electron_energy,
            max_cutoff_angle_cosine,
            tabular_evaluation_tol );
        }
        else if ( two_d_interp == MonteCarlo::LINLINLOG_INTERPOLATION )
        {
          MonteCarlo::ElasticElectronScatteringDistributionNativeFactory::getAngularGridAndPDF<Utility::LinLinLog,Utility::Direct>(
            angles,
            pdf,
            elastic_angle,
            elastic_pdf,
            max_electron_energy,
            max_cutoff_angle_cosine,
            tabular_evaluation_tol );
        }
        else
        {
          THROW_EXCEPTION( std::runtime_error,
                           "Error: the desired 2D interpolation policy " <<
                           two_d_interp <<
                           " is currently not supported!" );
        }
      }
      else if( two_d_grid == MonteCarlo::CORRELATED_GRID ||
               two_d_grid == MonteCarlo::UNIT_BASE_CORRELATED_GRID )
      {
        if ( two_d_interp == MonteCarlo::LOGLOGLOG_INTERPOLATION )
        {
          MonteCarlo::ElasticElectronScatteringDistributionNativeFactory::getAngularGridAndPDF<Utility::LogNudgedLogCosLog,Utility::Correlated>(
            angles,
            pdf,
            elastic_angle,
            elastic_pdf,
            max_electron_energy,
            max_cutoff_angle_cosine,
            tabular_evaluation_tol );
        }
        else if ( two_d_interp == MonteCarlo::LINLINLIN_INTERPOLATION )
        {
          MonteCarlo::ElasticElectronScatteringDistributionNativeFactory::getAngularGridAndPDF<Utility::LinLinLin,Utility::Correlated>(
            angles,
            pdf,
            elastic_angle,
            elastic_pdf,
            max_electron_energy,
            max_cutoff_angle_cosine,
            tabular_evaluation_tol );
        }
        else if ( two_d_interp == MonteCarlo::LINLINLOG_INTERPOLATION )
        {
          MonteCarlo::ElasticElectronScatteringDistributionNativeFactory::getAngularGridAndPDF<Utility::LinLinLog,Utility::Correlated>(
            angles,
            pdf,
            elastic_angle,
            elastic_pdf,
            max_electron_energy,
            max_cutoff_angle_cosine,
            tabular_evaluation_tol );
        }
        else
        {
          THROW_EXCEPTION( std::runtime_error,
                           "Error: the desired 2D interpolation policy " <<
                           two_d_interp <<
                           " is currently not supported!" );
        }
      }
      else
      {
        THROW_EXCEPTION( std::runtime_error,
                         "Error: the desired 2D grid policy " <<
                         two_d_grid <<
                         " is currently not supported!" );
      }

      elastic_angle[max_electron_energy] = angles;
      elastic_pdf[max_electron_energy] = pdf;
    }

    // Set new angular energy grid
    std::vector<double> new_energy_grid( angular_energy_grid.begin(), energy_bin );
    new_energy_grid.push_back( max_electron_energy );

    // Erase all distributions above the max electron energy
    for( energy_bin; energy_bin != angular_energy_grid.end(); ++energy_bin )
    {
      elastic_angle.erase( *energy_bin );
      elastic_pdf.erase( *energy_bin );
    }

    // Set the elastic cutoff data
    data_container.setElasticAngularEnergyGrid( new_energy_grid );
    data_container.setCutoffElasticPDF( elastic_pdf );
    data_container.setCutoffElasticAngles( elastic_angle );
  }
  else if ( max_electron_energy > angular_energy_grid.back() )
  {
    THROW_EXCEPTION( std::runtime_error,
                     "the desired max electron energy " <<
                     max_electron_energy <<
                     " is greater than the max tabulated energy " <<
                     angular_energy_grid.back() );
  }

  // Repopulate moment preserving data
  this->repopulateMomentPreservingDataImpl( cutoff_angle_cosine,
                                            tabular_evaluation_tol,
                                            number_of_moment_preserving_angles,
                                            two_d_interp );
}

// Repopulate the electron moment preserving data
void ENDLElectronPhotonRelaxationDataGenerator::repopulateMomentPreservingDataImpl(
                         const double cutoff_angle_cosine,
                         const double tabular_evaluation_tol,
                         const unsigned number_of_moment_preserving_angles,
                         const MonteCarlo::TwoDInterpolationType two_d_interp )
{
  TEST_FOR_EXCEPTION( cutoff_angle_cosine > 1.0,
                      std::runtime_error,
                      "The cutoff angle cosine must be between -1.0 and "
                      "1.0!" );
  TEST_FOR_EXCEPTION( cutoff_angle_cosine <= -1.0,
                      std::runtime_error,
                      "The cutoff angle cosine must be between -1.0 and "
                      "1.0!" );

  TEST_FOR_EXCEPTION( tabular_evaluation_tol <= 0.0,
                      std::runtime_error,
                      "The tabular evaluation tolerance must be between 0.0 "
                      "and 1.0!" );
  TEST_FOR_EXCEPTION( tabular_evaluation_tol >= 1.0,
                      std::runtime_error,
                      "The tabular evaluation tolerance must be between 0.0 "
                      "and 1.0!" );

  this->setCutoffAngleCosine( cutoff_angle_cosine );
  this->setNumberOfMomentPreservingAngles( number_of_moment_preserving_angles );
  this->setElectronTwoDInterpPolicy( two_d_interp );

  Data::ElectronPhotonRelaxationVolatileDataContainer& data_container =
    this->getVolatileDataContainer();

  std::vector<double> angular_energy_grid(
    data_container.getElasticAngularEnergyGrid() );

  FRENSIE_LOG_PARTIAL_NOTIFICATION( Utility::Bold( "Setting the moment preserving electron data ... " ) );
  FRENSIE_FLUSH_ALL_LOGS();

  // Clear the old moment preservinf data
  data_container.clearMomentPreservingData();

  // Check if moment preserving data can be generated
  if ( cutoff_angle_cosine > 0.999999 ||
       number_of_moment_preserving_angles < 1 )
  {
    FRENSIE_LOG_NOTIFICATION( Utility::BoldYellow( "done." ) );

    if( cutoff_angle_cosine > 0.999999 )
    {
      FRENSIE_LOG_NOTIFICATION( "\n" << Utility::BoldCyan( "  Note: " )
                                << "Moment preserving data was not generated "
                                << "because the cutoff angle cosine is "
                                << "greater than 0.999999 (mu = "
                                << cutoff_angle_cosine << ")." );
    }

    if( number_of_moment_preserving_angles < 1 )
    {
      FRENSIE_LOG_NOTIFICATION( "\n" << Utility::BoldCyan( "  Note: " )
                                << "Moment preserving data was not generated "
                                << "because the number of moment preserving "
                                << "angles is less than 1." );
    }
  }
  // Set the moment preserving data
  else
  {
    this->setMomentPreservingData( angular_energy_grid,
                                   tabular_evaluation_tol,
                                   two_d_interp );

    FRENSIE_LOG_NOTIFICATION( Utility::BoldGreen( "done." ) );
  }
}

// Set the atomic data
void ENDLElectronPhotonRelaxationDataGenerator::setRelaxationData()
{
  // Extract the subshell ENDF designators
  const std::set<unsigned>& subshells =
    d_endl_data_container->getSubshells();

  Data::ElectronPhotonRelaxationVolatileDataContainer& data_container =
    this->getVolatileDataContainer();

  // Assign the set of all subshells
  data_container.setSubshells( subshells );

  // Assign the subshell data
  for( auto&& subshell : subshells )
  {
    data_container.setSubshellOccupancy( subshell, d_endl_data_container->getSubshellOccupancy( subshell ) );
    data_container.setSubshellBindingEnergy( subshell, d_endl_data_container->getSubshellBindingEnergy( subshell ) );

    const std::map<unsigned,double>& radiative_transition_probs =
      d_endl_data_container->getRadiativeTransitionProbability( subshell );

    const std::map<unsigned,double>& radiative_transition_energies =
      d_endl_data_container->getRadiativeTransitionEnergy( subshell );

    const std::map<unsigned,std::map<unsigned,double> >& non_radiative_transition_probs =
      d_endl_data_container->getNonRadiativeTransitionProbability( subshell );

    const std::map<unsigned,std::map<unsigned,double> >& non_radiative_transition_energies =
      d_endl_data_container->getNonRadiativeTransitionEnergy( subshell );

    //std::cout << subshell << std::endl;

    std::vector<std::pair<unsigned,unsigned> > relaxation_vacancies;
    std::vector<double> relaxation_particle_energies, relaxation_probabilities;

    for( auto&& radiative_transition_probs_data : radiative_transition_probs )
    {
      relaxation_vacancies.push_back( std::make_pair( radiative_transition_probs_data.first, 0 ) );
      relaxation_probabilities.push_back( radiative_transition_probs_data.second );
      relaxation_particle_energies.push_back( radiative_transition_energies.find( radiative_transition_probs_data.first )->second );
    }

    for( auto&& non_radiative_transition_probs_data : non_radiative_transition_probs )
    {
      // std::cout << " " << non_radiative_transition_probs_data.first << std::endl;
      // if( subshell == 20 && non_radiative_transition_probs_data.first == 1 )
      //   std::cout << "map: " << Utility::toString(non_radiative_transition_energies) << std::endl;
      TEST_FOR_EXCEPTION( non_radiative_transition_energies.find( non_radiative_transition_probs_data.first ) == non_radiative_transition_energies.end(),
                            std::runtime_error,
                            "There are no non-radiative transition energies "
                            "for vacancy " << subshell << " transitioning to "
                            << non_radiative_transition_probs_data.first <<
                            "!" );
      
      std::map<unsigned,double>::const_iterator data_it =
        non_radiative_transition_probs_data.second.begin();

      std::map<unsigned,double>::const_iterator data_end =
        non_radiative_transition_probs_data.second.end();

      while( data_it != data_end )
      {
        relaxation_vacancies.push_back( std::make_pair( non_radiative_transition_probs_data.first, data_it->first ) );
        relaxation_probabilities.push_back( data_it->second );

        TEST_FOR_EXCEPTION( non_radiative_transition_energies.find( non_radiative_transition_probs_data.first )->second.find( data_it->first ) == non_radiative_transition_energies.find( non_radiative_transition_probs_data.first )->second.end(),
                            std::runtime_error,
                            "There are no non-radiative transition energies "
                            "for vacancy " << subshell << " transitioning to "
                            << non_radiative_transition_probs_data.first <<
                            " and " << data_it->first << "!" );
        
        relaxation_particle_energies.push_back( non_radiative_transition_energies.find( non_radiative_transition_probs_data.first )->second.find( data_it->first )->second );

        ++data_it;
      }
    }

    data_container.setSubshellRelaxationTransitions( subshell, relaxation_vacancies.size() );
    data_container.setSubshellRelaxationVacancies( subshell, relaxation_vacancies );
    data_container.setSubshellRelaxationParticleEnergies( subshell, relaxation_particle_energies );
    data_container.setSubshellRelaxationProbabilities( subshell, relaxation_probabilities );
  }
}

// Set the Compton profile data
void ENDLElectronPhotonRelaxationDataGenerator::setComptonProfileData()
{
  Data::ElectronPhotonRelaxationVolatileDataContainer& data_container =
    this->getVolatileDataContainer();

  const std::set<unsigned>& subshells = data_container.getSubshells();

  // Create the Compton profile generator
  std::shared_ptr<ComptonProfileGenerator> profile_generator =
    ComptonProfileGenerator::getDefaultGenerator( this->getAtomicNumber() );

  std::vector<double> half_momentum_grid;
  std::vector<std::vector<double> > half_profiles;

  profile_generator->generateHalfProfiles( half_profiles, half_momentum_grid );

  // Create the Compton profile subshell converter for this
  std::shared_ptr<const MonteCarlo::ComptonProfileSubshellConverter> converter;

  MonteCarlo::ComptonProfileSubshellConverterFactory::createConverter(
                                                     converter,
                                                     this->getAtomicNumber() );

  for( auto&& subshell : subshells )
  {
    unsigned compton_subshell_index = converter->convertSubshellToIndex(
                       Data::convertENDFDesignatorToSubshellEnum( subshell ) );

    std::vector<double> full_momentum_grid, full_profile;

    MonteCarlo::createFullProfileFromHalfProfile(
                                 half_momentum_grid.begin(),
                                 half_momentum_grid.end(),
                                 half_profiles[compton_subshell_index].begin(),
                                 half_profiles[compton_subshell_index].end(),
                                 full_momentum_grid,
                                 full_profile,
                                 true,
                                 false );

    // The evaluator
    std::shared_ptr<OccupationNumberEvaluator> evaluator;

  evaluator = OccupationNumberEvaluator::createEvaluator<Utility::LogLin,Utility::Units::AtomicMomentum>(
                              full_momentum_grid,
                              full_profile,
                              this->getOccupationNumberEvaluationTolerance() );

  std::function<double(double)> evaluation_wrapper =
      evaluator->getComptonProfileEvaluationWrapper();

    std::vector<double> optimized_momentum_grid( 5 ), evaluated_profile;
    optimized_momentum_grid[0] = -1.0;
    optimized_momentum_grid[1] = -0.5;
    optimized_momentum_grid[2] = 0.0;
    optimized_momentum_grid[3] = 0.5;
    optimized_momentum_grid[4] = 1.0;

    try{
      this->getDefaultPhotonGridGenerator().generateAndEvaluateInPlace(
                                                       optimized_momentum_grid,
                                                       evaluated_profile,
                                                       evaluation_wrapper );
    }
    EXCEPTION_CATCH_RETHROW(
                   std::runtime_error,
                   "Could not generate an optimized Compton "
                   "profile momentum grid for the "
                   << Data::convertENDFDesignatorToSubshellEnum( subshell ) <<
                   " subshell with the provided convergence parameters!" );

    data_container.setComptonProfileMomentumGrid( subshell,
                                                  optimized_momentum_grid );

    data_container.setComptonProfile( subshell, evaluated_profile );
  }
}

// Set the Occupation number data
void ENDLElectronPhotonRelaxationDataGenerator::setOccupationNumberData()
{
  Data::ElectronPhotonRelaxationVolatileDataContainer& data_container =
    this->getVolatileDataContainer();

  const std::set<unsigned>& subshells = data_container.getSubshells();

  std::set<unsigned>::const_iterator subshell = subshells.begin();

  // The evaluator
  std::shared_ptr<OccupationNumberEvaluator> evaluator;

  while( subshell != subshells.end() )
  {
    const std::vector<double>& momentum_grid =
      data_container.getComptonProfileMomentumGrid( *subshell );

    const std::vector<double>& compton_profile =
      data_container.getComptonProfile( *subshell );

    // Create the occupation number evaluator
    evaluator = OccupationNumberEvaluator::createEvaluator<Utility::LinLin,Utility::Units::MeCMomentum>(
                              momentum_grid,
                              compton_profile,
                              this->getOccupationNumberEvaluationTolerance() );

    std::function<double(double)> evaluation_wrapper =
      evaluator->getOccupationNumberEvaluationWrapper(
                              this->getOccupationNumberEvaluationTolerance() );

    // Create the occupation number grid
    std::vector<double> occupation_number_momentum_grid( 5 ),
      occupation_number;
    occupation_number_momentum_grid[0] = -1.0;
    occupation_number_momentum_grid[1] = -0.5;
    occupation_number_momentum_grid[2] = 0.0;
    occupation_number_momentum_grid[3] = 0.5;
    occupation_number_momentum_grid[4] = 1.0;

    try{
      this->getDefaultPhotonGridGenerator().generateAndEvaluateInPlace(
                                               occupation_number_momentum_grid,
                                               occupation_number,
                                               evaluation_wrapper );
    }
    EXCEPTION_CATCH_RETHROW(
                   std::runtime_error,
                   "Could not generate an optimized occupation "
                   "number momentum grid for the "
                   << Data::convertENDFDesignatorToSubshellEnum( *subshell ) <<
                   " subshell with the provided convergence parameters!" );

    // Fix the grid rounding errors
    std::vector<double>::iterator unity_occupation =
      std::lower_bound( occupation_number.begin(),
                        occupation_number.end(),
                        1.0 );

    while( unity_occupation != occupation_number.end() )
    {
      *unity_occupation = 1.0;

      ++unity_occupation;
    }

    data_container.setOccupationNumberMomentumGrid(
                                             *subshell,
                                             occupation_number_momentum_grid );
    data_container.setOccupationNumber( *subshell, occupation_number );

    ++subshell;
  }
}

// Set the Waller-Hartree scattering function data
void ENDLElectronPhotonRelaxationDataGenerator::setWallerHartreeScatteringFunctionData()
{
  // Create the evaluator, initial recoil momentum grid
  std::shared_ptr<const ScatteringFunctionEvaluator> evaluator;

  std::list<double> recoil_momentum_grid;

  {
    std::vector<double> raw_recoil_momentum(
      ++d_endl_data_container->getIncoherentScatteringFunctionArgument().begin(),
      d_endl_data_container->getIncoherentScatteringFunctionArgument().end() );

    std::vector<double> raw_scattering_function_values(
            ++d_endl_data_container->getIncoherentScatteringFunction().begin(),
            d_endl_data_container->getIncoherentScatteringFunction().end() );

    recoil_momentum_grid.push_back( raw_recoil_momentum.front() );
    recoil_momentum_grid.push_back( raw_recoil_momentum.back() );

    evaluator = DataGen::ScatteringFunctionEvaluator::createEvaluator<Utility::LogLog,Utility::Units::InverseCentimeter>(
                                              raw_recoil_momentum,
                                              raw_scattering_function_values );
  }

  boost::function<double (double)> grid_function =
    evaluator->getScatteringFunctionEvaluationWrapper();

  std::list<double> scattering_function;

  try{
    this->getDefaultPhotonGridGenerator().generateAndEvaluateInPlace(
                                                          recoil_momentum_grid,
                                                          scattering_function,
                                                          grid_function );
  }
  EXCEPTION_CATCH_RETHROW( std::runtime_error,
                           "Could not generate an optimized scattering "
                           "function recoil momentum grid with the provided "
                           "convergence parameters!" );

  // This grid point was left off because both the recoil momentum and the
  // scattering function are 0.0, which is incompatible with log-log
  // interpolation
  recoil_momentum_grid.push_front( 0.0 );
  scattering_function.push_front( 0.0 );

  std::vector<double> refined_recoil_momentum, refined_scattering_function;

  refined_recoil_momentum.assign( recoil_momentum_grid.begin(),
                                  recoil_momentum_grid.end() );

  refined_scattering_function.assign( scattering_function.begin(),
                                      scattering_function.end() );

  Data::ElectronPhotonRelaxationVolatileDataContainer& data_container =
    this->getVolatileDataContainer();

  data_container.setWallerHartreeScatteringFunctionMomentumGrid(
                                                     refined_recoil_momentum );
  data_container.setWallerHartreeScatteringFunction(
                                                 refined_scattering_function );
}

// Set the Waller-Hartree atomic form factor data
void ENDLElectronPhotonRelaxationDataGenerator::setWallerHartreeAtomicFormFactorData()
{
  // Create the evaluator, initial recoil momentum grid
  std::shared_ptr<const FormFactorEvaluator> evaluator;

  std::list<double> recoil_momentum_grid;

  {
    std::vector<double> raw_recoil_momentum(
              ++d_endl_data_container->getCoherentFormFactorArgument().begin(),
              d_endl_data_container->getCoherentFormFactorArgument().end() );

    std::vector<double> raw_form_factor(
                      ++d_endl_data_container->getCoherentFormFactor().begin(),
                      d_endl_data_container->getCoherentFormFactor().end() );

    recoil_momentum_grid.push_back( raw_recoil_momentum.front() );
    recoil_momentum_grid.push_back( raw_recoil_momentum.back() );

    // Create a form factor evaluator
    evaluator = DataGen::FormFactorEvaluator::createEvaluator<Utility::LogLog,Utility::Units::InverseCentimeter>(
                                                           raw_recoil_momentum,
                                                           raw_form_factor );
  }

  // Create the evaluation wrapper
  std::function<double(double)> evaluation_wrapper =
    evaluator->getFormFactorEvaluationWrapper();

  // Generate the new optimized recoil momentum grid
  std::list<double> form_factor;

  try{
    this->getDefaultPhotonGridGenerator().generateAndEvaluateInPlace(
                                                          recoil_momentum_grid,
                                                          form_factor,
                                                          evaluation_wrapper );
  }
  EXCEPTION_CATCH_RETHROW( std::runtime_error,
                           "Could not generate an optimized form "
                           "factor recoil momentum grid with the "
                           "provided convergence parameters!" );

  // The first value was kept off of the grid because the 0.0
  // value of the recoil momentum, which is not compatible with Log-Log
  // interpolation
  recoil_momentum_grid.push_front( d_endl_data_container->getCoherentFormFactorArgument().front() );
  form_factor.push_front( d_endl_data_container->getCoherentFormFactor().front() );

  std::vector<double> refined_recoil_momentum, refined_form_factor;

  refined_recoil_momentum.assign( recoil_momentum_grid.begin(),
                                  recoil_momentum_grid.end() );

  refined_form_factor.assign( form_factor.begin(),
                              form_factor.end() );

  Data::ElectronPhotonRelaxationVolatileDataContainer& data_container =
    this->getVolatileDataContainer();

  data_container.setWallerHartreeAtomicFormFactorMomentumGrid(
                                                     refined_recoil_momentum );
  data_container.setWallerHartreeAtomicFormFactor( refined_form_factor );

  // Generate the squared form factor from the newly generated form factor
  evaluator = FormFactorEvaluator::createEvaluator<Utility::LinLin,Utility::Units::InverseCentimeter>(
                                                       refined_recoil_momentum,
                                                       refined_form_factor );

  // Create the evaluation wrapper
  evaluation_wrapper = evaluator->getFormFactorSquaredEvaluationWrapper();

  // Initialize the squared recoil momentum grid
  std::vector<double> squared_recoil_momentum_grid( 2 ), squared_form_factor;

  squared_recoil_momentum_grid[0] = refined_recoil_momentum.front();
  squared_recoil_momentum_grid[0] *= squared_recoil_momentum_grid[0];

  squared_recoil_momentum_grid[1] = refined_recoil_momentum.back();
  squared_recoil_momentum_grid[1] *= squared_recoil_momentum_grid[1];

  try{
    this->getDefaultPhotonGridGenerator().generateAndEvaluateInPlace(
                                                  squared_recoil_momentum_grid,
                                                  squared_form_factor,
                                                  evaluation_wrapper );
  }
  EXCEPTION_CATCH_RETHROW( std::runtime_error,
                           "Could not generate an optimized squared "
                           "form factor squared recoil momentum grid with the "
                           "provided convergence parameters!" );

  data_container.setWallerHartreeSquaredAtomicFormFactorSquaredMomentumGrid(
                                                squared_recoil_momentum_grid );
  data_container.setWallerHartreeSquaredAtomicFormFactor( squared_form_factor );
}

// Extract the photon heating numbers
void ENDLElectronPhotonRelaxationDataGenerator::extractPhotonHeatingNumbers(
                        std::shared_ptr<const Utility::UnivariateDistribution>&
                        heating_numbers ) const
{
  // This is a dummy distribution since heating numbers are not present in
  // the ENDL tables
  heating_numbers.reset( new Utility::UniformDistribution(
                                                    this->getMinPhotonEnergy(),
                                                    this->getMaxPhotonEnergy(),
                                                    1.0 ) );
}

// Extract the Waller-Hartree incoherent cross section
void ENDLElectronPhotonRelaxationDataGenerator::extractWallerHartreeIncoherentCrossSection(
                        std::shared_ptr<const Utility::UnivariateDistribution>&
                        waller_hartree_incoherent_cs ) const
{
  waller_hartree_incoherent_cs.reset(
        new Utility::TabularDistribution<Utility::LogLog>(
                  d_endl_data_container->getIncoherentCrossSectionEnergyGrid(),
                  d_endl_data_container->getIncoherentCrossSection() ) );
}

// Extract the Waller-Hartree coherent cross section
void ENDLElectronPhotonRelaxationDataGenerator::extractWallerHartreeCoherentCrossSection(
                        std::shared_ptr<const Utility::UnivariateDistribution>&
                        waller_hartree_coherent_cs ) const
{
  waller_hartree_coherent_cs.reset(
     new Utility::TabularDistribution<Utility::LogLog>(
                    d_endl_data_container->getCoherentCrossSectionEnergyGrid(),
                    d_endl_data_container->getCoherentCrossSection() ) );
}

// Extract the pair production cross section
void ENDLElectronPhotonRelaxationDataGenerator::extractPairProductionCrossSection(
                        std::shared_ptr<const Utility::UnivariateDistribution>&
                        pair_production_cs ) const
{
  this->extractPhotonCrossSection<Utility::LogLog>(
          Utility::arrayViewOfConst( d_endl_data_container->getPairProductionCrossSectionEnergyGrid() ),
          Utility::arrayViewOfConst( d_endl_data_container->getPairProductionCrossSection() ),
          pair_production_cs,
          false );
}

// Extract the triplet production cross section
void ENDLElectronPhotonRelaxationDataGenerator::extractTripletProductionCrossSection(
                        std::shared_ptr<const Utility::UnivariateDistribution>&
                        triplet_production_cs ) const
{
  this->extractPhotonCrossSection<Utility::LogLog>(
          Utility::arrayViewOfConst( d_endl_data_container->getTripletProductionCrossSectionEnergyGrid() ),
          Utility::arrayViewOfConst( d_endl_data_container->getTripletProductionCrossSection() ),
          triplet_production_cs,
          false );
}

// Extract the subshell photoelectric effect cross section
void ENDLElectronPhotonRelaxationDataGenerator::extractSubshellPhotoelectricEffectCrossSections(
          std::vector<std::pair<unsigned,std::shared_ptr<const Utility::UnivariateDistribution> > >& cross_sections ) const
{
  const std::set<unsigned>& subshells =
    d_endl_data_container->getSubshells();

  cross_sections.resize( subshells.size() );

  int i = 0;

  for( auto&& subshell : subshells )
  {
    const Utility::ArrayView<const double> energy_grid( d_endl_data_container->getPhotoelectricCrossSectionEnergyGrid( subshell ) );

    const Utility::ArrayView<const double> cross_section( d_endl_data_container->getPhotoelectricCrossSection( subshell ) );

    cross_sections[i].first = subshell;

    this->extractPhotonCrossSection<Utility::LogLog>( energy_grid,
                                                      cross_section,
                                                      cross_sections[i].second,
                                                      false );

    ++i;
  }
}

// Create the heating numbers on the union energy grid
void ENDLElectronPhotonRelaxationDataGenerator::createHeatingNumbersOnUnionEnergyGrid(
                  const std::list<double>& union_energy_grid,
                  const std::shared_ptr<const Utility::UnivariateDistribution>&
                  original_cross_section,
                  std::vector<double>& cross_section ) const
{
  FRENSIE_LOG_WARNING( "The average heating numbers for photons cannot "
                       "currently be generated using ENDL tables!" );

  cross_section.resize( union_energy_grid.size(), 0.0 );
}

// Set the photon data
void ENDLElectronPhotonRelaxationDataGenerator::setPhotonData()
{
  // Extract the heating numbers
  std::shared_ptr<const Utility::UnivariateDistribution> heating_numbers;

  this->extractPhotonHeatingNumbers( heating_numbers );

  // Extract the Waller-Hartree incoherent cross section
  std::shared_ptr<const Utility::UnivariateDistribution>
    waller_hartree_incoherent_cs;

  this->extractWallerHartreeIncoherentCrossSection( waller_hartree_incoherent_cs );

  // Extract the Waller-Hartree coherent cross section
  std::shared_ptr<const Utility::UnivariateDistribution>
    waller_hartree_coherent_cs;

  this->extractWallerHartreeCoherentCrossSection( waller_hartree_coherent_cs );

  // Extract the pair production cross section
  std::shared_ptr<const Utility::UnivariateDistribution>
    pair_production_cs;

  this->extractPairProductionCrossSection( pair_production_cs );

  // Extract the triplet production cross section
  std::shared_ptr<const Utility::UnivariateDistribution>
    triplet_production_cs;

  this->extractTripletProductionCrossSection( triplet_production_cs );

  // Extract the subshell photoelectric effect cross sections
  std::vector<std::pair<unsigned,std::shared_ptr<const Utility::UnivariateDistribution> > >
    subshell_photoelectric_effect_css;

  this->extractSubshellPhotoelectricEffectCrossSections( subshell_photoelectric_effect_css );

  // Create the impulse approx. incoherent cross section evaluators
  std::vector<std::pair<unsigned,std::shared_ptr<const MonteCarlo::SubshellIncoherentPhotonScatteringDistribution> > >
    impulse_approx_incoherent_cs_evaluators;

  this->createSubshellImpulseApproxIncoherentCrossSectionEvaluators(
                                     impulse_approx_incoherent_cs_evaluators );

  // Create the union energy grid
  FRENSIE_LOG_PARTIAL_NOTIFICATION( " Creating " <<
                                    Utility::Italicized( "union energy grid" ) );
  FRENSIE_FLUSH_ALL_LOGS();

  std::list<double> union_energy_grid;

  this->initializePhotonUnionEnergyGrid( union_energy_grid );

  // Calculate the union energy grid
  boost::function<double (double pz)> grid_function =
    boost::bind( &Utility::UnivariateDistribution::evaluate,
                 boost::cref( *heating_numbers ),
                 _1 );

  try{
    this->getDefaultPhotonGridGenerator().generateInPlace( union_energy_grid,
                                                           grid_function );
  }
  EXCEPTION_CATCH_RETHROW( std::runtime_error,
                           "Could not generate an optimized photon "
                           "energy grid for the heating numbers with the "
                           "provided convergence parameters!" );

  FRENSIE_LOG_PARTIAL_NOTIFICATION( "." );
  FRENSIE_FLUSH_ALL_LOGS();

  grid_function = boost::bind( &Utility::UnivariateDistribution::evaluate,
                               boost::cref( *waller_hartree_incoherent_cs ),
                               _1 );

  try{
    this->getDefaultPhotonGridGenerator().generateInPlace( union_energy_grid,
                                                           grid_function );
  }
  EXCEPTION_CATCH_RETHROW( std::runtime_error,
                           "Could not generate an optimized photon "
                           "energy grid for the Waller-Hartree incoherent "
                           "cross section with the provided convergence"
                           "parameters!" );

  FRENSIE_LOG_PARTIAL_NOTIFICATION( "." );
  FRENSIE_FLUSH_ALL_LOGS();

  grid_function = boost::bind( &Utility::UnivariateDistribution::evaluate,
                               boost::cref( *waller_hartree_coherent_cs ),
                               _1 );

  try{
    this->getDefaultPhotonGridGenerator().generateInPlace( union_energy_grid,
                                                           grid_function );
  }
  EXCEPTION_CATCH_RETHROW( std::runtime_error,
                           "Could not generate an optimized photon "
                           "energy grid for the Waller-Hartree coherent "
                           "cross section with the provided convergence "
                           "parameters!" );

  FRENSIE_LOG_PARTIAL_NOTIFICATION( "." );
  FRENSIE_FLUSH_ALL_LOGS();

  grid_function = boost::bind( &Utility::UnivariateDistribution::evaluate,
                               boost::cref( *pair_production_cs ),
                               _1 );

  try{
    this->getDefaultPhotonGridGenerator().generateInPlace( union_energy_grid,
                                                           grid_function );
  }
  EXCEPTION_CATCH_RETHROW( std::runtime_error,
                           "Could not generate an optimized photon "
                           "energy grid for the pair production cross section "
                           "with the provided convergence parameters!" );

  FRENSIE_LOG_PARTIAL_NOTIFICATION( "." );
  FRENSIE_FLUSH_ALL_LOGS();

  grid_function = boost::bind( &Utility::UnivariateDistribution::evaluate,
                               boost::cref( *triplet_production_cs ),
                               _1 );

  try{
    this->getDefaultPhotonGridGenerator().generateInPlace( union_energy_grid,
                                                     grid_function );
  }
  EXCEPTION_CATCH_RETHROW( std::runtime_error,
                           "Could not generate an optimized photon "
                           "energy grid for the triplet production cross "
                           "section with the provided convergence "
                           "parameters!" );

  FRENSIE_LOG_PARTIAL_NOTIFICATION( "." );
  FRENSIE_FLUSH_ALL_LOGS();

  for( unsigned i = 0; i < subshell_photoelectric_effect_css.size(); ++i )
  {
    grid_function = boost::bind(
                   &Utility::UnivariateDistribution::evaluate,
                   boost::cref( *subshell_photoelectric_effect_css[i].second ),
                   _1 );

    try{
      this->getDefaultPhotonGridGenerator().generateInPlace( union_energy_grid,
                                                       grid_function );
    }
    EXCEPTION_CATCH_RETHROW( std::runtime_error,
                             "Could not generate an optimized photon "
                             "energy grid for the subshell "
                             << Data::convertENDFDesignatorToSubshellEnum(
                                subshell_photoelectric_effect_css[i].first ) <<
                             " photoelectric cross section with the provided "
                             "convergence parameters!" );

    FRENSIE_LOG_PARTIAL_NOTIFICATION( "." );
    FRENSIE_FLUSH_ALL_LOGS();
  }

  for( unsigned i = 0; i < impulse_approx_incoherent_cs_evaluators.size(); ++i)
  {
    grid_function = boost::bind(
             &MonteCarlo::SubshellIncoherentPhotonScatteringDistribution::evaluateIntegratedCrossSection,
             boost::cref( *impulse_approx_incoherent_cs_evaluators[i].second ),
             _1,
             this->getSubshellIncoherentEvaluationTolerance() );

    try{
      this->getDefaultPhotonGridGenerator().generateInPlace( union_energy_grid,
                                                       grid_function );
    }
    EXCEPTION_CATCH_RETHROW(
                          std::runtime_error,
                          "Could not generate an optimized photon "
                          "energy grid for the subshell "
                          << Data::convertENDFDesignatorToSubshellEnum(
                          impulse_approx_incoherent_cs_evaluators[i].first ) <<
                          " impulse approx. cross section with the "
                          "provided convergence parameters!" );

    FRENSIE_LOG_PARTIAL_NOTIFICATION( "." );
    FRENSIE_FLUSH_ALL_LOGS();
  }

  FRENSIE_LOG_NOTIFICATION( " done." );
  FRENSIE_FLUSH_ALL_LOGS();

  Data::ElectronPhotonRelaxationVolatileDataContainer& data_container =
    this->getVolatileDataContainer();

  // Set the union energy grid
  std::vector<double> energy_grid;
  energy_grid.assign( union_energy_grid.begin(), union_energy_grid.end() );
  data_container.setPhotonEnergyGrid( energy_grid );

  // Create and set the cross sections
  std::vector<double> cross_section;
  unsigned threshold;

  this->createHeatingNumbersOnUnionEnergyGrid( union_energy_grid,
                                               heating_numbers,
                                               cross_section );

  data_container.setAveragePhotonHeatingNumbers( cross_section );

  FRENSIE_LOG_PARTIAL_NOTIFICATION( " Setting the " <<
                                    Utility::Italicized( "Waller-Hartree incoherent " )
                                    << "cross section ... " );
  FRENSIE_FLUSH_ALL_LOGS();

  this->createCrossSectionOnUnionEnergyGrid( union_energy_grid,
                                             waller_hartree_incoherent_cs,
                                             cross_section,
                                             threshold,
                                             0.0,
                                             false );

  data_container.setWallerHartreeIncoherentCrossSection( cross_section );
  data_container.setWallerHartreeIncoherentCrossSectionThresholdEnergyIndex(
                                                                   threshold );

  FRENSIE_LOG_NOTIFICATION( Utility::BoldGreen( "done." ) );

  FRENSIE_LOG_PARTIAL_NOTIFICATION( " Setting the " <<
                                    Utility::Italicized( "Waller-Hartree coherent " )
                                    << "cross section ... " );
  FRENSIE_FLUSH_ALL_LOGS();

  this->createCrossSectionOnUnionEnergyGrid( union_energy_grid,
                                             waller_hartree_coherent_cs,
                                             cross_section,
                                             threshold,
                                             0.0,
                                             false );

  data_container.setWallerHartreeCoherentCrossSection( cross_section );
  data_container.setWallerHartreeCoherentCrossSectionThresholdEnergyIndex(
                                                                   threshold );

  FRENSIE_LOG_NOTIFICATION( Utility::BoldGreen( "done." ) );

  FRENSIE_LOG_PARTIAL_NOTIFICATION(" Setting the " <<
                                   Utility::Italicized( "pair production " )
                                   << "cross section ... " );
  FRENSIE_FLUSH_ALL_LOGS();

  this->createCrossSectionOnUnionEnergyGrid(
                     union_energy_grid,
                     pair_production_cs,
                     cross_section,
                     threshold,
                     2*Utility::PhysicalConstants::electron_rest_mass_energy,
                     true );

  data_container.setPairProductionCrossSection( cross_section );
  data_container.setPairProductionCrossSectionThresholdEnergyIndex(
                                                                   threshold );

  FRENSIE_LOG_NOTIFICATION( Utility::BoldGreen( "done." ) );

  FRENSIE_LOG_PARTIAL_NOTIFICATION( " Setting the " <<
                                    Utility::Italicized( "triplet production " )
                                    << "cross section ... " );
  FRENSIE_FLUSH_ALL_LOGS();

  this->createCrossSectionOnUnionEnergyGrid(
                     union_energy_grid,
                     triplet_production_cs,
                     cross_section,
                     threshold,
                     4*Utility::PhysicalConstants::electron_rest_mass_energy,
                     true );

  data_container.setTripletProductionCrossSection( cross_section );
  data_container.setTripletProductionCrossSectionThresholdEnergyIndex(
                                                                   threshold );

  FRENSIE_LOG_NOTIFICATION( Utility::BoldGreen( "done." ) );

  for( unsigned i = 0; i < subshell_photoelectric_effect_css.size(); ++i )
  {
    FRENSIE_LOG_PARTIAL_NOTIFICATION( " Setting " <<
                                      Utility::Italicized( "subshell " )
                                      << Utility::Italicized(Data::convertENDFDesignatorToSubshellEnum( subshell_photoelectric_effect_css[i].first ) )
                                      << Utility::Italicized( " photoelectric " )
                                      << "cross section ... " );
    FRENSIE_FLUSH_ALL_LOGS();

    this->createCrossSectionOnUnionEnergyGrid(
                            union_energy_grid,
                            subshell_photoelectric_effect_css[i].second,
                            cross_section,
                            threshold,
                            data_container.getSubshellBindingEnergy(
                                  subshell_photoelectric_effect_css[i].first ),
                            false );

    data_container.setSubshellPhotoelectricCrossSection(
                                    subshell_photoelectric_effect_css[i].first,
                                    cross_section );
    data_container.setSubshellPhotoelectricCrossSectionThresholdEnergyIndex(
                                    subshell_photoelectric_effect_css[i].first,
                                    threshold );

    FRENSIE_LOG_NOTIFICATION( Utility::BoldGreen( "done." ) );
  }

  for( unsigned i = 0; i < impulse_approx_incoherent_cs_evaluators.size(); ++i)
  {
    FRENSIE_LOG_PARTIAL_NOTIFICATION( " Setting " <<
                                      Utility::Italicized( "subshell " )
                                      << Utility::Italicized(Data::convertENDFDesignatorToSubshellEnum( impulse_approx_incoherent_cs_evaluators[i].first ) )
                                      << Utility::Italicized(" impulse approx incoherent " )
                                      << "cross section ... " );
    FRENSIE_FLUSH_ALL_LOGS();

    this->createCrossSectionOnUnionEnergyGrid(
                             union_energy_grid,
                             impulse_approx_incoherent_cs_evaluators[i].second,
                             cross_section,
                             threshold );

    data_container.setImpulseApproxSubshellIncoherentCrossSection(
                              impulse_approx_incoherent_cs_evaluators[i].first,
                              cross_section );
    data_container.setImpulseApproxSubshellIncoherentCrossSectionThresholdEnergyIndex(
                              impulse_approx_incoherent_cs_evaluators[i].first,
                              threshold );

    FRENSIE_LOG_NOTIFICATION( Utility::BoldGreen( "done." ) );
  }

  FRENSIE_LOG_PARTIAL_NOTIFICATION( " Setting the " <<
                                    Utility::Italicized( "total photoelectric " )
                                    << "cross section ... " );
  FRENSIE_FLUSH_ALL_LOGS();

  this->calculateTotalPhotoelectricCrossSection();

  FRENSIE_LOG_NOTIFICATION( Utility::BoldGreen( "done." ) );

  FRENSIE_LOG_PARTIAL_NOTIFICATION( " Setting the " <<
                                    Utility::Italicized( "impulse approx total incoherent " )
                                    << "cross section ... " );
  FRENSIE_FLUSH_ALL_LOGS();

  this->calculateImpulseApproxTotalIncoherentCrossSection();

  FRENSIE_LOG_NOTIFICATION( Utility::BoldGreen( "done." ) );

  FRENSIE_LOG_PARTIAL_NOTIFICATION( " Setting the " <<
                                    Utility::Italicized( "Waller-Hartree total " )
                                    << "cross section ... " );
  FRENSIE_FLUSH_ALL_LOGS();

  this->calculatePhotonTotalCrossSection( true );

  FRENSIE_LOG_NOTIFICATION( Utility::BoldGreen( "done." ) );

  FRENSIE_LOG_PARTIAL_NOTIFICATION( " Setting the " <<
                                    Utility::Italicized( "impulse approx total " )
                                    << "cross section ... " );

  this->calculatePhotonTotalCrossSection( false );

  FRENSIE_LOG_NOTIFICATION( Utility::BoldGreen( "done." ) );
}

// Create the subshell impulse approx incoherent cross section evaluators
void ENDLElectronPhotonRelaxationDataGenerator::createSubshellImpulseApproxIncoherentCrossSectionEvaluators(
     std::vector<std::pair<unsigned,std::shared_ptr<const MonteCarlo::SubshellIncoherentPhotonScatteringDistribution> > >& evaluators ) const
{
  const Data::ElectronPhotonRelaxationDataContainer& data_container =
    this->getDataContainer();

  const std::set<unsigned>& subshells = data_container.getSubshells();

  evaluators.resize( subshells.size() );

  int i = 0;

  for( auto&& subshell : subshells )
  {
    evaluators[i].first = subshell;

    const std::vector<double>& momentum_grid =
      data_container.getOccupationNumberMomentumGrid( subshell );

    const std::vector<double>& occupation_number =
      data_container.getOccupationNumber( subshell );

    std::shared_ptr<const Utility::UnivariateDistribution> occupation_number_dist(
       new Utility::TabularDistribution<Utility::LinLin>( momentum_grid,
                                                          occupation_number ) );

    evaluators[i].second.reset( new MonteCarlo::SubshellIncoherentPhotonScatteringDistribution(
                   Data::convertENDFDesignatorToSubshellEnum( subshell ),
                   data_container.getSubshellOccupancy( subshell ),
                   data_container.getSubshellBindingEnergy( subshell ),
                   occupation_number_dist ) );

    ++i;
  }
}

// Set the electron data
void ENDLElectronPhotonRelaxationDataGenerator::setElectronData()
{
//---------------------------------------------------------------------------//
// Set Elastic Data
//---------------------------------------------------------------------------//
  FRENSIE_LOG_PARTIAL_NOTIFICATION( " Setting the " <<
                                    Utility::Italicized( "elastic cutoff " )
                                    << "data ... " );
  FRENSIE_FLUSH_ALL_LOGS();

  Data::ElectronPhotonRelaxationVolatileDataContainer& data_container =
    this->getVolatileDataContainer();

  // Set the cutoff elastic scattering interpolation policy
  data_container.setCutoffElasticInterpPolicy( "Lin-Lin" );

  // Set elastic angular distribution
  std::map<double,std::vector<double> > elastic_pdf, elastic_angle;

  std::vector<double> angular_energy_grid(
    d_endl_data_container->getCutoffElasticAngularEnergyGrid() );

  /* ENDL gives the angular distribution in terms of the change in angle cosine
   * (delta_angle_cosine = 1 - angle_cosine). For the native format it needs to
   * be in terms of angle_cosine. This for loop flips the distribution and
   * changes the variables to angle_cosine.
   */
  std::vector<double>::iterator energy = angular_energy_grid.begin();
  std::vector<double>::iterator end_energy =
    Utility::Search::binaryUpperBound( energy,
                                       angular_energy_grid.end(),
                                       this->getMaxElectronEnergy() );
  end_energy++;

  for ( energy; energy != end_energy; ++energy )
  {
    calculateElasticAngleCosine(
        d_endl_data_container->getCutoffElasticAnglesAtEnergy( *energy ),
        d_endl_data_container->getCutoffElasticPDFAtEnergy( *energy ),
        elastic_angle[*energy],
        elastic_pdf[*energy] );
  }

  end_energy--;
  if( *end_energy != this->getMaxElectronEnergy() )
  {
    double max_cutoff_angle_cosine = 1.0;

      // Get the angular grid and pdf at the max energy
      if( this->getElectronTwoDGridPolicy() == MonteCarlo::DIRECT_GRID ||
          this->getElectronTwoDGridPolicy() == MonteCarlo::UNIT_BASE_GRID )
      {
        if ( this->getElectronTwoDInterpPolicy() == MonteCarlo::LOGLOGLOG_INTERPOLATION )
        {
          MonteCarlo::ElasticElectronScatteringDistributionNativeFactory::getAngularGridAndPDF<Utility::LogNudgedLogCosLog,Utility::Direct>(
            elastic_angle[this->getMaxElectronEnergy()],
            elastic_pdf[this->getMaxElectronEnergy()] ,
            elastic_angle,
            elastic_pdf,
            this->getMaxElectronEnergy(),
            max_cutoff_angle_cosine,
            this->getTabularEvaluationTolerance() );
        }
        else if ( this->getElectronTwoDInterpPolicy() == MonteCarlo::LINLINLIN_INTERPOLATION )
        {
          MonteCarlo::ElasticElectronScatteringDistributionNativeFactory::getAngularGridAndPDF<Utility::LinLinLin,Utility::Direct>(
            elastic_angle[this->getMaxElectronEnergy()],
            elastic_pdf[this->getMaxElectronEnergy()] ,
            elastic_angle,
            elastic_pdf,
            this->getMaxElectronEnergy(),
            max_cutoff_angle_cosine,
            this->getTabularEvaluationTolerance() );
        }
        else if ( this->getElectronTwoDInterpPolicy() == MonteCarlo::LINLINLOG_INTERPOLATION )
        {
          MonteCarlo::ElasticElectronScatteringDistributionNativeFactory::getAngularGridAndPDF<Utility::LinLinLog,Utility::Direct>(
            elastic_angle[this->getMaxElectronEnergy()],
            elastic_pdf[this->getMaxElectronEnergy()] ,
            elastic_angle,
            elastic_pdf,
            this->getMaxElectronEnergy(),
            max_cutoff_angle_cosine,
            this->getTabularEvaluationTolerance() );
        }
        else
        {
          THROW_EXCEPTION( std::runtime_error,
                           "Error: the desired 2D interpolation policy " <<
                           this->getElectronTwoDInterpPolicy() <<
                           " is currently not supported!" );
        }
      }
      else if( this->getElectronTwoDGridPolicy() == MonteCarlo::CORRELATED_GRID ||
               this->getElectronTwoDGridPolicy() == MonteCarlo::UNIT_BASE_CORRELATED_GRID )
      {
        if ( this->getElectronTwoDInterpPolicy() == MonteCarlo::LOGLOGLOG_INTERPOLATION )
        {
          MonteCarlo::ElasticElectronScatteringDistributionNativeFactory::getAngularGridAndPDF<Utility::LogNudgedLogCosLog,Utility::Correlated>(
            elastic_angle[this->getMaxElectronEnergy()],
            elastic_pdf[this->getMaxElectronEnergy()] ,
            elastic_angle,
            elastic_pdf,
            this->getMaxElectronEnergy(),
            max_cutoff_angle_cosine,
            this->getTabularEvaluationTolerance() );
        }
        else if ( this->getElectronTwoDInterpPolicy() == MonteCarlo::LINLINLIN_INTERPOLATION )
        {
          MonteCarlo::ElasticElectronScatteringDistributionNativeFactory::getAngularGridAndPDF<Utility::LinLinLin,Utility::Correlated>(
            elastic_angle[this->getMaxElectronEnergy()],
            elastic_pdf[this->getMaxElectronEnergy()] ,
            elastic_angle,
            elastic_pdf,
            this->getMaxElectronEnergy(),
            max_cutoff_angle_cosine,
            this->getTabularEvaluationTolerance() );
        }
        else if ( this->getElectronTwoDInterpPolicy() == MonteCarlo::LINLINLOG_INTERPOLATION )
        {
          MonteCarlo::ElasticElectronScatteringDistributionNativeFactory::getAngularGridAndPDF<Utility::LinLinLog,Utility::Correlated >(
            elastic_angle[this->getMaxElectronEnergy()],
            elastic_pdf[this->getMaxElectronEnergy()] ,
            elastic_angle,
            elastic_pdf,
            this->getMaxElectronEnergy(),
            max_cutoff_angle_cosine,
            this->getTabularEvaluationTolerance() );
        }
        else
        {
          THROW_EXCEPTION( std::runtime_error,
                           "Error: the desired 2D interpolation policy " <<
                           this->getElectronTwoDInterpPolicy() <<
                           " is currently not supported!" );
        }
      }
      else
      {
        THROW_EXCEPTION( std::runtime_error,
                         "Error: the desired 2D grid policy " <<
                         this->getElectronTwoDGridPolicy() <<
                         " is currently not supported!" );
      }

    elastic_angle.erase( *end_energy );
    elastic_pdf.erase( *end_energy );
  }

  // Set new angular energy grid
  while( angular_energy_grid.back() >= this->getMaxElectronEnergy() )
  {
    angular_energy_grid.pop_back();
  }
  angular_energy_grid.push_back( this->getMaxElectronEnergy() );

  // Set the elastic cutoff data
  data_container.setElasticAngularEnergyGrid( angular_energy_grid );
  data_container.setCutoffElasticPDF( elastic_pdf );
  data_container.setCutoffElasticAngles( elastic_angle );

  FRENSIE_LOG_NOTIFICATION( Utility::BoldGreen( "done." ) );

//---------------------------------------------------------------------------//
// Set Electron Cross Section Data Data
//---------------------------------------------------------------------------//
/*! \details The cross section data is needed for calculating the
 *  moment preserving data and must be set first.
 */
  FRENSIE_LOG_NOTIFICATION( " Setting the electron cross section data:" );
  FRENSIE_FLUSH_ALL_LOGS();

  this->setElectronCrossSectionsData();

  FRENSIE_LOG_PARTIAL_NOTIFICATION( " Setting the " <<
                                    Utility::Italicized( "elastic moment preserving " )
                                    << "data ... " );
  FRENSIE_FLUSH_ALL_LOGS();

  // Check if moment preserving data can be generated
  if ( this->getCutoffAngleCosine() > 0.999999 ||
       this->getNumberOfMomentPreservingAngles() < 1 )
  {
    FRENSIE_LOG_NOTIFICATION( Utility::BoldYellow( "done." ) );

    if( this->getCutoffAngleCosine() > 0.999999 )
    {
      FRENSIE_LOG_NOTIFICATION( Utility::BoldCyan( "  Note: " )
                                << "Moment preserving data was not generated "
                                << "because the cutoff angle cosine is "
                                << "greater than 0.999999 (mu="
                                << this->getCutoffAngleCosine() << ")." );
    }

    if( this->getNumberOfMomentPreservingAngles() < 1 )
    {
      FRENSIE_LOG_NOTIFICATION( Utility::BoldCyan( "  Note: " )
                                << "Moment preserving data was not "
                                << "generated because the number of moment "
                                << "preserving angles is less than 1." );
    }
  }
  // Set the moment preserving data
  else
  {
    this->setMomentPreservingData( angular_energy_grid,
                                   this->getTabularEvaluationTolerance(),
                                   this->getElectronTwoDInterpPolicy() );

    FRENSIE_LOG_NOTIFICATION( Utility::BoldGreen( "done." ) );
  }

//---------------------------------------------------------------------------//
// Set Electroionization Data
//---------------------------------------------------------------------------//
  FRENSIE_LOG_PARTIAL_NOTIFICATION( " Setting the " <<
                                    Utility::Italicized( "electroionization " )
                                    << "data ... " );
  FRENSIE_FLUSH_ALL_LOGS();

  std::set<unsigned>::iterator shell = data_container.getSubshells().begin();

  // Set the interpolation policy
  data_container.setElectroionizationRecoilInterpPolicy( "Lin-Lin" );

  // Loop through electroionization data for every subshell
  for ( shell; shell != data_container.getSubshells().end(); ++shell )
  {
    data_container.setElectroionizationEnergyGrid(
        *shell,
        d_endl_data_container->getElectroionizationRecoilEnergyGrid( *shell ) );

    data_container.setElectroionizationRecoilEnergy(
        *shell,
        d_endl_data_container->getElectroionizationRecoilEnergy( *shell ) );

    data_container.setElectroionizationRecoilPDF(
        *shell,
        d_endl_data_container->getElectroionizationRecoilPDF( *shell ) );
  }

  FRENSIE_LOG_NOTIFICATION( Utility::BoldGreen( "done." ) );

//---------------------------------------------------------------------------//
// Set Bremsstrahlung Data
//---------------------------------------------------------------------------//
  FRENSIE_LOG_PARTIAL_NOTIFICATION( " Setting the " <<
                                    Utility::Italicized( "bremsstrahlung " )
                                    << "data ... " );
  FRENSIE_FLUSH_ALL_LOGS();

  // Set the interpolation policy
  data_container.setBremsstrahlungPhotonInterpPolicy( "Lin-Lin" );

  data_container.setBremsstrahlungEnergyGrid(
    d_endl_data_container->getBremsstrahlungPhotonEnergyGrid() );

  data_container.setBremsstrahlungPhotonEnergy(
    d_endl_data_container->getBremsstrahlungPhotonEnergy() );

  data_container.setBremsstrahlungPhotonPDF(
    d_endl_data_container->getBremsstrahlungPhotonPDF() );

  FRENSIE_LOG_NOTIFICATION( Utility::BoldGreen( "done." ) );

//---------------------------------------------------------------------------//
// Set Atomic Excitation Data
//---------------------------------------------------------------------------//
  FRENSIE_LOG_PARTIAL_NOTIFICATION( " Setting the " <<
                                    Utility::Italicized( "atomic excitation " )
                                    << "data ... " );
  FRENSIE_FLUSH_ALL_LOGS();

  // Set the interpolation policy
  data_container.setAtomicExcitationEnergyLossInterpPolicy( "Lin-Lin" );

  // Set atomic excitation energy loss
  data_container.setAtomicExcitationEnergyGrid(
    d_endl_data_container->getAtomicExcitationEnergyGrid() );
  data_container.setAtomicExcitationEnergyLoss(
    d_endl_data_container->getAtomicExcitationEnergyLoss() );

  FRENSIE_LOG_NOTIFICATION( Utility::BoldGreen( "done." ) );
}

// Set the electron cross section union energy grid
void ENDLElectronPhotonRelaxationDataGenerator::setElectronCrossSectionsData()
{
  Data::ElectronPhotonRelaxationVolatileDataContainer& data_container =
    this->getVolatileDataContainer();

  // cross sections in the file
  std::shared_ptr<const Utility::UnivariateDistribution>
        bremsstrahlung_cross_section, atomic_excitation_cross_section,
        cutoff_elastic_cross_section, total_elastic_cross_section;

  // Initialize union energy grid
  std::list<double> union_energy_grid;
  this->initializeElectronUnionEnergyGrid( union_energy_grid );

//---------------------------------------------------------------------------//
// Get Elastic Data Cross Section Data
//---------------------------------------------------------------------------//
  // Get cutoff elastic cross section to union energy grid
  std::vector<double> raw_energy_grid =
    d_endl_data_container->getElasticEnergyGrid();

  cutoff_elastic_cross_section.reset(
    new Utility::TabularDistribution<Utility::LogLog>(
    raw_energy_grid,
    d_endl_data_container->getCutoffElasticCrossSection() ) );

  // merge raw energy grid with the union energy grid
  this->mergeElectronUnionEnergyGrid( raw_energy_grid, union_energy_grid );

  // Calculate the electron total elastic cross section
  this->calculateElectronTotalElasticCrossSection( total_elastic_cross_section,
                                                   raw_energy_grid );

//---------------------------------------------------------------------------//
// Get Electroionization Data Cross Section Data
//---------------------------------------------------------------------------//

  std::set<unsigned>::iterator shell = data_container.getSubshells().begin();

  std::vector<std::pair<unsigned,std::shared_ptr<const Utility::UnivariateDistribution> > >
    electroionization_cross_section( data_container.getSubshells().size() );

  unsigned i = 0;
  // Loop through electroionization data for every subshell
  for ( shell; shell != data_container.getSubshells().end(); ++shell )
  {
    // Get the raw energy grid
    raw_energy_grid =
        d_endl_data_container->getElectroionizationCrossSectionEnergyGrid(*shell);

    // merge raw energy grid with the union energy grid
    this->mergeElectronUnionEnergyGrid( raw_energy_grid, union_energy_grid );

    /*! \details There is conflicting documentation on the proper interpolation
     *  of the electroionization cross section data. The endl data file interp flag
     *  specifies lin-lin, but the documentation write-up says to use log-log
     *  interpolation on all cross sections. It was decided to match MCNP which
     *  uses log-log interpolation for electroionization.
     */
    this->extractElectronCrossSection<Utility::LogLog>(
                raw_energy_grid,
                d_endl_data_container->getElectroionizationCrossSection(*shell),
                electroionization_cross_section[i].second );

    // Set the shell identifier
    electroionization_cross_section[i].first = *shell;
    ++i;
  }

//---------------------------------------------------------------------------//
// Get Bremsstrahlung Cross Section Data
//---------------------------------------------------------------------------//

  raw_energy_grid =
    d_endl_data_container->getBremsstrahlungCrossSectionEnergyGrid();

  bremsstrahlung_cross_section.reset(
    new Utility::TabularDistribution<Utility::LogLog>(
    raw_energy_grid,
    d_endl_data_container->getBremsstrahlungCrossSection() ) );

  // merge raw energy grid with the union energy grid
  this->mergeElectronUnionEnergyGrid( raw_energy_grid, union_energy_grid );

//---------------------------------------------------------------------------//
// Get Atomic Excitation Data Cross Section Data
//---------------------------------------------------------------------------//

  raw_energy_grid = d_endl_data_container->getAtomicExcitationEnergyGrid();
  double true_atomic_excitation_threshold_energy = raw_energy_grid[0];
  this->extractElectronCrossSection<Utility::LogLog>(
              raw_energy_grid,
              d_endl_data_container->getAtomicExcitationCrossSection(),
              atomic_excitation_cross_section );

  // merge raw energy grid with the union energy grid
  this->mergeElectronUnionEnergyGrid( raw_energy_grid, union_energy_grid );

//---------------------------------------------------------------------------//
// Create union energy grid and calculate cross sections
//---------------------------------------------------------------------------//

  // Create the union energy grid
  FRENSIE_LOG_PARTIAL_NOTIFICATION( "   Creating " <<
                                    Utility::Italicized( "union energy grid " ) );
  FRENSIE_FLUSH_ALL_LOGS();

  // Calculate the union energy grid
  boost::function<double (double pz)> grid_function =
    boost::bind( &Utility::UnivariateDistribution::evaluate,
                 boost::cref( *cutoff_elastic_cross_section ),
                 _1 );

  try{
    this->getDefaultElectronGridGenerator().generateInPlace( union_energy_grid,
                                                             grid_function );
  }
  EXCEPTION_CATCH_RETHROW( std::runtime_error,
                           "Could not generate an optimized electron "
                           "energy grid for the cutoff elastic cross section "
                           "with the provided convergence parameters!" );

  FRENSIE_LOG_PARTIAL_NOTIFICATION( "." );
  FRENSIE_FLUSH_ALL_LOGS();

  grid_function = boost::bind(
        &Utility::UnivariateDistribution::evaluate,
        boost::cref( *total_elastic_cross_section ),
        _1 );

  try{
    this->getDefaultElectronGridGenerator().generateInPlace( union_energy_grid,
                                                             grid_function );
  }
  EXCEPTION_CATCH_RETHROW( std::runtime_error,
                           "Could not generate an optimized electron "
                           "energy grid for the total elastic cross section "
                           "with the provided convergence parameters!" );

  FRENSIE_LOG_PARTIAL_NOTIFICATION( "." );
  FRENSIE_FLUSH_ALL_LOGS();

  grid_function = boost::bind( &Utility::UnivariateDistribution::evaluate,
                               boost::cref( *bremsstrahlung_cross_section ),
                               _1 );

  try{
    this->getDefaultElectronGridGenerator().generateInPlace( union_energy_grid,
                                                             grid_function );
  }
  EXCEPTION_CATCH_RETHROW( std::runtime_error,
                           "Could not generate an optimized electron "
                           "energy grid for the bremsstrahlung cross section "
                           "with the provided convergence parameters!" );

  FRENSIE_LOG_PARTIAL_NOTIFICATION( "." );
  FRENSIE_FLUSH_ALL_LOGS();

  grid_function = boost::bind( &Utility::UnivariateDistribution::evaluate,
                               boost::cref( *atomic_excitation_cross_section ),
                               _1 );

  try{
    this->getDefaultElectronGridGenerator().refineInPlace(
                union_energy_grid,
                grid_function,
                atomic_excitation_cross_section->getLowerBoundOfIndepVar(),
                atomic_excitation_cross_section->getUpperBoundOfIndepVar() );
  }
  EXCEPTION_CATCH_RETHROW( std::runtime_error,
                           "Could not generate an optimized electron "
                           "energy grid for the atomic excitation cross "
                           "section with the provided convergence "
                           "parameters!" );

  FRENSIE_LOG_PARTIAL_NOTIFICATION( "." );
  FRENSIE_FLUSH_ALL_LOGS();

  for( unsigned i = 0; i < electroionization_cross_section.size(); ++i )
  {
    grid_function = boost::bind(
                   &Utility::UnivariateDistribution::evaluate,
                   boost::cref( *electroionization_cross_section[i].second ),
                   _1 );

    try{
      this->getDefaultElectronGridGenerator().refineInPlace(
        union_energy_grid,
        grid_function,
        electroionization_cross_section[i].second->getLowerBoundOfIndepVar(),
        electroionization_cross_section[i].second->getUpperBoundOfIndepVar() );
    }
    EXCEPTION_CATCH_RETHROW(
                          std::runtime_error,
                          "Could not generate an optimized electron "
                          "energy grid for the subshell "
                          << Data::convertENDFDesignatorToSubshellEnum(
                                  electroionization_cross_section[i].first ) <<
                          " electroionization cross section with the "
                          "provided convergence parameters!" );

    FRENSIE_LOG_PARTIAL_NOTIFICATION( "." );
    FRENSIE_FLUSH_ALL_LOGS();
  }

  FRENSIE_LOG_NOTIFICATION( Utility::BoldGreen( " done." ) );

  // Set the union energy grid
  std::vector<double> energy_grid;
  energy_grid.assign( union_energy_grid.begin(),
                      union_energy_grid.end() );

  data_container.setElectronEnergyGrid( energy_grid );

  // Set the electron cross section interpolation policy
  data_container.setElectronCrossSectionInterpPolicy( "Log-Log" );

  // Create and set the cross sections
  std::vector<double> cross_section;
  unsigned threshold;

  {
    // Set Elastic cross section data
    std::vector<double> cutoff_cross_section, total_cross_section;

    FRENSIE_LOG_PARTIAL_NOTIFICATION( "   Setting the " <<
                                      Utility::Italicized( "cutoff elastic " )
                                      << "cross section ... " );
    FRENSIE_FLUSH_ALL_LOGS();

    this->createCrossSectionOnUnionEnergyGrid(
                     union_energy_grid,
                     cutoff_elastic_cross_section,
                     cutoff_cross_section,
                     threshold,
                     cutoff_elastic_cross_section->getLowerBoundOfIndepVar(),
                     false );

    data_container.setCutoffElasticCrossSection( cutoff_cross_section );
    data_container.setCutoffElasticCrossSectionThresholdEnergyIndex( threshold );

    FRENSIE_LOG_NOTIFICATION( Utility::BoldGreen( "done." ) );

    FRENSIE_LOG_PARTIAL_NOTIFICATION( "   Setting the " <<
                                      Utility::Italicized( "screened Rutherford elastic " )
                                      << "cross section ... " );
    FRENSIE_FLUSH_ALL_LOGS();

    this->createCrossSectionOnUnionEnergyGrid(
                      union_energy_grid,
                      total_elastic_cross_section,
                      total_cross_section,
                      threshold,
                      total_elastic_cross_section->getLowerBoundOfIndepVar(),
                      false );

    data_container.setTotalElasticCrossSection( total_cross_section );
    data_container.setTotalElasticCrossSectionThresholdEnergyIndex( threshold );

    std::vector<double> raw_cross_section( total_cross_section.size() );
    for( int i = 0; i < total_cross_section.size(); ++i )
    {
      raw_cross_section[i] = total_cross_section[i] - cutoff_cross_section[i];

      // Calculate the relative difference between the total and cutoff cross sections
      double relative_difference = raw_cross_section[i]/total_cross_section[i];

      // Check for roundoff error and reduce to zero if needed
      if( relative_difference < 1.0e-6 )
      {
        raw_cross_section[i] = 0.0;

        // Update threshold index
        threshold = i+1;
      }
    }

    std::vector<double>::iterator start = raw_cross_section.begin();
    std::advance( start, threshold );

    cross_section.assign( start, raw_cross_section.end() );

    data_container.setScreenedRutherfordElasticCrossSection( cross_section );
    data_container.setScreenedRutherfordElasticCrossSectionThresholdEnergyIndex(
                                                                                threshold );
    FRENSIE_LOG_NOTIFICATION( Utility::BoldGreen( "done." ) );
  }

  FRENSIE_LOG_PARTIAL_NOTIFICATION( "   Setting the " <<
                                    Utility::Italicized( "bremsstrahlung " )
                                    << "cross section ... " );
  FRENSIE_FLUSH_ALL_LOGS();

  this->createCrossSectionOnUnionEnergyGrid(
                     union_energy_grid,
                     bremsstrahlung_cross_section,
                     cross_section,
                     threshold,
                     bremsstrahlung_cross_section->getLowerBoundOfIndepVar(),
                     false );

  data_container.setBremsstrahlungCrossSection( cross_section );
  data_container.setBremsstrahlungCrossSectionThresholdEnergyIndex( threshold );

  FRENSIE_LOG_NOTIFICATION( Utility::BoldGreen( "done." ) );

  FRENSIE_LOG_PARTIAL_NOTIFICATION( "   Setting the " <<
                                    Utility::Italicized( "atomic excitation " )
                                    << "cross section ... " );
  FRENSIE_FLUSH_ALL_LOGS();

  this->createCrossSectionOnUnionEnergyGrid(
                  union_energy_grid,
                  atomic_excitation_cross_section,
                  cross_section,
                  threshold,
                  true_atomic_excitation_threshold_energy,
                  true );

  data_container.setAtomicExcitationCrossSection( cross_section );
  data_container.setAtomicExcitationCrossSectionThresholdEnergyIndex(
                                                                   threshold );

  FRENSIE_LOG_NOTIFICATION( Utility::BoldGreen( "done." ) );

  for( unsigned i = 0; i < electroionization_cross_section.size(); ++i )
  {
    FRENSIE_LOG_PARTIAL_NOTIFICATION( "   Setting " <<
                                      Utility::Italicized( "subshell " )
                                      << Utility::Italicized(Data::convertENDFDesignatorToSubshellEnum( electroionization_cross_section[i].first ) )
                                      << Utility::Italicized( " electroionization " )
                                      << "cross section ... " );
    FRENSIE_FLUSH_ALL_LOGS();

    this->createCrossSectionOnUnionEnergyGrid(
                    union_energy_grid,
                    electroionization_cross_section[i].second,
                    cross_section,
                    threshold,
                    data_container.getSubshellBindingEnergy(
                      electroionization_cross_section[i].first ),
                    true );

    data_container.setElectroionizationCrossSection(
                                    electroionization_cross_section[i].first,
                                    cross_section );
    data_container.setElectroionizationCrossSectionThresholdEnergyIndex(
                                    electroionization_cross_section[i].first,
                                    threshold );

    FRENSIE_LOG_NOTIFICATION( Utility::BoldGreen( "done." ) );
  }

  FRENSIE_LOG_PARTIAL_NOTIFICATION( "   Setting the " <<
                                    Utility::Italicized( "Total Electron " )
                                    << "cross section ... " );
  FRENSIE_FLUSH_ALL_LOGS();

  this->calculateElectronTotalCrossSection();

  FRENSIE_LOG_NOTIFICATION( Utility::BoldGreen( "done." ) );
}

// Set the moment preserving data
void ENDLElectronPhotonRelaxationDataGenerator::setMomentPreservingData(
                         const std::vector<double>& angular_energy_grid,
                         const double tabular_evaluation_tol,
                         const MonteCarlo::TwoDInterpolationType two_d_interp )
{
  // Make sure the tolerance is valid
  testPrecondition( tabular_evaluation_tol > 0.0 );
  testPrecondition( tabular_evaluation_tol < 1.0 );

  // Create the coupled elastic distribution (combined Cutoff and Screened Rutherford)
  std::shared_ptr<const MonteCarlo::CoupledElasticElectronScatteringDistribution>
        coupled_distribution;

  // Get the cutoff and total elastic cross sections and the energy grid
  std::shared_ptr<std::vector<double> > cutoff_cross_section(
       new std::vector<double>( this->getDataContainer().getCutoffElasticCrossSection() ) );
  std::shared_ptr<std::vector<double> > total_cross_section(
       new std::vector<double>( this->getDataContainer().getTotalElasticCrossSection() ) );
  std::shared_ptr<std::vector<double> > energy_grid(
       new std::vector<double>( this->getDataContainer().getElectronEnergyGrid() ) );

  if ( two_d_interp == MonteCarlo::LOGLOGLOG_INTERPOLATION )
  {
    MonteCarlo::ElasticElectronScatteringDistributionNativeFactory::createCoupledElasticDistribution<Utility::LogNudgedLogCosLog,Utility::Correlated>(
        coupled_distribution,
        cutoff_cross_section,
        total_cross_section,
        energy_grid,
        this->getDataContainer().getCutoffElasticAngles(),
        this->getDataContainer().getCutoffElasticPDF(),
        angular_energy_grid,
        this->getDataContainer().getAtomicNumber(),
        MonteCarlo::MODIFIED_TWO_D_UNION,
        tabular_evaluation_tol );
  }
  else if ( two_d_interp == MonteCarlo::LINLINLIN_INTERPOLATION )
  {
    MonteCarlo::ElasticElectronScatteringDistributionNativeFactory::createCoupledElasticDistribution<Utility::LinLinLin,Utility::Correlated>(
        coupled_distribution,
        cutoff_cross_section,
        total_cross_section,
        energy_grid,
        this->getDataContainer().getCutoffElasticAngles(),
        this->getDataContainer().getCutoffElasticPDF(),
        angular_energy_grid,
        this->getDataContainer().getAtomicNumber(),
        MonteCarlo::MODIFIED_TWO_D_UNION,
        tabular_evaluation_tol );
  }
  else if ( two_d_interp == MonteCarlo::LINLINLOG_INTERPOLATION )
  {
    MonteCarlo::ElasticElectronScatteringDistributionNativeFactory::createCoupledElasticDistribution<Utility::LinLinLog,Utility::Correlated>(
        coupled_distribution,
        cutoff_cross_section,
        total_cross_section,
        energy_grid,
        this->getDataContainer().getCutoffElasticAngles(),
        this->getDataContainer().getCutoffElasticPDF(),
        angular_energy_grid,
        this->getDataContainer().getAtomicNumber(),
        MonteCarlo::MODIFIED_TWO_D_UNION,
        tabular_evaluation_tol );
  }

  // Construct the hash-based grid searcher for this atom
  std::shared_ptr<Utility::HashBasedGridSearcher<double> > grid_searcher(
     new Utility::StandardHashBasedGridSearcher<std::vector<double>, false>(
             energy_grid,
             100u ) );

  // Create the elastic traits
  std::shared_ptr<MonteCarlo::ElasticElectronTraits> elastic_traits(
    new MonteCarlo::ElasticElectronTraits( this->getDataContainer().getAtomicNumber() ) );

  // Create the moment evaluator of the elastic scattering distribution
  std::shared_ptr<DataGen::ElasticElectronMomentsEvaluator> moments_evaluator;
  moments_evaluator.reset(
    new DataGen::ElasticElectronMomentsEvaluator(
        this->getDataContainer().getCutoffElasticAngles(),
        energy_grid,
        grid_searcher,
        cutoff_cross_section,
        total_cross_section,
        this->getDataContainer().getScreenedRutherfordElasticCrossSectionThresholdEnergyIndex(),
        coupled_distribution,
        elastic_traits,
        this->getDataContainer().getCutoffAngleCosine() ) );

  // Moment preserving discrete angles and weights
  std::vector<double> discrete_angles, weights;

  // weights for a discrete angle cosine = 1
  std::vector<double> cross_section_reduction( angular_energy_grid.size() );

  // iterate through all angular energy bins
  for ( unsigned i = 0; i < angular_energy_grid.size(); ++i )
  {
    ENDLElectronPhotonRelaxationDataGenerator::calculateDiscreteAnglesAndWeights(
        moments_evaluator,
        angular_energy_grid[i],
        this->getDataContainer().getNumberOfMomentPreservingAngles(),
        discrete_angles,
        weights,
        cross_section_reduction[i] );

    this->getVolatileDataContainer().setMomentPreservingElasticDiscreteAngles(
                                                        angular_energy_grid[i],
                                                        discrete_angles );
    this->getVolatileDataContainer().setMomentPreservingElasticWeights(
                                                        angular_energy_grid[i],
                                                        weights );
  }

  // Set the cross section reduction
  this->getVolatileDataContainer().setMomentPreservingCrossSectionReduction( cross_section_reduction );
}

// Initialize the photon union energy grid
void ENDLElectronPhotonRelaxationDataGenerator::initializePhotonUnionEnergyGrid(
                                   std::list<double>& union_energy_grid ) const
{
  const Data::ElectronPhotonRelaxationDataContainer& data_container =
    this->getDataContainer();

  // Add the min photon energy to the union energy grid
  union_energy_grid.push_back( this->getMinPhotonEnergy() );

  // Add the binding energies
  this->addBindingEnergiesToUnionEnergyGrid( this->getMinPhotonEnergy(),
                                             this->getMaxPhotonEnergy(),
                                             true,
                                             union_energy_grid );

  // Add the pair production threshold
  const double pp_threshold =
    2*Utility::PhysicalConstants::electron_rest_mass_energy;

  if( pp_threshold > this->getMinPhotonEnergy() &&
      pp_threshold < this->getMaxPhotonEnergy() )
  {
    union_energy_grid.push_back( pp_threshold );
  }

  const double nudged_pp_threshold =
    pp_threshold*this->getPhotonThresholdEnergyNudgeFactor();

  if( nudged_pp_threshold > this->getMinPhotonEnergy() &&
      nudged_pp_threshold < this->getMinPhotonEnergy() )
  {
    union_energy_grid.push_back( nudged_pp_threshold );
  }

  // Add the triplet production threshold
  const double tp_threshold =
    4*Utility::PhysicalConstants::electron_rest_mass_energy;

  if( tp_threshold > this->getMinPhotonEnergy() &&
      tp_threshold < this->getMaxPhotonEnergy() )
  {
    union_energy_grid.push_back( tp_threshold );
  }

  const double nudged_tp_threshold =
    tp_threshold*this->getPhotonThresholdEnergyNudgeFactor();

  if( nudged_tp_threshold > this->getMinPhotonEnergy() &&
      nudged_tp_threshold < this->getMaxPhotonEnergy() )
  {
    union_energy_grid.push_back( nudged_tp_threshold );
  }

  // Add the max photon energy
  union_energy_grid.push_back( this->getMaxPhotonEnergy() );

  // Sort the union energy grid
  union_energy_grid.sort();
}

// Initialize the electron union energy grid
void ENDLElectronPhotonRelaxationDataGenerator::initializeElectronUnionEnergyGrid(
                                   std::list<double>& union_energy_grid ) const
{
  // Add the min electron energy to the union energy grid
  union_energy_grid.push_back( this->getMinElectronEnergy() );

  // Add the binding energies
  this->addBindingEnergiesToUnionEnergyGrid( this->getMinElectronEnergy(),
                                             this->getMaxElectronEnergy(),
                                             true,
                                             union_energy_grid );

  // Add the max electron energy
  union_energy_grid.push_back( this->getMaxElectronEnergy() );

  // Sort the union energy grid
  union_energy_grid.sort();
}

// Add binding energies to union energy grid
void ENDLElectronPhotonRelaxationDataGenerator::addBindingEnergiesToUnionEnergyGrid(
                                   const double min_energy,
                                   const double max_energy,
                                   const bool add_nudged_values,
                                   std::list<double>& union_energy_grid ) const
{
  const Data::ElectronPhotonRelaxationDataContainer& data_container =
    this->getDataContainer();

  const std::set<unsigned>& subshells = data_container.getSubshells();

  std::set<unsigned>::const_iterator subshell = subshells.begin();

  // Add the subshell binding energies
  while( subshell != subshells.end() )
  {
    double binding_energy =
      data_container.getSubshellBindingEnergy( *subshell );

    // The subshell photoelectric cross sections have a jump discontinuity at
    // the binding energy of the subshell. Ideally we would add in the binding
    // energy twice and evaluate it to zero at the first instance and the cross
    // section value at the second instance. We will use a very small shift
    // value to mimic this behavior since our grid generation algorithm is not
    // set up to handle this special case
    if( binding_energy - 1e-9 > min_energy &&
        binding_energy - 1e-9 < max_energy )
    {
      union_energy_grid.push_back( binding_energy - 1e-9 );
    }

    if( binding_energy > min_energy &&
        binding_energy < max_energy )
    {
      union_energy_grid.push_back( binding_energy );
    }

    // Since the subshell incoherent, pair production and triplet production
    // cross sections go to zero at the binding energy we will add another
    // grid point, which is the binding energy plus a small shift value.
    // Avoiding the zero value of the cross section has been shown to improve
    // grid generation times. The true threshold energy (binding energy) will
    // be added back into the grid for the subshell cross section at the end
    if( add_nudged_values )
    {
      const double nudged_binding_energy =
        binding_energy*this->getPhotonThresholdEnergyNudgeFactor();

      if( nudged_binding_energy > this->getMinPhotonEnergy() &&
          nudged_binding_energy < this->getMaxPhotonEnergy() )
      {
        union_energy_grid.push_back( nudged_binding_energy );
      }
    }

    ++subshell;
  }
}

// Merge a secondary energy grid with the electron union energy grid
void ENDLElectronPhotonRelaxationDataGenerator::mergeElectronUnionEnergyGrid(
     const std::vector<double>& energy_grid,
     std::list<double>& union_energy_grid ) const
{
  // Assign the new grid to the union grid
  union_energy_grid.insert( union_energy_grid.begin(),
                            energy_grid.begin(),
                            energy_grid.end() );

  // Sort the union energy grid
  union_energy_grid.sort();

  // Remove all energies less than the min
  while ( union_energy_grid.front() < this->getMinElectronEnergy() )
  {
    union_energy_grid.pop_front();
  }

  // Remove all energies greater than the max
  while ( union_energy_grid.back() > this->getMaxElectronEnergy() )
  {
    union_energy_grid.pop_back();
  }

  // Make sure the union energy grid values are unique
  union_energy_grid.unique();
}

// Create the cross section on the union energy grid
void ENDLElectronPhotonRelaxationDataGenerator::createCrossSectionOnUnionEnergyGrid(
   const std::list<double>& union_energy_grid,
   const std::shared_ptr<const Utility::UnivariateDistribution>& original_cross_section,
   std::vector<double>& cross_section,
   unsigned& threshold_index,
   const double true_threshold_energy,
   const bool zero_at_threshold ) const
{
  // Make the threshold region cross section (this step is not necessary if
  // the true threshold energy == the threshold energy present in the cross
  // section distribution)
  std::shared_ptr<const Utility::UnivariateDistribution>
    threshold_region_cross_section;

  // Get the cross sections current (possibly fictitious) threshold energy
  double current_threshold_energy =
    original_cross_section->getLowerBoundOfIndepVar();

  // Note: The cross section in the threshold region will either be a
  //       line from zero to the cross section at the current threshold
  //       energy (if the cross section goes to zero at its threshold - e.g.
  //       pair production) or a uniform distribution with a value equal to
  //       the cross section at the current threshold energy (if the cross
  //       section has a jump discontinuity at its threshold - e.g.
  //       photoelectric).
  {
    // Get the cross section at the current threshold energy
    double current_threshold_cross_section =
      original_cross_section->evaluate( current_threshold_energy );

    // Create a distribution for the range between the true threshold and
    // the current threshold
    if( true_threshold_energy < current_threshold_energy )
    {
      std::vector<double> threshold_region_energy_range( 2 );
      threshold_region_energy_range[0] = true_threshold_energy;
      threshold_region_energy_range[1] = current_threshold_energy;

      std::vector<double> threshold_region_cross_section_values( 2 );

      if( zero_at_threshold )
        threshold_region_cross_section_values[0] = 0.0;
      else
      {
        threshold_region_cross_section_values[0] =
          current_threshold_cross_section;
      }

      threshold_region_cross_section_values[1] =
        current_threshold_cross_section;

      threshold_region_cross_section.reset(
                          new Utility::TabularDistribution<Utility::LinLin>(
                                     threshold_region_energy_range,
                                     threshold_region_cross_section_values ) );
    }
  }

  // Evaluate the cross section on the union energy grid
  std::vector<double> raw_cross_section( union_energy_grid.size() );

  std::list<double>::const_iterator energy_grid_pt = union_energy_grid.begin();

  unsigned index = 0u;

  while( energy_grid_pt != union_energy_grid.end() )
  {
    if( *energy_grid_pt < true_threshold_energy )
      raw_cross_section[index] = 0.0;

    else if( *energy_grid_pt >= true_threshold_energy &&
             *energy_grid_pt < current_threshold_energy )
    {
      raw_cross_section[index] =
        threshold_region_cross_section->evaluate( *energy_grid_pt );
    }
    else
    {
      raw_cross_section[index] =
        original_cross_section->evaluate( *energy_grid_pt );
    }

    ++energy_grid_pt;
    ++index;
  }

  this->populateCrossSection( raw_cross_section,
                              cross_section,
                              threshold_index,
                              zero_at_threshold );
}

// Create the cross section on the union energy grid
void ENDLElectronPhotonRelaxationDataGenerator::createCrossSectionOnUnionEnergyGrid(
             const std::list<double>& union_energy_grid,
             const std::shared_ptr<const MonteCarlo::SubshellIncoherentPhotonScatteringDistribution>&
             original_cross_section,
             std::vector<double>& cross_section,
             unsigned& threshold_index ) const
{
  std::vector<double> raw_cross_section( union_energy_grid.size() );

  std::list<double>::const_iterator energy_grid_pt = union_energy_grid.begin();

  unsigned index = 0u;

  while( energy_grid_pt != union_energy_grid.end() )
  {
    raw_cross_section[index] =
      original_cross_section->evaluateIntegratedCrossSection(
                            *energy_grid_pt,
                            this->getSubshellIncoherentEvaluationTolerance() );

    ++energy_grid_pt;
    ++index;
  }

  // The subshell incoherent cross section is zero at the threshold energy
  this->populateCrossSection( raw_cross_section,
                              cross_section,
                              threshold_index,
                              true );
}

// Populate a cross section using the raw cross section
void ENDLElectronPhotonRelaxationDataGenerator::populateCrossSection(
                                  const std::vector<double>& raw_cross_section,
                                  std::vector<double>& cross_section,
                                  unsigned& threshold_index,
                                  const bool zero_at_threshold ) const
{
  // Find the threshold energy index
  std::vector<double>::const_iterator start =
    std::find_if( raw_cross_section.begin(),
                  raw_cross_section.end(),
                  notEqualZero );

  // If the value of the cross section at the threshold energy is zero, set
  // the threshold energy of the cross section to the last energy where it
  // was equal to zero.
  if( start != raw_cross_section.begin() && zero_at_threshold )
    --start;

  cross_section.assign( start, raw_cross_section.end() );

  threshold_index = std::distance( raw_cross_section.begin(), start );
}

// Calculate the total photoelectric cross section
void ENDLElectronPhotonRelaxationDataGenerator::calculateTotalPhotoelectricCrossSection()
{
  Data::ElectronPhotonRelaxationVolatileDataContainer& data_container =
    this->getVolatileDataContainer();

  const std::vector<double>& energy_grid =
    data_container.getPhotonEnergyGrid();

  std::vector<double> raw_cross_section( energy_grid.size(), 0.0 );

  const std::set<unsigned>& subshells = data_container.getSubshells();

  std::set<unsigned>::const_iterator subshell = subshells.begin();

  while( subshell != subshells.end() )
  {
    const std::vector<double>& subshell_photoelectric_cs =
      data_container.getSubshellPhotoelectricCrossSection( *subshell );

    unsigned start_index =
      energy_grid.size() - subshell_photoelectric_cs.size();

    for( unsigned i = 0; i < subshell_photoelectric_cs.size(); ++i )
      raw_cross_section[start_index+i] += subshell_photoelectric_cs[i];

    ++subshell;
  }

  std::vector<double> cross_section;
  unsigned threshold;

  // Note: The cross section has jump discontinuity at the lowest binding
  //       energy. We will therefore not shift the start index back one (as
  //       is done when the cross section is zero at the threshold).
  this->populateCrossSection( raw_cross_section,
                              cross_section,
                              threshold,
                              false );

  data_container.setPhotoelectricCrossSection( cross_section );
  data_container.setPhotoelectricCrossSectionThresholdEnergyIndex( threshold );
}

// Calculate the total impulse approx. incoherent cross section
void ENDLElectronPhotonRelaxationDataGenerator::calculateImpulseApproxTotalIncoherentCrossSection()
{
  Data::ElectronPhotonRelaxationVolatileDataContainer& data_container =
    this->getVolatileDataContainer();

  const std::vector<double>& energy_grid =
    data_container.getPhotonEnergyGrid();

  std::vector<double> raw_cross_section( energy_grid.size(), 0.0 );

  const std::set<unsigned>& subshells = data_container.getSubshells();

  std::set<unsigned>::const_iterator subshell = subshells.begin();

  while( subshell != subshells.end() )
  {
    const std::vector<double>& subshell_incoherent_cs =
      data_container.getImpulseApproxSubshellIncoherentCrossSection(*subshell);

    unsigned start_index = energy_grid.size() - subshell_incoherent_cs.size();

    for( unsigned i = 0; i < subshell_incoherent_cs.size(); ++i )
      raw_cross_section[start_index+i] += subshell_incoherent_cs[i];

    ++subshell;
  }

  std::vector<double> cross_section;
  unsigned threshold;

  // The value of the cross section at the threshold energy is zero. Set
  // the threshold energy of the cross section to the last energy where it
  // was equal to zero.
  this->populateCrossSection( raw_cross_section,
                              cross_section,
                              threshold,
                              true );

  data_container.setImpulseApproxIncoherentCrossSection( cross_section );
  data_container.setImpulseApproxIncoherentCrossSectionThresholdEnergyIndex(
                                                                   threshold );
}

// Calculate the photon total cross section
void ENDLElectronPhotonRelaxationDataGenerator::calculatePhotonTotalCrossSection(
                                  const bool use_waller_hartree_incoherent_cs )
{
  Data::ElectronPhotonRelaxationVolatileDataContainer& data_container =
    this->getVolatileDataContainer();

  const std::vector<double>& energy_grid =
    data_container.getPhotonEnergyGrid();

  std::vector<double> total_cross_section( energy_grid.size(), 0.0 );

  // Add the incoherent cs
  if( use_waller_hartree_incoherent_cs )
  {
    this->addCrossSectionToTotalCrossSection(
                       energy_grid,
                       data_container.getWallerHartreeIncoherentCrossSection(),
                       total_cross_section );
  }
  else
  {
    this->addCrossSectionToTotalCrossSection(
                       energy_grid,
                       data_container.getImpulseApproxIncoherentCrossSection(),
                       total_cross_section );
  }

  // Add the coherent cs
  this->addCrossSectionToTotalCrossSection(
                       energy_grid,
                       data_container.getWallerHartreeCoherentCrossSection(),
                       total_cross_section );

  // Add the pair production cs
  this->addCrossSectionToTotalCrossSection(
                       energy_grid,
                       data_container.getPairProductionCrossSection(),
                       total_cross_section );

  // Add the triplet production cs
  this->addCrossSectionToTotalCrossSection(
                       energy_grid,
                       data_container.getTripletProductionCrossSection(),
                       total_cross_section );

  // Add the photoelectric cs
  this->addCrossSectionToTotalCrossSection(
                       energy_grid,
                       data_container.getPhotoelectricCrossSection(),
                       total_cross_section );

  if( use_waller_hartree_incoherent_cs )
    data_container.setWallerHartreeTotalCrossSection( total_cross_section );
  else
    data_container.setImpulseApproxTotalCrossSection( total_cross_section );
}

// Calculate the electron total cross section
void ENDLElectronPhotonRelaxationDataGenerator::calculateElectronTotalCrossSection()
{
  Data::ElectronPhotonRelaxationVolatileDataContainer& data_container =
    this->getVolatileDataContainer();

  const std::vector<double>& energy_grid =
    data_container.getElectronEnergyGrid();

  std::vector<double> total_cross_section( energy_grid.size(), 0.0 );

  // Add the total elastic cs
  this->addCrossSectionToTotalCrossSection(
                       energy_grid,
                       data_container.getTotalElasticCrossSection(),
                       total_cross_section );

  // Add the atomic excitation cs
  this->addCrossSectionToTotalCrossSection(
                       energy_grid,
                       data_container.getAtomicExcitationCrossSection(),
                       total_cross_section );

  // Add the bremsstrahlung cs
  this->addCrossSectionToTotalCrossSection(
                       energy_grid,
                       data_container.getBremsstrahlungCrossSection(),
                       total_cross_section );

  std::set<unsigned>::iterator shell = data_container.getSubshells().begin();

  // Loop through electroionization data for every subshell
  for ( shell; shell != data_container.getSubshells().end(); ++shell )
  {
    // Add the electroionization subshell cs
    this->addCrossSectionToTotalCrossSection(
                        energy_grid,
                        data_container.getElectroionizationCrossSection(*shell),
                        total_cross_section );
  }

  data_container.setTotalElectronCrossSection( total_cross_section );
}

// Add cross section to total cross section
void ENDLElectronPhotonRelaxationDataGenerator::addCrossSectionToTotalCrossSection(
                               const std::vector<double>& energy_grid,
                               const std::vector<double>& cross_section,
                               std::vector<double>& total_cross_section ) const
{
  unsigned start_index = energy_grid.size() - cross_section.size();

  for( unsigned i = 0; i < cross_section.size(); ++i )
    total_cross_section[start_index+i] += cross_section[i];
}

// Calculate the elastic angular distribution for the angle cosine
void ENDLElectronPhotonRelaxationDataGenerator::calculateElasticAngleCosine(
    const std::vector<double>& raw_elastic_angle,
    const std::vector<double>& raw_elastic_pdf,
    std::vector<double>& elastic_angle,
    std::vector<double>& elastic_pdf ) const
{
  int size = raw_elastic_angle.size();
  int r_bin = size - 1;

  elastic_angle.resize( size );
  elastic_pdf.resize( size );

  for ( unsigned bin = 0; bin < size; ++bin )
  {
    elastic_pdf[r_bin] = raw_elastic_pdf[bin];
    long double angle_cosine = 1.0L - raw_elastic_angle[bin];
    elastic_angle[r_bin] = angle_cosine;
    r_bin--;
  }
}

// Calculate the elastic moment preserving discrete angle cosines and weights
void ENDLElectronPhotonRelaxationDataGenerator::calculateDiscreteAnglesAndWeights(
    const std::shared_ptr<DataGen::ElasticElectronMomentsEvaluator>& moments_evaluator,
    const double& energy,
    const int& number_of_moment_preserving_angles,
    std::vector<double>& discrete_angles,
    std::vector<double>& weights,
    double& cross_section_reduction )
{
  std::vector<Utility::long_float> legendre_moments;
  double precision = 1e-13;
  int n = ( number_of_moment_preserving_angles+1 )*2 -2;

  // Get the elastic moments
  moments_evaluator->evaluateElasticMoment( legendre_moments,
                                            energy,
                                            n,
                                            precision );

  // Use radau quadrature to find the discrete angles and weights from the moments
  std::shared_ptr<Utility::SloanRadauQuadrature> radau_quadrature(
      new Utility::SloanRadauQuadrature( legendre_moments ) );

  radau_quadrature->getRadauNodesAndWeights( discrete_angles,
                                             weights,
                                             number_of_moment_preserving_angles+1 );

  // Eliminate the forward discrete angle (mu = 1)
  discrete_angles.pop_back();
  weights.pop_back();

  // Re-normalize weights and set the cross_section_reduction to the sum of the weights
  cross_section_reduction = 0.0;
  for( int i = 0; i < weights.size(); ++i )
  {
    cross_section_reduction += weights[i];
  }

  for( int i = 0; i < weights.size(); ++i )
  {
    weights[i] /= cross_section_reduction;
  }
}

// Calculate the electron total elastic cross section
/*! \details The ENDL tables have a total elastic cross section alternatively
 * the total elastic cross section can also be calculated from integrating over
 * the cutoff distribution and the normalized screened Rutherford distribution.
 * The two versions of the cross section will usually be of the same magnitude,
 * but will not agree due to the course 2D grid and roundoff error.
 */
void ENDLElectronPhotonRelaxationDataGenerator::calculateElectronTotalElasticCrossSection(
                        std::shared_ptr<const Utility::UnivariateDistribution>&
                        total_elastic_cross_section,
                        const std::vector<double>& raw_energy_grid )
{
  Data::ElectronPhotonRelaxationVolatileDataContainer& data_container =
    this->getVolatileDataContainer();

  // Check to see if the ENDL or integrated cross section is wanted
  if ( !this->isElectronTotalElasticIntegratedCrossSectionModeOn() )
  {
    // Get total elastic cross section provided with ENDL
    total_elastic_cross_section.reset(
      new Utility::TabularDistribution<Utility::LogLog>(
        raw_energy_grid,
        d_endl_data_container->getTotalElasticCrossSection() ) );
  }
  else
  {
    std::shared_ptr<const MonteCarlo::CutoffElasticElectronScatteringDistribution>
        cutoff_endl_distribution;

    if ( this->getElectronTwoDInterpPolicy() == MonteCarlo::LOGLOGLOG_INTERPOLATION )
    {
      MonteCarlo::ElasticElectronScatteringDistributionNativeFactory::createCutoffElasticDistribution<Utility::LogLogCosLog,Utility::Correlated>(
            cutoff_endl_distribution,
            data_container.getCutoffElasticAngles(),
            data_container.getCutoffElasticPDF(),
            data_container.getElasticAngularEnergyGrid(),
            MonteCarlo::ElasticElectronTraits::mu_peak,
            this->getTabularEvaluationTolerance() );
    }
    else if ( this->getElectronTwoDInterpPolicy() == MonteCarlo::LINLINLIN_INTERPOLATION )
    {
      MonteCarlo::ElasticElectronScatteringDistributionNativeFactory::createCutoffElasticDistribution<Utility::LinLinLin,Utility::Correlated>(
            cutoff_endl_distribution,
            data_container.getCutoffElasticAngles(),
            data_container.getCutoffElasticPDF(),
            data_container.getElasticAngularEnergyGrid(),
            MonteCarlo::ElasticElectronTraits::mu_peak,
            this->getTabularEvaluationTolerance() );
    }
    else if ( this->getElectronTwoDInterpPolicy() == MonteCarlo::LINLINLOG_INTERPOLATION )
    {
      MonteCarlo::ElasticElectronScatteringDistributionNativeFactory::createCutoffElasticDistribution<Utility::LinLinLog,Utility::Correlated>(
            cutoff_endl_distribution,
            data_container.getCutoffElasticAngles(),
            data_container.getCutoffElasticPDF(),
            data_container.getElasticAngularEnergyGrid(),
            MonteCarlo::ElasticElectronTraits::mu_peak,
            this->getTabularEvaluationTolerance() );
    }

    // Create the elastic traits
    std::shared_ptr<MonteCarlo::ElasticElectronTraits> elastic_traits(
      new MonteCarlo::ElasticElectronTraits( data_container.getAtomicNumber() ) );

    std::vector<double> raw_elastic_cross_section =
        d_endl_data_container->getCutoffElasticCrossSection();

    // Calculate the total elastic cross section
    for (unsigned n = 0; n < raw_elastic_cross_section.size(); ++n)
    {
      // Get the energy
      double energy = raw_energy_grid[n];

      // Get the cutoff pdf
      double cutoff_pdf = cutoff_endl_distribution->evaluate(
                            energy, MonteCarlo::ElasticElectronTraits::mu_peak );

      // Get Moliere's Screening Constant
      double eta = elastic_traits->evaluateMoliereScreeningConstant( energy );

      // Calculate the integrated PDF value at the cutoff angle cosine
      double integrated_pdf =
        cutoff_pdf*MonteCarlo::ElasticElectronTraits::delta_mu_peak*
            ( MonteCarlo::ElasticElectronTraits::delta_mu_peak + eta )/( eta );

      // Evaluate the total coupled cross section at the incoming energy
      raw_elastic_cross_section[n] *= ( 1.0 +  integrated_pdf );
    }

    // Get total elastic cross section (same energy grid as cutoff)
    total_elastic_cross_section.reset(
      new Utility::TabularDistribution<Utility::LogLog>(
        raw_energy_grid,
        raw_elastic_cross_section ) );
  }
}

// Return the endl data container
const Data::ENDLDataContainer& ENDLElectronPhotonRelaxationDataGenerator::getENDLDataContainer() const
{
  return *d_endl_data_container;
}

} // end DataGen namespace

//---------------------------------------------------------------------------//
// end DataGen_ENDLElectronPhotonRelaxationDataGenerator.hpp
//---------------------------------------------------------------------------//
