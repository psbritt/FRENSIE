//---------------------------------------------------------------------------//
//!
//! \file   DataGen_StandardENDLDataGenerator.cpp
//! \author Luke Kersting
//! \brief  The standard endl data generator class def.
//!
//---------------------------------------------------------------------------//

// FRENSIE Includes
#include "DataGen_StandardENDLDataGenerator.hpp"
#include "Data_ENDLFileHandler.hpp"
#include "Data_ENDLHelperWrappers.hpp"
#include "Data_SubshellType.hpp"
#include "Utility_LoggingMacros.hpp"
#include "Utility_ExceptionTestMacros.hpp"
#include "Utility_DesignByContract.hpp"

namespace DataGen{

// Constructor
StandardENDLDataGenerator::StandardENDLDataGenerator(
                                const boost::filesystem::path& eadl_file_name,
                                const boost::filesystem::path& epdl_file_name,
                                const boost::filesystem::path& eedl_file_name )
  : ENDLDataGenerator(),
    d_eadl_file_name( eadl_file_name ),
    d_epdl_file_name( epdl_file_name ),
    d_eedl_file_name( eedl_file_name )
{
  d_eadl_file_name.make_preferred();
  
  TEST_FOR_EXCEPTION( !boost::filesystem::exists( eadl_file_name ),
                      std::runtime_error,
                      "The requested eadl file does not exist!" );

  d_epdl_file_name.make_preferred();
  
  TEST_FOR_EXCEPTION( !boost::filesystem::exists( epdl_file_name ),
                      std::runtime_error,
                      "The requested epdl file does not exist!" );

  d_eedl_file_name.make_preferred();
  
  TEST_FOR_EXCEPTION( !boost::filesystem::exists( eedl_file_name ),
                      std::runtime_error,
                      "The requested eedl file does not exist!" );
}

// Populate the ENDL data container
void StandardENDLDataGenerator::populateENDLDataContainer()
{
  // Populate the relaxation data
  this->populateEADLDataContainer();

  // Populate the photon data
  this->populateEPDLDataContainer();

  // Populate the electron data
  this->populateEEDLDataContainer();
}

// Populate the EADL data container
void StandardENDLDataGenerator::populateEADLDataContainer()
{
  // Set the relaxation data
  FRENSIE_LOG_PARTIAL_NOTIFICATION( "Setting the relaxation data " );
  FRENSIE_FLUSH_ALL_LOGS();

  // Open eadl file
  std::shared_ptr<Data::ENDLFileHandler>
    eadl_file_handler( new Data::ENDLFileHandler( d_eadl_file_name.string() ) );

  // Information in first header of the EADL file
  int atomic_number_in_table,
    outgoing_particle_designator,
    interpolation_flag;
  double atomic_weight;
  
  bool atomic_number_set = false;
  bool atomic_weight_set = false;
  
  // Information in the second header of the EADL file
  int reaction_type, electron_shell;
  
  // subshells and subshell data
  bool convert_subshell = true;
  std::vector<unsigned> subshells;
  std::set<unsigned> endf_subshells;
  unsigned endf_subshell;
  
  Data::ENDLVolatileDataContainer& data_container =
    this->getVolatileDataContainer();
  
  // Process every table in the EADL file
  while( eadl_file_handler->validFile() && !eadl_file_handler->endOfFile() )
  {
    // Read first table header and determine which element is being processed
    eadl_file_handler->readFirstTableHeader( atomic_number_in_table,
                                             outgoing_particle_designator,
                                             atomic_weight,
                                             interpolation_flag );
    
    if( !atomic_number_set )
    {
      data_container.setAtomicNumber( atomic_number_in_table );
      atomic_number_set = true;
    }

    if( !atomic_weight_set )
    {
      data_container.setAtomicWeight( atomic_weight );
      atomic_weight_set = true;
    }
    
    // Check that the EADL file is still valid (eof has not been reached)
    if( eadl_file_handler->endOfFile() )
      continue;
    
    TEST_FOR_EXCEPTION( atomic_number_in_table !=
                        data_container.getAtomicNumber(),
                        std::runtime_error,
                        "The atomic number in the table ("
                        << atomic_number_in_table << ") does not match the "
                        "atomic number set in the data container ("
                        << data_container.getAtomicNumber() <<
                        ")!" );
    
    // Read second table header and determine the reaction type
    eadl_file_handler->readSecondTableHeader( reaction_type,
                                              electron_shell );
    
    if( electron_shell > 0 )
    {
      // Convert subshell number to endf number
      endf_subshell =
        Data::convertEADLDesignatorToENDFDesignator( electron_shell );
      
      // insert subshell to set
      endf_subshells.insert( endf_subshell );
    }
    
    // Read and process the data in the current table, then store in the HDF5
    // file
    switch( reaction_type )
    {
      
      // Number of electrons in subshell
      case 91912:
      {
        // Interpolation should always be LinLin = 0
        testInvariant( interpolation_flag == 0 )
          
        std::map<unsigned,double> subshell_data;
        eadl_file_handler->mapTwoColumnTable( subshells,
                                              subshell_data,
                                              convert_subshell );

        // convert subshells to a set
        for( int i = 0; i < subshells.size(); ++i )
          endf_subshells.insert( subshells[i] );

        // set the subshells
        data_container.setSubshells( endf_subshells );

        // set the subshell data
        data_container.setSubshellOccupancy( subshell_data );

        FRENSIE_LOG_PARTIAL_NOTIFICATION( "." );
        FRENSIE_FLUSH_ALL_LOGS();
        
        break;
      }

      // Binding energy of a subshell
      case 91913:
      {
        // Interpolation should always be LinLin = 0
        testInvariant( interpolation_flag == 0 )

        std::map<unsigned,double> subshell_data;
        eadl_file_handler->mapTwoColumnTable( subshells,
                                              subshell_data,
                                              convert_subshell );

        testInvariant( subshells.size() == endf_subshells.size() );

        // set the subshell data
        data_container.setSubshellBindingEnergy( subshell_data );

        FRENSIE_LOG_PARTIAL_NOTIFICATION( "." );
        FRENSIE_FLUSH_ALL_LOGS();

        break;
      }

      // Kinetic energy of a subshell
      case 91914:
      {
        // Interpolation should always be LinLin = 0
        testInvariant( interpolation_flag == 0 )

        std::map<unsigned,double> subshell_data;
        eadl_file_handler->mapTwoColumnTable( subshells,
                                              subshell_data,
                                              convert_subshell );

        testInvariant( subshells.size() == endf_subshells.size() );

        // set the subshell data
        data_container.setSubshellKineticEnergy( subshell_data );

        FRENSIE_LOG_PARTIAL_NOTIFICATION( "." );
        FRENSIE_FLUSH_ALL_LOGS();

        break;
      }

      // Average radius of a subshell
      case 91915:
      {
        // Interpolation should always be LinLin = 0
        testInvariant( interpolation_flag == 0 )

        std::map<unsigned,double> subshell_data;
        eadl_file_handler->mapTwoColumnTable( subshells,
                                              subshell_data,
                                              convert_subshell );

        testInvariant( subshells.size() == endf_subshells.size() );

        // set the subshell data
        data_container.setSubshellAverageRadius( subshell_data );

        FRENSIE_LOG_PARTIAL_NOTIFICATION( "." );
        FRENSIE_FLUSH_ALL_LOGS();

        break;
      }

      // Radiative level width of a subshell
      case 91921:
      {
        // Interpolation should always be LinLin = 0
        testInvariant( interpolation_flag == 0 )

        std::map<unsigned,double> subshell_data;
        eadl_file_handler->mapTwoColumnTable( subshells,
                                              subshell_data,
                                              convert_subshell );

        // set the subshell data
        data_container.setSubshellRadiativeLevel( subshell_data );

        FRENSIE_LOG_PARTIAL_NOTIFICATION( "." );
        FRENSIE_FLUSH_ALL_LOGS();
        
        break;
      }

      // Non radiative level of a subshell
      case 91922:
      {
        // Interpolation should always be LinLin = 0
        testInvariant( interpolation_flag == 0 )

        std::map<unsigned,double> subshell_data;
        eadl_file_handler->mapTwoColumnTable( subshells,
                                              subshell_data,
                                              convert_subshell );

        // set the subshell data
        data_container.setSubshellNonRadiativeLevel( subshell_data );

        FRENSIE_LOG_PARTIAL_NOTIFICATION( "." );
        FRENSIE_FLUSH_ALL_LOGS();

        break;
      }

      // Radiative transition probability per subshell
      case 92931:
      {
        std::map<unsigned,double> indep_subshell_data, dep_subshell_data;
        eadl_file_handler->mapThreeColumnTable( subshells,
                                                indep_subshell_data,
                                                dep_subshell_data,
                                                true );

        data_container.setRadiativeTransitionProbability( endf_subshell,
                                                          indep_subshell_data );
        data_container.setRadiativeTransitionEnergy( endf_subshell,
                                                     dep_subshell_data );

        FRENSIE_LOG_PARTIAL_NOTIFICATION( "." );
        FRENSIE_FLUSH_ALL_LOGS();

        break;
      }

      // Nonradiative transition probability per subshell
      case 92932:
      {
        std::map<unsigned,std::vector<unsigned> > secondary_subshells;
        std::map<unsigned,std::map<unsigned,double> >
            indep_subshell_data, dep_subshell_data;
        eadl_file_handler->mapFourColumnTable( subshells,
                                               secondary_subshells,
                                               indep_subshell_data,
                                               dep_subshell_data,
                                               true );

        data_container.setNonRadiativeTransitionProbability( endf_subshell,
                                                             indep_subshell_data );
        data_container.setNonRadiativeTransitionEnergy( endf_subshell,
                                                        dep_subshell_data );

        FRENSIE_LOG_PARTIAL_NOTIFICATION( "." );
        FRENSIE_FLUSH_ALL_LOGS();

        break;
      }

      // The average number of particles per initial vacancy of a subshell
      case 92933:
      {
        // Interpolation should always be LinLin = 0
        testInvariant( interpolation_flag == 0 )

        std::map<unsigned,double> subshell_data;
        eadl_file_handler->mapTwoColumnTable( subshells,
                                              subshell_data,
                                              convert_subshell );

        // Average number of photons emitted per initial vacancy ( Yo == 7 )
        if( outgoing_particle_designator == 7 )
        {
          // set the subshell data
          data_container.setAveragePhotonsPerInitialVacancy( subshell_data );
        }
        // Average number of electrons emitted per initial vacancy ( Yo == 9 )
        else if( outgoing_particle_designator = 9 )
        {
          // set the subshell data
          data_container.setAverageElectronsPerInitialVacancy( subshell_data );
        }
        else
        {
          THROW_EXCEPTION( std::runtime_error,
                           "Unsupported outgoing particle designator "
                           "encountered in EADL table (Yo = "
                           << outgoing_particle_designator << ")!" );
        }

        FRENSIE_LOG_PARTIAL_NOTIFICATION( "." );
        FRENSIE_FLUSH_ALL_LOGS();

        break;
      }

      // The average energy of particles per initial vacancy
      case 92934:
      {
        // Interpolation should always be LinLin = 0
        testInvariant( interpolation_flag == 0 )

        std::map<unsigned,double> subshell_data;
        eadl_file_handler->mapTwoColumnTable( subshells,
                                              subshell_data,
                                              convert_subshell );

        // Average energy of photons emitted per initial vacancy ( Yo == 7 )
        if( outgoing_particle_designator == 7 )
        {
          // set the subshell data
          data_container.setAveragePhotonEnergyPerInitialVacancy( subshell_data );
        }
        
        // Average energy of electrons emitted per initial vacancy ( Yo == 9 )
        else if( outgoing_particle_designator == 9 )
        {
          // set the subshell data
          data_container.setAverageElectronEnergyPerInitialVacancy( subshell_data );
        }
        else
        {
          THROW_EXCEPTION( std::runtime_error,
                           "Unsupported outgoing particle designator "
                           "encountered in EADL table (Yo = "
                           << outgoing_particle_designator << ")!" );
        }

        FRENSIE_LOG_PARTIAL_NOTIFICATION( "." );
        FRENSIE_FLUSH_ALL_LOGS();

        break;
      }

      // The local deposition per initial vacancy of a subshell
      case 92935:
      {
        // Interpolation should always be LinLin = 0
        testInvariant( interpolation_flag == 0 )

        std::map<unsigned,double> subshell_data;
        eadl_file_handler->mapTwoColumnTable( subshells,
                                              subshell_data,
                                              convert_subshell );

        testInvariant( subshells.size() == endf_subshells.size() );

        // set the subshell data
        data_container.setLocalDepositionPerInitialVacancy( subshell_data );

        FRENSIE_LOG_PARTIAL_NOTIFICATION( "." );
        FRENSIE_FLUSH_ALL_LOGS();

        break;
      }

      // Unknown reaction type found
      default:
      { 
        THROW_EXCEPTION( std::runtime_error,
                         "An unknown reaction type was encountered while "
                         "processing the EADL file!" );
      }
    }
  }

  // Close the EADL file
  eadl_file_handler->closeENDLFile();
  
  FRENSIE_LOG_NOTIFICATION( " done." );
}

// Populate the EPDL data container
void StandardENDLDataGenerator::populateEPDLDataContainer()
{
  // Set the photon data
  FRENSIE_LOG_PARTIAL_NOTIFICATION( "Setting the photon data " );
  FRENSIE_FLUSH_ALL_LOGS();
  
  // Open epdl file
  std::shared_ptr<Data::ENDLFileHandler>
    epdl_file_handler( new Data::ENDLFileHandler( d_epdl_file_name.string() ) );

  // Information in first header of the EPDL file
  int atomic_number_in_table,
    outgoing_particle_designator,
    interpolation_flag;
  double atomic_weight;

  // Information in the second header of the EPDL file
  int reaction_type, electron_shell;

  // array of all the subshells read
  unsigned endf_subshell;
  std::set<unsigned> endf_subshells;

  Data::ENDLVolatileDataContainer& data_container =
    this->getVolatileDataContainer();

  // Process every table in the EPDL file
  while( epdl_file_handler->validFile() && !epdl_file_handler->endOfFile() )
  {
    // Read first table header and determine which element is being processed
    epdl_file_handler->readFirstTableHeader( atomic_number_in_table,
                                             outgoing_particle_designator,
                                             atomic_weight,
                                             interpolation_flag );

    // Check that the EPDL file is still valid (eof has not been reached)
    if( epdl_file_handler->endOfFile() )
      continue;

    TEST_FOR_EXCEPTION( atomic_number_in_table !=
                        data_container.getAtomicNumber(),
                        std::runtime_error,
                        "The atomic number in the table ("
                        << atomic_number_in_table << ") does not match the "
                        "atomic number set in the data container ("
                        << data_container.getAtomicNumber() <<
                        ")!" );

    // Read second table header and determine the reaction type
    epdl_file_handler->readSecondTableHeader( reaction_type, electron_shell );

    if( electron_shell > 0 )
    {
      // Convert subshell number to endf number
      endf_subshell =
        Data::convertEADLDesignatorToENDFDesignator( electron_shell );

      // insert subshell to set
      endf_subshells.insert( endf_subshell );
    }

    // Read and process the data in the current table, then store in the HDF5
    // file
    switch( reaction_type )
    {
      // Read in the integrated coherent cross section data
      case 71000:
      {
        std::vector<double> indep_data, dep_data;
        epdl_file_handler->processTwoColumnTable( indep_data, dep_data );

        data_container.setCoherentCrossSectionEnergyGrid( indep_data );
        data_container.setCoherentCrossSection( dep_data );

        FRENSIE_LOG_PARTIAL_NOTIFICATION( "." );
        FRENSIE_FLUSH_ALL_LOGS();
        
        break;
      }
      
      // Average energy of scattered photon from coherent scattering ignored
      case 71010:
      {
        std::vector<double> indep_data, dep_data;
        epdl_file_handler->processTwoColumnTable( indep_data, dep_data );

        data_container.setCoherentAveragePhotonIncidentEnergy( indep_data );
        data_container.setCoherentAveragePhotonEnergy( dep_data );

        FRENSIE_LOG_PARTIAL_NOTIFICATION( "." );
        FRENSIE_FLUSH_ALL_LOGS();
        
        break;
      }

      // Read in the integrated incoherent cross section data
      case 72000:
      {
        // Interpolation should always be LogLog = 5
        testInvariant( interpolation_flag == 5 )

        std::vector<double> indep_data, dep_data;
        epdl_file_handler->processTwoColumnTable( indep_data, dep_data );

        data_container.setIncoherentCrossSectionEnergyGrid( indep_data );
        data_container.setIncoherentCrossSection( dep_data );

        FRENSIE_LOG_PARTIAL_NOTIFICATION( "." );
        FRENSIE_FLUSH_ALL_LOGS();
        
        break;
      }

      // Average energy of scattered particle from incoherent scattering
      case 72010:
      {
        std::vector<double> indep_data, dep_data;
        epdl_file_handler->processTwoColumnTable( indep_data, dep_data );

        // Average energy of scattered photon from incoherent scattering ( Yo == 7 )
        if( outgoing_particle_designator == 7 )
        {
          data_container.setIncoherentAveragePhotonIncidentEnergy( indep_data );
          data_container.setIncoherentAveragePhotonEnergy( dep_data );
        }
        // Average energy of scattered electron from incoherent scattering ( Yo == 9 )
        else if( outgoing_particle_designator == 9 )
        {
          data_container.setIncoherentAverageElectronIncidentEnergy( indep_data );
          data_container.setIncoherentAverageElectronEnergy( dep_data );
        }
        else
        {
          THROW_EXCEPTION( std::runtime_error,
                           "Unsupported outgoing particle designator "
                           "encountered in EPDL table (Yo = "
                           << outgoing_particle_designator << ")!" );
        }

        FRENSIE_LOG_PARTIAL_NOTIFICATION( "." );
        FRENSIE_FLUSH_ALL_LOGS();
        
        break;
      }

      // Read in the integrated photoelectric cross section data
      case 73000:
      {
        // Interpolation should always be LogLog = 5
        testInvariant( interpolation_flag == 5 )

        std::vector<double> indep_data, dep_data;
        epdl_file_handler->processTwoColumnTable( indep_data, dep_data );

        // Read the total integrated photoelectric cross section
        if( electron_shell == 0 )
        {
          data_container.setPhotoelectricCrossSectionEnergyGrid( indep_data );
          data_container.setPhotoelectricCrossSection( dep_data );
        }
        else
        {
          data_container.setPhotoelectricCrossSectionEnergyGrid( endf_subshell,
                                                                 indep_data );
          data_container.setPhotoelectricCrossSection( endf_subshell,
                                                       dep_data );
        }

        FRENSIE_LOG_PARTIAL_NOTIFICATION( "." );
        FRENSIE_FLUSH_ALL_LOGS();
        
        break;
      }

      // Average energy of secondary particle from photoelectric effect
      case 73010:
      {
        std::vector<double> indep_data, dep_data;
        epdl_file_handler->processTwoColumnTable( indep_data, dep_data );

        if( electron_shell == 0 )
        {
          // Average energy of scattered photon from photoelectric effect ( Yo == 7 )
          if( outgoing_particle_designator == 7 )
          {
            data_container.setPhotoelectricAveragePhotonsIncidentEnergy( indep_data );
            data_container.setPhotoelectricAveragePhotonsEnergy( dep_data );
          }
          // Average energy of scattered electron from photoelectric effect ( Yo == 9 )
          else if( outgoing_particle_designator == 9 )
          {
            data_container.setPhotoelectricAverageElectronsIncidentEnergy( indep_data );
            data_container.setPhotoelectricAverageElectronsEnergy( dep_data );
          }
          else
          {
            THROW_EXCEPTION( std::runtime_error,
                             "Unsupported outgoing particle designator "
                             "encountered in EPDL table (Yo = "
                             << outgoing_particle_designator << ")!" );
          }
        }
        else
        {
          // Average energy of scattered photon from photoelectric effect ( Yo == 7 )
          if ( outgoing_particle_designator == 7 )
          {
            data_container.setPhotoelectricAveragePhotonsIncidentEnergy(
                endf_subshell,
                indep_data );
            data_container.setPhotoelectricAveragePhotonsEnergy(
                endf_subshell,
                dep_data );
          }
          // Average energy of scattered electron from photoelectric effect ( Yo == 9 )
          else if( outgoing_particle_designator == 9 )
          {
            data_container.setPhotoelectricAverageElectronsIncidentEnergy(
                endf_subshell,
                indep_data );
            data_container.setPhotoelectricAverageElectronsEnergy(
                endf_subshell,
                dep_data );
          }
          else
          {
            THROW_EXCEPTION( std::runtime_error,
                             "Unsupported outgoing particle designator "
                             "encountered in EPDL table (Yo = "
                             << outgoing_particle_designator << ")!" );
          }
        }

        FRENSIE_LOG_PARTIAL_NOTIFICATION( "." );
        FRENSIE_FLUSH_ALL_LOGS();
        
        break;
      }

      // Average energy to residual atom from photoelectric effect
      case 73011:
      {
        std::vector<double> indep_data, dep_data;
        epdl_file_handler->processTwoColumnTable( indep_data, dep_data );

        if( electron_shell == 0 )
        {
          data_container.setPhotoelectricAverageResidualIncidentEnergy(
            indep_data );
          data_container.setPhotoelectricAverageResidualEnergy( dep_data );
        }
        else
        {
          data_container.setPhotoelectricAverageResidualIncidentEnergy(
            endf_subshell,
            indep_data );
          data_container.setPhotoelectricAverageResidualEnergy(
            endf_subshell,
            dep_data );
        }

        FRENSIE_LOG_PARTIAL_NOTIFICATION( "." );
        FRENSIE_FLUSH_ALL_LOGS();
        
        break;
      }

      // Read the integrated pair production cross section
      case 74000:
      {
        // Interpolation should always be LogLog = 5
        testInvariant( interpolation_flag == 5 )

        std::vector<double> indep_data, dep_data;
        epdl_file_handler->processTwoColumnTable(
            indep_data,
            dep_data );

        data_container.setPairProductionCrossSectionEnergyGrid( indep_data );
        data_container.setPairProductionCrossSection( dep_data );

        FRENSIE_LOG_PARTIAL_NOTIFICATION( "." );
        FRENSIE_FLUSH_ALL_LOGS();
        
        break;
      }

      // Average energy of secondary particle from pair production
      case 74010:
      {
        std::vector<double> indep_data, dep_data;
        epdl_file_handler->processTwoColumnTable( indep_data, dep_data );

        // Average energy of secondary positron from pair production ( Yo == 8 )
        if( outgoing_particle_designator == 8 )
        {
          data_container.setPairProductionAveragePositronIncidentEnergy(
            indep_data );
          data_container.setPairProductionAveragePositronEnergy( dep_data );
        }
        // Average energy of secondary electron from pair production ( Yo == 9 )
        else if( outgoing_particle_designator == 9 )
        {
          data_container.setPairProductionAverageElectronIncidentEnergy(
            indep_data );
          data_container.setPairProductionAverageElectronEnergy( dep_data );
        }
        else
        {
          THROW_EXCEPTION( std::runtime_error,
                           "Unsupported outgoing particle designator "
                           "encountered in EPDL table (Yo = "
                           << outgoing_particle_designator << ")!" );
        }

        FRENSIE_LOG_PARTIAL_NOTIFICATION( "." );
        FRENSIE_FLUSH_ALL_LOGS();
        
        break;
      }

      // Read the integrated triplet production cross section
      case 75000:
      {
        // Interpolation should always be LogLog = 5
        testInvariant( interpolation_flag == 5 )

        std::vector<double> indep_data, dep_data;
        epdl_file_handler->processTwoColumnTable( indep_data,
                                                  dep_data );

        data_container.setTripletProductionCrossSectionEnergyGrid( indep_data );
        data_container.setTripletProductionCrossSection( dep_data );

        FRENSIE_LOG_PARTIAL_NOTIFICATION( "." );
        FRENSIE_FLUSH_ALL_LOGS();
        
        break;
      }

      // Average energy of secondary particle from triplet production
      case 75010:
      {
        std::vector<double> indep_data, dep_data;
        epdl_file_handler->processTwoColumnTable( indep_data, dep_data );

        // Average energy of secondary positron from triplet production ( Yo == 8 )
        if( outgoing_particle_designator == 8 )
        {
          data_container.setTripletProductionAveragePositronIncidentEnergy(
            indep_data );
          data_container.setTripletProductionAveragePositronEnergy( dep_data );
        }
        // Average energy of secondary electron from triplet production ( Yo == 9 )
        else if( outgoing_particle_designator == 9 )
        {
          data_container.setTripletProductionAverageElectronIncidentEnergy(
            indep_data );
          data_container.setTripletProductionAverageElectronEnergy( dep_data );
        }
        else
        {
          THROW_EXCEPTION( std::runtime_error,
                           "Unsupported outgoing particle designator "
                           "encountered in EPDL table (Yo = "
                           << outgoing_particle_designator << ")!" );
        }

        FRENSIE_LOG_PARTIAL_NOTIFICATION( "." );
        FRENSIE_FLUSH_ALL_LOGS();
        
        break;
      }

      // Read the atomic form factor
      case 93941:
      {
        // Interpolation should always be LogLog = 5
        testInvariant( interpolation_flag == 5 )

        std::vector<double> indep_data, dep_data;
        epdl_file_handler->processTwoColumnTable( indep_data, dep_data );

        data_container.setCoherentFormFactorArgument( indep_data );
        data_container.setCoherentFormFactor( dep_data );

        FRENSIE_LOG_PARTIAL_NOTIFICATION( "." );
        FRENSIE_FLUSH_ALL_LOGS();
        
        break;
      }

      // Read the scattering function
      case 93942:
      {
        // Interpolation should always be LogLog = 5
        testInvariant( interpolation_flag == 5 )

        std::vector<double> indep_data, dep_data;
        epdl_file_handler->processTwoColumnTable( indep_data, dep_data );

        data_container.setIncoherentScatteringFunctionArgument( indep_data );
        data_container.setIncoherentScatteringFunction( dep_data );

        FRENSIE_LOG_PARTIAL_NOTIFICATION( "." );
        FRENSIE_FLUSH_ALL_LOGS();
        
        break;
      }
      
      // Imaginary anomalous scattering factor
      case 93943:
      {
        // Interpolation should always be LinLin = 2
        testInvariant( interpolation_flag == 2 )

        std::vector<double> indep_data, dep_data;
        epdl_file_handler->processTwoColumnTable( indep_data, dep_data );

        data_container.setCoherentImaginaryAnomalousFactorIncidentEnergy(
            indep_data );
        data_container.setCoherentImaginaryAnomalousFactor( dep_data );

        FRENSIE_LOG_PARTIAL_NOTIFICATION( "." );
        FRENSIE_FLUSH_ALL_LOGS();
        
        break;
      }

      // Real anomalous scattering factor
      case 93944:
      {
        // Interpolation should always be LinLin = 2
        testInvariant( interpolation_flag == 2 )

        std::vector<double> indep_data, dep_data;
        epdl_file_handler->processTwoColumnTable( indep_data, dep_data );

        data_container.setCoherentRealAnomalousFactorIncidentEnergy(
            indep_data );
        data_container.setCoherentRealAnomalousFactor( dep_data );

        FRENSIE_LOG_PARTIAL_NOTIFICATION( "." );
        FRENSIE_FLUSH_ALL_LOGS();
        
        break;
      }
      
      // Unknown reaction type found
      default:
      {
        THROW_EXCEPTION( std::runtime_error,
                         "An unknown reaction type was encountered while "
                         "processing the EPDL file!" );
      }
    }
  }

  // Close the EPDL file
  epdl_file_handler->closeENDLFile();

  // Set the subshells
  if ( data_container.getSubshells().empty() )
    data_container.setSubshells( endf_subshells );

  FRENSIE_LOG_NOTIFICATION( "done." );
}

// Populate the EEDL data container
void StandardENDLDataGenerator::populateEEDLDataContainer()
{
  // Set the electron data
  FRENSIE_LOG_PARTIAL_NOTIFICATION( "Setting the electron data " );
  FRENSIE_FLUSH_ALL_LOGS();
  
  // Open eedl file
  std::shared_ptr<Data::ENDLFileHandler>
    eedl_file_handler( new Data::ENDLFileHandler( d_eedl_file_name.string() ) );

  // Information in first header of the EEDL file
  int atomic_number_in_table,
    outgoing_particle_designator,
    interpolation_flag;
  double atomic_weight;

  // Information in the second header of the EEDL file
  int reaction_type, electron_shell;
  unsigned endf_subshell;

  // array of all the subshells read
  std::set<unsigned> endf_subshells;

  // The elastic angular distribution energy grid and pdf
  std::vector<double> elastic_angular_energy_grid;
  std::map<double,std::vector<double> > elastic_pdf;

  Data::ENDLVolatileDataContainer& data_container =
    this->getVolatileDataContainer();

  // Process every table in the EEDL file
  while( eedl_file_handler->validFile() && !eedl_file_handler->endOfFile() )
  {
    // Read first table header and determine which element is being processed
    eedl_file_handler->readFirstTableHeader( atomic_number_in_table,
                                             outgoing_particle_designator,
                                             atomic_weight,
                                             interpolation_flag );

    // Check that the EEDL file is still valid (eof has not been reached)
    if( eedl_file_handler->endOfFile() )
      continue;

    TEST_FOR_EXCEPTION( atomic_number_in_table !=
                        data_container.getAtomicNumber(),
                        std::runtime_error,
                        "The atomic number in the table ("
                        << atomic_number_in_table << ") does not match the "
                        "atomic number set in the data container ("
                        << data_container.getAtomicNumber() <<
                        ")!" );

    // Read second table header and determine the reaction type
    eedl_file_handler->readSecondTableHeader( reaction_type,
                                              electron_shell );

    if( electron_shell > 0 )
    {
      // Convert subshell number to endf number
      endf_subshell =
        Data::convertEADLDesignatorToENDFDesignator( electron_shell );

      // insert subshell to set
      endf_subshells.insert( endf_subshell );
    }

    // Read and process the data in the current table, then store in the HDF5
    // file
    switch( reaction_type )
    {
      // Integrated elastic transport cross section data
      case 7000:
      {
        // Interpolation should always be LogLog = 5
        testInvariant( interpolation_flag == 5 )

        std::vector<double> energy_grid, cross_section;
        eedl_file_handler->processTwoColumnTable( energy_grid,
                                                  cross_section );

        data_container.setElasticEnergyGrid( energy_grid );
        data_container.setElasticTransportCrossSection( cross_section );

        FRENSIE_LOG_PARTIAL_NOTIFICATION( "." );
        FRENSIE_FLUSH_ALL_LOGS();
        
        break;
      }

      // Integrated large angle scattering cross section data
      case 8000:
      {
        // Interpolation should always be LogLog = 5
        testInvariant( interpolation_flag == 5 )

        std::vector<double> energy_grid, cross_section;
        eedl_file_handler->processTwoColumnTable( energy_grid,
                                                  cross_section );

        // Test that the cutoff energy grid is the same as the transport
        testInvariant( energy_grid.size() ==
                       data_container.getElasticEnergyGrid().size() );

        data_container.setCutoffElasticCrossSection( cross_section );

        FRENSIE_LOG_PARTIAL_NOTIFICATION( "." );
        FRENSIE_FLUSH_ALL_LOGS();
        
        break;
      }

      // Average energy to residual atom from elastic scattering
      case 8011:
      {
        // Interpolation should always be LinLin = 0
        testInvariant( interpolation_flag == 0 )

        std::vector<double> residual_incident_energy, residual_energy;
        eedl_file_handler->processTwoColumnTable( residual_incident_energy,
                                                  residual_energy );

        data_container.setCutoffElasticResidualIncidentEnergy( residual_incident_energy );
        data_container.setCutoffElasticResidualEnergy( residual_energy );

        FRENSIE_LOG_PARTIAL_NOTIFICATION( "." );
        FRENSIE_FLUSH_ALL_LOGS();
        
        break;
      }

      // Average energy of scattered electron from elastic scattering
      case 8010:
      {
        // Interpolation should always be LinLin = 0
        testInvariant( interpolation_flag == 0 )

        std::vector<double> scattered_incident_energy, scattered_energy;
        eedl_file_handler->processTwoColumnTable( scattered_incident_energy,
                                                  scattered_energy );

        data_container.setCutoffElasticScatteredElectronIncidentEnergy(
          scattered_incident_energy );
        data_container.setCutoffElasticScatteredElectronEnergy(
          scattered_energy );

        FRENSIE_LOG_PARTIAL_NOTIFICATION( "." );
        FRENSIE_FLUSH_ALL_LOGS();
        
        break;
      }

      // Elastic angular distribution of the scattered electron data
      case 8022:
      {
        // Interpolation should always be LinLin = 0
        testInvariant( interpolation_flag == 0 )

        std::map<double,std::vector<double> > elastic_angle;

        eedl_file_handler->mapThreeColumnTable( elastic_angular_energy_grid,
                                                elastic_angle,
                                                elastic_pdf );

        data_container.setCutoffElasticAngularEnergyGrid(
          elastic_angular_energy_grid );
        data_container.setCutoffElasticAngles( elastic_angle );
        data_container.setCutoffElasticPDF( elastic_pdf );

        FRENSIE_LOG_PARTIAL_NOTIFICATION( "." );
        FRENSIE_FLUSH_ALL_LOGS();
        
        break;
      }

      // Integrated total elastic cross section data
      case 10000:
      {
        // Interpolation should always be LogLog = 5
        testInvariant( interpolation_flag == 5 )

        std::vector<double> energy_grid, cross_section;
        eedl_file_handler->processTwoColumnTable( energy_grid,
                                                  cross_section );

        // Test that the energy grid is the same as the transport and cutoff
        testPostcondition( energy_grid.size() ==
                           data_container.getElasticEnergyGrid().size() );

        data_container.setTotalElasticCrossSection( cross_section );

        FRENSIE_LOG_PARTIAL_NOTIFICATION( "." );
        FRENSIE_FLUSH_ALL_LOGS();
        
        break;
      }

      // Extract the integrated ionization (electroionization) cross section
      case 81000:
      {
        // Interpolation should always be LinLin = 0
        testInvariant( interpolation_flag == 0 )

        std::vector<double> energy_grid, cross_section;
        eedl_file_handler->processTwoColumnTable( energy_grid,
                                                  cross_section );

        data_container.setElectroionizationCrossSectionEnergyGrid(
          endf_subshell,
          energy_grid );
        data_container.setElectroionizationCrossSection(
          endf_subshell,
          cross_section );

        FRENSIE_LOG_PARTIAL_NOTIFICATION( "." );
        FRENSIE_FLUSH_ALL_LOGS();
        
        break;
      }

      // Average energy of primary and secondary electrons from ionization
      case 81010:
      {
        // Interpolation should always be LinLin = 0
        testInvariant( interpolation_flag == 0 )

        std::vector<double> incident_energy, average_outgoing_energy;
        eedl_file_handler->processTwoColumnTable( incident_energy,
                                                  average_outgoing_energy );

        // Average energy of electron from ionization
        if( outgoing_particle_designator == 9 ) // 9 = electron
        {
          data_container.setElectroionizationAverageScatteredElectronIncidentEnergy(
              endf_subshell,
              incident_energy );
          data_container.setElectroionizationAverageScatteredElectronEnergy(
              endf_subshell,
              average_outgoing_energy );
        }
        // Average energy of secondary electron from ionization
        else if( outgoing_particle_designator == 19 ) // 19 = recoil electron
        {
          data_container.setElectroionizationAverageRecoilElectronIncidentEnergy(
              endf_subshell,
              incident_energy );
          data_container.setElectroionizationAverageRecoilElectronEnergy(
              endf_subshell,
              average_outgoing_energy );
        }
        else
        {
          THROW_EXCEPTION( std::runtime_error,
                           "Unsupported outgoing particle designator "
                           "encountered in EEDL table (Yo = "
                           << outgoing_particle_designator << ")!" );
        }

        FRENSIE_LOG_PARTIAL_NOTIFICATION( "." );
        FRENSIE_FLUSH_ALL_LOGS();
        
        break;
      }
      
      case 81021:
      {
        // Interpolation should always be LinLin = 0
        testInvariant( interpolation_flag == 0 )
        // The outgoing particle designator should be electron as recoil (19)
        testInvariant( outgoing_particle_designator == 19 );

        std::vector<double> electron_energy_grid;
        std::map<double,std::vector<double> > electroionization_recoil_energy,
                                              electroionization_recoil_pdf;

        // Read the recoil electron spectrum from ionization for a subshell
        // If electron_shell == 0 then no subshell data only total
        eedl_file_handler->mapThreeColumnTable( electron_energy_grid,
                                                electroionization_recoil_energy,
                                                electroionization_recoil_pdf );

        data_container.setElectroionizationRecoilEnergyGrid(
                          endf_subshell,
                          electron_energy_grid );

        data_container.setElectroionizationRecoilEnergy(
                          endf_subshell,
                          electroionization_recoil_energy );

        data_container.setElectroionizationRecoilPDF(
                          endf_subshell,
                          electroionization_recoil_pdf );

        FRENSIE_LOG_PARTIAL_NOTIFICATION( "." );
        FRENSIE_FLUSH_ALL_LOGS();
        
        break;
      }
      
      case 82000:
      {
        // Interpolation should always be LinLin = 0
        testInvariant( interpolation_flag == 0 )

        // Extract the integrated bremsstrahlung cross section
        std::vector<double> energy_grid, cross_section;
        eedl_file_handler->processTwoColumnTable( energy_grid,
                                                  cross_section );

        data_container.setBremsstrahlungCrossSectionEnergyGrid( energy_grid );
        data_container.setBremsstrahlungCrossSection( cross_section );

        FRENSIE_LOG_PARTIAL_NOTIFICATION( "." );
        FRENSIE_FLUSH_ALL_LOGS();
        
        break;
      }

      // Average energy of secondary particles from bremsstrahlung
      case 82010:
      {
        // Interpolation should always be LinLin = 0
        testInvariant( interpolation_flag == 0 )

        std::vector<double> incident_energy, average_outgoing_energy;
        eedl_file_handler->processTwoColumnTable( incident_energy,
                                                  average_outgoing_energy );

        // Average energy of secondary photon from bremsstrahlung
        if ( outgoing_particle_designator == 7 ) // 7 = photon
        {
          data_container.setBremsstrahlungAveragePhotonIncidentEnergy(
              incident_energy );
          data_container.setBremsstrahlungAveragePhotonEnergy(
              average_outgoing_energy );
        }
        // Average energy of secondary electron from bremsstrahlung
        else if( outgoing_particle_designator == 9 )
        {
          data_container.setBremsstrahlungAverageElectronIncidentEnergy(
              incident_energy );
          data_container.setBremsstrahlungAverageElectronEnergy(
              average_outgoing_energy );
        }
        else
        {
          THROW_EXCEPTION( std::runtime_error,
                           "Unsupported outgoing particle designator "
                           "encountered in EEDL table (Yo = "
                           << outgoing_particle_designator << ")!" );
        }

        FRENSIE_LOG_PARTIAL_NOTIFICATION( "." );
        FRENSIE_FLUSH_ALL_LOGS();
        
        break;
      }

      // Read the spectrum of the secondary photon from bremsstrahlung
      case 82021:
      {
        // Interpolation should always be LinLin = 0
        testInvariant( interpolation_flag == 0 )
        // The outgoing particle designator should be photon (7)
        testInvariant( outgoing_particle_designator == 7 );

        std::vector<double> electron_energy_grid;
        std::map<double,std::vector<double> > bremsstrahlung_photon_energy,
                                              bremsstrahlung_photon_pdf;

        eedl_file_handler->mapThreeColumnTable( electron_energy_grid,
                                                bremsstrahlung_photon_energy,
                                                bremsstrahlung_photon_pdf );

        data_container.setBremsstrahlungPhotonEnergyGrid( electron_energy_grid );

        data_container.setBremsstrahlungPhotonEnergy( bremsstrahlung_photon_energy );

        data_container.setBremsstrahlungPhotonPDF( bremsstrahlung_photon_pdf );

        FRENSIE_LOG_PARTIAL_NOTIFICATION( "." );
        FRENSIE_FLUSH_ALL_LOGS();
        
        break;
      }

      // Extract the integrated (atomic) excitation cross section
      case 83000:
      {
        // Interpolation should always be LinLin = 0
        testInvariant( interpolation_flag == 0 )

        std::vector<double> energy_grid, cross_section;
        eedl_file_handler->processTwoColumnTable( energy_grid, cross_section );

        data_container.setAtomicExcitationEnergyGrid( energy_grid );
        data_container.setAtomicExcitationCrossSection( cross_section );

        FRENSIE_LOG_PARTIAL_NOTIFICATION( "." );
        FRENSIE_FLUSH_ALL_LOGS();
        
        break;
      }
      
      case 83011:
      {
        // Read the average energy loss from excitation
        testInvariant( interpolation_flag == 0 );

        std::vector<double> atomic_excitation_energy_grid,
                            atomic_excitation_energy_loss;

        eedl_file_handler->processTwoColumnTable( atomic_excitation_energy_grid,
                                                  atomic_excitation_energy_loss );

        testPostcondition( atomic_excitation_energy_grid.size() ==
                           data_container.getAtomicExcitationEnergyGrid().size() );

        data_container.setAtomicExcitationEnergyLoss( atomic_excitation_energy_loss );

        FRENSIE_LOG_PARTIAL_NOTIFICATION( "." );
        FRENSIE_FLUSH_ALL_LOGS();
        
        break;
      }

      // Unknown reaction type found
      default:
      {
        THROW_EXCEPTION( std::runtime_error,
                         "An unknown reaction type was encountered while "
                         "processing the EPDL file!" );
      }
    }
  }

  // Close the EEDL file
  eedl_file_handler->closeENDLFile();

  // Set the subshells
  if( data_container.getSubshells().empty() )
    data_container.setSubshells( endf_subshells );

  FRENSIE_LOG_NOTIFICATION( " done." );
}

} // end DataGen namespace

//---------------------------------------------------------------------------//
// end DataGen_StandardENDLDataGenerator.cpp
//---------------------------------------------------------------------------//
