//---------------------------------------------------------------------------//
//!
//! \file   Utility_OArchivableObject_def.hpp
//! \author Alex Robinson
//! \brief  OArchivable object base class definition
//!
//---------------------------------------------------------------------------//

#ifndef UTILITY_OARCHIVABLE_OBJECT_DEF_HPP
#define UTILITY_OARCHIVABLE_OBJECT_DEF_HPP

// Std Lib Includes
#include <fstream>

// Boost Includes
#include <boost/filesystem.hpp>
#include <boost/serialization/nvp.hpp>

// FRENSIE Includes
#include "FRENSIE_Archives.hpp" // Must include first
#include "Utility_ExceptionTestMacros.hpp"
#include "Utility_ExceptionCatchMacros.hpp"

namespace Utility{

// Archive the object
/*! \details The file extension will be used to determine the archive type
 * (e.g. .xml, .txt, .bin, .h5fa)
 */
template<typename DerivedType>
void OArchivableObject<DerivedType>::saveToFile(
                         const boost::filesystem::path& archive_name_with_path,
                         const bool overwrite ) const
{
  // Check if overwriting is permitted
  if( !overwrite )
  {
    TEST_FOR_EXCEPTION( boost::filesystem::exists( archive_name_with_path ),
                        std::runtime_error,
                        "A file with the specified path and name already "
                        "exists!" );
  }

  // Verify that the parent directory exists
  if( archive_name_with_path.has_parent_path() )
  {
    TEST_FOR_EXCEPTION( !boost::filesystem::exists( archive_name_with_path.parent_path() ),
                        std::runtime_error,
                        "Cannot create the output archive "
                        << archive_name_with_path.string() <<
                        " because the parent directory does not exist!" );
  }
  
  // Initialize the archive ostream here to ensure that it gets deleted
  // after the archive
  std::unique_ptr<std::ostream> oarchive_stream;

  // Get the file extension
  std::string extension = archive_name_with_path.extension().string();

  // Create the oarchive
  if( extension == ".xml" || extension == ".txt" )
  {
    // Create the oarchive file stream
    oarchive_stream.reset(
                        new std::ofstream( archive_name_with_path.string() ) );
    
    if( extension == ".xml" )
    {
      boost::archive::xml_oarchive archive( *oarchive_stream );

      this->saveToArchive( archive );
    }
    else
    {
      boost::archive::text_oarchive archive( *oarchive_stream );

      this->saveToArchive( archive );
    }
  }
  else if( extension == ".bin" )
  {
    // Create the oarchive file stream
    oarchive_stream.reset( new std::ofstream( archive_name_with_path.string(),
                                              std::ofstream::binary ) );

    boost::archive::binary_oarchive archive( *oarchive_stream );

    this->saveToArchive( archive );
  }
  else if( extension == ".h5fa" )
  {
    Utility::HDF5OArchive archive( archive_name_with_path.string(),
                                   Utility::HDF5OArchiveFlags::OVERWRITE_EXISTING_ARCHIVE );

    this->saveToArchive( archive );
  }
  else
  {
    THROW_EXCEPTION( std::runtime_error,
                     "Cannot create the output archive "
                     << archive_name_with_path.string() <<
                     " because the extension type is not supported!" );
  }
}

// Archive the object using the required archive
template<typename DerivedType>
template<typename Archive>
void OArchivableObject<DerivedType>::saveToArchive( Archive& archive ) const
{
  try{
    archive << boost::serialization::make_nvp( this->getOArchiveName(), *dynamic_cast<const DerivedType*>(this) );
  }
  EXCEPTION_CATCH_RETHROW_AS( std::exception,
                              std::runtime_error,
                              "Unable to save the object to the desired "
                              "archive!" );
}
  
} // end Utility namespace

#endif // end UTILITY_OARCHIVABLE_OBJECT_DEF_HPP

//---------------------------------------------------------------------------//
// end Utility_OArchivableObject_def.hpp
//---------------------------------------------------------------------------//
