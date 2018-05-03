//---------------------------------------------------------------------------//
//!
//! \file   Data_ScatteringCenterPropertiesDatabase.i
//! \author Luke Kersting
//! \brief  The ScatteringCenterPropertiesDatabase classes interface file
//!
//---------------------------------------------------------------------------//

%{
// FRENSIE Includes
#include "PyFrensie_PythonTypeTraits.hpp"
#include "Data_ScatteringCenterPropertiesDatabaseImpl.hpp"
#include "Data_ScatteringCenterPropertiesDatabase.hpp"
#include "Data_ScatteringCenterPropertiesHelper.hpp"

// Add the Data namespace to the global lookup scope
using namespace Data;
%}

// C++ STL support
%include <stl.i>
%include <std_except.i>

// Include typemaps support
%include <typemaps.i>

// Include the data property helpers
%include "Data_PropertyHelpers.i"

// Standard exception handling
%include "exception.i"

// General exception handling
%exception
{
  try{
    $action;
    if( PyErr_Occurred() )
      SWIG_fail;
  }
  catch( Utility::ContractException& e )
  {
    SWIG_exception( SWIG_ValueError, e.what() );
  }
  catch( Data::InvalidScatteringCenterPropertiesData& e )
  {
    SWIG_exception( SWIG_RuntimeError, e.what() );
  }
  catch( std::runtime_error& e )
  {
    SWIG_exception( SWIG_RuntimeError, e.what() );
  }
  catch( ... )
  {
    SWIG_exception( SWIG_UnknownError, "Unknown C++ exception" );
  }
}

// Add some general typemaps

// Add typemaps for converting file_path to and from Python string
%typemap(in) const boost::filesystem::path& ( boost::filesystem::path temp ){
  temp = PyFrensie::convertFromPython<std::string>( $input );
  $1 = &temp;
}

%typemap(out) boost::filesystem::path {
  %append_output(PyFrensie::convertToPython( $1.string() ) );
}

%typemap(out) const boost::filesystem::path& {
  %append_output(PyFrensie::convertToPython( $1->string() ) );
}

%typemap(typecheck, precedence=1140) (const boost::filesystem::path&) {
  $1 = (PyString_Check($input)) ? 1 : 0;
}

//---------------------------------------------------------------------------//
// Add support for the ScatteringCenterPropertiesDatabaseImpl
//---------------------------------------------------------------------------//

%database_interface_setup(ScatteringCenterPropertiesDatabaseImpl)

// Include the ScatteringCenterPropertiesDatabaseImpl
%include "Data_ScatteringCenterPropertiesDatabaseImpl.hpp"

//---------------------------------------------------------------------------//
// Add support for the ScatteringCenterPropertiesDatabase
//---------------------------------------------------------------------------//

%extend Data::ScatteringCenterPropertiesDatabase
{
  void Data::ScatteringCenterPropertiesDatabase::saveToFile(
    const boost::filesystem::path& archive_name_with_path,
    const bool overwrite ) const
  {
    $self->saveToFile( archive_name_with_path, overwrite );

  }
}

%database_interface_setup(ScatteringCenterPropertiesDatabase)

// Include the ScatteringCenterPropertiesDatabase
%include "Data_ScatteringCenterPropertiesDatabase.hpp"

//---------------------------------------------------------------------------//
// end Data_ScatteringCenterPropertiesDatabase.i
//---------------------------------------------------------------------------//e
