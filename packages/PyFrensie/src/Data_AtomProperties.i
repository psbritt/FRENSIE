//---------------------------------------------------------------------------//
//!
//! \file   Data_AtomProperties.i
//! \author Luke Kersting
//! \brief  The AtomProperties classes interface file
//!
//---------------------------------------------------------------------------//

%{
// Std Lib Includes
#include <stdexcept>
#include <sstream>
#include <memory>

// FRENSIE Includes
#include "PyFrensie_PythonTypeTraits.hpp"
#include "Data_AtomType.hpp"
#include "Data_ScatteringCenterPropertiesHelper.hpp"
#include "Data_PhotoatomicDataProperties.hpp"
#include "Data_AdjointPhotoatomicDataProperties.hpp"
#include "Data_ElectroatomicDataProperties.hpp"
#include "Data_AdjointElectroatomicDataProperties.hpp"
#include "Data_AtomProperties.hpp"
// #include "Data_ZAID.hpp"
#include "Data_ExplicitTemplateInstantiationMacros.hpp"
// #include "Utility_Set.hpp"
#include "Utility_ContractException.hpp"
#include "Utility_SerializationHelpers.hpp"
#include "Utility_ToStringTraitsDecl.hpp"

// Add the Data namespace to the global lookup scope
using namespace Data;
%}

// C++ STL support
%include <stl.i>
%include <std_set.i>
%include <std_shared_ptr.i>
%include <std_except.i>

// Include typemaps support
%include <typemaps.i>

// Import the ToStringTraitsDecl
%import "Utility_ToStringTraitsDecl.hpp"

// Include the serialization helpers for handling macros
%include "Utility_SerializationHelpers.hpp"

// Include the explicit template instantiation helpers
%include "Data_ExplicitTemplateInstantiationMacros.hpp"

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

// General ignore directives
%ignore *::operator<<;

// Add general templates
%template(IntSet) std::set< unsigned int>;

//---------------------------------------------------------------------------//
// Add support for the AtomType
//---------------------------------------------------------------------------//
// Import the AtomType
%include "Data_AtomType.hpp"

//---------------------------------------------------------------------------//
// Add support for the PhotoatomicDataProperties
//---------------------------------------------------------------------------//

%atomic_properties_interface_setup( PhotoatomicDataProperties );

// Import the PhotoatomicDataProperties
%include "Data_PhotoatomicDataProperties.hpp"

//---------------------------------------------------------------------------//
// Add support for the AdjointPhotoatomicDataProperties
//---------------------------------------------------------------------------//

%atomic_properties_interface_setup( AdjointPhotoatomicDataProperties );

// Import the AdjointPhotoatomicDataProperties
%include "Data_AdjointPhotoatomicDataProperties.hpp"

//---------------------------------------------------------------------------//
// Add support for the ElectroatomicDataProperties
//---------------------------------------------------------------------------//

%atomic_properties_interface_setup( ElectroatomicDataProperties );

// Import the ElectroatomicDataProperties
%include "Data_ElectroatomicDataProperties.hpp"

//---------------------------------------------------------------------------//
// Add support for the AdjointElectroatomicDataProperties
//---------------------------------------------------------------------------//

%atomic_properties_interface_setup( AdjointElectroatomicDataProperties );

// Import the AdjointElectroatomicDataProperties
%include "Data_AdjointElectroatomicDataProperties.hpp"

//---------------------------------------------------------------------------//
// Add support for the AtomProperties
//---------------------------------------------------------------------------//

%ignore *::AtomicWeight;
%ignore Data::AtomProperties::AtomProperties( const Data::AtomType, const AtomicWeight );

%feature("docstring") Data::AtomProperties
"The AtomProperties class stores a atomic data properties. It can be used for
querying atomic data properties and for creating atomic data extractors or
container, which can be used to read atomic data."

%feature("autodoc", "atom(AtomProperties self) -> AtomType")
Data::AtomProperties::atom;

%feature("autodoc", "atomicNumber(AtomProperties self) -> unsigned")
Data::AtomProperties::atomicNumber;

%feature("autodoc", "atomicWeight(AtomProperties self) -> AtomicWeight")
Data::AtomProperties::atomicWeight;

%feature("autodoc", "atomicWeightRatio(AtomProperties self) -> double")
Data::AtomProperties::atomicWeightRatio;

// Allow shared pointers of AtomProperties objects
%shared_ptr( Data::AtomProperties );

// Rename the overloaded getDataFileVersions functions
%atom_properties_interface_setup( Photoatomic, photoatomic )
%atom_properties_interface_setup( AdjointPhotoatomic, adjointPhotoatomic )
%atom_properties_interface_setup( Electroatomic, electroatomic )
%atom_properties_interface_setup( AdjointElectroatomic, adjointElectroatomic )

// Add typemaps for converting AtomicWeight to and from Python float
%apply const Data::PhotoatomicDataProperties::AtomicWeight {
  const Data::AtomProperties::AtomicWeight }
%apply Data::PhotoatomicDataProperties::AtomicWeight {
  Data::AtomProperties::AtomicWeight }

// Import the AtomProperties
%include "Data_AtomProperties.hpp"

//---------------------------------------------------------------------------//
// end Data_AtomProperties.i
//---------------------------------------------------------------------------//e
