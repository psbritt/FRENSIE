//---------------------------------------------------------------------------//
//!
//! \file   DataGen.ElectronPhoton.i
//! \author Luke Kersting
//! \brief  The DataGen.ElectronPhoton sub-module swig interface file
//!
//---------------------------------------------------------------------------//

%define %data_gen_electron_photon_docstring
"
PyFrensie.DataGen.ElectronPhoton is the python interface to the FRENSIE
data_gen/electron_photon subpackage.

The purpose of ElectronPhoton is to provide tools for generating Electron-Photon
data.
"
%enddef

%module(package   = "PyFrensie.DataGen",
        autodoc   = "1",
        docstring = %data_gen_electron_photon_docstring) ElectronPhoton

%{
// Std Lib Includes
#include <sstream>

// PyTrilinos Includes
#define NO_IMPORT_ARRAY
#include "numpy_include.hpp"

// FRENSIE Includes
#include "DataGen_EPRDataGeneratorHelpers.hpp"
#include "Utility_ArchivableObject.hpp"
#include "Utility_ContractException.hpp"
%}

// Include std::string support
%include <std_string.i>

// C++ STL support
%include <std_shared_ptr.i>
%include <std_vector.i>

// Include typemaps support
%include <typemaps.i>

/*// C++ STL support*/
%include <stl.i>
%include <std_string.i>
%include <std_except.i>
%include <std_set.i>
%include <std_pair.i>
/*%include <std_vector.i>*/


// Include the Teuchos::ArrayRCP support
%include "PyFrensie_Array.i"

// Standard exception handling
%include "exception.i"

// Global swig features
%feature("autodoc", "1");

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
  catch( std::runtime_error& e )
  {
    SWIG_exception( SWIG_RuntimeError, e.what() );
  }
  catch( ... )
  {
    SWIG_exception( SWIG_UnknownError, "Unknown C++ exception" );
  }
}

//---------------------------------------------------------------------------//
// Create aliases for common type found in native data tables
//---------------------------------------------------------------------------//

// Allow std::vector<double> output type
%template(DoubleVector) std::vector<double>;

// Add use of std::shared_ptr
%shared_ptr(DataGen::StandardMomentPreservingElectronDataGenerator)

// Ignore populateMomentPreservingData
%ignore *::populateMomentPreservingData( Data::MomentPreservingElectronVolatileDataContainer&, const int& ) const;
// Ignore setMomentPreservingElectronData
%ignore *::setMomentPreservingElectronData( Data::MomentPreservingElectronVolatileDataContainer&, const int& ) const;

//---------------------------------------------------------------------------//
// Add support for the StandardMomentPreservingElectronDataGenerator
//---------------------------------------------------------------------------//
// Add a more detailed docstring for the StandardMomentPreservingElectronDataGenerator
%feature("docstring")
DataGen::EPRDataGeneratorHelpers
"
The EPRDataGeneratorHelpers can be used to generate
/*elastic moment preserving data.*/
/*A brief usage tutorial for this class is shown below:*/

/*  import PyFrensie.DataGen.ElectronPhoton, PyTrilinos.Teuchos, numpy, matplotlib.pyplot*/

/*  source = PyTrilinos.Teuchos.FileInputSource( 'datadir/cross_sections.xml' )*/
/*  xml_obj = source.getObject()*/
/*  cs_list = PyTrilinos.Teuchos.XMLParameterListReader().toParameterList( xml_obj )*/

/*  h_data_list = cs_list.get( 'H-ENDL' )*/
/*  h_endl_file_name = 'datadir' + h_data_list.get( 'electroatomic_file_path' )*/

/*  h_native_data = PyFrensie.Data.ENDL.ENDLDataContainer( h_endl_file_name )*/

/*  matplotlib.pyplot.loglog( h_native_data.getElasticEnergyGrid(), h_native_data.getCutoffElasticCrossSection() )*/
/*  matplotlib.pyplot.loglog( h_native_data.getElasticEnergyGrid(), h_native_data.getElasticTransportCrossSection() )*/
/*  matplotlib.pyplot.show()*/
"

// Include the EPRDataGeneratorHelpers.hpp
%include "DataGen_EPRDataGeneratorHelpers.hpp"


//---------------------------------------------------------------------------//
// Add support for the StandardMomentPreservingElectronDataGenerator
//---------------------------------------------------------------------------//

%feature("docstring")
DataGen::StandardMomentPreservingElectronDataGenerator
"The StandardMomentPreservingElectronDataGenerator can be used to generate
elastic moment preserving data."

// Add a general typemap
%apply std::vector<double>& OUTPUT { std::vector<double>& discrete_angles };
%apply std::vector<double>& OUTPUT { std::vector<double>& weights };

%feature("autodoc",
"evaluateDiscreteAnglesAndWeights(StandardMomentPreservingElectronDataGenerator self, const double & energy, const unsigned int & number_of_discrete_angles ) -> DoubleVector, DoubleVector" )
DataGen::StandardMomentPreservingElectronDataGenerator::evaluateDiscreteAnglesAndWeights;

%include "DataGen_StandardMomentPreservingElectronDataGenerator.hpp"

/*  // Generate elastic discrete angle cosines and weights*/
/*  void evaluateDiscreteAnglesAndWeights(*/
/*    const double& energy,*/
/*    const int& number_of_discrete_angles,*/
/*    std::vector<double>& discrete_angles,*/
/*    std::vector<double>& weights ) const;*/



//---------------------------------------------------------------------------//
// Turn off the exception handling
//---------------------------------------------------------------------------//
%exception;

//---------------------------------------------------------------------------//
// end DataGen.ElectronPhoton.i
//---------------------------------------------------------------------------//