//---------------------------------------------------------------------------//
//!
//! \file   Utility.DirectionDiscretization.i
//! \author Philip Britt
//! \brief  The Utility.DirectionDiscretization sub-module swig interface file
//!
//---------------------------------------------------------------------------//

%define %utility_direction_discretization_docstring
"
PyFrensie.Utility.DirectionDiscretization is the python interface to the FRENSIE utility/direction_discretization
subpackage.

The purpose of the direction_discretization package is to provide methods of discretizing the unit sphere that
direction phase space exists on.
"
%enddef

%module(package   = "PyFrensie.Utility",
        autodoc   = "1",
        docstring = %utility_direction_discretization_docstring) DirectionDiscretization

%{
// FRENSIE Includes
#include "PyFrensie_PythonTypeTraits.hpp"
#include "Utility_SerializationHelpers.hpp"

#include "Utility_PQLAQuadrature.hpp"
// Add the Utility namespace to the global lookup scope
using namespace Utility;
%}

// Standard exception handling
%include "exception.i"

// Global swig features
%feature("autodoc", "1");

%include <std_shared_ptr.i>
%include <std_array.i>

// Include typemaps support
%include <typemaps.i>

// Include the serialization helpers for handling macros
%include "Utility_SerializationHelpers.hpp"

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

// Add typemaps for converting std::array to and from Python array
%typemap(in) const std::array<double,3>{
  $1 = PyFrensie::convertFromPython<std::array<double,3> >( $input );
}

%ignore *::getSphericalTriangleVector();
%ignore *::sampleIsotropicallyFromTriangle();


%template(DirectionArray) std::array<double, 3>;

%shared_ptr(Utility::PQLAQuadrature)
%include "Utility_PQLAQuadrature.hpp"

//---------------------------------------------------------------------------//
// Turn off the exception handling
//---------------------------------------------------------------------------//
%exception;

//---------------------------------------------------------------------------//
// end MonteCarlo.DirectionDiscretization.i
//---------------------------------------------------------------------------//