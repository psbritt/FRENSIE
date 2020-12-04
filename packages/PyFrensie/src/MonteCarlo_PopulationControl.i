//---------------------------------------------------------------------------//
//!
//! \file   MonteCarlo_PopulationControl.i
//! \author Philip Britt
//! \brief  The population control class's interface file
//!
//---------------------------------------------------------------------------//

%{
// FRENSIE Includes
#include "PyFrensie_PythonTypeTraits.hpp"

#include "Utility_Mesh.hpp"
#include "Utility_StructuredHexMesh.hpp"
#include "Utility_TetMesh.hpp"

#include "MonteCarlo_ParticleHistoryObserver.hpp"
#include "MonteCarlo_DiscretizableParticleHistoryObserver.hpp"
#include "MonteCarlo_PopulationControl.hpp"
#include "MonteCarlo_WeightWindow.hpp"
#include "MonteCarlo_WeightWindowMesh.hpp"
#include "MonteCarlo_Importance.hpp"
#include "MonteCarlo_ImportanceMesh.hpp"

using namespace MonteCarlo;
%}

// C++ STL support
%include <stl.i>
%include <std_shared_ptr.i>
%include <std_map.i>
%include <std_vector.i>

// Include typemaps support
%include <typemaps.i>

// Mesh handling
%import "Utility.Mesh.i"
%include "Utility_StructuredHexMesh.hpp"
%include "Utility_TetMesh.hpp"

%shared_ptr(MonteCarlo::PopulationControl)
%include "MonteCarlo_PopulationControl.hpp"

%shared_ptr(MonteCarlo::Importance)
%include "MonteCarlo_Importance.hpp"

%typemap(in,numinputs=0) std::unordered_map<Utility::Mesh::ElementHandle,std::vector<double> >& importance_map (std::unordered_map<Utility::Mesh::ElementHandle,std::vector<double> > temp) "$1 = &temp;"

%shared_ptr(MonteCarlo::ImportanceMesh)
%include "MonteCarlo_ImportanceMesh.hpp"

%shared_ptr(MonteCarlo::WeightWindowBase)
%shared_ptr(MonteCarlo::WeightWindow)
%include "MonteCarlo_WeightWindow.hpp"

%typemap(in,numinputs=0) std::unordered_map<Utility::Mesh::ElementHandle,std::vector<MonteCarlo::WeightWindow> >& importance_map (std::unordered_map<Utility::Mesh::ElementHandle,std::vector<MonteCarlo::WeightWindow> > temp) "$1 = &temp;"

%shared_ptr(MonteCarlo::WeightWindowMesh)
%include "MonteCarlo_WeightWindowMesh.hpp"

//---------------------------------------------------------------------------//
// end MonteCarlo_PopulationControl.i
//---------------------------------------------------------------------------//