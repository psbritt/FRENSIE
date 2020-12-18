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
#include "MonteCarlo_WeightImportance.hpp"
#include "MonteCarlo_WeightImportanceMesh.hpp"
#include "MonteCarlo_Importance.hpp"
#include "MonteCarlo_ImportanceMesh.hpp"

%}

// C++ STL support
%include <stl.i>
%include <std_shared_ptr.i>
%include <std_unordered_map.i>
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

%template(ImportanceMap) std::unordered_map<Utility::Mesh::ElementHandle, std::vector<double>>;

%shared_ptr(MonteCarlo::ImportanceMesh)
%include "MonteCarlo_ImportanceMesh.hpp"

%shared_ptr(MonteCarlo::WeightImportance)
%include "MonteCarlo_WeightImportance.hpp"

/* Technically already defined above, but redefining it so a separate python class exists that has a 
 * recognizable name is formed
 */

%shared_ptr(MonteCarlo::WeightImportanceMesh)
%include "MonteCarlo_WeightImportanceMesh.hpp"

%shared_ptr(MonteCarlo::WeightWindowBase)
%include "MonteCarlo_WeightWindow.hpp"

%shared_ptr(MonteCarlo::WeightWindowMesh)
%include "MonteCarlo_WeightWindowMesh.hpp"

%template(WeightWindowVector) std::vector<MonteCarlo::WeightWindow>;
%template(WeightWindowMap) std::unordered_map<Utility::Mesh::ElementHandle, std::vector<MonteCarlo::WeightWindow>>;


//---------------------------------------------------------------------------//
// end MonteCarlo_PopulationControl.i
//---------------------------------------------------------------------------//