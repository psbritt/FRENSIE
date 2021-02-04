//---------------------------------------------------------------------------//
//!
//! \file   MonteCarlo_ParticleDistribution.i
//! \author Luke Kersting
//! \brief  The ParticleDistribution classes interface file
//!
//---------------------------------------------------------------------------//

%{
// FRENSIE Includes
#include "PyFrensie_PythonTypeTraits.hpp"
#include "MonteCarlo_ParticleDistribution.hpp"
#include "MonteCarlo_StandardParticleDistribution.hpp"
#include "MonteCarlo_GenericHistogramImportanceParticleDistribution.hpp"

#include "Utility_ToStringTraits.hpp"
#include "Utility_ToStringTraitsDecl.hpp"

using namespace MonteCarlo;
%}

// Import the ToStringTraitsDecl
%import "Utility_ToStringTraitsDecl.hpp"

// Ignore certian functions that use DimensionCounterMap
%ignore *::initializeDimensionCounters;
%ignore *::sampleAndRecordTrials;
%ignore *::sampleWithDimensionValueAndRecordTrials;


// ---------------------------------------------------------------------------//
// Add ParticleDistribution support
// ---------------------------------------------------------------------------//

%shared_ptr( MonteCarlo::ParticleDistribution )
%include "MonteCarlo_ParticleDistribution.hpp"

// ---------------------------------------------------------------------------//
// Add StandardParticleDistribution support
// ---------------------------------------------------------------------------//

%shared_ptr( MonteCarlo::StandardParticleDistribution )
%include "MonteCarlo_StandardParticleDistribution.hpp"

// ---------------------------------------------------------------------------//
// Add GenericHistogramImportanceParticleDistribution support
// ---------------------------------------------------------------------------//

%template(dimensionOrderArray) std::vector<MonteCarlo::PhaseSpaceDimension>;
%template(importanceDistributionBoundaryMap) std::map< MonteCarlo::PhaseSpaceDimension, std::vector< double > >;
%template(importanceDistributionPointerVector) std::vector< std::shared_ptr< MonteCarlo::PhaseSpaceDimensionDistribution> >;
%template(importanceDistributionDimensionMap) std::map< MonteCarlo::PhaseSpaceDimension, std::vector< std::shared_ptr< MonteCarlo::PhaseSpaceDimensionDistribution> > >;
%shared_ptr( MonteCarlo::GenericHistogramImportanceParticleDistribution )
%include "MonteCarlo_GenericHistogramImportanceParticleDistribution.hpp"

//---------------------------------------------------------------------------//
// end MonteCarlo_ParticleDistribution.i
//---------------------------------------------------------------------------//
