//---------------------------------------------------------------------------//
//!
//! \file   MonteCarlo_ImportanceSampledIndependentPhaseSpaceDimensionDistribution.cpp
//! \author Alex Robinson
//! \brief  Importance sampled independent phase space dimension distribution
//!         class template instantiations
//!
//---------------------------------------------------------------------------//

// FRENSIE Includes
#include "FRENSIE_Archives.hpp"
#include "MonteCarlo_ImportanceSampledIndependentPhaseSpaceDimensionDistribution.hpp"

BOOST_CLASS_EXPORT_IMPLEMENT( MonteCarlo::ImportanceSampledIndependentPrimarySpatialDimensionDistribution );
EXPLICIT_TEMPLATE_CLASS_INST( MonteCarlo::ImportanceSampledIndependentPhaseSpaceDimensionDistribution<MonteCarlo::PRIMARY_SPATIAL_DIMENSION> );
EXPLICIT_CLASS_SAVE_LOAD_INST( MonteCarlo::ImportanceSampledIndependentPhaseSpaceDimensionDistribution<MonteCarlo::PRIMARY_SPATIAL_DIMENSION> );

BOOST_CLASS_EXPORT_IMPLEMENT( MonteCarlo::ImportanceSampledIndependentSecondarySpatialDimensionDistribution );
EXPLICIT_TEMPLATE_CLASS_INST( MonteCarlo::ImportanceSampledIndependentPhaseSpaceDimensionDistribution<MonteCarlo::SECONDARY_SPATIAL_DIMENSION> );
EXPLICIT_CLASS_SAVE_LOAD_INST( MonteCarlo::ImportanceSampledIndependentPhaseSpaceDimensionDistribution<MonteCarlo::SECONDARY_SPATIAL_DIMENSION> );

BOOST_CLASS_EXPORT_IMPLEMENT( MonteCarlo::ImportanceSampledIndependentTertiarySpatialDimensionDistribution );
EXPLICIT_TEMPLATE_CLASS_INST( MonteCarlo::ImportanceSampledIndependentPhaseSpaceDimensionDistribution<MonteCarlo::TERTIARY_SPATIAL_DIMENSION> );
EXPLICIT_CLASS_SAVE_LOAD_INST( MonteCarlo::ImportanceSampledIndependentPhaseSpaceDimensionDistribution<MonteCarlo::TERTIARY_SPATIAL_DIMENSION> );

BOOST_CLASS_EXPORT_IMPLEMENT( MonteCarlo::ImportanceSampledIndependentPrimaryDirectionalDimensionDistribution );
EXPLICIT_TEMPLATE_CLASS_INST( MonteCarlo::ImportanceSampledIndependentPhaseSpaceDimensionDistribution<MonteCarlo::PRIMARY_DIRECTIONAL_DIMENSION> );
EXPLICIT_CLASS_SAVE_LOAD_INST( MonteCarlo::ImportanceSampledIndependentPhaseSpaceDimensionDistribution<MonteCarlo::PRIMARY_DIRECTIONAL_DIMENSION> );

BOOST_CLASS_EXPORT_IMPLEMENT( MonteCarlo::ImportanceSampledIndependentSecondaryDirectionalDimensionDistribution );
EXPLICIT_TEMPLATE_CLASS_INST( MonteCarlo::ImportanceSampledIndependentPhaseSpaceDimensionDistribution<MonteCarlo::SECONDARY_DIRECTIONAL_DIMENSION> );
EXPLICIT_CLASS_SAVE_LOAD_INST( MonteCarlo::ImportanceSampledIndependentPhaseSpaceDimensionDistribution<MonteCarlo::SECONDARY_DIRECTIONAL_DIMENSION> );

BOOST_CLASS_EXPORT_IMPLEMENT( MonteCarlo::ImportanceSampledIndependentTertiaryDirectionalDimensionDistribution );
EXPLICIT_TEMPLATE_CLASS_INST( MonteCarlo::ImportanceSampledIndependentPhaseSpaceDimensionDistribution<MonteCarlo::TERTIARY_DIRECTIONAL_DIMENSION> );
EXPLICIT_CLASS_SAVE_LOAD_INST( MonteCarlo::ImportanceSampledIndependentPhaseSpaceDimensionDistribution<MonteCarlo::TERTIARY_DIRECTIONAL_DIMENSION> );

BOOST_CLASS_EXPORT_IMPLEMENT( MonteCarlo::ImportanceSampledIndependentEnergyDimensionDistribution );
EXPLICIT_TEMPLATE_CLASS_INST( MonteCarlo::ImportanceSampledIndependentPhaseSpaceDimensionDistribution<MonteCarlo::ENERGY_DIMENSION> );
EXPLICIT_CLASS_SAVE_LOAD_INST( MonteCarlo::ImportanceSampledIndependentPhaseSpaceDimensionDistribution<MonteCarlo::ENERGY_DIMENSION> );

BOOST_CLASS_EXPORT_IMPLEMENT( MonteCarlo::ImportanceSampledIndependentTimeDimensionDistribution );
EXPLICIT_TEMPLATE_CLASS_INST( MonteCarlo::ImportanceSampledIndependentPhaseSpaceDimensionDistribution<MonteCarlo::TIME_DIMENSION> );
EXPLICIT_CLASS_SAVE_LOAD_INST( MonteCarlo::ImportanceSampledIndependentPhaseSpaceDimensionDistribution<MonteCarlo::TIME_DIMENSION> );

BOOST_CLASS_EXPORT_IMPLEMENT( MonteCarlo::ImportanceSampledIndependentWeightDimensionDistribution );
EXPLICIT_TEMPLATE_CLASS_INST( MonteCarlo::ImportanceSampledIndependentPhaseSpaceDimensionDistribution<MonteCarlo::WEIGHT_DIMENSION> );
EXPLICIT_CLASS_SAVE_LOAD_INST( MonteCarlo::ImportanceSampledIndependentPhaseSpaceDimensionDistribution<MonteCarlo::WEIGHT_DIMENSION> );

BOOST_CLASS_EXPORT_IMPLEMENT( MonteCarlo::ImportanceSampledIndependentSpatialIndexDimensionDistribution );
EXPLICIT_TEMPLATE_CLASS_INST( MonteCarlo::ImportanceSampledIndependentPhaseSpaceDimensionDistribution<MonteCarlo::SPATIAL_INDEX_DIMENSION> );
EXPLICIT_CLASS_SAVE_LOAD_INST( MonteCarlo::ImportanceSampledIndependentPhaseSpaceDimensionDistribution<MonteCarlo::SPATIAL_INDEX_DIMENSION> );

BOOST_CLASS_EXPORT_IMPLEMENT( MonteCarlo::ImportanceSampledIndependentDirectionIndexDimensionDistribution );
EXPLICIT_TEMPLATE_CLASS_INST( MonteCarlo::ImportanceSampledIndependentPhaseSpaceDimensionDistribution<MonteCarlo::DIRECTION_INDEX_DIMENSION> );
EXPLICIT_CLASS_SAVE_LOAD_INST( MonteCarlo::ImportanceSampledIndependentPhaseSpaceDimensionDistribution<MonteCarlo::DIRECTION_INDEX_DIMENSION> );

//---------------------------------------------------------------------------//
// end MonteCarlo_ImportanceSampledIndependentPhaseSpaceDimensionDistribution.cpp
//---------------------------------------------------------------------------//
