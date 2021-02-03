//---------------------------------------------------------------------------//
//!
//! \file   MonteCarlo_GenericHistogramImportanceParticleDistribution.hpp
//! \author Philip Britt
//! \brief  Generic Importance particle distribution declaration
//!
//---------------------------------------------------------------------------//

#ifndef MONTE_CARLO_GENERIC_HISTOGRAM_IMPORTANCE_PARTICLE_DISTRIBUTION_HPP
#define MONTE_CARLO_GENERIC_HISTOGRAM_IMPORTANCE_PARTICLE_DISTRIBUTION_HPP

// Std Lib Includes
#include <memory>
#include <functional>

// FRENSIE Includes
#include "MonteCarlo_ParticleDistribution.hpp"
#include "Utility_SpatialCoordinateConversionPolicy.hpp"
#include "Utility_DirectionalCoordinateConversionPolicy.hpp"
#include "Utility_QuantityTraits.hpp"
#include "Utility_Map.hpp"
#include "Utility_Set.hpp"
#include "Utility_Vector.hpp"
#include "Utility_StructuredHexMesh.hpp"
#include "Utility_PQLAQuadrature.hpp"

namespace MonteCarlo{

//! The standard particle distribution class
class GenericHistogramImportanceParticleDistribution : public ParticleDistribution
{
  // Typedef for scalar traits
  typedef Utility::QuantityTraits<double> QT;

  // Typedef for the dimension sampling function
  typedef std::function<void(const PhaseSpaceDimensionDistribution&,PhaseSpacePoint&)> DimensionSamplingFunction;

public:

  //! Typedef for the phase space dimension set
  typedef ParticleDistribution::DimensionSet DimensionSet;

  //! The dimension trial counter map
  typedef ParticleDistribution::DimensionCounterMap DimensionCounterMap;

  //! Basic Constructor
  GenericHistogramImportanceParticleDistribution( const std::string& name );

  //! Destructor
  ~GenericHistogramImportanceParticleDistribution()
  { /* ... */ }

  //! Set a dimension distribution
  void setIndependentDimensionDistribution(
                        const std::shared_ptr< PhaseSpaceDimensionDistribution >& dimension_distribution,
                        const bool finished_with_independent);

  void setImportanceDimensionDistributions(                      
                        const std::map< PhaseSpaceDimension, std::vector< std::shared_ptr< PhaseSpaceDimensionDistribution> > > importance_sampled_distributions,
                        const std::map< PhaseSpaceDimension, std::vector< double > > importance_distribution_boundaries,
                        const std::vector< PhaseSpaceDimension > importance_distribution_order);

  //! Use these BEFORE you set a distribution for them
  void setMeshIndexDimensionDistributionObject( std::shared_ptr<  Utility::StructuredHexMesh> mesh_object );
  
  void setDirectionIndexDimensionDistributionObject( std::shared_ptr< Utility::PQLAQuadrature> quadrature_pointer );

  //! Return the dimension distribution type name
  std::string getDimensionDistributionTypeName(
                          const PhaseSpaceDimension dimension ) const override;

  //! Check if the distribution is spatially uniform (somewhere)
  bool isSpatiallyUniform() const override;

  //! Check if the distribution is directionally uniform (isotropic)
  bool isDirectionallyUniform() const override;

  //! Initialize dimension counter map
  void initializeDimensionCounters( DimensionCounterMap& trials ) const override;

  //! Evaluate the distribution at the desired phase space point
  double evaluate( const ParticleState& particle ) const override;

  //! Sample a particle state from the distribution
  void sample( ParticleState& particle ) const override;

  //! Sample a particle state from the dist. and record the number of trials
  void sampleAndRecordTrials( ParticleState& particle,
                              DimensionCounterMap& trials ) const override;

  //! Sample a particle state with the desired dimension value
  void sampleWithDimensionValue( ParticleState& particle,
                                 const PhaseSpaceDimension dimension,
                                 const double dimension_value ) const override;

  //! Sample a particle state with the desired dim. value and record trials
  void sampleWithDimensionValueAndRecordTrials(
                                 ParticleState& particle,
                                 DimensionCounterMap& trials,
                                 const PhaseSpaceDimension dimension,
                                 const double dimension_value ) const override;

protected:

  //! Default Constructor
  GenericHistogramImportanceParticleDistribution();

private:

  // Sample the particle state using the desired sampling functor
  template<typename DimensionSamplingFunctor>
  void sampleImpl( DimensionSamplingFunctor& dimension_sampling_function,
                   ParticleState& particle ) const;

  // Save the state to an archive
  template<typename Archive>
  void save( Archive& ar, const unsigned version ) const;

  // Load the data from an archive
  template<typename Archive>
  void load( Archive& ar, const unsigned version );

  BOOST_SERIALIZATION_SPLIT_MEMBER();

  // Declare the boost serialization access object as a friend
  friend class boost::serialization::access;

  // The spatial coordinate conversion policy
  std::shared_ptr< Utility::SpatialCoordinateConversionPolicy>
  d_spatial_coord_conversion_policy;

  // The directional coordinate conversion policy
  std::shared_ptr< Utility::DirectionalCoordinateConversionPolicy>
  d_directional_coord_conversion_policy;

  bool d_independent_finished;

  // The independent particle source dimensions
  std::vector< PhaseSpaceDimension > d_dimension_order;

  // The particle source dimensions
  std::map<PhaseSpaceDimension,std::vector< std::shared_ptr<PhaseSpaceDimensionDistribution> > > d_dimension_distributions;

  // The map of distribution boundaries
  std::map< PhaseSpaceDimension, std::vector< double > > d_dimension_bounds;

  // Cached values for archiving
  std::shared_ptr< Utility::StructuredHexMesh> d_spatial_dimension_mesh;

  std::shared_ptr< Utility::PQLAQuadrature> d_direction_dimension_discretization;
};

} // end MonteCarlo namespace

BOOST_SERIALIZATION_CLASS_VERSION( GenericHistogramImportanceParticleDistribution, MonteCarlo, 0 );
BOOST_SERIALIZATION_CLASS_EXPORT_STANDARD_KEY( GenericHistogramImportanceParticleDistribution, MonteCarlo );
EXTERN_EXPLICIT_CLASS_SAVE_LOAD_INST( MonteCarlo, GenericHistogramImportanceParticleDistribution );

//---------------------------------------------------------------------------//
// Template Includes
//---------------------------------------------------------------------------//

#include "MonteCarlo_GenericHistogramImportanceParticleDistribution_def.hpp"

//---------------------------------------------------------------------------//

#endif // end MONTE_CARLO_STANDARD_PARTICLE_DISTRIBUTION_HPP

//---------------------------------------------------------------------------//
// end MonteCarlo_GenericHistogramImportanceParticleDistribution.hpp
//---------------------------------------------------------------------------//
