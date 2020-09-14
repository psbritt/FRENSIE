//---------------------------------------------------------------------------//
//!
//! \file   MonteCarlo_DiscretizedImportanceSampledParticleDistribution.hpp
//! \author Philip Britt
//! \brief  Discretized Importance Sampled Particle distribution class declaration (strictly for my thesis).
//!         Conversion policies are always cartesian, must have space, energy, and direction defined to use.
//!         This class assumes independent space, energy, and direction source dimensions (for now)
//!
//---------------------------------------------------------------------------//

#ifndef MONTE_CARLO_DISCRETIZED_IMPORTANCE_SAMPLED_PARTICLE_DISTRIBUTION_HPP
#define MONTE_CARLO_DISCRETIZED_IMPORTANCE_SAMPLED_PARTICLE_DISTRIBUTION_HPP

// Std Lib Includes
#include <iostream>
#include <memory>

// FRENSIE Includes
#include "MonteCarlo_PhaseSpaceDimensionDistribution.hpp"
#include "MonteCarlo_ParticleState.hpp"
#include "Utility_TypeNameTraits.hpp"
#include "Utility_Set.hpp"
#include "Utility_Mesh.hpp"
#include "Utility_PQLAQuadrature.hpp"
#include "MonteCarlo_ParticleDistribution.hpp"

namespace MonteCarlo{

//! The particle distribution base class
class DiscretizedImportanceSampledParticleDistribution: public ParticleDistribution
{

  //! Typedef for the phase space dimension distribution map
  typedef std::map<PhaseSpaceDimension,std::shared_ptr<PhaseSpaceDimensionDistribution> > DimensionDistributionMap;

public:

  //! The dimension trial counter map
  typedef ParticleDistribution::DimensionCounterMap DimensionCounterMap;

  // (for now only use cartesian coordinate policies for space and direction)

  //! Constructor (without specifics of which mesh elements are to be considered source)
  DiscretizedImportanceSampledParticleDistribution( const std::string& name,
                                                    const std::shared_ptr<const Utility::Mesh>& mesh,
                                                    const std::vector<double>& energy_bin_bounds,
                                                    const std::shared_ptr<const Utility::PQLAQuadrature>& direction_discretization);

  //! Constructor (with specifics of which mesh elements are to be considered source)
  DiscretizedImportanceSampledParticleDistribution(const std::string& name,
                                                   const std::shared_ptr<const Utility::Mesh>& mesh,
                                                   const std::vector<Utility::Mesh::ElementHandle>& mesh_source_elements,
                                                   const std::vector<double>& energy_bin_bounds,
                                                   const std::shared_ptr<const Utility::PQLAQuadrature>& direction_discretization);

  //! Destructor
  ~DiscretizedImportanceSampledParticleDistribution()
  { /* ... */ }

  //! Return the dimension distribution type name
  std::string getDimensionDistributionTypeName( const PhaseSpaceDimension dimension ) const override;

  //! Check if the distribution is spatially uniform (somewhere) (Only apply to primary PDF)
  bool isSpatiallyUniform() const override;

  //! Check if the distribution is directionally uniform (isotropic) (only apply to primary PDF)
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

  private:

  //! Mesh that covers the source
  std::shared_ptr<const Utility::Mesh> d_mesh;

  //! Direction discretization
  std::shared_ptr<const Utility::PQLAQuadrature> d_direction_discretization;

  //! Energy boundaries
  std::vector<double> d_energy_discretization_boundaries;

};

} // end MonteCarlo namespace

#endif // end MONTE_CARLO_PARTICLE_DISTRIBUTION_HPP

//---------------------------------------------------------------------------//
// end MonteCarlo_ParticleDistribution.hpp
//---------------------------------------------------------------------------//