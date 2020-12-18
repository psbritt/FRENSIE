//---------------------------------------------------------------------------//
//!
//! \file   tstWeightImportanceMesh.cpp
//! \author Philip Britt
//! \brief  ImportanceMesh test
//!
//---------------------------------------------------------------------------//
#include <iostream>
// std includes
#include <memory>

// FRENSIE Includes
#include "MonteCarlo_WeightImportance.hpp"
#include "MonteCarlo_WeightImportanceMesh.hpp"
#include "Utility_UnitTestHarnessWithMain.hpp"
#include "Utility_StructuredHexMesh.hpp"
#include "MonteCarlo_PhotonState.hpp"
#include "MonteCarlo_ParticleBank.hpp"
#include "Utility_RandomNumberGenerator.hpp"
#include "ArchiveTestHelpers.hpp"

//---------------------------------------------------------------------------//
// Testing Types
//---------------------------------------------------------------------------//

typedef TestArchiveHelper::TestArchives TestArchives;

//---------------------------------------------------------------------------//
// Testing Variables
//---------------------------------------------------------------------------//
std::shared_ptr<MonteCarlo::WeightImportanceMesh> weight_importance_mesh;
// Mesh variables
std::vector<double> x_planes;
std::vector<double> y_planes;
std::vector<double> z_planes;
std::shared_ptr<Utility::StructuredHexMesh> mesh;
std::vector<double> energy_bin_boundaries( 3 );
// Weight importance variables
std::unordered_map<Utility::Mesh::ElementHandle, std::vector<double>> weight_importance_mesh_map;

//---------------------------------------------------------------------------//
// Tests.
//---------------------------------------------------------------------------//

FRENSIE_UNIT_TEST( WeightImportanceMesh, getWeightImportance )
{
  MonteCarlo::PhotonState photon(0);

  photon.setEnergy( 1.0 );
  photon.setPosition(0.5, 0.5, 0.5);

  double weight_importance = weight_importance_mesh->getWeightImportance(photon);

  FRENSIE_CHECK_EQUAL(2.0, weight_importance);
}

FRENSIE_UNIT_TEST( WeightImportanceMesh, checkParticleWithPopulationController_split)
{

  MonteCarlo::PhotonState photon( 0 );
  MonteCarlo::ParticleBank particle_bank;
  photon.setEnergy( 1.0 );
  photon.setPosition( 0.5, 0.5, 0.5 );
  photon.setWeight( 4.8 );

  // 60% chance of splitting into 2, 40% chance of splitting into 3. Weight importance here = 2.0
  std::vector<double> fake_stream = { 0.5999, 0.60001 };
  Utility::RandomNumberGenerator::setFakeStream( fake_stream );

  // Split into 2
  weight_importance_mesh->checkParticleWithPopulationController( photon, particle_bank );
  FRENSIE_CHECK_CLOSE( photon.getWeight(), 2.0, 1e-15 );
  FRENSIE_CHECK_CLOSE( particle_bank.top().getWeight(), 2.0, 1e-15 );
  FRENSIE_CHECK_EQUAL( particle_bank.size(), 1 );

  particle_bank.pop();
  photon.setWeight( 4.8 );

  weight_importance_mesh->checkParticleWithPopulationController( photon, particle_bank );

  FRENSIE_CHECK_CLOSE( photon.getWeight(), 2.0, 1e-15 );
  FRENSIE_CHECK_CLOSE( particle_bank.top().getWeight(), 2.0, 1e-15 );
  FRENSIE_CHECK_EQUAL( particle_bank.size(), 2 );
  particle_bank.pop();
  FRENSIE_CHECK_CLOSE( particle_bank.top().getWeight(), 2.0, 1e-15 );
  Utility::RandomNumberGenerator::unsetFakeStream();
}

FRENSIE_UNIT_TEST(WeightImportanceMesh, checkParticleWithPopulationController_terminate)
{
  MonteCarlo::PhotonState photon( 0 );
  MonteCarlo::ParticleBank particle_bank;
  photon.setEnergy( 1.0 );
  photon.setPosition( 0.5, 0.5, 0.5 );
  photon.setWeight( 1.8 );

  // 90% chance of surviving, 10% chance of termination
  std::vector<double> fake_stream = { 0.10001, 0.0999 };
  Utility::RandomNumberGenerator::setFakeStream( fake_stream );

  weight_importance_mesh->checkParticleWithPopulationController( photon, particle_bank );

  FRENSIE_CHECK_CLOSE( photon.getWeight(), 2.0, 1e-15);
  FRENSIE_CHECK_EQUAL( particle_bank.size(), 0);

  photon.setWeight( 1.8 );

  weight_importance_mesh->checkParticleWithPopulationController( photon, particle_bank );

  FRENSIE_CHECK(photon.isGone());
  FRENSIE_CHECK_EQUAL( particle_bank.size(), 0);
  Utility::RandomNumberGenerator::unsetFakeStream();
}

FRENSIE_UNIT_TEST( WeightImportanceMesh, setMaxSplit )
{

  MonteCarlo::PhotonState photon( 0 );
  MonteCarlo::ParticleBank particle_bank;
  photon.setEnergy( 1.0 );
  photon.setPosition( 0.5, 0.5, 0.5 );
  photon.setWeight( 20.0 );

  weight_importance_mesh->checkParticleWithPopulationController( photon, particle_bank );

  FRENSIE_CHECK_EQUAL( photon.getWeight(), 4.0 );
  FRENSIE_CHECK_EQUAL( particle_bank.size(), 4 );

  for(unsigned i = 0; i < 4; i ++)
  {
    particle_bank.pop();
  }

  photon.setWeight( 20.0 );

  weight_importance_mesh->setMaxSplit(6);

  weight_importance_mesh->checkParticleWithPopulationController( photon, particle_bank );

  FRENSIE_CHECK_CLOSE( photon.getWeight(), 20.0/6.0, 1e-15 );
  FRENSIE_CHECK_EQUAL( particle_bank.size(), 5 );

}

//---------------------------------------------------------------------------//
// Check that an estimator can be archived
FRENSIE_UNIT_TEST_TEMPLATE_EXPAND( WeightImportanceMesh,
                                   archive,
                                   TestArchives )
{
  FETCH_TEMPLATE_PARAM( 0, RawOArchive );
  FETCH_TEMPLATE_PARAM( 1, RawIArchive );

  typedef typename std::remove_pointer<RawOArchive>::type OArchive;
  typedef typename std::remove_pointer<RawIArchive>::type IArchive;

  std::string archive_base_name( "test_weight_importance_mesh" );
  std::ostringstream archive_ostream;

  {
    std::unique_ptr<OArchive> oarchive;

    createOArchive( archive_base_name, archive_ostream, oarchive );

    std::shared_ptr<MonteCarlo::WeightImportanceMesh> mesh_archive_test;
    std::shared_ptr<MonteCarlo::WeightImportance> base_archive_test;

    { 
      mesh_archive_test = std::make_shared<MonteCarlo::WeightImportanceMesh>();
      base_archive_test = mesh_archive_test;
      // Set up spatial plane vectors
      std::vector<double> x_planes = {0.0, 1.0};
      std::vector<double> y_planes = {0.0, 0.5, 1.0};
      std::vector<double> z_planes = {0.0, 1.0};

      // Set up hex mesh
      std::shared_ptr<Utility::StructuredHexMesh> mesh = std::make_shared<Utility::StructuredHexMesh>( x_planes, y_planes, z_planes );
      mesh_archive_test->setMesh(mesh);

      // Set up discretization
      std::vector<double> energy_bin_boundaries( 3 );
      energy_bin_boundaries[0] = 0.0;
      energy_bin_boundaries[1] = 1.0;
      energy_bin_boundaries[2] = 2.0;

      base_archive_test->setDiscretization<MonteCarlo::OBSERVER_ENERGY_DIMENSION>(energy_bin_boundaries);

      std::unordered_map<Utility::Mesh::ElementHandle, std::vector<double>> weight_importance_mesh_map;

      for(int spatial_element = 0; spatial_element < 2; ++spatial_element)
      {
        std::vector<double> weight_importance_vector;
        for(int energy_element = 0; energy_element < 2; ++energy_element)
        {
          weight_importance_vector.push_back(static_cast<double>( 2*spatial_element + energy_element ) + 1.0);
        }
        weight_importance_mesh_map.emplace( spatial_element, weight_importance_vector );
      }

      mesh_archive_test->setWeightImportanceMap(weight_importance_mesh_map);


    }

    FRENSIE_REQUIRE_NO_THROW( (*oarchive) << BOOST_SERIALIZATION_NVP( mesh_archive_test ) );
    FRENSIE_REQUIRE_NO_THROW( (*oarchive) << BOOST_SERIALIZATION_NVP( base_archive_test ) );
  }

  // Copy the archive ostream to an istream
  std::istringstream archive_istream( archive_ostream.str() );

  // Load the archived distributions
  std::unique_ptr<IArchive> iarchive;

  createIArchive( archive_istream, iarchive );

  std::shared_ptr<MonteCarlo::WeightImportanceMesh> mesh_archive_test;
  std::shared_ptr<MonteCarlo::WeightImportance> base_archive_test;

  FRENSIE_REQUIRE_NO_THROW( (*iarchive) >> BOOST_SERIALIZATION_NVP( mesh_archive_test ) );
  FRENSIE_REQUIRE_NO_THROW( (*iarchive) >> BOOST_SERIALIZATION_NVP( base_archive_test ) );

  iarchive.reset();
  {
    std::shared_ptr<const Utility::StructuredHexMesh> underlying_mesh = std::dynamic_pointer_cast<const Utility::StructuredHexMesh>(mesh_archive_test->getMesh());

    FRENSIE_CHECK_EQUAL(underlying_mesh->getNumberOfElements(), 2);

    FRENSIE_CHECK_EQUAL(underlying_mesh->getXPlaneLocation(0), 0.0);
    FRENSIE_CHECK_EQUAL(underlying_mesh->getXPlaneLocation(1), 1.0);

    FRENSIE_CHECK_EQUAL(underlying_mesh->getYPlaneLocation(0), 0.0);
    FRENSIE_CHECK_EQUAL(underlying_mesh->getYPlaneLocation(1), 0.5);
    FRENSIE_CHECK_EQUAL(underlying_mesh->getYPlaneLocation(2), 1.0);

    FRENSIE_CHECK_EQUAL(underlying_mesh->getXPlaneLocation(0), 0.0);
    FRENSIE_CHECK_EQUAL(underlying_mesh->getXPlaneLocation(1), 1.0);

    std::vector<MonteCarlo::ObserverPhaseSpaceDimension> dimensions_discretized;
    base_archive_test->getDiscretizedDimensions(dimensions_discretized);


    FRENSIE_CHECK_EQUAL(dimensions_discretized.size(), 1);
    FRENSIE_CHECK_EQUAL(dimensions_discretized[0], MonteCarlo::ObserverPhaseSpaceDimension::OBSERVER_ENERGY_DIMENSION);

    std::vector<double> energy_discretization_bounds;

    base_archive_test->getDiscretization<MonteCarlo::ObserverPhaseSpaceDimension::OBSERVER_ENERGY_DIMENSION>(energy_discretization_bounds);

    FRENSIE_CHECK_EQUAL(energy_discretization_bounds[0], 0.0);
    FRENSIE_CHECK_EQUAL(energy_discretization_bounds[1], 1.0);
    FRENSIE_CHECK_EQUAL(energy_discretization_bounds[2], 2.0);
    FRENSIE_CHECK( mesh_archive_test.get() == base_archive_test.get() );

    const std::unordered_map<Utility::Mesh::ElementHandle, std::vector<double>>& underlying_weight_importance_map = mesh_archive_test->getWeightImportanceMap();
    for(int spatial_element = 0; spatial_element < 2; ++spatial_element)
    {
      for(int energy_element = 0; energy_element < 2; ++energy_element)
      {
        FRENSIE_CHECK_EQUAL(underlying_weight_importance_map.at(spatial_element)[energy_element], static_cast<double>(2*spatial_element + energy_element) + 1.0);
      }
    }
  }


}

//---------------------------------------------------------------------------//
// Custom Setup
//---------------------------------------------------------------------------//
FRENSIE_CUSTOM_UNIT_TEST_SETUP_BEGIN();

FRENSIE_CUSTOM_UNIT_TEST_INIT()
{

  x_planes = {0, 1, 2};
  y_planes = {0, 1};
  z_planes = {0, 1};

  mesh = std::make_shared<Utility::StructuredHexMesh>(x_planes, y_planes, z_planes);

  weight_importance_mesh = std::make_shared<MonteCarlo::WeightImportanceMesh>();
  weight_importance_mesh->setMesh(mesh);
  
  energy_bin_boundaries[0] = 0.0;
  energy_bin_boundaries[1] = 1e-1;
  energy_bin_boundaries[2] = 20.0;

  weight_importance_mesh->setDiscretization<MonteCarlo::OBSERVER_ENERGY_DIMENSION>(energy_bin_boundaries);

  std::vector<double> weight_importance_vector_1;

  weight_importance_vector_1.push_back(1.0);
  weight_importance_vector_1.push_back(2.0);

  weight_importance_mesh_map.emplace(0, weight_importance_vector_1);

  std::vector<double> weight_importance_vector_2;

  weight_importance_vector_2.push_back(3.5);
  weight_importance_vector_2.push_back(4.2);

  weight_importance_mesh_map.emplace(1, weight_importance_vector_2);

  weight_importance_mesh->setWeightImportanceMap(weight_importance_mesh_map);

  // Initialize the random number generator
  Utility::RandomNumberGenerator::createStreams();

}

FRENSIE_CUSTOM_UNIT_TEST_SETUP_END();

//---------------------------------------------------------------------------//
// end tstDefaultPopulationController.cpp
//---------------------------------------------------------------------------//
