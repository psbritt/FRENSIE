//---------------------------------------------------------------------------//
//!
//! \file   tstImportanceMesh.cpp
//! \author Philip Britt
//! \brief  ImportanceMesh test
//!
//---------------------------------------------------------------------------//
#include <iostream>
// std includes
#include <memory>

// FRENSIE Includes
#include "MonteCarlo_Importance.hpp"
#include "MonteCarlo_ImportanceMesh.hpp"
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
std::shared_ptr<MonteCarlo::ImportanceMesh> importance_mesh;

//---------------------------------------------------------------------------//
// Tests.
//---------------------------------------------------------------------------//

FRENSIE_UNIT_TEST( ImportanceMesh, getImportance )
{
  MonteCarlo::PhotonState photon(0);

  photon.setEnergy( 1.0 );
  photon.setPosition(0.5, 0.5, 0.5);

  double importance = importance_mesh->getImportance(photon);

  FRENSIE_CHECK_EQUAL(2.0, importance);
}

FRENSIE_UNIT_TEST( ImportanceMesh, checkParticleWithPopulationController_split)
{
  // 25% probability of 1 particle emerging, 75% probability of 2 with threshold at 0.25
  std::vector<double> fake_stream = {0.249, 0.251};
  Utility::RandomNumberGenerator::setFakeStream(fake_stream);
  {
    MonteCarlo::PhotonState photon(0);
    MonteCarlo::ParticleBank particle_bank;
    photon.setEnergy( 1.0);
    photon.setPosition(0.5, 0.5, 0.5);
    photon.setWeight(14.2);

    importance_mesh->checkParticleWithPopulationController(photon, particle_bank);

    FRENSIE_CHECK_EQUAL(photon.getWeight(), 14.2);
    FRENSIE_CHECK_EQUAL(photon.getImportancePair().first, 2.0);
    FRENSIE_CHECK_EQUAL(photon.getCollisionNumber(), 0);

    photon.setEnergy(1e-2);
    photon.setPosition(1.5, 0.5, 0.5);
    photon.incrementCollisionNumber();

    importance_mesh->checkParticleWithPopulationController(photon, particle_bank);

    double updated_weight = 14.2*(2.0/3.5);

    FRENSIE_CHECK_CLOSE(photon.getWeight(), updated_weight, 1e-14);
    FRENSIE_CHECK_EQUAL(photon.getImportancePair().second, 3.5);
    FRENSIE_CHECK(particle_bank.isEmpty());
    FRENSIE_CHECK_EQUAL(photon.getCollisionNumber(), 1);
  }

  {
    MonteCarlo::PhotonState photon(0);
    MonteCarlo::ParticleBank particle_bank;

    photon.setEnergy( 1.0);
    photon.setPosition(0.5, 0.5, 0.5);
    photon.setWeight(14.2);

    importance_mesh->checkParticleWithPopulationController(photon, particle_bank);

    FRENSIE_CHECK_EQUAL(photon.getWeight(), 14.2);
    FRENSIE_CHECK_EQUAL(photon.getImportancePair().first, 2.0);
    FRENSIE_CHECK_EQUAL(photon.getCollisionNumber(), 0);

    photon.setEnergy(1e-2);
    photon.setPosition(1.5, 0.5, 0.5);
    photon.incrementCollisionNumber();

    importance_mesh->checkParticleWithPopulationController(photon, particle_bank);

    double updated_weight = 14.2*(2.0/3.5);

    FRENSIE_CHECK_CLOSE(photon.getWeight(), updated_weight, 1e-14);
    FRENSIE_CHECK_CLOSE(particle_bank.top().getWeight(), updated_weight, 1e-14);
    FRENSIE_CHECK_EQUAL(photon.getImportancePair().second, 3.5);
    FRENSIE_CHECK_EQUAL(particle_bank.size(), 1);
    FRENSIE_CHECK_EQUAL(photon.getCollisionNumber(), 1);
  }
  
  Utility::RandomNumberGenerator::unsetFakeStream(); 
}

FRENSIE_UNIT_TEST(ImportanceMesh, checkParticleWithPopulationController_terminate)
{
  
  // Fake probability stream
  std::vector<double> fake_stream = {0.44, 0.42};
  Utility::RandomNumberGenerator::setFakeStream(fake_stream);
  {
    MonteCarlo::PhotonState photon(0);
    MonteCarlo::ParticleBank particle_bank;
    photon.setEnergy( 1e-2);
    photon.setPosition(1.5, 0.5, 0.5);
    photon.setWeight(3.0);

    importance_mesh->checkParticleWithPopulationController(photon, particle_bank);

    photon.setEnergy( 1.0);
    photon.setPosition(0.5, 0.5, 0.5);
    photon.incrementCollisionNumber();
    FRENSIE_CHECK_EQUAL(photon.getImportancePair().first, 3.5);
 
    importance_mesh->checkParticleWithPopulationController(photon, particle_bank);

    FRENSIE_CHECK(!photon.isGone());
    FRENSIE_CHECK_CLOSE(photon.getWeight(), 3.0*(3.5/2.0), 1e-14);

  }

  {
    MonteCarlo::PhotonState photon(0);
    MonteCarlo::ParticleBank particle_bank;

    photon.setEnergy( 1e-2);
    photon.setPosition(1.5, 0.5, 0.5);
    photon.setWeight(3.0);

    importance_mesh->checkParticleWithPopulationController(photon, particle_bank);

    photon.setEnergy( 1.0);
    photon.setPosition(0.5, 0.5, 0.5);
    photon.incrementCollisionNumber();
 
    importance_mesh->checkParticleWithPopulationController(photon, particle_bank);

    FRENSIE_CHECK(photon.isGone());
  }

  
  Utility::RandomNumberGenerator::unsetFakeStream();
  
}

//---------------------------------------------------------------------------//
// Check that an estimator can be archived
FRENSIE_UNIT_TEST_TEMPLATE_EXPAND( ImportanceMesh,
                                   archive,
                                   TestArchives )
{
  FETCH_TEMPLATE_PARAM( 0, RawOArchive );
  FETCH_TEMPLATE_PARAM( 1, RawIArchive );

  typedef typename std::remove_pointer<RawOArchive>::type OArchive;
  typedef typename std::remove_pointer<RawIArchive>::type IArchive;

  std::string archive_base_name( "test_importance_mesh" );
  std::ostringstream archive_ostream;

  {
    std::unique_ptr<OArchive> oarchive;

    createOArchive( archive_base_name, archive_ostream, oarchive );

    std::shared_ptr<MonteCarlo::ImportanceMesh> mesh_archive_test;
    std::shared_ptr<MonteCarlo::Importance> base_archive_test;

    { 
      mesh_archive_test = std::make_shared<MonteCarlo::ImportanceMesh>();
      base_archive_test = mesh_archive_test;
      // Set up spatial plane vectors
      std::vector<double> x_planes = {0.0, 1.0};
      std::vector<double> y_planes = {0.0, 0.5, 1.0};
      std::vector<double> z_planes = {0.0, 1.0};

      // Set up hex mesh
      std::shared_ptr<Utility::StructuredHexMesh> mesh = std::make_shared<Utility::StructuredHexMesh>(x_planes, y_planes, z_planes);
      mesh_archive_test->setMesh(mesh);

      // Set up discretization
      std::vector<double> energy_bin_boundaries( 3 );
      energy_bin_boundaries[0] = 0.0;
      energy_bin_boundaries[1] = 1.0;
      energy_bin_boundaries[2] = 2.0;

      base_archive_test->setDiscretization<MonteCarlo::OBSERVER_ENERGY_DIMENSION>(energy_bin_boundaries);

      std::unordered_map<Utility::Mesh::ElementHandle, std::vector<double>> importance_mesh_map;

      for(int i = 0; i < 2; ++i)
      {
        std::vector<double> importance_vector;
        for(int j = 0; j < 2; ++j)
        {
          importance_vector.push_back(static_cast<double>(2*i + j) + 1.0);
          std::cout << "Importance value" <<  static_cast<double>(2*i + j) + 1.0 << std::endl;
        }
        importance_mesh_map.emplace(i, importance_vector);
      }

      mesh_archive_test->setImportanceMap(importance_mesh_map);


    }

    FRENSIE_REQUIRE_NO_THROW( (*oarchive) << BOOST_SERIALIZATION_NVP( mesh_archive_test ) );
    FRENSIE_REQUIRE_NO_THROW( (*oarchive) << BOOST_SERIALIZATION_NVP( base_archive_test ) );
  }

  // Copy the archive ostream to an istream
  std::istringstream archive_istream( archive_ostream.str() );

  // Load the archived distributions
  std::unique_ptr<IArchive> iarchive;

  createIArchive( archive_istream, iarchive );

  std::shared_ptr<MonteCarlo::ImportanceMesh> mesh_archive_test;
  std::shared_ptr<MonteCarlo::Importance> base_archive_test;

  FRENSIE_REQUIRE_NO_THROW( (*iarchive) >> BOOST_SERIALIZATION_NVP( mesh_archive_test ) );
  FRENSIE_REQUIRE_NO_THROW( (*iarchive) >> BOOST_SERIALIZATION_NVP( base_archive_test ) );

  iarchive.reset();
  {
    FRENSIE_CHECK( mesh_archive_test.get() == base_archive_test.get() );

    MonteCarlo::PhotonState photon(0);

    photon.setEnergy( 0.5 );
    photon.setPosition(0.25, 0.75, 0.25);

    FRENSIE_CHECK_EQUAL(base_archive_test->getImportance(photon), 3.0);
  }

}

//---------------------------------------------------------------------------//
// Custom Setup
//---------------------------------------------------------------------------//
FRENSIE_CUSTOM_UNIT_TEST_SETUP_BEGIN();

FRENSIE_CUSTOM_UNIT_TEST_INIT()
{

  std::vector<double> x_planes = {0, 1, 2};
  std::vector<double> y_planes = {0, 1};
  std::vector<double> z_planes = {0, 1};

  std::shared_ptr<Utility::StructuredHexMesh> mesh = std::make_shared<Utility::StructuredHexMesh>(x_planes, y_planes, z_planes);
  importance_mesh = std::make_shared<MonteCarlo::ImportanceMesh>();
  importance_mesh->setMesh(mesh);
  
  std::vector<double> energy_bin_boundaries( 3 );
  energy_bin_boundaries[0] = 0.0;
  energy_bin_boundaries[1] = 1e-1;
  energy_bin_boundaries[2] = 20.0;

  importance_mesh->setDiscretization<MonteCarlo::OBSERVER_ENERGY_DIMENSION>(energy_bin_boundaries);

  std::unordered_map<Utility::Mesh::ElementHandle, std::vector<double>> importance_mesh_map;

  std::vector<double> importance_vector_1;

  importance_vector_1.push_back(1.0);
  importance_vector_1.push_back(2.0);

  importance_mesh_map.emplace(0, importance_vector_1);

  std::vector<double> importance_vector_2;

  importance_vector_2.push_back(3.5);
  importance_vector_2.push_back(4.2);

  importance_mesh_map.emplace(1, importance_vector_2);

  importance_mesh->setImportanceMap(importance_mesh_map);

  // Initialize the random number generator
  Utility::RandomNumberGenerator::createStreams();

}

FRENSIE_CUSTOM_UNIT_TEST_SETUP_END();

//---------------------------------------------------------------------------//
// end tstDefaultPopulationController.cpp
//---------------------------------------------------------------------------//
