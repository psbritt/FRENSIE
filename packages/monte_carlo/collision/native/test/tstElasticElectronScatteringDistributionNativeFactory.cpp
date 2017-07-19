//---------------------------------------------------------------------------//
//!
//! \file   tstElasticElectronScatteringDistributionNativeFactory.cpp
//! \author Luke Kersting
//! \brief  Elastic electron scattering distribution unit tests
//!
//---------------------------------------------------------------------------//

// Std Lib Includes
#include <iostream>

// Trilinos Includes
#include <Teuchos_UnitTestHarness.hpp>
#include <Teuchos_VerboseObject.hpp>
#include <Teuchos_RCP.hpp>

// FRENSIE Includes
#include "MonteCarlo_ElasticElectronScatteringDistributionNativeFactory.hpp"
#include "Data_AdjointElectronPhotonRelaxationDataContainer.hpp"
#include "Data_ElectronPhotonRelaxationDataContainer.hpp"
#include "Utility_RandomNumberGenerator.hpp"
#include "Utility_UnitTestHarnessExtensions.hpp"

//---------------------------------------------------------------------------//
// Testing Structs.
//---------------------------------------------------------------------------//
class TestElasticElectronScatteringDistributionNativeFactory : public MonteCarlo::ElasticElectronScatteringDistributionNativeFactory
{
public:

  // Allow public access to the protected member functions
  using MonteCarlo::ElasticElectronScatteringDistributionNativeFactory::createCrossSectionRatios;
  using MonteCarlo::ElasticElectronScatteringDistributionNativeFactory::createScatteringFunction;
  using MonteCarlo::ElasticElectronScatteringDistributionNativeFactory::createScatteringFunctionInSubrange;
};

//---------------------------------------------------------------------------//
// Testing Variables.
//---------------------------------------------------------------------------//
  typedef Utility::FullyTabularTwoDDistribution TwoDDist;

std::shared_ptr<Data::ElectronPhotonRelaxationDataContainer> data_container;
std::shared_ptr<Data::AdjointElectronPhotonRelaxationDataContainer>
    adjoint_data_container;
std::shared_ptr< const MonteCarlo::CutoffElasticElectronScatteringDistribution>
  native_cutoff_elastic_distribution;
std::shared_ptr< const MonteCarlo::ScreenedRutherfordElasticElectronScatteringDistribution>
  native_sr_elastic_distribution;
std::shared_ptr< const MonteCarlo::MomentPreservingElasticElectronScatteringDistribution>
  native_mp_elastic_distribution;
std::shared_ptr< const MonteCarlo::AnalogElasticElectronScatteringDistribution>
  native_analog_elastic_distribution;

  bool correlated_sampling_mode_on = true;

//---------------------------------------------------------------------------//
// Tests
//---------------------------------------------------------------------------//
// Check that the angular grid can be returned
TEUCHOS_UNIT_TEST( ElasticElectronScatteringDistributionNativeFactory,
                   getAngularGridAboveCutoff )
{

  double cutoff_angle_cosine = -1.0;
  double energy = 1.0e-5;
  std::vector<double> angular_grid =
    MonteCarlo::ElasticElectronScatteringDistributionNativeFactory::getAngularGridAboveCutoff(
                    data_container->getCutoffElasticAngles(),
                    energy,
                    cutoff_angle_cosine );
  // Test
  TEST_EQUALITY_CONST( angular_grid.size(), 2 );
  TEST_EQUALITY_CONST( angular_grid.front(), -1.0 );
  TEST_EQUALITY_CONST( angular_grid.back(), 0.999999 );

  energy = 1.001e-5;
  angular_grid =
    MonteCarlo::ElasticElectronScatteringDistributionNativeFactory::getAngularGridAboveCutoff(
                    data_container->getCutoffElasticAngles(),
                    energy,
                    cutoff_angle_cosine );
  // Test
  TEST_EQUALITY_CONST( angular_grid.size(), 2 );
  TEST_EQUALITY_CONST( angular_grid.front(), -1.0 );
  TEST_EQUALITY_CONST( angular_grid.back(), 0.999999 );

  angular_grid.clear();

  cutoff_angle_cosine = 0.9;
  energy = 1.0e-5;
  angular_grid =
    MonteCarlo::ElasticElectronScatteringDistributionNativeFactory::getAngularGridAboveCutoff(
                    data_container->getCutoffElasticAngles(),
                    energy,
                    cutoff_angle_cosine );
  // Test
  TEST_EQUALITY_CONST( angular_grid.size(), 2 );
  TEST_EQUALITY_CONST( angular_grid.front(), 0.9 );
  TEST_EQUALITY_CONST( angular_grid.back(), 0.999999 );

  angular_grid.clear();

  cutoff_angle_cosine = -1.0;
  energy = 1.0e+5;
  angular_grid =
    MonteCarlo::ElasticElectronScatteringDistributionNativeFactory::getAngularGridAboveCutoff(
                    data_container->getCutoffElasticAngles(),
                    energy,
                    cutoff_angle_cosine );
  // Test
  TEST_EQUALITY_CONST( angular_grid.size(), 90 );
  TEST_EQUALITY_CONST( angular_grid.front(), -1.0 );
  TEST_EQUALITY_CONST( angular_grid.back(), 0.999999 );

  cutoff_angle_cosine = -1.0;
  energy = 9.0e+4;
  angular_grid =
    MonteCarlo::ElasticElectronScatteringDistributionNativeFactory::getAngularGridAboveCutoff(
                    data_container->getCutoffElasticAngles(),
                    energy,
                    cutoff_angle_cosine );
  // Test
  TEST_EQUALITY_CONST( angular_grid.size(), 90 );
  TEST_EQUALITY_CONST( angular_grid.front(), -1.0 );
  TEST_EQUALITY_CONST( angular_grid.back(), 0.999999 );


  cutoff_angle_cosine = -1.0;
  energy = 1.0e-5;
  std::vector<double> raw_grid = 
    data_container->getCutoffElasticAngles( energy );
  angular_grid =
    MonteCarlo::ElasticElectronScatteringDistributionNativeFactory::getAngularGridAboveCutoff(
                                                raw_grid,
                                                cutoff_angle_cosine );
  // Test
  TEST_EQUALITY_CONST( angular_grid.size(), 2 );
  TEST_EQUALITY_CONST( angular_grid.front(), -1.0 );
  TEST_EQUALITY_CONST( angular_grid.back(), 0.999999 );

  angular_grid.clear();

  cutoff_angle_cosine = 0.9;
  angular_grid =
    MonteCarlo::ElasticElectronScatteringDistributionNativeFactory::getAngularGridAboveCutoff(
                                                raw_grid,
                                                cutoff_angle_cosine );
  // Test
  TEST_EQUALITY_CONST( angular_grid.size(), 2 );
  TEST_EQUALITY_CONST( angular_grid.front(), 0.9 );
  TEST_EQUALITY_CONST( angular_grid.back(), 0.999999 );

  angular_grid.clear();

  cutoff_angle_cosine = -1.0;
  energy = 1.0e+5;
  raw_grid = data_container->getCutoffElasticAngles( energy );
  angular_grid =
    MonteCarlo::ElasticElectronScatteringDistributionNativeFactory::getAngularGridAboveCutoff(
                                                raw_grid,
                                                cutoff_angle_cosine );
  // Test
  TEST_EQUALITY_CONST( angular_grid.size(), 90 );
  TEST_EQUALITY_CONST( angular_grid.front(), -1.0 );
  TEST_EQUALITY_CONST( angular_grid.back(), 0.999999 );

}

//---------------------------------------------------------------------------//
// Check that the angular grid can be returned
TEUCHOS_UNIT_TEST( ElasticElectronScatteringDistributionNativeFactory,
                   getAngularGrid )
{
  std::vector<double> angular_grid;
  double cutoff_angle_cosine = 0.9;

  // Test lowerest energy bin
  double energy = 1.0e-5;
  angular_grid =
    MonteCarlo::ElasticElectronScatteringDistributionNativeFactory::getAngularGrid(
        data_container->getCutoffElasticAngles(energy),
        cutoff_angle_cosine );

  TEST_EQUALITY_CONST( angular_grid.size(), 2 );
  TEST_EQUALITY_CONST( angular_grid.front(), -1.0 );
  TEST_EQUALITY_CONST( angular_grid.back(), 0.9 );

  // Test highest energy bin
  energy = 1.0e+5;
  angular_grid =
    MonteCarlo::ElasticElectronScatteringDistributionNativeFactory::getAngularGrid(
        data_container->getCutoffElasticAngles(energy),
        cutoff_angle_cosine );

  TEST_EQUALITY_CONST( angular_grid.size(), 21 );
  TEST_EQUALITY_CONST( angular_grid.front(), -1.0 );
  TEST_EQUALITY_CONST( angular_grid.back(), 0.9 );
}

//---------------------------------------------------------------------------//
// Check that the angular grid can be returned
TEUCHOS_UNIT_TEST( ElasticElectronScatteringDistributionNativeFactory,
                   getAngularGridAndPDF )
{
  std::vector<double> angular_grid, evaluated_pdf;
  double cutoff_angle_cosine = 0.9;

  // Test lowerest energy bin
  double energy = 1.0e-5;
  MonteCarlo::ElasticElectronScatteringDistributionNativeFactory::getAngularGridAndPDF(
    angular_grid,
    evaluated_pdf,
    data_container->getCutoffElasticAngles(energy),
    data_container->getCutoffElasticPDF(energy),
    cutoff_angle_cosine );

  TEST_EQUALITY_CONST( angular_grid.size(), 2 );
  TEST_EQUALITY_CONST( angular_grid.front(), -1.0 );
  TEST_EQUALITY_CONST( angular_grid.back(), 0.9 );
  TEST_EQUALITY_CONST( evaluated_pdf.size(), 2 );
  TEST_EQUALITY_CONST( evaluated_pdf.front(), 0.5 );
  TEST_FLOATING_EQUALITY( evaluated_pdf.back(), 5.49999549999550030e-01, 1e-12 );

  // Test highest energy bin
  energy = 1.0e+5;
  MonteCarlo::ElasticElectronScatteringDistributionNativeFactory::getAngularGridAndPDF(
    angular_grid,
    evaluated_pdf,
    data_container->getCutoffElasticAngles(energy),
    data_container->getCutoffElasticPDF(energy),
    cutoff_angle_cosine );

  TEST_EQUALITY_CONST( angular_grid.size(), 21 );
  TEST_EQUALITY_CONST( angular_grid.front(), -1.0 );
  TEST_EQUALITY_CONST( angular_grid.back(), 0.9 );
  TEST_EQUALITY_CONST( evaluated_pdf.size(), 21 );
  TEST_EQUALITY_CONST( evaluated_pdf.front(), 1.76576e-8 );
  TEST_FLOATING_EQUALITY( evaluated_pdf.back(), 1.70457187499999980e-04, 1e-12 );
}

//---------------------------------------------------------------------------//
// Check that the angular grid can be returned
TEUCHOS_UNIT_TEST( ElasticElectronScatteringDistributionNativeFactory,
                   getAngularGridAndPDF_LinLinLog )
{
  std::vector<double> angular_grid, evaluated_pdf;
  double evaluation_tol = 1e-7;
  double cutoff_angle_cosine = 1.0;

  // Test lowerest energy bin
  double energy = 1.0e-5;
  MonteCarlo::ElasticElectronScatteringDistributionNativeFactory::getAngularGridAndPDF<Utility::LinLinLog>(
    angular_grid,
    evaluated_pdf,
    data_container->getCutoffElasticAngles(),
    data_container->getCutoffElasticPDF(),
    energy,
    cutoff_angle_cosine,
    evaluation_tol );

  TEST_EQUALITY_CONST( angular_grid.size(), 2 );
  TEST_EQUALITY_CONST( angular_grid.front(), -1.0 );
  TEST_EQUALITY_CONST( angular_grid.back(), 0.999999 );
  TEST_EQUALITY_CONST( evaluated_pdf.size(), 2 );
  TEST_EQUALITY_CONST( evaluated_pdf.front(), 0.5 );
  TEST_EQUALITY_CONST( evaluated_pdf.back(), 0.5 );

  // Test inbetween energy bins
  energy = 20.0;
  MonteCarlo::ElasticElectronScatteringDistributionNativeFactory::getAngularGridAndPDF<Utility::LinLinLog>(
    angular_grid,
    evaluated_pdf,
    data_container->getCutoffElasticAngles(),
    data_container->getCutoffElasticPDF(),
    energy,
    cutoff_angle_cosine,
    evaluation_tol );

  TEST_EQUALITY_CONST( angular_grid.size(), 79 );
  TEST_EQUALITY_CONST( angular_grid.front(), -1.0 );
  TEST_EQUALITY_CONST( angular_grid.back(), 0.999999 );
  TEST_EQUALITY_CONST( evaluated_pdf.size(), 79 );
  TEST_EQUALITY_CONST( evaluated_pdf.front(), 5.19221275245657547e-08 );
  TEST_EQUALITY_CONST( evaluated_pdf.back(), 5.06129104868643801e+05 );

  // Test highest energy bin
  energy = 1.0e+5;
  MonteCarlo::ElasticElectronScatteringDistributionNativeFactory::getAngularGridAndPDF<Utility::LinLinLog>(
    angular_grid,
    evaluated_pdf,
    data_container->getCutoffElasticAngles(),
    data_container->getCutoffElasticPDF(),
    energy,
    cutoff_angle_cosine,
    evaluation_tol );

  TEST_EQUALITY_CONST( angular_grid.size(), 90 );
  TEST_EQUALITY_CONST( angular_grid.front(), -1.0 );
  TEST_EQUALITY_CONST( angular_grid.back(), 0.999999 );
  TEST_EQUALITY_CONST( evaluated_pdf.size(), 90 );
  TEST_EQUALITY_CONST( evaluated_pdf.front(), 1.76576e-8 );
  TEST_EQUALITY_CONST( evaluated_pdf.back(), 9.86374e5 );
}

//---------------------------------------------------------------------------//
// Check that the angular grid can be returned
TEUCHOS_UNIT_TEST( ElasticElectronScatteringDistributionNativeFactory,
                   getAngularGridAndPDF_LinLinLin )
{
  std::vector<double> angular_grid, evaluated_pdf;
  double evaluation_tol = 1e-7;
  double cutoff_angle_cosine = 1.0;

  // Test lowerest energy bin
  double energy = 1.0e-5;
  MonteCarlo::ElasticElectronScatteringDistributionNativeFactory::getAngularGridAndPDF<Utility::LinLinLin>(
    angular_grid,
    evaluated_pdf,
    data_container->getCutoffElasticAngles(),
    data_container->getCutoffElasticPDF(),
    energy,
    cutoff_angle_cosine,
    evaluation_tol );

  TEST_EQUALITY_CONST( angular_grid.size(), 2 );
  TEST_EQUALITY_CONST( angular_grid.front(), -1.0 );
  TEST_EQUALITY_CONST( angular_grid.back(), 0.999999 );
  TEST_EQUALITY_CONST( evaluated_pdf.size(), 2 );
  TEST_EQUALITY_CONST( evaluated_pdf.front(), 0.5 );
  TEST_EQUALITY_CONST( evaluated_pdf.back(), 0.5 );

  // Test inbetween energy bins
  energy = 20.0;
  MonteCarlo::ElasticElectronScatteringDistributionNativeFactory::getAngularGridAndPDF<Utility::LinLinLin>(
    angular_grid,
    evaluated_pdf,
    data_container->getCutoffElasticAngles(),
    data_container->getCutoffElasticPDF(),
    energy,
    cutoff_angle_cosine,
    evaluation_tol );

  TEST_EQUALITY_CONST( angular_grid.size(), 79 );
  TEST_EQUALITY_CONST( angular_grid.front(), -1.0 );
  TEST_EQUALITY_CONST( angular_grid.back(), 0.999999 );
  TEST_EQUALITY_CONST( evaluated_pdf.size(), 79 );
  TEST_EQUALITY_CONST( evaluated_pdf.front(), 6.14602802628500258e-08 );
  TEST_EQUALITY_CONST( evaluated_pdf.back(), 4.33429072791244427e+05 );


  // Test highest energy bin
  energy = 1.0e5;
  MonteCarlo::ElasticElectronScatteringDistributionNativeFactory::getAngularGridAndPDF<Utility::LinLinLin>(
    angular_grid,
    evaluated_pdf,
    data_container->getCutoffElasticAngles(),
    data_container->getCutoffElasticPDF(),
    energy,
    cutoff_angle_cosine,
    evaluation_tol );

  TEST_EQUALITY_CONST( angular_grid.size(), 90 );
  TEST_EQUALITY_CONST( angular_grid.front(), -1.0 );
  TEST_EQUALITY_CONST( angular_grid.back(), 0.999999 );
  TEST_EQUALITY_CONST( evaluated_pdf.size(), 90 );
  TEST_EQUALITY_CONST( evaluated_pdf.front(), 1.76576e-8 );
  TEST_EQUALITY_CONST( evaluated_pdf.back(), 9.86374e5 );
}

//---------------------------------------------------------------------------//
// Check that the cutoff distribution can be created
TEUCHOS_UNIT_TEST( ElasticElectronScatteringDistributionNativeFactory,
                   createCutoffElasticDistribution_LinLinLog )
{
  double cutoff_angle_cosine = 1.0;
  double evaluation_tol = 1e-7;

  MonteCarlo::ElasticElectronScatteringDistributionNativeFactory::createCutoffElasticDistribution<Utility::LinLinLog>(
        native_cutoff_elastic_distribution,
        *data_container,
        cutoff_angle_cosine,
        correlated_sampling_mode_on,
        evaluation_tol );

  // Set fake random number stream
  std::vector<double> fake_stream( 2 );
  // Cutoff
  fake_stream[0] = 0.5; // sample angle cosine = 1.249161208881750E-02
  fake_stream[1] = 0.5; // sample angle cosine = 1.249161208881750E-02

  Utility::RandomNumberGenerator::setFakeStream( fake_stream );

  double incoming_energy = 1.0e-3;
  double scattering_angle_cosine, outgoing_energy;

  // sample cutoff
  native_cutoff_elastic_distribution->sample( incoming_energy,
                                              outgoing_energy,
                                              scattering_angle_cosine );

  // Test
  TEST_FLOATING_EQUALITY( scattering_angle_cosine,
                          1.0 - 1.249161208881750E-02,
                          1e-12 );
  TEST_FLOATING_EQUALITY( outgoing_energy, 1.0e-3, 1e-12 );


  incoming_energy = 1.0e-4;

  // sample cutoff
  native_cutoff_elastic_distribution->sample( incoming_energy,
                                              outgoing_energy,
                                              scattering_angle_cosine );

  // Test
  TEST_FLOATING_EQUALITY( scattering_angle_cosine,
                          0.49375394395559125,
                          1e-12 );
  TEST_FLOATING_EQUALITY( outgoing_energy, 1.0e-4, 1e-12 );
}

//---------------------------------------------------------------------------//
// Check that the cutoff distribution can be created
TEUCHOS_UNIT_TEST( ElasticElectronScatteringDistributionNativeFactory,
                   createCutoffElasticDistribution_LinLinLin )
{
  double cutoff_angle_cosine = 1.0;
  double evaluation_tol = 1e-7;

  MonteCarlo::ElasticElectronScatteringDistributionNativeFactory::createCutoffElasticDistribution<Utility::LinLinLin>(
        native_cutoff_elastic_distribution,
        *data_container,
        cutoff_angle_cosine,
        correlated_sampling_mode_on,
        evaluation_tol );

  // Set fake random number stream
  std::vector<double> fake_stream( 2 );
  // Cutoff
  fake_stream[0] = 0.5; // sample angle cosine = 1.249161208881750E-02
  fake_stream[1] = 0.5; // sample angle cosine = 1.249161208881750E-02

  Utility::RandomNumberGenerator::setFakeStream( fake_stream );

  double incoming_energy = 1.0e-3;
  double scattering_angle_cosine, outgoing_energy;

  // sample cutoff
  native_cutoff_elastic_distribution->sample( incoming_energy,
                                              outgoing_energy,
                                              scattering_angle_cosine );

  // Test
  TEST_FLOATING_EQUALITY( scattering_angle_cosine,
                          1.0 - 1.249161208881750E-02,
                          1e-12 );
  TEST_FLOATING_EQUALITY( outgoing_energy, 1.0e-3, 1e-12 );


  incoming_energy = 1.0e-4;

  // sample cutoff
  native_cutoff_elastic_distribution->sample( incoming_energy,
                                              outgoing_energy,
                                              scattering_angle_cosine );

  // Test
  TEST_FLOATING_EQUALITY( scattering_angle_cosine,
                          0.0897730352646529,
                          1e-12 );
  TEST_FLOATING_EQUALITY( outgoing_energy, 1.0e-4, 1e-12 );
}

//---------------------------------------------------------------------------//
// Check that the cutoff distribution can be created
TEUCHOS_UNIT_TEST( ElasticElectronScatteringDistributionNativeFactory,
                   createCutoffElasticDistribution_LinLinLog_adjoint )
{
  double cutoff_angle_cosine = 1.0;
  double evaluation_tol = 1e-7;

  MonteCarlo::ElasticElectronScatteringDistributionNativeFactory::createCutoffElasticDistribution<Utility::LinLinLog>(
        native_cutoff_elastic_distribution,
        *adjoint_data_container,
        cutoff_angle_cosine,
        correlated_sampling_mode_on,
        evaluation_tol );

  // Set fake random number stream
  std::vector<double> fake_stream( 1 );
  fake_stream[0] = 0.5; // sample mu = 0.9875083879111824503
  Utility::RandomNumberGenerator::setFakeStream( fake_stream );

  double incoming_energy = 1.0e-3;
  double scattering_angle_cosine, outgoing_energy;

  // sample from distribution
  native_cutoff_elastic_distribution->sample( incoming_energy,
                                              outgoing_energy,
                                              scattering_angle_cosine );

  // test
  TEST_FLOATING_EQUALITY( scattering_angle_cosine, 0.98549702576858821956, 1e-12 );
  TEST_FLOATING_EQUALITY( outgoing_energy, 1.0e-3, 1e-12 );


  incoming_energy = 1.0e-4;
  // sample from distribution
  native_cutoff_elastic_distribution->sample( incoming_energy,
                                              outgoing_energy,
                                              scattering_angle_cosine );

  // Test
  TEST_FLOATING_EQUALITY( scattering_angle_cosine,
                          0.49274826288429413,
                          1e-12 );
  TEST_FLOATING_EQUALITY( outgoing_energy, 1.0e-4, 1e-12 );
}

//---------------------------------------------------------------------------//
// Check that the cutoff distribution can be created
TEUCHOS_UNIT_TEST( ElasticElectronScatteringDistributionNativeFactory,
                   createCutoffElasticDistribution_LinLinLin_adjoint )
{
  double cutoff_angle_cosine = 1.0;
  double evaluation_tol = 1e-7;

  MonteCarlo::ElasticElectronScatteringDistributionNativeFactory::createCutoffElasticDistribution<Utility::LinLinLin>(
        native_cutoff_elastic_distribution,
        *adjoint_data_container,
        cutoff_angle_cosine,
        correlated_sampling_mode_on,
        evaluation_tol );

  // Set fake random number stream
  std::vector<double> fake_stream( 1 );
  fake_stream[0] = 0.5; // sample mu = 0.9875083879111824503
  Utility::RandomNumberGenerator::setFakeStream( fake_stream );

  double incoming_energy = 1.0e-3;
  double scattering_angle_cosine, outgoing_energy;

  // sample from distribution
  native_cutoff_elastic_distribution->sample( incoming_energy,
                                              outgoing_energy,
                                              scattering_angle_cosine );

  // test
  TEST_FLOATING_EQUALITY( scattering_angle_cosine, 0.98549702576858821956, 1e-12 );
  TEST_FLOATING_EQUALITY( outgoing_energy, 1.0e-3, 1e-12 );


  incoming_energy = 1.0e-4;
  // sample from distribution
  native_cutoff_elastic_distribution->sample( incoming_energy,
                                              outgoing_energy,
                                              scattering_angle_cosine );

  // Test
  TEST_FLOATING_EQUALITY( scattering_angle_cosine,
                          0.0895901841607807,
                          1e-12 );
  TEST_FLOATING_EQUALITY( outgoing_energy, 1.0e-4, 1e-12 );
}

//---------------------------------------------------------------------------//
// Check that the screened Rutherford distribution can be created
TEUCHOS_UNIT_TEST( ElasticElectronScatteringDistributionNativeFactory,
                   createScreenedRutherfordElasticDistribution )
{
  double cutoff_angle_cosine = 1.0;
  double evaluation_tol = 1e-7;

  MonteCarlo::ElasticElectronScatteringDistributionNativeFactory::createCutoffElasticDistribution<Utility::LinLinLog>(
        native_cutoff_elastic_distribution,
        *data_container,
        cutoff_angle_cosine,
        correlated_sampling_mode_on,
        evaluation_tol );

  MonteCarlo::ElasticElectronScatteringDistributionNativeFactory::createScreenedRutherfordElasticDistribution(
        native_sr_elastic_distribution,
        native_cutoff_elastic_distribution,
        data_container->getAtomicNumber() );

  // Set fake random number stream
  std::vector<double> fake_stream( 3 );
  // Screened Rutherford
  fake_stream[0] = 0.0; // sample angle cosine = 1.0
  fake_stream[1] = 0.5; // sample angle cosine = 0.9999995
  fake_stream[2] = 1.0 - 1.0e-15; // sample angle cosine = 0.999999

  Utility::RandomNumberGenerator::setFakeStream( fake_stream );

  double incoming_energy = 1.0e-3;
  double scattering_angle_cosine, outgoing_energy;

  // sample screened rutherford
  native_sr_elastic_distribution->sample( incoming_energy,
                                          outgoing_energy,
                                          scattering_angle_cosine );

  // Test
  TEST_FLOATING_EQUALITY( scattering_angle_cosine, 1.0, 1e-12 );
  TEST_FLOATING_EQUALITY( outgoing_energy, 1.0e-3, 1e-12 );

  // sample screened rutherford
  native_sr_elastic_distribution->sample( incoming_energy,
                                          outgoing_energy,
                                          scattering_angle_cosine );

  // Test
  TEST_FLOATING_EQUALITY( scattering_angle_cosine, 0.9999995, 1e-12 );
  TEST_FLOATING_EQUALITY( outgoing_energy, 1.0e-3, 1e-12 );

  // sample screened rutherford
  native_sr_elastic_distribution->sample( incoming_energy,
                                          outgoing_energy,
                                          scattering_angle_cosine );

  // Test
  TEST_FLOATING_EQUALITY( scattering_angle_cosine, 1.0-1.0e-6, 1e-12 );
  TEST_FLOATING_EQUALITY( outgoing_energy, 1.0e-3, 1e-12 );
}

//---------------------------------------------------------------------------//
// Check that the screened Rutherford distribution can be created
TEUCHOS_UNIT_TEST( ElasticElectronScatteringDistributionNativeFactory,
                   createScreenedRutherfordElasticDistribution_adjoint )
{
  double cutoff_angle_cosine = 1.0;
  double evaluation_tol = 1e-7;

  MonteCarlo::ElasticElectronScatteringDistributionNativeFactory::createCutoffElasticDistribution<Utility::LinLinLog>(
        native_cutoff_elastic_distribution,
        *adjoint_data_container,
        cutoff_angle_cosine,
        correlated_sampling_mode_on,
        evaluation_tol );

  MonteCarlo::ElasticElectronScatteringDistributionNativeFactory::createScreenedRutherfordElasticDistribution(
        native_sr_elastic_distribution,
        native_cutoff_elastic_distribution,
        data_container->getAtomicNumber() );

  // Set fake random number stream
  std::vector<double> fake_stream( 3 );
  // Screened Rutherford
  fake_stream[0] = 0.0; // sample angle cosine = 1.0
  fake_stream[1] = 0.5; // sample angle cosine = 0.9999995
  fake_stream[2] = 1.0 - 1.0e-15; // sample angle cosine = 0.999999

  Utility::RandomNumberGenerator::setFakeStream( fake_stream );

  double incoming_energy = 1.0e-3;
  double scattering_angle_cosine, outgoing_energy;

  // sample screened rutherford
  native_sr_elastic_distribution->sample( incoming_energy,
                                          outgoing_energy,
                                          scattering_angle_cosine );

  // Test
  TEST_FLOATING_EQUALITY( scattering_angle_cosine, 1.0, 1e-12 );
  TEST_FLOATING_EQUALITY( outgoing_energy, 1.0e-3, 1e-12 );

  // sample screened rutherford
  native_sr_elastic_distribution->sample( incoming_energy,
                                          outgoing_energy,
                                          scattering_angle_cosine );

  // Test
  TEST_FLOATING_EQUALITY( scattering_angle_cosine, 0.9999995, 1e-12 );
  TEST_FLOATING_EQUALITY( outgoing_energy, 1.0e-3, 1e-12 );

  // sample screened rutherford
  native_sr_elastic_distribution->sample( incoming_energy,
                                          outgoing_energy,
                                          scattering_angle_cosine );

  // Test
  TEST_FLOATING_EQUALITY( scattering_angle_cosine, 1.0-1.0e-6, 1e-12 );
  TEST_FLOATING_EQUALITY( outgoing_energy, 1.0e-3, 1e-12 );
}

//---------------------------------------------------------------------------//
// Check that the moment preserving distribution can be created
TEUCHOS_UNIT_TEST( ElasticElectronScatteringDistributionNativeFactory,
                   createMomentPreservingElasticDistribution_LinLinLog )
{
  double cutoff_angle_cosine = 0.9;
  double evaluation_tol = 1e-7;

  MonteCarlo::ElasticElectronScatteringDistributionNativeFactory::createMomentPreservingElasticDistribution<Utility::LinLinLog>(
        native_mp_elastic_distribution,
        *data_container,
        cutoff_angle_cosine,
        correlated_sampling_mode_on,
        evaluation_tol );

  // Set fake random number stream
  std::vector<double> fake_stream( 2 );
  fake_stream[0] = 1.44375258484736535e-01; // sample mu = 9.23978505045084053e-01
  fake_stream[1] = 0.02;

  Utility::RandomNumberGenerator::setFakeStream( fake_stream );

  double incoming_energy = 1.0e-3;
  double scattering_angle_cosine, outgoing_energy;

  // sample moment preserving
  native_mp_elastic_distribution->sample( incoming_energy,
                                          outgoing_energy,
                                          scattering_angle_cosine );

  // Test
  TEST_FLOATING_EQUALITY( scattering_angle_cosine, 9.23978505045084053e-01, 1e-12 );
  TEST_FLOATING_EQUALITY( outgoing_energy, 1.0e-3, 1e-12 );


  incoming_energy = 1.0e-4;

  // sample moment preserving
  native_mp_elastic_distribution->sample( incoming_energy,
                                          outgoing_energy,
                                          scattering_angle_cosine );

  // Test
  TEST_FLOATING_EQUALITY( scattering_angle_cosine, 0.9197418038052813, 1e-12 );
  TEST_FLOATING_EQUALITY( outgoing_energy, 1.0e-4, 1e-12 );
}

//---------------------------------------------------------------------------//
// Check that the moment preserving distribution can be created
TEUCHOS_UNIT_TEST( ElasticElectronScatteringDistributionNativeFactory,
                   createMomentPreservingElasticDistribution_LinLinLin )
{
  double cutoff_angle_cosine = 0.9;
  double evaluation_tol = 1e-7;

  MonteCarlo::ElasticElectronScatteringDistributionNativeFactory::createMomentPreservingElasticDistribution<Utility::LinLinLin>(
        native_mp_elastic_distribution,
        *data_container,
        cutoff_angle_cosine,
        correlated_sampling_mode_on,
        evaluation_tol );

  // Set fake random number stream
  std::vector<double> fake_stream( 2 );
  fake_stream[0] = 1.44375258484736535e-01; // sample mu = 9.23978505045084053e-01
  fake_stream[1] = 0.02;

  Utility::RandomNumberGenerator::setFakeStream( fake_stream );

  double incoming_energy = 1.0e-3;
  double scattering_angle_cosine, outgoing_energy;

  // sample moment preserving
  native_mp_elastic_distribution->sample( incoming_energy,
                                          outgoing_energy,
                                          scattering_angle_cosine );

  // Test
  TEST_FLOATING_EQUALITY( scattering_angle_cosine, 9.23978505045084053e-01, 1e-12 );
  TEST_FLOATING_EQUALITY( outgoing_energy, 1.0e-3, 1e-12 );


  incoming_energy = 1.0e-4;

  // sample moment preserving
  native_mp_elastic_distribution->sample( incoming_energy,
                                          outgoing_energy,
                                          scattering_angle_cosine );

  // Test
  TEST_FLOATING_EQUALITY( scattering_angle_cosine, 0.9162754118818063, 1e-12 );
  TEST_FLOATING_EQUALITY( outgoing_energy, 1.0e-4, 1e-12 );
}

//---------------------------------------------------------------------------//
// Check that the moment preserving distribution can be created
TEUCHOS_UNIT_TEST( ElasticElectronScatteringDistributionNativeFactory,
                   createMomentPreservingElasticDistribution_LinLinLog_adjoint )
{
  double cutoff_angle_cosine = 0.9;
  double evaluation_tol = 1e-7;

  MonteCarlo::ElasticElectronScatteringDistributionNativeFactory::createMomentPreservingElasticDistribution<Utility::LinLinLog>(
        native_mp_elastic_distribution,
        *adjoint_data_container,
        cutoff_angle_cosine,
        correlated_sampling_mode_on,
        evaluation_tol );

  // Set fake random number stream
  std::vector<double> fake_stream( 2 );
  fake_stream[0] = 1.44375258484736535e-01; // sample mu = 9.23978505045084053e-01
  fake_stream[1] = 0.02;

  Utility::RandomNumberGenerator::setFakeStream( fake_stream );

  double incoming_energy = 1.0e-3;
  double scattering_angle_cosine, outgoing_energy;

  // sample moment preserving
  native_mp_elastic_distribution->sample( incoming_energy,
                                          outgoing_energy,
                                          scattering_angle_cosine );

  // Test
  TEST_FLOATING_EQUALITY( scattering_angle_cosine, 0.92398608900202416905, 1e-12 );
  TEST_FLOATING_EQUALITY( outgoing_energy, 1.0e-3, 1e-12 );


  incoming_energy = 1.0e-4;

  // sample moment preserving
  native_mp_elastic_distribution->sample( incoming_energy,
                                          outgoing_energy,
                                          scattering_angle_cosine );

  // Test
  TEST_FLOATING_EQUALITY( scattering_angle_cosine, 0.9197455957837513, 1e-12 );
  TEST_FLOATING_EQUALITY( outgoing_energy, 1.0e-4, 1e-12 );
}

//---------------------------------------------------------------------------//
// Check that the moment preserving distribution can be created
TEUCHOS_UNIT_TEST( ElasticElectronScatteringDistributionNativeFactory,
                   createMomentPreservingElasticDistribution_LinLinLin_adjoint )
{
  double cutoff_angle_cosine = 0.9;
  double evaluation_tol = 1e-7;

  MonteCarlo::ElasticElectronScatteringDistributionNativeFactory::createMomentPreservingElasticDistribution<Utility::LinLinLin>(
        native_mp_elastic_distribution,
        *adjoint_data_container,
        cutoff_angle_cosine,
        correlated_sampling_mode_on,
        evaluation_tol );

  // Set fake random number stream
  std::vector<double> fake_stream( 2 );
  fake_stream[0] = 1.44375258484736535e-01; // sample mu = 9.23978505045084053e-01
  fake_stream[1] = 0.02;

  Utility::RandomNumberGenerator::setFakeStream( fake_stream );

  double incoming_energy = 1.0e-3;
  double scattering_angle_cosine, outgoing_energy;

  // sample moment preserving
  native_mp_elastic_distribution->sample( incoming_energy,
                                          outgoing_energy,
                                          scattering_angle_cosine );

  // Test
  TEST_FLOATING_EQUALITY( scattering_angle_cosine, 0.92398608900202416905, 1e-12 );
  TEST_FLOATING_EQUALITY( outgoing_energy, 1.0e-3, 1e-12 );


  incoming_energy = 1.0e-4;

  // sample moment preserving
  native_mp_elastic_distribution->sample( incoming_energy,
                                          outgoing_energy,
                                          scattering_angle_cosine );

  // Test
  TEST_FLOATING_EQUALITY( scattering_angle_cosine, 0.9162761013324372, 1e-12 );
  TEST_FLOATING_EQUALITY( outgoing_energy, 1.0e-4, 1e-12 );
}

//---------------------------------------------------------------------------//
// Check that the analog distribution can be created
TEUCHOS_UNIT_TEST( ElasticElectronScatteringDistributionNativeFactory,
                   createAnalogElasticDistribution_LinLinLog )
{
  double evaluation_tol = 1e-7;

  MonteCarlo::ElasticElectronScatteringDistributionNativeFactory::createAnalogElasticDistribution<Utility::LinLinLog>(
        native_analog_elastic_distribution,
        data_container->getCutoffElasticAngles(),
        data_container->getCutoffElasticPDF(),
        data_container->getElasticAngularEnergyGrid(),
        data_container->getAtomicNumber(),
        correlated_sampling_mode_on,
        evaluation_tol );

  // Set fake random number stream
  std::vector<double> fake_stream( 4 );
  fake_stream[0] = 0.0; // sample angle cosine = -1.0
  fake_stream[1] = 0.5; // sample angle cosine = 0.98751141888536664304
  fake_stream[2] = 0.9; // sample angle cosine = 0.99879780594178391162
  fake_stream[3] = 1.0 - 1.0e-15; // sample angle cosine = 1.0

  Utility::RandomNumberGenerator::setFakeStream( fake_stream );

  double incoming_energy = 1.0e-3;
  double scattering_angle_cosine, outgoing_energy;

  // Test 1
  native_analog_elastic_distribution->sample( incoming_energy,
                                              outgoing_energy,
                                              scattering_angle_cosine );
  TEST_FLOATING_EQUALITY( scattering_angle_cosine, -1.0, 1e-12 );
  TEST_FLOATING_EQUALITY( outgoing_energy, 1.0e-3, 1e-12 );

  // Test 2
  native_analog_elastic_distribution->sample( incoming_energy,
                                              outgoing_energy,
                                              scattering_angle_cosine );
  TEST_FLOATING_EQUALITY( scattering_angle_cosine, 0.98751141888536664304, 1e-12 );
  TEST_FLOATING_EQUALITY( outgoing_energy, 1.0e-3, 1e-12 );

  // Test 3
  native_analog_elastic_distribution->sample( incoming_energy,
                                              outgoing_energy,
                                              scattering_angle_cosine );
  TEST_FLOATING_EQUALITY( scattering_angle_cosine, 0.99879780594178391162, 1e-12 );
  TEST_FLOATING_EQUALITY( outgoing_energy, 1.0e-3, 1e-12 );

  // Test 4
  native_analog_elastic_distribution->sample( incoming_energy,
                                              outgoing_energy,
                                              scattering_angle_cosine );
  TEST_FLOATING_EQUALITY( scattering_angle_cosine, 1.0, 1e-12 );
  TEST_FLOATING_EQUALITY( outgoing_energy, 1.0e-3, 1e-12 );

  // Set fake random number stream
  Utility::RandomNumberGenerator::setFakeStream( fake_stream );

  // Test with an energy inbetween bins
  incoming_energy = 1.0e-4;

  // Test 1
  native_analog_elastic_distribution->sample( incoming_energy,
                                              outgoing_energy,
                                              scattering_angle_cosine );
  TEST_FLOATING_EQUALITY( scattering_angle_cosine, -1.0, 1e-12 );
  TEST_FLOATING_EQUALITY( outgoing_energy, incoming_energy, 1e-12 );

  // Test 2
  native_analog_elastic_distribution->sample( incoming_energy,
                                              outgoing_energy,
                                              scattering_angle_cosine );
  TEST_FLOATING_EQUALITY( scattering_angle_cosine, 0.49375570944268343, 1e-12 );
  TEST_FLOATING_EQUALITY( outgoing_energy, incoming_energy, 1e-12 );

  // Test 3
  native_analog_elastic_distribution->sample( incoming_energy,
                                              outgoing_energy,
                                              scattering_angle_cosine );
  TEST_FLOATING_EQUALITY( scattering_angle_cosine, 0.8993989029708922, 1e-12 );
  TEST_FLOATING_EQUALITY( outgoing_energy, incoming_energy, 1e-12 );

  // Test 4
  native_analog_elastic_distribution->sample( incoming_energy,
                                              outgoing_energy,
                                              scattering_angle_cosine );
  TEST_FLOATING_EQUALITY( scattering_angle_cosine, 1.0, 1e-12 );
  TEST_FLOATING_EQUALITY( outgoing_energy, incoming_energy, 1e-12 );
}

//---------------------------------------------------------------------------//
// Check that the analog distribution can be created
TEUCHOS_UNIT_TEST( ElasticElectronScatteringDistributionNativeFactory,
                   createAnalogElasticDistribution_LinLinLin )
{
  double evaluation_tol = 1e-7;

  MonteCarlo::ElasticElectronScatteringDistributionNativeFactory::createAnalogElasticDistribution<Utility::LinLinLin>(
        native_analog_elastic_distribution,
        data_container->getCutoffElasticAngles(),
        data_container->getCutoffElasticPDF(),
        data_container->getElasticAngularEnergyGrid(),
        data_container->getAtomicNumber(),
        correlated_sampling_mode_on,
        evaluation_tol );

  // Set fake random number stream
  std::vector<double> fake_stream( 4 );
  fake_stream[0] = 0.0; // sample angle cosine = -1.0
  fake_stream[1] = 0.5; // sample angle cosine = 0.98751141888536664304
  fake_stream[2] = 0.9; // sample angle cosine = 0.99879780594178391162
  fake_stream[3] = 1.0 - 1.0e-15; // sample angle cosine = 1.0

  Utility::RandomNumberGenerator::setFakeStream( fake_stream );

  double incoming_energy = 1.0e-3;
  double scattering_angle_cosine, outgoing_energy;

  // Test 1
  native_analog_elastic_distribution->sample( incoming_energy,
                                              outgoing_energy,
                                              scattering_angle_cosine );
  TEST_FLOATING_EQUALITY( scattering_angle_cosine, -1.0, 1e-12 );
  TEST_FLOATING_EQUALITY( outgoing_energy, 1.0e-3, 1e-12 );

  // Test 2
  native_analog_elastic_distribution->sample( incoming_energy,
                                              outgoing_energy,
                                              scattering_angle_cosine );
  TEST_FLOATING_EQUALITY( scattering_angle_cosine, 0.98751141888536664304, 1e-12 );
  TEST_FLOATING_EQUALITY( outgoing_energy, 1.0e-3, 1e-12 );

  // Test 3
  native_analog_elastic_distribution->sample( incoming_energy,
                                              outgoing_energy,
                                              scattering_angle_cosine );
  TEST_FLOATING_EQUALITY( scattering_angle_cosine, 0.99879780594178391162, 1e-12 );
  TEST_FLOATING_EQUALITY( outgoing_energy, 1.0e-3, 1e-12 );

  // Test 4
  native_analog_elastic_distribution->sample( incoming_energy,
                                              outgoing_energy,
                                              scattering_angle_cosine );
  TEST_FLOATING_EQUALITY( scattering_angle_cosine, 1.0, 1e-12 );
  TEST_FLOATING_EQUALITY( outgoing_energy, 1.0e-3, 1e-12 );

  // Set fake random number stream
  Utility::RandomNumberGenerator::setFakeStream( fake_stream );

  // Test with an energy inbetween bins
  incoming_energy = 1.0e-4;

  // Test 1
  native_analog_elastic_distribution->sample( incoming_energy,
                                              outgoing_energy,
                                              scattering_angle_cosine );
  TEST_FLOATING_EQUALITY( scattering_angle_cosine, -1.0, 1e-12 );
  TEST_FLOATING_EQUALITY( outgoing_energy, incoming_energy, 1e-12 );

  // Test 2
  native_analog_elastic_distribution->sample( incoming_energy,
                                              outgoing_energy,
                                              scattering_angle_cosine );
  TEST_FLOATING_EQUALITY( scattering_angle_cosine, 8.9773765352948329e-02, 1e-12 );
  TEST_FLOATING_EQUALITY( outgoing_energy, incoming_energy, 1e-12 );

  // Test 3
  native_analog_elastic_distribution->sample( incoming_energy,
                                              outgoing_energy,
                                              scattering_angle_cosine );
  TEST_FLOATING_EQUALITY( scattering_angle_cosine, 0.8180725278128899, 1e-12 );
  TEST_FLOATING_EQUALITY( outgoing_energy, incoming_energy, 1e-12 );

  // Test 4
  native_analog_elastic_distribution->sample( incoming_energy,
                                              outgoing_energy,
                                              scattering_angle_cosine );
  TEST_FLOATING_EQUALITY( scattering_angle_cosine, 1.0, 1e-12 );
  TEST_FLOATING_EQUALITY( outgoing_energy, incoming_energy, 1e-12 );
}

//---------------------------------------------------------------------------//
// Check that the analog distribution can be created
TEUCHOS_UNIT_TEST( ElasticElectronScatteringDistributionNativeFactory,
                   createAnalogElasticDistribution_LinLinLog_adjoint )
{
  double evaluation_tol = 1e-7;

  MonteCarlo::ElasticElectronScatteringDistributionNativeFactory::createAnalogElasticDistribution<Utility::LinLinLog>(
        native_analog_elastic_distribution,
        adjoint_data_container->getAdjointCutoffElasticAngles(),
        adjoint_data_container->getAdjointCutoffElasticPDF(),
        adjoint_data_container->getAdjointElasticAngularEnergyGrid(),
        adjoint_data_container->getAtomicNumber(),
        correlated_sampling_mode_on,
        evaluation_tol );

  // Set fake random number stream
  std::vector<double> fake_stream( 4 );
  fake_stream[0] = 0.0; // sample angle cosine = -1.0
  fake_stream[1] = 0.5; // sample angle cosine = 0.98751141888536664304
  fake_stream[2] = 0.9; // sample angle cosine = 0.99879780594178391162
  fake_stream[3] = 1.0 - 1.0e-15; // sample angle cosine = 1.0
  Utility::RandomNumberGenerator::setFakeStream( fake_stream );

  double incoming_energy = 1.0e-3;
  double scattering_angle_cosine, outgoing_energy;

  // Test 1
  native_analog_elastic_distribution->sample( incoming_energy,
                                              outgoing_energy,
                                              scattering_angle_cosine );
  TEST_FLOATING_EQUALITY( scattering_angle_cosine, -1.0, 1e-12 );
  TEST_FLOATING_EQUALITY( outgoing_energy, 1.0e-3, 1e-12 );

  // Test 2
  native_analog_elastic_distribution->sample( incoming_energy,
                                              outgoing_energy,
                                              scattering_angle_cosine );
  TEST_FLOATING_EQUALITY( scattering_angle_cosine, 0.98549871388816456808, 1e-12 );
  TEST_FLOATING_EQUALITY( outgoing_energy, 1.0e-3, 1e-12 );

  // Test 3
  native_analog_elastic_distribution->sample( incoming_energy,
                                              outgoing_energy,
                                              scattering_angle_cosine );
  TEST_FLOATING_EQUALITY( scattering_angle_cosine, 0.99826383658664175069, 1e-12 );
  TEST_FLOATING_EQUALITY( outgoing_energy, 1.0e-3, 1e-12 );

  // Test 4
  native_analog_elastic_distribution->sample( incoming_energy,
                                              outgoing_energy,
                                              scattering_angle_cosine );
  TEST_FLOATING_EQUALITY( scattering_angle_cosine, 1.0, 1e-12 );
  TEST_FLOATING_EQUALITY( outgoing_energy, 1.0e-3, 1e-12 );


  // Reset fake random number stream
  Utility::RandomNumberGenerator::setFakeStream( fake_stream );

  // Test with an energy inbetween bins
  incoming_energy = 1.0e-4;

  // Test 1
  native_analog_elastic_distribution->sample( incoming_energy,
                                              outgoing_energy,
                                              scattering_angle_cosine );
  TEST_FLOATING_EQUALITY( scattering_angle_cosine,
                          -1.0,
                          1e-12 );
  TEST_FLOATING_EQUALITY( outgoing_energy, incoming_energy, 1e-12 );

  // Test 2
  native_analog_elastic_distribution->sample( incoming_energy,
                                              outgoing_energy,
                                              scattering_angle_cosine );
  TEST_FLOATING_EQUALITY( scattering_angle_cosine, 0.4927493569443319, 1e-12 );
  TEST_FLOATING_EQUALITY( outgoing_energy, incoming_energy, 1e-12 );

  // Test 3
  native_analog_elastic_distribution->sample( incoming_energy,
                                              outgoing_energy,
                                              scattering_angle_cosine );
  TEST_FLOATING_EQUALITY( scattering_angle_cosine, 0.8991319182937703, 1e-12 );
  TEST_FLOATING_EQUALITY( outgoing_energy, incoming_energy, 1e-12 );

  // Test 4
  native_analog_elastic_distribution->sample( incoming_energy,
                                              outgoing_energy,
                                              scattering_angle_cosine );
  TEST_FLOATING_EQUALITY( scattering_angle_cosine, 1.0, 1e-12 );
  TEST_FLOATING_EQUALITY( outgoing_energy, incoming_energy, 1e-12 );
}

//---------------------------------------------------------------------------//
// Check that the analog distribution can be created
TEUCHOS_UNIT_TEST( ElasticElectronScatteringDistributionNativeFactory,
                   createAnalogElasticDistribution_LinLinLin_adjoint )
{
  double evaluation_tol = 1e-7;

  MonteCarlo::ElasticElectronScatteringDistributionNativeFactory::createAnalogElasticDistribution<Utility::LinLinLin>(
        native_analog_elastic_distribution,
        adjoint_data_container->getAdjointCutoffElasticAngles(),
        adjoint_data_container->getAdjointCutoffElasticPDF(),
        adjoint_data_container->getAdjointElasticAngularEnergyGrid(),
        adjoint_data_container->getAtomicNumber(),
        correlated_sampling_mode_on,
        evaluation_tol );

  // Set fake random number stream
  std::vector<double> fake_stream( 4 );
  fake_stream[0] = 0.0; // sample angle cosine = -1.0
  fake_stream[1] = 0.5; // sample angle cosine = 0.98751141888536664304
  fake_stream[2] = 0.9; // sample angle cosine = 0.99879780594178391162
  fake_stream[3] = 1.0 - 1.0e-15; // sample angle cosine = 1.0
  Utility::RandomNumberGenerator::setFakeStream( fake_stream );

  double incoming_energy = 1.0e-3;
  double scattering_angle_cosine, outgoing_energy;

  // Test 1
  native_analog_elastic_distribution->sample( incoming_energy,
                                              outgoing_energy,
                                              scattering_angle_cosine );
  TEST_FLOATING_EQUALITY( scattering_angle_cosine, -1.0, 1e-12 );
  TEST_FLOATING_EQUALITY( outgoing_energy, 1.0e-3, 1e-12 );

  // Test 2
  native_analog_elastic_distribution->sample( incoming_energy,
                                              outgoing_energy,
                                              scattering_angle_cosine );
  TEST_FLOATING_EQUALITY( scattering_angle_cosine, 0.98549871388816456808, 1e-12 );
  TEST_FLOATING_EQUALITY( outgoing_energy, 1.0e-3, 1e-12 );

  // Test 3
  native_analog_elastic_distribution->sample( incoming_energy,
                                              outgoing_energy,
                                              scattering_angle_cosine );
  TEST_FLOATING_EQUALITY( scattering_angle_cosine, 0.99826383658664175069, 1e-12 );
  TEST_FLOATING_EQUALITY( outgoing_energy, 1.0e-3, 1e-12 );

  // Test 4
  native_analog_elastic_distribution->sample( incoming_energy,
                                              outgoing_energy,
                                              scattering_angle_cosine );
  TEST_FLOATING_EQUALITY( scattering_angle_cosine, 1.0, 1e-12 );
  TEST_FLOATING_EQUALITY( outgoing_energy, 1.0e-3, 1e-12 );


  // Reset fake random number stream
  Utility::RandomNumberGenerator::setFakeStream( fake_stream );

  // Test with an energy inbetween bins
  incoming_energy = 1.0e-4;

  // Test 1
  native_analog_elastic_distribution->sample( incoming_energy,
                                              outgoing_energy,
                                              scattering_angle_cosine );
  TEST_FLOATING_EQUALITY( scattering_angle_cosine, -1.0, 1e-12 );
  TEST_FLOATING_EQUALITY( outgoing_energy, incoming_energy, 1e-12 );

  // Test 2
  native_analog_elastic_distribution->sample( incoming_energy,
                                              outgoing_energy,
                                              scattering_angle_cosine );
  TEST_FLOATING_EQUALITY( scattering_angle_cosine, 8.9590792171902925e-02, 1e-12 );
  TEST_FLOATING_EQUALITY( outgoing_energy, incoming_energy, 1e-12 );

  // Test 3
  native_analog_elastic_distribution->sample( incoming_energy,
                                              outgoing_energy,
                                              scattering_angle_cosine );
  TEST_FLOATING_EQUALITY( scattering_angle_cosine, 0.8180239851450574, 1e-12 );
  TEST_FLOATING_EQUALITY( outgoing_energy, incoming_energy, 1e-12 );

  // Test 4
  native_analog_elastic_distribution->sample( incoming_energy,
                                              outgoing_energy,
                                              scattering_angle_cosine );
  TEST_FLOATING_EQUALITY( scattering_angle_cosine, 1.0, 1e-12 );
  TEST_FLOATING_EQUALITY( outgoing_energy, incoming_energy, 1e-12 );
}

//---------------------------------------------------------------------------//
// Check that the hybrid distribution can be created
TEUCHOS_UNIT_TEST( ElasticElectronScatteringDistributionNativeFactory,
                   createHybridElasticDistribution_LinLinLog )
{
  double cutoff_angle_cosine = 0.9;
  double evaluation_tol = 1e-7;

  Teuchos::ArrayRCP<double> energy_grid;
  energy_grid.assign(
    data_container->getElectronEnergyGrid().begin(),
    data_container->getElectronEnergyGrid().end() );

  Teuchos::RCP<Utility::HashBasedGridSearcher> grid_searcher(   
    new Utility::StandardHashBasedGridSearcher<Teuchos::ArrayRCP<const double>,false>(
         energy_grid,
         energy_grid[0],
         energy_grid[energy_grid.size()-1],
         100 ) );

  Teuchos::ArrayRCP<double> cutoff_cross_section;
  cutoff_cross_section.assign(
    data_container->getCutoffElasticCrossSection().begin(),
    data_container->getCutoffElasticCrossSection().end() );

  Teuchos::ArrayRCP<double> mp_cross_section;
  mp_cross_section.assign(
    data_container->getMomentPreservingCrossSection().begin(),
    data_container->getMomentPreservingCrossSection().end() );

  std::shared_ptr< const MonteCarlo::HybridElasticElectronScatteringDistribution>
    native_hybrid_elastic_distribution;

  MonteCarlo::ElasticElectronScatteringDistributionNativeFactory::createHybridElasticDistribution<Utility::LinLinLog>(
        native_hybrid_elastic_distribution,
        grid_searcher,
        energy_grid,
        cutoff_cross_section,
        mp_cross_section,
        *data_container,
        cutoff_angle_cosine,
        correlated_sampling_mode_on,
        evaluation_tol );

  double incoming_energy = 1.0e-3;

  // Set fake random number stream
  std::vector<double> fake_stream( 6 );
  fake_stream[0] = 1.5069669466805288e-01; // sample mu = 0.1 (cutoff)
  fake_stream[1] = 3.6013179335135637e-01; // sample mu = 0.9 (cutoff)
  fake_stream[2] = 3.6013179335136e-01; // sample mu = 9.2397850504508405e-01 (discrete)
  fake_stream[3] = 4.5251293108241e-01; // sample mu = 9.2397850504508405e-01 (discrete)
  fake_stream[4] = 4.5251293108242e-01; // sample mu = 9.8171108128432372e-01 (discrete)
  fake_stream[5] = 1.0-1e-15; // sample mu = 9.8171108128432372e-01 (discrete)

  Utility::RandomNumberGenerator::setFakeStream( fake_stream );

  double scattering_angle_cosine, outgoing_energy;

  // Test 1
  native_hybrid_elastic_distribution->sample( incoming_energy,
                                              outgoing_energy,
                                              scattering_angle_cosine );
  TEST_FLOATING_EQUALITY( scattering_angle_cosine, 0.1, 1e-12 );
  TEST_FLOATING_EQUALITY( outgoing_energy, incoming_energy, 1e-12 );

  // Test 2
  native_hybrid_elastic_distribution->sample( incoming_energy,
                                              outgoing_energy,
                                              scattering_angle_cosine );
  TEST_FLOATING_EQUALITY( scattering_angle_cosine, 0.9, 1e-12 );
  TEST_FLOATING_EQUALITY( outgoing_energy, incoming_energy, 1e-12 );

  // Test 3
  native_hybrid_elastic_distribution->sample( incoming_energy,
                                              outgoing_energy,
                                              scattering_angle_cosine );
  TEST_FLOATING_EQUALITY( scattering_angle_cosine, 9.2397850504508405e-01, 1e-12 );
  TEST_FLOATING_EQUALITY( outgoing_energy, incoming_energy, 1e-12 );

  // Test 4
  native_hybrid_elastic_distribution->sample( incoming_energy,
                                              outgoing_energy,
                                              scattering_angle_cosine );
  TEST_FLOATING_EQUALITY( scattering_angle_cosine, 9.2397850504508405e-01, 1e-12 );
  TEST_FLOATING_EQUALITY( outgoing_energy, incoming_energy, 1e-12 );

  // Test 5
  native_hybrid_elastic_distribution->sample( incoming_energy,
                                              outgoing_energy,
                                              scattering_angle_cosine );
  TEST_FLOATING_EQUALITY( scattering_angle_cosine, 9.8171108128432372e-01, 1e-12 );
  TEST_FLOATING_EQUALITY( outgoing_energy, incoming_energy, 1e-12 );

  // Test 6
  native_hybrid_elastic_distribution->sample( incoming_energy,
                                              outgoing_energy,
                                              scattering_angle_cosine );
  TEST_FLOATING_EQUALITY( scattering_angle_cosine, 9.8171108128432372e-01, 1e-12 );
  TEST_FLOATING_EQUALITY( outgoing_energy, incoming_energy, 1e-12 );


  // Test with an energy inbetween bins
  incoming_energy = 1.0e-4;

  // Set fake random number stream
  fake_stream.resize( 8 );
  fake_stream[0] = 4.5141090802130296e-01; // sample mu = 0.1 (cutoff)
  fake_stream[1] = 9.2532589918954045e-01; // sample mu = 0.9 (cutoff)
  fake_stream[2] = 9.2532589918955e-01; // sample mu = 9.197418038052812550e-01 (discrete)
  fake_stream[3] = 9.3610699179616e-01; // sample mu = 9.197418038052812550e-01 (discrete)
  fake_stream[4] = 9.3610699179617e-01; // sample mu = 9.486080919249011423e-01 (discrete)
  fake_stream[5] = 9.5694690447057e-01; // sample mu = 9.486080919249011423e-01 (discrete)
  fake_stream[6] = 9.5694690447058e-01; // sample mu = 9.731029893418076115e-01 (discrete)
  fake_stream[7] = 1.0-1e-15; // sample mu = 9.731029893418076115e-01 (discrete)

  Utility::RandomNumberGenerator::setFakeStream( fake_stream );

  // Test 1
  native_hybrid_elastic_distribution->sample( incoming_energy,
                                              outgoing_energy,
                                              scattering_angle_cosine );
  TEST_FLOATING_EQUALITY( scattering_angle_cosine, 0.1, 1e-12 );
  TEST_FLOATING_EQUALITY( outgoing_energy, incoming_energy, 1e-12 );

  // Test 2
  native_hybrid_elastic_distribution->sample( incoming_energy,
                                              outgoing_energy,
                                              scattering_angle_cosine );
  TEST_FLOATING_EQUALITY( scattering_angle_cosine, 0.9, 1e-12 );
  TEST_FLOATING_EQUALITY( outgoing_energy, incoming_energy, 1e-12 );

  // Test 3
  native_hybrid_elastic_distribution->sample( incoming_energy,
                                              outgoing_energy,
                                              scattering_angle_cosine );
  TEST_FLOATING_EQUALITY( scattering_angle_cosine, 9.197418038052812550e-01, 1e-12 );
  TEST_FLOATING_EQUALITY( outgoing_energy, incoming_energy, 1e-12 );

  // Test 4
  native_hybrid_elastic_distribution->sample( incoming_energy,
                                              outgoing_energy,
                                              scattering_angle_cosine );
  TEST_FLOATING_EQUALITY( scattering_angle_cosine, 9.197418038052812550e-01, 1e-12 );
  TEST_FLOATING_EQUALITY( outgoing_energy, incoming_energy, 1e-12 );

  // Test 5
  native_hybrid_elastic_distribution->sample( incoming_energy,
                                              outgoing_energy,
                                              scattering_angle_cosine );
  TEST_FLOATING_EQUALITY( scattering_angle_cosine, 9.486080919249011423e-01, 1e-12 );
  TEST_FLOATING_EQUALITY( outgoing_energy, incoming_energy, 1e-12 );

  // Test 6
  native_hybrid_elastic_distribution->sample( incoming_energy,
                                              outgoing_energy,
                                              scattering_angle_cosine );
  TEST_FLOATING_EQUALITY( scattering_angle_cosine, 9.486080919249011423e-01, 1e-12 );
  TEST_FLOATING_EQUALITY( outgoing_energy, incoming_energy, 1e-12 );

  // Test 7
  native_hybrid_elastic_distribution->sample( incoming_energy,
                                              outgoing_energy,
                                              scattering_angle_cosine );
  TEST_FLOATING_EQUALITY( scattering_angle_cosine, 9.731029893418076115e-01, 1e-12 );
  TEST_FLOATING_EQUALITY( outgoing_energy, incoming_energy, 1e-12 );

  // Test 8
  native_hybrid_elastic_distribution->sample( incoming_energy,
                                              outgoing_energy,
                                              scattering_angle_cosine );
  TEST_FLOATING_EQUALITY( scattering_angle_cosine, 9.731029893418076115e-01, 1e-12 );
  TEST_FLOATING_EQUALITY( outgoing_energy, incoming_energy, 1e-12 );
}

//---------------------------------------------------------------------------//
// Check that the hybrid distribution can be created
TEUCHOS_UNIT_TEST( ElasticElectronScatteringDistributionNativeFactory,
                   createHybridElasticDistribution_LinLinLin )
{
  double cutoff_angle_cosine = 0.9;
  double evaluation_tol = 1e-7;

  Teuchos::ArrayRCP<double> energy_grid;
  energy_grid.assign(
    data_container->getElectronEnergyGrid().begin(),
    data_container->getElectronEnergyGrid().end() );

  Teuchos::RCP<Utility::HashBasedGridSearcher> grid_searcher(   
    new Utility::StandardHashBasedGridSearcher<Teuchos::ArrayRCP<const double>,false>(
         energy_grid,
         energy_grid[0],
         energy_grid[energy_grid.size()-1],
         100 ) );

  Teuchos::ArrayRCP<double> cutoff_cross_section;
  cutoff_cross_section.assign(
    data_container->getCutoffElasticCrossSection().begin(),
    data_container->getCutoffElasticCrossSection().end() );

  Teuchos::ArrayRCP<double> mp_cross_section;
  mp_cross_section.assign(
    data_container->getMomentPreservingCrossSection().begin(),
    data_container->getMomentPreservingCrossSection().end() );

  std::shared_ptr< const MonteCarlo::HybridElasticElectronScatteringDistribution>
    native_hybrid_elastic_distribution;

  MonteCarlo::ElasticElectronScatteringDistributionNativeFactory::createHybridElasticDistribution<Utility::LinLinLin>(
        native_hybrid_elastic_distribution,
        grid_searcher,
        energy_grid,
        cutoff_cross_section,
        mp_cross_section,
        *data_container,
        cutoff_angle_cosine,
        correlated_sampling_mode_on,
        evaluation_tol );

  double incoming_energy = 1.0e-3;

  // Set fake random number stream
  std::vector<double> fake_stream( 6 );
  fake_stream[0] = 1.5069669466805288e-01; // sample mu = 0.1 (cutoff)
  fake_stream[1] = 3.6013179335135637e-01; // sample mu = 0.9 (cutoff)
  fake_stream[2] = 3.6013179335136e-01; // sample mu = 9.2397850504508405e-01 (discrete)
  fake_stream[3] = 4.5251293108241e-01; // sample mu = 9.2397850504508405e-01 (discrete)
  fake_stream[4] = 4.5251293108242e-01; // sample mu = 9.8171108128432372e-01 (discrete)
  fake_stream[5] = 1.0-1e-15; // sample mu = 9.8171108128432372e-01 (discrete)

  Utility::RandomNumberGenerator::setFakeStream( fake_stream );

  double scattering_angle_cosine, outgoing_energy;

  // Test 1
  native_hybrid_elastic_distribution->sample( incoming_energy,
                                              outgoing_energy,
                                              scattering_angle_cosine );
  TEST_FLOATING_EQUALITY( scattering_angle_cosine, 0.1, 1e-12 );
  TEST_FLOATING_EQUALITY( outgoing_energy, incoming_energy, 1e-12 );

  // Test 2
  native_hybrid_elastic_distribution->sample( incoming_energy,
                                              outgoing_energy,
                                              scattering_angle_cosine );
  TEST_FLOATING_EQUALITY( scattering_angle_cosine, 0.9, 1e-12 );
  TEST_FLOATING_EQUALITY( outgoing_energy, incoming_energy, 1e-12 );

  // Test 3
  native_hybrid_elastic_distribution->sample( incoming_energy,
                                              outgoing_energy,
                                              scattering_angle_cosine );
  TEST_FLOATING_EQUALITY( scattering_angle_cosine, 9.2397850504508405e-01, 1e-12 );
  TEST_FLOATING_EQUALITY( outgoing_energy, incoming_energy, 1e-12 );

  // Test 4
  native_hybrid_elastic_distribution->sample( incoming_energy,
                                              outgoing_energy,
                                              scattering_angle_cosine );
  TEST_FLOATING_EQUALITY( scattering_angle_cosine, 9.2397850504508405e-01, 1e-12 );
  TEST_FLOATING_EQUALITY( outgoing_energy, incoming_energy, 1e-12 );

  // Test 5
  native_hybrid_elastic_distribution->sample( incoming_energy,
                                              outgoing_energy,
                                              scattering_angle_cosine );
  TEST_FLOATING_EQUALITY( scattering_angle_cosine, 9.8171108128432372e-01, 1e-12 );
  TEST_FLOATING_EQUALITY( outgoing_energy, incoming_energy, 1e-12 );

  // Test 6
  native_hybrid_elastic_distribution->sample( incoming_energy,
                                              outgoing_energy,
                                              scattering_angle_cosine );
  TEST_FLOATING_EQUALITY( scattering_angle_cosine, 9.8171108128432372e-01, 1e-12 );
  TEST_FLOATING_EQUALITY( outgoing_energy, incoming_energy, 1e-12 );


  // Test with an energy inbetween bins
  incoming_energy = 1.0e-4;

  // Set fake random number stream
  fake_stream.resize( 8 );
  fake_stream[0] = 5.1354078318308127e-01; // sample mu = 0.1 (cutoff)
  fake_stream[1] = 9.2858622885104125e-01; // sample mu = 0.9 (cutoff)
  fake_stream[2] = 9.2858622885105e-01; // sample mu = 9.162754118818062787e-01 (discrete)
  fake_stream[3] = 9.3889661052004e-01; // sample mu = 9.162754118818062787e-01 (discrete)
  fake_stream[4] = 9.3889661052005e-01; // sample mu = 9.215238279035552482e-01 (discrete)
  fake_stream[5] = 9.5882663630330e-01; // sample mu = 9.215238279035552482e-01 (discrete)
  fake_stream[6] = 9.5882663630331e-01; // sample mu = 9.660600050252035054e-01 (discrete)
  fake_stream[7] = 1.0-1e-15; // sample mu = 9.660600050252035054e-01 (discrete)

  Utility::RandomNumberGenerator::setFakeStream( fake_stream );

  // Test 1
  native_hybrid_elastic_distribution->sample( incoming_energy,
                                              outgoing_energy,
                                              scattering_angle_cosine );
  TEST_FLOATING_EQUALITY( scattering_angle_cosine, 0.1, 1e-12 );
  TEST_FLOATING_EQUALITY( outgoing_energy, incoming_energy, 1e-12 );

  // Test 2
  native_hybrid_elastic_distribution->sample( incoming_energy,
                                              outgoing_energy,
                                              scattering_angle_cosine );
  TEST_FLOATING_EQUALITY( scattering_angle_cosine, 0.9, 1e-12 );
  TEST_FLOATING_EQUALITY( outgoing_energy, incoming_energy, 1e-12 );

  // Test 3
  native_hybrid_elastic_distribution->sample( incoming_energy,
                                              outgoing_energy,
                                              scattering_angle_cosine );
  TEST_FLOATING_EQUALITY( scattering_angle_cosine, 9.162754118818062787e-01, 1e-12 );
  TEST_FLOATING_EQUALITY( outgoing_energy, incoming_energy, 1e-12 );

  // Test 4
  native_hybrid_elastic_distribution->sample( incoming_energy,
                                              outgoing_energy,
                                              scattering_angle_cosine );
  TEST_FLOATING_EQUALITY( scattering_angle_cosine, 9.162754118818062787e-01, 1e-12 );
  TEST_FLOATING_EQUALITY( outgoing_energy, incoming_energy, 1e-12 );

  // Test 5
  native_hybrid_elastic_distribution->sample( incoming_energy,
                                              outgoing_energy,
                                              scattering_angle_cosine );
  TEST_FLOATING_EQUALITY( scattering_angle_cosine, 9.215238279035552482e-01, 1e-12 );
  TEST_FLOATING_EQUALITY( outgoing_energy, incoming_energy, 1e-12 );

  // Test 6
  native_hybrid_elastic_distribution->sample( incoming_energy,
                                              outgoing_energy,
                                              scattering_angle_cosine );
  TEST_FLOATING_EQUALITY( scattering_angle_cosine, 9.215238279035552482e-01, 1e-12 );
  TEST_FLOATING_EQUALITY( outgoing_energy, incoming_energy, 1e-12 );

  // Test 7
  native_hybrid_elastic_distribution->sample( incoming_energy,
                                              outgoing_energy,
                                              scattering_angle_cosine );
  TEST_FLOATING_EQUALITY( scattering_angle_cosine, 9.660600050252035054e-01, 1e-12 );
  TEST_FLOATING_EQUALITY( outgoing_energy, incoming_energy, 1e-12 );

  // Test 8
  native_hybrid_elastic_distribution->sample( incoming_energy,
                                              outgoing_energy,
                                              scattering_angle_cosine );
  TEST_FLOATING_EQUALITY( scattering_angle_cosine, 9.660600050252035054e-01, 1e-12 );
  TEST_FLOATING_EQUALITY( outgoing_energy, incoming_energy, 1e-12 );
}

//---------------------------------------------------------------------------//
// Check that the hybrid distribution can be created
TEUCHOS_UNIT_TEST( ElasticElectronScatteringDistributionNativeFactory,
                   createHybridElasticDistribution_LinLinLog_adjoint )
{
  double cutoff_angle_cosine = 0.9;
  double evaluation_tol = 1e-15;

  Teuchos::ArrayRCP<double> energy_grid;
  energy_grid.assign(
    adjoint_data_container->getAdjointElectronEnergyGrid().begin(),
    adjoint_data_container->getAdjointElectronEnergyGrid().end() );

  Teuchos::RCP<Utility::HashBasedGridSearcher> grid_searcher(   
    new Utility::StandardHashBasedGridSearcher<Teuchos::ArrayRCP<const double>,false>(
         energy_grid,
         energy_grid[0],
         energy_grid[energy_grid.size()-1],
         100 ) );

  Teuchos::ArrayRCP<double> cutoff_cross_section;
  cutoff_cross_section.assign(
    adjoint_data_container->getAdjointCutoffElasticCrossSection().begin(),
    adjoint_data_container->getAdjointCutoffElasticCrossSection().end() );

  Teuchos::ArrayRCP<double> mp_cross_section;
  mp_cross_section.assign(
    adjoint_data_container->getAdjointMomentPreservingCrossSection().begin(),
    adjoint_data_container->getAdjointMomentPreservingCrossSection().end() );

  std::shared_ptr< const MonteCarlo::HybridElasticElectronScatteringDistribution>
    native_hybrid_elastic_distribution;

  
  MonteCarlo::ElasticElectronScatteringDistributionNativeFactory::createHybridElasticDistribution<Utility::LinLinLog>(
        native_hybrid_elastic_distribution,
        grid_searcher,
        energy_grid,
        cutoff_cross_section,
        mp_cross_section,
        *adjoint_data_container,
        cutoff_angle_cosine,
        correlated_sampling_mode_on,
        evaluation_tol );

  // Set fake random number stream
  std::vector<double> fake_stream( 6 );
  fake_stream[0] = 1.0239714515463915e-02; // sample mu = 0.1 (cutoff)
  fake_stream[1] = 1.5586513237760702e-01; // sample mu = 0.9 (cutoff)
  fake_stream[2] = 1.55865132377608e-01; // sample mu = 9.23986089002024E-01 (discrete)
  fake_stream[3] = 2.92999320187278e-01; // sample mu = 9.23986089002024E-01 (discrete)
  fake_stream[4] = 2.92999320187279e-01; // sample mu = 9.78892622475528E-01 (discrete)
  fake_stream[5] = 1.0-1e-15; // sample mu = 9.78892622475528E-01 (discrete)


  Utility::RandomNumberGenerator::setFakeStream( fake_stream );

  double incoming_energy = 1.0e-3;
  double scattering_angle_cosine, outgoing_energy;

  // Test 1
  native_hybrid_elastic_distribution->sample( incoming_energy,
                                              outgoing_energy,
                                              scattering_angle_cosine );
  TEST_FLOATING_EQUALITY( scattering_angle_cosine, 0.1, 1e-12 );
  TEST_FLOATING_EQUALITY( outgoing_energy, 1.0e-3, 1e-12 );

  // Test 2
  native_hybrid_elastic_distribution->sample( incoming_energy,
                                              outgoing_energy,
                                              scattering_angle_cosine );
  TEST_FLOATING_EQUALITY( scattering_angle_cosine, 0.9, 1e-12 );
  TEST_FLOATING_EQUALITY( outgoing_energy, 1.0e-3, 1e-12 );

 // Test 3
  native_hybrid_elastic_distribution->sample( incoming_energy,
                                              outgoing_energy,
                                              scattering_angle_cosine );
  TEST_FLOATING_EQUALITY( scattering_angle_cosine, 9.23986089002024E-01, 1e-12 );
  TEST_FLOATING_EQUALITY( outgoing_energy, 1.0e-3, 1e-12 );

 // Test 4
  native_hybrid_elastic_distribution->sample( incoming_energy,
                                              outgoing_energy,
                                              scattering_angle_cosine );
  TEST_FLOATING_EQUALITY( scattering_angle_cosine, 9.23986089002024E-01, 1e-12 );
  TEST_FLOATING_EQUALITY( outgoing_energy, 1.0e-3, 1e-12 );

  // Test 5
  native_hybrid_elastic_distribution->sample( incoming_energy,
                                              outgoing_energy,
                                              scattering_angle_cosine );
  TEST_FLOATING_EQUALITY( scattering_angle_cosine, 9.78892622475528E-01, 1e-12 );
  TEST_FLOATING_EQUALITY( outgoing_energy, 1.0e-3, 1e-12 );

  // Test 6
  native_hybrid_elastic_distribution->sample( incoming_energy,
                                              outgoing_energy,
                                              scattering_angle_cosine );
  TEST_FLOATING_EQUALITY( scattering_angle_cosine, 9.78892622475528E-01, 1e-12 );
  TEST_FLOATING_EQUALITY( outgoing_energy, 1.0e-3, 1e-12 );
}

//---------------------------------------------------------------------------//
// Check that the distribution can be created
TEUCHOS_UNIT_TEST( ElasticElectronScatteringDistributionNativeFactory,
                   createScatteringFunction_LinLinLog )
{
  double evaluation_tol = 1e-7;
  double cutoff_angle_cosine = 1.0;

  std::shared_ptr<Utility::FullyTabularTwoDDistribution> scattering_function;

  TestElasticElectronScatteringDistributionNativeFactory::createScatteringFunction<Utility::LinLinLog>(
    data_container->getCutoffElasticAngles(),
    data_container->getCutoffElasticPDF(),
    data_container->getElasticAngularEnergyGrid(),
    scattering_function,
    cutoff_angle_cosine,
    evaluation_tol );

  TEST_FLOATING_EQUALITY( scattering_function->getLowerBoundOfPrimaryIndepVar(), 1e-5, 1e-12 );
  TEST_FLOATING_EQUALITY( scattering_function->getUpperBoundOfPrimaryIndepVar(), 1e5, 1e-12 );
}

//---------------------------------------------------------------------------//
// Check that the distribution can be created
TEUCHOS_UNIT_TEST( ElasticElectronScatteringDistributionNativeFactory,
                   createScatteringFunction_LinLinLin )
{
  double evaluation_tol = 1e-7;
  double cutoff_angle_cosine = 1.0;

  std::shared_ptr<Utility::FullyTabularTwoDDistribution> scattering_function;

  TestElasticElectronScatteringDistributionNativeFactory::createScatteringFunction<Utility::LinLinLin>(
    data_container->getCutoffElasticAngles(),
    data_container->getCutoffElasticPDF(),
    data_container->getElasticAngularEnergyGrid(),
    scattering_function,
    cutoff_angle_cosine,
    evaluation_tol );

  TEST_FLOATING_EQUALITY( scattering_function->getLowerBoundOfPrimaryIndepVar(), 1e-5, 1e-12 );
  TEST_FLOATING_EQUALITY( scattering_function->getUpperBoundOfPrimaryIndepVar(), 1e5, 1e-12 );
}

//---------------------------------------------------------------------------//
// Check that the distribution can be created
TEUCHOS_UNIT_TEST( ElasticElectronScatteringDistributionNativeFactory,
                   createScatteringFunctionInSubrange )
{
  double cutoff_angle_cosine = 0.9;
  std::vector<double> energy_grid = data_container->getElasticAngularEnergyGrid();

  Utility::Pair<double, std::shared_ptr<const Utility::UnitAwareTabularOneDDistribution<void, void> > > scattering_function;

  TestElasticElectronScatteringDistributionNativeFactory::createScatteringFunctionInSubrange(
    data_container->getCutoffElasticAngles( energy_grid[0] ),
    data_container->getCutoffElasticPDF( energy_grid[0] ),
    energy_grid[0],
    cutoff_angle_cosine,
    scattering_function );

  TEST_FLOATING_EQUALITY( scattering_function.first, 1e-5, 1e-12 );
  TEST_FLOATING_EQUALITY( scattering_function.second->getLowerBoundOfIndepVar(),
                          -1.0,
                          1e-12 );
  TEST_FLOATING_EQUALITY( scattering_function.second->getUpperBoundOfIndepVar(),
                          0.9,
                          1e-12 );
}

//---------------------------------------------------------------------------//
// Check that the hybrid distribution can be created
TEUCHOS_UNIT_TEST( ElasticElectronScatteringDistributionNativeFactory,
                   createCrossSectionRatios )
{
  double cutoff_angle_cosine = 0.9;
  double evaluation_tol = 1e-7;

  Teuchos::ArrayRCP<double> energy_grid;
  energy_grid.assign(
    data_container->getElectronEnergyGrid().begin(),
    data_container->getElectronEnergyGrid().end() );

  Teuchos::ArrayRCP<double> cutoff_cross_section;
  cutoff_cross_section.assign(
    data_container->getCutoffElasticCrossSection().begin(),
    data_container->getCutoffElasticCrossSection().end() );

  Teuchos::ArrayRCP<double> mp_cross_section;
  mp_cross_section.assign(
    data_container->getMomentPreservingCrossSection().begin(),
    data_container->getMomentPreservingCrossSection().end() );

  std::shared_ptr<TwoDDist> cutoff_scattering_function;
  TestElasticElectronScatteringDistributionNativeFactory::createScatteringFunction<Utility::LinLinLin>(
    data_container->getCutoffElasticAngles(),
    data_container->getCutoffElasticPDF(),
    data_container->getElasticAngularEnergyGrid(),
    cutoff_scattering_function,
    1.0,
    evaluation_tol );

  std::shared_ptr<const Utility::OneDDistribution> cross_section_ratios;
  TestElasticElectronScatteringDistributionNativeFactory::createCrossSectionRatios<Utility::LinLinLog>(
    energy_grid,
    cutoff_cross_section,
    mp_cross_section,
    cutoff_scattering_function,
    cutoff_angle_cosine,
    cross_section_ratios );

  TEST_FLOATING_EQUALITY( cross_section_ratios->getLowerBoundOfIndepVar(),
                          1e-5,
                          1e-12 );
  TEST_FLOATING_EQUALITY( cross_section_ratios->getUpperBoundOfIndepVar(),
                          1e5,
                          1e-12 );
}

//---------------------------------------------------------------------------//
// Custom setup
//---------------------------------------------------------------------------//
UTILITY_CUSTOM_TEUCHOS_UNIT_TEST_SETUP_BEGIN();

std::string test_native_file_name, test_adjoint_file_name;

UTILITY_CUSTOM_TEUCHOS_UNIT_TEST_COMMAND_LINE_OPTIONS()
{
  clp().setOption( "test_native_file",
                   &test_native_file_name,
                   "Test Native file name" );
  clp().setOption( "test_adjoint_file",
                   &test_adjoint_file_name,
                   "Test Adjoint file name" );
}

UTILITY_CUSTOM_TEUCHOS_UNIT_TEST_DATA_INITIALIZATION()
{
  // Create the native data file container
  data_container.reset( new Data::ElectronPhotonRelaxationDataContainer(
                        test_native_file_name ) );

  // Create the native adjoint data file container
  adjoint_data_container.reset( new Data::AdjointElectronPhotonRelaxationDataContainer(
                                test_adjoint_file_name ) );

  // Initialize the random number generator
  Utility::RandomNumberGenerator::createStreams();
}

UTILITY_CUSTOM_TEUCHOS_UNIT_TEST_SETUP_END();

//---------------------------------------------------------------------------//
// end tstElasticElectronScatteringDistribution.cpp
//---------------------------------------------------------------------------//