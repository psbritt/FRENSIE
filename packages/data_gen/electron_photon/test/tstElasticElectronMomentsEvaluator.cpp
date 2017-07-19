//---------------------------------------------------------------------------//
//!
//! \file   tstElasticElectronMomentsEvaluator.cpp
//! \author Luke Kersting
//! \brief  Elastic electron cross section moments evaluator unit tests
//!
//---------------------------------------------------------------------------//

// Std Lib Includes
#include <iostream>

// Trilinos Includes
#include <Teuchos_UnitTestHarness.hpp>
#include <Teuchos_VerboseObject.hpp>
#include <Teuchos_RCP.hpp>

// FRENSIE Includes
#include "DataGen_ElasticElectronMomentsEvaluator.hpp"
#include "Data_ElectronPhotonRelaxationDataContainer.hpp"
#include "Utility_UnitTestHarnessExtensions.hpp"
#include "Utility_GaussKronrodIntegrator.hpp"
#include "Utility_StandardHashBasedGridSearcher.hpp"
#include "MonteCarlo_ElasticElectronScatteringDistributionNativeFactory.hpp"
#include "MonteCarlo_ScreenedRutherfordElasticElectroatomicReaction.hpp"
#include "MonteCarlo_CutoffElasticElectroatomicReaction.hpp"

//---------------------------------------------------------------------------//
// Testing Structs.
//---------------------------------------------------------------------------//
class TestElasticElectronMomentsEvaluator : public DataGen::ElasticElectronMomentsEvaluator
{
public:
  TestElasticElectronMomentsEvaluator(
    const Data::ElectronPhotonRelaxationDataContainer& data_container )
    : ElasticElectronMomentsEvaluator( data_container, -1.0, 1e-7, true )
  { /* ... */ }

  ~TestElasticElectronMomentsEvaluator()
  { /* ... */ }

  // Allow public access to the ElasticElectronMomentsEvaluator protected member functions
  using DataGen::ElasticElectronMomentsEvaluator::evaluateCutoffPDFMoment;
  using DataGen::ElasticElectronMomentsEvaluator::evaluateScreenedRutherfordPDFMoment;
  using DataGen::ElasticElectronMomentsEvaluator::evaluateScreenedRutherfordPDFMomentByRecursion;
  using DataGen::ElasticElectronMomentsEvaluator::evaluateScreenedRutherfordPDFMomentByNumericalIntegration;
};

//---------------------------------------------------------------------------//
// Testing Variables
//---------------------------------------------------------------------------//

std::shared_ptr<Data::ElectronPhotonRelaxationDataContainer> pb_data, al_data;
std::shared_ptr<DataGen::ElasticElectronMomentsEvaluator>
  pb_evaluator, pb_lin_evaluator, al_evaluator;

std::shared_ptr<TestElasticElectronMomentsEvaluator>
  test_al_evaluator;

//---------------------------------------------------------------------------//
// Testing Functions.
//---------------------------------------------------------------------------//
bool notEqualZero( double value )
{
  return value != 0.0;
}

//---------------------------------------------------------------------------//
// Tests
//---------------------------------------------------------------------------//
// Check that the legendre expanded screened rutherford pdf can be evaluated
TEUCHOS_UNIT_TEST( ElasticElectronMomentsEvaluator,
                   evaluateLegendreExpandedRutherford )
{
  int n = 4;

  double eta = 2.51317958942017e3;
  double expanded_pdf =
    pb_evaluator->evaluateLegendreExpandedRutherford( 0.999999, 1.0e-5, eta, n );
  double pdf = 4.9999999999999983e-01;
  double legendre_expansion = 0.99999000002250015839;

  UTILITY_TEST_FLOATING_EQUALITY( expanded_pdf, pdf*legendre_expansion, 1e-12 );


  eta = 2.68213671998009;
  expanded_pdf =
    pb_evaluator->evaluateLegendreExpandedRutherford( 1.0, 1.0e-4, eta, n );
  pdf = 45.640157554246108873;
  legendre_expansion = 1.0;


  UTILITY_TEST_FLOATING_EQUALITY( expanded_pdf, pdf*legendre_expansion, 1e-12 );

  expanded_pdf =
    pb_evaluator->evaluateLegendreExpandedRutherford( 1.0, 5.5e-4, eta, n );
  pdf = 79.057620326347489481;
  legendre_expansion = 1.0;

  UTILITY_TEST_FLOATING_EQUALITY( expanded_pdf, pdf*legendre_expansion, 1e-12 );


  expanded_pdf =
    pb_evaluator->evaluateLegendreExpandedRutherford( 1.0, 1.0e-3, eta, n );
  pdf = 90.776239651506159589;
  legendre_expansion = 1.0;

  UTILITY_TEST_FLOATING_EQUALITY( expanded_pdf, pdf*legendre_expansion, 1e-12 );


  eta = 4.14887699806239e-14;
  expanded_pdf =
    pb_evaluator->evaluateLegendreExpandedRutherford( 0.999999, 1.0e5, eta, n );
  pdf = 0.041488766511823843564;
  legendre_expansion = 0.99999000002250015839;

  UTILITY_TEST_FLOATING_EQUALITY( expanded_pdf, pdf*legendre_expansion, 1e-10 );


n = 2;

  eta = 2.51317958942017e3;
  expanded_pdf =
    pb_evaluator->evaluateLegendreExpandedRutherford( 1.0, 1.0e-5, eta, n );
  pdf = 0.50000000039790215567;
  legendre_expansion = 1.0;

  UTILITY_TEST_FLOATING_EQUALITY( expanded_pdf, pdf*legendre_expansion, 1e-12 );
}

//---------------------------------------------------------------------------//
// Check that the legendre expanded screened rutherford pdf can be evaluated
TEUCHOS_UNIT_TEST( ElasticElectronMomentsEvaluator,
                   evaluateLegendreExpandedRutherford_lin )
{
  int n = 4;

  double eta = 2.51317958942017e3;
  double expanded_pdf =
    pb_lin_evaluator->evaluateLegendreExpandedRutherford( 0.999999, 1.0e-5, eta, n );
  double pdf = 0.49999999999999983347;
  double legendre_expansion = 0.99999000002250015839;

  UTILITY_TEST_FLOATING_EQUALITY( expanded_pdf, pdf*legendre_expansion, 1e-12 );


  eta = 2.68213671998009;
  expanded_pdf =
    pb_lin_evaluator->evaluateLegendreExpandedRutherford( 1.0, 1.0e-4, eta, n );
  pdf =  8.7076048088580542128;
  legendre_expansion = 1.0;

  UTILITY_TEST_FLOATING_EQUALITY( expanded_pdf, pdf*legendre_expansion, 1e-12 );

  expanded_pdf =
    pb_lin_evaluator->evaluateLegendreExpandedRutherford( 1.0, 5.5e-4, eta, n );
  pdf = 49.7436061285208595;
  legendre_expansion = 1.0;

  UTILITY_TEST_FLOATING_EQUALITY( expanded_pdf, pdf*legendre_expansion, 1e-12 );

  expanded_pdf =
    pb_lin_evaluator->evaluateLegendreExpandedRutherford( 1.0, 1.0e-3, eta, n );
  pdf = 90.776239651506159589;
  legendre_expansion = 1.0;

  UTILITY_TEST_FLOATING_EQUALITY( expanded_pdf, pdf*legendre_expansion, 1e-12 );


  eta = 4.14887699806239e-14;
  expanded_pdf =
    pb_lin_evaluator->evaluateLegendreExpandedRutherford( 0.999999, 1.0e5, eta, n );
  pdf =  0.041488766511823843564;
  legendre_expansion = 0.99999000002250015839;

  UTILITY_TEST_FLOATING_EQUALITY( expanded_pdf, pdf*legendre_expansion, 1e-10 );


n = 2;

  eta = 2.51317958942017e3;
  expanded_pdf =
    pb_evaluator->evaluateLegendreExpandedRutherford( 1.0, 1.0e-5, eta, n );
  pdf = 0.50000000039790215567;
  legendre_expansion = 1.0;

  UTILITY_TEST_FLOATING_EQUALITY( expanded_pdf, pdf*legendre_expansion, 1e-12 );
}

//---------------------------------------------------------------------------//
// Check that the legendre expanded differential cross section can be evaluated
TEUCHOS_UNIT_TEST( ElasticElectronMomentsEvaluator,
                   evaluateLegendreExpandedPDF )
{
  int n = 4;

  double expanded_pdf = pb_evaluator->evaluateLegendreExpandedPDF( 0.999999, 1.0e-5, n );
  double pdf = 0.50000025000012504606;
  double legendre_expansion = 0.99999000002250015839;
  UTILITY_TEST_FLOATING_EQUALITY( expanded_pdf, pdf*legendre_expansion, 1e-12 );

  expanded_pdf = pb_evaluator->evaluateLegendreExpandedPDF( -1.0, 1.0e-4, n );
  pdf = 0.54842166782042123252;
  legendre_expansion = 1.0;
  UTILITY_TEST_FLOATING_EQUALITY( expanded_pdf, pdf*legendre_expansion, 1e-12 );

  expanded_pdf = pb_evaluator->evaluateLegendreExpandedPDF( 0.999999, 1.0e-4, n );
  pdf = 45.642206638305196975;
  legendre_expansion = 0.99999000002250015839;
  UTILITY_TEST_FLOATING_EQUALITY( expanded_pdf, pdf*legendre_expansion, 1e-12 );

  expanded_pdf = pb_evaluator->evaluateLegendreExpandedPDF( 0.999999, 5.5e-4, n );
  pdf = 79.063811969654963718;
  legendre_expansion = 0.99999000002250015839;
  UTILITY_TEST_FLOATING_EQUALITY( expanded_pdf, pdf*legendre_expansion, 1e-12 );

  expanded_pdf = pb_evaluator->evaluateLegendreExpandedPDF( -1.0, 1.0e-3, n );
  pdf = 0.59684308564071753;
  legendre_expansion = 1.0;
  UTILITY_TEST_FLOATING_EQUALITY( expanded_pdf, pdf*legendre_expansion, 1e-12 );

  expanded_pdf = pb_evaluator->evaluateLegendreExpandedPDF( 0.999999, 1.0e-3, n );
  pdf = 90.784413026610266684;
  legendre_expansion = 0.99999000002250015839;
  UTILITY_TEST_FLOATING_EQUALITY( expanded_pdf, pdf*legendre_expansion, 1e-12 );


  expanded_pdf = pb_evaluator->evaluateLegendreExpandedPDF( 0.9999979, 1.0e5, n );
  pdf = 223727.0620065057301;
  legendre_expansion = 0.99997900009922502385;
  UTILITY_TEST_FLOATING_EQUALITY( expanded_pdf, pdf*legendre_expansion, 1e-10 );

n = 2;

  expanded_pdf = pb_evaluator->evaluateLegendreExpandedPDF( -1.0, 1.0e-5, n );
  pdf = 0.50000025000012504606;
  legendre_expansion = 1.0;
  UTILITY_TEST_FLOATING_EQUALITY( expanded_pdf, pdf*legendre_expansion, 1e-12 );

  expanded_pdf = pb_evaluator->evaluateLegendreExpandedPDF( 0.999999, 1.0e-5, n );
  pdf =  0.50000025000012504606;
  legendre_expansion = 0.99999700000149993606;
  UTILITY_TEST_FLOATING_EQUALITY( expanded_pdf, pdf*legendre_expansion, 1e-12 );
}


//---------------------------------------------------------------------------//
// Check that the legendre expanded differential cross section can be evaluated
TEUCHOS_UNIT_TEST( ElasticElectronMomentsEvaluator,
                   evaluateLegendreExpandedPDF_lin )
{
  int n = 4;

  double expanded_pdf =
    pb_lin_evaluator->evaluateLegendreExpandedPDF( 0.999999, 1.0e-5, n );
  double pdf = 0.50000025000012504606;
  double legendre_expansion = 0.99999000002250015839;

  UTILITY_TEST_FLOATING_EQUALITY( expanded_pdf, pdf*legendre_expansion, 1e-12 );

  expanded_pdf =
    pb_lin_evaluator->evaluateLegendreExpandedPDF( -1.0, 1.0e-4, n );
  pdf = 0.50880414414926977695;
  legendre_expansion = 1.0;

  UTILITY_TEST_FLOATING_EQUALITY( expanded_pdf, pdf*legendre_expansion, 1e-12 );

  expanded_pdf =
    pb_lin_evaluator->evaluateLegendreExpandedPDF( 0.999999, 1.0e-4, n );
  pdf = 8.707674138782865114;
  legendre_expansion = 0.99999000002250015839;

  UTILITY_TEST_FLOATING_EQUALITY( expanded_pdf, pdf*legendre_expansion, 1e-12 );

  expanded_pdf =
    pb_lin_evaluator->evaluateLegendreExpandedPDF( 0.999999, 5.5e-4, n );
  pdf = 49.746043582696565011;
  legendre_expansion = 0.99999000002250015839;

  UTILITY_TEST_FLOATING_EQUALITY( expanded_pdf, pdf*legendre_expansion, 1e-12 );

  expanded_pdf =
    pb_lin_evaluator->evaluateLegendreExpandedPDF( -1.0, 1.0e-3, n );

  pdf = 0.59684308564071753;
  legendre_expansion = 1.0;

  UTILITY_TEST_FLOATING_EQUALITY( expanded_pdf, pdf*legendre_expansion, 1e-12 );

  expanded_pdf =
    pb_lin_evaluator->evaluateLegendreExpandedPDF( 0.999999, 1.0e-3, n );

  pdf = 90.784413026610266684;
  legendre_expansion = 0.99999000002250015839;

  UTILITY_TEST_FLOATING_EQUALITY( expanded_pdf, pdf*legendre_expansion, 1e-12 );


  expanded_pdf =
    pb_lin_evaluator->evaluateLegendreExpandedPDF( 0.9999979, 1.0e5, n );

  pdf =  223727.0620065057301;
  legendre_expansion = 0.99997900009922502385;

  UTILITY_TEST_FLOATING_EQUALITY( expanded_pdf, pdf*legendre_expansion, 1e-10 );

n = 2;

  expanded_pdf =
    pb_lin_evaluator->evaluateLegendreExpandedPDF( -1.0, 1.0e-5, n );

  pdf = 0.50000025000012504606;
  legendre_expansion = 1.0;

  UTILITY_TEST_FLOATING_EQUALITY( expanded_pdf, pdf*legendre_expansion, 1e-12 );

  expanded_pdf =
    pb_lin_evaluator->evaluateLegendreExpandedPDF( 0.999999, 1.0e-5, n );

  pdf = 0.50000025000012504606;
  legendre_expansion = 0.99999700000149993606;

  UTILITY_TEST_FLOATING_EQUALITY( expanded_pdf, pdf*legendre_expansion, 1e-12 );
}

//---------------------------------------------------------------------------//
// Check that the screened rutherford moments can be evaluated
TEUCHOS_UNIT_TEST( ElasticElectronMomentsEvaluator,
                   evaluateScreenedRutherfordPDFMomentByRecursion )
{
  unsigned n = 2;
  double energy = 1.0e5;
  double tol = 1e-14;
  std::vector<Utility::long_float> moments(n+1);
  Utility::long_float rutherford_moment;
  Utility::long_float eta;

  eta = (Utility::long_float)2.855935652488E-15L;
  eta *= 2;

  moments.clear();
  moments[0] = 1.0L;
  moments[1] = 0.99999999999989729645L;
  moments[2] = 0.99999999999969188935L;

  test_al_evaluator->
    evaluateScreenedRutherfordPDFMomentByRecursion(
        rutherford_moment,
        eta,
        0 );

  UTILITY_TEST_FLOATING_EQUALITY( rutherford_moment.convert_to<double>(),
                                  moments[0].convert_to<double>(),
                                  tol );

  test_al_evaluator->
    evaluateScreenedRutherfordPDFMomentByRecursion(
        rutherford_moment,
        eta,
        1 );

  UTILITY_TEST_FLOATING_EQUALITY( rutherford_moment.convert_to<double>(),
                                  moments[1].convert_to<double>(),
                                  tol );

  test_al_evaluator->
    evaluateScreenedRutherfordPDFMomentByRecursion(
        rutherford_moment,
        eta,
        2 );

  UTILITY_TEST_FLOATING_EQUALITY( rutherford_moment.convert_to<double>(),
                                  moments[2].convert_to<double>(),
                                  tol );


  eta = (Utility::long_float)6.407445779740e-9L;
  eta *= 2;
  energy = 6.625e1L;

  moments.clear();
  moments[0] = 1.0L;
  moments[1] = 0.99999995609771409924L;
  moments[2] = 0.99999986829315983224L;

  test_al_evaluator->
    evaluateScreenedRutherfordPDFMomentByRecursion(
        rutherford_moment,
        eta,
        0 );

  UTILITY_TEST_FLOATING_EQUALITY( rutherford_moment.convert_to<double>(),
                                  moments[0].convert_to<double>(),
                                  tol );

  test_al_evaluator->
    evaluateScreenedRutherfordPDFMomentByRecursion(
        rutherford_moment,
        eta,
        1 );

  UTILITY_TEST_FLOATING_EQUALITY( rutherford_moment.convert_to<double>(),
                                  moments[1].convert_to<double>(),
                                  tol );

  test_al_evaluator->
    evaluateScreenedRutherfordPDFMomentByRecursion(
        rutherford_moment,
        eta,
        2 );

  UTILITY_TEST_FLOATING_EQUALITY( rutherford_moment.convert_to<double>(),
                                  moments[2].convert_to<double>(),
                                  tol );

  eta = (Utility::long_float)1.457796275283e-8L;
  eta *= 2;
  energy = 4.3750e1L;

  moments.clear();
  moments[0] = 1.0L;
  moments[1] = 0.99999992221948305306L;
  moments[2] = 0.99999976665848608979L;

  test_al_evaluator->
    evaluateScreenedRutherfordPDFMomentByRecursion(
        rutherford_moment,
        eta,
        0 );

  UTILITY_TEST_FLOATING_EQUALITY( rutherford_moment.convert_to<double>(),
                                  moments[0].convert_to<double>(),
                                  tol );

  test_al_evaluator->
    evaluateScreenedRutherfordPDFMomentByRecursion(
        rutherford_moment,
        eta,
        1 );

  UTILITY_TEST_FLOATING_EQUALITY( rutherford_moment.convert_to<double>(),
                                  moments[1].convert_to<double>(),
                                  tol );

  test_al_evaluator->
    evaluateScreenedRutherfordPDFMomentByRecursion(
        rutherford_moment,
        eta,
        2 );

  UTILITY_TEST_FLOATING_EQUALITY( rutherford_moment.convert_to<double>(),
                                  moments[2].convert_to<double>(),
                                  tol );

  energy = 1.0e-5;
  eta = (Utility::long_float)1.189646992674e1L;
  eta *= 2;

  moments.clear();
  moments[0] = 1.0L;
  moments[1] = 0.99999950000000700488L;
  moments[2] = 0.99999850000052101463L;

  test_al_evaluator->
    evaluateScreenedRutherfordPDFMomentByRecursion(
        rutherford_moment,
        eta,
        0 );

  UTILITY_TEST_FLOATING_EQUALITY( rutherford_moment.convert_to<double>(),
                                  moments[0].convert_to<double>(),
                                  tol );

  test_al_evaluator->
    evaluateScreenedRutherfordPDFMomentByRecursion(
        rutherford_moment,
        eta,
        1 );

  UTILITY_TEST_FLOATING_EQUALITY( rutherford_moment.convert_to<double>(),
                                  moments[1].convert_to<double>(),
                                  tol );

  test_al_evaluator->
    evaluateScreenedRutherfordPDFMomentByRecursion(
        rutherford_moment,
        eta,
        2 );

  UTILITY_TEST_FLOATING_EQUALITY( rutherford_moment.convert_to<double>(),
                                  moments[2].convert_to<double>(),
                                  tol );
}

//---------------------------------------------------------------------------//
// Check that the screened rutherford moments can be evaluated
TEUCHOS_UNIT_TEST( ElasticElectronMomentsEvaluator,
                   evaluateScreenedRutherfordPDFMomentByNumericalIntegration )
{
  double energy = 1.0e5;
  double tol = 1e-11;
  std::vector<Utility::long_float> moments(3);
  Utility::long_float rutherford_moment;
  Utility::long_float eta;

  moments.clear();
  moments[0] = 1.0L;
  moments[1] = 0.99999999999989729645L;
  moments[2] = 0.99999999999969188935L;

  test_al_evaluator->
    evaluateScreenedRutherfordPDFMomentByNumericalIntegration(
        rutherford_moment,
        energy,
        0,
        tol,
        100 );

  UTILITY_TEST_FLOATING_EQUALITY( rutherford_moment.convert_to<double>(),
                                  moments[0].convert_to<double>(),
                                  tol );

  test_al_evaluator->
    evaluateScreenedRutherfordPDFMomentByNumericalIntegration(
        rutherford_moment,
        energy,
        1,
        tol,
        100 );

  UTILITY_TEST_FLOATING_EQUALITY( rutherford_moment.convert_to<double>(),
                                  moments[1].convert_to<double>(),
                                  tol );

  test_al_evaluator->
    evaluateScreenedRutherfordPDFMomentByNumericalIntegration(
        rutherford_moment,
        energy,
        2,
        tol,
        100 );

  UTILITY_TEST_FLOATING_EQUALITY( rutherford_moment.convert_to<double>(),
                                  moments[2].convert_to<double>(),
                                  tol );

/*  NOTE: These tests can be added but the numerical integration is time
    consuming and the other tests cover the spread of energies.
 */

//  energy = 6.625e1L;

//  moments.clear();
//  moments[0] = 1.0L;
//  moments[1] = 0.99999995609771409924L;
//  moments[2] = 0.99999986829315983224L;

//  test_al_evaluator->
//    evaluateScreenedRutherfordPDFMomentByNumericalIntegration(
//        rutherford_moment,
//        energy,
//        0,
//        tol,
//        100 );

//  UTILITY_TEST_FLOATING_EQUALITY( rutherford_moment.convert_to<double>(),
//                                  moments[0].convert_to<double>(),
//                                  tol );

//  test_al_evaluator->
//    evaluateScreenedRutherfordPDFMomentByNumericalIntegration(
//        rutherford_moment,
//        energy,
//        1,
//        tol,
//        100 );

//  UTILITY_TEST_FLOATING_EQUALITY( rutherford_moment.convert_to<double>(),
//                                  moments[1].convert_to<double>(),
//                                  tol );

//  test_al_evaluator->
//    evaluateScreenedRutherfordPDFMomentByNumericalIntegration(
//        rutherford_moment,
//        energy,
//        2,
//        tol,
//        100 );

//  UTILITY_TEST_FLOATING_EQUALITY( rutherford_moment.convert_to<double>(),
//                                  moments[2].convert_to<double>(),
//                                  tol );

//  energy = 4.3750e1L;

//  moments.clear();
//  moments[0] = 1.0L;
//  moments[1] = 0.99999992221948305306L;
//  moments[2] = 0.99999976665848608979L;

//  test_al_evaluator->
//    evaluateScreenedRutherfordPDFMomentByNumericalIntegration(
//        rutherford_moment,
//        energy,
//        0,
//        tol,
//        100 );

//  UTILITY_TEST_FLOATING_EQUALITY( rutherford_moment.convert_to<double>(),
//                                  moments[0].convert_to<double>(),
//                                  tol );

//  test_al_evaluator->
//    evaluateScreenedRutherfordPDFMomentByNumericalIntegration(
//        rutherford_moment,
//        energy,
//        1,
//        tol,
//        100 );

//  UTILITY_TEST_FLOATING_EQUALITY( rutherford_moment.convert_to<double>(),
//                                  moments[1].convert_to<double>(),
//                                  tol );

//  test_al_evaluator->
//    evaluateScreenedRutherfordPDFMomentByNumericalIntegration(
//        rutherford_moment,
//        energy,
//        2,
//        tol,
//        100 );

//  UTILITY_TEST_FLOATING_EQUALITY( rutherford_moment.convert_to<double>(),
//                                  moments[2].convert_to<double>(),
//                                  tol );

  energy = 1.0e-5;

  moments.clear();
  moments[0] = 1.0L;
  moments[1] = 0.99999950000000700488L;
  moments[2] = 0.99999850000052101463L;

  test_al_evaluator->
    evaluateScreenedRutherfordPDFMomentByNumericalIntegration(
        rutherford_moment,
        energy,
        0,
        tol,
        100 );

  UTILITY_TEST_FLOATING_EQUALITY( rutherford_moment.convert_to<double>(),
                                  moments[0].convert_to<double>(),
                                  tol );

  test_al_evaluator->
    evaluateScreenedRutherfordPDFMomentByNumericalIntegration(
        rutherford_moment,
        energy,
        1,
        tol,
        100 );

  UTILITY_TEST_FLOATING_EQUALITY( rutherford_moment.convert_to<double>(),
                                  moments[1].convert_to<double>(),
                                  tol );

  test_al_evaluator->
    evaluateScreenedRutherfordPDFMomentByNumericalIntegration(
        rutherford_moment,
        energy,
        2,
        tol,
        100 );

  UTILITY_TEST_FLOATING_EQUALITY( rutherford_moment.convert_to<double>(),
                                  moments[2].convert_to<double>(),
                                  tol );
}

//---------------------------------------------------------------------------//
// Check that the screened rutherford moments can be evaluated
TEUCHOS_UNIT_TEST( ElasticElectronMomentsEvaluator,
                   evaluateScreenedRutherfordPDFMoment )
{
  unsigned n = 2;
  double energy = 1.0e5;
  double tol = 1e-14;
  std::vector<Utility::long_float> moments(n+1);
  Utility::long_float rutherford_moment;

  moments[0] = 1.0;
  moments[1] = 9.999999999999e-1L;
  moments[2] = 9.999999999997e-1L;

   test_al_evaluator->
    evaluateScreenedRutherfordPDFMoment( rutherford_moment,
                                         energy,
                                         0 );

  UTILITY_TEST_FLOATING_EQUALITY( rutherford_moment.convert_to<double>(),
                                  moments[0].convert_to<double>(),
                                  tol );

   test_al_evaluator->
    evaluateScreenedRutherfordPDFMoment( rutherford_moment,
                                         energy,
                                         1 );

  UTILITY_TEST_FLOATING_EQUALITY( rutherford_moment.convert_to<double>(),
                                  moments[1].convert_to<double>(),
                                  tol );

   test_al_evaluator->
    evaluateScreenedRutherfordPDFMoment( rutherford_moment,
                                         energy,
                                         2 );

  UTILITY_TEST_FLOATING_EQUALITY( rutherford_moment.convert_to<double>(),
                                  moments[2].convert_to<double>(),
                                  tol );

  energy = 4.3750e1;

  tol = 1e-16;
  moments.clear();
  moments[0] = 1.0L;
  moments[1] = 0.99999992222034858054L;
  moments[2] = 0.99999976666108267169L;

   test_al_evaluator->
    evaluateScreenedRutherfordPDFMoment( rutherford_moment,
                                         energy,
                                         0 );

  UTILITY_TEST_FLOATING_EQUALITY( rutherford_moment.convert_to<double>(),
                                  moments[0].convert_to<double>(),
                                  tol );

   test_al_evaluator->
    evaluateScreenedRutherfordPDFMoment( rutherford_moment,
                                         energy,
                                         1 );

  UTILITY_TEST_FLOATING_EQUALITY( rutherford_moment.convert_to<double>(),
                                  moments[1].convert_to<double>(),
                                  tol );

   test_al_evaluator->
    evaluateScreenedRutherfordPDFMoment( rutherford_moment,
                                         energy,
                                         2 );

  UTILITY_TEST_FLOATING_EQUALITY( rutherford_moment.convert_to<double>(),
                                  moments[2].convert_to<double>(),
                                  tol );

  energy = 1.0e-5;

  moments.clear();
  moments[0] = 1.0L;
  moments[1] = 0.99999950000000700488L;
  moments[2] = 0.99999850000052101463L;

  test_al_evaluator->
    evaluateScreenedRutherfordPDFMoment( rutherford_moment,
                                         energy,
                                         0 );

  UTILITY_TEST_FLOATING_EQUALITY( rutherford_moment.convert_to<double>(),
                                  moments[0].convert_to<double>(),
                                  tol );

  test_al_evaluator->
    evaluateScreenedRutherfordPDFMoment( rutherford_moment,
                                         energy,
                                         1 );

  UTILITY_TEST_FLOATING_EQUALITY( rutherford_moment.convert_to<double>(),
                                  moments[1].convert_to<double>(),
                                  tol );

  test_al_evaluator->
    evaluateScreenedRutherfordPDFMoment( rutherford_moment,
                                         energy,
                                         2 );

  UTILITY_TEST_FLOATING_EQUALITY( rutherford_moment.convert_to<double>(),
                                  moments[2].convert_to<double>(),
                                  tol );
}

//---------------------------------------------------------------------------//
// Check that the screened rutherford moments can be evaluated
TEUCHOS_UNIT_TEST( ElasticElectronMomentsEvaluator,
                   evaluateScreenedRutherfordMoment )
{
  unsigned n = 6;
  double energy = 1.0e5;
  double tol = 1e-15;
  std::vector<Utility::long_float> moments(n+1);
  Utility::long_float rutherford_moment;

  long double cross_section = 2.11160991169490e6L;

  moments[0] = cross_section;
  moments[1] = 9.9999999999933626914e-1L*cross_section;
  moments[2] = 9.9999999999800880748e-1L*cross_section;
  moments[3] = 9.9999999999601761514e-1L*cross_section;
  moments[4] = 9.9999999999336269232e-1L*cross_section;
  moments[5] = 9.9999999999004403926e-1L*cross_section;
  moments[6] = 9.9999999998606165627e-1L*cross_section;

  pb_evaluator->
    evaluateScreenedRutherfordMoment( rutherford_moment,
                                      energy,
                                      0 );

  UTILITY_TEST_FLOATING_EQUALITY( rutherford_moment.convert_to<double>(),
                                  moments[0].convert_to<double>(),
                                  tol );

  pb_evaluator->
    evaluateScreenedRutherfordMoment( rutherford_moment,
                                      energy,
                                      1 );

  UTILITY_TEST_FLOATING_EQUALITY( rutherford_moment.convert_to<double>(),
                                  moments[1].convert_to<double>(),
                                  tol );

  pb_evaluator->
    evaluateScreenedRutherfordMoment( rutherford_moment,
                                      energy,
                                      2 );

  UTILITY_TEST_FLOATING_EQUALITY( rutherford_moment.convert_to<double>(),
                                  moments[2].convert_to<double>(),
                                  tol );

  pb_evaluator->
    evaluateScreenedRutherfordMoment( rutherford_moment,
                                      energy,
                                      3 );

  UTILITY_TEST_FLOATING_EQUALITY( rutherford_moment.convert_to<double>(),
                                  moments[3].convert_to<double>(),
                                  tol );

  pb_evaluator->
    evaluateScreenedRutherfordMoment( rutherford_moment,
                                      energy,
                                      4 );

  UTILITY_TEST_FLOATING_EQUALITY( rutherford_moment.convert_to<double>(),
                                  moments[4].convert_to<double>(),
                                  tol );

  pb_evaluator->
    evaluateScreenedRutherfordMoment( rutherford_moment,
                                      energy,
                                      5 );

  UTILITY_TEST_FLOATING_EQUALITY( rutherford_moment.convert_to<double>(),
                                  moments[5].convert_to<double>(),
                                  tol );

  pb_evaluator->
    evaluateScreenedRutherfordMoment( rutherford_moment,
                                      energy,
                                      n );

  UTILITY_TEST_FLOATING_EQUALITY( rutherford_moment.convert_to<double>(),
                                  moments[6].convert_to<double>(),
                                  tol );

  tol = 1e-15;
  n = 0;
  moments.resize(n+1);
  moments[0] = 0.0;

  energy = 1.0e-5;
  pb_evaluator->
    evaluateScreenedRutherfordMoment( rutherford_moment,
                                      energy,
                                      0 );

  UTILITY_TEST_FLOATING_EQUALITY( rutherford_moment.convert_to<double>(),
                                  moments[0].convert_to<double>(),
                                  tol );

  energy = 6.30957;
  pb_evaluator->
    evaluateScreenedRutherfordMoment( rutherford_moment,
                                      energy,
                                      0 );

  UTILITY_TEST_FLOATING_EQUALITY( rutherford_moment.convert_to<double>(),
                                  moments[0].convert_to<double>(),
                                  tol );

  n = 2;
  energy = 1.0e5;
  moments.resize(n+1);

  cross_section = 3.8550699778229999356e5L;
  moments[0] = cross_section;
  moments[1] = 0.99999999999989729809L*cross_section;
  moments[2] = 0.99999999999969189428L*cross_section;

  al_evaluator->
    evaluateScreenedRutherfordMoment( rutherford_moment,
                                      energy,
                                      0 );

  UTILITY_TEST_FLOATING_EQUALITY( rutherford_moment.convert_to<double>(),
                                  moments[0].convert_to<double>(),
                                  tol );

  al_evaluator->
    evaluateScreenedRutherfordMoment( rutherford_moment,
                                      energy,
                                      1 );

  UTILITY_TEST_FLOATING_EQUALITY( rutherford_moment.convert_to<double>(),
                                  moments[1].convert_to<double>(),
                                  tol );

  al_evaluator->
    evaluateScreenedRutherfordMoment( rutherford_moment,
                                      energy,
                                      2 );

  UTILITY_TEST_FLOATING_EQUALITY( rutherford_moment.convert_to<double>(),
                                  moments[2].convert_to<double>(),
                                  tol );

  energy = 4.3750e1;

  tol = 1e-15;
  cross_section = 3.7500648271852370817e5L;
  moments.clear();
  moments[0] = cross_section;
  moments[1] = 0.99999992222034858054L*cross_section;
  moments[2] = 0.99999976666108267169L*cross_section;

  al_evaluator->
    evaluateScreenedRutherfordMoment( rutherford_moment,
                                      energy,
                                      0 );

  UTILITY_TEST_FLOATING_EQUALITY( rutherford_moment.convert_to<double>(),
                                  moments[0].convert_to<double>(),
                                  tol );

  al_evaluator->
    evaluateScreenedRutherfordMoment( rutherford_moment,
                                      energy,
                                      1 );

  UTILITY_TEST_FLOATING_EQUALITY( rutherford_moment.convert_to<double>(),
                                  moments[1].convert_to<double>(),
                                  tol );

  al_evaluator->
    evaluateScreenedRutherfordMoment( rutherford_moment,
                                      energy,
                                      2 );

  UTILITY_TEST_FLOATING_EQUALITY( rutherford_moment.convert_to<double>(),
                                  moments[2].convert_to<double>(),
                                  tol );

  energy = 1.0e-5;

  cross_section = 0.0L;
  moments.clear();
  moments[0] = cross_section;
  moments[1] = 9.50e-1L*cross_section;
  moments[2] = 8.550000000001e-1L*cross_section;

  al_evaluator->
    evaluateScreenedRutherfordMoment( rutherford_moment,
                                      energy,
                                      0 );

  UTILITY_TEST_FLOATING_EQUALITY( rutherford_moment.convert_to<double>(),
                                  moments[0].convert_to<double>(),
                                  tol );

  al_evaluator->
    evaluateScreenedRutherfordMoment( rutherford_moment,
                                      energy,
                                      1 );

  UTILITY_TEST_FLOATING_EQUALITY( rutherford_moment.convert_to<double>(),
                                  moments[1].convert_to<double>(),
                                  tol );

  al_evaluator->
    evaluateScreenedRutherfordMoment( rutherford_moment,
                                      energy,
                                      2 );

  UTILITY_TEST_FLOATING_EQUALITY( rutherford_moment.convert_to<double>(),
                                  moments[2].convert_to<double>(),
                                  tol );

}

//---------------------------------------------------------------------------//
// Check that the moments can be evaluated
TEUCHOS_UNIT_TEST( ElasticElectronMomentsEvaluator,
                   evaluateCutoffPDFMoment )
{
  unsigned n = 2;
  double energy = 1.0e5;
  double tol = 1e-13;
  std::vector<Utility::long_float> moments(n+1);
  Utility::long_float cutoff_moment;

  moments[0] = 0.99999999999990909993L;
  moments[1] = 0.99998613909145189036L;
  moments[2] = 0.99996009701780611083L;

  std::vector<double> angular_grid =
        MonteCarlo::ElasticElectronScatteringDistributionNativeFactory::getAngularGridAboveCutoff(
            al_data->getCutoffElasticAngles(),
            energy,
            -1.0 );

  Utility::GaussKronrodIntegrator<Utility::long_float> integrator( 1e-13 );

  test_al_evaluator->evaluateCutoffPDFMoment(
                        cutoff_moment,
                        angular_grid,
                        integrator,
                        energy,
                        0 );

  UTILITY_TEST_FLOATING_EQUALITY( cutoff_moment.convert_to<double>(),
                                  moments[0].convert_to<double>(),
                                  tol );

  test_al_evaluator->evaluateCutoffPDFMoment(
                        cutoff_moment,
                        angular_grid,
                        integrator,
                        energy,
                        1 );

  UTILITY_TEST_FLOATING_EQUALITY( cutoff_moment.convert_to<double>(),
                                  moments[1].convert_to<double>(),
                                  tol );

  test_al_evaluator->evaluateCutoffPDFMoment(
                        cutoff_moment,
                        angular_grid,
                        integrator,
                        energy,
                        2 );

  UTILITY_TEST_FLOATING_EQUALITY( cutoff_moment.convert_to<double>(),
                                  moments[2].convert_to<double>(),
                                  tol );
}

//---------------------------------------------------------------------------//
// Check that the moments can be evaluated
TEUCHOS_UNIT_TEST( ElasticElectronMomentsEvaluator,
                   evaluateCutoffMoment )
{
  unsigned n = 2;
  double energy = 1.0e5;
  double tol = 1e-13;
  std::vector<Utility::long_float> moments(n+1);
  Utility::long_float cutoff_moment;

  long double cross_section = 2.21769999999999997e-3L;

  moments[0] = 0.99999999999990909993L*cross_section;
  moments[1] = 0.99998613909145189036L*cross_section;
  moments[2] = 0.99996009701780611083L*cross_section;

  std::vector<double> angular_grid =
    MonteCarlo::ElasticElectronScatteringDistributionNativeFactory::getAngularGridAboveCutoff(
            al_data->getCutoffElasticAngles(),
            energy,
            -1.0 );

  Utility::GaussKronrodIntegrator<Utility::long_float> integrator( 1e-13 );

  test_al_evaluator->evaluateCutoffMoment(
                        cutoff_moment,
                        angular_grid,
                        integrator,
                        energy,
                        0 );

  UTILITY_TEST_FLOATING_EQUALITY( cutoff_moment.convert_to<double>(),
                                  moments[0].convert_to<double>(),
                                  tol );

  test_al_evaluator->evaluateCutoffMoment(
                        cutoff_moment,
                        angular_grid,
                        integrator,
                        energy,
                        1 );

  UTILITY_TEST_FLOATING_EQUALITY( cutoff_moment.convert_to<double>(),
                                  moments[1].convert_to<double>(),
                                  tol );

  test_al_evaluator->evaluateCutoffMoment(
                        cutoff_moment,
                        angular_grid,
                        integrator,
                        energy,
                        2 );

  UTILITY_TEST_FLOATING_EQUALITY( cutoff_moment.convert_to<double>(),
                                  moments[2].convert_to<double>(),
                                  tol );
}


//---------------------------------------------------------------------------//
// Check that the moments can be evaluated
TEUCHOS_UNIT_TEST( ElasticElectronMomentsEvaluator,
                   evaluateElasticMoment )
{
  double precision = 1e-10;
  double n = 0;
  double energy = 1.0e-5;
  double tol = 1e-10;

  std::shared_ptr<DataGen::ElasticElectronMomentsEvaluator> full_evaluator;

  full_evaluator.reset( new DataGen::ElasticElectronMomentsEvaluator(
                                    *pb_data,
                                    -1.0,
                                    1e-7,
                                    true ) );

  std::vector<Utility::long_float> total_moments(n+1);

  full_evaluator->evaluateElasticMoment( total_moments, energy, n, precision );

  double cutoff_cross_section = 2.48924e9;
  double sr_cross_section = 0.0;

  UTILITY_TEST_FLOATING_EQUALITY( total_moments[0].convert_to<double>(),
                                  1.0,
                                  tol );

  energy = 0.001;
  full_evaluator->evaluateElasticMoment( total_moments, energy, n, precision );


  cutoff_cross_section = 2.90281e8;
  sr_cross_section = 0.0;

  UTILITY_TEST_FLOATING_EQUALITY( total_moments[0].convert_to<double>(),
                                  1.0,
                                  tol );

  energy = 5.5e1;
  full_evaluator->evaluateElasticMoment( total_moments, energy, n, precision );

  cutoff_cross_section = 2.51776896157372e5;
  sr_cross_section = 1.89076556817226e6;

  UTILITY_TEST_FLOATING_EQUALITY( total_moments[0].convert_to<double>(),
                                  1.0,
                                  tol );
}

//---------------------------------------------------------------------------//
// Custom setup
//---------------------------------------------------------------------------//
UTILITY_CUSTOM_TEUCHOS_UNIT_TEST_SETUP_BEGIN();

  std::string test_native_pb_file_name, test_native_al_file_name;

UTILITY_CUSTOM_TEUCHOS_UNIT_TEST_COMMAND_LINE_OPTIONS()
{
  clp().setOption( "test_native_pb_file",
                   &test_native_pb_file_name,
                   "Test Native Pb file name" );
  clp().setOption( "test_native_al_file",
                   &test_native_al_file_name,
                   "Test Native Al file name" );
}

UTILITY_CUSTOM_TEUCHOS_UNIT_TEST_DATA_INITIALIZATION()
{
  double cutoff_angle_cosine = 0.9;
  double tabular_evaluation_tol = 1e-15;

  {
  // Create the native data file container
  al_data.reset( new Data::ElectronPhotonRelaxationDataContainer(
                                                   test_native_al_file_name ) );

  // Create the analog elastic distributions (combined Cutoff and Screened Rutherford)
  std::shared_ptr<const MonteCarlo::AnalogElasticElectronScatteringDistribution>
    analog_distribution;

  MonteCarlo::ElasticElectronScatteringDistributionNativeFactory::createAnalogElasticDistribution<Utility::LinLinLog>(
        analog_distribution,
        al_data->getCutoffElasticAngles(),
        al_data->getCutoffElasticPDF(),
        al_data->getElasticAngularEnergyGrid(),
        al_data->getAtomicNumber(),
        true,
        tabular_evaluation_tol );

  // Construct the hash-based grid searcher for this atom
  Teuchos::ArrayRCP<double> energy_grid;
  energy_grid.assign( al_data->getElectronEnergyGrid().begin(),
                      al_data->getElectronEnergyGrid().end() );

  Teuchos::RCP<Utility::HashBasedGridSearcher> grid_searcher(
     new Utility::StandardHashBasedGridSearcher<Teuchos::ArrayRCP<const double>, false>(
                                 energy_grid,
                                 100u ) );

  // Construct the cutoff reaction
  Teuchos::ArrayRCP<double> cutoff_cross_section;
  cutoff_cross_section.assign(
    al_data->getCutoffElasticCrossSection().begin(),
    al_data->getCutoffElasticCrossSection().end() );

  // Construct the screened Rutherford reaction
  Teuchos::ArrayRCP<double> rutherford_cross_section;
  rutherford_cross_section.assign(
    al_data->getScreenedRutherfordElasticCrossSection().begin(),
    al_data->getScreenedRutherfordElasticCrossSection().end() );

  // Create the moment evaluator
  al_evaluator.reset( new DataGen::ElasticElectronMomentsEvaluator(
        al_data->getCutoffElasticAngles(),
        energy_grid,
        grid_searcher,
        cutoff_cross_section,
        rutherford_cross_section,
        al_data->getCutoffElasticCrossSectionThresholdEnergyIndex(),
        al_data->getScreenedRutherfordElasticCrossSectionThresholdEnergyIndex(),
        analog_distribution,
        cutoff_angle_cosine ) );


  // Create the test moment evaluator
  test_al_evaluator.reset( new TestElasticElectronMomentsEvaluator(
                                    *al_data ) );
  }

  {
  // Create the native data file container
  pb_data.reset( new Data::ElectronPhotonRelaxationDataContainer(
                                    test_native_pb_file_name ) );

  bool linlinlog_interpolation_mode_on = true;

  // Create the moment evaluator
  pb_evaluator.reset(
    new DataGen::ElasticElectronMomentsEvaluator( *pb_data,
                                                  cutoff_angle_cosine,
                                                  tabular_evaluation_tol,
                                                  linlinlog_interpolation_mode_on ) );

  linlinlog_interpolation_mode_on = false;

  // Create the moment evaluator
  pb_lin_evaluator.reset(
    new DataGen::ElasticElectronMomentsEvaluator( *pb_data,
                                                  cutoff_angle_cosine,
                                                  tabular_evaluation_tol,
                                                  linlinlog_interpolation_mode_on ) );
  }
}

UTILITY_CUSTOM_TEUCHOS_UNIT_TEST_SETUP_END();

//---------------------------------------------------------------------------//
// end tstElasticElectronMomentsEvaluator.cpp
//---------------------------------------------------------------------------//