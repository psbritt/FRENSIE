//---------------------------------------------------------------------------//
//!
//! \file   tstAnalogElasticDistribution.cpp
//! \author Luke Kersting
//! \brief  Analog elastic distribution unit tests.
//!
//---------------------------------------------------------------------------//

// Std Lib Includes
#include <iostream>

// Boost Includes
#include <boost/units/systems/si.hpp>
#include <boost/units/systems/si/dimensionless.hpp>
#include <boost/units/systems/cgs.hpp>
#include <boost/units/io.hpp>

// Trilinos Includes
#include <Teuchos_UnitTestHarness.hpp>
#include <Teuchos_RCP.hpp>
#include <Teuchos_Array.hpp>
#include <Teuchos_ParameterList.hpp>
#include <Teuchos_XMLParameterListCoreHelpers.hpp>
#include <Teuchos_VerboseObject.hpp>

// FRENSIE Includes
#include "Utility_UnitTestHarnessExtensions.hpp"
#include "Utility_TabularOneDDistribution.hpp"
#include "Utility_AnalogElasticDistribution.hpp"
#include "Utility_RandomNumberGenerator.hpp"
#include "Utility_PhysicalConstants.hpp"
#include "Utility_UnitTraits.hpp"
#include "Utility_QuantityTraits.hpp"
#include "Utility_ElectronVoltUnit.hpp"
#include "Utility_AngleCosineUnit.hpp"

using boost::units::quantity;
using namespace Utility::Units;
namespace si = boost::units::si;
namespace cgs = boost::units::cgs;

typedef quantity<si::dimensionless> dl;

//---------------------------------------------------------------------------//
// Testing Variables
//---------------------------------------------------------------------------//

Teuchos::RCP<Teuchos::ParameterList> test_dists_list;

Teuchos::RCP<Utility::OneDDistribution> distribution;
Teuchos::RCP<Utility::TabularOneDDistribution> tab_distribution;

Teuchos::RCP<Utility::UnitAwareOneDDistribution<si::dimensionless,si::amount> >
  unit_aware_distribution;
Teuchos::RCP<Utility::UnitAwareTabularOneDDistribution<void,si::amount> >
unit_aware_tab_distribution;

unsigned atomic_number = 1u;

//---------------------------------------------------------------------------//
// Instantiation Macros.
//---------------------------------------------------------------------------//
#define UNIT_TEST_INSTANTIATION( type, name )                                \
  typedef Utility::LinLin LinLin;                                        \
  typedef Utility::LogLin LogLin;                                        \
  typedef Utility::LinLog LinLog;                                        \
  typedef Utility::LogLog LogLog;                                        \
  TEUCHOS_UNIT_TEST_TEMPLATE_1_INSTANT( type, name, LinLin )                \
  TEUCHOS_UNIT_TEST_TEMPLATE_1_INSTANT( type, name, LogLin )                \
  TEUCHOS_UNIT_TEST_TEMPLATE_1_INSTANT( type, name, LinLog )                \
  TEUCHOS_UNIT_TEST_TEMPLATE_1_INSTANT( type, name, LogLog )

//---------------------------------------------------------------------------//
// Testing Functions.
//---------------------------------------------------------------------------//
// Initialize the distribution
template<typename InterpolationPolicy, typename BaseDistribution>
void initialize( Teuchos::RCP<BaseDistribution>& dist )
{
  // Use the basic constructor
  Teuchos::Array<typename BaseDistribution::IndepQuantity>
    independent_values( 4 );
  Utility::setQuantity( independent_values[0], 1e-3 );
  Utility::setQuantity( independent_values[1], 1e-2 );
  Utility::setQuantity( independent_values[2], 1e-1 );
  Utility::setQuantity( independent_values[3], 1.0 );

  Teuchos::Array<typename BaseDistribution::DepQuantity> dependent_values( 4 );
  Utility::setQuantity( dependent_values[0], 1e2 );
  Utility::setQuantity( dependent_values[1], 1e1 );
  Utility::setQuantity( dependent_values[2], 1.0 );
  Utility::setQuantity( dependent_values[3], 1e-1 );

  dist.reset(new Utility::UnitAwareAnalogElasticDistribution<InterpolationPolicy,typename BaseDistribution::IndepUnit, typename BaseDistribution::DepUnit>(
                                                          independent_values,
                                                          dependent_values,
                                                          atomic_number ) );
}


// Initialize the distribution with a max CDF value
template<typename InterpolationPolicy, typename BaseDistribution>
void initializeWithMaxCDF( Teuchos::RCP<BaseDistribution>& dist )
{
  // Use the constructor with max CDF specified
  Teuchos::Array<typename BaseDistribution::IndepQuantity>
    independent_values( 4 );
  Utility::setQuantity( independent_values[0], 1e-3 );
  Utility::setQuantity( independent_values[1], 1e-2 );
  Utility::setQuantity( independent_values[2], 1e-1 );
  Utility::setQuantity( independent_values[3], 1.0 );

  Teuchos::Array<typename BaseDistribution::DepQuantity> dependent_values( 4 );
  Utility::setQuantity( dependent_values[0], 1e2 );
  Utility::setQuantity( dependent_values[1], 1e1 );
  Utility::setQuantity( dependent_values[2], 1.0 );
  Utility::setQuantity( dependent_values[3], 1e-1 );

////  UnnormCDFQuantity
//  double max_cdf = 10.0;

//  dist.reset(new Utility::UnitAwareAnalogElasticDistribution<InterpolationPolicy,typename BaseDistribution::IndepUnit, typename BaseDistribution::DepUnit>(
//                                                          independent_values,
//                                                          dependent_values,
//                                                          max_cdf ) );

  dist.reset(new Utility::UnitAwareAnalogElasticDistribution<InterpolationPolicy,typename BaseDistribution::IndepUnit, typename BaseDistribution::DepUnit>(
                                                          independent_values,
                                                          dependent_values,
                                                          atomic_number ) );
}

//---------------------------------------------------------------------------//
// Tests.
//---------------------------------------------------------------------------//
// Check that the distribution can be evaluated
TEUCHOS_UNIT_TEST_TEMPLATE_1_DECL( AnalogElasticDistribution,
                                   evaluate,
                                   InterpolationPolicy )
{
  initialize<InterpolationPolicy>( distribution );

  TEST_EQUALITY_CONST( distribution->evaluate( 0.0 ), 0.0 );
  TEST_EQUALITY_CONST( distribution->evaluate( 1e-3 ), 1e2 );
  TEST_EQUALITY_CONST( distribution->evaluate( 1e-2 ), 1e1 );
  TEST_EQUALITY_CONST( distribution->evaluate( 1e-1 ), 1.0 );
  TEST_EQUALITY_CONST( distribution->evaluate( 1.0 ), 1e-1 );
  TEST_EQUALITY_CONST( distribution->evaluate( 2.0 ), 0.0 );

  // Initialize the distribution with max CDF specified
  initializeWithMaxCDF<InterpolationPolicy>( distribution );

  TEST_EQUALITY_CONST( distribution->evaluate( 0.0 ), 0.0 );
  TEST_EQUALITY_CONST( distribution->evaluate( 1e-3 ), 1e2 );
  TEST_EQUALITY_CONST( distribution->evaluate( 1e-2 ), 1e1 );
  TEST_EQUALITY_CONST( distribution->evaluate( 1e-1 ), 1.0 );
  TEST_EQUALITY_CONST( distribution->evaluate( 1.0 ), 1e-1 );
  TEST_EQUALITY_CONST( distribution->evaluate( 2.0 ), 0.0 );

}

UNIT_TEST_INSTANTIATION( AnalogElasticDistribution, evaluate );

//---------------------------------------------------------------------------//
// Check that the unit-aware distribution can be evaluated
TEUCHOS_UNIT_TEST_TEMPLATE_1_DECL( UnitAwareAnalogElasticDistribution,
                                   evaluate,
                                   InterpolationPolicy )
{
  initialize<InterpolationPolicy>( unit_aware_distribution );

  TEST_EQUALITY_CONST( unit_aware_distribution->evaluate( 0.0 ),
                       0.0*si::mole );
  TEST_EQUALITY_CONST( unit_aware_distribution->evaluate( 1e-3 ),
                       1e2*si::mole );
  TEST_EQUALITY_CONST( unit_aware_distribution->evaluate( 1e-2 ),
                       1e1*si::mole );
  TEST_EQUALITY_CONST( unit_aware_distribution->evaluate( 1e-1 ),
                       1.0*si::mole );
  TEST_EQUALITY_CONST( unit_aware_distribution->evaluate( 1.0 ),
                       1e-1*si::mole );
  TEST_EQUALITY_CONST( unit_aware_distribution->evaluate( 2.0 ),
                       0.0*si::mole );
}

UNIT_TEST_INSTANTIATION( UnitAwareAnalogElasticDistribution, evaluate );

//---------------------------------------------------------------------------//
// Check that the PDF can be evaluated
TEUCHOS_UNIT_TEST_TEMPLATE_1_DECL( AnalogElasticDistribution,
                                   evaluatePDF,
                                   InterpolationPolicy )
{
  initialize<InterpolationPolicy>( distribution );

  TEST_EQUALITY_CONST( distribution->evaluatePDF( 0.0 ), 0.0 );
  TEST_FLOATING_EQUALITY( distribution->evaluatePDF( 1e-3 ),
                          67.340006734,
                          1e-6 );
  TEST_FLOATING_EQUALITY( distribution->evaluatePDF( 1e-2 ),
                          6.7340006734,
                          1e-6 );
  TEST_FLOATING_EQUALITY( distribution->evaluatePDF( 1e-1 ),
                          0.67340006734,
                          1e-6 );
  TEST_FLOATING_EQUALITY( distribution->evaluatePDF( 1.0 ),
                          0.067340006734,
                          1e-6 );
  TEST_EQUALITY_CONST( distribution->evaluatePDF( 2.0 ), 0.0 );

  // Initialize the distribution with max CDF specified
  initializeWithMaxCDF<InterpolationPolicy>( distribution );

  TEST_EQUALITY_CONST( distribution->evaluatePDF( 0.0 ), 0.0 );
  TEST_FLOATING_EQUALITY( distribution->evaluatePDF( 1e-3 ),
                          10.0,
                          1e-6 );
  TEST_FLOATING_EQUALITY( distribution->evaluatePDF( 1e-2 ),
                          1.0,
                          1e-6 );
  TEST_FLOATING_EQUALITY( distribution->evaluatePDF( 1e-1 ),
                          0.1,
                          1e-6 );
  TEST_FLOATING_EQUALITY( distribution->evaluatePDF( 1.0 ),
                          0.01,
                          1e-6 );
  TEST_EQUALITY_CONST( distribution->evaluatePDF( 2.0 ), 0.0 );
}

UNIT_TEST_INSTANTIATION( AnalogElasticDistribution, evaluatePDF );

//---------------------------------------------------------------------------//
// Check that the unit-aware PDF can be evaluated
TEUCHOS_UNIT_TEST_TEMPLATE_1_DECL( UnitAwareAnalogElasticDistribution,
                                   evaluatePDF,
                                   InterpolationPolicy )
{
  initialize<InterpolationPolicy>( unit_aware_distribution );

  TEST_EQUALITY_CONST( unit_aware_distribution->evaluatePDF( 0.0*si::dimensionless() ),
                       0.0/si::dimensionless() );
  UTILITY_TEST_FLOATING_EQUALITY(
                             unit_aware_distribution->evaluatePDF( 1e-3*si::dimensionless() ),
                             67.340006734/si::dimensionless(),
                             1e-6 );
  UTILITY_TEST_FLOATING_EQUALITY(
                             unit_aware_distribution->evaluatePDF( 1e-2*si::dimensionless() ),
                             6.7340006734/si::dimensionless(),
                             1e-6 );
  UTILITY_TEST_FLOATING_EQUALITY(
                             unit_aware_distribution->evaluatePDF( 1e-1*si::dimensionless() ),
                             0.67340006734/si::dimensionless(),
                             1e-6 );
  UTILITY_TEST_FLOATING_EQUALITY(
                              unit_aware_distribution->evaluatePDF( 1.0*si::dimensionless() ),
                              0.067340006734/si::dimensionless(),
                              1e-6 );
  TEST_EQUALITY_CONST( unit_aware_distribution->evaluatePDF( 2.0*si::dimensionless() ),
                       0.0/si::dimensionless() );
}

UNIT_TEST_INSTANTIATION( UnitAwareAnalogElasticDistribution, evaluatePDF );

//---------------------------------------------------------------------------//
// Check that the CDF can be evaluated
TEUCHOS_UNIT_TEST_TEMPLATE_1_DECL( AnalogElasticDistribution,
                                   evaluateCDF,
                                   InterpolationPolicy )
{
  initialize<InterpolationPolicy>( tab_distribution );

  TEST_EQUALITY_CONST( tab_distribution->evaluateCDF( 0.0 ), 0.0 );
  TEST_FLOATING_EQUALITY( tab_distribution->evaluateCDF( 1e-3 ),
                          0.0000000000,
                          1e-10 );
  TEST_FLOATING_EQUALITY( tab_distribution->evaluateCDF( 1e-2 ),
                          0.33333333333,
                          1e-10 );
  TEST_FLOATING_EQUALITY( tab_distribution->evaluateCDF( 1e-1 ),
                          0.66666666667,
                          1e-10 );
  TEST_FLOATING_EQUALITY( tab_distribution->evaluateCDF( 1.0 ),
                          1.0000000000,
                          1e-10 );
  TEST_EQUALITY_CONST( tab_distribution->evaluateCDF( 2.0 ), 1.0 );

  // Initialize the distribution with max CDF specified
  initializeWithMaxCDF<InterpolationPolicy>( tab_distribution );

  TEST_EQUALITY_CONST( tab_distribution->evaluateCDF( 0.0 ), 0.0 );
  TEST_FLOATING_EQUALITY( tab_distribution->evaluateCDF( 1e-3 ),
                          0.0000000000,
                          1e-10 );
  TEST_FLOATING_EQUALITY( tab_distribution->evaluateCDF( 1e-2 ),
                          0.0495,
                          1e-10 );
  TEST_FLOATING_EQUALITY( tab_distribution->evaluateCDF( 1e-1 ),
                          0.099,
                          1e-10 );
  TEST_FLOATING_EQUALITY( tab_distribution->evaluateCDF( 1.0 ),
                          1.0000000000,
                          1e-10 );
  TEST_EQUALITY_CONST( tab_distribution->evaluateCDF( 2.0 ), 1.0 );
}

UNIT_TEST_INSTANTIATION( AnalogElasticDistribution, evaluateCDF );

//---------------------------------------------------------------------------//
// Check that the unit-aware CDF can be evaluated
TEUCHOS_UNIT_TEST_TEMPLATE_1_DECL( UnitAwareAnalogElasticDistribution,
                                   evaluateCDF,
                                   InterpolationPolicy )
{
  initialize<InterpolationPolicy>( unit_aware_tab_distribution );

  TEST_EQUALITY_CONST( unit_aware_tab_distribution->evaluateCDF( 0.0*si::dimensionless() ),
                       0.0 );
  TEST_FLOATING_EQUALITY( unit_aware_tab_distribution->evaluateCDF( 1e-3*si::dimensionless() ),
                          0.0000000000,
                          1e-10 );
  TEST_FLOATING_EQUALITY( unit_aware_tab_distribution->evaluateCDF( 1e-2*si::dimensionless() ),
                          0.33333333333,
                          1e-10 );
  TEST_FLOATING_EQUALITY( unit_aware_tab_distribution->evaluateCDF( 1e-1*si::dimensionless() ),
                          0.66666666667,
                          1e-10 );
  TEST_FLOATING_EQUALITY( unit_aware_tab_distribution->evaluateCDF( 1.0*si::dimensionless() ),
                          1.0000000000,
                          1e-10 );
  TEST_EQUALITY_CONST( unit_aware_tab_distribution->evaluateCDF( 2.0*si::dimensionless() ),
                       1.0 );
}

UNIT_TEST_INSTANTIATION( UnitAwareAnalogElasticDistribution, evaluateCDF );

//---------------------------------------------------------------------------//
// Check that the distribution can be sampled
TEUCHOS_UNIT_TEST_TEMPLATE_1_DECL( AnalogElasticDistribution,
                                   sample,
                                   InterpolationPolicy )
{
  initialize<InterpolationPolicy>( distribution );

  std::vector<double> fake_stream( 2 );
  fake_stream[0] = 0.0;
  fake_stream[1] = 1.0 - 1e-15;

  Utility::RandomNumberGenerator::setFakeStream( fake_stream );

  double sample = distribution->sample();
  TEST_EQUALITY_CONST( sample, 1e-3 );

  sample = distribution->sample();
  TEST_FLOATING_EQUALITY( sample, 1.0, 1e-12 );

  Utility::RandomNumberGenerator::unsetFakeStream();
  Utility::RandomNumberGenerator::initialize();

  sample = distribution->sample();
  TEST_COMPARE( sample, >=, 1e-3 );
  TEST_COMPARE( sample, <=, 1.0 );


  // Initialize the distribution with max CDF specified
  initializeWithMaxCDF<InterpolationPolicy>( distribution );

  fake_stream.resize( 3 );
  fake_stream[0] = 0.0;
  fake_stream[1] = 0.1;
  fake_stream[2] = 1.0 - 1e-15;

  Utility::RandomNumberGenerator::setFakeStream( fake_stream );

  sample = distribution->sample();
  TEST_EQUALITY_CONST( sample, 1e-3 );

  sample = distribution->sample();
  TEST_FLOATING_EQUALITY( sample, 0.11005050633883220468, 1e-12 );

  sample = distribution->sample();
  TEST_EQUALITY_CONST( sample, 0.0 );

  Utility::RandomNumberGenerator::unsetFakeStream();
  Utility::RandomNumberGenerator::initialize();

  sample = distribution->sample();
  TEST_COMPARE( sample, >=, 0.0 );
  TEST_COMPARE( sample, <=, 0.11005050633883220468 );
}

UNIT_TEST_INSTANTIATION( AnalogElasticDistribution, sample );

//---------------------------------------------------------------------------//
// Check that the unit-aware distribution can be sampled
TEUCHOS_UNIT_TEST_TEMPLATE_1_DECL( UnitAwareAnalogElasticDistribution,
                                   sample,
                                   InterpolationPolicy )
{
  initialize<InterpolationPolicy>( unit_aware_distribution );

  std::vector<double> fake_stream( 2 );
  fake_stream[0] = 0.0;
  fake_stream[1] = 1.0 - 1e-15;

  Utility::RandomNumberGenerator::setFakeStream( fake_stream );

  quantity<si::dimensionless> sample = unit_aware_distribution->sample();
  TEST_EQUALITY_CONST( sample, 1e-3*si::dimensionless() );

  sample = unit_aware_distribution->sample();
  UTILITY_TEST_FLOATING_EQUALITY( sample, 1.0*si::dimensionless(), 1e-12 );

  Utility::RandomNumberGenerator::unsetFakeStream();
  Utility::RandomNumberGenerator::initialize();

  sample = unit_aware_distribution->sample();
  TEST_COMPARE( sample, >=, 1e-3*si::dimensionless() );
  TEST_COMPARE( sample, <=, 1.0*si::dimensionless() );
}

UNIT_TEST_INSTANTIATION( UnitAwareAnalogElasticDistribution, sample );

//---------------------------------------------------------------------------//
// Check that the distribution can be sampled
TEUCHOS_UNIT_TEST_TEMPLATE_1_DECL( AnalogElasticDistribution,
                                   sampleAndRecordTrials,
                                   InterpolationPolicy )
{
  initialize<InterpolationPolicy>( distribution );

  std::vector<double> fake_stream( 2 );
  fake_stream[0] = 0.0;
  fake_stream[1] = 1.0 - 1e-15;

  Utility::RandomNumberGenerator::setFakeStream( fake_stream );

  unsigned trials = 0;

  double sample = distribution->sampleAndRecordTrials( trials );
  TEST_EQUALITY_CONST( sample, 1e-3 );
  TEST_EQUALITY_CONST( 1.0/trials, 1.0 );

  sample = distribution->sampleAndRecordTrials( trials );
  TEST_FLOATING_EQUALITY( sample, 1.0, 1e-12 );
  TEST_EQUALITY_CONST( 2.0/trials, 1.0 );

  Utility::RandomNumberGenerator::unsetFakeStream();
  Utility::RandomNumberGenerator::initialize();

  sample = distribution->sampleAndRecordTrials( trials );
  TEST_COMPARE( sample, >=, 1e-3 );
  TEST_COMPARE( sample, <=, 1.0 );
  TEST_EQUALITY_CONST( 3.0/trials, 1.0 );

  // Initialize the distribution with max CDF specified
  initializeWithMaxCDF<InterpolationPolicy>( distribution );

  fake_stream.resize( 3 );
  fake_stream[0] = 0.0;
  fake_stream[1] = 0.1;
  fake_stream[2] = 1.0 - 1e-15;

  Utility::RandomNumberGenerator::setFakeStream( fake_stream );

  sample = distribution->sampleAndRecordTrials( trials );
  TEST_EQUALITY_CONST( sample, 1e-3 );
  TEST_EQUALITY_CONST( 4.0/trials, 1.0 );

  sample = distribution->sampleAndRecordTrials( trials );
  TEST_FLOATING_EQUALITY( sample, 0.11005050633883220468, 1e-12 );
  TEST_EQUALITY_CONST( 5.0/trials, 1.0 );

  sample = distribution->sampleAndRecordTrials( trials );
  TEST_EQUALITY_CONST( sample, 0.0 );
  TEST_EQUALITY_CONST( 6.0/trials, 1.0 );

  Utility::RandomNumberGenerator::unsetFakeStream();
  Utility::RandomNumberGenerator::initialize();

  sample = distribution->sampleAndRecordTrials( trials );
  TEST_COMPARE( sample, >=, 0.0 );
  TEST_COMPARE( sample, <=, 0.11005050633883220468 );
  TEST_EQUALITY_CONST( 7.0/trials, 1.0 );
}

UNIT_TEST_INSTANTIATION( AnalogElasticDistribution, sampleAndRecordTrials );

//---------------------------------------------------------------------------//
// Check that the unit-aware distribution can be sampled
TEUCHOS_UNIT_TEST_TEMPLATE_1_DECL( UnitAwareAnalogElasticDistribution,
                                   sampleAndRecordTrials,
                                   InterpolationPolicy )
{
  initialize<InterpolationPolicy>( unit_aware_distribution );

  std::vector<double> fake_stream( 2 );
  fake_stream[0] = 0.0;
  fake_stream[1] = 1.0 - 1e-15;

  Utility::RandomNumberGenerator::setFakeStream( fake_stream );

  unsigned trials = 0;

  quantity<si::dimensionless> sample =
    unit_aware_distribution->sampleAndRecordTrials( trials );
  TEST_EQUALITY_CONST( sample, 1e-3*si::dimensionless() );
  TEST_EQUALITY_CONST( 1.0/trials, 1.0 );

  sample = unit_aware_distribution->sampleAndRecordTrials( trials );
  UTILITY_TEST_FLOATING_EQUALITY( sample, 1.0*si::dimensionless(), 1e-12 );
  TEST_EQUALITY_CONST( 2.0/trials, 1.0 );

  Utility::RandomNumberGenerator::unsetFakeStream();
  Utility::RandomNumberGenerator::initialize();

  sample = unit_aware_distribution->sampleAndRecordTrials( trials );
  TEST_COMPARE( sample, >=, 1e-3*si::dimensionless() );
  TEST_COMPARE( sample, <=, 1.0*si::dimensionless() );
  TEST_EQUALITY_CONST( 3.0/trials, 1.0 );
}

UNIT_TEST_INSTANTIATION( UnitAwareAnalogElasticDistribution, sampleAndRecordTrials );

//---------------------------------------------------------------------------//
// Check that the distribution can be sampled
TEUCHOS_UNIT_TEST_TEMPLATE_1_DECL( AnalogElasticDistribution,
                                   sampleAndRecordBinIndex,
                                   InterpolationPolicy )
{
  initialize<InterpolationPolicy>( tab_distribution );

  std::vector<double> fake_stream( 2 );
  fake_stream[0] = 0.0;
  fake_stream[1] = 1.0 - 1e-15;

  Utility::RandomNumberGenerator::setFakeStream( fake_stream );

  unsigned bin_index;

  double sample = tab_distribution->sampleAndRecordBinIndex( bin_index );
  TEST_EQUALITY_CONST( sample, 1e-3 );
  TEST_EQUALITY_CONST( bin_index, 0u );

  sample = tab_distribution->sampleAndRecordBinIndex( bin_index );
  TEST_FLOATING_EQUALITY( sample, 1.0, 1e-12 );
  TEST_EQUALITY_CONST( bin_index, 2u );

  Utility::RandomNumberGenerator::unsetFakeStream();
  Utility::RandomNumberGenerator::initialize();

  sample = tab_distribution->sampleAndRecordBinIndex( bin_index );
  TEST_COMPARE( sample, >=, 1e-3 );
  TEST_COMPARE( sample, <=, 1.0 );

  // Initialize the distribution with max CDF specified
  initializeWithMaxCDF<InterpolationPolicy>( tab_distribution );

  fake_stream.resize( 3 );
  fake_stream[0] = 0.0;
  fake_stream[1] = 0.1;
  fake_stream[2] = 1.0 - 1e-15;

  Utility::RandomNumberGenerator::setFakeStream( fake_stream );

  sample = tab_distribution->sampleAndRecordBinIndex( bin_index );
  TEST_EQUALITY_CONST( sample, 1e-3 );
  TEST_EQUALITY_CONST( bin_index, 0u );

  sample = tab_distribution->sampleAndRecordBinIndex( bin_index );
  TEST_FLOATING_EQUALITY( sample, 0.11005050633883220468, 1e-12 );
  TEST_EQUALITY_CONST( bin_index, 2u );

  sample = tab_distribution->sampleAndRecordBinIndex( bin_index );
  TEST_FLOATING_EQUALITY( sample, 0.0, 1e-12 );
  TEST_EQUALITY_CONST( bin_index, 2u );

  Utility::RandomNumberGenerator::unsetFakeStream();
  Utility::RandomNumberGenerator::initialize();

  sample = tab_distribution->sampleAndRecordBinIndex( bin_index );
  TEST_COMPARE( sample, >=, 0.0 );
  TEST_COMPARE( sample, <=, 0.11005050633883220468 );
}

UNIT_TEST_INSTANTIATION( AnalogElasticDistribution, sampleAndRecordBinIndex );

//---------------------------------------------------------------------------//
// Check that the unit-aware distribution can be sampled
TEUCHOS_UNIT_TEST_TEMPLATE_1_DECL( UnitAwareAnalogElasticDistribution,
                                   sampleAndRecordBinIndex,
                                   InterpolationPolicy )
{
  initialize<InterpolationPolicy>( unit_aware_tab_distribution );

  std::vector<double> fake_stream( 2 );
  fake_stream[0] = 0.0;
  fake_stream[1] = 1.0 - 1e-15;

  Utility::RandomNumberGenerator::setFakeStream( fake_stream );

  unsigned bin_index;

  quantity<si::dimensionless> sample =
    unit_aware_tab_distribution->sampleAndRecordBinIndex( bin_index );
  TEST_EQUALITY_CONST( sample, 1e-3*si::dimensionless() );
  TEST_EQUALITY_CONST( bin_index, 0u );

  sample = unit_aware_tab_distribution->sampleAndRecordBinIndex( bin_index );
  UTILITY_TEST_FLOATING_EQUALITY( sample, 1.0*si::dimensionless(), 1e-12 );
  TEST_EQUALITY_CONST( bin_index, 2u );

  Utility::RandomNumberGenerator::unsetFakeStream();
  Utility::RandomNumberGenerator::initialize();

  sample = unit_aware_tab_distribution->sampleAndRecordBinIndex( bin_index );
  TEST_COMPARE( sample, >=, 1e-3*si::dimensionless() );
  TEST_COMPARE( sample, <=, 1.0*si::dimensionless() );
}

UNIT_TEST_INSTANTIATION( UnitAwareAnalogElasticDistribution,
                         sampleAndRecordBinIndex );

//---------------------------------------------------------------------------//
// Check that the distribution can be sampled
TEUCHOS_UNIT_TEST_TEMPLATE_1_DECL( AnalogElasticDistribution,
                                   sampleWithRandomNumber,
                                   InterpolationPolicy )
{
  initialize<InterpolationPolicy>( tab_distribution );

  double sample = tab_distribution->sampleWithRandomNumber( 0.0 );
  TEST_EQUALITY_CONST( sample, 1e-3 );

  sample = tab_distribution->sampleWithRandomNumber( 1.0 - 1e-15 );
  TEST_FLOATING_EQUALITY( sample, 1.0, 1e-12 );


  // Initialize the distribution with max CDF specified
  initializeWithMaxCDF<InterpolationPolicy>( tab_distribution );

  sample = tab_distribution->sampleWithRandomNumber( 0.0 );
  TEST_EQUALITY_CONST( sample, 1e-3 );

  sample = tab_distribution->sampleWithRandomNumber( 0.1 );
  TEST_FLOATING_EQUALITY( sample, 0.11005050633883220468, 1e-12 );

  sample = tab_distribution->sampleWithRandomNumber( 1.0 - 1e-15 );
  TEST_FLOATING_EQUALITY( sample, 0.0, 1e-12 );
}

UNIT_TEST_INSTANTIATION( AnalogElasticDistribution, sampleWithRandomNumber );

//---------------------------------------------------------------------------//
// Check that the unit-aware distribution can be sampled
TEUCHOS_UNIT_TEST_TEMPLATE_1_DECL( UnitAwareAnalogElasticDistribution,
                                   sampleWithRandomNumber,
                                   InterpolationPolicy )
{
  initialize<InterpolationPolicy>( unit_aware_tab_distribution );

  quantity<si::dimensionless> sample =
    unit_aware_tab_distribution->sampleWithRandomNumber( 0.0 );
  TEST_EQUALITY_CONST( sample, 1e-3*si::dimensionless() );

  sample = unit_aware_tab_distribution->sampleWithRandomNumber( 1.0 - 1e-15 );
  UTILITY_TEST_FLOATING_EQUALITY( sample, 1.0*si::dimensionless(), 1e-12 );
}

UNIT_TEST_INSTANTIATION( UnitAwareAnalogElasticDistribution,
                         sampleWithRandomNumber );

//---------------------------------------------------------------------------//
// Check that the distribution can be sampled from a subrange
TEUCHOS_UNIT_TEST_TEMPLATE_1_DECL( AnalogElasticDistribution,
                                   sampleInSubrange,
                                   InterpolationPolicy )
{
  initialize<InterpolationPolicy>( tab_distribution );

  std::vector<double> fake_stream( 2 );
  fake_stream[0] = 0.0;
  fake_stream[1] = 1.0 - 1e-15;

  Utility::RandomNumberGenerator::setFakeStream( fake_stream );

  double sample = tab_distribution->sampleInSubrange( 1e-1  );
  TEST_EQUALITY_CONST( sample, 1e-3 );

  sample = tab_distribution->sampleInSubrange( 1e-1 );
  TEST_FLOATING_EQUALITY( sample, 1e-1, 1e-12 );

  Utility::RandomNumberGenerator::unsetFakeStream();
  Utility::RandomNumberGenerator::initialize();

  sample = tab_distribution->sampleInSubrange( 1e-1 );
  TEST_COMPARE( sample, >=, 1e-3 );
  TEST_COMPARE( sample, <=, 1e-1 );


  // Initialize the distribution with max CDF specified
  initializeWithMaxCDF<InterpolationPolicy>( tab_distribution );

  Utility::RandomNumberGenerator::setFakeStream( fake_stream );

  sample = tab_distribution->sampleInSubrange( 1e-1 );
  TEST_EQUALITY_CONST( sample, 1e-3 );

  sample = tab_distribution->sampleInSubrange( 1e-1 );
  TEST_FLOATING_EQUALITY( sample, 1e-1, 1e-12 );

  Utility::RandomNumberGenerator::unsetFakeStream();
  Utility::RandomNumberGenerator::initialize();

  sample = tab_distribution->sampleInSubrange( 1e-1 );
  TEST_COMPARE( sample, >=, 1e-3 );
  TEST_COMPARE( sample, <=, 1e-1 );
}

UNIT_TEST_INSTANTIATION( AnalogElasticDistribution, sampleInSubrange );

//---------------------------------------------------------------------------//
// Check that the unit-aware distribution can be sampled from a subrange
TEUCHOS_UNIT_TEST_TEMPLATE_1_DECL( UnitAwareAnalogElasticDistribution,
                                   sampleInSubrange,
                                   InterpolationPolicy )
{
  initialize<InterpolationPolicy>( unit_aware_tab_distribution );

  std::vector<double> fake_stream( 2 );
  fake_stream[0] = 0.0;
  fake_stream[1] = 1.0 - 1e-15;

  Utility::RandomNumberGenerator::setFakeStream( fake_stream );

  quantity<si::dimensionless> sample =
    unit_aware_tab_distribution->sampleInSubrange( 1e-1*si::dimensionless()  );
  TEST_EQUALITY_CONST( sample, 1e-3*si::dimensionless() );

  sample = unit_aware_tab_distribution->sampleInSubrange( 1e-1*si::dimensionless() );
  UTILITY_TEST_FLOATING_EQUALITY( sample, 1e-1*si::dimensionless(), 1e-12 );

  Utility::RandomNumberGenerator::unsetFakeStream();
  Utility::RandomNumberGenerator::initialize();

  sample = unit_aware_tab_distribution->sampleInSubrange( 1e-1*si::dimensionless() );
  TEST_COMPARE( sample, >=, 1e-3*si::dimensionless() );
  TEST_COMPARE( sample, <=, 1e-1*si::dimensionless() );
}

UNIT_TEST_INSTANTIATION( UnitAwareAnalogElasticDistribution, sampleInSubrange );

//---------------------------------------------------------------------------//
// Check that the distribution can be sampled from a subrange
TEUCHOS_UNIT_TEST_TEMPLATE_1_DECL( AnalogElasticDistribution,
                                   sampleWithRandomNumberInSubrange,
                                   InterpolationPolicy )
{
  initialize<InterpolationPolicy>( tab_distribution );

  double sample =
    tab_distribution->sampleWithRandomNumberInSubrange( 0.0, 1e-1  );
  TEST_EQUALITY_CONST( sample, 1e-3 );

  sample = tab_distribution->sampleWithRandomNumberInSubrange( 1.0, 1e-1 );
  TEST_FLOATING_EQUALITY( sample, 1e-1, 1e-12 );

  // Initialize the distribution with max CDF specified
  initializeWithMaxCDF<InterpolationPolicy>( tab_distribution );

  sample = tab_distribution->sampleWithRandomNumberInSubrange( 0.0, 1e-1  );
  TEST_EQUALITY_CONST( sample, 1e-3 );

  sample = tab_distribution->sampleWithRandomNumberInSubrange( 1.0, 1e-1 );
  TEST_FLOATING_EQUALITY( sample, 1e-1, 1e-12 );

}

UNIT_TEST_INSTANTIATION( AnalogElasticDistribution,
                         sampleWithRandomNumberInSubrange );

//---------------------------------------------------------------------------//
// Check that the unit-aware distribution can be sampled from a subrange
TEUCHOS_UNIT_TEST_TEMPLATE_1_DECL( UnitAwareAnalogElasticDistribution,
                                   sampleWithRandomNumberInSubrange,
                                   InterpolationPolicy )
{
  initialize<InterpolationPolicy>( unit_aware_tab_distribution );

  quantity<si::dimensionless> sample =
    unit_aware_tab_distribution->sampleWithRandomNumberInSubrange(
                                                               0.0, 1e-1*si::dimensionless() );
  TEST_EQUALITY_CONST( sample, 1e-3*si::dimensionless() );

  sample = unit_aware_tab_distribution->sampleWithRandomNumberInSubrange(
                                                               1.0, 1e-1*si::dimensionless() );
  UTILITY_TEST_FLOATING_EQUALITY( sample, 1e-1*si::dimensionless(), 1e-12 );
}

UNIT_TEST_INSTANTIATION( UnitAwareAnalogElasticDistribution,
                         sampleWithRandomNumberInSubrange );

//---------------------------------------------------------------------------//
// Check that the upper bound of the distribution independent variable can be
// returned
TEUCHOS_UNIT_TEST_TEMPLATE_1_DECL( AnalogElasticDistribution,
                                   getUpperBoundOfIndepVar,
                                   InterpolationPolicy )
{
  initialize<InterpolationPolicy>( distribution );

  TEST_EQUALITY_CONST( distribution->getUpperBoundOfIndepVar(), 1.0 );


  // Initialize the distribution with max CDF specified
  initializeWithMaxCDF<InterpolationPolicy>( tab_distribution );

  TEST_EQUALITY_CONST( distribution->getUpperBoundOfIndepVar(), 1.0 );
}

UNIT_TEST_INSTANTIATION( AnalogElasticDistribution, getUpperBoundOfIndepVar );

//---------------------------------------------------------------------------//
// Check that the upper bound of the unit-aware distribution independent
// variable can be returned
TEUCHOS_UNIT_TEST_TEMPLATE_1_DECL( UnitAwareAnalogElasticDistribution,
                                   getUpperBoundOfIndepVar,
                                   InterpolationPolicy )
{
  initialize<InterpolationPolicy>( unit_aware_distribution );

  TEST_EQUALITY_CONST( unit_aware_distribution->getUpperBoundOfIndepVar(),
                       1.0*si::dimensionless() );
}

UNIT_TEST_INSTANTIATION( UnitAwareAnalogElasticDistribution,
                         getUpperBoundOfIndepVar );

//---------------------------------------------------------------------------//
// Check that the lower bound of the distribution independent variable can be
// returned
TEUCHOS_UNIT_TEST_TEMPLATE_1_DECL( AnalogElasticDistribution,
                                   getLowerBoundOfIndepVar,
                                   InterpolationPolicy )
{
  initialize<InterpolationPolicy>( distribution );

  TEST_EQUALITY_CONST( distribution->getLowerBoundOfIndepVar(), 1e-3 );


  // Initialize the distribution with max CDF specified
  initializeWithMaxCDF<InterpolationPolicy>( tab_distribution );

  TEST_EQUALITY_CONST( distribution->getLowerBoundOfIndepVar(), 1e-3 );
}

UNIT_TEST_INSTANTIATION( AnalogElasticDistribution, getLowerBoundOfIndepVar );

//---------------------------------------------------------------------------//
// Check that the lower bound of the unit-aware distribution independent
// variable can be returned
TEUCHOS_UNIT_TEST_TEMPLATE_1_DECL( UnitAwareAnalogElasticDistribution,
                                   getLowerBoundOfIndepVar,
                                   InterpolationPolicy )
{
  initialize<InterpolationPolicy>( unit_aware_distribution );

  TEST_EQUALITY_CONST( unit_aware_distribution->getLowerBoundOfIndepVar(),
                       1e-3*si::dimensionless() );
}

UNIT_TEST_INSTANTIATION( UnitAwareAnalogElasticDistribution,
                         getLowerBoundOfIndepVar );

//---------------------------------------------------------------------------//
// Check that the distribution type can be returned
TEUCHOS_UNIT_TEST_TEMPLATE_1_DECL( AnalogElasticDistribution,
                                   getDistributionType,
                                   InterpolationPolicy )
{
  initialize<InterpolationPolicy>( distribution );

  TEST_EQUALITY_CONST( distribution->getDistributionType(),
                       Utility::ANALOG_ELASTIC_DISTRIBUTION );
}

UNIT_TEST_INSTANTIATION( AnalogElasticDistribution, getDistributionType );

//---------------------------------------------------------------------------//
// Check that the unit-aware distribution type can be returned
TEUCHOS_UNIT_TEST_TEMPLATE_1_DECL( UnitAwareAnalogElasticDistribution,
                                   getDistributionType,
                                   InterpolationPolicy )
{
  initialize<InterpolationPolicy>( unit_aware_distribution );

  TEST_EQUALITY_CONST( unit_aware_distribution->getDistributionType(),
                       Utility::ANALOG_ELASTIC_DISTRIBUTION );
}

UNIT_TEST_INSTANTIATION( UnitAwareAnalogElasticDistribution, getDistributionType );

//---------------------------------------------------------------------------//
// Check if the distribution is tabular
TEUCHOS_UNIT_TEST_TEMPLATE_1_DECL( AnalogElasticDistribution,
                                   isTabular,
                                   InterpolationPolicy )
{
  initialize<InterpolationPolicy>( distribution );

  TEST_ASSERT( distribution->isTabular() );
}

UNIT_TEST_INSTANTIATION( AnalogElasticDistribution, isTabular );

//---------------------------------------------------------------------------//
// Check if the unit-aware distribution is tabular
TEUCHOS_UNIT_TEST_TEMPLATE_1_DECL( UnitAwareAnalogElasticDistribution,
                                   isTabular,
                                   InterpolationPolicy )
{
  initialize<InterpolationPolicy>( unit_aware_distribution );

  TEST_ASSERT( unit_aware_distribution->isTabular() );
}

UNIT_TEST_INSTANTIATION( UnitAwareAnalogElasticDistribution, isTabular );

//---------------------------------------------------------------------------//
// Check that the distribution is continuous
TEUCHOS_UNIT_TEST_TEMPLATE_1_DECL( AnalogElasticDistribution,
                                   isContinuous,
                                   InterpolationPolicy )
{
  initialize<InterpolationPolicy>( distribution );

  TEST_ASSERT( distribution->isContinuous() );
}

UNIT_TEST_INSTANTIATION( AnalogElasticDistribution, isContinuous );

//---------------------------------------------------------------------------//
// Check that the unit-aware distribution is continuous
TEUCHOS_UNIT_TEST_TEMPLATE_1_DECL( UnitAwareAnalogElasticDistribution,
                                   isContinuous,
                                   InterpolationPolicy )
{
  initialize<InterpolationPolicy>( unit_aware_distribution );

  TEST_ASSERT( unit_aware_distribution->isContinuous() );
}

UNIT_TEST_INSTANTIATION( UnitAwareAnalogElasticDistribution, isContinuous );

//---------------------------------------------------------------------------//
// Check if the distribution is compatible with the interpolation type
TEUCHOS_UNIT_TEST_TEMPLATE_1_DECL( AnalogElasticDistribution,
                                   isCompatibleWithInterpType,
                                   InterpolationPolicy )
{
  initialize<InterpolationPolicy>( distribution );

  if( boost::is_same<InterpolationPolicy,Utility::LinLin>::value )
  {
    TEST_ASSERT( distribution->isCompatibleWithInterpType<Utility::LinLin>() );
  }
  else
  {  
    TEST_ASSERT( !distribution->isCompatibleWithInterpType<Utility::LinLin>() );
  }

  if( boost::is_same<InterpolationPolicy,Utility::LinLog>::value )
  {
    TEST_ASSERT( distribution->isCompatibleWithInterpType<Utility::LinLog>() );
  }
  else
  {
    TEST_ASSERT( !distribution->isCompatibleWithInterpType<Utility::LinLog>() );
  }

  if( boost::is_same<InterpolationPolicy,Utility::LogLin>::value )
  {
    TEST_ASSERT( distribution->isCompatibleWithInterpType<Utility::LogLin>() );
  }
  else
  {
    TEST_ASSERT( !distribution->isCompatibleWithInterpType<Utility::LogLin>() );
  }

  if( boost::is_same<InterpolationPolicy,Utility::LogLog>::value )
  {
    TEST_ASSERT( distribution->isCompatibleWithInterpType<Utility::LogLog>() );
  }
  else
  {
    TEST_ASSERT( !distribution->isCompatibleWithInterpType<Utility::LogLog>() );
  }
}

UNIT_TEST_INSTANTIATION( AnalogElasticDistribution, isCompatibleWithInterpType );

//---------------------------------------------------------------------------//
// Check if the unit-aware distribution is compatible with the interp type
TEUCHOS_UNIT_TEST_TEMPLATE_1_DECL( UnitAwareAnalogElasticDistribution,
                                   isCompatibleWithInterpType,
                                   InterpolationPolicy )
{
  initialize<InterpolationPolicy>( unit_aware_distribution );

  if( boost::is_same<InterpolationPolicy,Utility::LinLin>::value )
  {
    TEST_ASSERT( unit_aware_distribution->isCompatibleWithInterpType<Utility::LinLin>() );
  }
  else
  {  
    TEST_ASSERT( !unit_aware_distribution->isCompatibleWithInterpType<Utility::LinLin>() );
  }

  if( boost::is_same<InterpolationPolicy,Utility::LinLog>::value )
  {
    TEST_ASSERT( unit_aware_distribution->isCompatibleWithInterpType<Utility::LinLog>() );
  }
  else
  {
    TEST_ASSERT( !unit_aware_distribution->isCompatibleWithInterpType<Utility::LinLog>() );
  }

  if( boost::is_same<InterpolationPolicy,Utility::LogLin>::value )
  {
    TEST_ASSERT( unit_aware_distribution->isCompatibleWithInterpType<Utility::LogLin>() );
  }
  else
  {
    TEST_ASSERT( !unit_aware_distribution->isCompatibleWithInterpType<Utility::LogLin>() );
  }

  if( boost::is_same<InterpolationPolicy,Utility::LogLog>::value )
  {
    TEST_ASSERT( unit_aware_distribution->isCompatibleWithInterpType<Utility::LogLog>() );
  }
  else
  {
    TEST_ASSERT( !unit_aware_distribution->isCompatibleWithInterpType<Utility::LogLog>() );
  }
}

UNIT_TEST_INSTANTIATION( UnitAwareAnalogElasticDistribution,
                         isCompatibleWithInterpType );

//---------------------------------------------------------------------------//
// Check that the distribution can be written to an xml file
TEUCHOS_UNIT_TEST_TEMPLATE_1_DECL( AnalogElasticDistribution,
                                   toParameterList,
                                   InterpolationPolicy )
{
  initialize<InterpolationPolicy>( distribution );

  typedef Utility::AnalogElasticDistribution<InterpolationPolicy> Distribution;

  Teuchos::RCP<Distribution> true_distribution =
    Teuchos::rcp_dynamic_cast<Distribution>( distribution );

  Teuchos::ParameterList parameter_list;

  parameter_list.set<Distribution>( "test distribution",
                                      *true_distribution );

  std::ostringstream xml_file_name;
  xml_file_name << "analog_elastic_" << InterpolationPolicy::name()
                << "_dist_test_list.xml";

  Teuchos::writeParameterListToXmlFile( parameter_list,
                                        xml_file_name.str() );

  Teuchos::RCP<Teuchos::ParameterList> read_parameter_list =
    Teuchos::getParametersFromXmlFile( xml_file_name.str() );

  TEST_EQUALITY( parameter_list, *read_parameter_list );

  Teuchos::RCP<Distribution>
    copy_distribution( new Distribution );

  *copy_distribution = read_parameter_list->get<Distribution>(
                                                          "test distribution");

  TEST_EQUALITY( *copy_distribution, *true_distribution );
}

UNIT_TEST_INSTANTIATION( AnalogElasticDistribution, toParameterList );

//---------------------------------------------------------------------------//
// Check that the unit-aware distribution can be written to an xml file
TEUCHOS_UNIT_TEST_TEMPLATE_1_DECL( UnitAwareAnalogElasticDistribution,
                                   toParameterList,
                                   InterpolationPolicy )
{
  initialize<InterpolationPolicy>( unit_aware_distribution );

  typedef Utility::UnitAwareAnalogElasticDistribution<InterpolationPolicy,si::dimensionless,si::amount> Distribution;

  Teuchos::RCP<Distribution> true_distribution =
    Teuchos::rcp_dynamic_cast<Distribution>( unit_aware_distribution );

  Teuchos::ParameterList parameter_list;

  parameter_list.set<Distribution>( "test distribution",
                                    *true_distribution );

  std::ostringstream xml_file_name;
  xml_file_name << "unit_aware_tabular_" << InterpolationPolicy::name()
                << "_dist_test_list.xml";

  Teuchos::writeParameterListToXmlFile( parameter_list,
                                        xml_file_name.str() );

  Teuchos::RCP<Teuchos::ParameterList> read_parameter_list =
    Teuchos::getParametersFromXmlFile( xml_file_name.str() );

  TEST_EQUALITY( parameter_list, *read_parameter_list );

  Teuchos::RCP<Distribution>
    copy_distribution( new Distribution );

  *copy_distribution = read_parameter_list->get<Distribution>(
                                                          "test distribution");

  TEST_EQUALITY( *copy_distribution, *true_distribution );
}

UNIT_TEST_INSTANTIATION( UnitAwareAnalogElasticDistribution, toParameterList );

//---------------------------------------------------------------------------//
// Check that the distribution can be read from an xml file
TEUCHOS_UNIT_TEST( AnalogElasticDistribution, fromParameterList )
{
  Utility::AnalogElasticDistribution<Utility::LinLin> distribution_1 =
    test_dists_list->get<Utility::AnalogElasticDistribution<Utility::LinLin> >( "Analog Elastic Distribution A" );

  TEST_EQUALITY_CONST( distribution_1.getLowerBoundOfIndepVar(), 0.001 );
  TEST_EQUALITY_CONST( distribution_1.getUpperBoundOfIndepVar(),
                       Utility::PhysicalConstants::pi );

  distribution_1 =
    test_dists_list->get<Utility::AnalogElasticDistribution<Utility::LinLin> >( "Analog Elastic Distribution B" );

  TEST_EQUALITY_CONST( distribution_1.getLowerBoundOfIndepVar(), 0.001 );
  TEST_EQUALITY_CONST( distribution_1.getUpperBoundOfIndepVar(), 1.0 );

  Utility::AnalogElasticDistribution<Utility::LogLog> distribution_2 =
    test_dists_list->get<Utility::AnalogElasticDistribution<Utility::LogLog> >( "Analog Elastic Distribution C" );

  TEST_EQUALITY_CONST( distribution_2.getLowerBoundOfIndepVar(), 0.001 );
  TEST_EQUALITY_CONST( distribution_2.getUpperBoundOfIndepVar(), 10.0 );
}

//---------------------------------------------------------------------------//
// Check that the unit-aware distribution can be read from an xml file
TEUCHOS_UNIT_TEST( UnitAwareAnalogElasticDistribution, fromParameterList )
{
  Utility::UnitAwareAnalogElasticDistribution<Utility::LinLin,si::dimensionless,si::amount>
    distribution_1 =
    test_dists_list->get<Utility::UnitAwareAnalogElasticDistribution<Utility::LinLin,si::dimensionless,si::amount> >( "Unit-Aware Analog Elastic Distribution A" );

  TEST_EQUALITY_CONST( distribution_1.getLowerBoundOfIndepVar(), 0.001 );
  TEST_EQUALITY_CONST( distribution_1.getUpperBoundOfIndepVar(),
                       Utility::PhysicalConstants::pi );

  distribution_1 =
    test_dists_list->get<Utility::UnitAwareAnalogElasticDistribution<Utility::LinLin,si::dimensionless,si::amount> >( "Unit-Aware Analog Elastic Distribution B" );

  TEST_EQUALITY_CONST( distribution_1.getLowerBoundOfIndepVar(), 0.001 );
  TEST_EQUALITY_CONST( distribution_1.getUpperBoundOfIndepVar(), 1.0 );

  Utility::UnitAwareAnalogElasticDistribution<Utility::LogLog,si::dimensionless,si::amount>
    distribution_2 =
    test_dists_list->get<Utility::UnitAwareAnalogElasticDistribution<Utility::LogLog,si::dimensionless,si::amount> >( "Unit-Aware Analog Elastic Distribution C" );

  TEST_EQUALITY_CONST( distribution_2.getLowerBoundOfIndepVar(), 0.001 );
  TEST_EQUALITY_CONST( distribution_2.getUpperBoundOfIndepVar(), 10.0 );
}

//---------------------------------------------------------------------------//
// Check that distributions can be scaled
TEUCHOS_UNIT_TEST_TEMPLATE_4_DECL( UnitAwareAnalogElasticDistribution,
                                   explicit_conversion,
                                   IndepUnitA,
                                   DepUnitA,
                                   IndepUnitB,
                                   DepUnitB )
{
  typedef typename Utility::UnitTraits<IndepUnitA>::template GetQuantityType<double>::type IndepQuantityA;
  typedef typename Utility::UnitTraits<typename Utility::UnitTraits<IndepUnitA>::InverseUnit>::template GetQuantityType<double>::type InverseIndepQuantityA;

  typedef typename Utility::UnitTraits<IndepUnitB>::template GetQuantityType<double>::type IndepQuantityB;
  typedef typename Utility::UnitTraits<typename Utility::UnitTraits<IndepUnitB>::InverseUnit>::template GetQuantityType<double>::type InverseIndepQuantityB;

  typedef typename Utility::UnitTraits<DepUnitA>::template GetQuantityType<double>::type DepQuantityA;
  typedef typename Utility::UnitTraits<DepUnitB>::template GetQuantityType<double>::type DepQuantityB;

  initialize<Utility::LinLin>( distribution );

  // Copy from unitless distribution to distribution type A
  Utility::UnitAwareAnalogElasticDistribution<Utility::LinLin,IndepUnitA,DepUnitA>
    unit_aware_dist_a_copy = Utility::UnitAwareAnalogElasticDistribution<Utility::LinLin,IndepUnitA,DepUnitA>::fromUnitlessDistribution( *Teuchos::rcp_dynamic_cast<Utility::AnalogElasticDistribution<Utility::LinLin> >( distribution ) );

  // Copy from distribution type A to distribution type B
  Utility::UnitAwareAnalogElasticDistribution<Utility::LinLin,IndepUnitB,DepUnitB>
    unit_aware_dist_b_copy( unit_aware_dist_a_copy );

  IndepQuantityA indep_quantity_a =
    Utility::QuantityTraits<IndepQuantityA>::initializeQuantity( 0.0 );
  InverseIndepQuantityA inv_indep_quantity_a =
    Utility::QuantityTraits<InverseIndepQuantityA>::initializeQuantity( 0.0 );
  DepQuantityA dep_quantity_a =
    Utility::QuantityTraits<DepQuantityA>::initializeQuantity( 0.0 );

  IndepQuantityB indep_quantity_b( indep_quantity_a );
  InverseIndepQuantityB inv_indep_quantity_b( inv_indep_quantity_a );
  DepQuantityB dep_quantity_b( dep_quantity_a );

  UTILITY_TEST_FLOATING_EQUALITY(
                           unit_aware_dist_a_copy.evaluate( indep_quantity_a ),
                           dep_quantity_a,
                           1e-15 );
  UTILITY_TEST_FLOATING_EQUALITY(
                        unit_aware_dist_a_copy.evaluatePDF( indep_quantity_a ),
                        inv_indep_quantity_a,
                        1e-15 );
  UTILITY_TEST_FLOATING_EQUALITY(
                           unit_aware_dist_b_copy.evaluate( indep_quantity_b ),
                           dep_quantity_b,
                           1e-15 );
  UTILITY_TEST_FLOATING_EQUALITY(
                        unit_aware_dist_b_copy.evaluatePDF( indep_quantity_b ),
                        inv_indep_quantity_b,
                        1e-15 );

  Utility::setQuantity( indep_quantity_a, 0.1 );
  Utility::setQuantity( inv_indep_quantity_a, 0.67340006734 );
  Utility::setQuantity( dep_quantity_a, 1.0 );

  indep_quantity_b = IndepQuantityB( indep_quantity_a );
  inv_indep_quantity_b = InverseIndepQuantityB( inv_indep_quantity_a );
  dep_quantity_b = DepQuantityB( dep_quantity_a );

  UTILITY_TEST_FLOATING_EQUALITY(
                           unit_aware_dist_a_copy.evaluate( indep_quantity_a ),
                           dep_quantity_a,
                           1e-15 );
  UTILITY_TEST_FLOATING_EQUALITY(
                        unit_aware_dist_a_copy.evaluatePDF( indep_quantity_a ),
                        inv_indep_quantity_a,
                        1e-6 );
  UTILITY_TEST_FLOATING_EQUALITY(
                           unit_aware_dist_b_copy.evaluate( indep_quantity_b ),
                           dep_quantity_b,
                           1e-15 );
  UTILITY_TEST_FLOATING_EQUALITY(
                        unit_aware_dist_b_copy.evaluatePDF( indep_quantity_b ),
                        inv_indep_quantity_b,
                        1e-6 );
}

typedef si::energy si_energy;
typedef cgs::energy cgs_energy;
typedef si::amount si_amount;
typedef si::length si_length;
typedef cgs::length cgs_length;
typedef si::mass si_mass;
typedef cgs::mass cgs_mass;
typedef si::dimensionless si_dimensionless;
typedef cgs::dimensionless cgs_dimensionless;

TEUCHOS_UNIT_TEST_TEMPLATE_4_INSTANT( UnitAwareAnalogElasticDistribution,
                                      explicit_conversion,
                                      si_dimensionless,
                                      si_amount,
                                      cgs_dimensionless,
                                      si_amount );
TEUCHOS_UNIT_TEST_TEMPLATE_4_INSTANT( UnitAwareAnalogElasticDistribution,
                                      explicit_conversion,
                                      cgs_dimensionless,
                                      si_amount,
                                      si_dimensionless,
                                      si_amount );
TEUCHOS_UNIT_TEST_TEMPLATE_4_INSTANT( UnitAwareAnalogElasticDistribution,
                                      explicit_conversion,
                                      si_dimensionless,
                                      si_length,
                                      cgs_dimensionless,
                                      cgs_length );
TEUCHOS_UNIT_TEST_TEMPLATE_4_INSTANT( UnitAwareAnalogElasticDistribution,
                                      explicit_conversion,
                                      cgs_dimensionless,
                                      cgs_length,
                                      si_dimensionless,
                                      si_length );
TEUCHOS_UNIT_TEST_TEMPLATE_4_INSTANT( UnitAwareAnalogElasticDistribution,
                                      explicit_conversion,
                                      si_dimensionless,
                                      si_mass,
                                      cgs_dimensionless,
                                      cgs_mass );
TEUCHOS_UNIT_TEST_TEMPLATE_4_INSTANT( UnitAwareAnalogElasticDistribution,
                                      explicit_conversion,
                                      cgs_dimensionless,
                                      cgs_mass,
                                      si_dimensionless,
                                      si_mass );
TEUCHOS_UNIT_TEST_TEMPLATE_4_INSTANT( UnitAwareAnalogElasticDistribution,
                                      explicit_conversion,
                                      si_dimensionless,
                                      si_dimensionless,
                                      cgs_dimensionless,
                                      cgs_dimensionless );
TEUCHOS_UNIT_TEST_TEMPLATE_4_INSTANT( UnitAwareAnalogElasticDistribution,
                                      explicit_conversion,
                                      cgs_dimensionless,
                                      cgs_dimensionless,
                                      si_dimensionless,
                                      si_dimensionless );
TEUCHOS_UNIT_TEST_TEMPLATE_4_INSTANT( UnitAwareAnalogElasticDistribution,
                                      explicit_conversion,
                                      si_dimensionless,
                                      void,
                                      cgs_dimensionless,
                                      void );
TEUCHOS_UNIT_TEST_TEMPLATE_4_INSTANT( UnitAwareAnalogElasticDistribution,
                                      explicit_conversion,
                                      cgs_dimensionless,
                                      void,
                                      si_dimensionless,
                                      void );
//TEUCHOS_UNIT_TEST_TEMPLATE_4_INSTANT( UnitAwareAnalogElasticDistribution,
//                                      explicit_conversion,
//                                      ElectronVolt,
//                                      si_amount,
//                                      si_energy,
//                                      si_amount );
//TEUCHOS_UNIT_TEST_TEMPLATE_4_INSTANT( UnitAwareAnalogElasticDistribution,
//                                      explicit_conversion,
//                                      ElectronVolt,
//                                      si_amount,
//                                      cgs_energy,
//                                      si_amount );
//TEUCHOS_UNIT_TEST_TEMPLATE_4_INSTANT( UnitAwareAnalogElasticDistribution,
//                                      explicit_conversion,
//                                      ElectronVolt,
//                                      si_amount,
//                                      KiloElectronVolt,
//                                      si_amount );
//TEUCHOS_UNIT_TEST_TEMPLATE_4_INSTANT( UnitAwareAnalogElasticDistribution,
//                                      explicit_conversion,
//                                      ElectronVolt,
//                                      si_amount,
//                                      MegaElectronVolt,
//                                      si_amount );
//TEUCHOS_UNIT_TEST_TEMPLATE_4_INSTANT( UnitAwareAnalogElasticDistribution,
//                                      explicit_conversion,
//                                      KiloElectronVolt,
//                                      si_amount,
//                                      si_energy,
//                                      si_amount );
//TEUCHOS_UNIT_TEST_TEMPLATE_4_INSTANT( UnitAwareAnalogElasticDistribution,
//                                      explicit_conversion,
//                                      KiloElectronVolt,
//                                      si_amount,
//                                      cgs_energy,
//                                      si_amount );
//TEUCHOS_UNIT_TEST_TEMPLATE_4_INSTANT( UnitAwareAnalogElasticDistribution,
//                                      explicit_conversion,
//                                      KiloElectronVolt,
//                                      si_amount,
//                                      ElectronVolt,
//                                      si_amount );
//TEUCHOS_UNIT_TEST_TEMPLATE_4_INSTANT( UnitAwareAnalogElasticDistribution,
//                                      explicit_conversion,
//                                      KiloElectronVolt,
//                                      si_amount,
//                                      MegaElectronVolt,
//                                      si_amount );
//TEUCHOS_UNIT_TEST_TEMPLATE_4_INSTANT( UnitAwareAnalogElasticDistribution,
//                                      explicit_conversion,
//                                      MegaElectronVolt,
//                                      si_amount,
//                                      si_energy,
//                                      si_amount );
//TEUCHOS_UNIT_TEST_TEMPLATE_4_INSTANT( UnitAwareAnalogElasticDistribution,
//                                      explicit_conversion,
//                                      MegaElectronVolt,
//                                      si_amount,
//                                      cgs_energy,
//                                      si_amount );
//TEUCHOS_UNIT_TEST_TEMPLATE_4_INSTANT( UnitAwareAnalogElasticDistribution,
//                                      explicit_conversion,
//                                      MegaElectronVolt,
//                                      si_amount,
//                                      ElectronVolt,
//                                      si_amount );
//TEUCHOS_UNIT_TEST_TEMPLATE_4_INSTANT( UnitAwareAnalogElasticDistribution,
//                                      explicit_conversion,
//                                      MegaElectronVolt,
//                                      si_amount,
//                                      KiloElectronVolt,
//                                      si_amount );
//TEUCHOS_UNIT_TEST_TEMPLATE_4_INSTANT( UnitAwareAnalogElasticDistribution,
//                                      explicit_conversion,
//                                      void,
//                                      MegaElectronVolt,
//                                      void,
//                                      KiloElectronVolt );

//---------------------------------------------------------------------------//
// Custom setup
//---------------------------------------------------------------------------//
UTILITY_CUSTOM_TEUCHOS_UNIT_TEST_SETUP_BEGIN();

std::string test_dists_xml_file;

UTILITY_CUSTOM_TEUCHOS_UNIT_TEST_COMMAND_LINE_OPTIONS()
{
  clp().setOption( "test_dists_xml_file",
                   &test_dists_xml_file,
                   "Test distributions xml file name" );
}

UTILITY_CUSTOM_TEUCHOS_UNIT_TEST_DATA_INITIALIZATION()
{
  TEUCHOS_ADD_TYPE_CONVERTER( Utility::AnalogElasticDistribution<Utility::LinLin> );
  TEUCHOS_ADD_TYPE_CONVERTER( Utility::AnalogElasticDistribution<Utility::LogLin> );
  TEUCHOS_ADD_TYPE_CONVERTER( Utility::AnalogElasticDistribution<Utility::LinLog> );
  TEUCHOS_ADD_TYPE_CONVERTER( Utility::AnalogElasticDistribution<Utility::LogLog> );
  typedef Utility::UnitAwareAnalogElasticDistribution<Utility::LinLin,si::dimensionless,si::amount> UnitAwareAnalogElasticLinLinDist;
  TEUCHOS_ADD_TYPE_CONVERTER( UnitAwareAnalogElasticLinLinDist );
  typedef Utility::UnitAwareAnalogElasticDistribution<Utility::LogLin,si::dimensionless,si::amount> UnitAwareAnalogElasticLogLinDist;
  TEUCHOS_ADD_TYPE_CONVERTER( UnitAwareAnalogElasticLogLinDist );
  typedef Utility::UnitAwareAnalogElasticDistribution<Utility::LinLog,si::dimensionless,si::amount> UnitAwareAnalogElasticLinLogDist;
  TEUCHOS_ADD_TYPE_CONVERTER( UnitAwareAnalogElasticLinLogDist );
  typedef Utility::UnitAwareAnalogElasticDistribution<Utility::LogLog,si::dimensionless,si::amount> UnitAwareAnalogElasticLogLogDist;
  TEUCHOS_ADD_TYPE_CONVERTER( UnitAwareAnalogElasticLogLogDist );

  test_dists_list = Teuchos::getParametersFromXmlFile( test_dists_xml_file );

  // Initialize the random number generator
  Utility::RandomNumberGenerator::createStreams();
}

UTILITY_CUSTOM_TEUCHOS_UNIT_TEST_SETUP_END();

//---------------------------------------------------------------------------//
// end tstAnalogElasticDistribution.cpp
//---------------------------------------------------------------------------//