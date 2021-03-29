//---------------------------------------------------------------------------//
//!
//! \file   tstPQLAQuadrature.cpp
//! \author Philip Britt
//! \brief  PQLA Quadrature class unit tests
//!
//---------------------------------------------------------------------------//

// Std Lib Includes
#include <memory>
#include <array>
#include <utility>
#include <cmath>

// FRENSIE Includes
#include "Utility_3DCartesianVectorHelpers.hpp"
#include "Utility_PQLAQuadrature.hpp"
#include "Utility_UnitTestHarnessWithMain.hpp"
#include "FRENSIE_config.hpp"
#include "ArchiveTestHelpers.hpp"
#include "Utility_RandomNumberGenerator.hpp"

// Add unit test here - one octant, true random number stream (n=1) check that w is distributed uniformly on poles, and arctan(v/u) distributed uniformly
// Next test is testing if it samples from each triangle uniformly

//---------------------------------------------------------------------------//
// Test n = 1 quadrature order
//---------------------------------------------------------------------------//
/*
FRENSIE_UNIT_TEST(PQLAQuadrature, sampleIsotropically_n1)
{

  size_t quadrature_order = 1;
  std::shared_ptr<Utility::PQLAQuadrature> direction_discretization = std::make_shared<Utility::PQLAQuadrature>(quadrature_order);
  std::vector<double> z_bin_bounds;
  std::vector<double> theta_bin_bounds;
  // Initialize the random number generator
  Utility::RandomNumberGenerator::createStreams();
  // 1000 bins in each dimension
  size_t number_of_bins = 1000;
  for(unsigned i = 0; i <= number_of_bins; ++i)
  {
    double frac = static_cast<double>(i)/number_of_bins;
    z_bin_bounds.push_back(frac);
    theta_bin_bounds.push_back((M_PI/2)*frac);
  }
  std::vector<double> z_results;
  z_results.resize(number_of_bins);
  std::vector<double> theta_results;
  theta_results.resize(number_of_bins);
  // Take 10 million samples from triangle 0 and bin
  size_t number_of_samples = 300000000;
  for(unsigned i = 0; i < number_of_samples; ++i)
  {
    std::array<double, 3> direction;
    direction_discretization->sampleIsotropicallyFromTriangle(direction, 0);

    double z_value = direction[2];
    double theta_value = atan(direction[2]/direction[1]);
    for( size_t z_bin = 0; z_bin < z_bin_bounds.size()-1; ++z_bin )
    {
      if( z_bin_bounds[z_bin] <= z_value && z_value < z_bin_bounds[z_bin+1] )
      {
        z_results[z_bin] += 1.0;
        break;
      }
    }

    for( size_t theta_bin = 0; theta_bin < theta_bin_bounds.size()-1; ++theta_bin)
    {
      if( theta_bin_bounds[theta_bin] <= theta_value && theta_value < theta_bin_bounds[theta_bin+1])
      {
        theta_results[theta_bin] += 1.0;
        break;
      }
    }

  }

  double expectation_value = 1.0/static_cast<double>(number_of_bins);
  double theta_sum = 0;
  double z_sum = 0;

  for( size_t bin = 0; bin < z_results.size(); ++bin)
  {
    theta_results[bin] /= number_of_samples;
    theta_sum += theta_results[bin];

    FRENSIE_CHECK_FLOATING_EQUALITY(theta_results[bin], expectation_value, 1e-2);

    z_results[bin] /= number_of_samples;
    z_sum += z_results[bin];

    FRENSIE_CHECK_FLOATING_EQUALITY(z_results[bin], expectation_value, 1e-2);
  }

  FRENSIE_CHECK_FLOATING_EQUALITY(theta_sum, 1.0, 1e-14);
  FRENSIE_CHECK_FLOATING_EQUALITY(z_sum, 1.0, 1e-14)

}
*/
//---------------------------------------------------------------------------//
// Test n = 2 quadrature order
//---------------------------------------------------------------------------//
FRENSIE_UNIT_TEST(PQLAQuadrature, sampleIsotropically_n2)
{

  size_t quadrature_order = 2;
  std::shared_ptr<Utility::PQLAQuadrature> direction_discretization = std::make_shared<Utility::PQLAQuadrature>(quadrature_order);
  std::vector<double> z_bin_bounds;
  std::vector<double> theta_bin_bounds;
  // Initialize the random number generator
  Utility::RandomNumberGenerator::createStreams();
  // 1000 bins in each dimension
  size_t number_of_bins = 1000;
  for(unsigned i = 0; i <= number_of_bins; ++i)
  {
    double frac = static_cast<double>(i)/number_of_bins;
    z_bin_bounds.push_back(2*frac-1);
    theta_bin_bounds.push_back((M_PI)*frac-(M_PI/2));
  }

  // Form CDF for sampling from separate triangles
  std::vector<double> octant_triangle_PDF;
  double area_sum = 0;
  std::vector<double> octant_triangle_CDF;
  octant_triangle_CDF.push_back(0);

  for(size_t i = 0; i < 32; ++i)
  {
    octant_triangle_PDF.push_back(direction_discretization->getTriangleArea(i)/(4*M_PI));
    area_sum += octant_triangle_PDF[i];
    octant_triangle_CDF.push_back(area_sum);
  }

  FRENSIE_CHECK_FLOATING_EQUALITY(area_sum, 1.0, 1e-14);
  FRENSIE_CHECK_FLOATING_EQUALITY(octant_triangle_CDF.back(), 1.0, 1e-14);

  std::vector<double> z_results;
  z_results.resize(number_of_bins);
  std::vector<double> theta_results;
  theta_results.resize(number_of_bins);
  // Take 10 million samples from triangle 0 and bin
  size_t number_of_samples = 200000000;
  for(unsigned i = 0; i < number_of_samples; ++i)
  {
    std::array<double, 3> direction;

    double random_number = Utility::RandomNumberGenerator::getRandomNumber<double>();
    unsigned triangle_index;
    for(size_t i = 0; i < octant_triangle_CDF.size()-1; ++i)
    {
      if(octant_triangle_CDF[i] <= random_number && random_number < octant_triangle_CDF[i+1])
      {
        triangle_index = i;
        break;
      }
    }

    direction_discretization->sampleIsotropicallyFromTriangle(direction, triangle_index);

    double z_value = direction[2];
    double theta_value = atan(direction[2]/direction[1]);
    for( size_t z_bin = 0; z_bin < z_bin_bounds.size()-1; ++z_bin )
    {
      if( z_bin_bounds[z_bin] <= z_value && z_value < z_bin_bounds[z_bin+1] )
      {
        z_results[z_bin] += 1.0;
        break;
      }
    }

    for( size_t theta_bin = 0; theta_bin < theta_bin_bounds.size()-1; ++theta_bin)
    {
      if( theta_bin_bounds[theta_bin] <= theta_value && theta_value < theta_bin_bounds[theta_bin+1])
      {
        theta_results[theta_bin] += 1.0;
        break;
      }
    }

  }

  double expectation_value = 1.0/static_cast<double>(number_of_bins);
  double theta_sum = 0;
  double z_sum = 0;

  for( size_t bin = 0; bin < z_results.size(); ++bin)
  {
    theta_results[bin] /= number_of_samples;
    theta_sum += theta_results[bin];

    FRENSIE_CHECK_FLOATING_EQUALITY(theta_results[bin], expectation_value, 1e-2);

    z_results[bin] /= number_of_samples;
    z_sum += z_results[bin];

    FRENSIE_CHECK_FLOATING_EQUALITY(z_results[bin], expectation_value, 1e-2);
  }

  FRENSIE_CHECK_FLOATING_EQUALITY(theta_sum, 1.0, 1e-14);
  FRENSIE_CHECK_FLOATING_EQUALITY(z_sum, 1.0, 1e-14)


}

//---------------------------------------------------------------------------//
// end tstPQLAQuadrature.cpp
//---------------------------------------------------------------------------//
