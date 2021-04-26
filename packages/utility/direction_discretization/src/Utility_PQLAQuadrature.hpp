//---------------------------------------------------------------------------//
//!
//! \file   Utility_PQLAQuadrature.hpp
//! \author Philip Britt
//! \brief  PQLA Direction Quadrature handler declaration
//!
//---------------------------------------------------------------------------//

#ifndef UTILITY_PQLA_QUADRATURE
#define UTILITY_PQLA_QUADRATURE

// Boost Includes
#include <boost/serialization/version.hpp>
#include <boost/serialization/export.hpp>
#include <boost/serialization/shared_ptr.hpp>

// FRENSIE includes
#include "Utility_Vector.hpp"
#include "Utility_Tuple.hpp"
#include "Utility_Array.hpp"
#include "Utility_ExplicitSerializationTemplateInstantiationMacros.hpp"
#include "Utility_SerializationHelpers.hpp"
#include "Utility_3DCartesianVectorHelpers.hpp"

namespace Utility{

struct SphericalTriangle
{
  /*! Vector that contains a tuple representing spherical triangle parameters.
   * \details First element of the tuple is an array that contains the
   * 2-norm direction representing a vertex of the triangle. 
   * Second is the length of the spherical triangle side opposite from that vertex
   * (or angle that the 2 other vertices of the triangle make with each other).
   * Third is the angle made from the sides of the spherical triangle from that vertex.
   * Note their order does NOT matter so long as they are consistent with the above definition
   */
  std::vector<std::tuple<std::array<double, 3>, double, double>> triangle_parameter_vector;

  //! Area of the triangle
  double area;

  void computeAndStoreTriangleParameters(std::vector<std::array<double, 3>>& vertex_vector)
  {
    // Put methods in struct to simplify this part.
    
    // calculate cosine of length of side of spherical triangle opposite from respective vertex (for use later, not kept as member data)
    std::vector<double> opposite_cos_vector {calculateCosineOfAngleBetweenUnitVectors(vertex_vector[1].data(), vertex_vector[2].data()),
                                              calculateCosineOfAngleBetweenUnitVectors(vertex_vector[0].data(), vertex_vector[2].data()),
                                              calculateCosineOfAngleBetweenUnitVectors(vertex_vector[0].data(), vertex_vector[1].data())};

    // calculate length of side of spherical triangle opposite from respective vertex (in radians b/c unit sphere)
    std::vector<double> opposite_side_length_vector{acos(opposite_cos_vector[0]), acos(opposite_cos_vector[1]), acos(opposite_cos_vector[2])};

    std::vector<double> angle_vector{acos((opposite_cos_vector[0] - opposite_cos_vector[1]*opposite_cos_vector[2])/(sin(opposite_side_length_vector[1])*sin(opposite_side_length_vector[2]))),
                                    acos((opposite_cos_vector[1] - opposite_cos_vector[0]*opposite_cos_vector[2])/(sin(opposite_side_length_vector[0])*sin(opposite_side_length_vector[2]))),
                                    acos((opposite_cos_vector[2] - opposite_cos_vector[0]*opposite_cos_vector[1])/(sin(opposite_side_length_vector[0])*sin(opposite_side_length_vector[1])))};

    for(size_t vert = 0; vert < 3; ++vert)
    {
      triangle_parameter_vector.push_back( std::make_tuple(vertex_vector[vert], opposite_side_length_vector[vert], angle_vector[vert]));
    }

    // Store triangle area
    area = angle_vector[0]+angle_vector[1]+angle_vector[2]-M_PI;
  }

  // Serialize the data
  template<typename Archive>
  void serialize( Archive& ar, const unsigned version )
  { 
    ar & BOOST_SERIALIZATION_NVP( triangle_parameter_vector );
    ar & BOOST_SERIALIZATION_NVP( area );
  }
};

class PQLAQuadrature
{

  public:

  //! Constructor
  PQLAQuadrature(unsigned quadrature_order);

  //! Destructor
  ~PQLAQuadrature()
  { /* ... */ }

  //! Find whether or not a triangle id is valid
  bool isTriangleIDValid( const size_t triangle_id) const;

  //! Find which triangle bin a direction vector is in
  size_t findTriangleBin( const std::array<double, 3>& direction) const;

  //! Find which triangle bin a direction vector is in
  size_t findTriangleBin( const double x_direction, const double y_direction, const double z_direction) const;

  //! Return the order of the quadrature
  unsigned getQuadratureOrder() const;

  //! Get the total number of triangles
  size_t getNumberOfTriangles() const;
  
  //! Get the area of a specific spherical triangle
  double getTriangleArea(const size_t triangle_index) const;

  /*! Get a random direction from within a spherical triangle (evenly distributed probability) - reference here
   * \details reference: Stratified Sampling of Spherical Triangles, James Arvo, SIGGRAPH '95
   */
  void sampleIsotropicallyFromTriangle(std::array<double, 3>& direction_vector, 
                                       const size_t triangle_index) const;

  //! Used for archive testing
  const std::vector<SphericalTriangle>& getSphericalTriangleVector() const;

  std::array<double, 3> getTriangleCentroid(const size_t triangle_index) const;

  private:

  //! Default constructor (for archiving)
  PQLAQuadrature()
  { /* ... */ }

  //! Vector operation for the purpose of sampleIsotropicallyFromTriangle
  void isotropicSamplingVectorOperation(const std::array<double, 3>& vertex_1,
                                        const std::array<double, 3>& vertex_2,
                                        std::array<double, 3>& result_vector) const;

  //! Converts direction vector to 1-norm normalized vector
  void normalizeVectorToOneNorm( const std::array<double, 3>& direction_2_norm,
                                 std::array<double, 3>& direction_1_norm) const;
  
  //! Converts direction vector to 1-norm normalized vector
  void normalizeVectorToOneNorm(  const double x_direction, 
                                  const double y_direction, 
                                  const double z_direction,
                                  std::array<double, 3>& direction_1_norm) const;

  //! Take lower bounding plane indices of direction vector to form triangle index
  size_t calculatePositiveTriangleBinIndex(const unsigned i_x, const unsigned i_y, const unsigned i_z) const;

  //! Take direction signs to calculate secondary index
  size_t findSecondaryIndex(const bool x_sign, const bool y_sign, const bool z_sign) const;

  //! Get a specific reference to a spherical triangle
  void getSphericalTriangle(const size_t triangle_index,
                            SphericalTriangle& triangle) const ;

  //! Quadrature order
  unsigned d_quadrature_order;

  //! Vector that stores POSITIVE DOMAIN spherical triangles
  std::vector<SphericalTriangle> d_spherical_triangle_vector;

  //! Serialize the data
  template<typename Archive>
  void serialize( Archive& ar, const unsigned version );

  //! Declare the boost serialization access object as a friend
  friend class boost::serialization::access;

};

// Serialize the data
template<typename Archive>
void PQLAQuadrature::serialize( Archive& ar, const unsigned version )
{
  // Serialize the member data
  ar & BOOST_SERIALIZATION_NVP( d_quadrature_order );
  ar & BOOST_SERIALIZATION_NVP( d_spherical_triangle_vector );
}

} // end Utility namespace

BOOST_SERIALIZATION_CLASS_VERSION( PQLAQuadrature, Utility, 0 );
BOOST_SERIALIZATION_CLASS_EXPORT_STANDARD_KEY( PQLAQuadrature, Utility );
EXTERN_EXPLICIT_CLASS_SERIALIZE_INST( Utility, PQLAQuadrature );

#endif // end UTILITY_PQLA_QUADRATURE

//---------------------------------------------------------------------------//
// end Utility_PQLADiscetization.hpp
//---------------------------------------------------------------------------//