//---------------------------------------------------------------------------//
//!
//! \file   CellFactory.hpp
//! \author Alex Robinson
//! \brief  Cell factory class declaration
//!
//---------------------------------------------------------------------------//

#ifndef CELL_FACTORY_HPP
#define CELL_FACTORY_HPP

// Std Lib Includes
#include <string>
#include <map>

// Trilinos Includes
#include <Teuchos_ArrayView.hpp>
#include <Teuchos_RCP.hpp>

// FACEMC Includes
#include "Cell.hpp"
#include "Surface.hpp"
#include "Polygon.hpp"
#include "Tuple.hpp"

namespace FACEMC{

/*! \brief Enumeration to specify the cell attribute(s) to calc. numerically
 */
enum NumericalCellIntegrationType{
  CELL_VOLUME,
  CELL_SURFACE_AREAS,
  CELL_VOLUME_AND_SURFACE_AREAS
};  

/*! A factory class for creating smart pointers to cell objects and possibly
 * calculating the volume and surface areas.
 */
template<typename Cell, typename SurfaceMap>
class CellFactory 
{
public:

  //@{
  //! Typedefs
  //! Typedef for cell pointer
  typedef typename Teuchos::RCP<Cell> CellPtr;
  //! Typedef for bounding box
  typedef Quad<Pair<typename Cell::ordinalType, NumericalCellIntegrationType>,
	       Pair<typename Cell::scalarType,typename Cell::scalarType>,
	       Pair<typename Cell::scalarType,typename Cell::scalarType>,
	       Pair<typename Cell::scalarType,typename Cell::scalarType> >
	       BoudingBox;
  //@}

private:
  
  //! Typedef for cell ordinal type
  typedef typename Cell::ordinalType cellOrdinalType;
  //! Typedef for surface ordinal type
  typedef typename Cell::surfaceOrdinalType surfaceOrdinalType;
  //! Typedef for scalar type
  typedef typename Cell::scalarType scalarType;
  //! Typedef for suface-sense pairs iterator
  typedef typename Cell::SurfaceSensePairsIterator SurfaceSensePairsIterator;
  //! Typedef for Boolean array
  typedef typename Cell::BooleanArray BooleanArray;
  //! Typedef for scalar traits
  typedef Teuchos::ScalarTraits<scalarType> ST;
  //! Typedef for cell ordinal traits
  typedef Teuchos::OrdinalTraits<cellOrdinalType> CellOT;
  //! Typedef for surface ordinal traits
  typedef Teuchos::OrdinalTraits<surfaceOrdinalType> SurfaceOT;
  //! Typedef for intersection point
  typedef Pair<Trip<scalarType,scalarType,scalarType>,
	       Trip<surfaceOrdinalType,surfaceOrdinalType,surfaceOrdinalType> >
  IntersectionPoint;
  //!Typedef for polygon
  typedef Polygon<surfaceOrdinalType,scalarType> Polygon;

public:

  //! Constructor
  CellFactory( const SurfaceMap &global_surface_map );

  //! Destructor
  ~CellFactory()
  { /* ... */ }

  //! Create the desired cell
  CellPtr create( const typename Cell::ordinalType id,
		  std::string &cell_definition,
		  const bool calculate_volume_and_area = false ) const;

  //! Get the cells needing Monte Carlo integration
  Teuchos::Array<typename Cell::ordinalType> 
  getCellsNeedingMonteCarloIntegration() const;

protected:

  //! Calculate the volume and area of a polyhedral cell
  void calculatePolyhedralCellVolumeAndArea( CellPtr &cell );

  //! Create the polygons bounding the cell
  static void createBoundingPolygons( 
				  const Teuchos::Array<Polygon> &cell_polygons,
				  CellPtr &cell );

  //! Calculate the volume of a polyhedral cell using bounding polygons
  static scalarType calculatePolyhedralCellVolumeFromPolygons( 
				const Teuchos::Array<Polygon> &cell_polygons );

  //! Create a bounding box for a polyhedral cell from bounding polygons
  static void createBoundingBoxFromPolygons( 
				BoundingBox &cell_bounding_box,
				const Teuchos::Array<Polygon> &cell_polygons );

  //! Calculate the intersection points of planes with a plane of interest
  static void calculateIntersectionPointsOnPlane(
			    std::list<IntersectionPoint> &intersection_points,
			    const SurfaceSensePairsIterator &plane_of_polygon,
			    const CellPtr &cell ) const;

  //! Calculate the volume and area of a rotationally symmetric cell
  void calculateRotationallySymmetricCellVolumeAndArea( CellPtr &cell ); 

private:

  // Stored copy of the global surface map (weak pointer)
  Teuchos::RCP<const SurfaceMap> d_global_surface_map;

  // Bounding boxes of cells needing MC integration to determine vol (and area)
  Teuchos::Array<BoundingBox> d_bounding_boxes;
};

} // end FACEMC namespace

//---------------------------------------------------------------------------//
// Template includes.
//---------------------------------------------------------------------------//

#include CellFactory_def.hpp

//---------------------------------------------------------------------------//

#endif // end CELL_FACTORY_HPP

//---------------------------------------------------------------------------//
// end CellFactory.hpp
//---------------------------------------------------------------------------//

