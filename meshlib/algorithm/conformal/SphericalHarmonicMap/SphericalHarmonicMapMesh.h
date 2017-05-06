/*!
*      \file SphericalHarmonicMap.h
*      \brief Mesh for Spherical Harmonic Map
*	   \author David Gu
*      \date Document 12/15/2010
*
*/
#ifndef  _SPHERICAL_HARMONIC_MAP_H_
#define  _SPHERICAL_HARMONIC_MAP_H_

#include <map>
#include <vector>

#include "Mesh/BaseMesh.h"
#include "Mesh/Vertex.h"
#include "Mesh/HalfEdge.h"
#include "Mesh/Edge.h"
#include "Mesh/Face.h"
#include "mesh/iterators.h"
#include "mesh/boundary.h"
#include "Parser/parser.h"
#include "Parser/traits_io.h"

namespace MeshLib
{
namespace Holomorphy
{
/*!
*	\brief CSHMVertex class
*
*	Vertex class for spherical harmonic map
*/
class CSHMVertex : public CVertex
{
public:

	/*!
	 *	CSHMVertex constructor
	 */
	CSHMVertex() {};
	/*!
	 *	CSHMVertex destructor
	 */
	~CSHMVertex() {};

	/*!	Vertex spherical harmonic map image coordinates
	 */
	CPoint & u() { return m_u; };

	/*!	Vertex normal
	 */
	CPoint & normal() { return m_normal; };

	/*!	Vertex Laplacian
	 */
	CPoint & L() { return m_L; };

	/*! write vertex traits from the vertex string */
	void  _to_string();

protected:	//output
	/*! vertex harmonic map image coordinates */
	CPoint  m_u;	
	/*! normal */
	CPoint  m_normal;
	/*! Laplacian */
	CPoint  m_L;
};


/*! write vertex traits m_u to the string */
inline	void  CSHMVertex::_to_string()
	{
		CParser parser( m_string );
		parser._removeToken( "u" );

		parser._toString( m_string );
		
		std::string line;
		std::stringstream iss(line);
		iss << "u=(" << m_u[0] << " " << m_u[1] << " " << m_u[2] << ")";

		if( m_string.length() > 0 )
		{
			m_string += " ";
		}
		m_string += iss.str();

	};


/*! \brief CSHMEdge class
*
* Edge class for computing spherical harmonioc maps
*/

class CSHMEdge : public  CEdge
  {
  public:

  /*! CSHMEdge constructor
  */
	 CSHMEdge() { m_weight = 0; };
  /*! CSHMEdge destructor
  */
    ~CSHMEdge(){};
	/*! Edge weight */
	double & weight() { return m_weight; };
	/*! Edge Length */
	double & length() { return m_length; };

  protected: //output
	/*! edge weight */
	double   m_weight;
	/*! edge length */
	double   m_length;
};

/*! \brief CSHMFace class
*
* Face class for computing spherical harmonioc maps
*/

class CSHMFace : public  CFace
  {
  public:

  /*! CSHMFace constructor
  */
	 CSHMFace() { m_area = 0; };
  /*! CSHMFace destructor
  */
    ~CSHMFace(){};
	/*! Face area   */
	double & area() { return m_area; };
	/*! Face normal */
	CPoint & normal() { return m_normal; };

  protected: 
	/* Face area */
	double   m_area;
	/* Face normal */
	CPoint   m_normal;
};

/*! \brief CSHMHalfEdge class
*
*  Halfedge class for computing spherical harmonic maps
*/
class CSHMHalfEdge : public  CHalfEdge
  {
  public:
	/*! CSHMHalfEdge constructor */
	 CSHMHalfEdge() { m_angle = 0; };
	/*! CHCFHalfEdge destructor */
    ~CSHMHalfEdge(){};
	  /*! corner angle */
	double &   angle() { return m_angle; };
	
  //no output
  protected: 
	  /*! corner angle */
	  double   m_angle;

};

/*! \brief CSphericalHarmonicMapMesh class
*
*  Mesh class for computing spherical harmonic maps 
*/
template<typename V, typename E, typename F, typename H>
class CSphericalHarmonicMapMesh : public CBaseMesh<V,E,F,H>
{
public:
	typedef V CVertex;
	typedef E CEdge;
	typedef F CFace;
	typedef H CHalfEdge;

	typedef CBoundary<V,E,F,H>  CBoundary;
	typedef CLoop<V,E,F,H> CLoop;
	typedef MeshVertexIterator<V,E,F,H> MeshVertexIterator;
	typedef MeshFaceIterator<V,E,F,H> MeshFaceIterator;
	typedef MeshEdgeIterator<V,E,F,H> MeshEdgeIterator;
	typedef FaceVertexIterator<V,E,F,H> FaceVertexIterator;
	typedef VertexVertexIterator<V,E,F,H> VertexVertexIterator;
	typedef VertexEdgeIterator<V,E,F,H> VertexEdgeIterator;
	typedef VertexFaceIterator<V,E,F,H> VertexFaceIterator;
	typedef VertexOutHalfedgeIterator<V,E,F,H> VertexOutHalfedgeIterator;
public:
};


typedef CSphericalHarmonicMapMesh<CSHMVertex, CSHMEdge, CSHMFace, CSHMHalfEdge> CSHMMesh;

} //namespace Holomorphy

unsigned long long Holomorphy::CSHMMesh::m_input_traits  = 0;
unsigned long long Holomorphy::CSHMMesh::m_output_traits = VERTEX_U;

} //namespace MeshLib

#endif  _SPHERICAL_HARMONIC_MAP_MESH_