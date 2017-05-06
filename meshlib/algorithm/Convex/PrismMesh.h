/*!
*      \file PrismMesh.h
*      \brief Mesh for Prism
*	   \author David Gu
*      \date Documented 03/21/2013
*
*/

#ifndef  _PRISM_MESH_H_
#define  _PRISM_MESH_H_

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
/*!
*	\brief CPrismVertex class
*
*	Vertex class for Convex Cap
*   adding vertex height, vertex curvature traits 
*/
class CPrismVertex : public CVertex
{

public:
	/*!
	 *	CPrismVertex constructor
	 */
	CPrismVertex() {};
	/*!
	 *	CPrismVertex destructor
	 */
	~CPrismVertex() {};

};


/*!
*	\brief CPrismEdge class
*
*	Edge class for prism
*/
class CPrismEdge : public  CEdge
  {
  public:
    /*!	CPrismEdge constructor
	 */
	 CPrismEdge() { m_l=0; };
    /*!	CPrismEdge destructor
	 */
    ~CPrismEdge(){};
	/*!	edge target length
	 */
	double & l() { return m_l; };

	/*!	edge dihedral angle
	 */
	double & angle() { return m_angle; };

  protected:
	/*! edge length */
	double   m_l;
	
	/*! edge dihedral angle */
	double   m_angle;
};



/*!
*	\brief CPrismFace class
*
*	Face class for convex cap
*   adding face normal trait
*/
class CPrismFace : public  CFace
  {
  public:
    /*!	CPrismFace constructor
	 */
	 CPrismFace() {};
    /*!	CPrismFace destructor
	 */
    ~CPrismFace(){};
	/*!	Face normal
	 */
	CPoint & normal() { return m_normal; };

  protected:
	CPoint m_normal;
};


/*!
*	\brief CPrismHalfEdge class
*
*	HalfEdge class for Prism
*   adding Corner Angle trait
*/
class CPrismHalfEdge : public  CHalfEdge
  {
  public:
    /*!	CPrismHalfEdge constructor
	 */
	 CPrismHalfEdge() {};
    /*!	CPrismHalfEdge destructor
	 */
    ~CPrismHalfEdge(){};
	/*!	Corner angle
	 */
	double & angle() { return m_angle; };

  protected:
	  double m_angle;

};

/*-------------------------------------------------------------------------------------------------------------------------------------

	Prism Mesh Class

--------------------------------------------------------------------------------------------------------------------------------------*/
/*!
 *	\brief CPrismMesh class
 *
 *	Mesh class for prism
 */
template<typename V, typename E, typename F, typename H>
class CPrismMesh : public CBaseMesh<V,E,F,H>
{
public:
	typedef CBoundary<V,E,F,H> CBoundary;
	typedef CLoop<V,E,F,H> CLoop;
	typedef MeshVertexIterator<V,E,F,H> MeshVertexIterator;
	typedef MeshEdgeIterator<V,E,F,H> MeshEdgeIterator;
	typedef VertexVertexIterator<V,E,F,H> VertexVertexIterator;
	typedef VertexEdgeIterator<V,E,F,H> VertexEdgeIterator;
	typedef VertexInHalfedgeIterator<V,E,F,H> VertexInHalfedgeIterator;
	typedef MeshFaceIterator<V,E,F,H> MeshFaceIterator;
	typedef FaceVertexIterator<V,E,F,H> FaceVertexIterator;
	typedef VertexFaceIterator<V,E,F,H> VertexFaceIterator;
	typedef FaceHalfedgeIterator<V,E,F,H> FaceHalfedgeIterator;
};

/*! Mesh class for CHarmonicMapper class, Abbreviated as 'CHMMesh'
 */
typedef CPrismMesh<CPrismVertex, CPrismEdge, CPrismFace, CPrismHalfEdge> CPrMesh;	
/*! CHMMesh has no input traits, and has VERTEX_UV output traits
 */
unsigned long long CPrMesh::m_input_traits  = EDGE_LENGTH;
unsigned long long CPrMesh::m_output_traits = 0;
};
#endif  _PRISM_MESH_H_