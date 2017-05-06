/*!
*      \file ConvexCapMesh.h
*      \brief Mesh for Convex Cap Embedding 
*	   \author David Gu
*      \date Documented 03/21/2013
*
*/

#ifndef  _CONVEX_CAP_MESH_H_
#define  _CONVEX_CAP_MESH_H_

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
*	\brief CConvexCapVertex class
*
*	Vertex class for Convex Cap
*   adding vertex height, vertex curvature traits 
*/
class CConvexCapVertex : public CVertex
{

public:
	/*!
	 *	CConvexCapVertex constructor
	 */
	CConvexCapVertex() {};
	/*!
	 *	CConvexCapVertex destructor
	 */
	~CConvexCapVertex() {};

	/*!
	 *	Write vertex uv to vertex string
	 */
	void _to_string();
	
	double & h() { return m_h; };
	double & k() { return m_k; };

protected:
	/*! Vertex height */
	double m_h;
	/*! Vertex curvature */
	double m_k;
};


//converting vertex uv trait to string

/*! write v->h to the string 
*/
inline void CConvexCapVertex::_to_string()
{
	CParser parser( m_string );
	parser._removeToken( "h" );
	parser._toString( m_string );
	std::stringstream iss;
	iss << "h=(" << m_h << ")";

	if( m_string.length() > 0 )
	{
		m_string += " ";
	}
	m_string += iss.str();
}

/*!
*	\brief CConvexCapEdge class
*
*	Edge class for Convex Cap
*   adding edge weight trait
*/
class CConvexCapEdge : public  CEdge
  {
  public:
    /*!	CConvexCapEdge constructor
	 */
	 CConvexCapEdge() { m_l=0; };
    /*!	CConvexCapEdge destructor
	 */
    ~CConvexCapEdge(){};
	/*!	edge target length
	 */
	double & l() { return m_l; };

	/*!	edge dihedral angle
	 */
	double & angle() { return m_angle; };

	/*! read edge traits from string, such as edge length
	 */
	void _from_string();

	/*!	edge source length
	 */
	double & source_l() { return m_source_l; };

	/*!	edge target length
	 */
	double & target_l() { return m_target_l; };

  protected:
	/*! edge target length trait */
	double   m_l;
	double   m_source_l;
	double   m_target_l;
	
	/*! edge dihedral angle */
	double   m_angle;
};

//read edge length from string
/*!	Read edge->length from the edge->string
 *
 */
inline void CConvexCapEdge::_from_string()
{

  CParser parser( m_string );

  for( std::list<CToken*>::iterator iter = parser.tokens().begin() ; iter != parser.tokens().end(); ++ iter )
  {
	  CToken * token = *iter;

	  if( token->m_key == "l" )
	  {
		 std::string line = strutil::trim( token->m_value, "()");
		 m_l = strutil::parseString<double>(line) ;		
	  }
  }

}


/*!
*	\brief CConvexCapFace class
*
*	Face class for convex cap
*   adding face normal trait
*/
class CConvexCapFace : public  CFace
  {
  public:
    /*!	CConvexCapFace constructor
	 */
	 CConvexCapFace() {};
    /*!	CConvexCapFace destructor
	 */
    ~CConvexCapFace(){};
	/*!	Face normal
	 */
	CPoint & normal() { return m_normal; };

  protected:
	CPoint m_normal;
};


/*!
*	\brief CConvexCapHalfEdge class
*
*	HalfEdge class for Convex Cap
*   adding Corner Angle trait
*/
class CConvexCapHalfEdge : public  CHalfEdge
  {
  public:
    /*!	CConvexCapHalfEdge constructor
	 */
	 CConvexCapHalfEdge() {};
    /*!	CConvexCapHalfEdge destructor
	 */
    ~CConvexCapHalfEdge(){};

	/*!	omega angle trait
	 */
	double & omega() { return m_omega; };

	/*!	alpha angle trait
	 */
	double & alpha() { return m_alpha; };

  protected:

	  /*! omega angle trait */
	double m_omega;

	  /*! alpha angle trait */
	double m_alpha;
};

/*-------------------------------------------------------------------------------------------------------------------------------------

	Convex Cap Mesh Class

--------------------------------------------------------------------------------------------------------------------------------------*/
/*!
 *	\brief CConvexCapMesh class
 *
 *	Mesh class for convex cap purpose
 */
template<typename V, typename E, typename F, typename H>
class CConvexCapMesh : public CBaseMesh<V,E,F,H>
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
typedef CConvexCapMesh<CConvexCapVertex, CConvexCapEdge, CConvexCapFace, CConvexCapHalfEdge> CCCMesh;	
/*! CHMMesh has no input traits, and has VERTEX_UV output traits
 */
unsigned long long CCCMesh::m_input_traits  = EDGE_LENGTH;
unsigned long long CCCMesh::m_output_traits = 0;
};
#endif  _CONVEX_CAP_MESH_H_