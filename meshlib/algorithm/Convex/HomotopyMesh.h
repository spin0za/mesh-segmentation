/*!
*      \file HomotopyMesh.h
*      \brief Mesh for Convex Polyhedron Embedding 
*	   \author David Gu
*      \date Documented 03/21/2013
*
*/

#ifndef  _HOMOTOPY_MESH_H_
#define  _HOMOTOPY_MESH_H_

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
class CHomotopyVertex : public CVertex
{

public:
	/*!
	 *	CHomotopyVertex constructor
	 */
	CHomotopyVertex() { m_index = 0; m_touched = false; };
	/*!
	 *	CHomotopyVertex destructor
	 */
	~CHomotopyVertex() {};

	/*!
	 *	Write vertex uv to vertex string
	 */
	void _to_string();
	
	int & idx() { return m_index; };

	bool & touched() { return m_touched; };

protected:
	/*! vertex index */
	int m_index;
	/*! touched */
	bool m_touched;
};

/*! write v->index to the string 
*/
inline void CHomotopyVertex::_to_string()
{
	CParser parser( m_string );
	parser._removeToken( "h" );
	parser._toString( m_string );
	std::stringstream iss;
	iss << "h=(" << m_index << ")";

	if( m_string.length() > 0 )
	{
		m_string += " ";
	}
	m_string += iss.str();
};

/*!
*	\brief CHomotopyEdge class
*
*	Edge class for Convex Embedding
*   adding edge weight trait
*/
class CHomotopyEdge : public  CEdge
  {
  public:
    /*!	CConvexCapEdge constructor
	 */
	 CHomotopyEdge() { m_l=0; m_index = 0; m_touched = false; };
    /*!	CConvexCapEdge destructor
	 */
    ~CHomotopyEdge(){};
	/*!	edge target length
	 */
	double & l() { return m_l; };
	/*!	edge index
	 */
	int & idx() { return m_index; };

	/*! read edge traits from string, such as edge length
	 */
	void _from_string();

	bool& touched() { return m_touched; };
	
  protected:
	/*! edge target length trait */
	double   m_l;
	int		 m_index;
	bool     m_touched;
};

//read edge length from string
/*!	Read edge->length from the edge->string
 *
 */
inline void CHomotopyEdge::_from_string()
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

};


/*!
*	\brief CHomotopyFace class
*
*	Face class for Convex Embedding
*
*/
class CHomotopyFace : public  CFace
  {
  public:
    /*!	CHomotopyFace constructor
	 */
	 CHomotopyFace() { m_index = 0; m_touched = false; };
    /*!	CHomotopyFace destructor
	 */
    ~CHomotopyFace(){};
	/*!	Face index
	 */
	int & idx() { return m_index; };

	/*!	If the Face has been touched
	 */
	bool& touched() { return m_touched; };
	
  protected:
	int		 m_index;
	bool     m_touched;
};



/*-------------------------------------------------------------------------------------------------------------------------------------

	Homotopy Mesh Class

--------------------------------------------------------------------------------------------------------------------------------------*/
/*!
 *	\brief CHomotopyMesh class
 *
 *	Mesh class for convex embedding
 */
template<typename V, typename E, typename F, typename H>
class CHomotopyMesh : public CBaseMesh<V,E,F,H>
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
typedef CHomotopyMesh<CHomotopyVertex, CHomotopyEdge, CHomotopyFace, CHalfEdge> CHTMesh;	
/*! CHMMesh has no input traits, and has VERTEX_UV output traits
 */
unsigned long long CHTMesh::m_input_traits  = EDGE_LENGTH;
unsigned long long CHTMesh::m_output_traits = 0;

}
#endif  _HOMOTOPY_MESH_H_