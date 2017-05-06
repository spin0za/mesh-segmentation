/*!
*      \file ShortestPathMesh.h
*      \brief Mesh for Shortest Path
*	   \author David Gu
*      \date Document 10/11/2010
*
*/

#ifndef  _SHORTEST_PATH_MESH_H_
#define  _SHORTEST_PATH_MESH_H_

#include <map>
#include <vector>

#include "Mesh/BaseMesh.h"
#include "Mesh/Vertex.h"
#include "Mesh/HalfEdge.h"
#include "Mesh/Edge.h"
#include "Mesh/Face.h"
#include "mesh/iterators.h"
#include "mesh/boundary.h"
#include "Parser/traits_io.h"

namespace MeshLib
{

/*!	\brief CSPEdge class
 *
 *  Edge class for shortest path 
 *  adding edge sharp trait to  label shortest paths 
 */
class CSPEdge : public  CEdge
  {
  public:
	  /*! whether the edge is on the shortest path */
	bool     m_sharp;
  public:
	  /*! CSPEdge constructure */
    CSPEdge() { m_sharp = false; };
	  /*! CSPEdge destructure */
    ~CSPEdge(){};
	/*! whether the edge is on the shorest path */
	bool & sharp() { return m_sharp; };

	/*! write to string */
	void _to_string();
};

/*! write vertex traits m_sharp to the string */
inline	void  CSPEdge::_to_string()
	{
		CParser parser( m_string );
		parser._removeToken( "sharp" );
		parser._toString( m_string );
		
		std::string line;
		std::stringstream iss(line);
		if( m_sharp )
		{
			iss << "sharp";
		}
		if( m_string.length() > 0 )
		{
			m_string += " ";
		}
		m_string += iss.str();

	};


/*!	\brief CSPVertex class
 *
 *	Vertex class for shortest path
 *  adding the following traits: vertex boundary index, parent in breadth first searching tree, touched in searching process 
 */
class CSPVertex : public CVertex
{
public:
	/*! CSPVertex constructor
	*/
	CSPVertex() { m_index = 0; m_touched = false; m_parent = NULL; };
	/*! CSPVertex destructor
	*/
	~CSPVertex() {};
	/*!
	 *  if the vertex is on the boundary, the boundary loop index 
	 */
	int & idx() { return m_index; };
	/*!
	 *	whether the vertex has been touched in the breadth first searching process
	 */
	bool& touched() { return m_touched; };
	/*!
	 *	the parent of the vertex in the breadth first searching tree
	 */
	CSPVertex * &parent() { return m_parent;};
	/*!
	 *	the edge pointing to the parent
	 */
	CSPEdge   * &bridge() { return m_bridge; };

protected:
	/*! vertex index */
	int m_index;
	/*! whether the vertex has been accessed */
	bool m_touched;
	/*! the parent of the current vertex */
	CSPVertex * m_parent;
	/*! the edge pointing to the parent  */
	CSPEdge   * m_bridge;
	 
};

/*!	\brief CShortestPathMesh class
 *
 *	Mesh class for computing shortest paths between two boundary loops
 */

template<typename V, typename E, typename F, typename H>
class CShortestPathMesh : public CBaseMesh<V,E,F,H>
{
public:
	typedef CBoundary<V,E,F,H> CBoundary;
	typedef CLoop<V,E,F,H> CLoop;
	typedef MeshVertexIterator<V,E,F,H> MeshVertexIterator;
	typedef MeshEdgeIterator<V,E,F,H> MeshEdgeIterator;
	typedef VertexVertexIterator<V,E,F,H> VertexVertexIterator;
	typedef VertexEdgeIterator<V,E,F,H> VertexEdgeIterator;
public:
};


typedef CShortestPathMesh<CSPVertex, CSPEdge, CFace, CHalfEdge> CSPMesh;

unsigned long long CSPMesh::m_input_traits  = 0;
unsigned long long CSPMesh::m_output_traits = EDGE_SHARP;
};
#endif  _SHORTEST_PATH_MESH_H_