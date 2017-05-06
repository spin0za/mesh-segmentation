/*!
*      \file DynamicWedgeMesh.h
*      \brief Dynamic Mesh for Wedge mesh 
*	   \author David Gu
*      \date Document 10/08/2013
*
*/


#ifndef  _DYNAMIC_WEDGE_MESH_H_
#define  _DYNAMIC_WEDGE_MESH_H_

#include <map>
#include <vector>

#include "Mesh/BaseMesh.h"
#include "Mesh/DynamicMesh.h"
#include "Mesh/Boundary.h"
#include "Mesh/Vertex.h"
#include "Mesh/HalfEdge.h"
#include "Mesh/Edge.h"
#include "Mesh/Face.h"
#include "mesh/iterators.h"
#include "Parser/parser.h"
#include "Parser/traits_io.h"

namespace MeshLib
{
namespace Topology
{
  class CBaseWedge;

/*-------------------------------------------------------------------------------------------

 Vertex Trait 

--------------------------------------------------------------------------------------------*/
  /*!	\brief CWedgeVertex class
   *
   *   Vertex class for Wedge mesh
   *   adding two traits, the topological valence and the father vertex id.
   */
  class CDynamicWedgeVertex : public  CVertex
  {
  public:
    /*! CWedgeVertex constructor. */
	CDynamicWedgeVertex() { m_valence = 0; m_father = 0; };  
	 /*! CWedgeVertex destructor. */
	~CDynamicWedgeVertex() {};
	
	/*! convert traits to vertex string. */
	void _to_string();
	
	/*! topological valence. */
	int & valence() { return m_valence; };
	/*! the father vertex of current vertex. */
	int & father()  { return m_father;  };

  protected:
	/*! topological valence. */
	int     m_valence;
	/*! vertex father id. */ 
	int     m_father;
  };

/*! Add vertex father trait to vertex string. */
inline void CDynamicWedgeVertex::_to_string()
{
  CParser parser( m_string );
  parser._removeToken( "father" );

  parser._toString( m_string );
  std::stringstream iss;
  iss << "father=(" << m_father << ")";

  if( m_string.length() > 0 ) m_string += " ";
  m_string += iss.str();
};


/*-------------------------------------------------------------------------------------------

 Edge Trait 

--------------------------------------------------------------------------------------------*/
/*! \brief CDynamicWedgeEdge class
 *
 * Edge class for Wedge Mesh
 * Add edge sharp field.
 */
class CDynamicWedgeEdge : public  CEdge
  {
  public:
	/*! CWedgeEdge constructor. */
    CDynamicWedgeEdge() { m_sharp = false; };
	/*! CWedgeEdge destructor. */
    ~CDynamicWedgeEdge(){};
	/*! read sharp trait from edge string.*/
	void _from_string();
	/*! whether the current edge is sharp.*/
	bool & sharp() { return m_sharp; };

  protected:
	  /*! whether current edge is sharp. */
	  bool m_sharp;
  
};

/*!	Read edge sharp trait from the edge string.
 */
inline void CDynamicWedgeEdge::_from_string( )
{
	CParser parser( m_string );

	  for( std::list<CToken*>::iterator iter = parser.tokens().begin() ; iter != parser.tokens().end(); ++ iter )
	  {
		  CToken * token = *iter;
		  if( token->m_key == "sharp" )
		  {
			  m_sharp = true;
			  break;
		  }
	  }
};


/*-------------------------------------------------------------------------------------------

 Half Edge Trait 

--------------------------------------------------------------------------------------------*/
/*! \brief CWedgeHalfEdge class
 *
*    HalfEdge class for WedgeMesh
*	Add wedge trait, each halfedge belongs to a unique wedge.
*/
class CDynamicWedgeHalfEdge : public  CHalfEdge
{
public:
	/*! CWedgeHalfEdge constructor */
  CDynamicWedgeHalfEdge() { m_wedge = NULL; m_child = NULL; };
  /*! CWedgeHalfEdge destructor */
  ~CDynamicWedgeHalfEdge() {};
  /*! the wedge current corner belongs to.*/
  CBaseWedge* & wedge() { return m_wedge; };
  /*! the child halfedge */
  CDynamicWedgeHalfEdge * & child() { return m_child; };

protected:
  /*! the wedge current corner belongs to.*/
  CBaseWedge * m_wedge;
  /*! the child of the current halfedge */
  CDynamicWedgeHalfEdge *  m_child;
};

/*! \brief CSliceMesh class
 *
 *  Mesh class for Wedge Mesh, for slciing a mesh with sharp edges to an open mesh.
 */
template<typename V, typename E, typename F, typename H>
class CDynamicSliceMesh : public CDynamicMesh<V,E,F,H>
{
public:
	typedef V CVertex;
	typedef E CEdge;
	typedef F CFace;
	typedef H CHalfEdge;

	typedef CBoundary<V,E,F,H>             CBoundary;
	typedef CLoop<V,E,F,H>				  CLoop;

	typedef MeshVertexIterator<V,E,F,H>   MeshVertexIterator;
	typedef MeshFaceIterator<V,E,F,H>	  MeshFaceIterator;
	typedef MeshEdgeIterator<V,E,F,H>     MeshEdgeIterator;
	typedef VertexVertexIterator<V,E,F,H> VertexVertexIterator;
	typedef FaceHalfedgeIterator<V,E,F,H> FaceHalfedgeIterator;
	typedef FaceVertexIterator<V,E,F,H>   FaceVertexIterator;
	typedef VertexEdgeIterator<V,E,F,H>   VertexEdgeIterator;
	typedef VertexInHalfedgeIterator<V,E,F,H>   VertexInHalfedgeIterator;
};


typedef CDynamicSliceMesh<CDynamicWedgeVertex, CDynamicWedgeEdge, CFace, CDynamicWedgeHalfEdge> CDynamicSMesh;

} //namespace Topology

unsigned long long Topology::CDynamicSMesh::m_input_traits  = EDGE_SHARP;
unsigned long long Topology::CDynamicSMesh::m_output_traits = VERTEX_FATHER;

}//namespace MeshLib
#endif  _DYNAMIC_WEDGE_MESH_H_