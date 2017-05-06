/*!
*      \file DynamicCutGraphMesh.h
*      \brief Mesh for Computing Cut Graph on Dynamic Mesh
*	   \author David Gu
*      \date Document 10/08/2013
*
*/
#ifndef  _DYNAMIC_CUT_GRAPH_MESH_H_
#define  _DYNAMIC_CUT_GRAPH_MESH_H_

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

namespace Topology
{

/*!
*	\brief CSHMVertex class
*
*	Vertex class for spherical harmonic map
*/
class GCutGraphVertex : public CVertex
{
public:

	/*!
	 *	GCutGraphVertex constructor
	 */
	GCutGraphVertex() { m_touched = false; m_valence = 0; m_father = NULL;};
	/*!
	 *	GCutGraphVertex destructor
	 */
	~GCutGraphVertex() {};

	/*!	Vertex cg_touched trait
	 */
	bool &touched() { return m_touched; };
	/*!	Vertex topological cg_valence
	 */
	int  &valence() { return m_valence; };
	/*! father in the tree
	 */
	GCutGraphVertex * &father() { return m_father; };

protected:	
	/*! cg_touched */
	bool  m_touched;
	/*! topological cg_valence */
	int   m_valence;
	/*! vertex cg_father in the tree */
	GCutGraphVertex * m_father; 
};



/*! \brief GCutGraphEdge class
*
* Edge class for computing cut graph
*/
class GCutGraphEdge : public  CEdge
  {
  public:

  /*! GCutGraphEdge constructor
  */
	 GCutGraphEdge() { m_sharp = false; };
  /*! GCutGraphEdge destructor
  */
    ~GCutGraphEdge(){};

	/*! Edge sharp */
	bool & sharp() { return m_sharp; };
	/*! Edge mark */
	int  & mark()  { return m_mark; };

	/*! output the string
	 */
	void _to_string();

  protected: 
	/* edge sharp flat */
	bool   m_sharp;
	/* mark */
	int m_mark;
};

/*! write vertex traits m_sharp to the string */
inline	void  GCutGraphEdge::_to_string()
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


/*! \brief GCutGraphFace class
*
* Face class for computing spherical harmonioc maps
*/
class GCutGraphFace : public  CFace
  {
  public:

  /*! GCutGraphFace constructor
  */
	 GCutGraphFace() { m_touched = false; };
  /*! GCutGraphFace destructor
  */
    ~GCutGraphFace(){};
	/*!	Face cg_touched flat
	*/
	bool& touched() { return m_touched; };

  protected: 
	/* Face touched flag */
	bool m_touched;
};

class GCutGraphHalfEdge : public  CHalfEdge
  {
  public:

  /*! GCutGraphFace constructor
  */
	 GCutGraphHalfEdge() { };
  /*! GCutGraphFace destructor
  */
    ~GCutGraphHalfEdge(){};
};



/*! \brief CDynamicCutGraphMesh class
*
*  Mesh class for computing cut graph on dynamic mesh
*/


template<typename V, typename E, typename F, typename H>
class CDynamicCutGraphMesh : public CDynamicMesh<V,E,F,H>
{
public:
	
	typedef	V	CVertex;
	typedef	E	CEdge;
	typedef	F	CFace;
	typedef	H	CHalfEdge;

	typedef CBoundary<V,E,F,H>  CBoundary;
	typedef CLoop<V,E,F,H> CLoop;
	typedef MeshVertexIterator<V,E,F,H> MeshVertexIterator;
	typedef MeshFaceIterator<V,E,F,H> MeshFaceIterator;
	typedef MeshEdgeIterator<V,E,F,H> MeshEdgeIterator;
	typedef VertexVertexIterator<V,E,F,H> VertexVertexIterator;
	typedef VertexEdgeIterator<V,E,F,H> VertexEdgeIterator;
	typedef VertexFaceIterator<V,E,F,H> VertexFaceIterator;
	typedef VertexOutHalfedgeIterator<V,E,F,H> VertexOutHalfedgeIterator;
	typedef FaceHalfedgeIterator<V,E,F,H> FaceHalfedgeIterator;
public:
	CDynamicCutGraphMesh()		{	}
	~CDynamicCutGraphMesh()	{	}
};



typedef CDynamicCutGraphMesh<GCutGraphVertex, GCutGraphEdge, GCutGraphFace, CHalfEdge> DynamicCutGraphMesh;

} //namespace Topology

unsigned long long Topology::DynamicCutGraphMesh::m_input_traits  = 0;
unsigned long long Topology::DynamicCutGraphMesh::m_output_traits = EDGE_SHARP;


} //namespace MeshLib

#endif  _DYNAMIC_CUT_GRAPH_MESH_H_