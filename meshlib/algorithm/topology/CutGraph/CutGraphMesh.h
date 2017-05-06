/*!
*      \file CutGraphMesh.h
*      \brief Mesh for Computing Cut Graph
*	   \author David Gu
*      \date Document 12/25/2010
*
*/
#ifndef  _CUT_GRAPH_MESH_H_
#define  _CUT_GRAPH_MESH_H_

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
class CCutGraphVertex : public CVertex
{
public:

	/*!
	 *	CCutGraphVertex constructor
	 */
	CCutGraphVertex() { m_touched = false; m_valence = 0; m_father = NULL;};
	/*!
	 *	CCutGraphVertex destructor
	 */
	~CCutGraphVertex() {};

	/*!	Vertex touched trait
	 */
	bool &touched() { return m_touched; };
	/*!	Vertex topological valence
	 */
	int     & valence() { return m_valence; };
	/*! father in the tree
	 */
	CCutGraphVertex * &father() { return m_father; };

protected:	
	/*! touched */
	bool  m_touched;
	/*! topological valence */
	int   m_valence;
	/*! vertex father in the tree */
	CCutGraphVertex * m_father; 
};


/*! \brief CCutGraphEdge class
*
* Edge class for computing cut graph
*/

class CCutGraphEdge : public  CEdge
  {
  public:

  /*! CCutGraphEdge constructor
  */
	 CCutGraphEdge() { m_sharp = false; };
  /*! CCutGraphEdge destructor
  */
    ~CCutGraphEdge(){};

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
inline	void  CCutGraphEdge::_to_string()
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


/*! \brief CCutGraphFace class
*
* Face class for computing spherical harmonioc maps
*/

class CCutGraphFace : public  CFace
  {
  public:

  /*! CCutGraphFace constructor
  */
	 CCutGraphFace() { m_touched = false; };
  /*! CCutGraphFace destructor
  */
    ~CCutGraphFace(){};
	/*!	Face touched flat
	*/
	bool& touched() { return m_touched; };

  protected: 
	/* Face touched flag */
	bool m_touched;
};

/*! \brief CCutGraphMesh class
*
*  Mesh class for computing cut graph 
*/
template<typename V, typename E, typename F, typename H>
class CCutGraphMesh : public CBaseMesh<V,E,F,H>
{
public:
	typedef V	CVertex;
	typedef E	CEdge;
	typedef F	CFace;
	typedef H	CHalfEdge;

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
};


typedef CCutGraphMesh<CCutGraphVertex, CCutGraphEdge, CCutGraphFace, CHalfEdge> CutGraphMesh;

} //namespace Topology

unsigned long long Topology::CutGraphMesh::m_input_traits  = 0;
unsigned long long Topology::CutGraphMesh::m_output_traits = EDGE_SHARP;

} //namespace MeshLib

#endif  _CUT_GRAPH_MESH_