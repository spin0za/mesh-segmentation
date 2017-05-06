/*!
*      \file ExtremalLengthMesh.h
*      \brief Mesh for Computing Extremal Length 
*	   \author David Gu
*      \date Document 12/29/2010
*
*/
#ifndef  _EXTREMAL_LENGTH_MESH_H_
#define  _EXTREMAL_LENGTH_MESH_H_

#include <map>
#include <vector>
#include <queue>

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
*	\brief CELVertex class
*
*	Vertex class for extremal length
*/
class CELVertex : public CVertex
{
public:

	/*!
	 *	CELVertex constructor
	 */
	CELVertex() { m_index = 0; m_valence = 0; };
	/*!
	 *	CELVertex destructor
	 */
	~CELVertex() {};

	/*! 
	 *	Vertex index
	 */
	int &     idx()    { return m_index;  };
	
	/*! Topological valence
	 */
	int & valence()    { return m_valence; };

	/*! whether the vertex is fixed
	 */
	bool & fixed()    { return m_fixed; };

	/*!	Vertex texture uv coordinates
	 */
	CPoint2 & uv() { return m_uv; };
	/*!	Vertex function u
	*/
	double  &  u() { return m_u;  };
	/*! convert vertex traits to the vertex string */
	void  _to_string();

protected:	//output
	/*! vertex texture coordinates */
	CPoint2     m_uv;	
protected:
	/*! vertex index */
	int         m_index;
	/*! vertex topological valence */
	int         m_valence;
	/*! whether the current vertex is fiexed */
	bool        m_fixed; 
	/*! vertex function */
    double      m_u;
	
};


/*! write vertex traits uv to the string */
inline	void  CELVertex::_to_string()
	{
		CParser parser( m_string );
		parser._removeToken( "uv" );

		parser._toString( m_string );
		
		std::string line;
		std::stringstream iss(line);
		iss << "uv=(" << m_uv[0] << " " << m_uv[1] << ")";

		if( m_string.length() > 0 )
		{
			m_string += " ";
		}
		m_string += iss.str();

	};
/*! \brief CELEdge class
*
* Edge class for computing extremal length
*/

class CELEdge : public  CEdge
  {
  public:
  /*! CELEdge constructor
  */
	 CELEdge() { m_weight = 0; m_du = 0; m_sharp = false; };
  /*! CELEdge destructor
  */
    ~CELEdge(){};
	/*! Edge 1-form */
	double & du() { return m_du; };
	/*! Edge weight */
	double & weight() { return m_weight; };
	/*! Edge length */
	double & length() { return m_length; };

	/*! convert edge traits to string */
	void _to_string();

	/*! read from string */
	void _from_string();

	/*! sharp flag for the edge */
	bool & sharp() { return m_sharp; };

  protected: //output
	 /*! edge 1-form */
    double   m_du;
	/*! edge weight */
	double   m_weight;
	/*! edge sharp flat */
	bool     m_sharp;
	/*! edge length */
	double   m_length;
};

/*!	Read edge sharp trait from the edge string.
 */
inline void CELEdge::_from_string( )
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


/*! convert edge 1-form trait to edge string */
inline void CELEdge::_to_string()
{
		CParser parser( m_string );
		parser._removeToken( "du" );

		parser._toString( m_string );
		
		std::string line;
		std::stringstream iss(line);
		iss << "du=(" << m_du << ")";
		m_string = iss.str();
};

/*! \brief CELHalfEdge class
*
*  Halfedge class for computing extremal length
*/
class CELHalfEdge : public  CHalfEdge
  {
  public:
	/*! CHCFHalfEdge constructor */
	 CELHalfEdge() { m_angle = 0; };
	/*! CHCFHalfEdge destructor */
    ~CELHalfEdge(){};
	  /*! corner angle */
	double &   angle() { return m_angle; };
	
  //no output
  protected: 
	  /*! corner angle */
	  double   m_angle;

};

/*! \brief CExtremalLengthMesh class
*
*  Mesh class for computing extremal length 
*/
template<typename V, typename E, typename F, typename H>
class CExtremalLengthMesh : public CBaseMesh<V,E,F,H>
{
public:
	
	typedef V CVertex;
	typedef E CEdge;
	typedef F CFace;
	typedef H CHalfEdge;

	typedef CBoundary<V,E,F,H>  CBoundary;
	typedef CLoop<V,E,F,H> CLoop;
	typedef CLoopSegment<V,E,F,H> CLoopSegment;

	typedef MeshVertexIterator<V,E,F,H> MeshVertexIterator;
	typedef MeshFaceIterator<V,E,F,H> MeshFaceIterator;
	typedef MeshEdgeIterator<V,E,F,H> MeshEdgeIterator;
	typedef VertexVertexIterator<V,E,F,H> VertexVertexIterator;
	typedef VertexEdgeIterator<V,E,F,H> VertexEdgeIterator;
	typedef VertexOutHalfedgeIterator<V,E,F,H> VertexOutHalfedgeIterator;
public:
};

typedef CExtremalLengthMesh<CELVertex, CELEdge, CFace, CELHalfEdge> CELMesh;

} //namespace Holomorphy

unsigned long long Holomorphy::CELMesh::m_input_traits  = EDGE_SHARP;
unsigned long long Holomorphy::CELMesh::m_output_traits = EDGE_DU;

} //namespace MeshLib

#endif  _EXTREMAL_LENGTH_MESH_H_