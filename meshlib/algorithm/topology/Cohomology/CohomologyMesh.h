/*!
*      \file CohomologyMesh.h
*      \brief Mesh for Cohomology
*	   \author David Gu
*      \date Document 10/11/2010
*
*/
#ifndef  _COHOMOLOGY_MESH_H_
#define  _COHOMOLOGY_MESH_H_

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

namespace Topology
{

/*!
*	\brief CHCFVertex class
*
*	Vertex class for harmonic closed form
*/
class CCHVertex : public CVertex
{
public:

	/*!
	 *	CHCFVertex constructor
	 */
	CCHVertex() { m_index = 0; };
	/*!
	 *	CHCFVertex destructor
	 */
	~CCHVertex() {};
	/*!
	 *	Vertex index
	 */
	int & idx() { return m_index; };
	/*!
	 *	Vertex father vertex
	 */
	int &     father() { return m_father; };
	/*!	Vertex texture uv coordinates
	 */
	CPoint2 & uv() { return m_uv; };
	/*!	Vertex function u
	*/
	double  &  u() { return m_u;  };
	/*!	Vertex topological valence
	 */
	int     & valence() { return m_valence; };
	/*!	whether the vertex has been accessed
	*/
	bool    & touched() { return m_touched; };
	/*! convert vertex traits to the vertex string */
	void  _to_string();
	/*! read vertex traits from the string */
	void  _from_string();

protected:	//output
	/*! vertex texture coordinates */
	CPoint2     m_uv;	
protected:
	/*! vertex's father vertex id */
	int         m_father;
	/*! whether the current vertex has been accessed */
	bool        m_touched; //for integration
	/*! vertex topological valence */
	int         m_valence; //topological valence
	/*! vertex function */
    double      m_u;
	/*  vertex index */
	int			m_index;
	
};

/*! Read vertex father from vertex string */
inline	void CCHVertex::_from_string()
{
	  CParser parser( m_string );

    std::list<CToken*> & tokens = parser.tokens();
    for( std::list<CToken*>::iterator titer = tokens.begin(); titer != tokens.end(); titer ++ )
    {
        CToken * pT = *titer;
		    if( pT->m_key == "father" )
		    {
				std::string str = strutil::trim( pT->m_value, "()");
				m_father = strutil::parseString<int>( str );
		    }
    }
};

/*! write vertex traits uv to the string */
inline	void  CCHVertex::_to_string()
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
/*! \brief CCHEdge class
*
* Edge class for computing harmonioc closed form
*/

class CCHEdge : public  CEdge
  {
  public:
  /*! CCHEdge constructor
  */
	 CCHEdge() { m_du = 0; };
  /*! CCHEdge destructor
  */
    ~CCHEdge(){};
	/*! Edge 1-form */
	double & du() { return m_du; };

	/*! convert edge traits to string */
	void _to_string();

  protected: //output
	 /*! edge 1-form */
    double   m_du;
};

/*! convert edge 1-form trait to edge string */
inline void CCHEdge::_to_string()
{
		CParser parser( m_string );
		parser._removeToken( "du" );

		parser._toString( m_string );
		
		std::string line;
		std::stringstream iss(line);
		iss << "du=(" << m_du << ")";

		if( m_string.length() > 0 )
		{
			m_string += " ";
			m_string += iss.str();
		}
		else
		{
			m_string = iss.str();
		}
};


/*! \brief CChomologyMesh class
*
*  Mesh class for computing closed harmonic 1-form 
*/
template<typename V, typename E, typename F, typename H>
class CCohomologyMesh : public CBaseMesh<V,E,F,H>
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
	typedef VertexVertexIterator<V,E,F,H> VertexVertexIterator;
	typedef VertexEdgeIterator<V,E,F,H> VertexEdgeIterator;
	typedef VertexOutHalfedgeIterator<V,E,F,H> VertexOutHalfedgeIterator;
public:
};


typedef CCohomologyMesh<CCHVertex, CCHEdge, CFace, CHalfEdge> CCHMesh;

} //namespace Topology

unsigned long long Topology::CCHMesh::m_input_traits  = VERTEX_FATHER;
unsigned long long Topology::CCHMesh::m_output_traits = VERTEX_UV | EDGE_DU;

} //namespace MeshLib

#endif  _COHOMOLOGY_MESH_H_