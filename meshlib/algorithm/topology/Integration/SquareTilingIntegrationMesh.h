/*! \file SquareTilingIntegrationMesh.h
*  \brief Mesh class for square tiling
*  \author David Gu
*  \data Documented on 10/12/2010
*/

/*******************************************************************************
*      Square Tiling Integration Class
*
*       Copyright (c) Stony Brook University
*
*    Purpose:
*
*       Mesh Integration Trait
* 
*       David Gu October 25, 2014,  gu@cs.stonybrook.edu
*
*******************************************************************************/


#ifndef  _SQUARE_TILING_INTEGRATION_MESH_H_
#define  _SQUARE_TILING_INTEGRATION_MESH_H_

#include <map>
#include <vector>
#include <list>

#include "Mesh/Vertex.h"
#include "Mesh/Edge.h"
#include "Mesh/Face.h"
#include "Mesh/HalfEdge.h"

#include "Mesh/BaseMesh.h"
#include "Mesh/iterators.h"
#include "Parser/parser.h"
#include "Parser/traits_io.h"

namespace MeshLib
{
namespace Topology
{
/*---------------------------------------------------------------------------------------------------------------------------------------

	Integration Vertex Trait

----------------------------------------------------------------------------------------------------------------------------------------*/
/*!
* \brief CSquareTilingIntegrationVertex class
*
*  Vertex class for computing integration
*  The following traits are added
*  Vertex texture coordinates m_uv
*  Wheter the vertex has been accessed in the integration process m_touched
*  The father vertex id m_father
*/
class CSquareTilingIntegrationVertex : public  CVertex
{

  public:
	  /*! Wheter the vertex has been accessed in the integration process m_touched */
	  bool & touched()  { return m_touched; };
	  /*! Vertex texture coordinates */
	  double & u()    { return m_u; };
	  /*! Vertex father id */
	  int     & father() { return m_father; };
protected:
   /*! Wheter the vertex has been accessed in the integration process m_touched */
    bool        m_touched;
	 /*! Vertex texture coordinates */
    double		m_u;
	/*! Vertex father id */
	int         m_father;
  
  public:
	/*! CIntegrationVertex constructor */
    CSquareTilingIntegrationVertex() { m_father = 0; m_touched = false;m_u = 0;};
	/*! CIntegrationVertex destructor */
    ~CSquareTilingIntegrationVertex(){};
	/*! Save vertex uv coordinates to the vertex string */
	void _to_string();
	/*! Get vertex father from the vertex string */
	void _from_string();

 };

// Get vertex father from the vertex string 

inline void CSquareTilingIntegrationVertex::_from_string( )
{
	  CParser parser( m_string );

	  for( std::list<CToken*>::iterator iter = parser.tokens().begin() ; iter != parser.tokens().end(); ++ iter )
	  {
		  CToken * token = *iter;
		  if( token->m_key == "father" )
		  {
			  std::string line = strutil::trim( token->m_value, "()");
			  m_father = strutil::parseString<int>( line );	
		  }
		  if( token->m_key == "u" )
		  {
				std::string str = strutil::trim( token->m_value, "()");
				m_u = strutil::parseString<double>( str );
		  }
	  }
};

// Save vertex uv coordinates to the vertex string 
inline void CSquareTilingIntegrationVertex::_to_string()
{
	CParser parser( m_string );
	parser._removeToken( "u" );

	parser._toString( m_string );

	std::string line;
	std::stringstream iss( line );

	iss << "u=(" << m_u << ")";

	 if( m_string.size() > 0 )
    {
		  m_string += " ";
    }
    m_string += iss.str();
};



/*---------------------------------------------------------------------------------------------------------------------------------------

	Square Tiling Integration Edge Trait

----------------------------------------------------------------------------------------------------------------------------------------*/
/*! \brief CSquareTilingIntegrationEdge class
 *
 *	Edge class for computing integration
 *  The edge harmonic 1-form trait is added : m_du
 */
class CSquareTilingIntegrationEdge : public  CEdge
{
  public:
	/*!  The edge combinatorial harmonic 1-form trait */
	  double  & du() { return m_du; };
protected:
	/*!  The edge combinatorial harmonic 1-form trait */
    double		m_du;

 public:
	/*! Read edge combinatorial harmonic 1-form m_du from the string with the key token "du" */
	void _from_string();
	/*! write edge combinatorial harmonic 1-form m_du to the string with the key token "du" */
	void _to_string();
};

// Read edge combinatorial harmonic 1-form m_du from the string with the key token "du" 
inline void CSquareTilingIntegrationEdge::_from_string()
{
	  CParser parser( m_string );
	
	  for( std::list<CToken*>::iterator iter = parser.tokens().begin() ; iter != parser.tokens().end(); ++ iter )
	  {
		  CToken * token = *iter;
		  if( token->m_key == "du" )
		  {
				std::string str = strutil::trim( token->m_value, "()");
				m_du = strutil::parseString<double>( str );
		 }
	  }
};

//write holomorphic 1-form trait m_duv to the string "duv"

inline void CSquareTilingIntegrationEdge::_to_string()
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

/*---------------------------------------------------------------------------------------------------------------------------------------

	Square Tiling Integration Face Trait

----------------------------------------------------------------------------------------------------------------------------------------*/
/*! \brief CSquareTilingIntegrationFace class
 *
 *	Face class for computing integration
 *  The face u trait : m_u
 */
class CSquareTilingIntegrationFace : public  CFace
{
  public:
	  /*!  The face u coordinate */
	  double  & u() { return m_u; };
	  /*! whether the face has been touched */
	  bool    & touched() { return m_touched; };

protected:
	/*!  The face coordinate trait */
    double		m_u;
	/*!  The face has been touched */
	bool        m_touched;

 public:
	/*! Read face u-coordinates from the string with the key token "u" */
	void _from_string();
	/*! write face u-coordinate to the string with the key token "u" */
	void _to_string();
};

// Read face u-coordinates from the string with the key token "u" 

inline void CSquareTilingIntegrationFace::_from_string()
{
	  CParser parser( m_string );
	
	  for( std::list<CToken*>::iterator iter = parser.tokens().begin() ; iter != parser.tokens().end(); ++ iter )
	  {
		CToken * token = *iter;
		if( token->m_key == "u" )
		{
			std::string str = strutil::trim( token->m_value, "()");
			m_u = strutil::parseString<double>( str );
		}
	  }
};

// write face u-coordinate to the string with the key token "u" 

inline void CSquareTilingIntegrationFace::_to_string()
{
	CParser parser( m_string );
	parser._removeToken( "u" );

	parser._toString( m_string );
	
	std::string line;
	std::stringstream iss(line);
	iss << "u=(" << m_u << ")";

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


/*---------------------------------------------------------------------------------------------------------------------------------------

	Integration Trait

----------------------------------------------------------------------------------------------------------------------------------------*/
/*! \brief CSquareTilingIntegrationMesh class
 *  Mesh class for integrating holomorphic 1-form on a mesh
*/
template<typename V, typename E, typename F, typename H>
class CSquareTilingIntegrationMesh : public CBaseMesh<V,E,F,H>
{
public:
	typedef V CVertex;
	typedef E CEdge;
	typedef F CFace;
	typedef H CHalfEdge;

	typedef MeshVertexIterator<V,E,F,H> MeshVertexIterator;
	typedef MeshFaceIterator<V,E,F,H> MeshFaceIterator;
	typedef MeshEdgeIterator<V,E,F,H> MeshEdgeIterator;
	typedef VertexVertexIterator<V,E,F,H> VertexVertexIterator;
	typedef VertexEdgeIterator<V,E,F,H> VertexEdgeIterator;
	typedef FaceHalfedgeIterator<V,E,F,H> FaceHalfedgeIterator;
	typedef FaceVertexIterator<V,E,F,H> FaceVertexIterator;
public:
};


typedef CSquareTilingIntegrationMesh<CSquareTilingIntegrationVertex, CSquareTilingIntegrationEdge, CSquareTilingIntegrationFace, CHalfEdge> CSTIMesh;

} //namespace Topology

unsigned long long Topology::CSTIMesh::m_input_traits  = VERTEX_FATHER | EDGE_DU;
unsigned long long Topology::CSTIMesh::m_output_traits = VERTEX_U | FACE_U;

} //namespace MeshLib

#endif  _SQUARE_TILING_INTEGRATION_TRAIT_H_