/*! \file IntegrationMesh.h
*  \brief Mesh class for Integrating 1-forms on meshes
*  \author David Gu
*  \data Documented on 10/12/2010
*/

/*******************************************************************************
*      Integration Trait Class
*
*       Copyright (c) Stony Brook University
*
*    Purpose:
*
*       Mesh Integration Trait
* 
*       David Gu June 27, 2008,  gu@cs.stonybrook.edu
*
*******************************************************************************/


#ifndef  _INTEGRATION_MESH_H_
#define  _INTEGRATION_MESH_H_

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
* \brief CIntegrationVertex class
*
*  Vertex class for computing integration
*  The following traits are added
*  Vertex texture coordinates m_uv
*  Wheter the vertex has been accessed in the integration process m_touched
*  The father vertex id m_father
*/
class CIntegrationVertex : public  CVertex
{

  public:
	  /*! Wheter the vertex has been accessed in the integration process m_touched */
	  bool & touched()  { return m_touched; };
	  /*! Vertex texture coordinates */
	  CPoint2 & uv()    { return m_uv; };
	  /*! Vertex father id */
	  int     & father() { return m_father; };
protected:
   /*! Wheter the vertex has been accessed in the integration process m_touched */
    bool        m_touched;
	 /*! Vertex texture coordinates */
    CPoint2		m_uv;
	/*! Vertex father id */
	int         m_father;
  
  public:
	/*! CIntegrationVertex constructor */
    CIntegrationVertex() { m_father = 0; m_touched = false;};
	/*! CIntegrationVertex destructor */
    ~CIntegrationVertex(){};
	/*! Save vertex uv coordinates to the vertex string */
	void _to_string();
	/*! Get vertex father from the vertex string */
	void _from_string();

 };

// Get vertex father from the vertex string 

inline void CIntegrationVertex::_from_string( )
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
	  }
};

// Save vertex uv coordinates to the vertex string 
inline void CIntegrationVertex::_to_string()
{
	CParser parser( m_string );
	parser._removeToken( "uv" );

	parser._toString( m_string );

	std::string line;
	std::stringstream iss( line );

	iss << "uv=(" << m_uv[0] << " " << m_uv[1] << ")";

	 if( m_string.size() > 0 )
    {
		  m_string += " ";
    }
    m_string += iss.str();
};



/*---------------------------------------------------------------------------------------------------------------------------------------

	Integration Edge Trait

----------------------------------------------------------------------------------------------------------------------------------------*/
/*! \brief CIntegrationEdge class
 *
 *	Edge class for computing integration
 *  The edge holomorphic 1-form trait is added : m_duv
 */
class CIntegrationEdge : public  CEdge
{
  public:
	/*!  The edge holomorphic 1-form trait */
	  CPoint2  & duv() { return m_duv; };
protected:
	/*!  The edge holomorphic 1-form trait */
    CPoint2		m_duv;

 public:
	/*! Read edge holomorphic 1-form m_duv from the string with the key token "duv" */
	void _from_string();
	/*! write edge holomorphic 1-form m_duv to the string with the key token "duv" */
	void _to_string();
};

// Read edge holomorphic 1-form m_duv from the string with the key token "duv" 
inline void CIntegrationEdge::_from_string()
{
	  CParser parser( m_string );
	
	  for( std::list<CToken*>::iterator iter = parser.tokens().begin() ; iter != parser.tokens().end(); ++ iter )
	  {
		  CToken * token = *iter;
		  if( token->m_key == "duv" )
		  {
			  token->m_value >> m_duv;
		  }
	  }
};

//write holomorphic 1-form trait m_duv to the string "duv"

inline void CIntegrationEdge::_to_string()
{
	CParser parser( m_string );
	parser._removeToken( "duv" );

	parser._toString( m_string );
	
	std::string line;
	std::stringstream iss(line);
	iss << "duv=(" << m_duv[0] << " " << m_duv[1] << ")";

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
/*! \brief CIntegrationMesh class
 *  Mesh class for integrating holomorphic 1-form on a mesh
*/
template<typename V, typename E, typename F, typename H>
class CIntegrationMesh : public CBaseMesh<V,E,F,H>
{
public:
	typedef V CVertex;
	typedef E CEdge;
	typedef F CFace;
	typedef H CHalfEdge;

	typedef MeshVertexIterator<V,E,F,H> MeshVertexIterator;
	typedef MeshEdgeIterator<V,E,F,H> MeshEdgeIterator;
	typedef VertexVertexIterator<V,E,F,H> VertexVertexIterator;
	typedef VertexEdgeIterator<V,E,F,H> VertexEdgeIterator;
public:
};


typedef CIntegrationMesh<CIntegrationVertex, CIntegrationEdge, CFace, CHalfEdge> CIMesh;

} //namespace Topology

unsigned long long Topology::CIMesh::m_input_traits  = VERTEX_FATHER | EDGE_DUV;
unsigned long long Topology::CIMesh::m_output_traits = VERTEX_UV;

} //namespace MeshLib

#endif  _INTEGRATION_TRAIT_H_