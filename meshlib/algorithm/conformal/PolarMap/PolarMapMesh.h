/*!
*      \file PolarMapMesh.h
*      \brief Mesh for Computing Exponential map \f$ z \to e^z \f$
*	   \author David Gu
*      \date Documented on 10/12/2010
*
*/
/*******************************************************************************
*      Polar Map Mesh
*
*       Copyright (c) Stony Brook University
*
*    Purpose:
*
*       Polar Map Mesh
* 
*       David Gu June 27, 2008,  gu@cs.stonybrook.edu
*
*******************************************************************************/


#ifndef  _POLAR_MAP_MESH_H_
#define  _POLAR_MAP_MESH_H_

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

namespace Holomorphy
{

/*---------------------------------------------------------------------------------------------------------------------------------------

	Polar Vertex Trait

----------------------------------------------------------------------------------------------------------------------------------------*/
/*! \brief CPolarMapVertex class
 *
 *   Vertex class for computing exponential map
 *   Traits:
 *   vertex texture coordinates m_uv
 *   vertex father id  m_father
 */

class CPolarMapVertex : public  CVertex
{

  public:
	/*!   vertex texture coordinates m_uv */
	CPoint2 & uv()    { return m_uv; };
	/*!   vertex father id  m_father */
	int     & father() { return m_father; };

protected:
	/*! vertex texture coordinates m_uv */
    CPoint2		m_uv;
	/*! vertex father id */
	int         m_father;
  
  public:
	/*! CPolarMapVertex constructor */
    CPolarMapVertex() { m_father = 0;};
	/*! CPolarMapVertex destructor */
    ~CPolarMapVertex(){};
	
	/*! write vertex uv to vertex string*/
	void _to_string();
	/*! read vertex uv from vertex string*/
	void _from_string();

 };

//read in the vertex father information

inline void CPolarMapVertex::_from_string( )
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
		  if( token->m_key == "uv" )
		  {
			  token->m_value >> m_uv;
		  }
	  }
};

//write vertex uv coordinates to the string

inline void CPolarMapVertex::_to_string()
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

	PolarMap Mesh

----------------------------------------------------------------------------------------------------------------------------------------*/

/*! \brief CPolarMapMesh class
*
*	Mesh class for computing exponential map \f$z \to e^{z}\f$
*/
template<typename V, typename E, typename F, typename H>
class CPolarMapMesh : public CBaseMesh<V,E,F,H>
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


typedef CPolarMapMesh<CPolarMapVertex, CEdge, CFace, CHalfEdge> CPMMesh;

} //namespace Holomorphy

unsigned long long Holomorphy::CPMMesh::m_input_traits  = VERTEX_FATHER | VERTEX_UV;
unsigned long long Holomorphy::CPMMesh::m_output_traits = VERTEX_UV;

} //namespace MeshLib

#endif  _POLAR_MAP_MESH_H_