/*!
*      \file SlitMapMesh.h
*      \brief Mesh for Computing Slit mappings
*	   \author David Gu
*      \date Documented on 10/12/2010
*
*/

#ifndef  _SLIT_MAP_MESH_H_
#define  _SLIT_MAP_MESH_H_


#include <map>
#include <vector>
#include <list>

#include "Mesh/Vertex.h"
#include "Mesh/Edge.h"
#include "Mesh/Face.h"
#include "Mesh/HalfEdge.h"

#include "Mesh/BaseMesh.h"
#include "Mesh/boundary.h"
#include "Mesh/iterators.h"
#include "Parser/parser.h"
#include "Parser/traits_io.h"
#include "Operator/Operator.h"


namespace MeshLib
{

namespace Holomorphy
{

/*-------------------------------------------------------------------------------------------------------------------------------------

	Slit Map Edge Trait Class

--------------------------------------------------------------------------------------------------------------------------------------*/
/*! \brief CSlitMapEdge class
* 
*	Edge class for computing Slit maps
*   Triat:
*   holomorphic 1-form m_duv
*/

class CSlitMapEdge : public  CEdge
  {
  protected:
	/*!   holomorphic 1-form defined on edge */
	CPoint2  m_duv;
  public:
	/*!   holomorphic 1-form defined on edge */
	  CPoint2 & duv() { return m_duv; };

  public:
	  /*! CSlitMapEdge constructor */
    CSlitMapEdge() {};
	  /*! CSlitMapEdge destructor */
    ~CSlitMapEdge(){};

  public:
	  /*! read holomorphic form from edge string*/
	void _from_string();
	  /*! write holomorphic form to edge string*/
	void _to_string();
  };

// read holomorphic form from edge string
inline void CSlitMapEdge::_from_string()
{
	  CParser parser( m_string );

		std::list<CToken*> & tokens = parser.tokens();
		for( std::list<CToken*>::iterator titer = tokens.begin(); titer != tokens.end(); titer ++ )
		{
			CToken * pT = *titer;
			if( pT->m_key == "duv" )
			{
					pT->m_value >> m_duv;
			}
		}
};

// write holomorphic form to edge string
inline void CSlitMapEdge::_to_string()
{
	CParser parser( m_string );
	parser._removeToken( "duv" );
	parser._toString( m_string );

	std::stringstream iss;
	iss << "duv=("<< m_duv[0] << " " << m_duv[1] << ")";

	if( m_string.length() > 0 ) m_string += " ";
	m_string += iss.str();
};


/*-------------------------------------------------------------------------------------------------------------------------------------

	Slit Map Mesh Class

--------------------------------------------------------------------------------------------------------------------------------------*/
/*! \brief CSplitMapMesh class
*
*	Mesh for computing slit maps
*/

template<typename V, typename E, typename F, typename H>
class CSlitMapMesh : public CBaseMesh<V,E,F,H>
{
public:
	typedef V CVertex;
	typedef E CEdge;
	typedef F CFace;
	typedef H CHalfEdge;

	typedef CBoundary<V,E,F,H> CBoundary;
	typedef CLoop<V,E,F,H> CLoop;
	typedef MeshVertexIterator<V,E,F,H> MeshVertexIterator;
	typedef MeshEdgeIterator<V,E,F,H> MeshEdgeIterator;
	typedef VertexVertexIterator<V,E,F,H> VertexVertexIterator;
	typedef VertexEdgeIterator<V,E,F,H> VertexEdgeIterator;
public:
};


typedef CSlitMapMesh<CVertex, CSlitMapEdge, CFace, CHalfEdge> CSMMesh;

} //namespace Holomorphy

unsigned long long Holomorphy::CSMMesh::m_input_traits  = EDGE_DUV;
unsigned long long Holomorphy::CSMMesh::m_output_traits = EDGE_DUV;

} //namespace MeshLib

#endif  _SLIT_MAP_MESH_H_