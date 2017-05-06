/*! \file CanonicalBasisMesh.h
*  \brief Mesh class for Computing Canonical Basis for fundamental group
*  \author David Gu
*  \data Documented on 03/20/2011
*/


#ifndef  _CANONICAL_BASIS_MESH_H_
#define  _CANONICAL_BASIS_MESH_H_

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

namespace MeshLib
{

namespace Topology
{
/*!
* \brief CCanonicalBasisVertex class
*
*/
class CCanonicalBasisVertex : public  CVertex
{

  public:
	  /*! Vertex ancesotr id */
	  int     & ancestor() { return m_ancestor; };

protected:
	/*! Vertex ancestor id */
	int         m_ancestor;

  public:
	/*! CIntegrationVertex constructor */
    CCanonicalBasisVertex() { m_ancestor = 0; };
	/*! CIntegrationVertex destructor */
    ~CCanonicalBasisVertex(){};

	/*! Get vertex father from the vertex string */
	void _from_string();

 };

/*! Compuare two Dijkstra Vertices according to their father traits */

class CompareVertexAncestor 
{
public:
    bool operator()(CCanonicalBasisVertex * pV1, CCanonicalBasisVertex * pV2)
    {
       return (pV1->ancestor() > pV2->ancestor() );
    }
};

// Get vertex ancestor from the vertex string 

inline void CCanonicalBasisVertex::_from_string( )
{
	  CParser parser( m_string );

	  for( std::list<CToken*>::iterator iter = parser.tokens().begin() ; iter != parser.tokens().end(); ++ iter )
	  {
		  CToken * token = *iter;
		  if( token->m_key == "ancestor" )
		  {
			  std::string line = strutil::trim( token->m_value, "()");
			  m_ancestor = strutil::parseString<int>( line );	
		  }
	  }
};

/*! \brief CCanonicalBasisEdge class
 *
 *	Edge class for computing the shortest path
 *  The edge trait is added : m_sharp
 */
class CCanonicalBasisEdge : public  CEdge
{
  public:
	  /*! CCanonicalBasisEdge constructor */
	  CCanonicalBasisEdge() { m_sharp = false; };
	  /*! CDijkstraEdge destructor */
	  ~CCanonicalBasisEdge() {};
	  /*!  The edge sharp trait */
	  bool  & sharp() { return m_sharp; };
      static void * m_pMesh;

protected:
	/*!  The edge sharp trait */
    bool		m_sharp;

public:
	/*! Read edge sharp from the string with the key token "sharp" */
	void _from_string();
	/*! write edge sharp to the string with the key token "sharp" */
	void _to_string();
};

// Read edge sharp from the string with the key token "duv" 
inline void CCanonicalBasisEdge::_from_string()
{
	  CParser parser( m_string );
	
	  for( std::list<CToken*>::iterator iter = parser.tokens().begin() ; iter != parser.tokens().end(); ++ iter )
	  {
		  CToken * token = *iter;
		  if( token->m_key == "sharp" )
		  {
			  m_sharp = true;
			  continue;
		  }
	  }
};

//write sharp trait m_duv to the string "duv"

inline void CCanonicalBasisEdge::_to_string()
{
	CParser parser( m_string );
	parser._removeToken( "sharp" );
	parser._toString( m_string );

	if( !m_sharp ) return;

	
	
	std::string line;
	std::stringstream iss(line);
	iss << "sharp";

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





/*! \brief CDijkstraFace class
 *
 *	Face class for computing the shortest path
 *  The face trait is added : m_topo_valence
 */
class CCanonicalBasisFace : public  CFace
{
  public:

	/*!  If the face has been touched */
	  bool& touched() { return m_touched; };

protected:
	/*!  If the face has been touched */
	bool    m_touched;
};


/*! \brief CDijkstraMesh class
 *  Mesh class for computing the shortest path
*/
template<typename V, typename E, typename F, typename H>
class CCanonicalBasisMesh : public CBaseMesh<V,E,F,H>
{
public:
	typedef CBoundary<V,E,F,H> CBoundary;
	typedef CLoop<V,E,F,H> CLoop;
	typedef MeshVertexIterator<V,E,F,H> MeshVertexIterator;
	typedef MeshFaceIterator<V,E,F,H> MeshFaceIterator;
	typedef MeshEdgeIterator<V,E,F,H> MeshEdgeIterator;
	typedef VertexVertexIterator<V,E,F,H> VertexVertexIterator;
	typedef VertexEdgeIterator<V,E,F,H> VertexEdgeIterator;
	typedef VertexFaceIterator<V,E,F,H> VertexFaceIterator;
	typedef FaceEdgeIterator<V,E,F,H> FaceEdgeIterator;
	typedef FaceHalfedgeIterator<V,E,F,H> FaceHalfedgeIterator;
	typedef FaceVertexIterator<V,E,F,H> FaceVertexIterator;
public:
};

/*!
 *	Declare the static data member of CCanonicalBasisEdge::m_pMesh
 */
void * CCanonicalBasisEdge::m_pMesh = NULL;

typedef CCanonicalBasisMesh<CCanonicalBasisVertex, CCanonicalBasisEdge, CCanonicalBasisFace, CHalfEdge> CCBMesh;

} //namespace Topology
} //namespace MeshLib

#endif  _CANONICAL_BASIS_MESH_H_