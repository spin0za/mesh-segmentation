/*! \file DijkstraMesh.h
*  \brief Mesh class for Computing Shortest Paths
*  \author David Gu
*  \data Documented on 03/20/2011
*/


#ifndef  _DIJKSTRA_MESH_H_
#define  _DIJKSTRA_MESH_H_

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

/*!
* \brief CDijkstraVertex class
*
*  Vertex class for computing shortest path
*  The following traits are added
*  Vertex distance to the source
*  Wheter the vertex has been accessed 
*  The father vertex id m_father
*/
class CDijkstraVertex : public  CVertex
{

  public:
	  /*! Wheter the vertex has been accessed in the integration process m_touched */
	  bool & touched()  { return m_touched; };
	  /*! Vertex distance */
	  double & distance()    { return m_distance; };
	  /*! Vertex father id */
	  int     & father() { return m_father; };
	  /*! Vertex parent id */
	  int     & parent() { return m_parent; };
	  /*! Vertex ancesotr id */
	  int     & ancestor() { return m_ancestor; };


	  /*! Compare two DijkstraVertex */

	  const bool operator <(const CDijkstraVertex & pV) const 
	  {     
		  return (m_distance < pV.m_distance); 
	  }
	  /*! Previous vertex */
	  CDijkstraVertex * & previous() { return m_previous; };
	  
	  /*! root vertex */
	  CDijkstraVertex * & root() { return m_root; };

	  /*! UV coordinates */
	  CPoint2  & uv() { return m_uv; };

protected:
   /*! Wheter the vertex has been accessed in the integration process m_touched */
    bool        m_touched;
	 /*! Vertex distance */
    double m_distance;
	/*! Vertex parent id */
	int         m_parent;
	/*! Vertex ancestor id */
	int         m_ancestor;
	/*! Vertex father id */
	int         m_father;

	/*! Previous vertex */
	CDijkstraVertex * m_previous;
	/*! Root vertex */
	CDijkstraVertex * m_root;
	/*! uv coordinates */
	CPoint2           m_uv;

  public:
	/*! CIntegrationVertex constructor */
    CDijkstraVertex() { m_parent = 0; m_touched = false; m_distance = 1e+20;};
	/*! CIntegrationVertex destructor */
    ~CDijkstraVertex(){};

	/*! Get vertex father from the vertex string */
	void _from_string();

 };

/*! Compuare two Dijkstra Vertices according to their distance traits */

class CompareVertex 
{
public:
    bool operator()(CDijkstraVertex * pV1, CDijkstraVertex * pV2)
    {
       return (pV1->distance() > pV2->distance() );
    }
};


/*! Compuare two Dijkstra Vertices according to their father traits */

class CompareVertexParent 
{
public:
    bool operator()(CDijkstraVertex * pV1, CDijkstraVertex * pV2)
    {
       return (pV1->parent() > pV2->parent() );
    }
};

class CompareVertexFather 
{
public:
    bool operator()(CDijkstraVertex * pV1, CDijkstraVertex * pV2)
    {
       return (pV1->father() > pV2->father() );
    }
};

// Get vertex father from the vertex string 

inline void CDijkstraVertex::_from_string( )
{
	  CParser parser( m_string );

	  for( std::list<CToken*>::iterator iter = parser.tokens().begin() ; iter != parser.tokens().end(); ++ iter )
	  {
		  CToken * token = *iter;
		  if( token->m_key == "parent" )
		  {
			  std::string line = strutil::trim( token->m_value, "()");
			  m_parent = strutil::parseString<int>( line );	
		  }
		  else if( token->m_key == "father" )
		  {
			  std::string line = strutil::trim( token->m_value, "()");
			  m_father = strutil::parseString<int>( line );	
		  }
		  else if( token->m_key == "ancestor" )
		  {
			  std::string line = strutil::trim( token->m_value, "()");
			  m_ancestor = strutil::parseString<int>( line );	
		  }
	  }
};


/*! \brief CDijkstraEdge class
 *
 *	Edge class for computing the shortest path
 *  The edge trait is added : m_sharp
 */
class CDijkstraEdge : public  CEdge
{
  public:
	  /*! CDijkstraEdge constructor */
	  CDijkstraEdge() { m_length = -1; m_sharp = false; };
	  /*! CDijkstraEdge destructor */
	  ~CDijkstraEdge() {};
	  /*!  The edge sharp trait */
	  bool  & sharp() { return m_sharp; };
      /*!  Edge length */
	  double & length(){ return m_length;};

protected:
	/*!  The edge sharp trait */
    bool		m_sharp;
	/*! Edge length */
	double m_length;
 public:
	/*! Read edge sharp from the string with the key token "sharp" */
	void _from_string();
	/*! write edge sharp to the string with the key token "sharp" */
	void _to_string();
};

// Read edge sharp from the string with the key token "duv" 
inline void CDijkstraEdge::_from_string()
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
		  if( token->m_key == "l" )
		  {
			  std::string t = token->m_value;
			  t.erase(0, t.find_first_not_of("()") );
		      t.erase(t.find_last_not_of("()") + 1);
			  
			  m_length = strutil::parseString<double>( t );
		  }
	  }
};

//write sharp trait m_duv to the string "duv"

inline void CDijkstraEdge::_to_string()
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

/*! Compuare two Dijkstra Edges according to the loop length */

class CompareEdge 
{
public:
    bool operator()(CDijkstraEdge * pE1, CDijkstraEdge * pE2)
    {
		CDijkstraVertex * pV1 = (CDijkstraVertex *)(pE1->halfedge(0)->target());
		CDijkstraVertex * pV2 = (CDijkstraVertex *)(pE1->halfedge(0)->source());
		double d1 = pV1->distance() + pV2->distance() + pE1->length();

		pV1 = (CDijkstraVertex *)(pE2->halfedge(0)->target());
		pV2 = (CDijkstraVertex *)(pE2->halfedge(0)->source());
		double d2 = pV1->distance() + pV2->distance() + pE2->length();

		return ( d1 > d2 );
    }
};


/*! \brief CDijkstraFace class
 *
 *	Face class for computing the shortest path
 *  The face trait is added : m_topo_valence
 */
class CDijkstraFace : public  CFace
{
  public:
    /*!  Topological valence */
	  int & _topo_valence(){ return m_topo_valence;};
	/*!  The segment id */
	  int & _id() { return m_id; };
	/*!  If the face has been touched */
	  bool& touched() { return m_touched; };

protected:
	/*!  The topological valence */
    int		m_topo_valence;
	/*!  The segment id */
	int     m_id;
	/*!  If the face has been touched */
	bool    m_touched;
};


/*! \brief CDijkstraMesh class
 *  Mesh class for computing the shortest path
*/
template<typename V, typename E, typename F, typename H>
class CDijkstraMesh : public CBaseMesh<V,E,F,H>
{
public:
	typedef CBoundary<V,E,F,H> CBoundary;
	typedef CLoop<V,E,F,H> CLoop;
	typedef MeshVertexIterator<V,E,F,H> MeshVertexIterator;
	typedef MeshFaceIterator<V,E,F,H> MeshFaceIterator;
	typedef MeshEdgeIterator<V,E,F,H> MeshEdgeIterator;
	typedef VertexVertexIterator<V,E,F,H> VertexVertexIterator;
	typedef VertexFaceIterator<V,E,F,H> VertexFaceIterator;
	typedef VertexEdgeIterator<V,E,F,H> VertexEdgeIterator;
	typedef VertexOutHalfedgeIterator<V,E,F,H> VertexOutHalfedgeIterator;
	typedef FaceEdgeIterator<V,E,F,H> FaceEdgeIterator;
	typedef FaceHalfedgeIterator<V,E,F,H> FaceHalfedgeIterator;
	typedef FaceVertexIterator<V,E,F,H> FaceVertexIterator;
public:
};


typedef CDijkstraMesh<CDijkstraVertex, CDijkstraEdge, CDijkstraFace, CHalfEdge> CDKMesh;


};

#endif  _DIJKSTRA_MESH_H_