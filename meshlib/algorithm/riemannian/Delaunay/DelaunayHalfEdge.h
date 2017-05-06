/*! \file DelaunayHalfEdge.h
 *  \brief Vertex, Edge classes for Delaunay Algrithm
 *  \date Documented on 10/14/2010
 */

#ifndef  _DELAUNAY_HALFEDGE_H_
#define  _DELAUNAY_HALFEDGE_H_

#include <map>
#include <vector>
#include <queue>
#include <list>

#include "Mesh/Vertex.h"
#include "Mesh/Edge.h"
#include "Mesh/Face.h"
#include "Mesh/HalfEdge.h"
#include "Parser/Parser.h"

namespace MeshLib
{

	/*!	\brief CDelaunayVertex class
	 *  
	 *	Vertex class for Delaunay triangulation
	 *	Traits, vertex uv coordinates m_uv
	 */
class CDelaunayVertex : public  CVertex
  {
  protected:
    CPoint2  m_uv;
	int      m_father;
	bool     m_touched;
  public:
	  
	  int     & father() { return m_father; };
	  CPoint2 & uv() { return m_uv; };
	  bool    & touched() { return m_touched; };
	 
	  void	  _from_string();
	  void    _to_string();
  };



 inline void CDelaunayVertex::_from_string()
 {
		  CParser parser( m_string );
		
		  for( std::list<CToken*>::iterator iter = parser.tokens().begin() ; iter != parser.tokens().end(); ++ iter )
		  {
			  CToken * token = *iter;
			  if( token->m_key == "uv" )
			  {
				  token->m_value >> m_uv;
			  }
		  }
  };

  inline void CDelaunayVertex::_to_string()
  {
		  CParser parser( m_string );
		  parser._removeToken( "uv" );
			
		  parser._toString( m_string );

		  std::stringstream iss;
		  iss << "uv=(" << m_uv[0] << " " << m_uv[1] << ")";

		  if( m_string.length() > 0 )
			  m_string += " ";
		  m_string += iss.str();
  };

	/*!	\brief CDelaunayEdge class
	 *  
	 *	Edge class for Delaunay triangulation
	 *	Traits, frozen if the current edge is frozen
	 */

 class CDelaunayEdge : public  CEdge
 {
 public:
	 CDelaunayEdge() { m_frozen = false; };
	 bool & frozen() { return m_frozen; };
	 void _to_string();
 protected:
	 bool m_frozen;
 };

 inline void CDelaunayEdge::_to_string()
 {
	 if( m_frozen ) m_string="sharp";
 }

 /*!	\brief CDelaunayHalfEdge class
	 *  
	 *	HalfEdge class for Delaunay triangulation
	 */
 class CDelaunayHalfEdge : public CHalfEdge
 {
 };

 	/*!	\brief CDelaunayFace class
	 *  
	 *	Face class for Delaunay triangulation
	 *	Traits, 
	 *  m_inside whether the current face is inside or outside
	 *  m_minAngle minimal angle
	 *  m_quality  quality of the triangle
	 */

class CDelaunayFace : public CFace
 {
 public:
	 
	 CDelaunayFace() { m_inside = true;  }; //by default all faces are inside

	 bool   & inside() {   return m_inside;  };
	 double & minAngle() { return m_minAngle; };
	 double & quality()  { return m_quality;  };
	 void    _to_string();
	 bool   & touched() { return m_touched; };
 protected:
	 bool m_touched;
	 bool m_inside;
	 double m_minAngle;
	 double m_quality;
 };

inline void CDelaunayFace::_to_string()
{
	if( !m_inside )
	{
		m_string = "rgb=(1 0 0)";
	}
};

//segment class
/*! \brief CSegment class
 *
 *  Segment in PLSG
 */
template<typename CVertex>
class CSegment
{
 public:
	  CSegment(){};
	 ~CSegment(){};
	  CVertex * & vertex( int id ) { return m_vertices[id]; };
 protected:
	  CVertex * m_vertices[2]; //v[0] source, v[1] target
 };


}
#endif