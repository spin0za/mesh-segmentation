/*!
*      \file StitchMesh.h
*      \brief Mesh for Stiching two meshes
*	   \author David Gu
*      \date Documented on 10/12/2010
*
*/
#ifndef  _STITCH_MESH_H_
#define  _STITCH_MESH_H_

#include <map>
#include <vector>

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

/*--------------------------------------------------------------------------------------------------------------------------------------

	Stitch Face Trait

---------------------------------------------------------------------------------------------------------------------------------------*/
/*!	\brief CStitchFace class
 *
 *	Face class for stitching two meshes together
 *
 *	Traits
 *	segment id for the current face m_segment_id
 *  face index for the current face m_idx
 */
  class CStitchFace: public CFace
  {
  protected:
  
	  /*!	segment id for the current face m_segment_id */
	  int m_segment_id;
	  /*!   face index */
	  int m_idx;
	
  public:

	  /*!	segment id for the current face m_segment_id */
	  int & segment() { return m_segment_id; };
	  /*!   face index */
	  int & idx()     { return m_idx;        };
	
  public:
	/*! read segment id from the string */
	void _from_string( );
	/*! write segment id to the string */
	void _to_string( );
	/* CStitchFace constructor */
    CStitchFace() { m_segment_id = 0; m_idx = 0; };
	/* CStitchFace destructor  */
    ~CStitchFace() {};
  };


  /*! read segment id from the string */
  inline void CStitchFace::_from_string()
	{
		  CParser parser( m_string );
		
		  for( std::list<CToken*>::iterator iter = parser.tokens().begin() ; iter != parser.tokens().end(); ++ iter )
		  {
			  CToken * token = *iter;
			  if( token->m_key == "sid" )
			  {
				  std::string line = strutil::trim( token->m_value, "()");
				  m_segment_id = strutil::parseString<int>( line );
			  }
		  }
	 };
//write face traits to string
  inline void CStitchFace::_to_string()
  {
	  CParser parser( m_string );
	  parser._removeToken( "sid" );
	  parser._toString( m_string );

	  std::stringstream iss;
	  iss << "sid=(" << m_segment_id << ")";

	  if( m_string.length() > 0 )
	  {
		  m_string += " ";
	  }
	  m_string += iss.str();
   };

/*--------------------------------------------------------------------------------------------------------------------------------------

	Stitch Vertex Trait

---------------------------------------------------------------------------------------------------------------------------------------*/
 /*! \brief CStitchVertex class
  *
  *  Vertex class for stitching two meshes
  *  Traits:
  *  vertex index : m_index
  *  vertex father: m_father
  *  vertex uv coordinates: m_uv
  */
  class CStitchVertex : public  CVertex
  {
  protected:
	/*!  vertex index */
    int         m_index;
	/*! vertex father */
	int         m_father;
	/*! vertex texture uv coordinates */
	CPoint2     m_uv;
	/*! vertex rgb color */
	CPoint      m_rgb;
	/*! if the vertex is touched */
	bool        m_touched;

  public:
	/*!  vertex index */
	  int & idx()   { return m_index; };
	/*! vertex father */
	  int &father() { return m_father; };
	/*! vertex texture uv coordinates */
	  CPoint2 & uv(){ return m_uv; };
	/*!   rgb  */
	  CPoint &rgb() { return m_rgb;};
	/*! touched */
	  bool &touched() { return m_touched; };

  public:
	  /*! CStitchVertex Constructor */
    CStitchVertex() { m_index = 0; m_father =0;};
	/*! CStitchVertex destructor */
    ~CStitchVertex() {};
	/*! read father from vertex string */
	void _from_string();
	/*! write father and uv to vertex string */
	void _to_string();

  };

 //write vertex traits to string
  inline void CStitchVertex::_to_string()
  {
	  CParser parser( m_string );
	  parser._removeToken( "father" );
	  parser._removeToken( "rgb" );
	  parser._toString( m_string );

	  std::stringstream iss;
	  iss << "father=(" << m_father << ") ";
	  iss << "rgb=(" << m_rgb[0]<<" "<< m_rgb[1]<<" "<< m_rgb[2] << ")";

	  if( m_string.length() > 0 )
	  {
		  m_string += " ";
	  }


	  m_string += iss.str();
   };

//read vertex traits from string 
  inline void CStitchVertex::_from_string()
 {
		  CParser parser( m_string );
		
		  for( std::list<CToken*>::iterator iter = parser.tokens().begin() ; iter != parser.tokens().end(); ++ iter )
		  {
			  CToken * token = *iter;
			  if( token->m_key == "father" )
			  {
				  std::string line = strutil::trim( token->m_value, "()");
				  m_father = strutil::parseString<int>(line);
			  }
			  if( token->m_key == "uv" )
			  {
				  token->m_value >> m_uv;
			  }
			  if( token->m_key == "rgb" )
			  {
				 std::string line = strutil::trim( token->m_value, "()");
				 line >> m_rgb;		
			  }
		  }
	  };




/*--------------------------------------------------------------------------------------------------------------------------------------

	Stitch Edge Trait

---------------------------------------------------------------------------------------------------------------------------------------*/
/*! \brief CStichEdge class
* 
*	Edge class for stitching two meshes 
*   Traits:
*   whether the edge is sharp : m_sharp
*/
class CStitchEdge : public  CEdge
  {
  protected:
	  /*! whether the edge is sharp */
    bool    m_sharp;
 
  public:
	  /*! CStitchEdge constructor */
    CStitchEdge() { m_sharp = false;};
	/*! CStitchEdge destructor */
    ~CStitchEdge(){};
	/*! read sharp trait from edge string */
	void _from_string();
	/*! write sharp trait from edge string */
	void _to_string();
	/*! whether the edge is labeled as sharp */
	bool & sharp() { return m_sharp; };
};

/*! read sharp trait from edge string */
inline void CStitchEdge::_from_string()
{
	CParser parser( m_string );

	for( std::list<CToken*>::iterator iter = parser.tokens().begin() ; iter != parser.tokens().end(); ++ iter )
	{
			CToken * token = *iter;

			if( token->m_key == "sharp" )
			{
				m_sharp = true;		
			}
	};
};

/*! write sharp trait from edge string */
inline void CStitchEdge::_to_string()
{
	if( !m_sharp ) return;

	if( m_string.length() > 0 )
	{
		m_string += " ";
	}
	m_string += "sharp";
};

/*!\brief CStitchHalfEdge class
*
*	HalfEdge class for stiching two meshes
*   trait: corner angle
*/
class CStitchHalfEdge : public  CHalfEdge
  {
  public:
	  /*! corner angle */
	double &   angle() { return m_angle; };

  public:
	/*! CStitchHalfEdge constructor */
	 CStitchHalfEdge() { m_angle = 0; };
	/*! CStitchHalfEdge destructor */
    ~CStitchHalfEdge(){};
	/*! read corner angle from the halfedge string */
	void _from_string();
	/*! write corner angle to the halfedge string */
	void _to_string(); 

  //no output
  protected: 
	  double   m_angle;

};

// read corner angle from the halfedge string 
inline void CStitchHalfEdge::_from_string()
{
	CParser parser( m_string );

	for( std::list<CToken*>::iterator titer = parser.tokens().begin(); titer != parser.tokens().end(); titer ++ )
	{
		CToken* pT = *titer;
		if( pT->m_key == "a" )
		{
			std::string line = strutil::trim( pT->m_value,"()");
			m_angle = strutil::parseString<double>( line );
		}
	}
};

//write corner angle to the halfedge string 
inline void CStitchHalfEdge::_to_string()
{
	CParser parser( m_string );
	parser._removeToken( "a" );

	parser._toString( m_string );
	
	std::string line;
	std::stringstream iss;
	iss << "a=(" << m_angle << ")";

	if( m_string.length() > 0 )
	{
		m_string += " ";
	}
	m_string += iss.str();
};

/*! \brief CStitchMesh class
*
*	Mesh class for stitching two meshes
*/
template<typename V, typename E, typename F, typename H>
class CStitchMesh : public CBaseMesh<V,E,F,H>
{
public:
	typedef V CVertex;
	typedef E CEdge;
	typedef F CFace;
	typedef H CHalfEdge;

	typedef CBoundary<V,E,F,H>            CBoundary;
	typedef CLoop<V,E,F,H>                CLoop;
	typedef MeshVertexIterator<V,E,F,H>   MeshVertexIterator;
	typedef MeshFaceIterator<V,E,F,H>	  MeshFaceIterator;
	typedef MeshEdgeIterator<V,E,F,H>     MeshEdgeIterator;
	typedef VertexVertexIterator<V,E,F,H> VertexVertexIterator;
	typedef FaceHalfedgeIterator<V,E,F,H> FaceHalfedgeIterator;
	typedef FaceVertexIterator<V,E,F,H>   FaceVertexIterator;
};


typedef CStitchMesh<CStitchVertex, CStitchEdge, CStitchFace, CStitchHalfEdge> CSTMesh;

} //namespace Topology

unsigned long long Topology::CSTMesh::m_input_traits = VERTEX_FATHER;
unsigned long long Topology::CSTMesh::m_output_traits = VERTEX_FATHER;

} //namespace MeshLib
#endif  _STITCH_MESH_H_