/*!
*      \file UnifiedRicciFlowMesh.h
*      \brief Mesh for Unified Ricci Flow
*	   \author David Gu
*      \date Documented 01/10/2014
*
*/

#ifndef  _GRAPH_EMBEDDING_RICCI_FLOW_MESH_H_
#define  _GRAPH_EMBEDDING_RICCI_FLOW_MESH_H_

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
namespace RicciFlow
{

/*!
*	\brief CGraphEmbeddingRicciFlowVertex class
*
*	Vertex class for Unified Ricci flow
*   Traits:
*   index
*	father
*   log_radius
*   curvature
*   vertex uv
*   target_curvature
*   touched
*/
class CGraphEmbeddingRicciFlowVertex : public CVertex
{

public:

	/*! 
	 *	Each bit in the traits indicates whether the vertex class has the correpsonding trait.
	 *  e.g. if ( traits & TRAIT_UV ), then vertex->huv() needs to be stored int the vertex string.
	 */
	static unsigned int traits;

	/*!
	 *	CUnifiedRicciFlowVertex constructor
	 */
	CGraphEmbeddingRicciFlowVertex() { m_index = 0; };
	/*!
	 *	CRicciFlowVertex destructor
	 */
	~CGraphEmbeddingRicciFlowVertex() {};

	/*!
	 *	Vertex uv trait
	 */
	CPoint2 & huv() { return m_huv; };
	/*!
	 *	Vertex index trait
	 */
	int &     idx() { return m_index; };
	/*!
	 *  Vertex father trait
	 */
	int &	  father(){ return m_father;};
	/*!
	 *	Vertex log radius
	 */
	double &   u()    { return m_log_radius; };
	/*!
	 *	Vertex curvature trait
	 */
	double &   k()    { return m_curvature; };
	/*!
	 *	Vertex target curvature
	 */
	double & target_k() { return m_target_curvature; };
	/*!
	 *	whether the vertex has been touched
	 */
	bool  & touched()  { return m_touched; };

	/*!
	 *	Read vertex uv from vertex string
	 */
	void _from_string();
	/*!
	 *	Write vertex uv to vertex string
	 */
	void _to_string();
	/*!
	 * Topological valence of the vertex
	 */
	int & valence() { return m_valence; };
	/*!
	 *	Vertex color
	 */
	//CPoint & rgb() { return m_rgb; };

	/*!
	 *	vertex radius
	 */
	double & r() { return m_r; };
	/*!
	 *	circle packing scheme
	 */
	double & e() { return m_e; };
	/*!
	 *  vertex b
	 */
	double & b() { return m_b; };
	/*!
	 *	type
	 */
	int    & type() { return m_type; };
protected:
	/*! Vertex uv */
	CPoint2 m_huv;
	/*! Vertex index */
	int     m_index;
	//father
	int     m_father;
	//log radius
	double  m_log_radius;
	//curvature
	double  m_curvature;
	//target curvature
	double  m_target_curvature;
	//touched
	bool    m_touched;
	//topological valence
	int     m_valence;
	//vertex radius
	double  m_r;
	//scheme
	double  m_e;
	//b value
	double m_b;
	//type
	int    m_type;
	//vertex color
	//CPoint  m_rgb;
};

//read father from string
inline void CGraphEmbeddingRicciFlowVertex::_from_string()
{
  CParser parser( m_string );

  for( std::list<CToken*>::iterator iter = parser.tokens().begin() ; iter != parser.tokens().end(); ++ iter )
  {
	  CToken * token = *iter;
	  if( token->m_key == "type" )
	  {
		 std::string line = strutil::trim( token->m_value, "()");
		 m_type = strutil::parseString<int>( line );	
		 continue;	
	  }
	  if( token->m_key == "K" )
	  {
		 std::string line = strutil::trim( token->m_value, "()");
		 m_target_curvature = strutil::parseString<double>( line );	
		 continue;
		
	  }

	  if( token->m_key == "r" )
	  {
		 std::string line = strutil::trim( token->m_value, "()");
		 m_r = strutil::parseString<double>( line );	
		 continue;
		
	  }
	  /*
	  if( token->m_key == "rgb" )
	  {
		 std::string line = strutil::trim( token->m_value, "()");
		 line >> m_rgb;
		 continue;
	  }
	  */
	  /*
	  if( token->m_key == "uv" )
	  {
		 std::string line = strutil::trim( token->m_value, "()");
		 line >> m_huv;		
		 continue;
	  }
	  */
  }
}

//converting vertex uv trait to string

inline void CGraphEmbeddingRicciFlowVertex::_to_string()
{
	CParser parser( m_string );
	parser._removeToken( "type" );
	parser._removeToken( "K" );

	if( traits & TRAIT_UV )
	{
		parser._removeToken( "uv" );
	}
	if( traits & TRAIT_RADIUS )
	{
		parser._removeToken( "r" );
	}
	//if( traits & TRAIT_RGB )
	//{
	//	parser._removeToken( "rgb" );
	//}
	parser._toString( m_string );
	std::stringstream iss;
	iss << "type=(" << m_type <<") ";
	iss << "K=(" << m_target_curvature << ") ";

	if( traits & TRAIT_UV )
	{
		iss.precision(10);
		iss << "uv=(" << m_huv[0] << " " << m_huv[1] << ") ";
	}
	//if( traits & TRAIT_RGB )
	//{
	//	iss << "rgb=(" << m_rgb[0] << " " << m_rgb[1] << " " << m_rgb[2] << ") ";
	//}
	if( traits & TRAIT_RADIUS )
	{
		iss.precision(10);
		double r = exp( m_log_radius );
		iss << "r=(" << r << ") ";
	}

	if( m_string.length() > 0 )
	{
		m_string += " ";
	}
	m_string += iss.str();


}

/*!
*	\brief CGraphEmbeddingRicciFlowEdge class
*
*	Edge class for Ricci flow
*   trait:
*   edge length
*   edge weight \$f \frac{\partial \theta_i}{\partial u_j} = w_k \f$
*   edge inversive distance \f$ \cos \phi \f$, intersection angle is \f$\phi\f$
*/
class CGraphEmbeddingRicciFlowEdge : public  CEdge
  {
  public:
    /*!	CUnifiedRicciFlowEdge constructor
	 */
	 CGraphEmbeddingRicciFlowEdge() { m_length = 0; m_weight = 0; m_sharp = false; };
    /*!	CRicciFlowEdge destructor
	 */
    ~CGraphEmbeddingRicciFlowEdge(){};
	/*!	edge weight trait
	 */
	double & weight() { return m_weight; };
	/*!	edge length trait
	 */
	double & length() { return m_length; };
	/*!	initial edge length trait
	 */
	double & initial_length() { return m_initial_length; };
	/*! edge inversive distance
	*/
	double & inversive_distance() { return m_inversive_distance; };
	/*! edge sharp
	 */
	bool sharp() { return m_sharp; };
	/*!	eta
	 */
	double & eta() { return m_eta; };

    /*!
	 *	read sharp trait from the edge string
	 */
	void _from_string() 
	{ 
		CEdge::_from_string(); 
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
	/*!
	 *	write no traits to edge string
	 */
	void _to_string() { CEdge::_to_string(); };
	int & idx() { return m_index; };

  protected:
	  /*! edge weight trait */
	double   m_weight;
	  //edge length
	double   m_length;
	//initial edge length
	double   m_initial_length;
	//inversive distance
	double   m_inversive_distance;
	/*!
	 *	if edge is sharp
	 */
	bool     m_sharp;
	/*! edge eta - conformal structure coefficient
	 */
	double   m_eta;
	/*! index
	 */
	int      m_index;
};

/*!
*	\brief CGraphEmbeddingRicciFlowHalfEdge class
*
*	HalfEdge class for Ricci flow
*   trait:
*   corner angle
*/

class CGraphEmbeddingRicciFlowHalfEdge : public  CHalfEdge
{
public:
  CGraphEmbeddingRicciFlowHalfEdge() { m_angle = 0.0; };
  ~CGraphEmbeddingRicciFlowHalfEdge() {};
  /*! corner angle */
  double & angle() { return m_angle; };
  /*! on half edge [v_i,v_j], \f$ \frac{\parital \theat_i}{\partial u_j} */		
  double & dtheta_du() { return m_theta_u; };
  /*! on half edge [v_i,v_j], \f$ \frac{\parital l_{ij}}{\partial u_j} */		
  double & dl_du()     { return m_l_u; };

  double & c_w_ii()    { return m_c_w_ii; };
  double & c_w_ij()    { return m_c_w_ij; };

protected:
  double m_angle;
  CPoint m_s;
  double m_theta_u;
  double m_l_u;
  double m_c_w_ii;
  double m_c_w_ij;

};



/*!
*	\brief CGraphEmbeddingRicciFlowFace class
*
*	Face class for Ricci flow
*   trait:
*	whether the face has been touched
*/

class CGraphEmbeddingRicciFlowFace: public CFace
{
public:
	CGraphEmbeddingRicciFlowFace() { m_touched = false; };
	~CGraphEmbeddingRicciFlowFace() {};
	bool & touched() { return m_touched; };
	CPoint & normal(){ return m_normal; };
	/*!
	 *	Write vertex uv to vertex string
	 */
	void _to_string();
	int & idx() { return m_index; };

protected:
	int  m_index;
	bool m_touched;
	CPoint m_normal;
};


inline void CGraphEmbeddingRicciFlowFace::_to_string()
{
	CParser parser( m_string );
	parser._removeToken( "type" );
	parser._toString( m_string );
	std::stringstream iss;
	int type = (m_touched)?1:0;
	iss << "type=(" << type <<") ";
	if( m_string.length() > 0 )
	{
		m_string += " ";
	}
	m_string += iss.str();


}


/*-------------------------------------------------------------------------------------------------------------------------------------

	Ricci flow Mesh Class

--------------------------------------------------------------------------------------------------------------------------------------*/
/*!
 *	\brief CUnifiedRicciFlowMesh class
 *
 *	Mesh class for unified Ricci flow
 */
template<typename V, typename E, typename F, typename H>
class CGraphEmbeddingRicciFlowMesh : public CBaseMesh<V,E,F,H>
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
	typedef VertexFaceIterator<V,E,F,H> VertexFaceIterator;
	typedef VertexInHalfedgeIterator<V,E,F,H> VertexInHalfedgeIterator;
	typedef VertexOutHalfedgeIterator<V,E,F,H> VertexOutHalfedgeIterator;
	typedef FaceHalfedgeIterator<V,E,F,H> FaceHalfedgeIterator;
	typedef FaceEdgeIterator<V,E,F,H> FaceEdgeIterator;
	typedef MeshFaceIterator<V,E,F,H> MeshFaceIterator;
	typedef FaceVertexIterator<V,E,F,H> FaceVertexIterator;
	typedef MeshHalfEdgeIterator<V,E,F,H> MeshHalfEdgeIterator;
public:
};


typedef CGraphEmbeddingRicciFlowMesh<CGraphEmbeddingRicciFlowVertex, 
CGraphEmbeddingRicciFlowEdge, 
CGraphEmbeddingRicciFlowFace, 
CGraphEmbeddingRicciFlowHalfEdge> CGERFMesh;	

};	//namespace RicciFlow

};	//namespace MeshLib

#endif  _GRAPH_EMBEDDING_RICCI_FLOW_MESH_H_