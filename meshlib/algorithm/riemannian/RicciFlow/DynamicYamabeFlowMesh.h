/*!
*      \file DynamicYamabeFlowMesh.h
*      \brief Mesh for Dynamic Yamabe Flow
*	   \author David Gu
*      \date Documented 10/16/2010
*
*/

#ifndef  _DYNAMIC_YAMABE_FLOW_MESH_H_
#define  _DYNAMIC_YAMABE_FLOW_MESH_H_

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

#define TRAIT NORMAL  1
#define TRAIT_FATHER  2
#define TRAIT_UV      4
#define TRAIT_RGB     8
#define TRAIT_PARENT  16

/*!
*	\brief CRicciFlowVertex class
*
*	Vertex class for Ricci flow
*   Traits:
*   index
*	father
*   log_radius
*   curvature
*   vertex uv
*   target_curvature
*   touched
*/
class CDynamicYamabeFlowVertex : public CVertex
{

public:

	/*! 
	 *	Each bit in the traits indicates whether the vertex class has the correpsonding trait.
	 *  e.g. if ( traits & TRAIT_UV ), then vertex->huv() needs to be stored int the vertex string.
	 */
	static unsigned int traits;

	/*!
	 *	CRicciFlowVertex constructor
	 */
	CDynamicYamabeFlowVertex() { m_index = 0; };
	/*!
	 *	CRicciFlowVertex destructor
	 */
	~CDynamicYamabeFlowVertex() {};

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

	CDynamicYamabeFlowVertex* & dual() { return m_dual; };

	CDynamicYamabeFlowVertex* & ancestor() { return m_ancestor; };

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
	//vertex color
	//CPoint  m_rgb;

	//dual vertex
	CDynamicYamabeFlowVertex * m_dual;
	//ancestor
	CDynamicYamabeFlowVertex * m_ancestor;
};

//read father from string
inline void CDynamicYamabeFlowVertex::_from_string()
{
  CParser parser( m_string );

  for( std::list<CToken*>::iterator iter = parser.tokens().begin() ; iter != parser.tokens().end(); ++ iter )
  {
	  CToken * token = *iter;
	  if( token->m_key == "father" )
	  {
		 std::string line = strutil::trim( token->m_value, "()");
		 m_father = strutil::parseString<int>( line );	
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

inline void CDynamicYamabeFlowVertex::_to_string()
{
	CParser parser( m_string );
	
	if( traits & TRAIT_UV )
	{
		parser._removeToken( "uv" );
	}
	//if( traits & TRAIT_RGB )
	//{
	//	parser._removeToken( "rgb" );
	//}
	parser._toString( m_string );
	std::stringstream iss;
	
	if( traits & TRAIT_UV )
	{
		iss << "uv=(" << m_huv[0] << " " << m_huv[1] << ") ";
	}
	//if( traits & TRAIT_RGB )
	//{
	//	iss << "rgb=(" << m_rgb[0] << " " << m_rgb[1] << " " << m_rgb[2] << ") ";
	//}
	if( m_string.length() > 0 )
	{
		m_string += " ";
	}
	m_string += iss.str();
}

/*!
*	\brief CDynamicYamabeFlowEdge class
*
*	Edge class for Dynamic Yamabe Flow
*   trait:
*   edge length
*   edge weight \$f \frac{\partial \theta_i}{\partial u_j} = w_k \f$
*   edge inversive distance \f$ \cos \phi \f$, intersection angle is \f$\phi\f$
*/
class CDynamicYamabeFlowEdge : public  CEdge
  {
  public:
    /*!	CRicciFlowEdge constructor
	 */
	 CDynamicYamabeFlowEdge() { m_length = 0; m_weight = 0; m_sharp = false; };
    /*!	CRicciFlowEdge destructor
	 */
    ~CDynamicYamabeFlowEdge(){};
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
};

/*!
*	\brief CDynamicYamabeFlowHalfEdge class
*
*	HalfEdge class for Ricci flow
*   trait:
*   corner angle
*/

class CDynamicYamabeFlowHalfEdge : public  CHalfEdge
{
public:
  CDynamicYamabeFlowHalfEdge() { m_angle = 0.0; };
  ~CDynamicYamabeFlowHalfEdge() {};
  /*! corner angle */
  double & angle() { return m_angle; };
  /*! on half edge [v_i,v_j], \f$ \frac{\parital \theat_i}{\partial u_j} */		
  double & dtheta_du() { return m_theta_u; };
  /*! on half edge [v_i,v_j], \f$ \frac{\parital l_{ij}}{\partial u_j} */		
  double & dl_du()     { return m_l_u; };

protected:
  double m_angle;
  CPoint m_s;
  double m_theta_u;
  double m_l_u;

};

/*!
*	\brief CDynamicYamabeFlowFace class
*
*	Face class for Dynamic Yamabe Flow
*   trait:
*	whether the face has been touched
*/

class CDynamicYamabeFlowFace: public CFace
{
public:
	CDynamicYamabeFlowFace() { m_touched = false; };
	~CDynamicYamabeFlowFace() {};
	bool & touched() { return m_touched; };
	CPoint & normal(){ return m_normal; };
protected:
	bool m_touched;
	CPoint m_normal;
};



/*
template<typename V, typename E, typename F, typename H>
class CDynamicYamabeFlowMesh : public CDynamicMesh<V,E,F,H>
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
*/

typedef CDynamicRicciFlowMesh<CDynamicYamabeFlowVertex, CDynamicYamabeFlowEdge, CDynamicYamabeFlowFace, CDynamicYamabeFlowHalfEdge> CDYFMesh;	

} //namespace RicciFlow
} //namespace MeshLib

#endif  _DYNAMIC_YAMABE_FLOW_MESH_H_