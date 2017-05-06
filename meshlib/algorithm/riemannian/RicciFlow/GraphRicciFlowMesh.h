/*!
*      \file GraphRicciFlowMesh.h
*      \brief Mesh for Ricci Flow for Graph visualization
*	   \author David Gu
*      \date Documented 03/23/2011
*
*/

#ifndef  _GRAPH_RICCI_FLOW_MESH_H_
#define  _GRAPH_RICCI_FLOW_MESH_H_

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

namespace MeshLib
{

namespace RicciFlow
{
/*!
*	\brief CGraphRicciFlowVertex class
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
class CGraphRicciFlowVertex : public CVertex
{

public:
	/*!
	 *	CGraphRicciFlowVertex constructor
	 */
	CGraphRicciFlowVertex() { m_index = 0; };
	/*!
	 *	CGraphRicciFlowVertex destructor
	 */
	~CGraphRicciFlowVertex() {};

	int &     idx() { return m_index; };

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


protected:
	/*! Vertex index */
	int     m_index;
	//log radius
	double  m_log_radius;
	//curvature
	double  m_curvature;
	//target curvature
	double  m_target_curvature;
};


/*!
*	\brief CGraphRicciFlowEdge class
*
*	Edge class for Ricci flow
*   trait:
*   edge length
*   edge weight \$f \frac{\partial \theta_i}{\partial u_j} = w_k \f$
*   edge inversive distance \f$ \cos \phi \f$, intersection angle is \f$\phi\f$
*/
class CGraphRicciFlowEdge : public  CEdge
  {
  public:
    /*!	CRicciFlowEdge constructor
	 */
	 CGraphRicciFlowEdge() { m_length = 0; m_weight = 0; m_sharp = false; };
    /*!	CRicciFlowEdge destructor
	 */
    ~CGraphRicciFlowEdge(){};
	/*!	edge weight trait
	 */
	double & weight() { return m_weight; };
	/*!	edge length trait
	 */
	double & length() { return m_length; };
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
				  continue;
			  }
			  if( token->m_key == "l" )
			  {
				  m_length = strutil::parseString<double>( token->m_value );
				  continue;
			  }
		  }	
	};
	/*!
	 *	write no traits to edge string
	 */
	void _to_string() 
	{ 
		CParser parser( m_string );
		parser._removeToken( "l" );
		parser._toString( m_string );
		
		
		std::string line;
		std::stringstream iss(line);
		iss << "l=(" << m_length << ")";

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

  protected:
	  /*! edge weight trait */
	double   m_weight;
	  //edge length
	double   m_length;
	//inversive distance
	double   m_inversive_distance;
	/*!
	 *	if edge is sharp
	 */
	bool     m_sharp;
};

/*!
*	\brief CGraphRicciFlowHalfEdge class
*
*	HalfEdge class for Ricci flow
*   trait:
*   corner angle
*/

class CGraphRicciFlowHalfEdge : public  CHalfEdge
{
public:
  CGraphRicciFlowHalfEdge() { m_angle = 0.0; };
  ~CGraphRicciFlowHalfEdge() {};
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
*	\brief CGraphRicciFlowFace class
*
*	Face class for Ricci flow
*   trait:
*	whether the face has been touched
*/

class CGraphRicciFlowFace: public CFace
{
public:
	CGraphRicciFlowFace() { m_touched = false; };
	~CGraphRicciFlowFace() {};
	bool & touched() { return m_touched; };
	CPoint & normal(){ return m_normal; };
protected:
	bool m_touched;
	CPoint m_normal;
};

/*-------------------------------------------------------------------------------------------------------------------------------------

	Harmonic Mapper Mesh Class

--------------------------------------------------------------------------------------------------------------------------------------*/
/*!
 *	\brief CHarmonicMapperMesh class
 *
 *	Mesh class for harmonic mapping purpose
 */
template<typename V, typename E, typename F, typename H>
class CGraphRicciFlowMesh : public CBaseMesh<V,E,F,H>
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


typedef CGraphRicciFlowMesh<CGraphRicciFlowVertex, CGraphRicciFlowEdge, CGraphRicciFlowFace, CGraphRicciFlowHalfEdge> CGRFMesh;	

} //namespace RicciFlow
} //namespace MeshLib
#endif  _GRAPH_RICCI_FLOW_MESH_H_