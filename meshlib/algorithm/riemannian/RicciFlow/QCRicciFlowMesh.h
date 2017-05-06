/*!
*      \file QCRicciFlowMesh.h
*      \brief Mesh for Quasi-Conformal Ricci Flow
*	   \author David Gu
*      \date Documented 03/08/2011
*
*/

#ifndef  _QC_RICCI_FLOW_MESH_H_
#define  _QC_RICCI_FLOW_MESH_H_

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
#include "RicciFlowMesh.h"
#include <complex>

namespace MeshLib
{

namespace RicciFlow
{
/*!
*	\brief CQCRicciFlowVertex class
*
*	Vertex class for Quasi-Conformal Ricci flow
*   Traits:
*   index
*	father
*   log_radius
*   curvature
*   vertex uv
*   target_curvature
*   touched
*   mu
*/
class CQCRicciFlowVertex : public CVertex
{

public:
	/*!
	 *	CQCRicciFlowVertex constructor
	 */
	CQCRicciFlowVertex() { m_index = 0; };
	/*!
	 *	CRicciFlowVertex destructor
	 */
	~CQCRicciFlowVertex() {};

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
	CPoint & rgb() { return m_rgb; };
	/*!
	 *	Vertex Beltrami coefficient
	 */ 
	std::complex<double> & mu() { return m_mu; };
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
	CPoint  m_rgb;
	//vertex Beltrami coefficient
	std::complex<double> m_mu;
};

//read father from string
inline void CQCRicciFlowVertex::_from_string()
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
	  if( token->m_key == "rgb" )
	  {
		 std::string line = strutil::trim( token->m_value, "()");
		 line >> m_rgb;		
	  }
	  if( token->m_key == "uv" )
	  {
		 std::string line = strutil::trim( token->m_value, "()");
		 line >> m_huv;		
	  }
	  if( token->m_key == "mu" )
	  {
	     CPoint2 t;
		 std::string line = strutil::trim( token->m_value, "()");
		 line >> t;
		 m_mu = std::complex<double>( t[0], t[1] );
	  }
  }
}

//converting vertex uv trait to string

inline void CQCRicciFlowVertex::_to_string()
{
	CParser parser( m_string );
	parser._removeToken( "uv" );
	parser._removeToken( "mu" );
	parser._toString( m_string );
	std::stringstream iss;
	iss << "uv=(" << m_huv[0] << " " << m_huv[1] << ")";
	iss << " ";
	iss << "mu=(" << m_mu.real() << " " << m_mu.imag() << ")";

	if( m_string.length() > 0 )
	{
		m_string += " ";
	}
	m_string += iss.str();
}

/*!
*	\brief CQCRicciFlowEdge class
*
*	Edge class for Quasi-Conformal Ricci flow
*   trait:
*   edge length
*   edge weight \$f \frac{\partial \theta_i}{\partial u_j} = w_k \f$
*   edge inversive distance \f$ \cos \phi \f$, intersection angle is \f$\phi\f$
*/
class CQCRicciFlowEdge : public  CEdge
  {
  public:
    /*!	CRicciFlowEdge constructor
	 */
	 CQCRicciFlowEdge() { m_length = 0; m_weight = 0; m_sharp = false; };
    /*!	CRicciFlowEdge destructor
	 */
    ~CQCRicciFlowEdge(){};
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
	/*!	edge length trait using Auxillary metric
	 */
	double & mu_length() { return m_mu_length; };
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
	  //edge auxillary metric length
	double   m_mu_length;
	//inversive distance
	double   m_inversive_distance;
	/*!
	 *	if edge is sharp
	 */
	bool     m_sharp;
};

/*!
*	\brief CQCRicciFlowHalfEdge class
*
*	HalfEdge class for Quasi-Conformal Ricci flow
*   trait:
*   corner angle
*/

class CQCRicciFlowHalfEdge : public  CHalfEdge
{
public:
  CQCRicciFlowHalfEdge() { m_angle = 0.0; };
  ~CQCRicciFlowHalfEdge() {};
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
*	\brief CQCRicciFlowFace class
*
*	Face class for Quasi-Conformal Ricci flow
*   trait:
*	whether the face has been touched
*/

class CQCRicciFlowFace: public CFace
{
public:
	CQCRicciFlowFace() { m_touched = false; };
	~CQCRicciFlowFace() {};
	bool & touched() { return m_touched; };
	CPoint & normal(){ return m_normal; };
protected:
	bool m_touched;
	CPoint m_normal;
};

/*!
 *	\brief CQCRicciFlowMesh class
 *
 *	Mesh class for harmonic mapping purpose
 */
/*
template<typename V, typename E, typename F, typename H>
class CRicciFlowMesh : public CBaseMesh<V,E,F,H>
{
public:
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

typedef CRicciFlowMesh<CQCRicciFlowVertex, CQCRicciFlowEdge, CQCRicciFlowFace, CQCRicciFlowHalfEdge> CQCRFMesh;	

}; //namespace RicciFlow
}; //namespace MeshLib

#endif  _QC_RICCI_FLOW_MESH_H_