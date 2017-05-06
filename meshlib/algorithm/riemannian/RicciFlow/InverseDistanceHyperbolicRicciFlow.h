/*! \file InversiveDistanceHyperbolicRicciFlow.h
 *  \brief Inversive Distance Hyperbolic Ricci flow algorithm
 *  \author David Gu
 *  \date   documented on 10/17/2010
 *
 *	Algorithm for Inversive Distance Hyperbolic Ricci Flow
 */

#ifndef _INVERSIVE_DISTANCE_HYPERBOLIC_RICCI_FLOW_H_
#define _INVERSIVE_DISTANCE_HYPERBOLIC_RICCI_FLOW_H_

#include <map>
#include <vector>
#include <Eigen/Dense>
#include <Eigen/Sparse>
#include "Mesh/BaseMesh.h"
#include "Mesh/Vertex.h"
#include "Mesh/HalfEdge.h"
#include "Mesh/Edge.h"
#include "Mesh/Face.h"
#include "mesh/iterators.h"
#include "mesh/boundary.h"
#include "Parser/parser.h"
#include "RicciFlowMesh.h"
#include "BaseRicciFlow.h"


namespace MeshLib
{
	
namespace RicciFlow
{
/*! \brief Class CInversiveDistanceHyperbolicRicciFlow
*
*	Algorithm for computing Hyperbolic Ricci flow using inversive distance circle packing metric
*/
template<typename M>
class CInversiveDistanceHyperbolicRicciFlow : public CBaseHyperbolicRicciFlow<M>
  {
  public:
    /*! \brief CTangentialRicciFlow constructor
	 *  \param pMesh the input mesh
	 */
	  CInversiveDistanceHyperbolicRicciFlow( M * pMesh );
    /*! \brief CTangentialRicciFlow destructor
	 */
	  ~CInversiveDistanceHyperbolicRicciFlow(){};
	

  protected:

	  /*!
	   *	Calculate each edge length, has to be defined in the derivated classes
	   */
	  void _length( double u1, double u2, typename M::CEdge * e );

	/*!
	 *	Cosine law, has to be defined in the derivated classes
	 */
	//double _cosine_law( double a, double b, double c ) ;

	/*!
	 *	inverse cosine law, has to be defined in the derivated classes
	 *  \param a edge length a
	 *  \param b edge length b
	 *  \param iv inversive distance, corresponding to the cosine value of the angle between a and b
	 *  return the length of the 3rd edge
	 */
	double _inverse_cosine_law( double a, double b, double iv ) ;


	/*!
	 *	Calculate the edge weight
	 */
    void _calculate_edge_weight();

	/*!
	 *	Set the target curvature on each vertex
	 */
    void    _set_target_curvature();

	 /*! compute initial vertex radius and edge length, edge inversive distance */
	void _calculate_initial_inversive_distance();
	/*! computing the Jacobi matrix \f$ \frac{\partial \theta_i}{\partial u_j} \f$ in the triangle f
	* \param f the input triangle
	*/
	void _face_jacobi( typename M::CFace * f );
	/*! computing the Jacobi matrix \f$ \frac{\partial \theta_i}{\partial u_j} \f$ in the triangle 
	* \param a three angles
	* \param l three edge lengths
	* \param r three radii
	* \output CMatrix the Jacobi matrix
	*/
	
	Eigen::Matrix3d _jacobi( double a[3], double l[3], double r[3] );

  };


template<typename M>
CInversiveDistanceHyperbolicRicciFlow<M>::CInversiveDistanceHyperbolicRicciFlow( M * pMesh ) : CBaseHyperbolicRicciFlow<M>( pMesh )
{
	//m_pMesh = pMesh;
	int vid = 0;
	for( M::MeshVertexIterator viter(m_pMesh); !viter.end(); viter ++ )
	{
		M::CVertex * pV = *viter;
		pV->idx() = vid ++;
	}
	_calculate_initial_inversive_distance();
};

//Compute the edge length
template<typename M>
void CInversiveDistanceHyperbolicRicciFlow<M>::_length( double u1, double u2, typename M::CEdge * e )
{
	double r1 = _u2r( u1 );
    double r2 = _u2r( u2 );
	double iv = e->inversive_distance();

    double l = _inverse_cosine_law( r1, r2, iv );
	e->length() = l;
};

//given edge lengths a,b and angle C, return edge length c

template<typename M>
double CInversiveDistanceHyperbolicRicciFlow<M>::_inverse_cosine_law( double a, double b, double iv )
{
  double c;
  c = cosh(a) * cosh(b) + sinh(a)*sinh(b) * iv ;
  c = log( c+sqrt(c*c-1) );
  return c;
}


//Calculate edge weight

template<typename M>
void CInversiveDistanceHyperbolicRicciFlow<M>::_calculate_edge_weight()
{
  for( M::MeshFaceIterator fiter( m_pMesh ); !fiter.end(); ++ fiter)
  {
	  M::CFace * f = *fiter;
	  _face_jacobi( f );
  }
};

//set target curvature

template<typename M>
void CInversiveDistanceHyperbolicRicciFlow<M>::_set_target_curvature()
{
  for( M::MeshVertexIterator viter( m_pMesh ); !viter.end(); viter ++ )
  {
	  M::CVertex * v = *viter;
    v->target_k() = 0;
  }
};


template<typename M>
void CInversiveDistanceHyperbolicRicciFlow<M>::_calculate_initial_inversive_distance()
{
	// set initial inversive distance and radius

	for( M::MeshVertexIterator viter( m_pMesh );  !viter.end(); viter ++ )
	{
		M::CVertex * v = *viter;
		v->u()   = INT_MAX;
	}
	for( M::MeshEdgeIterator eiter( m_pMesh ); !eiter.end(); eiter ++ )
	{
		M::CEdge * e = *eiter;
		e->length() = m_pMesh->edgeLength( e );
	};

	for( M::MeshFaceIterator fiter( m_pMesh ); !fiter.end(); fiter ++ )
	{
		M::CFace * f = *fiter;

		double edgeLen[3];

		M::CHalfEdge * he = m_pMesh->faceMostCcwHalfEdge( f );
		for( int i = 0; i < 3; i ++ )
		{
			M::CEdge * e = m_pMesh->halfedgeEdge( he );
			edgeLen[i] = e->length();			
			he = m_pMesh->faceNextCcwHalfEdge( he );
		}

		he = m_pMesh->faceMostCcwHalfEdge( f );
		for( int i=0; i<3; i++ ) {
			M::CVertex * v = m_pMesh->halfedgeSource( he );
			double r = (edgeLen[i] + edgeLen[(i+2)%3] - edgeLen[(i+1)%3])/2;
			double u = _r2u( r );
			if( v->u()>u )
				v->u() = u;
			he = m_pMesh->faceNextCcwHalfEdge( he );
		}
	}

	for( M::MeshEdgeIterator eiter( m_pMesh ); !eiter.end();  eiter++ )
	{
		M::CEdge * e = *eiter;

		M::CVertex * v1 = m_pMesh->edgeVertex1( e );
		M::CVertex * v2 = m_pMesh->edgeVertex2( e );

		double l = (v1->point()-v2->point()).norm();

		double u1 = v1->u();
		double u2 = v2->u();

		double r1 = _u2r( u1 );
		double r2 = _u2r( u2 );
		e->inversive_distance() = ( cosh( l ) - cosh( r1 ) * cosh( r2 ) )/( sinh(r1) * sinh(r2) );
	}

};

/*
#define a1 (a[0])
#define a2 (a[1])
#define a3 (a[2])

#define l1 (l[0])
#define l2 (l[1])
#define l3 (l[2])

#define r1 (r[0])
#define r2 (r[1])
#define r3 (r[2])
*/
template<typename M>
Eigen::Matrix3d CInversiveDistanceHyperbolicRicciFlow<M>::_jacobi( double a[3], double l[3], double r[3] )
{
	Eigen::Matrix3d R, M1, M2,M3, M4;

	R.setZero(); M1.setZero();M2.setZero();M3.setZero();M4.setZero();

	//M1
	for( int k = 0; k < 3; k ++ ) M1(k,k) = sinh( l[k] );
	
	//M2
	for( int k = 0; k < 3; k ++ ) M2(k,k) = -1;
	M2(0,1) = cos ( a[2] );
	M2(0,2) = cos( a[1] );
	M2(1,2) = cos( a[0] );

	M2(1,0) = cos( a[2] );
	M2(2,0) = cos( a[1] );
	M2(2,1) = cos( a[0] );

	//M4
	for( int k = 0; k < 3; k ++ ) M4(k,k) = sinh( r[k] );

	//M3

	M3(0,1) = (-cosh(r[2])  + cosh(l[0]) *cosh(r[1]) )/(sinh(l[0]) *sinh(r[1]) );
	M3(0,2) = (-cosh(r[1]) + cosh(l[0]) *cosh(r[2])  )/(sinh(l[0]) *sinh(r[2]) );

	M3(1,0) = (-cosh(r[2]) + cosh(l[1]) *cosh(r[0])  )/(sinh(l[1]) *sinh(r[0]) );
	M3(1,2) = (-cosh(r[0]) + cosh(l[1]) *cosh(r[2])  )/(sinh(l[1]) *sinh(r[2]) );

	M3(2,0) = (-cosh(r[1]) + cosh(l[2]) *cosh(r[0])  )/(sinh(l[2]) *sinh(r[0]) );
	M3(2,1) = (-cosh(r[0]) + cosh(l[2]) *cosh(r[1])  )/(sinh(l[2]) *sinh(r[1]) );

	R = M1 *  M2;
	R = R  *  M3;
	R = R  *  M4;

	//area
	double A = -1/( sin(a[0])  * sinh(l[1]) * sinh(l[2]) );

	R = R * A;

	return R;
}

template<typename M>
void CInversiveDistanceHyperbolicRicciFlow<M>::_face_jacobi( typename M::CFace * f )
{
	M::CHalfEdge * h[3];
	int k = 0;

	for( M::FaceHalfedgeIterator fhiter( f ); !fhiter.end(); fhiter ++ )
	{
		M::CHalfEdge * he = *fhiter;
		h[k++] = he;
	}
/*
	M:CVertex * v[3];

	for( int k = 0; k < 3; k ++ )
	{
		v[k] = m_pMesh->halfedgeTarget( h[(k+1)%3] );
	}
*/
	double a[3];
	for( int k = 0; k < 3; k ++ )
	{
		a[k] = h[(k+1)%3]->angle();
	}

	double r[3];
	for( int k = 0; k < 3; k ++ )
	{
		M::CVertex  * pV = m_pMesh->halfedgeTarget( h[(k+1)%3] );
		r[k] = _u2r(  pV->u() );
	}

	double l[3];
	for( int k = 0; k < 3; k ++ )
	{
		M::CEdge * e = m_pMesh->halfedgeEdge( h[k] );
		l[k] = e->length();
	}

	Eigen::Matrix3d M;
	M = _jacobi( a,l,r);
//	M.output();

	h[0]->c_w_ii() = -M(2,2);
	h[1]->c_w_ii() = -M(0,0);
	h[2]->c_w_ii() = -M(1,1);

	h[0]->c_w_ij() = -M(1,2);
	h[1]->c_w_ij() = -M(2,0);
	h[2]->c_w_ij() = -M(0,1);

}



} //namespace RicciFlow
} //namespace MeshLib

#endif  _INVERSIVE_DISTANCE_HYPERBOLIC_RICCI_FLOW_H_