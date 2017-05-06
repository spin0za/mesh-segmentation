/*! \file UnifiedRicciFlow.h
 *  \brief Unified Ricci flow algorithm
 *  \author David Gu
 *  \date   documented on 01/01/2014
 *
 *	Algorithm for general Ricci Flow
 */

#ifndef _UNIFIED_RICCI_FLOW_H_
#define _UNIFIED_RICCI_FLOW_H_

#include <cmath>
#include <map>
#include <vector>
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
#include "Operator/Operator.h"
#include "Laplace/Poisson.h"



namespace MeshLib
{

namespace RicciFlow
{

enum Geometry { EUCLIDEAN, SPHERICAL, HYPERBOLIC };
enum Scheme   { TANGENTIAL, INVERSIVE_DISTANCE, YAMABE, VIRTUAL_RADIUS };


/*! \brief CUnifiedRicciFlow class
*
*	Algorithm for computing general Ricci flow
*/
template<typename M>
  class CUnifiedRicciFlow
  {
  public:
    /*! \brief CUnifiedRicciFlow constructor
	 *  \param pMesh the input mesh
	 */
    CUnifiedRicciFlow( M * pMesh, Geometry geometry, Scheme scheme );
    /*! \brief CBaseRicciFlow destructor
	 */
    ~CUnifiedRicciFlow();
	/*!	Computing the metric
	 */
	virtual void _calculate_metric( double threshold = 5e-4, double step_length = 5e-2 );

	/*!	Computing the metric
	 */
	virtual void _calculate_metric_torus_with_holes( double threshold = 5e-4, double step_length = 5e-2 );

	Geometry & geometry() { return m_geometry; };
	Scheme   & scheme()   { return m_scheme;   };

	virtual void set_target_curvature();
	virtual double update_target_curvature();

  protected:
    /*!
     *	the input mesh
	 */
    M	  * m_pMesh;
	/*!
	 *	the boundary of the mesh
	 */
	typename M::CBoundary m_boundary;

  protected:

	 void _vertex_radius();
	 void _edge_length();
	 void _corner_angle();
	 void _vertex_curvature();
	 void _face_Hessian();

	 Geometry m_geometry;
	 Scheme   m_scheme;

  private:

	 /*!
	  *	from conformal factor to compute circle radii
	  */
	  double __u2r( double u );
	  /*!
	   * from circle radius to conformal factor
	   */
	  double __r2u( double r );
	  /*!
	   *	compute edge length
	   */
	  double __l( double u1, double e1, double u2, double e2, double eta );
	  /*!
	   *	inverse cosine law, return angle A
	   */
	  double __inverse_cosine_law( double a, double b, double c ); 
	  /*!
	   *	tau(i,j,k)
	   */
	 double __tau( double li, double e_j, double r_j, double e_k, double r_k);
	/*!
	 *	s(x)
	 */
	double __s( double x );
	/*!
	 *	A( theta_k, l_i, l_j );
	 */
	double __A( double theta_k, double l_i, double l_j );
	/*!
	 *	Jacobian matrix
	 */
	Eigen::Matrix3d _jacobian( double a[3], double l[3], double r[3], double e[3] );
	void _face_jacobian( typename M::CFace * f );

	/*!
	 *	normalize discrete conformal factor
	 */
	void __normalize();

  public:

	/*!
	 *	Calculate vertex curvature error
	 */
   virtual double   _calculate_curvature_error();


    /*!
	 *	Curvature flow 
	 */
     virtual bool   _flow( double );
	 /*!
	  *	Newton's method to optimize the entropy energy
	  * \param threshold err bound
	  * \param step_length step length
	  */
     virtual void   _Newton( double threshold, double step_length );

  private:
	  /*!
	   *	set up initial conditions
	   */
	  void __tangential_CP_E2();
	  void __tangential_CP_S2();
	  void __tangential_CP_H2();

	  void __inversive_distance_CP_E2();
	  void __inversive_distance_CP_S2();
	  void __inversive_distance_CP_H2();

	  void __Yamabe_flow_E2();
	  void __Yamabe_flow_S2();
	  void __Yamabe_flow_H2();

	  void __virtual_radius_CP_E2();
	  void __virtual_radius_CP_S2();
	  void __virtual_radius_CP_H2();
};

template<typename M>
void CUnifiedRicciFlow<M>::set_target_curvature()
{
	for( M::MeshVertexIterator viter( m_pMesh ); !viter.end(); viter ++ )
	{
		M::CVertex * pV = *viter;
		pV->target_k() = 0;
	}
};

template<typename M>
double CUnifiedRicciFlow<M>::update_target_curvature()
{
  for( M::MeshVertexIterator viter( m_pMesh ); !viter.end(); viter ++ )
  {
	  M::CVertex * v = *viter;
	  if( v->boundary() ) continue;
	  v->target_k() = 0;
  }


  std::vector<M::CLoop*> & loops = m_boundary.loops();

  int i = 0;

  double total_error = 0;

  for( size_t i = 0; i < loops.size(); i ++ )
  {
	  //if( i < 2 ) continue;
	  M::CLoop * pL = loops[i];
		
	  std::vector<M::CHalfEdge*> hes;

	  for( std::list<M::CHalfEdge*>::iterator hiter = pL->halfedges().begin();
		  hiter != pL->halfedges().end(); hiter ++ )
	  {
		  M::CHalfEdge * pH = *hiter;
		  hes.push_back( pH );
	  }
		
	  double length = 0;
	  for( size_t i = 0; i < hes.size(); i ++ )
	  {
		  M::CHalfEdge * pH = hes[i];
		  M::CEdge     * e = m_pMesh->halfedgeEdge( pH );
		  length += e->length();
	  }
	
	  double error = 0;

	  for( size_t i = 0; i < hes.size(); i ++ )
	  {
		  M::CHalfEdge * pH = hes[i];
		  M::CHalfEdge * nH = hes[(i+1)%hes.size()];
		  M::CVertex   * pV = m_pMesh->halfedgeTarget( pH );

		  M::CEdge * pE = m_pMesh->halfedgeEdge( pH );
		  M::CEdge * nE = m_pMesh->halfedgeEdge( nH );
		  double target_k = -2 * PI * ( pE->length() + nE->length() )/( 2 * length );
		  error += fabs( pV->target_k()  - target_k );
		  pV->target_k() = target_k;
	  }
	  total_error += error;
	  std::cout << "Error is " << error <<  " ";

  }

  std::cout << std::endl;

  return total_error;
};


//Constructor
template<typename M>
CUnifiedRicciFlow<M>::CUnifiedRicciFlow( M * pMesh, Geometry geometry, Scheme scheme ): m_pMesh( pMesh), m_boundary( pMesh ), m_geometry( geometry), m_scheme( scheme )
{

	//set vertex e
	for( M::MeshVertexIterator viter( m_pMesh ); !viter.end(); viter ++ )
	{
		M::CVertex * pV = * viter;
		switch( m_scheme )
		{
		case TANGENTIAL:
		case INVERSIVE_DISTANCE:
			pV->e() = 1;
			break;
		case YAMABE:
			pV->e() = 0;
			break;
		case VIRTUAL_RADIUS:
			pV->e() = -1;
			break;
		}
	}
	//set target curvature
	for( M::MeshVertexIterator viter( m_pMesh ); !viter.end(); viter ++ )
	{
		M::CVertex * pV = *viter;
		pV->target_k() = 0;
	}

	//compute edge length
	for( M::MeshEdgeIterator eiter( m_pMesh ); !eiter.end(); eiter ++ )
	{
		M::CEdge * pE = *eiter;
		pE->length() = m_pMesh->edgeLength( pE );
	}

	set_target_curvature();

	switch( m_geometry )
	{
	case EUCLIDEAN:
		switch( m_scheme )
		{
		case TANGENTIAL:
			__tangential_CP_E2();
			break;
		case INVERSIVE_DISTANCE:
			__inversive_distance_CP_E2();
			break;
		case YAMABE:
			__Yamabe_flow_E2();
			break;
		case VIRTUAL_RADIUS:
			__virtual_radius_CP_E2();
			break;
		}
		break;
	case HYPERBOLIC:
		switch( m_scheme )
		{
		case TANGENTIAL:
			__tangential_CP_H2();
			break;
		case INVERSIVE_DISTANCE:
			__inversive_distance_CP_H2();
			break;
		case YAMABE:
			__Yamabe_flow_H2();
			break;
		case VIRTUAL_RADIUS:
			__virtual_radius_CP_H2();
			break;
		}
		break;
	case SPHERICAL:
		switch( m_scheme )
		{
		case TANGENTIAL:
			__tangential_CP_S2();
			break;
		case INVERSIVE_DISTANCE:
			__inversive_distance_CP_S2();
			break;
		case YAMABE:
			__Yamabe_flow_S2();
			break;
		case VIRTUAL_RADIUS:
			__virtual_radius_CP_S2();
			break;
		}
		break;
	};

};
//Destructor
template<typename M>
CUnifiedRicciFlow<M>::~CUnifiedRicciFlow( )
{
};

template<typename M>
void CUnifiedRicciFlow<M>::_vertex_radius()
{
	for( M::MeshVertexIterator viter( m_pMesh ); !viter.end(); viter ++ )
	{
		M::CVertex * v = *viter;
		v->r() = __u2r( v->u() );
	}
};

template<typename M>
void CUnifiedRicciFlow<M>::_edge_length( )
{

	for( M::MeshEdgeIterator eiter(m_pMesh);  !eiter.end();  eiter++ )
  {
	  M::CEdge * e = *eiter;

	  M::CVertex * v1 = m_pMesh->edgeVertex1( e );
	  M::CVertex * v2 = m_pMesh->edgeVertex2( e );

      double u1  = v1->u();
	  double e1  = v1->e();
      double u2  = v2->u();
	  double e2  = v2->e();
	  double eta = e->eta();
	  
	  e->length() = __l( u1, e1, u2, e2, eta );
  }
};


//Calculate corner angle
template<typename M>
void CUnifiedRicciFlow<M>::_corner_angle( )
{
	for ( M::MeshFaceIterator fiter( m_pMesh); ! fiter.end(); fiter ++ )
	{
	  M::CFace * f = *fiter;

	  M::CHalfEdge * he[3];

      he[0] = m_pMesh->faceMostCcwHalfEdge( f ); 

	  for( int i = 0; i < 3; i ++ )
      {
        he[(i+1)%3] = m_pMesh->	faceNextCcwHalfEdge(he[i]);
      }

      double l[3];
      for(int i = 0; i < 3; i ++ )
      {
		  M::CEdge * e = m_pMesh->halfedgeEdge( he[i] );
          l[i] = e->length();
      }

      for(int i = 0; i < 3; i ++ )
      {
          he[(i+1)%3]->angle() = __inverse_cosine_law( l[i], l[(i+1)%3] ,l[(i+2)%3] );
      }
	}

};


//Calculate vertex curvature
template<typename M>
void CUnifiedRicciFlow<M>::_vertex_curvature()
{
	for( M::MeshVertexIterator viter( m_pMesh ); !viter.end(); viter ++ )
	{
	  M::CVertex * v = *viter;
      double k  = (v->boundary() )? PI: PI * 2;
	  for( M::VertexInHalfedgeIterator vhiter( m_pMesh, v ); !vhiter.end();  ++vhiter )
      {
		  M::CHalfEdge * he = *vhiter;
          k -= he->angle();
      }
      v->k() = k;
	}
};




template<typename M>
Eigen::Matrix3d CUnifiedRicciFlow<M>::_jacobian( double a[3], double l[3], double r[3], double e[3] )
{
	Eigen::Matrix3d L, Theta, D, IL, J;

	L.setZero(); IL.setZero(); Theta.setZero(); D.setZero(); J.setZero();

	//L
	for( int k = 0; k < 3; k ++ )  L(k,k) = __s( l[k] );
	//IL
	for( int k = 0; k < 3; k ++ ) IL(k,k) = 1.0/L(k,k);

	//Theta
	for( int k = 0; k < 3; k ++ ) Theta(k,k) = -1;
	Theta(0,1) = cos ( a[2] );
	Theta(0,2) = cos( a[1] );
	Theta(1,2) = cos( a[0] );

	Theta(1,0) = cos( a[2] );
	Theta(2,0) = cos( a[1] );
	Theta(2,1) = cos( a[0] );

	//D
	for( int i = 0; i < 3; i ++ )
	for( int j = 0; j < 3; j ++ )
	{
		if( i == j ) continue;
		int k = 3 - i - j;
		D(i,j) = __tau( l[i], e[j], r[j], e[k], r[k] );
	}

	J = L *  Theta * IL * D;

	//area
	double A = -1/( 2 * __A(a[0], l[1], l[2]) );

	J = J * A;

	return J;
}

template<typename M>
void CUnifiedRicciFlow<M>::_face_jacobian( typename M::CFace * f )
{
	M::CHalfEdge * h[3];
	int k = 0;

	for( M::FaceHalfedgeIterator fhiter( f ); !fhiter.end(); fhiter ++ )
	{
		M::CHalfEdge * he = *fhiter;
		h[k++] = he;
	}

	double a[3];
	for( int k = 0; k < 3; k ++ )
	{
		a[k] = h[(k+1)%3]->angle();
	}

	double r[3];
	double e[3];
	for( int k = 0; k < 3; k ++ )
	{
		M::CVertex  * pV = m_pMesh->halfedgeTarget( h[(k+1)%3] );
		r[k] = pV->r();
		e[k] = pV->e();
	}

	double l[3];
	for( int k = 0; k < 3; k ++ )
	{
		M::CEdge * e = m_pMesh->halfedgeEdge( h[k] );
		l[k] = e->length();
	}

	Eigen::Matrix3d J;
	J = _jacobian( a,l,r, e );

	h[0]->c_w_ii() = -J(2,2);
	h[1]->c_w_ii() = -J(0,0);
	h[2]->c_w_ii() = -J(1,1);

	h[0]->c_w_ij() = -J(1,2);
	h[1]->c_w_ij() = -J(2,0);
	h[2]->c_w_ij() = -J(0,1);
};


template<typename M>
void CUnifiedRicciFlow<M>::_face_Hessian()
{
	for( M::MeshFaceIterator fiter( m_pMesh ); !fiter.end(); fiter ++ )
	{
		M::CFace * f = *fiter;
		_face_jacobian( f );
	}
};

template<typename M>
double CUnifiedRicciFlow<M>::__r2u(double r )
{
	double u = 0;
	switch( m_geometry )
	{
	case EUCLIDEAN:
		u = log( r );
		break;
	case HYPERBOLIC:
		u = log( tanh( r/2.0 ) );
		break;
	case SPHERICAL:
		u = log( tan( r/2.0 ) );
		break;
	};
	return u;
};


template<typename M>
double CUnifiedRicciFlow<M>::__u2r( double u )
{

	double r = 0;
	switch( m_geometry )
	{
	case EUCLIDEAN:
		r = exp( u );
		break;
	case HYPERBOLIC:
		if( m_scheme == YAMABE )
			r = 1.0;
		else
			r = log( (1+exp(u) )/(1-exp(u) ) );
		break;
	case SPHERICAL:
		if( m_scheme == YAMABE )
			r = 1.0;
		else
			r = 2 * atan( exp(u) );
		break;
	}
	return r;
};

template<typename M>
double CUnifiedRicciFlow<M>::__l( double u1, double e1, double u2, double e2, double eta )
{

	double l;

	switch( m_geometry )
	{
	case EUCLIDEAN:
		l = sqrt( 2 * eta * exp(u1) * exp(u2) + e1 * exp( 2*u1 ) + e2 * exp( 2 * u2 ) );
		break;
	case HYPERBOLIC:
		l = (4 * eta * exp(u1+u2) + ( 1 + e1 * exp( 2*u1 ) ) * ( 1 + e2 * exp( 2 * u2 ) ))/ (( 1 - e1 * exp( 2 * u1 ) ) * ( 1 - e2 * exp(2 * u2)));
		l = log( l + sqrt( l * l - 1) ); // acosh( l );
		break;
	case SPHERICAL:
		l = (-4 * eta * exp(u1+u2) + ( 1 - e1 * exp( 2*u1 ) ) * ( 1 - e2 * exp( 2 * u2 ) ))/ (( 1 + e1 * exp( 2 * u1 ) ) * ( 1 + e2 * exp(2 * u2)));
		l = acos(l); 
		break;
	}
	return l;
};


template<typename M>
double CUnifiedRicciFlow<M>::__tau( double l_i, double e_j, double r_j, double e_k, double r_k )
{
	double s;

	switch( m_geometry )
	{
	case EUCLIDEAN:
		s = (l_i * l_i + e_j * r_j * r_j - e_k * r_k * r_k)/2.0;
		break;
	case HYPERBOLIC:
		s = cosh(l_i) * pow( cosh(r_j), e_j ) - pow( cosh( r_k), e_k );
		break;
	case SPHERICAL:
		s = cos(l_i) * pow( cos(r_j), e_j ) - pow( cos( r_k), e_k );
		break;
	}
	return s;
};

template<typename M>
double CUnifiedRicciFlow<M>::__s( double x )
{
	double r;
	switch( m_geometry )
	{
	case EUCLIDEAN:
		r = x;
		break;
	case HYPERBOLIC:
		r = sinh( x );
		break;
	case SPHERICAL:
		r = sin(x );
		break;
	};
	return r;
};

template<typename M> 
double CUnifiedRicciFlow<M>::__inverse_cosine_law( double a, double b, double c )
{
	double A;

	switch( m_geometry )
	{
	case EUCLIDEAN:
		{
			double t = (b*b + c*c - a*a) /(2.0 * b * c);
			
			if( t > 1 ) t  = 1.0;
			if( t < -1) t = -1.0;

			A = acos(t);
		}
		break;
	case HYPERBOLIC:
		A = acos( (cosh(b) * cosh( c ) - cosh(a) )/( sinh(b) * sinh( c ) ));
		break;
	case SPHERICAL:
		A = acos( (-cos(b) * cos( c )  + cos( a ) )/ ( sin( b ) * sin( c ) ));
		break;
	};
	return A;
};

template<typename M>
double CUnifiedRicciFlow<M>::__A( double theta_k, double l_i, double l_j )
{
	return sin( theta_k ) * __s( l_i ) * __s( l_j )/2.0;
}

//compute curvature error

template<typename M>
double CUnifiedRicciFlow<M>::_calculate_curvature_error()
{
  double max_error = -1;
  M::CVertex * vert = NULL;

  for(M::MeshVertexIterator viter( m_pMesh); !viter.end() ; viter ++ )
   {
	   M::CVertex * v = *viter;
	   double k = v->target_k() - v->k();      
       k = fabs( k );
	   if( k > max_error )
	   {
		max_error = k;
		vert = v;
	   }
   }
   printf("Vertex id is %d\n", vert->id() );
   return max_error; 
};

//compute metric

template<typename M>
void CUnifiedRicciFlow<M>::_calculate_metric( double threshold, double step_length )
{
  if( m_geometry == SPHERICAL )
  {
		_flow(step_length);
  }
  else
  {
	_Newton( threshold, step_length );     
  }

};

//compute metric

template<typename M>
void CUnifiedRicciFlow<M>::_calculate_metric_torus_with_holes( double threshold, double step_length )
{
  if( m_geometry == SPHERICAL )
  {
		_flow(step_length);
  }
  else
  {
	  while(true)
	  {
		 if( update_target_curvature() < 0.01 ) break;
		_Newton( threshold * 128, step_length );     
	  }
	  _Newton( threshold, step_length );     
  }

};

//gradient flow method

template<typename M>
bool CUnifiedRicciFlow<M>::_flow( double error_threshold )
{
  int num = m_pMesh->numVertices();

  for( int k = 0; k < 64000; k ++  )
	  {
 		  _vertex_radius();
		  _edge_length();
		  _corner_angle();
		  _vertex_curvature();

		  double error =  _calculate_curvature_error();
		  printf("Current error is %f\r\n", error );

		  if( error < error_threshold)  return true;
  	  

	  //set b 
		  for(M::MeshVertexIterator viter( m_pMesh ); !viter.end(); ++ viter )
		  {
			  M::CVertex * v = *viter;
			  double dif = v->target_k() - v->k();
			  v->u() += dif * 2e-2;
		  }
    }
    return false;
};

/*!
 * In Euclidean case, normalize the discrete conformal factor
 */
template<typename M>
void CUnifiedRicciFlow<M>::__normalize()
{
	if( m_geometry != EUCLIDEAN ) return;

	double sum = 0;
	for( M::MeshVertexIterator viter( m_pMesh ); !viter.end(); viter ++  )
	{
		M::CVertex * v = *viter;
		sum  += v->u();
	}
	sum  /= m_pMesh->numVertices();

	for( M::MeshVertexIterator viter( m_pMesh ); !viter.end(); viter ++  )
	{
		M::CVertex * v = *viter;
		v->u() -= sum;
	}

}
//Newton's method for optimizing entropy energy

template<typename M>
void CUnifiedRicciFlow<M>::_Newton( double threshold, double step_length )
{
	const clock_t begin_time = std::clock();

	int num = m_pMesh->numVertices();

	int round = 0;

  	while( true )
	{
		//the order of the following functions really matters
		_vertex_radius();
 		_edge_length();
		_corner_angle();
		_vertex_curvature();

		_face_Hessian();

		double error =  _calculate_curvature_error();
		//printf("Newton's Method: Current error is %f\r\n", error );
		if( error < threshold) break;
	

		for( M::MeshVertexIterator viter( m_pMesh ); !viter.end(); viter ++  )
		{
			M::CVertex * v = *viter;
			v->huv()[0] = v->u();	
			v->u() = v->target_k() - v->k();
			v->b() = v->target_k() - v->k();
		}

		CPoisson<M> P(m_pMesh );
		P.hyperbolic_solve();

		__normalize();

		for( M::MeshVertexIterator viter( m_pMesh ); !viter.end(); viter ++  )
		{
			M::CVertex * v = *viter;
			v->u() = v->huv()[0] + v->u() * step_length;
		}

		round ++;
  }
	std::cout << "Total Time " << float( std::clock () - begin_time ) /  CLOCKS_PER_SEC << std::endl;

	std::cout << "Newton's Method : Total Rounds " << round << std::endl;
};


template<typename M>
void CUnifiedRicciFlow<M>::__tangential_CP_E2()
{
	for( M::MeshVertexIterator viter( m_pMesh ); !viter.end(); viter ++ )
	{
		M::CVertex * pV = * viter;
		pV->r() = 1.0;
		pV->u() = 0;
		pV->e() = 1;
	}

	for( M::MeshEdgeIterator eiter( m_pMesh ); !eiter.end(); eiter ++ )
	{
		M::CEdge * pE = *eiter;
		pE->eta() = 1.0;
	}
};

template<typename M>
void CUnifiedRicciFlow<M>::__inversive_distance_CP_E2()
{
	for( M::MeshVertexIterator viter( m_pMesh ); !viter.end(); viter ++ )
	{
		CURFMesh::CVertex * pV = * viter;
		pV->e() = 1;
	}

	for( M::MeshVertexIterator viter( m_pMesh );  !viter.end(); viter ++ )
	{
		CURFMesh::CVertex * v = *viter;
		v->u()   = 1e+10;
	}
	for( M::MeshEdgeIterator eiter( m_pMesh ); !eiter.end(); eiter ++ )
	{
		CURFMesh::CEdge * e = *eiter;
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
			double u = log( r );
			if( v->u()>u )
			{
				v->u() = u;
				v->r() = r;
			}
			he = m_pMesh->faceNextCcwHalfEdge( he );
		}
	}

	for( M::MeshEdgeIterator eiter( m_pMesh ); !eiter.end();  eiter++ )
	{
		M::CEdge * e = *eiter;

		M::CVertex * v1 = m_pMesh->edgeVertex1( e );
		M::CVertex * v2 = m_pMesh->edgeVertex2( e );

		double l = e->length();
		double r1 = v1->r();
		double r2 = v2->r();
		e->eta() = ( l*l - r1*r1 - r2*r2  ) / ( 2*r1*r2 );
	}
};

template<typename M>
void CUnifiedRicciFlow<M>::__Yamabe_flow_E2()
{
	for( M::MeshVertexIterator viter( m_pMesh ); !viter.end(); viter ++ )
	{
		M::CVertex * pV = * viter;
		pV->e() = 0;
	}

	for( M::MeshEdgeIterator eiter( m_pMesh ); !eiter.end(); eiter ++ )
	{
		M::CEdge * e = *eiter;
		e->length() = m_pMesh->edgeLength( e );
	};

	for( M::MeshVertexIterator viter( m_pMesh );  !viter.end(); viter ++ )
	{
		M::CVertex * v = *viter;
		v->u()   = 0;
	}

	for( M::MeshEdgeIterator eiter( m_pMesh ); !eiter.end();  eiter++ )
	{
		M::CEdge * e = *eiter;
		double l = e->length();
		e->eta() = l * l/2.0;
	}
};

template<typename M>
void CUnifiedRicciFlow<M>::__virtual_radius_CP_E2()
{
	for( M::MeshVertexIterator viter( m_pMesh ); !viter.end(); viter ++ )
	{
		M::CVertex * pV = * viter;
		pV->e() = -1;
	}

	for( M::MeshEdgeIterator eiter( m_pMesh ); !eiter.end(); eiter ++ )
	{
		M::CEdge * e = *eiter;
		e->length() = m_pMesh->edgeLength( e );
	};


	for( M::MeshVertexIterator viter( m_pMesh );  !viter.end(); viter ++ )
	{
		M::CVertex * v = *viter;
		v->u()   = 1e+10;
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
			double u = log( r );
			if( v->u()>u )
			{
				v->u() = u;
				v->r() = r;
			}
			he = m_pMesh->faceNextCcwHalfEdge( he );
		}
	}

	for( M::MeshVertexIterator viter( m_pMesh );  !viter.end(); viter ++ )
	{
		M::CVertex * v = *viter;
		v->u()   = v->u() - 2.0;
	};

	for( M::MeshEdgeIterator eiter( m_pMesh ); !eiter.end();  eiter++ )
	{
		M::CEdge * e = *eiter;

		M::CVertex * v1 = m_pMesh->edgeVertex1( e );
		M::CVertex * v2 = m_pMesh->edgeVertex2( e );

		double l = e->length();
		
		double u1 = v1->u();
		double u2 = v2->u();

		e->eta() = ( l*l + exp(u1)*exp(u1) + exp(u2)*exp(u2)  ) / ( 2*exp(u1)*exp(u2) );
	};

};

template<typename M>
void CUnifiedRicciFlow<M>::__tangential_CP_H2()
{
	for( M::MeshVertexIterator viter( m_pMesh ); !viter.end(); viter ++ )
	{
		M::CVertex * pV = * viter;
		pV->target_k() = 0;
		pV->r() = 1.0;
		pV->u() = log(tanh(pV->r()/2));
		pV->e() = 1;
	}
	for( M::MeshEdgeIterator eiter( m_pMesh ); !eiter.end(); eiter ++ )
	{
		M::CEdge * pE = *eiter;
		pE->eta() = 1.0;
	}
};

template<typename M>
void CUnifiedRicciFlow<M>::__inversive_distance_CP_H2()
{
	for( M::MeshVertexIterator viter( m_pMesh );  !viter.end(); viter ++ )
	{
		M::CVertex * v = *viter;
		v->u()   = 1e+10;
		v->r()   = 1e+10;
		v->e()   = 1.0;
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
			double u = log( tanh(r/2.0) );
			if( v->r()>r )
			{
				v->r() = r;
				v->u() = u;
			}
			he = m_pMesh->faceNextCcwHalfEdge( he );
		}
	}

	for( M::MeshEdgeIterator eiter( m_pMesh ); !eiter.end();  eiter++ )
	{
		M::CEdge * e = *eiter;

		M::CVertex * v1 = m_pMesh->edgeVertex1( e );
		M::CVertex * v2 = m_pMesh->edgeVertex2( e );
		double l = e->length();
		double u1 = v1->u();
		double u2 = v2->u();
		e->eta() = ( cosh(l) * ( 1- exp(2*u1))*(1-exp(2*u2)) - (1+exp(2*u1))*(1+exp(2*u2)) )/( 4 * exp( u1+ u2) ); 
	}
};

template<typename M>
void CUnifiedRicciFlow<M>::__Yamabe_flow_H2()
{
	for( M::MeshVertexIterator viter( m_pMesh );  !viter.end(); viter ++ )
	{
		M::CVertex * v = *viter;
		v->e()   = 0.0;
		v->u()   = 0.0;
		v->r()   = 1.0;
	}

	for( M::MeshEdgeIterator eiter( m_pMesh ); !eiter.end(); eiter ++ )
	{
		M::CEdge * e = *eiter;
		e->length() = m_pMesh->edgeLength( e );
	};

	for( M::MeshEdgeIterator eiter( m_pMesh ); !eiter.end();  eiter++ )
	{
		M::CEdge * e = *eiter;

		M::CVertex * v1 = m_pMesh->edgeVertex1( e );
		M::CVertex * v2 = m_pMesh->edgeVertex2( e );

		double l = e->length();

		double u1 = v1->u();
		double u2 = v2->u();

		e->eta() = ( cosh(l) - 1 )/( 4 * exp( u1+ u2) ); 
	}

};

template<typename M>
void CUnifiedRicciFlow<M>::__virtual_radius_CP_H2()
{
	for( M::MeshVertexIterator viter( m_pMesh );  !viter.end(); viter ++ )
	{
		M::CVertex * v = *viter;
		v->u()   = 1e+10;
		v->r()   = 1e+10;
		v->e()   = -1.0;
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
			double u = log( tanh(r/2.0) );
			if( v->r()>r )
			{
				v->r() = r;
				v->u() = u;
			}
			he = m_pMesh->faceNextCcwHalfEdge( he );
		}
	}

	for( M::MeshEdgeIterator eiter( m_pMesh ); !eiter.end();  eiter++ )
	{
		M::CEdge * e = *eiter;

		M::CVertex * v1 = m_pMesh->edgeVertex1( e );
		M::CVertex * v2 = m_pMesh->edgeVertex2( e );

		double l = e->length();

		double u1 = v1->u();
		double u2 = v2->u();

		double e1 = v1->e();
		double e2 = v2->e();

		e->eta() = ( cosh(l) * ( 1 + exp(2*u1))*(1 + exp(2*u2)) - (1-exp(2*u1))*(1-exp(2*u2)) )/( 4 * exp( u1+ u2) ); 
	}
};


template<typename M>
void CUnifiedRicciFlow<M>::__tangential_CP_S2()
{

	for( M::MeshVertexIterator viter( m_pMesh ); !viter.end(); viter ++ )
	{
		M::CVertex * pV = * viter;
		pV->r() = 1.0/1024;
		pV->u() = log(tan(pV->r()/2));
		pV->e() = 1;
	}

	for( M::MeshEdgeIterator eiter( m_pMesh ); !eiter.end(); eiter ++ )
	{
		M::CEdge * pE = *eiter;
		pE->eta() = 1.0;
	}
};

template<typename M>
void CUnifiedRicciFlow<M>::__inversive_distance_CP_S2()
{
	for( M::MeshVertexIterator viter( m_pMesh ); !viter.end(); viter ++ )
	{
		M::CVertex * pV = * viter;
		pV->e() = 1;
	}

	for( M::MeshVertexIterator viter( m_pMesh );  !viter.end(); viter ++ )
	{
		M::CVertex * v = *viter;
		v->u()   = 1e+10;
	}
	
	double max_length = -1e+10;

	for( M::MeshEdgeIterator eiter( m_pMesh ); !eiter.end(); eiter ++ )
	{
		M::CEdge * e = *eiter;
		e->length() = m_pMesh->edgeLength( e );
		max_length = ( max_length > e->length() )? max_length: e->length();
	};

	for( M::MeshEdgeIterator eiter( m_pMesh ); !eiter.end(); eiter ++ )
	{
		M::CEdge * e = *eiter;
		e->length() /= (max_length * 512 );
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
			double u = log( tan(r/2) );
			if( v->u()>u )
			{
				v->u() = u;
				v->r() = r;
			}
			he = m_pMesh->faceNextCcwHalfEdge( he );
		}
	}

	for( M::MeshEdgeIterator eiter( m_pMesh ); !eiter.end();  eiter++ )
	{
		M::CEdge * e = *eiter;

		M::CVertex * v1 = m_pMesh->edgeVertex1( e );
		M::CVertex * v2 = m_pMesh->edgeVertex2( e );

		double l = e->length();
		double u1 = v1->u();
		double u2 = v2->u();
		e->eta() = ( cos(l) * ( 1+ exp(2*u1)) * (1+exp(2*u2)) - ( 1-exp(2*u1) ) * ( 1 - exp( 2 * u2) ) )/( -4*exp(u1+u2));
	}
};

//tangential 
template<typename M>
void CUnifiedRicciFlow<M>::__Yamabe_flow_S2()
{
	for( M::MeshVertexIterator viter( m_pMesh );  !viter.end(); viter ++ )
	{
		M::CVertex * v = *viter;
		v->e()   = 0.0;
		v->u()   = 0.0;
		v->r()   = 1.0;
	}

	double max_length = -1e+10;
	for( M::MeshEdgeIterator eiter( m_pMesh ); !eiter.end(); eiter ++ )
	{
		M::CEdge * e = *eiter;
		e->length() = m_pMesh->edgeLength( e );
		max_length = ( max_length > e->length() )? max_length:  e->length();
	};

	for( M::MeshEdgeIterator eiter( m_pMesh ); !eiter.end(); eiter ++ )
	{
		M::CEdge * e = *eiter;
		e->length() /= (max_length * 3.0);
		e->length() *=  PI;
	};

	for( M::MeshEdgeIterator eiter( m_pMesh ); !eiter.end();  eiter++ )
	{
		M::CEdge * e = *eiter;

		M::CVertex * v1 = m_pMesh->edgeVertex1( e );
		M::CVertex * v2 = m_pMesh->edgeVertex2( e );

		double l = e->length();

		double u1 = v1->u();
		double u2 = v2->u();

		e->eta() = ( cos(l) - 1 )/( -4 * exp( u1+ u2) ); 
	}
};

template<typename M>
void CUnifiedRicciFlow<M>::__virtual_radius_CP_S2()
{
	for( M::MeshVertexIterator viter( m_pMesh );  !viter.end(); viter ++ )
	{
		M::CVertex * v = *viter;
		v->u()   = 1e+10;
		v->r()   = 1e+10;
		v->e()   = -1.0;
	}

	double max_length = -1e+10;
	for( M::MeshEdgeIterator eiter( m_pMesh ); !eiter.end(); eiter ++ )
	{
		M::CEdge * e = *eiter;
		e->length() = m_pMesh->edgeLength( e );
		max_length = ( max_length > e->length() )? max_length: e->length();
	};

	for( M::MeshEdgeIterator eiter( m_pMesh ); !eiter.end(); eiter ++ )
	{
		M::CEdge * e = *eiter;
		e->length() /= ( max_length * 512.0);
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
			double u = log( tan(r/2.0) );
			if( v->r()>r )
			{
				v->r() = r;
				v->u() = u;
			}
			he = m_pMesh->faceNextCcwHalfEdge( he );
		}
	}

	for( M::MeshEdgeIterator eiter( m_pMesh ); !eiter.end();  eiter++ )
	{
		M::CEdge * e = *eiter;
		M::CVertex * v1 = m_pMesh->edgeVertex1( e );
		M::CVertex * v2 = m_pMesh->edgeVertex2( e );
		double l = e->length();
		double u1 = v1->u();
		double u2 = v2->u();
		e->eta() = ( cos(l) * ( 1 + exp(2*u1))*(1 + exp(2*u2)) - (1-exp(2*u1))*(1-exp(2*u2)) )/( -4 * exp( u1+ u2) ); 
	}

};


} //namespace RicciFlow
} //namespace MeshLib

#endif  _UNIFIED_RICCI_FLOW_H_