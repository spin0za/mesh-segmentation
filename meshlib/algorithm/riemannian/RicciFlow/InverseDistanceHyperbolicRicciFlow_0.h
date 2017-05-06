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
class CInversiveDistanceHyperbolicRicciFlow 
  {
  public:
    /*! \brief CTangentialRicciFlow constructor
	 *  \param pMesh the input mesh
	 */
	  CInversiveDistanceHyperbolicRicciFlow( M * pMesh );
    /*! \brief CTangentialRicciFlow destructor
	 */
	  ~CInversiveDistanceHyperbolicRicciFlow(){};
	/*!	Computing the metric
	 */
	void _calculate_metric();
	
	void _calculate_corner_angle( );
	double _calculate_curvature_error();
	void _calculate_vertex_curvature( );
	void _calculate_edge_length( );

  protected:

	void _Newton( double threshold, double step_length );

	  /*!
	   *	Calculate each edge length, has to be defined in the derivated classes
	   */
	  void _length( double u1, double u2, typename M::CEdge * e );

	/*!
	 *	Cosine law, has to be defined in the derivated classes
	 */
	double _cosine_law( double a, double b, double c ) ;

	/*!
	 *	inverse cosine law, has to be defined in the derivated classes
	 *  \param a edge length a
	 *  \param b edge length b
	 *  \param iv inversive distance, corresponding to the cosine value of the angle between a and b
	 *  return the length of the 3rd edge
	 */
	double _inverse_cosine_law( double a, double b, double iv ) ;
	 /*!
	  *	Normalization
	  * \param du the du vector
	  * \param n dimension of the du vector
	  *  do nothing
	  */
	void _normalization( Eigen::VectorXd & du, int n );


	/*!
	 *	Calculate the edge weight
	 */
    void _calculate_edge_weight();
	/*!
	 *	Calculate the edge weight for a halfedge
	 *  \param pH input halfedge
	 */
	void _calculate_edge_weight( typename M::CHalfEdge * pH );

	/*!
	 *	Set the target curvature on each vertex
	 */
    void    _set_target_curvature();
	/*!
	 *	Curvature flow 
	 */
    bool   _flow( double );

	 /*!
	  *	calculate hessian matrix Hessain 
	  * \param SparseMatrix
	  */
	void _calculate_Hessain( Eigen::SparseMatrix<double> & M );

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

	M * m_pMesh;

  private:
	  /*!
	   *	Converting u to radius r
	   */
	  double _u2r( double u );
	  /*!
	   *	Converting radius r to u
	   */
	  double _r2u( double r );
  };


template<typename M>
CInversiveDistanceHyperbolicRicciFlow<M>::CInversiveDistanceHyperbolicRicciFlow( M * pMesh )
{
	m_pMesh = pMesh;
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


//Calculate corner angle
template<typename M>
double CInversiveDistanceHyperbolicRicciFlow<M>::_cosine_law( double a, double b, double c )
{
	double C;
	C = acos( (cosh(a) * cosh(b)-cosh(c))/(sinh(a)*sinh(b)) );
	return C;
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

// u = log tanh (r/2)
template<typename M>
double CInversiveDistanceHyperbolicRicciFlow<M>::_u2r( double u )
{  
    double e = exp(u);
    double r;
    r = log(  (1+e)/(1 - e ) );
    return r;
}

// u = log tanh (r/2)
template<typename M>
double CInversiveDistanceHyperbolicRicciFlow<M>::_r2u( double r )
{  
    return log(tanh(r/2));
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

//compute metric

template<typename M>
void CInversiveDistanceHyperbolicRicciFlow<M>::_calculate_metric()
{

	_set_target_curvature();
	double error_threshold = 1e-6;
	//double step_length = 1.0;
	double step_length = 0.5;
	_Newton( error_threshold, step_length );

};

//gradient flow method

template<typename M>
bool CInversiveDistanceHyperbolicRicciFlow<M>::_flow( double error_threshold )
{
  int num = m_pMesh->numVertices();

  for( int k = 0; k < 64; k ++  )
	  {
 		  _calculate_edge_length();
	      _set_target_curvature();
		  _calculate_edge_weight();

		  _calculate_corner_angle();
		  _calculate_vertex_curvature();

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

//normalization

template<typename M>
void CInversiveDistanceHyperbolicRicciFlow<M>::_normalization( Eigen::VectorXd & x, int num )
{
	//for hyperbolic Ricci flow, do nothing
};

template<typename M>
void CInversiveDistanceHyperbolicRicciFlow<M>::_calculate_Hessain( Eigen::SparseMatrix<double> & M )
{

	std::vector<Eigen::Triplet<double> > M_coefficients;

  for( M::MeshVertexIterator viter( m_pMesh ); !viter.end(); ++ viter )
  {
	  M::CVertex * v = *viter;

      double a_ii = 0;

	  for( M::VertexInHalfedgeIterator hiter( m_pMesh, v ); !hiter.end(); ++ hiter )
	  {
		  M::CHalfEdge * he = *hiter;
	  	a_ii -= he->c_w_ii();
	  }
	  M_coefficients.push_back( Eigen::Triplet<double> ( v->idx(), v->idx(), a_ii ) );
      //pMatrix->AddElementTail( v->idx(), v->idx(), a_ii);   
  }

  for( M::MeshEdgeIterator eiter( m_pMesh ); !eiter.end(); ++ eiter )
  {
	  M::CEdge * e = *eiter;

      double a_ij = 0;

	  M::CHalfEdge * he = m_pMesh->edgeHalfedge(e,0);

      a_ij -= he->c_w_ij();

	  M::CHalfEdge *  twin = m_pMesh->halfedgeSym( he );
      if( twin != NULL )
      {
            a_ij -= twin->c_w_ij();
      }

	  M::CVertex * v0 = m_pMesh->edgeVertex1( e );
	  M::CVertex * v1 = m_pMesh->edgeVertex2( e );

      //pMatrix->AddElementTail( v0->idx(), v1->idx(),  a_ij );   
      //pMatrix->AddElementTail( v1->idx(), v0->idx(),  a_ij );   
	  M_coefficients.push_back( Eigen::Triplet<double> ( v0->idx(), v1->idx(), a_ij ) );
	  M_coefficients.push_back( Eigen::Triplet<double> ( v1->idx(), v0->idx(), a_ij ) );

  }

  M.setFromTriplets(M_coefficients.begin(), M_coefficients.end());
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

	h[0]->c_w_ii() = M(2,2);
	h[1]->c_w_ii() = M(0,0);
	h[2]->c_w_ii() = M(1,1);

	h[0]->c_w_ij() = M(1,2);
	h[1]->c_w_ij() = M(2,0);
	h[2]->c_w_ij() = M(0,1);

}



//Newton's method for optimizing entropy energy

template<typename M>
void CInversiveDistanceHyperbolicRicciFlow<M>::_Newton( double threshold, double step_length )
{
	int num = m_pMesh->numVertices();

	
  	while( true )
	{
		//the order of the following functions really matters

 		_calculate_edge_length();
		_calculate_corner_angle();
		_calculate_vertex_curvature();
		_calculate_edge_weight();

		double error =  _calculate_curvature_error();
		printf("Newton's Method: Current error is %f\r\n", error );

		if( error < threshold) break;
	

		Eigen::SparseMatrix<double>  M( num, num );
		M.setZero();
		_calculate_Hessain( M );

		Eigen::VectorXd b(num);

		for( M::MeshVertexIterator viter( m_pMesh ); !viter.end(); viter ++  )
		{
			M::CVertex * v = *viter;
			int idx = v->idx();
			b(idx) = v->target_k() - v->k();
		}

		//std::cout << M;
		//Eigen::ConjugateGradient<Eigen::SparseMatrix<double>> solver;
		Eigen::SimplicialLDLT<Eigen::SparseMatrix<double>> solver;
		std::cerr << "Eigen Decomposition" << std::endl;
		solver.compute(M);
		std::cerr << "Eigen Decomposition Finished" << std::endl;
	
		if( solver.info() != Eigen::Success )
		{
			std::cerr << "Waring: Eigen decomposition failed" << std::endl;
		}

		Eigen::VectorXd x = solver.solve(b);
		if( solver.info() != Eigen::Success )
		{
			std::cerr << "Waring: Eigen decomposition failed" << std::endl;
		}
		//std::cout << x.transpose();

		_normalization( x, num );

		for( M::MeshVertexIterator viter( m_pMesh ); !viter.end(); viter ++  )
		{
			M::CVertex * v = *viter;
			int idx = v->idx();
			v->u() += x(idx) * step_length;
		}

  }

};




//Calculate corner angle
template<typename M>
void CInversiveDistanceHyperbolicRicciFlow<M>::_calculate_corner_angle( )
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
          he[(i+1)%3]->angle() = _cosine_law( l[(i+1)%3] ,l[(i+2)%3] ,l[i] );
      }
  }

};

//Calculate vertex curvature
template<typename M>
void CInversiveDistanceHyperbolicRicciFlow<M>::_calculate_vertex_curvature()
{

  for( M::MeshVertexIterator viter( m_pMesh ); !viter.end(); viter ++ )
  {
	  M::CVertex * v = *viter;
      double k  = (v->boundary() )? PI: PI * 2;
	  for( M::VertexInHalfedgeIterator vh( m_pMesh, v ); !vh.end();  ++vh )
      {
		  M::CHalfEdge * he = *vh;
          k -= he->angle();
      }
      v->k() = k;
  }

};


//compute curvature error

template<typename M>
double CInversiveDistanceHyperbolicRicciFlow<M>::_calculate_curvature_error()
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

//Compute the edge length
template<typename M>
void CInversiveDistanceHyperbolicRicciFlow<M>::_calculate_edge_length( )
{

	for( M::MeshEdgeIterator eiter(m_pMesh);  !eiter.end();  eiter++ )
  {
	  M::CEdge * e = *eiter;

	  M::CVertex * v1 = m_pMesh->edgeVertex1( e );
	  M::CVertex * v2 = m_pMesh->edgeVertex2( e );

      double u1 = v1->u();
      double u2 = v2->u();
	  
	  _length( u1, u2, e );
  }

};

} //namespace RicciFlow
} //namespace MeshLib

#endif  _INVERSIVE_DISTANCE_HYPERBOLIC_RICCI_FLOW_H_