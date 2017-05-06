/*! \file QCVirtualInversiveDistanceRicciFlow.h
 * \brief Quasi-Conformal Euclidean inversive distance Ricci flow with virtual radius
 *  \author David Gu
 *  \date   documented on 02/11/2011
 *
 *	Algorithm for Inversive distance  Ricci Flow with virtual radius
 */

#ifndef _QC_VIRTUAL_INVERSIVE_DISTANCE_RICCI_FLOW_H_
#define _QC_VIRTUAL_INVERSIVE_DISTANCE_RICCI_FLOW_H_

#include <map>
#include <vector>
#include <Eigen/Sparse>
#include <Eigen/Dense>

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
#include "Geometry/Circle.h"

namespace MeshLib
{

namespace RicciFlow
{
	/*! \brief CVirtualInversiveDistanceRicciFlow class
	 *  
	 *  Algorithm for Quasi-Conformal Inversive distance Ricci flow with virtual radius
	 */
	 
 template<typename M>
 class CQCVirtualInversiveDistanceRicciFlow : public CBaseRicciFlow<M>
 {
  public:
	  /*! \brief CQCVirtualInversiveDistanceRicciFlow constructor
	   *
	   *  call base class constructor 
	   */
	  CQCVirtualInversiveDistanceRicciFlow( M * pMesh );
	  /*! override the virtual function calculate_metric */
	  void _calculate_metric( double step_length );
	

  protected:

	  /*!
	   *	Calculate each edge length, has to be defined in the derivated classes
	   */
	  void _length( double u1, double u2, typename M::CEdge * e );

	/*!
	 *	Cosine law, has to be defined in the derivated classes
	 */
	double _cosine_law( double a, double b, double c ) ;
	 /*!
	  *	Normalization
	  * \param du the du vector
	  * \param n dimension of the du vector
	  */
	void _normalization( Eigen::VectorXd & du, int n );


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
	/*! calculate edge weight by power center */
	void _calculate_edge_weight_by_power_center();
	/*! calculate power center on a face head
	    \param head input face
	*/
	void _calculate_power_center( typename M::CFace* head );
	/*! normalize the mesh
	*/
	void _normalize_mesh();
 };

template<typename M>
CQCVirtualInversiveDistanceRicciFlow<M>::CQCVirtualInversiveDistanceRicciFlow( M * pMesh ): CBaseRicciFlow( pMesh)
{
	_normalize_mesh();

	for( M::MeshVertexIterator viter( m_pMesh ); !viter.end(); viter ++ )
	{
		M::CVertex * v = *viter;
	  CPoint2 p = v->huv();
	  std::complex<double> mu(p[0],p[1]);
	  double r = std::abs(mu);
	  double theta = std::arg(mu);
	  mu = std::polar(  r*(1-r), theta);
	  v->mu() = mu;
	}

	//set up the auxillary metric

	for( M::MeshEdgeIterator eiter( m_pMesh ); !eiter.end(); eiter ++ )
	{
		M::CEdge * e = *eiter;
		M::CVertex * v1 = m_pMesh->edgeVertex1(e);
		M::CVertex * v2 = m_pMesh->edgeVertex2(e);

	  CPoint2 p1 = v1->huv();
	  CPoint2 p2 = v2->huv();

	  std::complex<double> z1(p1[0], p1[1]);
	  std::complex<double> z2(p2[0], p2[1]);

	  std::complex<double> dz = z2 - z1;

	  std::complex<double> mu1 = v1->mu();
	  std::complex<double> mu2 = v2->mu();
	  std::complex<double> mu  = (mu1 + mu2)/2.0;

	  e->mu_length() = std::abs( dz + mu * std::conj( dz ));
	}

	_calculate_initial_inversive_distance();
};


 //calculate edge weight

template<typename M>
void CQCVirtualInversiveDistanceRicciFlow<M>::_calculate_edge_weight()
{
	_calculate_edge_weight_by_power_center();
}

template<typename M>
void CQCVirtualInversiveDistanceRicciFlow<M>::_calculate_metric( double step_length )
{

	  _set_target_curvature();
	  //_flow( 1e-5 );
	  double error_threshold = 1e-6;
	  //double step_length = 0.5;
	  _Newton( error_threshold, step_length );
		  
};      


template<typename M>
void CQCVirtualInversiveDistanceRicciFlow<M>::_calculate_initial_inversive_distance()
{
	// set initial inversive distance and radius

	for( M::MeshVertexIterator viter( m_pMesh );  !viter.end(); viter ++ )
	{
		M::CVertex * v = *viter;
		v->u()   = INT_MAX;
	};

	for( M::MeshEdgeIterator eiter( m_pMesh ); !eiter.end(); eiter ++ )
	{
		M::CEdge * e = *eiter;
		//e->length() = m_pMesh->edgeLength( e );
		e->length() = e->mu_length();
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
				v->u() = u;
			he = m_pMesh->faceNextCcwHalfEdge( he );
		}
	}

	for( M::MeshVertexIterator viter( m_pMesh );  !viter.end(); viter ++ )
	{
		M::CVertex * v = *viter;
		v->u()   = v->u() - 2.0;
		//v->u()   = v->u();
	};

	for( M::MeshEdgeIterator eiter( m_pMesh ); !eiter.end();  eiter++ )
	{
		M::CEdge * e = *eiter;

		M::CVertex * v1 = m_pMesh->edgeVertex1( e );
		M::CVertex * v2 = m_pMesh->edgeVertex2( e );

		double l = e->mu_length();
		//(v1->point()-v2->point()).norm();
		
		double u1 = v1->u();
		double u2 = v2->u();

		//double r1 = exp(u1);
		//double r2 = exp(u2);
		//e->inversive_distance() = ( l*l + r1*r1 + r2*r2  ) / ( 2*r1*r2 );
		e->inversive_distance() = ( l*l + exp(u1)*exp(u1) + exp(u2)*exp(u2)  ) / ( 2*exp(u1)*exp(u2) );
	};

};



template<typename M>
void CQCVirtualInversiveDistanceRicciFlow<M>::_calculate_power_center( typename M::CFace * head )
{

	M::CHalfEdge * he[3];

	he[0] = m_pMesh->faceMostCcwHalfEdge( head  );
	he[1] = m_pMesh->faceNextCcwHalfEdge( he[0] );
	he[2] = m_pMesh->faceNextCcwHalfEdge( he[1] );

	std::vector<M::CVertex*> av;

	for( int i = 0; i < 3; i ++ )
	{
		av.push_back( m_pMesh->halfedgeSource( he[i] ) );
	}

	M::CEdge * e[3];
	for( int i = 0; i < 3; i ++ )
	{
		e[i] =  m_pMesh->halfedgeEdge( he[i] );
	}
	
	av[0]->uv() = CPoint2(0,0);
	av[1]->uv() = CPoint2( e[0]->length(), 0 );

	CPoint2 c0,c1;

	_circle_circle_intersection( CCircle( av[0]->uv(), e[2]->length() ),
							     CCircle( av[1]->uv(), e[1]->length() ),
								 c0,c1);


	if( cross( av[1]->uv() - av[0]->uv(), c0- av[0]->uv() ) > 0 )
	{
		av[2]->uv() = c0;
	}
	else
	{
		av[2]->uv() = c1;
	}

	CPoint p[3];

	for( int i = 0; i < 3; i ++ )
	{
		p[i] = CPoint( av[i]->uv()[0], av[i]->uv()[1], exp( av[i]->u() ) );
	}

	double C[3];
	for( int i = 0; i < 3; i ++ )
	{
		C[i] = p[i]*p[i];
	}
	
	double A[2][2];

	A[0][0] = 2 * (p[0][0] - p[1][0]);
	A[0][1] = 2 * (p[0][1] - p[1][1]);
	A[1][0] = 2 * (p[0][0] - p[2][0]);
	A[1][1] = 2 * (p[0][1] - p[2][1]);

	double B[2];
	B[0] = C[0] - C[1];
	B[1] = C[0] - C[2];

	double IA[2][2];

	double D = A[0][0] * A[1][1] - A[0][1] * A[1][0];

	IA[0][0] =  A[1][1]/D;
	IA[1][1] =  A[0][0]/D;
	IA[0][1] = -A[0][1]/D;
	IA[1][0] = -A[1][0]/D;

	

	//compute the perpendicular foot from the orthogonal circle center 
	CPoint2 c( IA[0][0] * B[0] + IA[0][1] * B[1], IA[1][0] * B[0] + IA[1][1] * B[1] );
	

	for( int i = 0; i < 3; i ++ )
	{
		CPoint2 t = av[(i+1)%3]->uv();
		CPoint2 s = av[(i+0)%3]->uv();

		CPoint2 d = t-s;
		double lambda = ( c-s )*d/mag2(d);
		CPoint2 ft = t * lambda  + s * (1-lambda);
		
		CPoint2 n(-d[1],d[0]);
		n = n/n.norm();	

		//printf("%f \n", (c-ft)*d);
		double h = (c-ft)*n;	//frac{\partial \theta_s }{\partial u_t } = h/l
		double dt = (t-ft).norm();  //frac{\partial l }{\parital u_t } = dt
		double ds = (s-ft).norm();  //frac{\partial l }{\parital u_s } = ds
		
		he[i]->dtheta_du() = h/e[i]->length();
		he[i]->dl_du() = dt;
		//printf("%f ", h/e[i]->length());
	}
	//printf("\n");
};


 //calculate edge weight

template<typename M>
void CQCVirtualInversiveDistanceRicciFlow<M>::_calculate_edge_weight_by_power_center()
{
	for( M::MeshFaceIterator fiter( m_pMesh ); !fiter.end(); fiter ++ )
	{
		M::CFace * head = *fiter;
		_calculate_power_center( head );
	}

	for( M::MeshEdgeIterator eiter( m_pMesh ); !eiter.end(); eiter ++ )
	  {
		  M::CEdge * e = *eiter;
		  e->weight() = 0.0;
		  M::CHalfEdge * he = m_pMesh->edgeHalfedge( e, 0 );
		  e->weight() += he->dtheta_du();
		  he = m_pMesh->edgeHalfedge( e, 1 );
		  if( he != NULL)
			  e->weight() += he->dtheta_du();	  
	  }
};


//Euclidean cosine law
template<typename M>
double CQCVirtualInversiveDistanceRicciFlow<M>::_cosine_law( double a, double b, double c )
{
          double cs =  ( a*a + b * b  - c * c )/( 2.0 * a * b );
          assert( cs <= 1.0 && cs >= -1.0 );
          return acos( cs );    
};

//normalization

template<typename M>
void CQCVirtualInversiveDistanceRicciFlow<M>::_normalization( Eigen::VectorXd & x, int num )
{

	double s = 0;
	for(int i = 0; i < num; i ++ )
	{
		s += x(i);
	}
	s /= num;

	for (int i = 0; i < num; i++)
	{
	 x(i) -= s;
	}
};


//set target curvature

template<typename M>
void CQCVirtualInversiveDistanceRicciFlow<M>::_set_target_curvature()
{
  for( M::MeshVertexIterator viter( m_pMesh ); !viter.end(); viter ++ )
  {
	  M::CVertex * v = *viter;
    v->target_k() = 0;
  }
/*
  std::vector<M::CLoop*>& pLs = m_boundary.loops();
	
  int id = 0;

  for( std::vector<M::CLoop*>::iterator liter = pLs.begin(); liter != pLs.end(); liter++)
  {
	  M::CLoop * pL = *liter;
	  std::list<H*> & pHes = pL->halfedges();

	  double sum = 0;
	  double inv_sum = 0;
	
	  for( std::list<H*>::iterator hiter = pHes.begin(); hiter != pHes.end(); hiter ++ )
	  {
			H  * he = *hiter;
			E  * pe = m_pMesh->halfedgeEdge( he );

			sum  += pe->length();
			inv_sum += 1/pe->length();
 	  }

	  for( std::list<H*>::iterator hiter = pHes.begin(); hiter != pHes.end(); hiter ++ )
	  {
			H * ce = *hiter;
			V   * pv = m_pMesh->halfedgeTarget( ce );
			H * he = m_pMesh->vertexMostCcwInHalfEdge( pv );
			H * te =  m_pMesh->vertexMostClwOutHalfEdge( pv );

			double L = ( m_pMesh->halfedgeEdge(he)->length() + m_pMesh->halfedgeEdge( te )->length()  )/2.0;

		 // map all boundaries to circular holes
			if( id == 0 )
			{
				double tk = 2*PI*L/sum;
				pv->target_k() = tk;
			}
			else
			{
				double tk = -2*PI*L/sum;
				pv->target_k() = tk;
			}
	  }
	  id ++;
  }
*/
};


 //calculate edge length
template<typename M>
void CQCVirtualInversiveDistanceRicciFlow<M>::_normalize_mesh()
{
	CPoint pmax(-1e+5,-1e+5,-1e+5);
	CPoint pmin( 1e+5, 1e+5, 1e+5);

	for( M::MeshVertexIterator viter( m_pMesh );  !viter.end(); viter ++ )
	{
		M::CVertex * v = *viter;
		for( int i = 0; i < 3; i ++ )
		{
			pmax[i] = (pmax[i]>v->point()[i])?pmax[i]:v->point()[i];
			pmin[i] = (pmin[i]<v->point()[i])?pmin[i]:v->point()[i];
		}
	};

	CPoint d = pmax - pmin;
	double s = d[0];
	s = (s > d[1] )?s:d[1];
	s = (s > d[2] )?s:d[2];
	s = s/2.0;

	for( M::MeshVertexIterator viter( m_pMesh );  !viter.end(); viter ++ )
	{
		M::CVertex * v = *viter;
		v->point() = (v->point() - d/2 )/s;
	};

	//we also normalize the huv
	CPoint2 center(0,0);
	for( M::MeshVertexIterator viter( m_pMesh );  !viter.end(); viter ++ )
	{
		M::CVertex * v = *viter;
		center = center + v->huv();
	};
	
	center = center/m_pMesh->numVertices();

	for( M::MeshVertexIterator viter( m_pMesh );  !viter.end(); viter ++ )
	{
		M::CVertex * v = *viter;
		v->huv() = v->huv() - center;
	};

	double ds = -1;

	for( M::MeshVertexIterator viter( m_pMesh );  !viter.end(); viter ++ )
	{
		M::CVertex * v = *viter;
		for( int i = 0; i < 2; i ++ )
		{
			ds=(ds> fabs(v->huv()[i]))?ds: fabs(v->huv()[i]);
		}
	};

	for( M::MeshVertexIterator viter( m_pMesh );  !viter.end(); viter ++ )
	{
		M::CVertex * v = *viter;
		v->huv() /= ds;
	};


};


 //calculate edge length
template<typename M>
void CQCVirtualInversiveDistanceRicciFlow<M>::_length( double u1, double u2, typename M::CEdge * e )
{
	//double r1 = exp( u1 );
	//double r2 = exp( u2 );

	double inv_dis = e->inversive_distance();
	//assert(  2*r1*r2*inv_dis -r1*r1 - r2*r2 > 0 );
	//e->length() = sqrt( 2*r1*r2*inv_dis -r1*r1 - r2*r2 );
	e->length() = sqrt( 2*exp(u1)*exp(u2)*inv_dis -exp(u1)*exp(u1) - exp(u2)*exp(u2) );
	
};

} //namespace RicciFlow
} //namespace MeshLib

#endif  _QC_VIRTUAL_INVERSIVE_DISTANCE_RICCI_FLOW_H_