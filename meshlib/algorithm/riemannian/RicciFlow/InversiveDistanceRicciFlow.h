/*! \file InversiveDistanceRicciFlow.h
 * \brief Euclidean inversive distance Ricci flow
 *  \author David Gu
 *  \date   documented on 10/17/2010
 *
 *	Algorithm for Inversive distance  Ricci Flow
 */

#ifndef _INVERSIVE_DISTANCE_RICCI_FLOW_H_
#define _INVERSIVE_DISTANCE_RICCI_FLOW_H_

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

	/*! \brief CInversiveDistanceRicciFlow class
	 *  
	 *  Algorithm for Inversive distance Ricci flow
	 */
	 
 template<typename M>
 class CInversiveDistanceRicciFlow : public CBaseRicciFlow<M>
  {
  public:
	  /*! \brief CInversiveDistanceRicciFlow constructor
	   *
	   *  call base class constructor 
	   */
	  CInversiveDistanceRicciFlow( M * pMesh );
	  /*! override the virtual function calculate_metric */
	  virtual void _calculate_metric();
	

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
	 *	Calculate the edge weight
	 */
    void _calculate_edge_weight();

	/*!
	 *	Set the target curvature on each vertex
	 */
    virtual void    _set_target_curvature();


	/*! compute initial vertex radius and edge length, edge inversive distance */
	void _calculate_initial_inversive_distance();
	/*! calculate edge weight by power center */
	void _calculate_edge_weight_by_power_center();
	/*! calculate power center on a face head
	    \param head input face
	*/
	void _calculate_power_center( typename M::CFace* head );
  };




template<typename M>
CInversiveDistanceRicciFlow<M>::CInversiveDistanceRicciFlow( M * pMesh ): CBaseRicciFlow( pMesh)
{
	_calculate_initial_inversive_distance();
}

 //calculate edge length
template<typename M>
void CInversiveDistanceRicciFlow<M>::_length( double u1, double u2, typename M::CEdge * e )
{
	  double r1 = exp( u1 );
	  double r2 = exp( u2 );

	  double inv_dis = e->inversive_distance();

	  e->length() = sqrt( r1*r1 + r2*r2 + 2*r1*r2*inv_dis );
};

 //calculate edge weight

template<typename M>
void CInversiveDistanceRicciFlow<M>::_calculate_edge_weight()
{
	_calculate_edge_weight_by_power_center();
	return;

	for( M::MeshEdgeIterator eiter( m_pMesh ); !eiter.end(); eiter ++ )
	  {
		  M::CEdge * e = *eiter;
		  e->weight() = 0.0;
	  }

	for(  M::MeshFaceIterator fiter( m_pMesh ) ; !fiter.end(); fiter ++ )
	{
		M::CFace * f = *fiter;

		double r[3];
		double l[3];
		M::CHalfEdge * he = m_pMesh->faceMostCcwHalfEdge( f );
		for( int i = 0; i < 3; i ++ )
		{
			M::CVertex * v = m_pMesh->halfedgeSource( he );
			r[i] = exp( v->u() );
			M::CEdge * e = m_pMesh->halfedgeEdge( he );
			l[(i+2)%3] = e->length();
			he = m_pMesh->faceNextCcwHalfEdge( he );
		}

		double L1 = l[0];
		double L2 = l[1];
		double L3 = l[2];
		double r1 = r[0];
		double r2 = r[1];
		double r3 = r[2];
		
		double face_area = sqrt( (L1+L2+L3)*(L2+L3-L1)*(L3+L1-L2)*(L1+L2-L3) )/4.0;
		double ratio = -1.0/(face_area*2) ;
		
		double a = L1*L1;
		double b = L2*L2;
		double c = L3*L3;

		double x = (r2*r2-r3*r3)/a;
		double y = (r3*r3-r1*r1)/b;
		double z = (r1*r1-r2*r2)/c;

		double deri[3];
		deri[0] = ( c-a-b - 2*a*x - (c+a-b)*z )/4 * ratio;
		deri[1] = ( a-b-c - 2*b*y - (a+b-c)*x )/4 * ratio;
		deri[2] = ( b-c-a - 2*c*z - (b+c-a)*y )/4 * ratio;

		he = m_pMesh->faceMostCcwHalfEdge( f );
		for(int i = 0; i < 3; i ++ )
		{
			M::CEdge * e = m_pMesh->halfedgeEdge( he );
			e->weight() += deri[i];
			he = m_pMesh->faceNextCcwHalfEdge( he );
		}
	}
}

template<typename M>
void CInversiveDistanceRicciFlow<M>::_calculate_metric()
{

	  _set_target_curvature();
	  double error_threshold = 1e-6;
	  double step_length = 0.5;
	  _Newton( error_threshold, step_length );
		  
};      


template<typename M>
void CInversiveDistanceRicciFlow<M>::_calculate_initial_inversive_distance()
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
			double u = log( r );
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

		double r1 = exp( u1 );
		double r2 = exp( u2 );
		e->inversive_distance() = ( l*l - r1*r1 - r2*r2  ) / ( 2*r1*r2 );
		assert( e->inversive_distance() >= 0.9 );
	}

};



template<typename M>
void CInversiveDistanceRicciFlow<M>::_calculate_power_center( typename M::CFace* head )
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

	CCircle C[3];

	for( int i = 0; i < 3; i ++ )
	{
		C[i] = CCircle( av[i]->uv(), exp( av[i]->u() ) );
		//printf("Center (%f,%f) radius %f\n", av[i]->huv()[0], av[i]->huv()[1], exp( av[i]->u()) );
	}
	CCircle power = orthogonal( C );
	
	//compute the perpendicular foot from the orthogonal circle center 
	CPoint2 c = power.c();
	

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
		//double ds = (s-ft).norm();  //frac{\partial l }{\parital u_s } = ds
		
		he[i]->dtheta_du() = h/e[i]->length();
		he[i]->dl_du() = dt;
		//printf("%f ", h/e[i]->length());
	}
	//printf("\n");
};


 //calculate edge weight

template<typename M>
void CInversiveDistanceRicciFlow<M>::_calculate_edge_weight_by_power_center()
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
double CInversiveDistanceRicciFlow<M>::_cosine_law( double a, double b, double c )
{
          double cs =  ( a*a + b * b  - c * c )/( 2.0 * a * b );
          assert( cs <= 1.0 && cs >= -1.0 );
          return acos( cs );    
};



//set target curvature

template<typename M>
void CInversiveDistanceRicciFlow<M>::_set_target_curvature()
{
  for( M::MeshVertexIterator viter( m_pMesh ); !viter.end(); viter ++ )
  {
	  M::CVertex * v = *viter;
	  v->target_k() = 0;
  }
};

}	//namespace RicciFlow
}	//namespace MeshLib

#endif  _INVERSIVE_DISTANCE_RICCI_FLOW_H_