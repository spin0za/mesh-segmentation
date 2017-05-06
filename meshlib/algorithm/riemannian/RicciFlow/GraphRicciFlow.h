/*! \file GraphRicciFlow.h
 *  \brief General Euclidean Ricci flow algorithm for Graph Visualization
 *  \author David Gu
 *  \date   documented on 03/23/2011
 *
 *	Algorithm for general Ricci Flow for graph visualization
 */

#ifndef _GRAPH_RICCI_FLOW_H_
#define _GRAPH_RICCI_FLOW_H_

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
#include "Operator/Operator.h"
#include "GraphRicciFlowMesh.h"


namespace MeshLib
{

namespace RicciFlow
{

/*! \brief Class CGraphRicciFlow
*
*	Algorithm for computing Ricci flow
*/
template<typename M>
  class CGraphRicciFlow
  {
  public:
    /*! \brief CRicciFlow constructor
	 *  \param pMesh the input mesh
	 */
    CGraphRicciFlow( M * pMesh );
    /*! \brief CRicciFlow destructor
	 */
    ~CGraphRicciFlow();
	/*!	Computing the metric
	 */
	virtual void _calculate_metric();


  protected:
    /*!
     *	the input mesh
	 */
    M	  * m_pMesh;
	/*!
	 *	boundary of the input mesh
	 */
	typename M::CBoundary		  m_boundary;

  protected:
	  /*!
	   *	Calculate each edge length
	   */
    virtual void _calculate_edge_length();
	/*!
	 *	Calculate corner angle
	 */
    virtual void _calculate_corner_angle();
	/*!
	 *	Calculate vertex curvature
	 */
    void _calculate_vertex_curvature();
	/*!
	 *	Calculate the edge weight
	 */
    virtual void _calculate_edge_weight();

	/*!
	 *	Set the target curvature on each vertex
	 */
    void    _set_target_curvature();
	/*!
	 *	Calculate vertex curvature error
	 */
   double   _calculate_curvature_error();
    /*!
	 *	Curvature flow 
	 */
     bool   _flow( double );
	 /*!
	  *	Newton's method to optimize the entropy energy
	  */
	 void   _Newton(std::vector<typename M::CVertex* > &, double );
  };

//Constructor
template<typename M>
CGraphRicciFlow<M>::CGraphRicciFlow( M * pMesh ): m_pMesh( pMesh), m_boundary( pMesh )
{
};

//Destructor
template<typename M>
CGraphRicciFlow<M>::~CGraphRicciFlow( )
{
};

//Compute the edge length
template<typename M>
void CGraphRicciFlow<M>::_calculate_edge_length( )
{

	for( M::MeshEdgeIterator eiter(m_pMesh);  !eiter.end();  eiter++ )
  {
	  M::CEdge * e = *eiter;

	  M::CVertex * v1 = m_pMesh->edgeVertex1( e );
	  M::CVertex * v2 = m_pMesh->edgeVertex2( e );

      double u1 = v1->u();
      double u2 = v2->u();

	  e->length() = exp(u1) + exp(u2);
  }

};


//Calculate corner angle
template<typename M>
void CGraphRicciFlow<M>::_calculate_corner_angle( )
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
          double cs =  ( l[(i+1)%3] * l[(i+1)%3] + l[(i+2)%3] * l[(i+2)%3]  - l[i] * l[i] )/( 2.0 * l[(i+1)%3] * l[(i+2)%3] );
          assert( cs <= 1.0 && cs >= -1.0 );
          he[(i+1)%3]->angle() = acos( cs );
      }

  }
    
};

//Calculate vertex curvature
template<typename M>
void CGraphRicciFlow<M>::_calculate_vertex_curvature()
{
  double total_k = 0;

  for( M::MeshVertexIterator viter( m_pMesh ); !viter.end(); viter ++ )
  {
      V * v = *viter;
      double k  = (v->boundary() )? PI: PI * 2;
	  for( M::VertexInHalfedgeIterator vh( m_pMesh, v ); !vh.end();  ++vh )
      {
		  M::CHalfEdge * he = *vh;
          k -= he->angle();
      }
      v->k() = k;
      total_k += v->k();
  }

};

//Calculate edge weight

template<typename M>
void CGraphRicciFlow<M>::_calculate_edge_weight()
{
	for( M::MeshEdgeIterator eiter( m_pMesh ); !eiter.end(); eiter ++ )
  {
	  M::CEdge * e = *eiter;
      e->weight() = 0.0;
  }

	for(  M::MeshFaceIterator fiter( m_pMesh ) ; !fiter.end(); fiter ++ )
  {
	  M::CFace * f = *fiter;

      double r[3];
      int i = 0;
	  for( M::FaceHalfedgeIterator hiter(  f ); !hiter.end(); ++hiter )
      {
		  M::CHalfEdge * he = *hiter;
		  M::CVertex   * v = m_pMesh->halfedgeTarget(he);
         r[i++] = exp( v->u() );
      }

      double w = sqrt(r[0]*r[1]*r[2]/(r[0]+r[1]+r[2]));
      
	  for( M::FaceEdgeIterator eiter(f); !eiter.end(); ++ eiter )
      {
		  M::CEdge * e = * eiter;
          e->weight() += w/e->length();
      }
  }
};

//set target curvature

template<typename M>
void CGraphRicciFlow<M>::_set_target_curvature()
{
  int euler = m_pMesh->numVertices() + m_pMesh->numFaces() - m_pMesh->numEdges();
  for( M::MeshVertexIterator viter( m_pMesh ); !viter.end(); viter ++ )
  {
	  M::CVertex * v = *viter;
	  v->target_k() = euler * 2 * PI/(double)(m_pMesh->numVertices());
  }


};

//compute curvature error

template<typename M>
double CGraphRicciFlow<M>::_calculate_curvature_error()
{
  double max_error = -1;
  M::CVertex * vert = NULL;

  for( M::MeshVertexIterator viter( m_pMesh); !viter.end() ; viter ++ )
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
void CGraphRicciFlow<M>::_calculate_metric()
{

	std::vector<M::CVertex*> vertices;

  int idx = 0;
  for( M::MeshVertexIterator viter( m_pMesh ); !viter.end(); viter ++ )
	{
		M::CVertex * v = *viter;
		vertices.push_back(  v );
		v->u() = 0;
		v->idx() = idx ++; 
	}

  //double error = 1e-6;
  double error = 5e-4;

  _calculate_edge_length();
  _set_target_curvature();
  _Newton( vertices, error );
};      



//gradient flow method

template<typename M>
bool CGraphRicciFlow<M>::_flow( double error_threshold )
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
		    V * v = *viter;
			double dif = v->target_k() - v->k();
            v->u() += dif * 2e-2;
		  }
    }
    return false;
};

//Newton's method for optimizing entropy energy

template<typename M>
void CGraphRicciFlow<M>::_Newton( std::vector<typename M::CVertex* > & vertices, double error_threshold )
{
  int num = m_pMesh->numVertices();

  	while( true )
	{
 		_calculate_edge_length();
		_calculate_edge_weight();

		_calculate_corner_angle();
		_calculate_vertex_curvature();

		double error =  _calculate_curvature_error();
		printf("Current error is %f\r\n", error );

		if( error < error_threshold) break;
	  

		Eigen::VectorXd b(num);
		//set b 
		for(int i = 0; i < num; i ++ )
		{
		    V * v = vertices[i];
			b(i) = v->target_k() - v->k();
		}


		std::vector<Eigen::Triplet<double> > _coefficients;

		for( std::list<E*>::iterator eiter = m_pMesh->edges().begin(); eiter != m_pMesh->edges().end(); eiter ++ )
		{
		  E * e = *eiter;
		  V * v1 = m_pMesh->edgeVertex1( e );
		  V * v2 = m_pMesh->edgeVertex2( e );

		  _coefficients.push_back( Eigen::Triplet<double>( v1->idx(), v2->idx(), -e->weight() ) );
		  _coefficients.push_back( Eigen::Triplet<double>( v2->idx(), v1->idx(), -e->weight() ) );
		
		  //pMatrix->AddElementTail(v1->idx(), v2->idx(), -e->weight() );
		  //pMatrix->AddElementTail(v2->idx(), v1->idx(), -e->weight() );
		}
  

	  for(int i = 0; i < ( int ) vertices.size(); i ++ )
	  {
		  V * v = vertices[i];

		  double w = 0;
	      
		  for( CRicciFlowMesh<V,E,F,H>::VertexEdgeIterator veiter( v ); !veiter.end(); veiter ++ )
		  {
			  E * pE = *veiter;
			  w += pE->weight();
		  }
		   //pMatrix->AddElementTail( i, i, w); 
		   _coefficients.push_back( Eigen::Triplet<double>( i, i, w ) );
	  }

		Eigen::SparseMatrix<double> M(num,num);
		M.setZero();

		M.setFromTriplets( _coefficients.begin(), _coefficients.end() );

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

		//normalize x
		double s = 0;
		for(int i = 0; i < num; i ++ )
		{
			s += x(i);
		}
		s /= num;
		for(int i = 0; i < num; i ++ )
		{
			x(i) -= s;
		}
	

		for (int i = 0; i < num; i++)
		{
		  V * v = vertices[i];

		  assert( v->idx() == i ); 
		  v->u() += x[i];
		}

  }
};

} //namespace RicciFlow
} //namespace MeshLib

#endif  _GRAPH_RICCI_FLOW_H_