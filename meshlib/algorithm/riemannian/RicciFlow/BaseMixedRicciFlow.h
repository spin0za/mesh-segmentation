/*! \file BaseMixedRicciFlow.h
 *  \brief Base class for Mixed Ricci flow algorithm: interior vertices target curvatre are set, boundary conformal factors are set
 *  
 *  \author David Gu
 *  \date   documented on 08/03/2011
 *
 *	Algorithm for general Mixed Ricci Flow
 */

#ifndef _BASE_MIXED_RICCI_FLOW_H_
#define _BASE_MIXED_RICCI_FLOW_H_

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

namespace MeshLib
{

namespace RicciFlow
{
/*! \brief BaseClass CBaseRicciFlow
*
*	Algorithm for computing general Ricci flow
*/
template<typename M>
class CBaseMixedRicciFlow : public CBaseRicciFlow<M>
  {
  public:
    /*! \brief CBaseMixedRicciFlow constructor
	 *  \param pMesh the input mesh
	 */
    CBaseMixedRicciFlow( M * pMesh );
    /*! \brief CBaseRicciFlow destructor
	 */
    ~CBaseMixedRicciFlow();

  protected:
	  /*!
	   *	Calculate each edge length, has to be defined in the derivated classes
	   */
	  virtual void _length( double u1, double u2, typename M::CEdge * e )=0;

	/*!
	 *	Cosine law, has to be defined in the derivated classes
	 */
	virtual double _cosine_law( double a, double b, double c ) = 0;
	/*!
	 *	Calculate corner angle
	 */
    void _calculate_corner_angle();
	/*!
	 *	Calculate vertex curvature error
	 */
   double   _calculate_curvature_error();

	/*!
	 *	Calculate the edge weight
	 */
    virtual void _calculate_edge_weight() = 0;

	/*!
	 *	Set the target curvature on each interior vertex
	 */
    virtual void  _set_target_curvature() =  0;
	/*!
	 *	Set the target conformal factor on each boundary vertex
	 */
    virtual void  _set_target_conformal_factor() =  0;
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
	 /*!
	  *	Normalization
	  * \param du the du vector
	  * \param n dimension of the du vector
	  */
	 virtual void _normalization( Eigen::VectorXd & du, int n ) = 0;
	 /*!
	  *	calculate hessian matrix Hessain 
	  * \param SparseMatrix
	  */
	 virtual void _calculate_Hessain( Eigen::SparseMatrix<double> & M );

	 int m_num_vertices;
  };

//Constructor
template<typename M>
CBaseMixedRicciFlow<M>::CBaseMixedRicciFlow( M * pMesh ): CBaseRicciFlow<M>( pMesh )
{
  int idx = 0;
  for( M::MeshVertexIterator viter( m_pMesh ); !viter.end(); viter ++ )
	{
		M::CVertex * v = *viter;
	  if( m_pMesh->isBoundary( v ) ) continue;
      v->u() = 0;
      v->idx() = idx ++; 
	}
  m_num_vertices = idx;
};

//Destructor
template<typename M>
CBaseMixedRicciFlow<M>::~CBaseMixedRicciFlow( )
{
};

//Calculate corner angle
template<typename M>
void CBaseMixedRicciFlow<M>::_calculate_corner_angle( )
{
	//CStructure<M, V,E,F,H> pS( m_pMesh );
	//pS._metric_2_angle();
	
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


//compute curvature error

template<typename M>
double CBaseMixedRicciFlow<M>::_calculate_curvature_error()
{
  double max_error = -1;
  M::CVertex * vert = NULL;

  for(M::MeshVertexIterator viter( m_pMesh); !viter.end() ; viter ++ )
   {
	   M::CVertex * v = *viter;
	   if( m_pMesh->isBoundary( v ) ) continue;
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


//gradient flow method

template<typename M>
bool CBaseMixedRicciFlow<M>::_flow( double error_threshold )
{
  int num = m_pMesh->numVertices();

  for( int k = 0; k < 64000; k ++  )
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

//Newton's method for optimizing entropy energy

template<typename M>
void CBaseMixedRicciFlow<M>::_Newton( double threshold, double step_length )
{
	int num = m_num_vertices;
	double* b = new double[num];
	assert( b );
	memset(b,0,sizeof(double)*num);

	double* x = new double[num];
	assert( x != NULL );
	memset(x,0,sizeof(double)*num);


  	while( true )
	{
		//the order of the following functions really matters

 		_calculate_edge_length();
		_calculate_corner_angle();
		_calculate_vertex_curvature();
		_calculate_edge_weight();

		double error =  _calculate_curvature_error();
		printf("Current error is %f\r\n", error );

		if( error < threshold) break;
	  
		Eigen::SparseMatrix<double> M(num,num );
		M.setZero();


		//set b 
		Eigen::VectorXd b(num);
		for( M::MeshVertexIterator viter( m_pMesh ); !viter.end(); viter ++  )
		{
			M::CVertex * v = *viter;
			if( m_pMesh->isBoundary( v ) ) continue;

			int idx = v->idx();
			b(idx) = v->target_k() - v->k();
		}

		_calculate_Hessain( M );

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

		_normalization( x, num );

		for( M::MeshVertexIterator viter( m_pMesh ); !viter.end(); viter ++  )
		{
			M::CVertex * v = *viter;
			if( m_pMesh->isBoundary( v ) ) continue;
			int idx = v->idx();
			v->u() += x(idx) * step_length;
		}

  }

};


//Newton's method for optimizing entropy energy

template<typename M>
void CBaseMixedRicciFlow<M>::_calculate_Hessain( Eigen::SparseMatrix<double> & M )
{
	std::vector<Eigen::Triplet<double> > _coefficients;

	//set A
	for( M::MeshEdgeIterator eiter( m_pMesh ); !eiter.end(); eiter ++  )
	{
		M::CEdge * e = *eiter;
		M::CVertex * v1 = m_pMesh->edgeVertex1( e );
		M::CVertex * v2 = m_pMesh->edgeVertex2( e );
	  
	  if( m_pMesh->isBoundary( v1 ) || m_pMesh->isBoundary( v2 ) ) continue;

	  //pMatrix->AddElementTail(v1->idx(), v2->idx(), -e->weight() );
	  //pMatrix->AddElementTail(v2->idx(), v1->idx(), -e->weight() );

	  _coefficients.push_back( Eigen::Triplet<double>( v1->idx(), v2->idx(), -e->weight() ) );
	  _coefficients.push_back( Eigen::Triplet<double>( v2->idx(), v1->idx(), -e->weight() ) );
	}
  
	for( M::MeshVertexIterator viter( m_pMesh ); !viter.end(); viter ++  )
	{
		M::CVertex * v = *viter;
		if( m_pMesh->isBoundary( v ) ) continue;

		int idx = v->idx();
		double w = 0;
		for( M::VertexEdgeIterator veiter( v ); !veiter.end(); veiter ++ )
		{
			M::CEdge * pE = *veiter;
		  
			M::CVertex * v1 = m_pMesh->edgeVertex1( pE );
			M::CVertex * v2 = m_pMesh->edgeVertex2( pE );

		  if( m_pMesh->isBoundary( v1 ) || m_pMesh->isBoundary( v2 ) ) continue;

		  w += pE->weight();
		}
		//pMatrix->AddElementTail( idx, idx, w); 
		_coefficients.push_back( Eigen::Triplet<double>( idx, idx, w ) );
	}

	M.setFromTriplets( _coefficients.begin(), _coefficients.end() );
};

} //namespace Ricci Flow
} //namespace MeshLib

#endif  _BASE_MIXED_RICCI_FLOW_H_