/*! \file GenusOneInversiveDistanceRicciFlow.h
 * \brief Euclidean inversive distance Ricci flow for genus one surface with boundaries
 *  \author David Gu
 *  \date   documented on 10/17/2010
 *
 *	Algorithm for Inversive distance  Ricci Flow
 */

#ifndef _GENUS_ONE_INVERSIVE_DISTANCE_RICCI_FLOW_H_
#define _GENUS_ONE__INVERSIVE_DISTANCE_RICCI_FLOW_H_

#include "InversiveDistanceRicciFlow.h"


namespace MeshLib
{
namespace RicciFlow
{
	/*! \brief CGenusOneInversiveDistanceRicciFlow class
	 *  
	 *  Algorithm for Genus One Inversive distance Ricci flow
	 */
	 
 template<typename M>
 class CGenusOneInversiveDistanceRicciFlow : public CInversiveDistanceRicciFlow<M>
  {
  public:
	  /*! \brief CGenusOneInversiveDistanceRicciFlow constructor
	   *
	   *  call base class constructor 
	   */
	  CGenusOneInversiveDistanceRicciFlow( M * pMesh );
	  /*! override the virtual function calculate_metric */
	  void _calculate_metric();
	

  protected:

	/*!
	 *	Set the target curvature on each vertex
	 */
    void    _set_target_curvature();
	/*!
	 *	Gradient flow method
	 */
	bool _flow( double error_threshold );
	/*!
	 *	the difference between the old target curvature and the current target curvature
	 */
	double m_difference;
  };

template<typename M>
CGenusOneInversiveDistanceRicciFlow<M>::CGenusOneInversiveDistanceRicciFlow( M * pMesh ): CInversiveDistanceRicciFlow( pMesh)
{
}

template<typename M>
void CGenusOneInversiveDistanceRicciFlow<M>::_calculate_metric()
{
	  double error_threshold = 1e-6;
	  double step_length = 0.5;

	  _calculate_edge_length();
	  _set_target_curvature();
	  
	  //for( int k = 0; k < 8; k ++ )
	  int k = 0;
	  while( true )
	  {
		printf("\n\n Round %d, the target curvature difference is %f \n\n", k ++, m_difference );
	
		_Newton( error_threshold, step_length );
		_set_target_curvature();
		//if( m_difference < 2e-2 ) break;
	  	if( _flow( error_threshold ) ) break;

	  }      	  
};      


//set target curvature

template<typename M>
void CGenusOneInversiveDistanceRicciFlow<M>::_set_target_curvature()
{

  for( M::MeshVertexIterator viter( m_pMesh ); !viter.end(); viter ++ )
  {
	  M::CVertex * v = *viter;
	  if( m_pMesh->isBoundary( v ) ) continue;
	  v->target_k() = 0;
  }

  std::vector<M::CLoop*>& pLs = m_boundary.loops();
	
  m_difference = -1e+3;

  for( std::vector<M::CLoop*>::iterator liter = pLs.begin(); liter != pLs.end(); liter++)
  {
	  M::CLoop * pL = *liter;
	  std::list<M::CHalfEdge*> & pHes = pL->halfedges();

	  double sum = 0;

	  for( std::list<M::CHalfEdge*>::iterator hiter = pHes.begin(); hiter != pHes.end(); hiter ++ )
	  {
		  M::CHalfEdge  * he = *hiter;
		  M::CEdge  * pe = m_pMesh->halfedgeEdge( he );
			sum  += pe->length();
 	  }
	  

	  for( std::list<M::CHalfEdge*>::iterator hiter = pHes.begin(); hiter != pHes.end(); hiter ++ )
	  {
		  M::CHalfEdge * ce = *hiter;
		  M::CVertex   * pv = m_pMesh->halfedgeTarget( ce );

		  M::CHalfEdge * he = m_pMesh->vertexMostCcwInHalfEdge( pv );
		  M::CHalfEdge * te =  m_pMesh->vertexMostClwOutHalfEdge( pv );

		  M::CEdge * pE = m_pMesh->halfedgeEdge(he);
		  M::CEdge * nE = m_pMesh->halfedgeEdge( te );

			double L = ( pE->length() + nE->length()  )/2.0;
			double tk = -2*PI*L/sum;

			double difference = fabs( pv->target_k() - tk );
			pv->target_k() = tk;
			
			m_difference =( difference > m_difference )? difference:m_difference;
			
	  }

  }

  printf("Set_target_curvature: %f\n", m_difference );

};

template<typename M>
bool CGenusOneInversiveDistanceRicciFlow<M>::_flow( double error_threshold )
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
            v->u() += dif * 5e-2;
		  }
    }
    return false;
};

}	//namespace RicciFlow
}	//namespace MeshLib

#endif  _GENUS_ONE_INVERSIVE_DISTANCE_RICCI_FLOW_H_