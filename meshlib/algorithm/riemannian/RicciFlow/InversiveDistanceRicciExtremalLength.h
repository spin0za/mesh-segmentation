/*! \file InversiveDistanceRicciExtremalLength.h
 * \brief Euclidean Inversive Distance Ricci flow for Extremal Length
 *  \author David Gu
 *  \date   documented on 11/15/2010
 *
 *	Algorithm for Inversive Distance Ricci Flow for Extremal Length
 */

#ifndef _INVERSIVE_DISTANCE_RICCI_FLOW_EXTREMAL_LENGTH_H_
#define _INVERSIVE_DISTANCE_RICCI_FLOW_EXTREMAL_LENGTH_H_

#include <map>
#include <vector>

#include "InversiveDistanceRicciFlow.h"

namespace MeshLib
{

namespace RicciFlow
{

	/*! \brief CInversiveDistanceRicciFlowExtremalLength class
	 *  
	 *  Algorithm for Euclidean Inversive Distance Ricci flow for Extremal Length
	 */
	 
 template<typename M>
 class CInversiveDistanceRicciFlowExtremalLength : public CInversiveDistanceRicciFlow<M>
  {
  public:
	  /*! \brief CYamabeExtremalLength constructor
	   *
	   *  call base class constructor 
	   */
	  CInversiveDistanceRicciFlowExtremalLength( M * pMesh );

  protected:

	/*!
	 *	Set the target curvature on each vertex
	 */
    void    _set_target_curvature();
	/*!
	 *	compute the vertex topological valence
	 */
	void _calculate_topo_valence();
  };

template<typename M>
CInversiveDistanceRicciFlowExtremalLength<M>::CInversiveDistanceRicciFlowExtremalLength( M * pMesh ): CInversiveDistanceRicciFlow( pMesh)
{
	_calculate_topo_valence();
}

template<typename M>
void CInversiveDistanceRicciFlowExtremalLength<M>::_calculate_topo_valence()
{
	for( M::MeshVertexIterator viter( m_pMesh ); !viter.end(); viter ++ )
    {
		M::CVertex * v = *viter;
      v->valence() = 0;
      for(M::VertexEdgeIterator veiter( v ); !veiter.end(); ++ veiter )
      {
		  M::CEdge * e = *veiter;
        if( e->sharp() ) v->valence() ++;
      }
    }
}


//set target curvature

template<typename M>
void CInversiveDistanceRicciFlowExtremalLength<M>::_set_target_curvature()
{
	int count = 0;
  for( M::MeshVertexIterator viter( m_pMesh ); !viter.end(); viter ++ )
  {
	  M::CVertex * v = *viter;
    v->target_k() = 0;
	
	if( v->valence() > 2 )
	{
		count ++;
		v->target_k() = PI/2.0;
	}
  }
	printf("There are %d corners.\n", count );
};

} //namespace RicciFlow
} //namespace MeshLib

#endif  