/*! \file VirtualInversiveDistanceRicciExtremalLength.h
 * \brief Euclidean Inversive Distance Ricci flow for Extremal Length with virtual raidus
 *  \author David Gu
 *  \date   documented on 11/15/2010
 *
 *	Algorithm for Inversive Distance Ricci Flow for Extremal Length with virtual radius
 */

#ifndef _VIRTUAL_INVERSIVE_DISTANCE_RICCI_FLOW_EXTREMAL_LENGTH_H_
#define _VIRTUAL_INVERSIVE_DISTANCE_RICCI_FLOW_EXTREMAL_LENGTH_H_

#include <map>
#include <vector>

#include "VirtualInversiveDistanceRicciFlow.h"

namespace MeshLib
{
namespace RicciFlow
{
	/*! \brief CVirtualInversiveDistanceRicciFlowExtremalLength class
	 *  
	 *  Algorithm for Euclidean Inversive Distance Ricci flow for Extremal Length with virtual radius
	 */
	 
 template<typename M>
 class CVirtualInversiveDistanceRicciFlowExtremalLength : public CVirtualInversiveDistanceRicciFlow<M>
  {
  public:
	  /*! \brief CVirtualInversiveDistanceRicciFlowExtremalLength constructor
	   *
	   *  call base class constructor 
	   */
	  CVirtualInversiveDistanceRicciFlowExtremalLength( M * pMesh );

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
CVirtualInversiveDistanceRicciFlowExtremalLength<M>::CVirtualInversiveDistanceRicciFlowExtremalLength( M * pMesh ): CVirtualInversiveDistanceRicciFlow( pMesh)
{
	_calculate_topo_valence();
}

template<typename M>
void CVirtualInversiveDistanceRicciFlowExtremalLength<M>::_calculate_topo_valence()
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
void CVirtualInversiveDistanceRicciFlowExtremalLength<M>::_set_target_curvature()
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