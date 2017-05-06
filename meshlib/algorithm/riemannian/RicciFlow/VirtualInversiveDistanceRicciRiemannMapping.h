/*! \file InversiveDistanceRicciRiemannMapping.h
 * \brief Euclidean Inversive Distance Ricci flow for Riemann Mapping
 *  \author David Gu
 *  \date   documented on 11/15/2010
 *
 *	Algorithm for Euclidean Inversive Distance Ricci flow for Riemann Mapping
 */

#ifndef _VIRTURAL_INVERSIVE_DISTANCE_RICCI_FLOW_RIEMANN_MAPPING_H_
#define _VIRTURAL_INVERSIVE_DISTANCE_RICCI_FLOW_RIEMANN_MAPPING_H_

#include <map>
#include <vector>

#include "VirtualInversiveDistanceRicciFlow.h"

namespace MeshLib
{
namespace RicciFlow
{
	/*! \brief CVirtualInversiveDistanceRicciFlowRiemannMapping class
	 *  
	 *  Algorithm for Virtual Inversive Distance RicciFlow with virutal radius for Riemann Mapping
	 */
	 
 template<typename M>
 class CVirtualInversiveDistanceRicciFlowRiemannMapping : public CVirtualInversiveDistanceRicciFlow<M>
  {
  public:
	  /*! \brief CVirtualInversiveDistanceRicciFlowRiemannMapping constructor
	   *
	   *  call base class constructor 
	   */
	  CVirtualInversiveDistanceRicciFlowRiemannMapping( M * pMesh );

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
CVirtualInversiveDistanceRicciFlowRiemannMapping<M>::CVirtualInversiveDistanceRicciFlowRiemannMapping( M * pMesh ): CVirtualInversiveDistanceRicciFlow( pMesh)
{
	//_calculate_topo_valence();
}

template<typename M>
void CVirtualInversiveDistanceRicciFlowRiemannMapping<M>::_calculate_topo_valence()
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
void CVirtualInversiveDistanceRicciFlowRiemannMapping<M>::_set_target_curvature()
{
  for( M::MeshVertexIterator viter( m_pMesh ); !viter.end(); viter ++ )
  {
	  M::CVertex * v = *viter;
    v->target_k() = 0;
  }
};

} //namespace RicciFlow
} //namespace MeshLib

#endif  