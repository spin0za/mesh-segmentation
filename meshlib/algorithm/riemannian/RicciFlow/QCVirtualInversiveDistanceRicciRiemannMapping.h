/*! \file QCInversiveDistanceRicciRiemannMapping.h
 * \brief Quasi-Conformal Euclidean Inversive Distance Ricci flow for Riemann Mapping
 *  \author David Gu
 *  \date   documented on 03/09/2011
 *
 *	Algorithm for Quasi-Conformal Euclidean Inversive Distance Ricci flow for Riemann Mapping
 */

#ifndef _QC_VIRTURAL_INVERSIVE_DISTANCE_RICCI_FLOW_RIEMANN_MAPPING_H_
#define _QC_VIRTURAL_INVERSIVE_DISTANCE_RICCI_FLOW_RIEMANN_MAPPING_H_

#include <map>
#include <vector>

#include "QCVirtualInversiveDistanceRicciFlow.h"

namespace MeshLib
{

namespace RicciFlow
{
	/*! \brief CQCVirtualInversiveDistanceRicciFlowRiemannMapping class
	 *  
	 *  Algorithm for Quasi-Conformal Inversive Distance RicciFlow with virutal radius for Riemann Mapping
	 */
	 
 template<typename M>
 class CQCVirtualInversiveDistanceRicciFlowRiemannMapping : public CQCVirtualInversiveDistanceRicciFlow<M>
  {
  public:
	  /*! \brief CQCVirtualInversiveDistanceRicciFlowRiemannMapping constructor
	   *
	   *  call base class constructor 
	   */
	  CQCVirtualInversiveDistanceRicciFlowRiemannMapping( M * pMesh );

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
CQCVirtualInversiveDistanceRicciFlowRiemannMapping<M>::CQCVirtualInversiveDistanceRicciFlowRiemannMapping( M * pMesh ): CQCVirtualInversiveDistanceRicciFlow( pMesh)
{
	//_calculate_topo_valence();
}

template<typename M>
void CQCVirtualInversiveDistanceRicciFlowRiemannMapping<M>::_calculate_topo_valence()
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
void CQCVirtualInversiveDistanceRicciFlowRiemannMapping<M>::_set_target_curvature()
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