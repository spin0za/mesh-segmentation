/*! \file QCYamabeRiemannMapping.h
 * \brief Euclidean Quasi-Conformal Yamabe flow for Riemann Mapping
 *  \author David Gu
 *  \date   documented on 03/08/2011
 *
 *	Algorithm for Euclidean Yamabe Ricci flow for Riemann Mapping
 */

#ifndef _QC_YAMABE_FLOW_RIEMANN_MAPPING_H_
#define _QC_YAMABE_FLOW_RIEMANN_MAPPING_H_

#include <map>
#include <vector>

#include "QCYamabeFlow.h"

namespace MeshLib
{

namespace RicciFlow
{
	/*! \brief CQCYamabeFlowRiemannMapping class
	 *  
	 *  Algorithm for Quasi-Conformal Euclidean Yamabe flow for Riemann Mapping
	 */
	 
 template<typename M>
 class CQCYamabeFlowRiemannMapping : public CQCYamabeFlow<M>
  {
  public:
	  /*! \brief CQCYamabeFlowRiemannMapping constructor
	   *
	   *  call base class constructor 
	   */
	  CQCYamabeFlowRiemannMapping( M * pMesh );

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
CQCYamabeFlowRiemannMapping<M>::CQCYamabeFlowRiemannMapping( M * pMesh ): CQCYamabeFlow( pMesh)
{
}

template<typename M>
void CQCYamabeFlowRiemannMapping<M>::_calculate_topo_valence()
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
void CQCYamabeFlowRiemannMapping<M>::_set_target_curvature()
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