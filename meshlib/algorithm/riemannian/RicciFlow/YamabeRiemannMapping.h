/*! \file YamabeRiemannMapping.h
 * \brief Euclidean Yamabe flow for Riemann Mapping
 *  \author David Gu
 *  \date   documented on 11/15/2010
 *
 *	Algorithm for Euclidean Yamabe Ricci flow for Riemann Mapping
 */

#ifndef _YAMABE_FLOW_RIEMANN_MAPPING_H_
#define _YAMABE_FLOW_RIEMANN_MAPPING_H_

#include <map>
#include <vector>

#include "YamabeFlow.h"

namespace MeshLib
{

namespace RicciFlow
{

	/*! \brief CYamabeFlowRiemannMapping class
	 *  
	 *  Algorithm for Euclidean Yamabe flow for Riemann Mapping
	 */
	 
 template<typename M>
 class CYamabeFlowRiemannMapping : public CYamabeFlow<M>
  {
  public:
	  /*! \brief CYamabeFlowRiemannMapping constructor
	   *
	   *  call base class constructor 
	   */
	  CYamabeFlowRiemannMapping( M * pMesh );

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
CYamabeFlowRiemannMapping<M>::CYamabeFlowRiemannMapping( M * pMesh ): CYamabeFlow( pMesh)
{
	//_calculate_topo_valence();
}

template<typename M>
void CYamabeFlowRiemannMapping<M>::_calculate_topo_valence()
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
void CYamabeFlowRiemannMapping<M>::_set_target_curvature()
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