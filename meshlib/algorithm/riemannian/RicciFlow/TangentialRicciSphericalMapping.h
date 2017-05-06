/*! \file TangentialRicciSphericalMapping.h
 * \brief Euclidean Tangential Ricci flow for Spherical Mapping
 *  \author David Gu
 *  \date   documented on 07/31/2011
 *
 *	Algorithm for Tangential Ricci Flow for Spherical Mapping
 */

#ifndef _TANGENTIAL_RICCI_FLOW_SPHERICAL_MAPPING_H_
#define _TANGENTIAL_RICCI_FLOW_SPHERICAL_MAPPING_H_

#include <map>
#include <vector>

#include "TangentialRicciFlow.h"

namespace MeshLib
{

namespace RicciFlow
{

	/*! \brief CTangentialRicciFlowSphericalMapping class
	 *  
	 *  Algorithm for Tangential Ricci flow for Spherical Mapping
	 */
	 
 template<typename M>
 class CTangentialRicciFlowSphericalMapping : public CTangentialRicciFlow<M>
  {
  public:
	  /*! \brief CTangentialRicciFlowSphericalMapping constructor
	   *
	   *  call base class constructor 
	   */
	  CTangentialRicciFlowSphericalMapping( M * pMesh );

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
CTangentialRicciFlowSphericalMapping<M>::CTangentialRicciFlowSphericalMapping( M * pMesh ): CTangentialRicciFlow( pMesh)
{
	_calculate_topo_valence();
}

template<typename M>
void CTangentialRicciFlowSphericalMapping<M>::_calculate_topo_valence()
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
void CTangentialRicciFlowSphericalMapping<M>::_set_target_curvature()
{
	int count = 0;
  for( M::MeshVertexIterator viter( m_pMesh ); !viter.end(); viter ++ )
  {
	  M::CVertex * v = *viter;

	if( !v->boundary() )
	{
		v->target_k() = 0;
		count ++;
	}
	else
	{
		v->target_k() = 2.0 * PI/3.0;
	}
  }
	printf("There are %d corners.\n", count );
};

} //namespace RicciFlow
} //namespace MeshLib

#endif  