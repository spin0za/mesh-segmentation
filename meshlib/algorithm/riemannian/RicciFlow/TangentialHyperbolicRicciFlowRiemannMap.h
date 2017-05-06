/*! \file TangentialHyperbolicRicciFlowRiemannMap.h
 *  \brief Tangential Hyperbolic Ricci flow for Riemann Mapping algorithm
 *  \author David Gu
 *  \date   documented on 08/03/2011
 *
 *	Algorithm for general Ricci Flow
 */

#ifndef _TANGENTIAL_HYPERBOLIC_RICCI_FLOW_RIEMANN_MAP_H_
#define _TANGENTIAL_HYPERBOLIC_RICCI_FLOW_RIEMANN_MAPH_

#include <map>
#include <vector>

#include "TangentialMixedHyperbolicRicciFlow.h"

namespace MeshLib
{

namespace RicciFlow
{
/*! \brief Class CTangentialHyperbolicRicciFlow
*
*	Algorithm for computing Hyperbolic Ricci flow using tangential Circle packing metric
*/
template<typename M>
class CTangentialHyperbolicRicciFlowRiemannMap : public CTangentialMixedHyperbolicRicciFlow<M>
  {
  public:
    /*! \brief CTangentialRicciFlowRiemannMap constructor
	 *  \param pMesh the input mesh
	 */
	  CTangentialHyperbolicRicciFlowRiemannMap( M * pMesh, double boundary_u ): CTangentialMixedHyperbolicRicciFlow( pMesh ){ m_boundary_u = boundary_u; };
    /*! \brief CTangentialHyperbolicRicciFlowRiemannMap destructor
	 */
	  ~CTangentialHyperbolicRicciFlowRiemannMap() {};

  protected:

	/*!
	 *	Set the target curvature on each interior vertex
	 */
    void    _set_target_curvature();
	/*!
	 *	Set the target conformal factor on each boundary vertex
	 */
    void    _set_target_conformal_factor();
	/*!
	 *	the boundary u value
	 */
	double m_boundary_u;
  };




//set target curvature

template<typename M>
void CTangentialHyperbolicRicciFlowRiemannMap<M>::_set_target_curvature()
{
  for( M::MeshVertexIterator viter( m_pMesh ); !viter.end(); viter ++ )
  {
	  M::CVertex * v = *viter;
	if( m_pMesh->isBoundary( v ) ) continue;
    v->target_k() = 0;
  }
};


//set target conformal factor

template<typename M>
void CTangentialHyperbolicRicciFlowRiemannMap<M>::_set_target_conformal_factor()
{
  for( M::MeshVertexIterator viter( m_pMesh ); !viter.end(); viter ++ )
  {
	  M::CVertex * v = *viter;
	if( ! m_pMesh->isBoundary( v ) ) continue;
    v->u() = m_boundary_u;
  }
};

} //namespace RicciFlow
} //namespace MeshLib

#endif  _TANGENTIAL_HYPERBOLIC_RICCI_FLOW_RIEMANN_MAP_H_