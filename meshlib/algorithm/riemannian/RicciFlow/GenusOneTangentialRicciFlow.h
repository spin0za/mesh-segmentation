/*! \file GenusOneTangentialRicciFlow.h
 *  \brief General Euclidean Ricci flow algorithm for genus one surface with boundaries
 *  \author David Gu
 *  \date   documented on 10/13/2012
 *
 *	Algorithm for general Ricci Flow
 */

#ifndef _GENUS_ONE_TANGENTIAL_RICCI_FLOW_H_
#define _GENUS_ONE_TANGENTIAL_RICCI_FLOW_H_

#include "TangentialRicciFlow.h"


namespace MeshLib
{

namespace RicciFlow
{

/*! \brief Class CGenusOneTangentialRicciFlow
*
*	Algorithm for computing Ricci flow for genus one surface 
*/
template<typename M>
class CGenusOneTangentialRicciFlow : public CTangentialRicciFlow<M>
  {
  public:
    /*! \brief CGenusOneTangentialRicciFlow constructor
	 *  \param pMesh the input mesh
	 */
	  CGenusOneTangentialRicciFlow( M * pMesh ): CTangentialRicciFlow<M>( pMesh ) {};
    /*! \brief CGenusOneTangentialRicciFlow destructor
	 */
	  ~CGenusOneTangentialRicciFlow(){};

  protected:

	/*!
	 *	Set the target curvature on each vertex
	 */
    virtual void    _set_target_curvature();

  };



//set target curvature

template<typename M>
void CGenusOneTangentialRicciFlow<M>::_set_target_curvature()
{
  for( M::MeshVertexIterator viter( m_pMesh ); !viter.end(); viter ++ )
  {
	  M::CVertex * v = *viter;
	  v->target_k() = 0;
  }

  std::vector<M::CLoop*>& pLs = m_boundary.loops();
	

  for( std::vector<M::CLoop*>::iterator liter = pLs.begin(); liter != pLs.end(); liter++)
  {
	  M::CLoop * pL = *liter;
	  std::list<M::CHalfEdge*> & pHes = pL->halfedges();

	  double sum = 0;
	  double inv_sum = 0;
	
	  for( std::list<M::CHalfEdge*>::iterator hiter = pHes.begin(); hiter != pHes.end(); hiter ++ )
	  {
		  M::CHalfEdge  * he = *hiter;
		  M::CEdge      * pe = m_pMesh->halfedgeEdge( he );

			sum  += pe->length();
			inv_sum += 1/pe->length();
 	  }

	  for( std::list<M::CHalfEdge*>::iterator hiter = pHes.begin(); hiter != pHes.end(); hiter ++ )
	  {
		  M::CHalfEdge * ce = *hiter;
		  M::CVertex   * pv = m_pMesh->halfedgeTarget( ce );
		  M::CHalfEdge * he = m_pMesh->vertexMostCcwInHalfEdge( pv );
		  M::CHalfEdge * te =  m_pMesh->vertexMostClwOutHalfEdge( pv );

			double L = ( m_pMesh->halfedgeEdge(he)->length() + m_pMesh->halfedgeEdge( te )->length()  )/2.0;
			double tk = -2*PI*L/sum;
			pv->target_k() = tk;
	  }
  }

};

}	//namespace RicciFlow
	
}	//namespace MeshLib

#endif  _GENUS_ONE_TANGENTIAL_RICCI_FLOW_H_