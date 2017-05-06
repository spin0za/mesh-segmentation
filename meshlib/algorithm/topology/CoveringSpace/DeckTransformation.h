/*!
*      \file DeckTransformation.h
*      \brief Compute Deck Transformation Group
*	   \author David Gu
*      \date Documented 03/31/2011
*
*/

#ifndef  _DECK_TRANSFORMATION_H_
#define  _DECK_TRANSFORMATION_H_

#include <map>
#include <vector>

#include "Riemannian/RicciFlow/RicciFlowMesh.h"

namespace MeshLib
{
namespace Topology
{
/*!
*	\brief CDeckTransformation
*	
*	Deck Transformation Group
* 
*/
template<typename M>
class CDeckTransformation
{

public:
	/*!
	 *	CDeckTransformation constructor
	 *  \param pMesh    closed mesh
	 *  \param pDomain  foundamental domain
	 */
	CDeckTransformation( M * pMesh, M * pDomain ):m_pMesh( pMesh), m_pDomain(pDomain), m_boundary( pDomain) {};
	/*!
	 *	CDeckTransformation destructor
	 */
	~CDeckTransformation() {};
	virtual void _calculate_generators() = 0;

protected:
	
	M * m_pMesh;
	M * m_pDomain;
	typename M::CBoundary m_boundary;

};


}
};
#endif  _DECK_TRANSFORMATION_H_