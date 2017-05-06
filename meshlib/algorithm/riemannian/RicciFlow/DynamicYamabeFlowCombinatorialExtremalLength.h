/*! \file   DynamicYamabeFlowCombinatorialExtremalLength.h
 *  \brief  Dynamic Euclidean Yamabe flow for Extremal Length using combinatorial structure
 *  \author David Gu
 *  \date   documented on 10/27/2013
 *
 *	Algorithm for Dynamic Yamabe Flow for Extremal Length
 */

#ifndef _DYNAMIC_YAMABE_FLOW_COMBINATORIAL_EXTREMAL_LENGTH_H_
#define _DYNAMIC_YAMABE_FLOW_COMBINATORIAL_EXTREMAL_LENGTH_H_


#include "DynamicYamabeFlowExtremalLength.h"


namespace MeshLib
{
namespace RicciFlow
{
	/*! \brief CDynamicYamabeFlowCombinatorialExtremalLength class
	 *  
	 *  Algorithm for Dynamic Euclidean Yamabe flow for Extremal Length using Combinatorial Structure
	 */
	 
 template<typename M>
 class CDynamicYamabeFlowCombinatorialExtremalLength : public CDynamicYamabeFlowExtremalLength<M>
  {
  public:
	  /*! \brief CDynamicYamabeFlowCombinatorialExtremalLength constructor
	   *
	   *  call base class constructor 
	   */
	  CDynamicYamabeFlowCombinatorialExtremalLength( M * pMesh ): CDynamicYamabeFlowExtremalLength<M>( pMesh )
	  {
		  _calculate_initial_edge_length();

		 for( M::MeshEdgeIterator eiter( m_pMesh ); !eiter.end(); eiter ++ )
		 {
			 M::CEdge * pE = *eiter;
			 pE->length() = 1.0;
		 }
		  
	  };
	

  protected:

	 void _calculate_initial_edge_length()
	 {
		 for( M::MeshEdgeIterator eiter( m_pMesh ); !eiter.end(); eiter ++ )
		 {
			 M::CEdge * pE = *eiter;
			 pE->initial_length() = 1.0;
		 }
	 };
  
 };

} //namespace RicciFlow
} //namespace MeshLib

#endif  _DYNAMIC_YAMABE_FLOW_COMBINATORIAL_EXTREMAL_LENGTH_H_