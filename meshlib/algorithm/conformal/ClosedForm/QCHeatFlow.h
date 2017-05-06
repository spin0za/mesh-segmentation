/*!
*      \file QCHeatFlow.h
*      \brief Algorithm for base harmonic closed non-exact forms under quasi-conformal structure
*	   \author David Gu
*      \date Document 06/23/2011
*
*		Computing harmonic forms (for closed surfaces)
*
*/

#ifndef _QC_HEAT_FLOW_H_
#define _QC_HEAT_FLOW_H_

#include "BaseHeatFlow.h"

namespace MeshLib
{

namespace Holomorphy
{
	/*! \brief CQCHeatFlow class
	*  
	*  Algorithm for computing harmonic forms under quasi-conformal structure
	*/
	template<typename M>
	class CQCHeatFlow : public CBaseHeatFlow<M>
  {
  public:
	  /*! CQCHeatFlow constructor
	  *
	  *	\param pMesh input closed mesh
	  * \param pWMesh sliced open mesh
	  */

	  CQCHeatFlow( M * pMesh, M * pWMesh ) : CBaseHeatFlow<M>( pMesh, pWMesh ){};
	/*! CQCHeatFlow
	*/

	  ~CQCHeatFlow(){};

  protected:
	  void _angle_structure()
	  {
		COperator<M> pC( m_pMesh );
		pC._parameter_mu_2_metric();
		pC._metric_2_angle();
		pC._angle_2_Laplace();
	  };

  };

}	//namespace Holomorphy
}	//namespace MeshLib
#endif