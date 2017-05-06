/*!
*      \file HarmonicForm.h
*      \brief Algorithm for base harmonic forms (for closed surfaces)
*	   \author David Gu
*      \date Document 12/26/2010
*
*		Computing harmonic forms (for closed surfaces)
*
*/

/********************************************************************************************************************************
*
*      Harmonic Form Class
*
*       Copyright (c) Stony Brook University
*
*    Purpose:
*
*       Harmonic Form
* 
*       David Gu December 26, 2010,  gu@cs.stonybrook.edu
*
*
*      Input:
*         
*           Original mesh, the mesh cut through a homology basis
*
*      Output:
*
*           Closed non-exact Harmonic 1-form. and the mesh with UV coordinates.
*
*********************************************************************************************************************************/


#ifndef _HEAT_FLOW_H_
#define _HEAT_FLOW_H_

#include "BaseHeatFlow.h"

namespace MeshLib
{

namespace Holomorphy
{
	/*! \brief CBaseHarmonicForm class
	*  
	*  Algorithm for computing harmonic forms
	*/
	template<typename M>
	class CHeatFlow : public CBaseHeatFlow<M>
  {
  public:
	  /*! CHeatFlow constructor
	  *
	  *	\param pMesh input closed mesh
	  * \param pWMesh sliced open mesh
	  */

	  CHeatFlow( M * pMesh, M * pWMesh ) : CBaseHeatFlow<M>( pMesh, pWMesh ){};
	/*! CBaseHarmonicForm destructor
	*/

	  ~CHeatFlow(){};

  protected:
	  void _angle_structure()
	  {
		COperator<CHCFMesh> pC( m_pMesh );
		pC._embedding_2_metric();
		pC._metric_2_angle();
		pC._angle_2_Laplace();	//compute edge weight
	  };

  };
}	//namespace Holomorphy
}   //namespace MeshLib
#endif