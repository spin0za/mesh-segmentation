/*!
*      \file CombinatorialHarmonicForm.h
*      \brief Algorithm for base harmonic forms (for closed surfaces) assuming all the triangles are equilateral
*	   \author David Gu
*      \date Document 10/25/2014
*
*		Computing harmonic forms (for closed surfaces)
*
*/

/********************************************************************************************************************************
*
*      Combinatorial Harmonic Form Class
*
*       Copyright (c) Stony Brook University
*
*    Purpose:
*
*       Combinatorial Harmonic Form
* 
*       David Gu October 25, 2010,  gu@cs.stonybrook.edu
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


#ifndef _COMBINATORIAL_HEAT_FLOW_H_
#define _COMBINATORIAL_HEAT_FLOW_H_

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
	class CCombinatorialHeatFlow : public CBaseHeatFlow<M>
  {
  public:
	  /*! CCombinatorialHeatFlow constructor
	  *
	  *	\param pMesh input closed mesh
	  * \param pWMesh sliced open mesh
	  */

	  CCombinatorialHeatFlow( M * pMesh, M * pWMesh ) : CBaseHeatFlow<M>( pMesh, pWMesh ){};
	/*! CCombinatorialHarmonicForm destructor
	*/

	  ~CCombinatorialHeatFlow(){};

  protected:
	  void _angle_structure()
	  {
		COperator<CHCFMesh> pC( m_pMesh );
		pC._combinatorial_Laplace();	//compute edge weight
	  };

  };
}	//namespace Holomorphy
}   //namespace MeshLib
#endif