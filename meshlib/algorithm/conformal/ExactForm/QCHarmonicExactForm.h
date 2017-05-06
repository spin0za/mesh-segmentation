/*!
*      \file QCHarmonicExactForm.h
*      \brief Algorithm for computing exact harmonic forms using Quasi-Conformal structure
*	   \author David Gu
*      \date Document 06/23/2011
*
*		Compute harmonic exact 1-forms for genus zero surfaces with multiple boundary components.
*
*/



#ifndef _SLIT_MAP_QC_HARMONIC_EXACT_FORM_H_
#define _SLIT_MAP_QC_HARMONIC_EXACT_FORM_H_

#include  <math.h>
#include <queue>
#include "BaseHarmonicExactForm.h"

namespace MeshLib
{
namespace Holomorphy
{

/*! \brief CQCHarmonicExactForm class
*
*	Compute exact harmonic 1-forms for genus zero surfaces with multiple boundary components under 
*	quasi-conformal structure
*/
	template<typename M>
	class CQCHarmonicExactForm :public CBaseHarmonicExactForm<M>
  {
  public:
	  /*!	CQCHarmonicExactForm constructor
	  *     \param pMesh input mesh
	  */
	  CQCHarmonicExactForm( M * pMesh ):CBaseHarmonicExactForm( pMesh ){};
	/*!	CQCHarmonicExactForm destructor
	*/
	  ~CQCHarmonicExactForm(){};

  protected:
	/*!
	 *	Compute the angle structure
	 */
	void _angle_structure()
	{
		COperator<M> pC( m_pMesh );
		pC._parameter_mu_2_metric();
		pC._metric_2_angle();
		pC._angle_2_Laplace();
	}
	
  };

} //namespace Holomorphy
} //namespace MeshLib

#endif _SLIT_MAP_QC_HARMONIC_EXACT_FORM_H_