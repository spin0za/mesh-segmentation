/*!
*      \file HarmonicExactForm.h
*      \brief Algorithm for computing exact harmonic forms
*	   \author David Gu
*      \date Document 10/11/2010
*
*		Compute harmonic exact 1-forms for genus zero surfaces with multiple boundary components.
*
*/

/********************************************************************************************************************************
*
*      Harmonic Exact Form Class
*
*       Copyright (c) Stony Brook University
*
*    Purpose:
*
*       Harmonic Exact Form
* 
*       David Gu June 27, 2008,  gu@cs.stonybrook.edu
*
*
*      Input:
*         
*           Multiple Connected Annulus
*
*      Output:
*
*           Exact Harmonic 1-forms, which is -1 on one inner boundary, and 0 on all the other boundaries. The vertex uv is also set
*
*********************************************************************************************************************************/

/*---------------------------------------------------------------------------------------------------------------------------------
#include "ExactForm/HarmonicExactForm.h"
int main( int argc, char * argv[] )
{

	if( strcmp( argv[1], "-exact_form" ) == 0 )
	{
		CHarmonicMesh hm;
		he.read_m( argv[2] );

		CHarmonicExactForm exact_form( & hm );
		exact_form.calculate_harmonic_exact_form( argv[3] );
		return 0;
	}
}

}
----------------------------------------------------------------------------------------------------------------------------------*/

#ifndef _SLIT_MAP_HARMONIC_EXACT_FORM_H_
#define _SLIT_MAP_HARMONIC_EXACT_FORM_H_

#include  <math.h>
#include <queue>
#include "BaseHarmonicExactForm.h"

namespace MeshLib
{

namespace Holomorphy
{
/*! \brief CHarmonicExactForm class
*
*	Compute exact harmonic 1-forms for genus zero surfaces with multiple boundary components
*/
	template<typename M>
   class CHarmonicExactForm :public CBaseHarmonicExactForm<M>
  {
  public:
	  /*!	CHarmonicExactForm constructor
	  *     \param pMesh input mesh
	  */
	  CHarmonicExactForm( M * pMesh ):CBaseHarmonicExactForm( pMesh ){};
	/*!	CHarmonicExactForm destructor
	*/
	  ~CHarmonicExactForm(){};

  protected:
	/*!
	 *	Compute the angle structure
	 */
	void _angle_structure()
	{
		COperator<M> pC( m_pMesh );
		pC._embedding_2_metric();
		pC._metric_2_angle();
		pC._angle_2_Laplace();	//compute edge weight
	}
	
  };

} //namespace Holomorphy
} //namespace MeshLib
#endif _SLIT_MAP_HARMONIC_EXACT_FORM_H_