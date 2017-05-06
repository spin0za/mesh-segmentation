/*!
*      \file BaseHarmonicExactForm.h
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

#ifndef _SLIT_MAP_BASE_HARMONIC_EXACT_FORM_H_
#define _SLIT_MAP_BASE_HARMONIC_EXACT_FORM_H_

#include  <math.h>
#include <queue>
#include "HarmonicMesh.h"
#include "Operator/Operator.h"
#include "Laplace/Laplace.h"

namespace MeshLib
{

namespace Holomorphy
{

/*! \brief CBaseHarmonicExactForm class
*
*	Compute exact harmonic 1-forms for genus zero surfaces with multiple boundary components
*/
  template<typename M>
  class CBaseHarmonicExactForm
  {
  public:
	  /*!	CBaseHarmonicExactForm constructor
	  *     \param pMesh input mesh
	  */
    CBaseHarmonicExactForm( M * pMesh );
	/*!	CBaseHarmonicExactForm destructor
	*/
    ~CBaseHarmonicExactForm();
	/*!	Calculate harmonic exact forms
	*   \param prefix the prefix of output mesh file
	*/
    void calculate_harmonic_exact_form( const char * prefix);

  protected:
	  /*! the input mesh */
    M * m_pMesh;
		/*! the boundary of the input mesh */
	typename M::CBoundary m_boundary;
	
	/*! Compute harmonic exact form, the harmonic function equals to 1 on pL, and zero on other boundary components.
	 * \param pL the boundary loop, on which the function equals to one.
	 */
	void _harmonic_exact_form( typename M::CLoop * pL );
	/*!
	 *	Compute the angle structure
	 */
	virtual void _angle_structure()=0;

  };


//CBaseHarmonicExactForm constructor
// pMesh is the input mesh
template<typename M>
CBaseHarmonicExactForm<M>::CBaseHarmonicExactForm( M * pMesh ):m_pMesh( pMesh ), m_boundary( m_pMesh )
{
	for( M::MeshVertexIterator viter( m_pMesh ); !viter.end(); ++ viter )
	{
		M::CVertex * pV = *viter;
		pV->fixed() = pV->boundary();
	}
}


//CHarmonicExactForm destructor
template<typename M>
CBaseHarmonicExactForm<M>::~CBaseHarmonicExactForm()
{
}


//Compute the harmonic exact form, such that the harmonic function
//equals to one on the boundary loop pL
template<typename M>
void CBaseHarmonicExactForm<M>::_harmonic_exact_form( typename M::CLoop * pL )
{
	for( M::MeshVertexIterator viter( m_pMesh ); !viter.end(); ++ viter )
	{
		M::CVertex * pV = *viter;
		pV->u() = 0.0;
	}

	for( std::list<M::CHalfEdge*>::iterator hiter = pL->halfedges().begin(); hiter != pL->halfedges().end(); hiter ++ )
	{
		M::CHalfEdge * he = *hiter;
		M::CVertex   * pV = m_pMesh->halfedgeVertex( he );
		pV->u() = -1.0;
	}

	CLaplace<CHarmonicMesh> L( m_pMesh );

	L.solve();
	
	for( M::MeshVertexIterator viter( m_pMesh ); !viter.end(); ++ viter )
	{
		M::CVertex * pV = *viter;
		pV->uv() = CPoint2( pV->u(), 0.53 );
	}

	for( M::MeshEdgeIterator eiter( m_pMesh ); !eiter.end(); eiter ++ )
	{
		M::CEdge * e = *eiter;
		M::CVertex * v1 = m_pMesh->edgeVertex1( e );
		M::CVertex * v2 = m_pMesh->edgeVertex2( e );
		e->du() = v2->u() - v1->u();
	}
}

//Compute harmonic exact forms 
//output the results to the meshes "prefix_k.du.m"
template<typename M>
void CBaseHarmonicExactForm<M>::calculate_harmonic_exact_form(const char * prefix )
{

	_angle_structure();		

	std::vector<M::CLoop*>& loops = m_boundary.loops();

	for( size_t k = 1; k < loops.size(); k ++ )
	{
		printf("%f \n", loops[k]->length());
	}

	for( size_t k = 1; k < loops.size(); k ++ )
	{		
	  _harmonic_exact_form( loops[k] );

	  std::stringstream iss;
	  iss << prefix << "_" << k-1 << ".du.m ";
	  m_pMesh->write_m( iss.str().c_str() );
	}
}



} //namespace Holomorphy
} //namespace MeshLib
#endif _SLIT_MAP_BASE_HARMONIC_EXACT_FORM_H_