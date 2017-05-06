/*!
*      \file HolomorphicForm.h
*      \brief Algorithm for computing holomorphic 1-forms, Hodge star operator
*	   \author David Gu
*      \date Document 10/12/2010
*
*		Algorithm for computing holomorphic 1-forms, Hodge Star operator
*
*/
/********************************************************************************************************************************
*
*      Holomorphic 1-form Class
*
*       Copyright (c) Stony Brook University
*
*    Purpose:
*
*       Compute the holomorphic 1-forms
* 
*       David Gu June 27, 2008,  gu@cs.stonybrook.edu
*
*
*      Input:
*         
*           Original mesh, the mesh cut through a shortest path connecting one inner boundary and the exterior boundary
*
*      Output:
*
*           Closed non-exact Harmonic 1-form. and the mesh with UV coordinates.
*
*********************************************************************************************************************************/

/*---------------------------------------------------------------------------------------------------------------------------------
#include "HolomorphicForm/HolomorphicForm.h"

using namespace MeshLib;

int main( int argc, char * argv[] )
{
	if( strcmp( argv[1], "-holomorphic_form" ) == 0 )
	{
	std::list<CHoloFormMesh*> meshes;


	for( int i = 2; i < argc; i ++ )
	{
		CHoloFormMesh * pMesh = new CHoloFormMesh;
		assert( pMesh );
		pMesh->read_m( argv[i] );
		meshes.push_back( pMesh );
	}

	CHolomorphicForm form( meshes );
	form.conjugate();

	int id = 0;
	for( std::list<CHoloFormMesh*>::iterator miter = meshes.begin(); miter != meshes.end(); miter++)
	{
		CHoloFormMesh * pM = *miter;
		std::stringstream iss;
		iss << "holomorphic_form_" << id++ << ".m";
		pM->write_m( iss.str().c_str() );
	}

	for( std::list<CHoloFormMesh*>::iterator miter = meshes.begin(); miter != meshes.end(); miter++)
	{
		CHoloFormMesh * pM = *miter;
		delete pM;
	}
	
}
----------------------------------------------------------------------------------------------------------------------------------*/
#ifndef _HOLOMORPHIC_FORM_H_
#define _HOLOMORPHIC_FORM_H_

#include  <math.h>
#include <queue>
#include <list>
#include <vector>
#include "HolomorphicFormMesh.h"
#include "BaseHolomorphicForm.h"
#include "Structure/Structure.h"

namespace MeshLib
{

namespace Holomorphy
{
/*! \brief CHolomorphicForm class
 *
 *	Compute holomorphic forms on a mesh
 */
template<typename M>
class CHolomorphicForm : public  CBaseHolomorphicForm<M>
{
public:
	/*! CHolomorphicForm constructor
	*	\param meshes the list of meshes with stores the harmonic 1-forms, which form the basis 
	*   of the first cohomology group of the mesh
	*/
	CHolomorphicForm( std::list<M*> & meshes );
	/*! CHolomorphicForm destructor
	*/
	~CHolomorphicForm();
protected:
	void _angle_structure();
};

//CHolomorphicForm constructor
//\param meshes are the basis of harmonic 1-form group
template<typename M>
CHolomorphicForm<M>::CHolomorphicForm( std::list<M*> & meshes ) :CBaseHolomorphicForm<M>( meshes )
{
};

//CHolomorphicForm destructor
template<typename M>
CHolomorphicForm<M>::~CHolomorphicForm()
{
};

//Compute the angle structure
template<typename M>
void CHolomorphicForm<M>::_angle_structure()
{
	for( std::vector<M*>::iterator miter = m_meshes.begin(); miter != m_meshes.end(); miter ++ )
	{
		M * pM = *miter;

		COperator<M> pC( pM );
		pC._embedding_2_metric();
		pC._metric_2_angle();
	}

}


} //namespace Holomorphy
} //namespace MeshLib
#endif