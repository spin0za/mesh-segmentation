/*!
*      \file QCHolomorphicForm.h
*      \brief Algorithm for computing holomorphic 1-forms, Hodge star operator, using quasi-conformal structure
*	   \author David Gu
*      \date Document 10/12/2010
*
*		Algorithm for computing holomorphic 1-forms, Hodge Star operator
*
*/

#ifndef _QC_HOLOMORPHIC_FORM_H_
#define _QC_HOLOMORPHIC_FORM_H_

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
class CQCHolomorphicForm : public  CBaseHolomorphicForm<M>
{
public:
	/*! CQCHolomorphicForm constructor
	*	\param meshes the list of meshes with stores the harmonic 1-forms, which form the basis 
	*   of the first cohomology group of the mesh
	*/
	CQCHolomorphicForm( std::list<M*> & meshes );
	/*! CHolomorphicForm destructor
	*/
	~CQCHolomorphicForm();
protected:
	void _angle_structure();
};

//CQCHolomorphicForm constructor
//\param meshes are the basis of harmonic 1-form group
template<typename M>
CQCHolomorphicForm<M>::CQCHolomorphicForm( std::list<M*> & meshes ) :CBaseHolomorphicForm<M>( meshes )
{
};

//CHolomorphicForm destructor
template<typename M>
CQCHolomorphicForm<M>::~CQCHolomorphicForm()
{
};

//Compute the angle structure
template<typename M>
void CQCHolomorphicForm<M>::_angle_structure()
{
	for( std::vector<M*>::iterator miter = m_meshes.begin(); miter != m_meshes.end(); miter ++ )
	{
		M * pM = *miter;

		COperator<M> pC( pM );
		pC._parameter_mu_2_metric();
		pC._metric_2_angle();
	}

}


} //namespace Holomorphy
} //namespace MeshLib

#endif _QC_HOLOMORPHIC_FORM_H_