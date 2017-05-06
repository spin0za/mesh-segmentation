/*! \file LinearCombination.h
*   \brief Algorithm for lineary combining 1-forms on a mesh
*	\author David Gu
*   \date Documented on 10/26/2014
*/

/********************************************************************************************************************************
*
*      LinearCombinationOneForm  Class
*
*       Copyright (c) Stony Brook University
*
*    Purpose:
*
*       LinearCombination Class
* 
*       David Gu October 26, 2014,  gu@cs.stonybrook.edu
*
*
*      Input:
*         
*           Closed/harmonic 1-form meshes, linearCombinations Coefficients
*
*      Output:
*
*           Closed/harmonic 1-form mesh
*
*********************************************************************************************************************************/


#ifndef _LINEAR_COMBINATION_ONE_FORM_H_
#define _LINEAR_COMBINATION_ONE_FORM_H_

#include <queue>
#include "IntegrationMesh.h"

namespace MeshLib
{
namespace Topology
{
/*! \brief CLinearCombination class
*
*	Compute the linear combination of holomorphic 1-forms
*/ 
template<typename M>
class CLinearCombinationOneForm
{
public:
	/*! CLinearCombination constructor
	* \param pForm input holomorphic 1-form 
	* \param linear combination coefficient
	*/
	CLinearCombinationOneForm( M * pForm, double lambda );
	/*! CLinearCombination destructor */
	~CLinearCombinationOneForm();
	/*! linear combine */
	void _LinearCombine( M * pForm, double lambda);

private:
	/*! Holomorphic 1-form */
	M * m_pForm;
	/*! Linear Combination Coefficient */
	double m_coefficient;
};

//CLinearCombination constructor
//\param pForm the input holomorphic 1-form
//\param lambda the linear combination coefficient
template<typename M>
CLinearCombinationOneForm<M>::CLinearCombinationOneForm( M * pForm, double lambda )
{
	m_pForm = pForm;
	m_coefficient = lambda;


	for( M::MeshEdgeIterator eiter( m_pForm ); !eiter.end(); eiter ++ )
	{
		M::CEdge * e = *eiter;
		e->du() = e->du() * lambda;
	}	
};

//CLinearCombination destructor
template<typename M>
CLinearCombinationOneForm<M>::~CLinearCombinationOneForm()
{
};


//Linear combination
template<typename M>
void CLinearCombinationOneForm<M>::_LinearCombine(M * pForm, double lambda)
{

	for( M::MeshEdgeIterator eiter( pForm ); !eiter.end(); eiter ++ )
	{
		M::CEdge * e = *eiter;
		M::CVertex * v1 = pForm->edgeVertex1( e );
		M::CVertex * v2 = pForm->edgeVertex2( e );
		
		int id1 = v1->id();
		int id2 = v2->id();
		
		M::CVertex * w1 = m_pForm->idVertex( id1 );
		M::CVertex * w2 = m_pForm->idVertex( id2 );
		M::CEdge * we = m_pForm->vertexEdge( w1, w2 );
		
		we->du() = we->du() + (e->du() * lambda);
	}	
}



} //namespace Topology
} //namespace MeshLib

#endif