/*! \file LinearCombination.h
*   \brief Algorithm for lineary combining holomorphic 1-forms on a mesh
*	\author David Gu
*   \date Documented on 1/29/2011
*/

/********************************************************************************************************************************
*
*      LinearCombination  Class
*
*       Copyright (c) Stony Brook University
*
*    Purpose:
*
*       LinearCombination Class
* 
*       David Gu January 9, 2011,  gu@cs.stonybrook.edu
*
*
*      Input:
*         
*           Holomorphic 1-form meshes, linearCombinations Coefficients
*
*      Output:
*
*           Holomorphic 1-form mesh
*
*********************************************************************************************************************************/


#ifndef _LINEAR_COMBINATION_H_
#define _LINEAR_COMBINATION_H_

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
class CLinearCombination
{
public:
	/*! CLinearCombination constructor
	* \param pForm input holomorphic 1-form 
	* \param linear combination coefficient
	*/
	CLinearCombination( M * pForm, double lambda );
	/*! CLinearCombination destructor */
	~CLinearCombination();
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
CLinearCombination<M>::CLinearCombination( M * pForm, double lambda )
{
	m_pForm = pForm;
	m_coefficient = lambda;


	for( M::MeshEdgeIterator eiter( m_pForm ); !eiter.end(); eiter ++ )
	{
		M::CEdge * e = *eiter;
		e->duv() = e->duv() * lambda;
	}	
};

//CLinearCombination destructor
template<typename M>
CLinearCombination<M>::~CLinearCombination()
{
};


//Linear combination
template<typename M>
void CLinearCombination<M>::_LinearCombine(M * pForm, double lambda)
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
		
		we->duv() = we->duv() + (e->duv() * lambda);
	}	
}



} //namespace Topology
} //namespace MeshLib

#endif