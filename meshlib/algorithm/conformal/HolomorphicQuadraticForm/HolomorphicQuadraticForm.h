/*!
*      \file HolomorphicQuadraticForm.h
*      \brief Algorithm for computing holomorphic quadratic forms
*	   \author David Gu
*      \date Document 06/16/2013
*
*		Algorithm for computing holomorphic 1-forms, Hodge Star operator
*
*/

#ifndef _HOLOMORPHIC_QUADRATIC_FORM_H_
#define _HOLOMORPHIC_QUADRATIC_FORM_H_

#include  <math.h>
#include <queue>
#include <list>
#include <vector>
#include "HolomorphicQuadraticFormMesh.h"
#include "Structure/Structure.h"
#include <Eigen/Dense>
#include <Eigen/Sparse>


namespace MeshLib
{

namespace Holomorphy
{

/*! \brief CHolomorphicQuadraticForm class
 *
 *	Compute holomorphic quadratic forms on a mesh
 */
template<typename M>
class CHolomorphicQuadraticForm
{
public:
	/*! CHolomorphicQuadraticForm constructor
	*	\param meshes the list of meshes with stores the holomorphic quadratic forms, which form the basis 
	*   of the group of all holomorphic quadratic forms
	*/
	//CHolomorphicQuadraticForm( std::list<CHoloQuadFormMesh*> & meshes );
	CHolomorphicQuadraticForm(  M * pMesh );
	/*! CHolomorphicForm destructor
	*/
	~CHolomorphicQuadraticForm();
	/*!	The list of meshes storing holomorphic quadratic forms
	*/
	std::vector< M *> & meshes() { return m_meshes; };

	/*! linear combination
	 */
	void linear_combine( std::list<M *> & meshes, std::vector<double> & coefficients  );

	/*!
	 *	normalize the holomorphic quadratic form
	 */
	void normalize();
	/*!
	 *	convert holomorphic quadratic form to a flat metric with cone singularities
	 */
	void form_to_metric();
	/*!
	 *	pointer to the mesh
	 */
	M * mesh() { return m_pMesh; };

protected:
	/*!	The list of meshes storing harmonic form bases
	*/
	std::vector<M*> m_meshes;
	/*!	Pointer to the mesh
	*/
	M* m_pMesh;

	/*!
	 *	the norm of the holomorphic quadratic form
	 */
	double norm();
};


//CHolomorphicQuadraticForm constructor
//\param meshes are the basis of holomorphic quadratic forms
template<typename M>
CHolomorphicQuadraticForm<M>::CHolomorphicQuadraticForm( M *  pMesh )
{
	m_pMesh = pMesh;

	for( M::MeshEdgeIterator eiter( m_pMesh ); !eiter.end(); eiter ++ )
	{
		M::CEdge * pE = *eiter;
		
		M::CVertex * pV0 = m_pMesh->edgeVertex1( pE );
		M::CVertex * pV1 = m_pMesh->edgeVertex2( pE );

		std::complex<double> dz = pV0->z() - pV1->z();

		pE->dz() = dz * dz;

		pE->length() = std::abs( dz );
	}

};

//CHolomorphicQuadraticForm destructor
template<typename M>
CHolomorphicQuadraticForm<M>::~CHolomorphicQuadraticForm()
{
};




//the norm of the CHolomorphicQuadraticForm
template<typename M>
double CHolomorphicQuadraticForm<M>::norm()
{
	double A = 0;

	for( M::MeshFaceIterator fiter( m_pMesh ); !fiter.end(); fiter ++ )
	{
		M::CFace * pF = *fiter;
		std::vector<double> d;
		for( M::FaceEdgeIterator feiter( pF ); !feiter.end(); feiter ++ )
		{
			M::CEdge * pE = *feiter;
			d.push_back( pE->length() );
		}

		double s = ( d[0] + d[1] + d[2] )/2.0;

		pF->area() = sqrt( s * ( s - d[0] ) * ( s - d[1] ) * ( s - d[2] ) );
		
		A += pF->area();
	}
	std::cout << "Norm is " << A << std::endl;
	return A;
};

//normalize the holomorphic quadratic form
template<typename M>
void CHolomorphicQuadraticForm<M>::normalize()
{
	double A = norm();
	
	for( M::MeshEdgeIterator eiter( m_pMesh ); !eiter.end(); eiter ++ )
	{
		CHoloQuadFormEdge * pE = *eiter;
		pE->dz() /= A;

		pE->length() = std::sqrt( std::abs( pE->dz() ));
	}

	//verification
	norm();

}

//convert holomorphic quadratic form to a flat metric
template<typename M>
void CHolomorphicQuadraticForm<M>::form_to_metric()
{
	for( M::MeshEdgeIterator eiter( m_pMesh ); !eiter.end(); eiter ++ )
	{
		M::CEdge * pE = *eiter;
		pE->length() = std::sqrt( std::abs( pE->dz() ));
	}
}

/*! linear combination
*/
template<typename M>
void CHolomorphicQuadraticForm<M>::linear_combine( std::list<M*> & meshes, std::vector<double> & coefficients )
{
	std::vector<M*> pms;

	for( std::list<M*>::iterator miter = meshes.begin(); miter != meshes.end(); miter ++ )
	{
		M * pM = *miter;
		pms.push_back( pM );
	}

	for( M::MeshEdgeIterator eiter( m_pMesh ); !eiter.end(); eiter ++ )
	{
		M::CEdge * pE = *eiter;

		M::CVertex * pV0 = m_pMesh->edgeVertex1( pE );
		M::CVertex * pV1 = m_pMesh->edgeVertex2( pE );
		
		std::complex<double> dz = 0;

		for( size_t i = 0; i < meshes.size(); i ++ )
		{
			M * pM = pms[i];
			M::CVertex * pW0 = pM->idVertex( pV0->id());
			M::CVertex * pW1 = pM->idVertex( pV1->id());
			M::CEdge   * pWe = pM->vertexEdge( pW0, pW1 );
			dz += pWe->dz() * coefficients[i];
		}

		pE->dz() = dz;
	}

};

} //namespace Holomorphy
} //namespace MeshLib
#endif