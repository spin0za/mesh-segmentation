/*!
*      \file ConvexCap.h
*      \brief Algorithm for Convex Cap Embedding
*	   \author David Gu
*      \date Document 03/22/2013
*
*		Convex Cap Isomeetric Embedding
*/


#ifndef _PRISM_H_
#define _PRISM_H_

#include <vector>
#include "PrismMesh.h"


namespace MeshLib
{
/*!
 *	\brief CPrism class
 *
 */
	class CPrism
	{
	public:
		/*!	CConvexCap constructor
		 *	\param pMesh the input mesh
		 */
		CPrism();
		/*!	CHarmonicMapper destructor
		 */
		~CPrism();
		/*!  Compute the harmonic map using direct method
		 */
		void _map();

		CPrMesh & interior_prism() { return m_interior_prism; };
		CPrMesh & boundary_prism() { return m_boundary_prism; };

		/*!
		 * Embed the interior prism
		 */
		void _embed_interior_prism( std::vector<double> & ls, std::vector<double> & hs, CPrMesh & _prism );

		/*!
		 *	Compute face normal of prism
		 */
		void __compute_prism_face_normal( CPrMesh & );
		/*!
		 *	Compute edge dihedral angle of prism
		 */
		void __compute_prism_edge_dihedral_angle( CPrMesh & );
		/*!
		 *	Compute face corner angle
		 */
		void __compute_prism_face_corner_angle( CPrMesh & );


	protected:
		/*!	The input surface mesh
		 */
		CPrMesh m_interior_prism;
		CPrMesh m_boundary_prism;

	protected:
		/*!
		 *	create an interior prism
		 */
		void _create_interior_prism();
		/*!
		 *	create a boundary prism
		 */
		void _create_boundary_prism();
		/*!
		 * Embed the boundary prism
		 */
		//void _embed_boundary_prism( CConvexCapFace * pFace, CPrMesh & _prism );

	};
}

#endif _PRISM_

