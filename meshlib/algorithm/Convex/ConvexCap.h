/*!
*      \file ConvexCap.h
*      \brief Algorithm for Convex Cap Embedding
*	   \author David Gu
*      \date Document 03/22/2013
*
*		Convex Cap Isomeetric Embedding
*/


#ifndef _CONVEX_CAP_H_
#define _CONVEX_CAP_H_

#include <vector>
#include "ConvexCapMesh.h"


namespace MeshLib
{
/*!
 *	\brief CConvexCap class
 *
 */
	class CConvexCap
	{
	public:
		/*!	CConvexCap constructor
		 *	\param pMesh the input mesh
		 */
		CConvexCap(CCCMesh* pMesh);
		/*!	CHarmonicMapper destructor
		 */
		~CConvexCap();
		/*!  Compute the harmonic map using direct method
		 */
		void _map();

	protected:
		/*!	The input surface mesh
		 */
		CCCMesh * m_pMesh;
		/*!	The boundary of m_pMesh
		 */
		CCCMesh::CBoundary m_boundary;
		/*!
		 *  The interior prism
		 */
		CCCMesh m_interior_prism;
		/*!
		 *	The boundary prism
		 */
		CCCMesh m_boundary_prism;

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
		 * Embed the interior prism
		 */
		void _embed_interior_prism( CConvexCapFace * pFace, CCCMesh & _prism );
		/*!
		 * Embed the boundary prism
		 */
		void _embed_boundary_prism( CConvexCapFace * pFace, CCCMesh & _prism );

	private:
		/*!
		 *	Compute face normal of prism
		 */
		void __compute_prism_face_normal( CCCMesh & );
		/*!
		 *	Compute edge dihedral angle of prism
		 */
		void __compute_prism_edge_dihedral_angle( CCCMesh & );
	};
}

#endif _CONVEX_CAP_

/*-------------------------------------------------------------------------------------------
//Example using CHarmonicMapper

#include "HarmonicMapper/HarmonicMapper.h"

using namespace MeshLib;	


void help(char * exe )
{
	printf("Usage:\n");
	printf("%s -harmonic_map input_mesh output_mesh\n", exe);


}
int main( int argc, char * argv[] )
{
	if( argc < 3 )
	{
		help( argv[0] );
		return 0;
	}
	if( strcmp( argv[1] , "-harmonic_map") == 0 )
	{
		CHMMesh mesh;
		mesh.read_m( argv[2] );

		CHarmonicMapper mapper( & mesh );
		mapper._map();
		mesh.write_m( argv[3] );
		return 0;
	}

	help( argv[0] );
	return 0;
} 
-------------------------------------------------------------------------------------------*/