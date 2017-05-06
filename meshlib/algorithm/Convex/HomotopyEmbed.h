/*!
*      \file HomotopyEmbed.h
*      \brief Algorithm for Convex Embedding
*	   \author David Gu
*      \date Document 03/22/2013
*
*		Convex Isomeetric Embedding
*/


#ifndef _HOMOTOPY_EMBED_H_
#define _HOMOTOPY_EMBED_H_

#include <vector>
#include "HomotopyMesh.h"


namespace MeshLib
{
/*!
 *	\brief CHomotopyEmbed class
 *
 */
	class CHomotopyEmbed
	{
	public:
		/*!	CHomotopyEmbed constructor
		 *	\param pMesh the input mesh
		 */
		CHomotopyEmbed(CHTMesh* pMesh);
		/*!	CHarmonicMapper destructor
		 */
		~CHomotopyEmbed();
		/*!  Compute the harmonic map using direct method
		 */
		void _embed();

	protected:
		/*!	The input surface mesh
		 */
		CHTMesh * m_pMesh;
		/*!	The boundary of m_pMesh
		 */
		CHTMesh::CBoundary m_boundary;
	};

};

#endif _HOMOTOPY_EMBED_


