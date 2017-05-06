/*!
*      \file QuadrilateralHarmonicMapper.h
*      \brief Algorithm for harmonic mapping for topological quadrilaterals
*	   \author David Gu
*      \date Document 09/29/2013
*
*		Simple harmonic map, that maps a topological quadrilateral to topological quadrilateral
*		with minimal harmonic energy.
*/


#ifndef _QUADRILATERAL_HARMONIC_MAPPER_H_
#define _QUADRILATERAL_HARMONIC_MAPPER_H_

#include <vector>
#include "TMapMesh.h"
#include "BaseHarmonicMapper.h"

namespace MeshLib
{

namespace Holomorphy
{
/*!
 *	\brief CQuadrilateralHarmonicMapper class
 *
 */

	template<typename M>
	class CQuadrilateralHarmonicMapper : public CBaseHarmonicMapper<M>
	{
	public:
		/*!	CQuadrilateralHarmonicMapper constructor
		 *	\param pMesh the input mesh
		 */
		CQuadrilateralHarmonicMapper(M * pMesh ):CBaseHarmonicMapper(pMesh){};
		/*!	CBaseHarmonicMapper destructor
		 */
		~CQuadrilateralHarmonicMapper(){};

		/*! Set the boundary markers
		 */
		void set_boundary_markers( std::vector<typename M::CVertex*> & markers, std::vector<CPoint2> & positions )
		{
			m_boundary_markers = markers;
			m_boundary_positions= positions;
			_boundary_segmentation( m_boundary_markers );
		};

	protected:
		/*!
		 *	boundary segmentation
		 */
		void _boundary_segmentation( std::vector<typename M::CVertex*>& markers );
		/*!
		 *	boundary segments, each segment is an array of vertices
		 */
		std::vector<std::vector<typename M::CVertex*>> m_boundary_segments;
	
		virtual void _manipulate_mu();
		

		void set_boundary_constraints( int k );
	};


	template<typename M>
	void	CQuadrilateralHarmonicMapper<M>::_manipulate_mu()
	{
		double max = 0;
		//manupilate mu
		double sum = 0;
		double variance = 0;

		for( M::MeshFaceIterator fiter( m_pMesh ); !fiter.end(); fiter ++ )
		{
			M::CFace * f = *fiter;
			std::complex<double> mu = f->mu();
			double r = std::abs(mu);
			variance += r * r;
			max = ( max > r )? max:r;
			if( r > 1 ) r = 0.99;
			sum += r;
		}
		double mean = sum/m_pMesh->numFaces();
		variance -= mean * mean;

		for( M::MeshFaceIterator fiter( m_pMesh ); !fiter.end(); fiter ++ )
		{
			M::CFace * f = *fiter;
			std::complex<double> mu = f->mu();
			double r = std::abs(mu);
			mu /= r;
			mu *= mean;
			f->mu() = mu;
		}				

		std::cout << "Max mu norm " << max << " Mean " << mean << "Standard deviation " << sqrt( variance ) << std::endl;
	};


	template<typename M>
	void	CQuadrilateralHarmonicMapper<M>::_boundary_segmentation( std::vector<typename M::CVertex*>& markers )
	{
		//compute boundary loop and divide it into four segments

		M::CBoundary bound( m_pMesh );
		M::CLoop * pL = bound.loops()[0];

		std::vector<M::CVertex*> bnd_markers;
		for( int i = 0; i < 4; i ++ )
		{
			bnd_markers.push_back( markers[i] );
		}

		pL->divide( bnd_markers );

		for( size_t i = 0; i < pL->segments().size(); i ++ )
		{
			M::CLoopSegment * pS = pL->segments()[i];
			std::cout << pS->start()->id() << " " << pS->end()->id() << std::endl;

			std::vector<M::CVertex*> verts;

			for( size_t j = 0; j < pS->halfedges().size(); j ++ )
			{
				M::CHalfEdge * ph = pS->halfedges()[j];
				M::CVertex   * pv = m_pMesh->halfedgeSource( ph );
				pv->uv() = pS->start()->uv();
				verts.push_back( pv );

				std::cout << pv->point()[0] << " " << pv->point()[1] << " " << pv->point()[2] << std::endl;
			}
			verts.push_back( pS->end() );
			std::cout << pS->end()->point()[0] << " " << pS->end()->point()[1] << " " << pS->end()->point()[2] << std::endl;
			
			m_boundary_segments.push_back( verts );
		}

		
	}




	template<typename M>
	void	CQuadrilateralHarmonicMapper<M>::set_boundary_constraints( int k )
		{
			
			std::vector<M::CVertex*> * bottom; 
			std::vector<M::CVertex*> * top;    

			if( k == 0 ) 
			{
				bottom = &m_boundary_segments[3];
				top    = &m_boundary_segments[1];
				for( size_t i = 0; i < (*bottom).size(); i ++ )
				{
					M::CVertex * pv = (*bottom)[i];
					pv->u() = m_boundary_positions[0][0];
					pv->fixed() = true;
					//pv->u() = pv->uv()[0];
					//std::cout << pv->id() << " " << pv->u() << std::endl;
				}

				for( size_t i = 0; i < (*top).size(); i ++ )
				{
					M::CVertex * pv = (*top)[i];
					pv->u() = m_boundary_positions[1][0];
					pv->fixed() = true;
					//pv->u() = pv->uv()[0];
					//std::cout << pv->id() << " " << pv->u() << std::endl;
				}

			}
			else
			{
				bottom = &m_boundary_segments[0];
				top    = &m_boundary_segments[2];

				for( size_t i = 0; i < (*bottom).size(); i ++ )
				{
					M::CVertex * pv = (*bottom)[i];
					pv->u() = m_boundary_positions[0][1];
					pv->fixed() = true;
					//pv->u() = pv->uv()[1];
					//std::cout << pv->id() << " " << pv->u() << std::endl;
				}

				for( size_t i = 0; i < (*top).size(); i ++ )
				{
					M::CVertex * pv = (*top)[i];
					pv->u() = m_boundary_positions[2][1];
					pv->fixed() = true;
					//pv->u() = pv->uv()[1];
					//std::cout << pv->id() << " " << pv->u() << std::endl;
				}

			}

		};


} //namespace Holomorphy
} //namespace MeshLib

#endif

