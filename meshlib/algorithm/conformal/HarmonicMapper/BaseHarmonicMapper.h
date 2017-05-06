/*!
*      \file BaseHarmonicMapper.h
*      \brief Algorithm for harmonic mapping
*	   \author David Gu
*      \date Document 09/29/2013
*
*		Base harmonic map, that maps a topological quadrilateral to topological quadrilateral
*		with minimal harmonic energy.
*/


#ifndef _BASE_HARMONIC_MAPPER_H_
#define _BASE_HARMONIC_MAPPER_H_

#include <vector>
#include "TMapMesh.h"
#include "Operator/Operator.h"

namespace MeshLib
{

namespace Holomorphy
{

/*!
 *	\brief CBaseHarmonicMapper class
 *
 *  Compute the harmonic map by solving Dirichlet problem
 * \f[
 *		\left\{
 *		\begin{array}{lcc}
 *		 \Delta u &\equiv& 0\\
 *		 u|_{\partial \Omega} &=& f
 *		 \end{array}
 *		\right.
 * \f]
 */

	template<typename M>
	class CBaseHarmonicMapper
	{
	public:
		/*!	CBaseHarmonicMapper constructor
		 *	\param pMesh the input mesh
		 */
		CBaseHarmonicMapper(M * pMesh );
		/*!	CBaseHarmonicMapper destructor
		 */
		~CBaseHarmonicMapper();

		/*!	input the interior markers
		 */
		virtual void set_interior_markers( std::vector<typename M::CVertex*> & markers, std::vector<CPoint2> & positions )
		{
			m_interior_markers = markers;
			m_interior_positions = positions;
		};

		/*!	input the boundary markers
		 */
		virtual void set_boundary_markers( std::vector<typename M::CVertex*> & markers, std::vector<CPoint2> & positions )
		{
			m_boundary_markers = markers;
			m_boundary_positions = positions;
		};

		/*! Teichmuller map
		 */
		void _tmap();

		/*!	qc_harmonic_map
		 */
		void _qc_harmonic_map();

	protected:

		/*!	compute edge weight
		 */
		void _initialize_edge_weight();
		/*!	update edge weight
		 */

		void _update_edge_weight();
		virtual void _manipulate_mu();
		/*!	The input surface mesh
		 */
		M* m_pMesh;
		
		/*!
		 *	compute the harmonic functioon defined on the mesh, the input mesh 
		 *  1. vertex->fixed() field has been set, for anchors, this field is true, for unknowns, this field is false
		 *  2. for anchors, the vertex->u() is set as the Dirichlet condition
		 */
		void _harmonic_function();
		
		/*! The collection of interior markers
		 */
		std::vector<typename M::CVertex *> m_interior_markers;
		/*!	The corresponding positions of the interior markers
		 */
		std::vector<CPoint2> m_interior_positions;

		/*!	The boundary markers
		 */
		std::vector<typename M::CVertex *> m_boundary_markers;
		/*!	The corresponding positions of the boundary markers
		 */
		std::vector<CPoint2> m_boundary_positions;

	protected:

		virtual void set_boundary_constraints( int k ) = 0;
		virtual void set_interior_constraints( int k );
		/*!
		 *	compute a harmonic mapping using current \mu
		 */
		void _map();
	};

	template<typename M>
	CBaseHarmonicMapper<M>::CBaseHarmonicMapper( M * pMesh )
	{
		m_pMesh = pMesh;
	};

	template<typename M>
	CBaseHarmonicMapper<M>::~CBaseHarmonicMapper()
	{
	};

	template<typename M>
	void	CBaseHarmonicMapper<M>::_initialize_edge_weight()
	{
		//Compute cotangent edge weight
		COperator<M> pC( m_pMesh );
		pC._embedding_2_metric();  //convert embedding to metric
		pC._metric_2_angle();	   //convert metric to angle
		pC._angle_2_Laplace();	   //convert angle to cotangent edge weight
				
	};

	template<typename M>
	void	CBaseHarmonicMapper<M>::_update_edge_weight()
	{
		//Compute cotangent edge weight
		COperator<M> pC( m_pMesh );
		pC._parameter_2_mu();  //convert embedding to metric

		_manipulate_mu();

		pC._parameter_mu_2_angle();	   //convert metric to angle
		pC._angle_2_Laplace();	   //convert angle to cotangent edge weight
				
	};


	template<typename M>
	void	CBaseHarmonicMapper<M>::_manipulate_mu()
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
	void CBaseHarmonicMapper<M>::_harmonic_function()
	{
		CLaplace<M> L( m_pMesh );
		L.solve();
	};


	template<typename M>
	void CBaseHarmonicMapper<M>::set_interior_constraints( int k )
	{
		for( size_t i = 0; i < m_interior_markers.size(); i ++ )
		{
			M::CVertex * pV = m_interior_markers[i];
			pV->fixed() = true;
			pV->u() = m_interior_positions[i][k];
		}
	};


	template<typename M>
	void	CBaseHarmonicMapper<M>::_map()
	{
		for( int k = 0; k < 2; k ++ )
		{
			for( M::MeshVertexIterator viter( m_pMesh ); !viter.end(); viter ++ )
			{
				M::CVertex * pv = *viter;
				pv->fixed() = false;
			}
			set_interior_constraints(k);
			set_boundary_constraints(k);
			_harmonic_function();
			for( M::MeshVertexIterator viter( m_pMesh ); !viter.end(); viter ++ )
			{
				M::CVertex * pv = *viter;
				pv->uv()[k] = pv->u();
			}

		}

		for( M::MeshVertexIterator viter( m_pMesh ); !viter.end(); viter ++ )
		{
			M::CVertex * pv = *viter;
			CPoint2 uv = pv->uv();
			pv->w() = std::complex<double>( uv[0], uv[1] );
		}
	};
	
	template<typename M>
	void	CBaseHarmonicMapper<M>::_tmap()
	{
		_initialize_edge_weight();
		_map();
		
		for( int i = 0; i < 16; i ++ )
		{
			_update_edge_weight();
			_map();
		}
	};

	//assume the vertex->z() is set, the face->mu() is set, boundary constraints are set

	template<typename M>
	void	CBaseHarmonicMapper<M>::_qc_harmonic_map()
	{
		//Compute cotangent edge weight
		COperator<M> pC( m_pMesh );
		pC._parameter_mu_2_angle();	   //convert metric to angle
		pC._angle_2_Laplace();	   //convert angle to cotangent edge weight
		_map();

		for( M::MeshFaceIterator fiter( m_pMesh ); !fiter.end(); fiter ++ )
		{
			M::CFace * pF = *fiter;
			pF->nu() = pF->mu();
		}

		pC._parameter_2_mu();

		for( M::MeshFaceIterator fiter( m_pMesh ); !fiter.end(); fiter ++ )
		{
			M::CFace * pF = *fiter;
			double  err = std::abs( pF->nu() - pF->mu() );
			if( err < 0.1 ) continue;
			std::cout << pF->id() << " " << err << std::endl;
		}

	}

} //namespace Holomorphy
} //namespace MeshLib
#endif

