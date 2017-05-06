/*!
*      \file ConformalWelding.h
*      \brief Algorithm for conformal welding using Zipper Algorithm
*	   \author David Gu
*      \date Document 02/13/2014
*
*/


#ifndef _CONFORMAL_WELDING_H_
#define _CONFORMAL_WELDING_H_

#include <vector>
#include <deque>
#include "Mesh/iterators.h"
#include "ZipperMesh.h"

namespace MeshLib
{

	template<typename M>
	class CConformalWelding
	{
	public:
		/*!	
		 *  CConformalWelding constructor
		 *	\param pMesh0 - interior mesh
		 *  \param pMesh1 - exterior mesh
		 */
		CConformalWelding( M * pMesh0, M * pMesh1 );
		/*!	
		 *  CConformalWelding destructor
		 */
		~CConformalWelding();

		/*!
		 *	maps the complement of the arc (z_0, z_1, z_2) to the upper half plane
		 */
		void _map_arc_completment_to_upper_half_plane();
	    /*!
		 *	upper half plane to the unit disk
		 */
		std::complex<double> _map_upper_half_plane_to_unit_disk(std::complex<double> z);
		/*!
		 *	unit disk to upper half plane
		 */
		
		std::complex<double> _map_unit_disk_to_upper_half_plane(std::complex<double> z);
		/*!
		 *	normalize
		 */
		
		void _normalize();
		/*!
		 *	unzip one segment
		 */
		void _zip();

		/*!
		 *	boundary vertices
		 */
		typename M::CVertex * operator[]( size_t k ) { return m_pts0[k]; };
		typename M::CVertex * operator()( size_t k ) { return m_pts1[k]; };

		/*!
		 *	output the uv coordinates
		 */
		void output( const char * fileName );

		/*!
		 *	compute the conformal map
		 */
		void _map();
		/*!
		 *  map z0,z1,z2 and w0, w1, w2 to -\infty, -1, 0
		 */
		void initialize();
		/*!
		 *	final step
		 */
		void _step_2();
		/*!
		 *	step 3
		 */
		void _step_3();
		/*!
		 *	step 4
		 */
		void _step_4();
		/*!
		 *	step 5
		 */
		void _step_5();
		/*!
         *  size of boundary vertices
		 */
		size_t size() { return m_pts0.size(); };

	protected:
		/*!
		 *  The input interior mesh
		 */
		M* m_pMesh0;
		/*!
		 *	The input exterior mesh
		 */
		M* m_pMesh1;
		/*	
		 * interior boundary
		 */
		typename M::CBoundary m_boundary0;
		/*
		 *  exteriior boundary
		 */
		typename M::CBoundary m_boundary1;
		/*
		 * array of boundary points of interior mesh
		 */
		std::deque<typename M::CVertex*> m_pts0;
		/*
		 * array of boundary points of exterior mesh
		 */
		std::deque<typename M::CVertex*> m_pts1;

		/*
		 *	next zero point	
		 */
		int m_order;

	};


/*!	CHarmonicMapper constructor 
*	Count the number of interior vertices, boundary vertices and the edge weight
*
*/
template<typename M>
CConformalWelding<M>::CConformalWelding( M* pMesh0, M * pMesh1 ): m_pMesh0( pMesh0 ), m_pMesh1( pMesh1), m_boundary0( m_pMesh0 ), m_boundary1( m_pMesh1 )
{
	M::CLoop * pL0 = m_boundary0.loops()[0];

	for( std::list<M::CHalfEdge*>::iterator hiter = pL0->halfedges().begin(); hiter != pL0->halfedges().end(); hiter ++ )
	{
		M::CHalfEdge * pH = *hiter;
		M::CVertex   * pV = m_pMesh0->halfedgeTarget( pH );
		m_pts0.push_back( pV );
	}

	M::CLoop * pL1 = m_boundary1.loops()[0];

	for( std::list<M::CHalfEdge*>::iterator hiter = pL1->halfedges().begin(); hiter != pL1->halfedges().end(); hiter ++ )
	{
		M::CHalfEdge * pH = *hiter;
		M::CVertex   * pV = m_pMesh1->halfedgeTarget( pH );
		m_pts1.push_front( pV );
	}

	while( true )
	{
		M::CVertex * pV = m_pts1[0];
		M::CVertex * pW = m_pts0[0];
		if( pV->id() == pW->id() ) break;
		m_pts1.pop_front();
		m_pts1.push_back( pV );
	}

	for( size_t i = 0; i < m_pts0.size(); i ++ )
	{
		M::CVertex * pV = m_pts1[i];
		M::CVertex * pW = m_pts0[i];

		std::cout << pV->id() << " " << pW->id() << std::endl;
	}
};

/*!
 *	CConformalWelding destructor
 */
template<typename M>
CConformalWelding<M>::~CConformalWelding()
{
};



/*!
 *	maps the upper half plane to the unit disk
 */
template<typename M>
std::complex<double> CConformalWelding<M>::_map_upper_half_plane_to_unit_disk(std::complex<double> z)
{
	std::complex<double> zi(0,1);

	std::complex<double> w = (z-zi)/(z+zi);
	return w;
};

/*!
 *	maps the upper half plane to the unit disk
 */
template<typename M>
std::complex<double> CConformalWelding<M>::_map_unit_disk_to_upper_half_plane(std::complex<double> z)
{
	if( std::norm( z - 1.0 ) < 1e-8 ) return std::complex<double>(1e+20,0);

	std::complex<double> zi(0,1);

	std::complex<double> w = zi*(1.0+z)/(1.0-z);
	return w;
};

/*!
 *	initialize, maps (z0,z1,z2) and (w0,w1,w2) to (-\infty, -1, 0 )
 */
template<typename M>
void CConformalWelding<M>::initialize()
{
	size_t n = m_pts0.size();

	M::CVertex * pV2 = m_pts0[n-1];
	M::CVertex * pV1 = m_pts0[n-2];
	M::CVertex * pV0 = m_pts0[n-3];

	std::complex<double> z0 = pV0->z();
	std::complex<double> z1 = pV1->z();
	std::complex<double> z2 = pV2->z();

	for( size_t i = 0; i < n-3; i ++ )
	{
		M::CVertex * pV = m_pts0[i];
		std::complex<double> z = pV->z();
		z = (z2-z0)/(z0-z1)*(z-z1)/(z-z2);
		pV->z() = z;
	}
	pV0->z() = std::complex<double>(-1,0);
	pV1->z() = std::complex<double>( 0,0);
	pV2->z() = std::complex<double>(1e+20,0);

	for( M::MeshVertexIterator viter( m_pMesh0 ); !viter.end(); viter ++ )
	{
		M::CVertex * pV = *viter;
		if( pV->boundary() ) continue;
		std::complex<double> z = pV->z();
		z = (z2-z0)/(z0-z1)*(z-z1)/(z-z2);
		pV->z() = z;
	}

	M::CVertex * pW2 = m_pts1[n-1];
	M::CVertex * pW1 = m_pts1[n-2];
	M::CVertex * pW0 = m_pts1[n-3];

	std::complex<double> w0 = pW0->z();
	std::complex<double> w1 = pW1->z();
	std::complex<double> w2 = pW2->z();

	for( size_t i = 0; i < n-3; i ++ )
	{
		M::CVertex * pW = m_pts1[i];
		std::complex<double> w = pW->z();
		w = (w2-w0)/(w0-w1)*(w-w1)/(w-w2);
		pW->z() = w;
	}
	pW0->z() = std::complex<double>(-1,0);
	pW1->z() = std::complex<double>( 0,0);
	pW2->z() = std::complex<double>(1e+20,0);

	for( M::MeshVertexIterator viter( m_pMesh1 ); !viter.end(); viter ++ )
	{
		M::CVertex * pW = *viter;
		if( pW->boundary() ) continue;
		std::complex<double> w = pW->z();
		w = (w2-w0)/(w0-w1)*(w-w1)/(w-w2);
		pW->z() = w;
	}

	//
	m_order = n-3;

	for( M::MeshVertexIterator viter( m_pMesh0 ); !viter.end(); viter ++)
	{
		M::CVertex * pV = *viter;
		std::complex<double> z = pV->z();
		z += std::complex<double>(1,0);
		std::complex<double> zeta = std::sqrt( z );
		if( std::arg( zeta ) < 0 ) zeta = -zeta;
		z = std::complex<double>(0,1) * zeta;
		pV->z() = z;
	}

	for(size_t i = 0; i < m_order; i ++ )
	{
		M::CVertex * pV = m_pts0[i];
		std::complex<double> z = pV->z();
		double x = std::abs( z );
		pV->z() = std::complex<double>( -x, 0 );
	}

	for( M::MeshVertexIterator viter( m_pMesh1 ); !viter.end(); viter ++)
	{
		M::CVertex * pW = *viter;
		std::complex<double> w = pW->z();
		w += std::complex<double>(1,0);
		std::complex<double> zeta = std::sqrt( w );
		if( std::arg( zeta ) > 0 ) zeta = -zeta;
		w = std::complex<double>(0,1) * zeta;
		pW->z() = w;
	}

	for(size_t i = 0; i < m_order; i ++ )
	{
		M::CVertex * pW = m_pts1[i];
		std::complex<double> z = pW->z();
		double x = std::abs( z );
		pW->z() = std::complex<double>( x, 0 );
	}

	m_order --;
};

/*!
 *	initialize, maps (z0,z1,z2) and (w0,w1,w2) to (-\infty, -1, 0 )
 */
template<typename M>
void CConformalWelding<M>::_zip()
{
	if( m_order < 0 ) return;

	std::cout << m_pts0[m_order]->z() << " " << m_pts1[m_order]->z() << std::endl;
	
	//az/(cz+1)

	double alpha = m_pts0[m_order]->z().real();
	double beta  = m_pts1[m_order]->z().real();

	double m[2][2];

	m[0][0] = 1.0/4.0;
	m[0][1] = 1.0/4.0;
	m[1][0] = 1.0/2.0;
	m[1][1] =-1.0/2.0;

	double b[2];
	b[0] = -1.0/alpha;
	b[1] =  1.0/beta;

	double a[2];

	a[0] = m[0][0] * b[0] + m[0][1] * b[1];
	a[1] = m[1][0] * b[0] + m[1][1] * b[1];

	std::complex<double> z = std::complex<double>(alpha,0);
	std::cout << a[0] * z /( a[1] * z + 1.0 ) << " ";
	z = std::complex<double>(beta,0);
	std::cout << a[0] * z /( a[1] * z + 1.0 ) << std::endl;

/*
	for( M::MeshVertexIterator viter( m_pMesh0 ); !viter.end(); viter ++ )
	{
		M::CVertex * pV = *viter;
		std::complex<double> z = pV->z();
		z = a[0] * z /( a[1] * z + 1.0 );
		if( std::arg( pV->z() ) < 0 )
		{
			std::cout << std::arg( pV->z() );
		}
		pV->z() = z * 0.2;
	}
*/
	for( M::MeshVertexIterator viter( m_pMesh0 ); !viter.end(); viter ++ )
	{
		M::CVertex * pW = *viter;
		std::complex<double> w = pW->z();
		w = a[0] * w /( a[1] * w + 1.0 );
		/*
		if( std::arg( pW->z() ) < 0 )
		{
			std::cout << w << " " << std::arg(pW->z() ) << std::endl;
		}
		*/
		pW->z() = w * 0.2;
	}

	for( M::MeshVertexIterator viter( m_pMesh1 ); !viter.end(); viter ++ )
	{
		M::CVertex * pW = *viter;
		std::complex<double> w = pW->z();
		w = a[0] * w /( a[1] * w + 1.0 );
		/*
		if( std::arg( pW->z() ) < 0 )
		{
			std::cout << w << " " << std::arg(pW->z() ) << std::endl;
		}
		*/
		pW->z() = w * 0.2;
	}


/*

	for( M::MeshVertexIterator viter( m_pMesh0 ); !viter.end(); viter ++ )
	{
		M::CVertex * pV = *viter;

		std::complex<double> z0 = pV->z();
		std::complex<double> z = pV->z();

		z = z * z - 0.01;		
		z = std::complex<double>(0,1) * std::sqrt( -z );

		std::complex<double> w0 = this->_map_upper_half_plane_to_unit_disk( z0 );

		if( std::norm( w0 ) > 0.98 )
		{
			if( z0.real() * z.real() < 0 ) z = -z;
		}	

		pV->z() = z;
	}
*/
	for( M::MeshVertexIterator viter( m_pMesh0 ); !viter.end(); viter ++ )
	{
		M::CVertex * pV = *viter;

		std::complex<double> z0 = pV->z();

		double a0 = std::arg( z0 );
		std::complex<double> z = pV->z();
		z = z * z - 0.01;

		z = std::complex<double>(0,1) * std::sqrt( -z );

		std::complex<double> w0 = this->_map_upper_half_plane_to_unit_disk( z0 );
		if( std::norm( w0 ) > 0.98 )
		{
			if( z0.real() * z.real() < 0 ) z = -z;
		}	

		pV->z() = z;
	}


	for( M::MeshVertexIterator viter( m_pMesh1 ); !viter.end(); viter ++ )
	{
		M::CVertex * pV = *viter;

		std::complex<double> z0 = pV->z();

		double a0 = std::arg( z0 );
		std::complex<double> z = pV->z();
		z = z * z - 0.01;

		z = std::complex<double>(0,1) * std::sqrt( -z );

		std::complex<double> w0 = this->_map_upper_half_plane_to_unit_disk( z0 );
		if( std::norm( w0 ) > 0.98 )
		{
			if( z0.real() * z.real() < 0 ) z = -z;
		}	

		pV->z() = z;
	}


	for( size_t i = 0; i < m_order; i ++ )
	{
		m_pts0[i]->z() = std::complex<double>( m_pts0[i]->z().real(), 0);
		m_pts1[i]->z() = std::complex<double>( m_pts1[i]->z().real(), 0);
	}

	m_pts0[m_order]->z() = std::complex<double>(0,0);
	m_pts1[m_order]->z() = std::complex<double>(0,0);
	m_order --;
	std::cout << "Order " << m_order << std::endl;

	size_t n = m_pts0.size();
	std::complex<double> Zend = m_pts0[n-1]->z();
	std::cout << "Zend " << Zend << std::endl;

/*
	for( M::MeshVertexIterator viter( m_pMesh0 ); !viter.end(); viter ++ )
	{
		M::CVertex * pV = *viter;

		if( pV == m_pts0[n-1] )
		{
			pV->z() = 1e+20;
			continue;
		}

		std::complex<double> z0 = pV->z();

		std::complex<double> z = pV->z() /( pV->z() - Zend );

		pV->z() = z;
	}

	for( M::MeshVertexIterator viter( m_pMesh1 ); !viter.end(); viter ++ )
	{
		M::CVertex * pV = *viter;

		if( pV == m_pts1[n-1] )
		{
			pV->z() = 1e+20;
			continue;
		}

		std::complex<double> z0 = pV->z();

		std::complex<double> z = pV->z() /( pV->z() - Zend );

		pV->z() = z;
	}
*/

}


/*!
 *	maps the upper half plane to the unit disk
 */
template<typename M>
void CConformalWelding<M>::_normalize()
{
	std::complex<double> w(0,0);

	for( M::MeshVertexIterator viter( m_pMesh0 ); !viter.end(); viter ++ )
	{
		M::CVertex * pV = *viter;
		std::complex<double> zeta = _map_upper_half_plane_to_unit_disk( pV->z() );
		w += zeta;
	}
	for( M::MeshVertexIterator viter( m_pMesh1 ); !viter.end(); viter ++ )
	{
		M::CVertex * pV = *viter;
		std::complex<double> zeta = _map_upper_half_plane_to_unit_disk( pV->z() );
		w += zeta;
	}
	w /= (m_pMesh0->numVertices()+m_pMesh1->numVertices());

	w = _map_unit_disk_to_upper_half_plane( w );

	std::complex<double> r = std::complex<double>(0,1)/w;
	
	double a = r.real();
	double b = -r.imag();

	w = a * w /( b * w + 1.0);
	std::cout << w << std::endl;

	for( M::MeshVertexIterator viter( m_pMesh0 ); !viter.end(); viter ++ )
	{
		M::CVertex * pV = *viter;
		std::complex<double> z = pV->z();
		z = a * z /( b * z + 1.0 );
		pV->z() = z;
	}
	for( M::MeshVertexIterator viter( m_pMesh1 ); !viter.end(); viter ++ )
	{
		M::CVertex * pV = *viter;
		std::complex<double> z = pV->z();
		z = a * z /( b * z + 1.0 );
		pV->z() = z;
	}

};

/*!
 *	maps the upper half plane to the unit disk
 */
template<typename M>
void CConformalWelding<M>::_step_2()
{
	size_t n = m_pts0.size();
	std::complex<double> Zend = m_pts0[n-1]->z();

	for( M::MeshVertexIterator viter( m_pMesh0 ); !viter.end(); viter ++ )
	{
		M::CVertex * pV = *viter;
		if( pV->id() == m_pts0[n-1]->id() )
		{
			pV->z() = std::complex<double>( -1e+20, 0 );
			continue;
		}

		std::complex<double> z = pV->z();

		z = z/(z-Zend);
		//z = z * z;
		//z = 0.01 * z/(z-3.0);

		pV->z() = z;
	}

	for( M::MeshVertexIterator viter( m_pMesh1 ); !viter.end(); viter ++ )
	{
		M::CVertex * pV = *viter;
		if( pV->id() == m_pts1[n-1]->id() )
		{
			pV->z() = std::complex<double>( 1e+20, 0 );
			continue;
		}
		std::complex<double> z = pV->z();

		z = z/(z-Zend);
		//z = z * z;
		//z = 0.01 * z/(z-3.0);
		pV->z() = z;
	}
};

/*!
 *	maps the upper half plane to the unit disk
 */
template<typename M>
void CConformalWelding<M>::_step_3()
{
	size_t n = m_pts0.size();
	std::complex<double> Zend = m_pts0[n-1]->z();

	for( M::MeshVertexIterator viter( m_pMesh0 ); !viter.end(); viter ++ )
	{
		M::CVertex * pV = *viter;
		if( pV->id() == m_pts0[n-1]->id() )
		{
			pV->z() = std::complex<double>( 1e+20, 0 );
			continue;
		}

		std::complex<double> z = pV->z();
		pV->z() = z * z;
	}

	for( M::MeshVertexIterator viter( m_pMesh1 ); !viter.end(); viter ++ )
	{
		M::CVertex * pV = *viter;
		if( pV->id() == m_pts1[n-1]->id() )
		{
			pV->z() = std::complex<double>( 1e+20, 0 );
			continue;
		}
		std::complex<double> z = pV->z();
		pV->z() = z * z;
	}
};

/*!
 *	maps the upper half plane to the unit disk
 */
template<typename M>
void CConformalWelding<M>::_step_5()
{
	std::complex<double> w(0,0);

	for( M::MeshVertexIterator viter( m_pMesh0 ); !viter.end(); viter ++ )
	{
		M::CVertex * pV = *viter;
		w += pV->z();
	}
	w /= (m_pMesh0->numVertices());


	for( M::MeshVertexIterator viter( m_pMesh0 ); !viter.end(); viter ++ )
	{
		M::CVertex * pV = *viter;
		std::complex<double> z = pV->z();
		z = ( z - w)/( 1.0 - std::conj( w) * z );
		pV->z() = z;
	}
	for( M::MeshVertexIterator viter( m_pMesh1 ); !viter.end(); viter ++ )
	{
		M::CVertex * pV = *viter;
		std::complex<double> z = pV->z();
		z = ( z - w)/( 1.0 - std::conj( w) * z );
		pV->z() = z;
	}

};

template<typename M>
void CConformalWelding<M>::_step_4()
{
	size_t n = m_pts0.size();

	std::complex<double> z0 = m_pts0[0]->z();
	std::complex<double> z1 = m_pts0[n/2]->z();
	std::complex<double> z2 = m_pts0[n-1]->z();

	for( M::MeshVertexIterator viter( m_pMesh0 ); !viter.end(); viter ++ )
	{
		M::CVertex * pV = *viter;
		std::complex<double> z = pV->z();
		pV->z() = z/z1;
	}

	for( M::MeshVertexIterator viter( m_pMesh1 ); !viter.end(); viter ++ )
	{
		M::CVertex * pV = *viter;
		std::complex<double> z = pV->z();
		pV->z() = z/z1;
	}

	for( M::MeshVertexIterator viter( m_pMesh0 ); !viter.end(); viter ++ )
	{
		M::CVertex * pV = *viter;
		std::complex<double> z = pV->z();
		pV->z() = _map_upper_half_plane_to_unit_disk( z );
	}

	for( M::MeshVertexIterator viter( m_pMesh1 ); !viter.end(); viter ++ )
	{
		M::CVertex * pV = *viter;
		std::complex<double> z = pV->z();
		pV->z() = _map_upper_half_plane_to_unit_disk( z );
	}

};


template<typename M>
void CConformalWelding<M>::_map()
{
	initialize();
	while( m_order >= 0 )
		_zip();

	_step_2();
	_step_3();
	_step_4();
	_step_5();
}
};
#endif

