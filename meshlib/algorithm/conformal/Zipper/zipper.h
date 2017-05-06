/*!
*      \file Zipper.h
*      \brief Algorithm for Gauss Curvature
*	   \author David Gu
*      \date Document 02/13/2014
*
*/


#ifndef _ZIPPER_H_
#define _ZIPPER_H_

#include <vector>
#include "Mesh/iterators.h"
#include "ZipperMesh.h"

#ifndef PI
#define PI 3.14159265358979323846
#endif
namespace MeshLib
{

	template<typename M>
	class CZipper
	{
	public:
		/*!	CZipper constructor
		 *	\param pMesh the input mesh
		 */
		CZipper( M* pMesh);
		/*!	CZipper destructor
		 */
		~CZipper();

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
		typename M::CVertex * operator[]( size_t k ) { return m_pts[k]; };

		/*!
		 *	output the uv coordinates
		 */
		void output( const char * fileName );

		/*!
		 *	compute the conformal map
		 */
		void _map();

	protected:
		/*!	The input surface mesh
		 */
		M* m_pMesh;
		/*	boundary
		 */
		typename M::CBoundary m_boundary;
		/*	array of boundary points
		 */
		std::vector<typename M::CVertex*> m_pts;
		int m_order;
	};


/*!	CHarmonicMapper constructor 
*	Count the number of interior vertices, boundary vertices and the edge weight
*
*/
template<typename M>
CZipper<M>::CZipper( M* pMesh ): m_pMesh( pMesh ), m_boundary( pMesh )
{
	M::CLoop * pL = m_boundary.loops()[0];

	for( std::list<M::CHalfEdge*>::iterator hiter = pL->halfedges().begin(); hiter != pL->halfedges().end(); hiter ++ )
	{
		M::CHalfEdge * pH = *hiter;
		M::CVertex   * pV = m_pMesh->halfedgeTarget( pH );
		m_pts.push_back( pV );
	}

};

/*!
 *	CZipper destructor
 */
template<typename M>
CZipper<M>::~CZipper()
{
};

/*!
 *	maps the completment of the circular arc through z0,z1,z2 to upper half plane
 */
template<typename M>
void CZipper<M>::_map_arc_completment_to_upper_half_plane()
{
	std::complex<double> zi(0,1);

	std::complex<double> z0 = m_pts[0]->z();
	std::complex<double> z1 = m_pts[1]->z();
	std::complex<double> z2 = m_pts[2]->z();

	std::complex<double> c = (z0-z1)/(z1-z2);

	m_pts[0]->z() = std::complex<double>( -1e+20,0 );
	m_pts[1]->z() = std::complex<double>( -1,0 );
	m_pts[2]->z() = std::complex<double>(  0,0 );

	for( size_t i = 3; i < m_pts.size(); i ++ )
	{
		M::CVertex * pV = m_pts[i];
		std::complex<double> z = pV->z();
		std::complex<double> w = zi * std::sqrt( c*(z - z2)/(z-z0) );
		pV->z() = w;
	}

	for( M::MeshVertexIterator viter( m_pMesh ); !viter.end(); viter ++ )
	{
		M::CVertex * pV = *viter;
		if( pV->boundary() ) continue;
		std::complex<double> z = pV->z();
		std::complex<double> w = zi * std::sqrt( c*(z - z2)/(z-z0) );
		pV->z() = w;
	}

	m_order = 3;
};

/*!
 *	maps the completment of the circular arc through z0,z1,z2 to upper half plane
 */
template<typename M>
void CZipper<M>::_zip()
{
	if( (size_t) m_order >= m_pts.size() )
	{
		std::cout << "The process has been finished !" << std::endl;
		return;
	}
	std::vector<double> x_coords;

	for( size_t i = 0; i < m_pts.size(); i ++ )
	{
		x_coords.push_back( m_pts[i]->z().real() );
	}

	std::complex<double> a = m_pts[m_order]->z();

	std::cout << m_order << "/" << m_pts.size() << ": a is " << a << std::endl;

	if( fabs( a.real()) < 1e-5 )
	{
		double c = a.imag();

		for( M::MeshVertexIterator viter( m_pMesh ); !viter.end(); viter ++ )
		{
			M::CVertex * pV = *viter;
			std::complex<double> z = pV->z();
			std::complex<double> w = z*z+c*c;
			w = std::complex<double>(0,1) * std::sqrt( -w );
			pV->z() = w;
		}

		for( size_t i = 0; i < (size_t)( m_order-1 ); i ++ )
		{
			M::CVertex * pV = m_pts[i];
			std::complex<double> w = pV->z();

			if( x_coords[i] > 0 ) 
				w = std::complex<double>( fabs( w.real() ), 0 ); 
			else
				w = std::complex<double>(-fabs( w.real() ), 0 ); 

			pV->z() = w;
		}
		
		m_pts[m_order-1]->z() = std::complex<double>( -c, 0 );
		m_pts[m_order]->z() = std::complex<double>( 0, 0 );
				

		m_order++;
		return;
	}

	double b = std::norm(a)/a.real();
	double c = std::norm(a)/a.imag();

	std::cout << a/(std::complex<double>(1.0,0)-a/b) << " " << c << std::endl;


	for( M::MeshVertexIterator viter( m_pMesh ); !viter.end(); viter ++ )
	{
		M::CVertex * pV = *viter;
		if( pV->boundary() ) continue;

		std::complex<double> z = pV->z();
		z = z/(std::complex<double>(1.0,0)-z/b);
		std::complex<double> w = z*z+c*c;
		w = std::complex<double>(0,1) * std::sqrt( -w );
		pV->z() = w;
	}

	for( size_t i = 0; i < (size_t)(m_order-1); i ++ )
	{
		M::CVertex * pV = m_pts[i];
		std::complex<double> w;
		std::complex<double> z = pV->z();
		if( std::norm( std::complex<double>(1.0,0)-z/b ) < 1e-8 )
		{
			w = std::complex<double>(1e+20,0);
		}
		else
		{
			z = z/(std::complex<double>(1.0,0)-z/b);
			w = z*z+c*c;
			w = std::complex<double>(0,1) * std::sqrt( -w );
		}

		if( b < 0 )
		{
			if( x_coords[i] < b || x_coords[i] > 0 ) 
				w = std::complex<double>( fabs( w.real() ), 0 ); 
			else
				w = std::complex<double>(-fabs( w.real() ), 0 ); 
		}
		else
		{
			if( x_coords[i] < b && x_coords[i] > 0 ) 
				w  = std::complex<double>( fabs( w.real() ), 0 );
			else
				w = std::complex<double>( -fabs( w.real() ), 0 );
		}

		pV->z() = w;
	}
	
	m_pts[m_order-1]->z() = std::complex<double>( -c, 0 );
	m_pts[m_order]->z() = std::complex<double>( 0, 0 );

	for( size_t i = m_order + 1; i < m_pts.size(); i ++ )
	{
		M::CVertex * pV = m_pts[i];
		std::complex<double> z = pV->z();

		z = z/(std::complex<double>(1.0,0)-z/b);
		std::complex<double> w = z*z+c*c;
		w = std::complex<double>(0,1) * std::sqrt( -w );
		pV->z() = w;
	}

	m_order++;
};


/*!
 *	maps the upper half plane to the unit disk
 */
template<typename M>
std::complex<double> CZipper<M>::_map_upper_half_plane_to_unit_disk(std::complex<double> z)
{
	std::complex<double> zi(0,1);

	std::complex<double> w = (z-zi)/(z+zi);
	return w;
};

/*!
 *	maps the upper half plane to the unit disk
 */
template<typename M>
std::complex<double> CZipper<M>::_map_unit_disk_to_upper_half_plane(std::complex<double> z)
{
	std::complex<double> zi(0,1);

	std::complex<double> w = zi*(1.0+z)/(1.0-z);
	return w;
};

/*!
 *	maps the upper half plane to the unit disk
 */
template<typename M>
void CZipper<M>::_normalize()
{
	std::complex<double> w(0,0);

	for( M::MeshVertexIterator viter( m_pMesh ); !viter.end(); viter ++ )
	{
		M::CVertex * pV = *viter;
		std::complex<double> zeta = _map_upper_half_plane_to_unit_disk( pV->z() );
		w += zeta;
	}
	w /= m_pMesh->numVertices();

	w = _map_unit_disk_to_upper_half_plane( w );

	std::complex<double> r = std::complex<double>(0,1)/w;
	
	double a = r.real();
	double b = -r.imag();

	w = a * w /( b * w + 1.0);
	std::cout << w << std::endl;

	for( M::MeshVertexIterator viter( m_pMesh ); !viter.end(); viter ++ )
	{
		M::CVertex * pV = *viter;
		std::complex<double> z = pV->z();
		z = a * z /( b * z + 1.0 );
		pV->z() = z;
	}

};


/*!
 *	output mesh
 */
template<typename M>
void CZipper<M>::output( const char * file_name )
{
	for( M::MeshVertexIterator viter( m_pMesh ); !viter.end(); viter ++ )
	{
		M::CVertex * pV = *viter;
		pV->z() = _map_upper_half_plane_to_unit_disk( pV->z() );
	}

	m_pMesh->write_m( file_name );
};

/*!
 *	compute the conformal map
 */
template<typename M>
void CZipper<M>::_map()
{
	_map_arc_completment_to_upper_half_plane();
	while( m_order < m_pts.size() )
	{
		_zip();
		_normalize();
	}

	for( int i = 0;i < 8; i ++ )
	{
		_normalize();
	}
};

};
#endif

