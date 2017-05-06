/*! \file Octree.h
*   \brief Octree
*   \author David Gu
*   \date   documented on 02/11/2011
*
*   Mesh for viewer 
*/
#ifndef  _OCTREE_H_
#define  _OCTREE_H_

#include <map>
#include <vector>
#include <cmath>
#include <iostream>
#include <limits>
#include <cstdlib>
#include <ctime>

#include "Mesh/Vertex.h"
#include "Mesh/Edge.h"
#include "Mesh/Face.h"
#include "Mesh/HalfEdge.h"

#include "Mesh/BaseMesh.h"
#include "Mesh/boundary.h"
#include "Mesh/iterators.h"
#include "Parser/parser.h"
#include "TriangleCubeIntersect.h"

namespace MeshLib
{
	class CSample
	{
	public:
		CPoint  & point() { return m_point; };
		CPoint2 & uv()    { return m_uv;    };
		CPoint  & rgb()   { return m_rgb;   };
	protected:
		CPoint  m_point;
		CPoint2 m_uv;
		CPoint  m_rgb;
	};

	class CTriangle
	{
	public:
		CTriangle() {};
		CTriangle( const CTriangle & tri )
		{
			for( int i = 0; i < 3; i ++ )
			{
				v[i] = tri.v[i]; uv[i] = tri.uv[i]; rgb[i] = tri.rgb[i];
			}
		}
		~CTriangle(){};
	public:
		CPoint    v[3];
		CPoint2  uv[3];
		CPoint  rgb[3];
		//project point in, the output is out, 
		//if the projection is out of the triangle
		// return false; else retun true
		bool _project( CPoint & in, CPoint & out  );
		bool _project( CPoint & in, CSample & out );
	};

	/*!	COctreeNode class
	 *
	 */
	class COctreeNode
	{
	public:
		COctreeNode( CPoint p, CPoint q )
		{
			m_corner[0] = p;
			m_corner[1] = q;
			
			for( int i = 0; i < 2; i ++ )
			for( int j = 0; j < 2; j ++ )
			for( int k = 0; k < 2; k ++ )
			{
				m_child[i][j][k] = NULL;
			}
		};

		~COctreeNode()
		{
			for( int i = 0; i < 2; i ++ )
			for( int j = 0; j < 2; j ++ )
			for( int k = 0; k < 2; k ++ )
			{
				if( m_child[i][j][k] != NULL) delete m_child[i][j][k];
			}

		}

		void subdivide( std::vector<CPoint> & pts, int n )
		{
			if( n <= 0 || pts.empty() ) return;

			double len = ( m_corner[1][0] - m_corner[0][0] )/2.0;
			
			for( size_t id = 0; id < pts.size(); id ++ )
			{
				CPoint p = pts[id];
				bool inside = true;
				for( int i = 0; i < 3; i ++ )
				{
					if( p[i] < m_corner[0][i] || p[i] > m_corner[1][i] ) inside = false; 
				}
				if( inside )
					m_pts.push_back( p );
			}


			for( int i = 0; i < 2; i ++ )
			for( int j = 0; j < 2; j ++ )
			for( int k = 0; k < 2; k ++ )
			{
				CPoint p = m_corner[0] + CPoint( i, j, k ) * len;
				CPoint q = p + CPoint(len, len, len );
				m_child[i][j][k] = new COctreeNode(p,q);

				m_child[i][j][k]->subdivide( m_pts, n-1 );
			}
						
		}
		
		void subdivide( std::vector<CTriangle*> & trs, int n )
		{
			if( n < 0 || trs.empty() ) return;

			CTriangleCubeIntersect TC;

			double len = ( m_corner[1][0] - m_corner[0][0] )/2.0;

			CPoint center = (m_corner[0] + m_corner[1] )/2.0;
			
			for( size_t id = 0; id < trs.size(); id ++ )
			{
				CTriangle * pT = trs[id];
				
				CPoint a = ( pT->v[0] - center )/(len*2);
				CPoint b = ( pT->v[1] - center )/(len*2);
				CPoint c = ( pT->v[2] - center )/(len*2);

				long inside = TC.test( a, b, c );
				if( !inside )
					m_triangles.push_back( pT );
			}
			
			//current node is a leaf
			if( n == 0 ) return;
	
			for( int i = 0; i < 2; i ++ )
			for( int j = 0; j < 2; j ++ )
			for( int k = 0; k < 2; k ++ )
			{
				CPoint p = m_corner[0] + CPoint( i, j, k ) * len;
				CPoint q = p + CPoint(len, len, len );
				m_child[i][j][k] = new COctreeNode(p,q);

				m_child[i][j][k]->subdivide( m_triangles, n-1 );
			}
						
		}

		void _sample( std::vector<CPoint>  & samples );
		void _sample( std::vector<CSample> & samples );

	public:

		CPoint m_corner[2];
		COctreeNode * m_child[2][2][2];
		std::vector<CPoint> m_pts;
		std::vector<CTriangle*> m_triangles;	
	};

	/*!	COctree class
	 *
	 */
	class COctree
	{
	public:
		COctree();
		~COctree();
		COctreeNode * root() { return m_root; };
		void _construct( std::vector<CPoint> & pts, int n );
		void _construct( int n );
		void _insert_triangle( CPoint a, CPoint b, CPoint c );
		void _insert_triangle( const CTriangle & tri );
		void _sample( std::vector<CPoint> & samples );
		void _sample( std::vector<CSample> & samples );

	protected:
		COctreeNode * m_root;
		std::vector<CTriangle*> m_trs;
	};
};

#endif