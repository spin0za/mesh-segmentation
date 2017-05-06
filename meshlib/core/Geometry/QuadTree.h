/*! \file QuadTree.h
*   \brief Octree
*   \author David Gu
*   \date   documented on 02/14/2011
*
*   Mesh for viewer 
*/
#ifndef  _QUAD_TREE_H_
#define  _QUAD_TREE_H_

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

	/*!	CQuadTreeNode class
	 *
	 */
	template <typename T>
	class CQuadTreeNode
	{
	public:
		CQuadTreeNode<T>( CPoint2 p, CPoint2 q )
		{
			m_corner[0] = p;
			m_corner[1] = q;
			
			for( int i = 0; i < 2; i ++ )
			for( int j = 0; j < 2; j ++ )
			{
				m_child[i][j] = NULL;
			}
		};

		~CQuadTreeNode<T>()
		{
			for( int i = 0; i < 2; i ++ )
			for( int j = 0; j < 2; j ++ )
			{
				if( m_child[i][j] != NULL) delete m_child[i][j];
				m_child[i][j] = NULL;
			}

		}

		std::vector<T*> & samples() { return m_samples; };

		void subdivide( std::vector<T*> & samples, int n )
		{
			if( n < 0 || samples.empty() ) return;
			
			double len = ( m_corner[1][0] - m_corner[0][0] )/2.0;

			CPoint2 center = (m_corner[0] + m_corner[1] )/2.0;
			
			for( size_t id = 0; id < samples.size(); id ++ )
			{
				T * pT = samples[id];
				if( pT->intersect( m_corner[0], m_corner[1] ) ) 
					m_samples.push_back( pT ); 	
			}
			
#if 0
			//for debugging purpose
			if( n == 1 && m_samples.empty() )
			{
				int a = 0;

				for( size_t id = 0; id < samples.size(); id ++ )
				{
					T * pT = samples[id];
					if( pT->intersect( m_corner[0], m_corner[1] ) ) 
						m_samples.push_back( pT ); 	
				}
			}
#endif

			//current node is a leaf
			if( n == 0 || m_samples.empty() ) return;
	
			for( int i = 0; i < 2; i ++ )
			for( int j = 0; j < 2; j ++ )
			{
				CPoint2 p = m_corner[0] + CPoint2( i, j ) * len;
				CPoint2 q = p + CPoint2(len, len );
				m_child[i][j] = new CQuadTreeNode(p,q);
				if( m_child[i][j] == NULL )
				{
					std::cerr<< "Error in allocation" << std::endl;
				}
				m_child[i][j]->subdivide( m_samples, n-1 );
			}
		
			size_t c = 0;
			
			for( int i = 0; i < 2; i ++ )
			for( int j = 0; j < 2; j ++ )
			{
				c += m_child[i][j]->samples().size();
			}

			if( !m_samples.empty() && c == 0 )
			{
					std::cerr<< "Error in subdivision" << std::endl;				
			}

		};

		T * _locate( CPoint2 & pt )
		{
			if( m_child[0][0] == NULL )
			{
				for( size_t i = 0; i < m_samples.size(); i ++ )
				{
					T* pT = m_samples[i];
					if( pT->inside( pt ) ) 
						return pT;
				}

				return NULL;
			}

			double len = ( m_corner[1][0] - m_corner[0][0] )/2.0;
			CPoint2 center = (m_corner[0] + m_corner[1] )/2.0;
			
	
			for( int i = 0; i < 2; i ++ )
			for( int j = 0; j < 2; j ++ )
			{
				CPoint2 p = m_corner[0] + CPoint2( i, j ) * len;
				CPoint2 q = p + CPoint2(len, len );
				
				if( pt[0] >= p[0] && pt[1] >= p[1] && pt[0] <= q[0] && pt[1] <= q[1] )
				{
					return m_child[i][j]->_locate( pt );
				}
			}

			return NULL;

		};

	public:

		CPoint2 m_corner[2];
		CQuadTreeNode * m_child[2][2];
		std::vector<T*> m_samples;	
	};

	/*!	CQuadTree class
	 *
	 */
	template <typename T>
	class CQuadTree
	{
	public:
		CQuadTree<T>();
		~CQuadTree<T>();
		CQuadTreeNode<T> * root() { return m_root; };

		void _construct( std::vector<T*> & samples, int n );
		T* CQuadTree<T>::_locate( CPoint2 & pt );

	protected:

		CQuadTreeNode<T> * m_root;

	};


	template <typename T>
	CQuadTree<T>::CQuadTree()
	{
		m_root = new CQuadTreeNode<T>( CPoint2(-1,-1), CPoint2(1,1) );
		//m_root = new CQuadTreeNode<T>( CPoint2(0.8,0.8), CPoint2(0.85,0.85) ); for debugging
		if( m_root == NULL )
		{
			std::cerr<< "Error in allocation" << std::endl;
		}
	};

	template <typename T>
	CQuadTree<T>::~CQuadTree()
	{
		delete m_root;
	};
	
	template <typename T>
	void CQuadTree<T>::_construct( std::vector<T*> & samples, int n )
	{
		m_root->subdivide( samples, n );
	};

	template <typename T>
	T* CQuadTree<T>::_locate( CPoint2 & pt )
	{
		return m_root->_locate( pt );
	};

};

#endif