/*! \file ConvexHull.h
*   \brief ConvexHull
*   \author David Gu
*   \date   documented on 02/12/2011
*
*   Mesh for viewer 
*/
#ifndef  _CCONVEX_HULL_H_
#define  _CCONVEX_HULL_H_

#include <map>
#include <vector>
#include <cmath>
#include <iostream>
#include <fstream>
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
#include "Geometry/MemoryPool.h"

namespace MeshLib
{
	namespace ConvexHull
	{
	/****************************************************************************************
		Matrix Operations
	****************************************************************************************/
		class CMatrix
		{
		public:
			
			CMatrix(){
			for( int i = 0; i < 4; i ++ )
				for( int j = 0; j < 4; j ++ ) 
					m_v[i][j] = 0;
			
			};
			~CMatrix(){};
			double & operator()(int i, int j) { return m_v[i][j]; }
			double determinant() { 
				  double value =
				  m_v[0][3] * m_v[1][2] * m_v[2][1] * m_v[3][0]-m_v[0][2] * m_v[1][3] * m_v[2][1] * m_v[3][0]-m_v[0][3] * m_v[1][1] * m_v[2][2] * m_v[3][0]+m_v[0][1] * m_v[1][3] * m_v[2][2] * m_v[3][0]+
				  m_v[0][2] * m_v[1][1] * m_v[2][3] * m_v[3][0]-m_v[0][1] * m_v[1][2] * m_v[2][3] * m_v[3][0]-m_v[0][3] * m_v[1][2] * m_v[2][0] * m_v[3][1]+m_v[0][2] * m_v[1][3] * m_v[2][0] * m_v[3][1]+
				  m_v[0][3] * m_v[1][0] * m_v[2][2] * m_v[3][1]-m_v[0][0] * m_v[1][3] * m_v[2][2] * m_v[3][1]-m_v[0][2] * m_v[1][0] * m_v[2][3] * m_v[3][1]+m_v[0][0] * m_v[1][2] * m_v[2][3] * m_v[3][1]+
				  m_v[0][3] * m_v[1][1] * m_v[2][0] * m_v[3][2]-m_v[0][1] * m_v[1][3] * m_v[2][0] * m_v[3][2]-m_v[0][3] * m_v[1][0] * m_v[2][1] * m_v[3][2]+m_v[0][0] * m_v[1][3] * m_v[2][1] * m_v[3][2]+
				  m_v[0][1] * m_v[1][0] * m_v[2][3] * m_v[3][2]-m_v[0][0] * m_v[1][1] * m_v[2][3] * m_v[3][2]-m_v[0][2] * m_v[1][1] * m_v[2][0] * m_v[3][3]+m_v[0][1] * m_v[1][2] * m_v[2][0] * m_v[3][3]+
				  m_v[0][2] * m_v[1][0] * m_v[2][1] * m_v[3][3]-m_v[0][0] * m_v[1][2] * m_v[2][1] * m_v[3][3]-m_v[0][1] * m_v[1][0] * m_v[2][2] * m_v[3][3]+m_v[0][0] * m_v[1][1] * m_v[2][2] * m_v[3][3];
			   
				  return value;
			}

		protected:
			double   m_v[4][4];

		};

		std::ostream & operator<< ( std::ostream & is, CMatrix & m );

		/******************************************************************************************************
				
			CPoint4 in R4 class

		*******************************************************************************************************/
		class CVertex;
		class CEdge;
		class CFace;
		class CTet;
		class CSimplex4;

		class CPoint4
		{
		public:
			CPoint4(){};
			CPoint4( const CPoint4 & pt ) { for( int i = 0; i < 4; i ++ ) m_v[i] = pt.m_v[i]; };
			CPoint4( double x, double y, double z, double w) { m_v[0] = x; m_v[1] = y; m_v[2] = z; m_v[3] = w; };
			~CPoint4(){};
			
			double & operator[]( int i ) { return m_v[i]; };
			
			CPoint4 operator-( CPoint4 & p )
			{
				return CPoint4( m_v[0]-p[0], m_v[1]-p[1], m_v[2] - p[2], m_v[3]-p[3] );
			};

		protected:
			double m_v[4];
		};


		/******************************************************************************************************
				
			Simplex classes

		*******************************************************************************************************/

		class CVertex
		{
		public:
			CVertex( int id ) { m_id = id; };
			CVertex( CPoint4 p ) { m_point = p; };
			~CVertex() {};
			CPoint4 & point() { return m_point; };
			int & id() { return m_id; };
			std::vector<CFace*> & faces() { return m_faces; };
			
			CFace * find( CFace * pF );
			void remove( CFace * pF  );
			

		protected:
			CPoint4 m_point;
			int     m_id;
			std::vector<CFace*> m_faces;
		};

		class CEdge
		{
		public:
			CEdge( CVertex * v0, CVertex * v1 ) { m_v[0] = v0; m_v[1] = v1; };
			~CEdge(){};
			CVertex * & operator[](int i ) { return m_v[i]; };
		protected:
			CVertex * m_v[2];
		};

		class CFace
		{
		public:
			CFace() { m_tets[0] = NULL; m_tets[1] = NULL; };

			CFace( CVertex * v0, CVertex * v1, CVertex * v2 ) { m_v[0] = v0; m_v[1] = v1; m_v[2] = v2; m_tets[0] = NULL; m_tets[1]=NULL; };
			~CFace(){};
			CVertex *& operator[]( int i ) { return m_v[i]; };
			CVertex *  operator[](int i ) const { return m_v[i]; };

			//verify if two faces are equal
			bool operator== ( const CFace & face ) const
			{
				std::vector<int> ids;
				std::vector<int> fid;
				for( int i = 0; i < 3; i ++ ) { ids.push_back( m_v[i]->id() ); fid.push_back( face[i]->id() ); };

				std::sort( ids.begin(), ids.end() );
				std::sort( fid.begin(), fid.end() );

				for( int i = 0; i < 3; i ++ )
				{
					if( ids[i] != fid[i] ) return false;
				}

				return true;
			};
			
			CTet* & tets(int i) { return m_tets[i]; };
			CVertex * min_vert() 
			{
				for( int i = 0; i < 3; i ++ )
				{
					if( m_v[i]->id() <= m_v[(i+1)%3]->id() && m_v[i]->id() <= m_v[(i+2)%3]->id() ) return m_v[i];
				}
				return NULL;
			};

			CFace * & next() { return m_next; };
			CFace * & prev() { return m_prev; };

		protected:
			CVertex * m_v[3];
			//assumption:
			//[m_tets[0], this] = +1
			//[m_tets[1], this] = -1
			CTet * m_tets[2];
		protected:
			CFace * m_prev, * m_next;
		};

		class CTet
		{
		public:
			CTet() { m_visible = false; };
			CTet( CVertex * v0, CVertex * v1, CVertex * v2, CVertex * v3 ) { m_v[0] = v0; m_v[1] = v1; m_v[2] = v2; m_v[3] = v3; m_visible = false; };
			~CTet(){};
			CVertex *& operator[]( int i ) { return m_v[i]; };
			CVertex *  operator[]( int i ) const { return m_v[i]; };

			//verify if two faces are equal
			bool operator== ( const CTet & tet ) const
			{
				std::vector<int> ids;
				std::vector<int> fid;
				for( int i = 0; i < 4; i ++ ) { ids.push_back( m_v[i]->id() ); fid.push_back( tet[i]->id() ); };

				std::sort( ids.begin(), ids.end() );
				std::sort( fid.begin(), fid.end() );

				for( int i = 0; i < 4; i ++ )
				{
					if( ids[i] != fid[i] ) return false;
				}

				return true;
			};
			
			bool & visible() { return m_visible; };
			CTet * & next() { return m_next; };
			CTet * & prev() { return m_prev; };

		protected:
			CVertex * m_v[4];
			bool      m_visible;
		protected:
			CTet * m_prev, * m_next;
		};

		class CSimplex4
		{
		public:
			CSimplex4() {};

			CSimplex4( CVertex * v0, CVertex * v1, CVertex * v2, CVertex * v3, CVertex * v4 ) { m_v[0] = v0; m_v[1] = v1; m_v[2] = v2; m_v[3] = v3; m_v[4] = v4; };
			~CSimplex4(){};
			CVertex *& operator[]( int i ) { return m_v[i]; };
			CVertex *  operator[]( int i ) const { return m_v[i]; };
			double volume();
		protected:
			CVertex * m_v[5];
		};



		/*******************************************************************************************************
		
		Convex Hull Algorithm

		********************************************************************************************************/



		class CConvexHull
		{
		public:
			CConvexHull();
			~CConvexHull();
			void compute_convex_hull( std::vector<CPoint4> & pts );
			std::set<CFace*> & faces() { return m_faces; };
			void _output( const char * name );
			
		protected:
			
			void _initialize();
			void _insert_one_vertex( CVertex * pV );
			void _visibility( CVertex * pV );
			void _remove_upper_tets();
		
		protected:
			//boundary operator
			void boundary( CSimplex4 & Simplex, std::set<CTet*> & tets );
			//boundary operator
			void boundary( CTet & Simplex, std::set<CFace*> & faces );

			int contangency( CTet * pTet, CFace * pF );
			
			void consistency_check();
		
		protected:
			std::set<CTet*>    m_tets;
			std::set<CFace*>   m_faces;
			std::vector<CVertex*> m_verts;

			CMPool<CFace>    m_face_pool;
			CMPool<CTet>     m_tet_pool;
		};


	};

}

#endif