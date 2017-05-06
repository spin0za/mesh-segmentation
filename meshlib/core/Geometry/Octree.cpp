#include "Octree.h"

using namespace MeshLib;

COctree::COctree()
{
	m_root = new COctreeNode( CPoint(-1,-1,-1), CPoint(1,1,1) );
}

COctree::~COctree()
{
	delete m_root;

	for( size_t i = 0; i < m_trs.size(); i ++ )
	{
		delete m_trs[i];
	}
}

void COctree::_construct( std::vector<CPoint> & pts, int n )
{
	m_root->subdivide( pts, n );
}

void COctree::_construct( int n )
{
	m_root->subdivide( m_trs, n );
}

void COctree::_insert_triangle( CPoint a, CPoint b, CPoint c )
{
	CTriangle * pT = new CTriangle;
	pT->v[0] = a;
	pT->v[1] = b;
	pT->v[2] = c;

	m_trs.push_back( pT );
}

void COctree::_insert_triangle( const CTriangle & tri )
{
	CTriangle * pT = new CTriangle( tri );
	m_trs.push_back( pT );
}


void COctreeNode::_sample( std::vector<CPoint> & samples )
{
	if( m_triangles.empty() ) return;

	//leaf node
	if( m_child[0][0][0] == NULL )
	{
		CPoint rd;	
		for (int i =0; i<3;++i)
		{
			rd[i] =((double) rand() / (RAND_MAX+1)) * 2 - 1.0;
		}
		
		double len = m_corner[1][0] - m_corner[0][0];
		//ditter
		CPoint p = (m_corner[0] + m_corner[1])/2.0 + rd * len/4.0;

		for( size_t k = 0; k < m_triangles.size(); k ++ )
		{
			CTriangle * pT = m_triangles[k];
			CPoint q;
			if( pT->_project( p, q ) ) 
			{
				samples.push_back( q );
				return;
			}
		}

	}
	else
	{
		for( int i = 0; i < 2; i ++ )
		for( int j = 0; j < 2; j ++ )
		for( int k = 0; k < 2; k ++ )
		{
			m_child[i][j][k]->_sample( samples );
		}
	}
}

void COctreeNode::_sample( std::vector<CSample> & samples )
{
	if( m_triangles.empty() ) return;

	//leaf node
	if( m_child[0][0][0] == NULL )
	{
		CPoint rd;	
		for (int i =0; i<3;++i)
		{
			rd[i] =((double) rand() / (RAND_MAX+1)) * 2 - 1.0;
		}
		
		double len = m_corner[1][0] - m_corner[0][0];
		//ditter
		CPoint p = (m_corner[0] + m_corner[1])/2.0 + rd * len/4.0;

		for( size_t k = 0; k < m_triangles.size(); k ++ )
		{
			CTriangle * pT = m_triangles[k];
			CSample  q;
			if( pT->_project( p, q ) ) 
			{
				samples.push_back( q );
				return;
			}
		}

	}
	else
	{
		for( int i = 0; i < 2; i ++ )
		for( int j = 0; j < 2; j ++ )
		for( int k = 0; k < 2; k ++ )
		{
			m_child[i][j][k]->_sample( samples );
		}
	}
}

void COctree::_sample( std::vector<CPoint> & samples )
{
	srand((unsigned)time(NULL));
	if( m_root != NULL )
		m_root->_sample( samples );
};

void COctree::_sample( std::vector<CSample> & samples )
{
	srand((unsigned)time(NULL));
	if( m_root != NULL )
		m_root->_sample( samples );
};


bool CTriangle::_project( CPoint &in, CPoint &out )
{
	CPoint n = (v[1] - v[0])^(v[2]-v[0]);
	double d = n.norm();
	if( d < 1e-8 ) return false;
	n = n/d;

	double t = (v[0] - in ) * n;

	out = in + n * t;
	
	CPoint bary;

	for( int i = 0; i < 3; i ++ )
	{
		CPoint l = v[(i+0)%3] - out;
		CPoint r = v[(i+1)%3] - out;
		bary[i] = (l^r) * n;
	}

	return (bary[0]>=0 && bary[1] >= 0 && bary[2] >= 0 );
};


bool CTriangle::_project( CPoint &in, CSample & sample )
{
	CPoint n = (v[1] - v[0])^(v[2]-v[0]);
	double d = n.norm();
	if( d < 1e-8 ) return false;

	n = n/d;

	double t = (v[0] - in ) * n;

	CPoint out = in + n * t;
	
	CPoint bary;

	for( int i = 0; i < 3; i ++ )
	{
		CPoint l = v[(i+1)%3] - out;
		CPoint r = v[(i+2)%3] - out;
		bary[i] = (l^r) * n/d;
	}
	
	sample.point() = out;
	sample.uv()    = uv[0] * bary[0] + uv[1] * bary[1] + uv[2] * bary[2];
	sample.rgb()   = rgb[0] * bary[0] + rgb[1] * bary[1] + rgb[2] * bary[2];

	return (bary[0]>=0 && bary[1] >= 0 && bary[2] >= 0 );
};