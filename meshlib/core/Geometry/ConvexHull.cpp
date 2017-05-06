#include "ConvexHull.h"

using namespace MeshLib;


//boundary operator

void ConvexHull::CConvexHull::boundary( ConvexHull::CSimplex4 & Simplex, std::set<ConvexHull::CTet*> & tets )
{
	for( int i = 0; i < 5; i ++ )
	{
		std::vector<ConvexHull::CVertex*> vts;
		
		for( int j = 0; j < 5; j ++ )
		{
			if( j == i ) continue;
			vts.push_back( Simplex[j] );
		}
	
		//swap
		if( i%2 )
		{
			ConvexHull::CVertex * v = vts[0];
			vts[0] = vts[1];
			vts[1] = v;
		}


		ConvexHull::CTet * pT = m_tet_pool.allocate();

		(*pT)[0] = vts[0];
		(*pT)[1] = vts[1];
		(*pT)[2] = vts[2];
		(*pT)[3] = vts[3];


		tets.insert( pT );
	}
};


//boundary operator

void ConvexHull::CConvexHull::boundary( ConvexHull::CTet & Simplex, std::set<ConvexHull::CFace*> & faces )
{
	for( int i = 0; i < 4; i ++ )
	{
		std::vector<ConvexHull::CVertex*> vts;
		
		for( int j = 0; j < 4; j ++ )
		{
			if( j == i ) continue;
			vts.push_back( Simplex[j] );
		}
	
		//swap
		if( i%2 )
		{
			ConvexHull::CVertex * v = vts[0];
			vts[0] = vts[1];
			vts[1] = v;
		}

		ConvexHull::CFace * pF = m_face_pool.allocate();

		(*pF)[0] = vts[0];
		(*pF)[1] = vts[1];
		(*pF)[2] = vts[2];

		faces.insert( pF );
	}
};



bool predicate( const ConvexHull::CFace * & pF, const ConvexHull::CFace * & pW )
{
	return (*pF)==(*pW);
};

int swap( std::vector<int> & ids )
{
	int c = 0;
	for( size_t i = 0; i < ids.size() - 1; i ++ )
	{
		for( size_t j = 0; j < ids.size() - 1 - i; j ++ )
		{
			if( ids[j] > ids[j+1] ) 
			{
				int id = ids[j]; ids[j] = ids[j+1]; ids[j+1]=id; c++;
			}
		}
	}
	return c;
}

//verify if two faces share the same orientation
int  consistent( const ConvexHull::CFace * pF, const ConvexHull::CFace *  pW )
{
	assert( (*pF) == (*pW) );
	
	std::vector<int> ids;
	std::vector<int> fid;
	for( int i = 0; i < 3; i ++ ) { ids.push_back( (*pF)[i]->id() ); fid.push_back( (*pW)[i]->id() ); };
	
	int c1 = swap( ids );
	int c2 = swap( fid );

	if( (c1 + c2) % 2 ) return -1;
	return +1;

};

std::ostream & operator<< ( std::ostream & is, ConvexHull::CMatrix & m )
{
	for( int k = 0; k < 4; k ++ )
	{
		for( int j = 0; j < 4; j ++ )
		{
			is << m(k,j) << " ";
		}
		is << std::endl;
	}
	is << std::endl;
	return is;
};


//testing the visibility of all tets with respect to pV

void ConvexHull::CConvexHull::_visibility( ConvexHull::CVertex * pV )
{

	ConvexHull::CSimplex4 Simplex;
	
	std::set<CTet*> invisible_tets;
	//remove visible tets
	std::vector<ConvexHull::CTet*> visible_tets;

	//visibility testing
	Simplex[0] = pV;
	for( std::set<CTet*>::iterator iter = m_tets.begin(); iter != m_tets.end(); iter ++ )
	{
		ConvexHull::CTet * pT = *iter;

		for( int k = 0; k < 4; k ++ )
		{
			Simplex[k+1] = (*pT)[k];
		}

		if( Simplex.volume() < 0 ) 
		{
			pT->visible() = true;
			visible_tets.push_back( pT );
		}
		else
		{
			pT->visible() = false;
			invisible_tets.insert( pT );
		}
	}

	m_tets.clear();
	m_tets = invisible_tets;

	for( size_t i = 0; i < visible_tets.size(); i ++ )
	{	
		ConvexHull::CTet * pT = visible_tets[i];
		m_tet_pool.delocate( pT );
	}

};

//constructor

ConvexHull::CConvexHull::CConvexHull( )
{
	
};

//construct the initial convex hull

void ConvexHull::CConvexHull::_initialize()
{
	if( m_verts.size() < 5 ) 
	{
		std::cout << "The input points should be more than 5" << std::endl;
		return;
	}

	//initialize 
	ConvexHull::CSimplex4 Simplex;
	
	for( int i = 0; i < 5; i ++ )
	{
		Simplex[i] = m_verts[i];
	}

	//swap
	if( Simplex.volume() < 0 )
	{
		ConvexHull::CVertex* pV  = Simplex[0]; Simplex[0]= Simplex[1]; Simplex[1]   = pV;
	}

	//construct the initial tets
	boundary( Simplex, m_tets );
			
	//construct initial faces
	for( std::set<CTet*>::iterator iter = m_tets.begin(); iter != m_tets.end(); iter ++)
	{
		ConvexHull::CTet * pT = *iter;

		std::set<ConvexHull::CFace*> fs;
		boundary( *pT, fs );

		for( std::set<ConvexHull::CFace*>::iterator iter = fs.begin(); iter != fs.end(); iter++ )
		{
			ConvexHull::CFace * pF = *iter;
			ConvexHull::CVertex * pV = pF->min_vert();
			ConvexHull::CFace * pD = pV->find( pF );
			
			if( pD != NULL )
			{
				//[pT, pF] = -1
				pD->tets(1) =  pT;
				m_face_pool.delocate( pF );
			}
			else
			{
				//[pT, pF] = +1
				pF->tets(0) = pT; 
				pV->faces().push_back( pF );
				m_faces.insert( pF );
			}
		}
	};
};


void ConvexHull::CConvexHull::compute_convex_hull( std::vector<CPoint4> & pts )
{
	for( size_t i = 0; i < pts.size(); i ++ )
	{
		ConvexHull::CVertex * pV = new ConvexHull::CVertex((int)(i+1));
		pV->point() = pts[i];
		m_verts.push_back( pV );
	}

	_initialize();
	consistency_check();
	for( size_t k = 5; k < m_verts.size(); k ++ )
	{
		if( k % 128 == 0 )
		{
			std::cout << k << "/" << m_verts.size() << std::endl;
		}
		ConvexHull::CVertex * pV = m_verts[k];
		_insert_one_vertex( pV );
		consistency_check();
	}
	_remove_upper_tets();
}

int ConvexHull::CConvexHull::contangency( ConvexHull::CTet * pTet, ConvexHull::CFace * pF )
{
	std::set<ConvexHull::CFace*> fs;
	boundary( *pTet, fs );
	for( std::set<ConvexHull::CFace*>::iterator iter = fs.begin(); iter != fs.end(); iter ++ )
	{
		ConvexHull::CFace * pBF = *iter;
		if( (*pBF) == (*pF) ) 
		{
			return consistent( pBF, pF );
		}
	}

	return 0;

}


void ConvexHull::CConvexHull::consistency_check()
{
	return ;

	for( std::set<CFace*>::iterator iter = m_faces.begin(); iter != m_faces.end(); iter ++ )
	{
		ConvexHull::CFace * pF = *iter;
		ConvexHull::CTet  * pT0 = pF->tets(0);
		ConvexHull::CTet  * pT1 = pF->tets(1);

		if( contangency( pT0, pF ) != 1 )
		{
			std::cout << "Error" << std::endl;
		}
		if( contangency( pT1, pF ) != -1 )
		{
			std::cout << "Error" << std::endl;
		}

	}
}

void ConvexHull::CConvexHull::_insert_one_vertex( CVertex * pV )
{

	//test the visibility of each tetrahedron
	_visibility( pV );


	std::vector<ConvexHull::CFace*> visible_faces;
	std::vector<ConvexHull::CFace*> silhoutte_faces;
	std::set<ConvexHull::CFace*>   invisible_faces;

	for( std::set<CFace*>::iterator iter = m_faces.begin(); iter != m_faces.end(); iter ++ )
	{
		CFace * pF = *iter;
		
		CTet  * pT = pF->tets(0);
		CTet  * pW = pF->tets(1);
		
		if(   pT->visible() &&  pW->visible() )  visible_faces.push_back( pF );
		else
		{
			invisible_faces.insert( pF );
		}
		if( ( pT->visible() && !pW->visible() ) || (!pT->visible() &&  pW->visible()) )  silhoutte_faces.push_back( pF );
	}

	m_faces.clear();
	m_faces = invisible_faces;

	//generate new tets
	std::vector<ConvexHull::CTet*> new_tets;

	//Silhoutte face
	for( size_t i = 0; i < silhoutte_faces.size(); i ++ )
	{
		CFace * pF = silhoutte_faces[i];

		if( pF->tets(0)->visible() )
		{
			CTet * pTet = m_tet_pool.allocate();

			(*pTet)[0] = pV;
			(*pTet)[1] = (*pF)[0];
			(*pTet)[2] = (*pF)[1];
			(*pTet)[3] = (*pF)[2];
				
			pF->tets(0) = pTet;
			m_tets.insert( pTet );
			new_tets.push_back( pTet );
		}
		else //pF->tets()[1]->visible()
		{
			CTet * pTet = m_tet_pool.allocate();

			(*pTet)[0] = pV;
			(*pTet)[1] = (*pF)[1];
			(*pTet)[2] = (*pF)[0];
			(*pTet)[3] = (*pF)[2];

			pF->tets(1) = pTet;
			m_tets.insert( pTet );
			new_tets.push_back( pTet );
		}
	}

	//generate new faces

	for( size_t i = 0; i < new_tets.size(); i ++ )
	{
		ConvexHull::CTet * pT = new_tets[i];
		std::set<ConvexHull::CFace*> fs;
		boundary( *pT, fs );

		for( std::set<ConvexHull::CFace*>::iterator iter=fs.begin(); iter != fs.end(); iter ++ )
		{
			ConvexHull::CFace   * pF = *iter;
			ConvexHull::CVertex * pV = pF->min_vert();
			ConvexHull::CFace   * pD = pV->find( pF );

			if( pD != NULL )
			{
				int id = consistent( pD,pF );
				switch( id )
				{
				case +1:
					pD->tets(0) = pT;
					break;
				case -1:
					pD->tets(1) = pT;
					break;
				}
				m_face_pool.delocate( pF );
			}
			else
			{
				pF->tets(0) = pT;
				pV->faces().push_back( pF );
				m_faces.insert( pF );
			}
		}
	}



	//remove visible faces
	for( size_t i = 0; i < visible_faces.size(); i ++ )
	{
		ConvexHull::CFace   * pF = visible_faces[i];
		ConvexHull::CVertex * pV = pF->min_vert();
		pV->remove( pF );
		m_face_pool.delocate( pF );
	}

};


double ConvexHull::CSimplex4::volume()
{
	CMatrix m;

	for( int i = 1; i < 5; i ++ )
	{
		CVertex * pV = m_v[i];
		CPoint4 q = pV->point() - m_v[0]->point();
		for( int j = 0; j < 4; j ++ )
		{
			m(i-1,j) = q[j];
		}
	};

	double vol = m.determinant();

	//std::cout << "Volume is " << vol << std::endl;
	return vol;
};


ConvexHull::CFace * ConvexHull::CVertex::find( ConvexHull::CFace * pF )
{
	for( size_t k = 0; k < m_faces.size(); k ++ )
	{
		ConvexHull::CFace * pW = m_faces[k];
		if( *pF == *pW ) return pW;
	}
	return NULL;
};


void ConvexHull::CVertex::remove( ConvexHull::CFace * pF )
{
	std::vector<ConvexHull::CFace*>::iterator iter = std::find( m_faces.begin(), m_faces.end(), pF );
	if( iter != m_faces.end() ) 
	{
		m_faces.erase( iter );
	}
};


ConvexHull::CConvexHull::~CConvexHull()
{
	m_faces.clear();
	m_tets.clear();

	for( size_t k = 0; k < m_verts.size(); k ++ )
	{
		delete m_verts[k];
	}
	m_verts.clear();

};


void ConvexHull::CConvexHull::_remove_upper_tets()
{
	ConvexHull::CVertex v(-1);
	v.point() = ConvexHull::CPoint4(0,0,0, 1e+20);
	_visibility( &v );

	//remove visible tets
	std::vector<ConvexHull::CTet*> visible_tets;
	for( std::set<ConvexHull::CTet*>::iterator iter = m_tets.begin(); iter != m_tets.end(); iter ++ )
	{	
		ConvexHull::CTet * pT = *iter;
		if( pT->visible() ) visible_tets.push_back( pT );
	}

	std::vector<ConvexHull::CFace*> visible_faces;

	for( std::set<ConvexHull::CFace*>::iterator iter = m_faces.begin(); iter != m_faces.end(); iter ++ )
	{
		CFace * pF = *iter;
		
		CTet  * pT = pF->tets(0);
		CTet  * pW = pF->tets(1);
		
		if(   pT->visible() &&  pW->visible() )  visible_faces.push_back( pF );
	}



	for( size_t i = 0; i < visible_tets.size(); i ++ )
	{
		ConvexHull::CTet * pT = visible_tets[i];
		std::set<ConvexHull::CTet*>::iterator iter = std::find( m_tets.begin(), m_tets.end(), pT );
		assert( iter != m_tets.end() );
		m_tets.erase( iter );
		m_tet_pool.delocate( pT );
	}

	//remove visible faces
	for( size_t i = 0; i < visible_faces.size(); i ++ )
	{
		ConvexHull::CFace * pF = visible_faces[i];
		std::set<ConvexHull::CFace*>::iterator iter = std::find( m_faces.begin(), m_faces.end(), pF );
		m_faces.erase( iter );
		ConvexHull::CVertex * pV = pF->min_vert();
		pV->remove( pF );
		m_face_pool.delocate( pF );
	}
}

void ConvexHull::CConvexHull::_output( const char * name )
{
	std::ofstream of(name);
	for( size_t k = 0; k < m_verts.size(); k ++ )
	{
		ConvexHull::CVertex * pV = m_verts[k];
		of << "Vertex " << pV->id()<< " " << pV->point()[0]<< " " << pV->point()[1]<<" " << pV->point()[2] << std::endl;
	}
	size_t k = 0;
	for( std::set<ConvexHull::CTet*>::iterator iter = m_tets.begin(); iter != m_tets.end(); iter ++ )
	{
		ConvexHull::CTet * pT = *iter;
		of << "Tet " << ++k << " " << (*pT)[0]->id() << " " << (*pT)[1]->id() << " " << (*pT)[2]->id() << " " << (*pT)[3]->id()  << std::endl;
	}
	of.close();

}
