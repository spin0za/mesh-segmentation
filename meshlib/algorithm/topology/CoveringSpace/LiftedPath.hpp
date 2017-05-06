template<typename V, typename E, typename F, typename H, typename M>
CoveringSpace<V,E,F,H,M>::CLiftedPath::CLiftedPath( M* pM, M* pD, const char * input )
{
	m_pMesh   = pM;
	m_pDomain = pD;

	std::ifstream myfile(input);

	std::vector<H*> hedges;
	
	if (myfile.is_open())
	{
		while ( myfile.good() )
		{
			int id1, id2;
			myfile >> id1;
			myfile >> id2;

			V * pV1 = pM->idVertex( id1 );
			V * pV2 = pM->idVertex( id2 );

			E * pE = pM->vertexEdge( pV1, pV2 );
			H * pH = pM->edgeHalfedge( pE, 0 );
			
			if( pM->halfedgeTarget( pH ) == pV2 )
				hedges.push_back( pH );
			else
			{
				pH = pM->edgeHalfedge(pE,1);
				hedges.push_back( pH );
			}
		}
	}

	myfile.close();	

	H  * pH = hedges[0];
	V  * pV = pM->halfedgeSource( pH );

	CNeighborhood * pN  = new CNeighborhood( pV, pM );
	assert( pN );
	
	m_neighborhoods.push_back( pN );

	for( size_t i = 0; i < hedges.size(); i ++ )
	{
		H * pH = hedges[i];
		V   * pV = pM->halfedgeTarget( pH );
		CNeighborhood * pN  = new CNeighborhood( pV, pM );
		assert( pN );

		m_neighborhoods.push_back( pN );
	}
};



template<typename V, typename E, typename F, typename H, typename M>
void CoveringSpace<V,E,F,H,M>::CLiftedPath::_overlap( CNeighborhood * pN1, CNeighborhood * pN2, std::vector<std::pair<CFace*,CFace*> > & pairs )
{
	std::priority_queue<std::pair<CFace*,CNeighborhood*>, std::vector<std::pair<CFace*,CNeighborhood*>>, compareFaceInNeighborhood>  fqueue;

	for( size_t i = 0; i < pN1->faces().size(); i ++ )
	{
		CFace * pF = pN1->faces()[i];
		fqueue.push( std::pair<CFace*,CNeighborhood*>(pF,pN1) );
	}

	for( size_t i = 0; i < pN2->faces().size(); i ++ )
	{
		CFace * pF = pN2->faces()[i];
		fqueue.push( std::pair<CFace*,CNeighborhood*>(pF,pN2) );
	}
	
	while( !fqueue.empty() )
	{
		std::pair<CFace*,CNeighborhood*> pF = fqueue.top();
		fqueue.pop();

		if( fqueue.empty() ) continue;

		std::pair<CFace*,CNeighborhood*>  pD = fqueue.top();

		if( pF.first->base() != pD.first->base() ) continue;

		fqueue.pop();
		
		std::pair<CFace*,CFace*> pair_face;

		std::cout<< "Overlap " << pF.second->center()->base()->id()<< " " << pF.first->base()->id() << " ";
		std::cout<< pD.second->center()->base()->id()<< " " << pD.first->base()->id() << std::endl;

		pair_face.first  = (pF.second == pN1)? pF.first: pD.first;
		pair_face.second = (pF.second == pN2)? pF.first: pD.first;

		pairs.push_back( pair_face );
	}
};





template<typename V, typename E, typename F, typename H, typename M>
CMobius CoveringSpace<V,E,F,H,M>::CLiftedPath::MT( std::complex<double> za, std::complex<double> zb )
{
	CMobius  m(za,0); 
	zb = m * zb;
	Complex theta= std::conj( zb );

	return CMobius( za, std::arg( theta) );

}

template<typename V, typename E, typename F, typename H, typename M>
CMobius CoveringSpace<V,E,F,H,M>::CLiftedPath::_transformation( std::complex<double> A, std::complex<double> B, std::complex<double> IA, std::complex<double> IB )
{
	CMobius m0  = MT(A,B);

	std::cout << m0 * B << std::endl;

	CMobius m1  = MT(IA,IB);

	std::cout << m1 * IB << std::endl;

	CMobius im1 = inverse( m1 );
	CMobius  m  = im1 * m0;

	return m;
}

template<typename V, typename E, typename F, typename H, typename M>
CMobius CoveringSpace<V,E,F,H,M>::CLiftedPath::_transform( CFace * pF0, CFace * pF1 )
{
	CHalfEdge * pSH = pF1->halfedge()->next();
	CHalfEdge * pTH = pF0->halfedge()->next();

	H * pBSH = pSH->base();
	H * pBTH = pTH->base();

	E * pBSE = m_pMesh->halfedgeEdge( pBSH );
	E * pBTE = m_pMesh->halfedgeEdge( pBTH );
	
	std::cout << "Target: " ;
	std::cout << "Source " << pTH->source()->uv() << "End " << pTH->target()->uv() << std::endl;


	assert( pSH->base() == pTH->base() );

	CMobius mob = _transformation( pSH->source()->uv(), pSH->target()->uv(), pTH->source()->uv(), pTH->target()->uv()); 

	std::cout <<"Matched: ";
	std::cout <<"Source "<< pSH->source()->base()->id()<<":" << mob * pSH->source()->uv();
	std::cout <<"Target "<< pSH->target()->base()->id()<<":" << mob * pSH->target()->uv() << std::endl;

	return mob;
}

template<typename V, typename E, typename F, typename H, typename M>
void CoveringSpace<V,E,F,H,M>::CLiftedPath::_embed(  CFundamentalDomain * pD )
{
	for( size_t i = 0; i < m_neighborhoods.size(); i ++ )
	{
		CNeighborhood * pN = m_neighborhoods[i];
		pN->_embed();
	}

	//embed the first neighborhood
	CNeighborhood * pN = m_neighborhoods[0];
	
	CFace		  * pF = pN->faces()[0];
	F             * pBF = pF->base();
	F             * pOF = pD->dual_face( pBF );

	CHalfEdge	  * pH = pF->halfedge();
	std::vector<CVertex*> pVs;

	for( int i = 0; i < 3; i ++ )
	{
		pVs.push_back( pH->target());
		pH = pH->next();
	}

	std::vector<V*> pOVs;

	for( M::FaceVertexIterator fviter( pOF ); !fviter.end(); fviter ++ )
	{
		V * pV = *fviter;
		pOVs.push_back( pV );
	}
	
	int b = 0;
	for( int i = 0; i < 3; i ++ )
	{
		if( pOVs[i]->father() == pVs[0]->base()->id() )
		{
			b = i;
			break;
		}
	}

	std::vector<V*> pWs;

	for( int i = 0; i < 3 ; i ++ )
	{
		pWs.push_back( pOVs[(i+b)%3] );
	}

	CMobius mob = _transformation( pVs[0]->uv(), pVs[1]->uv(), 
		std::complex<double>(pWs[0]->uv()[0], pWs[0]->uv()[1]), 
		std::complex<double>(pWs[1]->uv()[0], pWs[1]->uv()[1]) );
	
	(*pN) *= mob;

	for( size_t i = 0; i < m_neighborhoods.size() - 1; i ++ )
	{
		CNeighborhood * pN = m_neighborhoods[i];
		CNeighborhood * pW = m_neighborhoods[i+1];

		std::vector<std::pair<CFace*,CFace*> > pairs;

		_overlap( pN, pW, pairs );

		std::pair<CFace*,CFace*> head = pairs[0];

		CMobius mob = _transform( head.first, head.second );

		(*pW) *= mob;		
	}

}



template<typename V, typename E, typename F, typename H, typename M>
void CoveringSpace<V,E,F,H,M>::CLiftedPath::_output( const char * output_file )
{
	std::ofstream os;
	os.open( output_file );

	int id = 1;

	for( size_t i = 0; i < m_neighborhoods.size(); i ++ )
	{
		CNeighborhood * pN = m_neighborhoods[i];

		for( size_t j = 0; j < pN->vertices().size(); j ++ )
		{
			CVertex * pV = pN->vertices()[j];
			pV->id() = id ++;
			os <<"Vertex " << pV->id()<< " " << pV->uv().real() << " " << pV->uv().imag() << " 0 {rgb=(0 1 0) father=(" << pV->base()->id()<<")}" << std::endl;
		}		
	}

	int fid = 1;

	for( size_t i = 0; i < m_neighborhoods.size(); i ++ )
	{
		CNeighborhood * pN = m_neighborhoods[i];
		for( size_t j = 0; j < pN->faces().size(); j ++ )
		{
			CFace * pF = pN->faces()[j];
			CHalfEdge * pH = pF->halfedge();

			os <<"Face "<< fid++ << " ";

			for( int k = 0; k < 3; k ++ )
			{
				os << pH->target()->id() << " ";
				pH = pH->next();
			}
			os << std::endl;
		}		
	}
	
	os.close();

}

template<typename V, typename E, typename F, typename H, typename M>
CMobius CoveringSpace<V,E,F,H,M>::CLiftedPath::_fuchsian_transformation()
{
	size_t n = m_neighborhoods.size();
	CNeighborhood * pS = m_neighborhoods[0];
	CNeighborhood * pT = m_neighborhoods[n-1];

	std::vector<std::pair<CFace*,CFace*> > pairs;

	_overlap( pS, pT, pairs );

	std::pair<CFace*,CFace*> head = pairs[0];

	CMobius mob = _transform( head.first, head.second );

	return mob;
}


template<typename V, typename E, typename F, typename H, typename M>
void CoveringSpace<V,E,F,H,M>::CLiftedPath::_straighten(  CFundamentalDomain * pD, const char * name, const char * plane_name )
{


	CMobius mob = _fuchsian_transformation();

	std::complex<double> T = mob * std::complex<double>(0,0);

	std::vector<std::pair<CPoint,CPoint2> >  Pts;
	_ray( pD, CPoint2(0,0), CPoint2(T.real(), T.imag()), Pts );

	CPoint2 pT = pD->origin()->uv();

	_read_vertex_uv<M, V, E, F, H>( m_pDomain );

	Pts.clear();
	_ray( pD, CPoint2(0,0), pT, Pts );


	std::ofstream os;
	os.open( name );

	for( size_t i = 0; i < Pts.size()-1; i ++ )
	{
		CPoint p = Pts[i+0].first;
		CPoint q = Pts[i+1].first;

		os <<"L 0 0 0 "<<std::endl;
		os <<"d 1 1 0 "<<std::endl;
		os <<"s 1 1 1 "<<std::endl;
		os <<"g 1 0 0 "<<std::endl;

		os <<"v "<< p[0] << " " << p[1] << " " << p[2] << std::endl;
		os <<"v "<< q[0] << " " << q[1] << " " << q[2] << std::endl;
		os <<"E 0 0 0 "<<std::endl << std::endl;

	}
	os.close();


	std::ofstream pos;
	pos.open( plane_name );

	for( size_t i = 0; i < Pts.size()-1; i ++ )
	{
		CPoint2 p = Pts[i+0].second;
		CPoint2 q = Pts[i+1].second;

		pos <<"L 0 0 0 "<<std::endl;
		pos <<"d 1 1 0 "<<std::endl;
		pos <<"s 1 1 1 "<<std::endl;
		pos <<"g 1 0 0 "<<std::endl;

		pos <<"v "<< p[0] << " " << p[1] << " " << 0 << std::endl;
		pos <<"v "<< q[0] << " " << q[1] << " " << 0 << std::endl;
		pos <<"E 0 0 0 "<<std::endl << std::endl;

	}
	pos.close();

}

template<typename V, typename E, typename F, typename H, typename M>
bool CoveringSpace<V,E,F,H,M>::CLiftedPath::_intersect( CPoint2 ptS, CPoint2 ptT, E * pE, CPoint2 & q , double & lambda )
{
		//starting
		double epsilon = 1e-12;

		CLine2D line( ptS, ptT );

		V * pS = m_pDomain->edgeVertex1( pE );
		V * pT = m_pDomain->edgeVertex2( pE );

		double s = line._evaluate( pS->uv() );
		double t = line._evaluate( pT->uv() );

		if( (s < -epsilon && t < -epsilon ) || (t > epsilon && s > epsilon ) ) return false;


		CLine2D S( pS->uv(), pT->uv() );

		MeshLib::_intersect( line, S, q );

		//verify if the intersection is on the line segment between pS and pT

		CPoint2 pt_d = ptT - ptS;
		CPoint2 pt_s =  q  - ptS;
		
		if( pt_s* pt_d < -epsilon ) return false;

		lambda = pt_s.norm()/pt_d.norm();

		if( lambda > 1+epsilon || lambda < -epsilon )  return false;

		return true;
}

template<typename V, typename E, typename F, typename H, typename M>
void CoveringSpace<V,E,F,H,M>::CLiftedPath::_ray(  CFundamentalDomain * pD, CPoint2 ptS, CPoint2 ptT, std::vector< std::pair<CPoint,CPoint2> > & pts )
{
	//input the base point coordinates
	pts.push_back( std::pair<CPoint,CPoint2>( pD->base()->point(), ptS) );


	//starting
	CPoint2 q;
	double lambda;

	CLine2D line( ptS, ptT );

	E * pCurrentEdge     = NULL;
	H * pCurrentHalfedge = NULL;

	for( M::VertexOutHalfedgeIterator vhiter( m_pDomain, pD->origin() ); !vhiter.end(); vhiter ++ )
	{
		H * pH = *vhiter;
		H * pN = m_pDomain->halfedgeNext( pH );
		E * pE =  m_pDomain->halfedgeEdge( pN );

		if( _intersect( ptS, ptT, pE, q, lambda ) )
		{
			pCurrentEdge = pE;
			pCurrentHalfedge = pN;
			break;
		}
	}


	V * pS = m_pDomain->edgeVertex1( pCurrentEdge );
	V * pT = m_pDomain->edgeVertex2( pCurrentEdge );

	V * pFS = m_pMesh->idVertex( pS->father() );
	V * pFT = m_pMesh->idVertex( pT->father() );

	double nu = (q - pS->uv()).norm() / ( pT->uv()-pS->uv() ).norm();
	CPoint Q = pFS->point() * (1-nu) + pFT->point() * nu;

	pts.push_back( std::pair<CPoint,CPoint2>(Q,q) );

	while( true )
	{
	
		pCurrentHalfedge = m_pDomain->halfedgeSym( pCurrentHalfedge );
	
		if( pCurrentHalfedge == NULL )
		{
			//break;
			std::cout << "Fuchsian Transform" << std::endl;

			E * pDE = pD->dual_edge(pCurrentEdge);

			_Fuchsian_Transform( pDE, pCurrentEdge );

			pCurrentEdge = pDE;
			pCurrentHalfedge = m_pDomain->edgeHalfedge( pDE, 0 );
		}

		std::cout << m_pDomain->halfedgeSource( pCurrentHalfedge )->id() << " " << m_pDomain->halfedgeTarget( pCurrentHalfedge )->id() << " ";

		H * pH = pCurrentHalfedge;
		for( int i = 0; i < 2; i ++ )
		{
			pH = m_pDomain->halfedgeNext( pH );
			E * pE = m_pDomain->halfedgeEdge( pH );
			if( _intersect( ptS, ptT, pE, q, lambda ) ) 
			{
				pCurrentEdge = pE;
				pCurrentHalfedge = pH;
				break;
			}
		}
		
		double mu = (q - ptS).norm()/(ptT - ptS).norm();
		std::cout << mu << std::endl;
		
		if( mu >= 1 ) break;

		{
			V * pS = m_pDomain->edgeVertex1( pCurrentEdge );
			V * pT = m_pDomain->edgeVertex2( pCurrentEdge );

			V * pFS = m_pMesh->idVertex( pS->father() );
			V * pFT = m_pMesh->idVertex( pT->father() );

			double nu = (q - pS->uv()).norm() / ( pT->uv()-pS->uv() ).norm();
			CPoint Q = pFS->point() * (1-nu) + pFT->point() * nu;
			pts.push_back( std::pair<CPoint,CPoint2>(Q,q) );
		}

	}
	pts.push_back( std::pair<CPoint,CPoint2>(pD->base()->point(),ptT) );

}


template<typename V, typename E, typename F, typename H, typename M>
void CoveringSpace<V,E,F,H,M>::CLiftedPath::_Fuchsian_Transform(  E* eS, E * eT )
{
	V * s1 = m_pDomain->edgeVertex1( eS );
	V * t1 = m_pDomain->edgeVertex2( eS );

	V * s2 = m_pDomain->edgeVertex1( eT );
	V * t2 = m_pDomain->edgeVertex2( eT );

	if( s1->father() != s2->father() ) __swap<V*>(s1,t1);

	assert( s1->father() == s2->father() );
	assert( t1->father() == t2->father() );


	std::complex<double>  A(s1->uv()[0], s1->uv()[1]);
	std::complex<double>  B(t1->uv()[0], t1->uv()[1]);

	std::complex<double> IA(s2->uv()[0], s2->uv()[1]);
	std::complex<double> IB(t2->uv()[0], t2->uv()[1]);

	CMobius mob = _transformation( A,B, IA, IB);

	for( M::MeshVertexIterator viter( m_pDomain ); !viter.end(); viter ++ )
	{
		V * pV = *viter;
		std::complex<double> z( pV->uv()[0], pV->uv()[1] );
		z = mob * z;
		pV->uv() = CPoint2( z.real(), z.imag() ); 
	}
}


template<typename V, typename E, typename F, typename H, typename M>
CoveringSpace<V,E,F,H,M>::CLiftedPath::~CLiftedPath()
{
	for( size_t i = 0; i < m_neighborhoods.size(); i ++ )
	{
		delete m_neighborhoods[i];
	}
	m_neighborhoods.clear();
};