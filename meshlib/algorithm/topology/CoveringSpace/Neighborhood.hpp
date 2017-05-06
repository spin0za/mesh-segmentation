template<typename V, typename E, typename F, typename H, typename M>
CoveringSpace<V,E,F,H,M>::CNeighborhood::CNeighborhood( V* pV, M* pM )
{
	m_pMesh = pM;

	std::map<V*,CVertex*> vmap;

	m_center = new CVertex;
	m_verts.push_back( m_center );
	assert( m_center );
	m_center->base() = pV;
	vmap[pV] = m_center;
	
	for( M::VertexVertexIterator vviter(pV ); !vviter.end(); vviter++ )
	{
		V * pBV = *vviter;
		CVertex * pW = new CVertex;
		assert( pW );
		pW->base() = pBV;
		m_verts.push_back( pW );
		vmap[pBV] = pW;
	}

	for( M::VertexFaceIterator vfiter(pV ); !vfiter.end(); vfiter++ )
	{
		F * pBF = *vfiter;
		
		std::vector<CVertex*> verts;

		for( M::FaceVertexIterator fviter( pBF ); !fviter.end(); fviter ++ )
		{
			V * pBV = *fviter;
			CVertex * pCV = vmap[pBV];
			verts.push_back( pCV );
		}

		CFace * pF = createFace( verts );
		assert( pF );
		pF->base() = pBF;
		
		//compute the base for each halfedge
		std::vector<H*> hes;

		for( M::FaceHalfedgeIterator hviter( pBF ); !hviter.end(); hviter ++ )
		{
			H * pH = *hviter;
			hes.push_back( pH );
		}

		for( size_t i = 0; i < 3; i ++ )
		{
			H * pH = hes[i];
			
			CHalfEdge * ph = pF->halfedge();
			for( int j = 0; j < 3; j ++ )
			{
				if( ph->target()->base() == m_pMesh->halfedgeTarget( pH ) )
				{
					ph->base() = pH;
					break;
				}

				ph = ph->next();
			}
		}
	}

	//glue the edges

	std::priority_queue<CHalfEdge*, std::vector<CHalfEdge*>, compareHalfEdge>  hqueue;

	for( std::list<CHalfEdge*>::iterator hiter = m_halfedges.begin(); hiter != m_halfedges.end(); hiter ++ )
	{
		CHalfEdge * pH = *hiter;
		hqueue.push( pH );
	}

	while( !hqueue.empty() )
	{
		CHalfEdge* pH = hqueue.top();
		hqueue.pop();

		if( hqueue.empty() ) break;

		CVertex * pS = pH->source();
		CVertex * pT = pH->target();

		std::cout << pS->base()->id() << " " << pT->base()->id() << std::endl;

		CHalfEdge* pD = hqueue.top();
		std::cout << pD->source()->base()->id() << " " << pD->target()->base()->id() << std::endl << std::endl;
		
		if( !compareHalfEdge()( pH, pD ) && !compareHalfEdge()( pD, pH ) ) 
		{
			CEdge * pE = new CEdge;
			assert( pE );
			m_edges.push_back( pE );

			pE->halfedge(0) = pH;
			pE->halfedge(1) = pD;

			pH->edge() = pE;
			pD->edge() = pE;

			H * pBH = pH->base();
			E * pBE = m_pMesh->halfedgeEdge( pBH );
			pE->base() = pBE;
		}
	}
};

template<typename V, typename E, typename F, typename H, typename M>
void CoveringSpace<V,E,F,H,M>::CNeighborhood::_embed()
{
	//copy metric

	for( std::list<CEdge*>::iterator eiter = m_edges.begin(); eiter != m_edges.end(); eiter++ )
	{
		CEdge * pE = *eiter;
		pE->length() = pE->base()->length();
	}

	//compute corner angle
	for( size_t i = 0; i < m_faces.size(); i ++ )
	{
		CFace * pF = m_faces[i];

		CHalfEdge * hes[3];

		CHalfEdge * pH = pF->halfedge();
		for( int i = 0; i < 3; i ++ )
		{
			hes[i] = pH;
			pH = pH->next();
		}

		double l[3];
		for( int i = 0; i < 3; i ++ )
		{
			H * pBH = hes[i]->base();
			E * pBE = m_pMesh->halfedgeEdge( pBH );
			l[i] = pBE->length();
		}

		for( int i = 0; i < 3; i ++ )
		{
			hes[(i+1)%3]->angle() = hyperbolic_cosine_law()( l[(i+1)%3], l[(i+2)%3], l[i]); 
		}	
	}

	//for half edges ccwly

	CFace * pF = m_faces[0];

	CHalfEdge * pHead = pF->halfedge();

	for( int i = 0; i < 3; i ++ )
	{
		if( pHead->target() == m_center ) break;
		pHead = pHead->next();
	}

	//find surrounding in halfedges CCWly

	std::vector<CHalfEdge*> halfedges; 
	
	CHalfEdge * pH = pHead;
	do{
		halfedges.push_back( pH );
		pH = pH->dual()->prev();
	}while( pH != pHead );

	m_center->uv() = std::complex<double>(0,0);

	double theta = 0;
	for( size_t i = 0; i < halfedges.size(); i ++ )
	{
		CHalfEdge * pH = halfedges[i];
		theta += pH->angle();

		CVertex   * pV = pH->source();
		CEdge     * pE = pH->edge();

		double   r = pE->length();
		r = (exp(r) - 1.0)/(exp(r)+1.0);
		pV->uv() = std::complex<double>( r * cos( theta), r * sin( theta ) );			
	}

};


template<typename V, typename E, typename F, typename H, typename M>
CoveringSpace<V,E,F,H,M>::CNeighborhood::~CNeighborhood()
{
	//delete m_center;
	for( size_t i = 0; i < m_verts.size(); i ++ )
	{
		CVertex * pV = m_verts[i];
		delete pV;
	}
	m_verts.clear();

	for( size_t i = 0; i < m_faces.size(); i ++ )
	{
		CFace * pF = m_faces[i];
		delete pF;
	}
	m_faces.clear();

	for( std::list<CHalfEdge*>::iterator hiter = m_halfedges.begin(); hiter != m_halfedges.end(); hiter  ++ )
	{
		CHalfEdge * pH = *hiter;
		delete pH;
	}
	m_halfedges.clear();

	for( std::list<CEdge*>::iterator eiter = m_edges.begin(); eiter != m_edges.end(); eiter  ++ )
	{
		CEdge * pE = *eiter;
		delete pE;
	}
	m_edges.clear();
};




template<typename V, typename E, typename F, typename H, typename M>
void CoveringSpace<V,E,F,H,M>::CNeighborhood::operator*=( CMobius & mob )
{
	for( size_t i = 0; i < m_verts.size(); i ++ )
	{
		CVertex * pV = m_verts[i];
		pV->uv() = mob * pV->uv();
		std::cout << pV->base()->id() << ":" << pV->uv() << std::endl;
	}
};


