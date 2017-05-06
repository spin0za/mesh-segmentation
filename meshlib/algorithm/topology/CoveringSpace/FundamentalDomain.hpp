//Constructor for fundamental domain
template<typename V, typename E, typename F, typename H, typename M>
CoveringSpace<V,E,F,H,M>::CFundamentalDomain::CFundamentalDomain( M* pM, M* pD )
{
	m_pMesh   = pM;
	m_pDomain = pD;
	
	_locate_base();
	_construct_face_map();
	_construct_edge_map();
};

//construct projection map from the closed mesh to the fundamental domain mesh
template<typename V, typename E, typename F, typename H, typename M>
void CoveringSpace<V,E,F,H,M>::CFundamentalDomain::_locate_base()
{

	m_origin = NULL;
	double min_d = 1e+10;

	for( M::MeshVertexIterator viter( m_pDomain ); !viter.end(); viter ++ )
	{
		 V * pV = *viter;
		 if( pV->uv().norm() < min_d )
		 {
			 m_origin = pV;
			 min_d = pV->uv().norm();
		 }
	}

	m_base = m_pMesh->idVertex( m_origin->father() );
}

//construct projection map from the closed mesh to the fundamental domain mesh
template<typename V, typename E, typename F, typename H, typename M>
void CoveringSpace<V,E,F,H,M>::CFundamentalDomain::_construct_face_map()
{
	//find covering face for m_pDomain, m_pMesh

	for( M::MeshFaceIterator fiter( m_pDomain ); !fiter.end(); fiter ++ )
	{
		F * pF = *fiter;
		std::vector<V*> pVs;
		for( M::FaceVertexIterator fviter( pF ); !fviter.end(); fviter ++ )
		{
			V * pV = *fviter;
			pVs.push_back( pV );
		}
		V * pW0 = m_pMesh->idVertex( pVs[0]->father() );
		V * pW1 = m_pMesh->idVertex( pVs[1]->father() );
		E * pWE = m_pMesh->vertexEdge( pW0, pW1 );
		H * pWH = m_pMesh->edgeHalfedge( pWE, 0 );
		H * pWN = m_pMesh->halfedgeNext( pWH );
		V * pW2 = m_pMesh->halfedgeTarget( pWN );
		F * pWF = m_pMesh->halfedgeFace( pWH );

		if( pW2->id() == pVs[2]->father() )
		{
			m_face_map[pWF] = pF;
			m_face_map[pF]  = pWF;
			continue;
		}

		pWH = m_pMesh->edgeHalfedge( pWE, 1 );
		pWN = m_pMesh->halfedgeNext( pWH );
		pW2 = m_pMesh->halfedgeTarget( pWN );
		pWF = m_pMesh->halfedgeFace( pWH );
		assert( pW2->id() == pVs[2]->father() );
		m_face_map[pWF] = pF;
		m_face_map[pF]  = pWF;
	}
};

//construct dual edge map

template<typename V, typename E, typename F, typename H, typename M>
void CoveringSpace<V,E,F,H,M>::CFundamentalDomain::_construct_edge_map()
{
	//find dual edge
	std::priority_queue< std::pair<E*,M*>, std::vector<std::pair<E*,M*>>, compareEdgeMeshPair>  queue;

	for( M::MeshEdgeIterator eiter( m_pDomain ); !eiter.end(); eiter ++ )
	{
		 E * pE = *eiter;
		 if( !pE->boundary() ) continue;
		 queue.push( std::pair<E*,M*>( pE,m_pDomain) );
	}

	while( ! queue.empty() )
	{
		std::pair<E*,M*> pE = queue.top();
		queue.pop();
		std::pair<E*,M*> pW = queue.top();
		if( !compareEdgeMeshPair()( pE, pW) && !compareEdgeMeshPair()( pW, pE) )
		{
			queue.pop();
			m_edge_map[pE.first] = pW.first;
			m_edge_map[pW.first] = pE.first;
		}
	}
};

