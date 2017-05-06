/*! \file Dijkstra.h
*   \brief Algorithm for computing the shortest paths on a mesh
*	\author David Gu
*   \date Documented on 03/20/2011
*/

#ifndef _DIJKSTRA_H_
#define _DIJKSTRA_H_

#include <queue>
#include <set>

#include "Mesh/boundary.h"
#include "DijkstraMesh.h"
#include "UndirectedGraph.h"

namespace MeshLib
{
/*! \brief CDijkstra class
*
*	Compute the tight fundamental group generators
*   1. Compute the shortest non-separating, essential loop through a given point
*   2. Cut the surface along the loop
*   3. Compute the shortest non-separating, essential loop, connecting two corresponding points on the boundary.

*/ 
class CDijkstra
{
public:
	/*! CDijkstra constructor
	* \param pMesh input closed mesh 
	* \param pDomain open mesh, pMesh sliced along a loop
	*/
	CDijkstra( CDKMesh * pForm, CDKMesh * pDomain );
	/*! CDijkstra destructor */
	~CDijkstra();

	/*! Shortest Basis */
	void _shortest_basis( CDijkstraVertex * root );

	/*! 3. \beta_k's are a_k's, then compute the conjugate b_k */
	//int shortest_path( const char * loop_file, const char * father_loop_file );
	int shortest_path( const char * loop_file, const char * father_loop_file, const char * ancestor_loop_file );

private:

	/*! Closed mesh */
	CDKMesh * m_pMesh;
	/*! The open mesh  */
	CDKMesh * m_pDomain;
	/*! boundary of the open mesh */
	CDKMesh::CBoundary m_boundary;
	
	/*! Check if the loop is a separation loop */
	bool  _separating( CDKMesh::CLoop & loop );


	/*! trace a loop, from root1, through pE, to root2, root1 and root2 may be the same, may be different */
	double _trace( CDijkstraEdge * pE,  CDKMesh::CLoop & loop );

	
	//find a shortest path from pStart to pEnd, less than the threshold

	double _shortest_path( CDijkstraVertex* pStart, CDijkstraVertex* pEnd, double threshold, CDKMesh::CLoop & loop  );

	/*! compute the distance field */
	void _tree(std::vector<CDijkstraVertex*> & roots );
	void _cotree( CUndirectedGraph<CDijkstraFace*,CDijkstraEdge*, CompareEdge, CompareVertex> &G );

};


//trace from pE->vertex1() to the root, and pE->vertex2() to the root, 
//store the loop in CLoop, root statisfies the condition that root->root() == root
inline double CDijkstra::_trace( CDijkstraEdge * pE, CDKMesh::CLoop & loop )
{
	double dis = 0;

	CHalfEdge * pH = m_pDomain->edgeHalfedge( pE, 0 );

	CDijkstraVertex *pS = m_pDomain->halfedgeSource( pH );
	CDijkstraVertex *pT = m_pDomain->halfedgeTarget( pH );

	loop.halfedges().push_back( pH );
	dis += pE->length();
	
	CDijkstraVertex *pV = pT;
	while( pV != pV->root() )
	{
		CDijkstraVertex * pW = pV->previous();
		CDijkstraEdge   * pE = m_pDomain->vertexEdge( pW,pV );
		dis += pE->length();
		CHalfEdge * pH = m_pDomain->edgeHalfedge( pE, 0 );
		pV = pW;
		loop.halfedges().push_back( pH );
	}
	pV = pS;
	while( pV != pV->root() )
	{
		CDijkstraVertex * pW = pV->previous();
		CDijkstraEdge   * pE = m_pDomain->vertexEdge( pW,pV );
		dis += pE->length();
		CHalfEdge * pH = m_pDomain->edgeHalfedge( pE, 0 );
		pV = pW;
		loop.halfedges().push_front( pH );
	}
	return dis;
};




//Verify if a loop is separating

inline bool  CDijkstra::_separating( CDKMesh::CLoop & loop )
{
	for( CDKMesh::MeshEdgeIterator eiter( m_pMesh ); !eiter.end(); eiter ++ )
	{
		CDijkstraEdge *pE = *eiter;
		pE->sharp() = false;
	}
	for( std::list<CHalfEdge*>::iterator hiter = loop.halfedges().begin(); hiter != loop.halfedges().end(); hiter ++ )
	{
		CHalfEdge * pH = *hiter;
		CDijkstraEdge * pE = m_pDomain->halfedgeEdge( pH );

		CDijkstraVertex * pS = m_pDomain->edgeVertex1( pE );
		CDijkstraVertex * pT = m_pDomain->edgeVertex2( pE );

		CDijkstraVertex * pWS = m_pMesh->idVertex( pS->father() );
		CDijkstraVertex * pWT = m_pMesh->idVertex( pT->father() );
		
		CDijkstraEdge   * pWE = m_pMesh->vertexEdge( pWS, pWT );
		pWE->sharp() = true;
	}

	CDijkstraFace * head = NULL;
	for( CDKMesh::MeshFaceIterator fiter( m_pMesh ); !fiter.end(); fiter ++ )
	{
		CDijkstraFace *pF = *fiter;
		pF->touched() = false;
		head = pF;
	}

	std::queue<CDijkstraFace*> queue;

	head->touched() = true;
	queue.push( head );

	while( !queue.empty() )
	{
		CDijkstraFace * pF = queue.front();
		queue.pop();

		for( CDKMesh::FaceHalfedgeIterator fhiter( pF ); !fhiter.end(); fhiter ++ )
		{
			CHalfEdge * pH = *fhiter;
			CHalfEdge * pS = m_pMesh->halfedgeSym( pH );
			if( pS == NULL ) continue;
			CDijkstraEdge * pE = m_pMesh->halfedgeEdge( pH );
			if( pE->sharp() ) continue;


			CDijkstraFace * pSF = m_pMesh->halfedgeFace( pS );
			if( pSF->touched() ) continue;
			pSF->touched() = true;
			queue.push( pSF );
		}
	}

	for( CDKMesh::MeshFaceIterator fiter( m_pMesh ); !fiter.end(); fiter ++ )
	{
		CDijkstraFace * pF = *fiter;
		if( !pF->touched() ) return true;
	}

	return false;
};




//called by the following method : _shortest_path( CDKMesh::CLoop  )

inline double CDijkstra::_shortest_path( CDijkstraVertex* pStart, CDijkstraVertex* pEnd, double threshold, CDKMesh::CLoop & loop  )
{

	std::priority_queue<CDijkstraVertex*, std::vector<CDijkstraVertex*>, CompareVertex> vqueue;

	for( CDKMesh::MeshVertexIterator viter( m_pDomain ); !viter.end(); viter ++ )
	{
		CDijkstraVertex * v = *viter;
		v->touched()  = false;
		v->distance() = 1e+30;
		v->previous() = NULL;
		v->root() = NULL;
	}

	for( CDKMesh::MeshEdgeIterator eiter( m_pDomain ); !eiter.end(); eiter ++ )
	{
		CDijkstraEdge * e = *eiter;
		if( e->length() < 0 )
			e->length() = m_pDomain->edgeLength( e );
	}

	CDijkstraVertex* root = pStart;
	root->touched() = true;
	root->distance()    = 0;
	root->root() = root;
	vqueue.push( root );
	
	while( !vqueue.empty() )
	{
		CDijkstraVertex * pV = vqueue.top();
		vqueue.pop();

		if( pV->distance() > threshold )
		{
			std::cout << "Give up" << std::endl;
			return 1e+30;
		}

		if( pV == pEnd )
		{
			break;
		}

		for( CDKMesh::VertexVertexIterator vviter( pV ); !vviter.end(); vviter ++ )
		{
			CDijkstraVertex * pW = *vviter;
			if( pW->touched() ) continue;

			CDijkstraEdge * pE = m_pDomain->vertexEdge( pV, pW );

			if( pV->distance() + pE->length() < pW->distance() )
			{
				pW->distance() = pV->distance() + pE->length();
				pW->previous() = pV;
				pW->root() = pV->root();
				pW->touched() = true;
				vqueue.push( pW );
			}
		}	
	}
	
	double dis = 0;
	CDijkstraVertex *pV = pEnd;
	while( pV != pV->root() )
	{
		CDijkstraVertex * pW = pV->previous();
		CDijkstraEdge   * pE = m_pDomain->vertexEdge( pW,pV );
		dis += pE->length();
		CHalfEdge * pH = m_pDomain->edgeHalfedge( pE, 0 );
		pV = pW;
		loop.halfedges().push_back( pH );
	}
	return dis;
};


//compute the distance field

inline void CDijkstra::_tree(std::vector<CDijkstraVertex*> & roots )
{

	std::priority_queue<CDijkstraVertex*, std::vector<CDijkstraVertex*>, CompareVertex> vqueue;

	for( CDKMesh::MeshVertexIterator viter( m_pDomain ); !viter.end(); viter ++ )
	{
		CDijkstraVertex * v = *viter;
		v->touched()  = false;
		v->distance() = 1e+30;
		v->previous() = NULL;
		v->root() = NULL;
	}

	for( CDKMesh::MeshEdgeIterator eiter( m_pDomain ); !eiter.end(); eiter ++ )
	{
		CDijkstraEdge * e = *eiter;
		if( e->length() < 0 )
			e->length() = m_pDomain->edgeLength( e );
	}

	for( unsigned int i = 0; i < roots.size(); i ++ )
	{
		CDijkstraVertex* root = roots[i];
		root->touched() = true;
		root->distance()    = 0;
		root->root() = root;
		vqueue.push( root );
	}


	while( !vqueue.empty() )
	{
		//print( vqueue );

		CDijkstraVertex * pV = vqueue.top();
		vqueue.pop();

		for( CDKMesh::VertexVertexIterator vviter( pV ); !vviter.end(); vviter ++ )
		{
			CDijkstraVertex * pW = *vviter;
			if( pW->touched() ) continue;

			CDijkstraEdge * pE = m_pDomain->vertexEdge( pV, pW );

			if( pV->distance() + pE->length() < pW->distance() )
			{
				pW->distance() = pV->distance() + pE->length();
				pW->previous() = pV;
				pW->root() = pV->root();
				pW->touched() = true;
				vqueue.push( pW );
			}
		}	
	}
};


inline void CDijkstra::_cotree(	CUndirectedGraph<CDijkstraFace*,CDijkstraEdge*, CompareEdge, CompareVertex> & G )
{
	for( CDKMesh::MeshEdgeIterator eiter( m_pDomain ); !eiter.end(); eiter ++ )
	{
		CDijkstraEdge * pE = *eiter;
		pE->sharp() = true;
	}	

	for( CDKMesh::MeshVertexIterator viter( m_pDomain ); !viter.end(); viter ++ )
	{
		CDijkstraVertex * pV = *viter;
		CDijkstraVertex * pW = pV->previous();
		if( pW == NULL ) continue;
		if( pV == pW   ) continue;

		CDijkstraEdge * pE = m_pDomain->vertexEdge( pV, pW );
		pE->sharp() = false;
	}	


	std::set<CDijkstraFace*> fset;
	std::set<CDijkstraEdge*> eset;

	for( CDKMesh::MeshEdgeIterator eiter( m_pDomain ); !eiter.end(); eiter ++ )
	{
		CDijkstraEdge * pE = *eiter;
		if( pE->boundary() ) continue;

		if( pE->sharp() )
		{
			eset.insert( pE );
			CDijkstraFace * pF = m_pDomain->edgeFace1( pE );
			fset.insert( pF );
			pF = m_pDomain->edgeFace2( pE );
			if( pF != NULL )
			{
				fset.insert( pF );
			}
		}
	}	

	for( std::set<CDijkstraFace*>::iterator fiter = fset.begin(); fiter != fset.end(); fiter ++ )
	{
		CDijkstraFace* pF = *fiter;
		G.insert_node( pF );
	}

	for( std::set<CDijkstraEdge*>::iterator eiter = eset.begin(); eiter != eset.end(); eiter ++ )
	{
		CDijkstraEdge * pE = *eiter;

		CDijkstraFace * pF1 = m_pDomain->edgeFace1( pE );
		CDijkstraFace * pF2 = m_pDomain->edgeFace2( pE );

		CNode<CDijkstraFace*,CDijkstraEdge*> * pN1 = G.find( pF1 );
		CNode<CDijkstraFace*,CDijkstraEdge*> * pN2 = G.find( pF2 );

		G.insert_link( pN1, pN2, pE );
	}
};

//Implement Erickson's algorithm, the shortest basis system is written on the mesh m_pDomain as 
//sharp edges

inline void CDijkstra::_shortest_basis( CDijkstraVertex * root )
{
	std::vector<CDijkstraVertex*> roots;
	roots.push_back( root );
	
	_tree( roots );
	CUndirectedGraph<CDijkstraFace*,CDijkstraEdge*, CompareEdge, CompareVertex> G;
	_cotree( G );
	G._prune();
	std::vector<CDijkstraEdge*> seeds;
	G._maximal_spanning_tree( seeds );

	for( CDKMesh::MeshEdgeIterator eiter( m_pDomain ); !eiter.end(); eiter ++ )
	{
		CDijkstraEdge * pE = *eiter;
		pE->sharp() = false;
	}


	for( unsigned int k = 0; k < seeds.size(); k ++ )
	{
		CDijkstraEdge * pE = seeds[k];
		CDKMesh::CLoop loop( m_pDomain );
		double dis = _trace(pE,loop  );

		for( std::list<CHalfEdge*>::iterator hiter = loop.halfedges().begin(); hiter != loop.halfedges().end(); hiter ++)
		{
			CHalfEdge* pH = *hiter;
			CDijkstraEdge* pE = m_pDomain->halfedgeEdge( pH );
			pE->sharp() = true;
		}
	}
/*
	for( unsigned int k = 0; k < G.links().size(); k ++ )
	{
		CDijkstraEdge * pE = G.links()[k]->key();

		CDijkstraVertex * pV1 = (CDijkstraVertex *)(pE->halfedge(0)->target());
		CDijkstraVertex * pV2 = (CDijkstraVertex *)(pE->halfedge(0)->source());
		double d = pV1->distance() + pV2->distance() + pE->length();
		std::cout << " " << d << std::endl;
	}
*/

	//m_pDomain->write_m( "test.m" );
};



/* Used to compute b_i */

inline int CDijkstra::shortest_path( const char * loop_file, const char * father_loop_file, const char * ancestor_loop_file )
{

	std::priority_queue<CDijkstraVertex*, std::vector<CDijkstraVertex*>, MeshLib::CompareVertexFather>  vqueue;
	
	for( unsigned int i = 0; i < m_boundary.loops().size(); i ++ )
	{
		CDKMesh::CLoop * pL = m_boundary.loops()[i];
		
		for( std::list<CHalfEdge*>::iterator hiter = pL->halfedges().begin(); hiter != pL->halfedges().end(); hiter ++ )
		{
			CHalfEdge * pH = *hiter;
			CDijkstraVertex * pT = m_pDomain->halfedgeTarget(pH);
			vqueue.push( pT );
		}
	}

	CDijkstraVertex * pStart = NULL;
	CDijkstraVertex * pEnd   = NULL;
	double min_dis = 1e+30;

	int n = vqueue.size();

	while( !vqueue.empty() )
	{
		CDijkstraVertex * pV = vqueue.top();
		vqueue.pop();
		CDijkstraVertex * pW = vqueue.top();
		vqueue.pop();
		std::cout << n - vqueue.size() << " " << n << std::endl;

		while(  !vqueue.empty() && pV->father() != pW->father() )
		{
			pV = pW;
			pW = vqueue.top();
			vqueue.pop();
		}
		
		if( pV->father() != pW->father() && vqueue.empty() )
		{
			break;
		}

		//crossing vertex with previous loop,skip it

		if( !vqueue.empty() && pW->father() == vqueue.top()->father() )
		{
			printf("Skipping\n");
			while( !vqueue.empty() && vqueue.top()->father() == pW->father() )
			{
				vqueue.pop();
			}
			continue;
		}

		CDKMesh::CLoop loop(m_pDomain); 
		double dis = _shortest_path( pV,pW, min_dis, loop );
		if( dis < min_dis && !_separating( loop ) )
		{
			min_dis = dis;
			pStart = pV;
			pEnd   = pW;
		}
	}


	if( pStart == NULL && pEnd == NULL )
	{
		printf("Erro: There is no non-separating loops\n");
		std::string line;
		std::cin >> line;
		return -1;
	}

	CDKMesh::CLoop loop(m_pDomain); 
	double dis = _shortest_path( pStart,pEnd, min_dis + 1, loop );

	//label shortest essential loop
	for( CDKMesh::MeshEdgeIterator eiter( m_pDomain ); !eiter.end(); eiter ++ )
	{
		CDijkstraEdge * pE = *eiter;
		pE->sharp() = false;
	}

	for( std::list<CHalfEdge*>::iterator hiter = loop.halfedges().begin(); hiter != loop.halfedges().end(); hiter ++ )
	{
		CHalfEdge * pH = *hiter;
		if( pH == NULL ) continue;
		CDijkstraEdge * pE = m_pDomain->halfedgeEdge( pH );
		pE->sharp() = true;
	}
	
	loop.write( loop_file );


	CDKMesh::CLoop father_loop( m_pMesh );
	for( std::list<CHalfEdge*>::iterator hiter = loop.halfedges().begin(); hiter != loop.halfedges().end(); hiter ++ )
	{
		CHalfEdge * pH = *hiter;

		CDijkstraVertex * pS = m_pDomain->halfedgeSource(pH);
		CDijkstraVertex * pT = m_pDomain->halfedgeTarget(pH);

		CDijkstraVertex * pWS = m_pMesh->idVertex( pS->father() );
		CDijkstraVertex * pWT = m_pMesh->idVertex( pT->father() );

		CDijkstraEdge   * pWE = m_pMesh->vertexEdge( pWS, pWT );

		CHalfEdge * pWH = m_pMesh->edgeHalfedge( pWE, 0);

		if( m_pMesh->halfedgeSource( pWH ) != pWS )
		{
			pWH = m_pMesh->edgeHalfedge( pWE, 1);
		}

		father_loop.halfedges().push_back( pWH );

	}

	father_loop.write(father_loop_file);

	//output to the file
	std::ofstream myfile;
	myfile.open (ancestor_loop_file);
	std::vector<CHalfEdge*> hes;
	for( std::list<CHalfEdge*>::iterator hiter = loop.halfedges().begin(); hiter != loop.halfedges().end(); hiter ++ )
	{
		CHalfEdge * pH = *hiter;
		hes.push_back( pH );
	}



	CHalfEdge * pC = hes[0];
	CHalfEdge * pN = hes[1];

	CDijkstraVertex * pCS = m_pMesh->halfedgeSource( pC );
	CDijkstraVertex * pCT = m_pMesh->halfedgeTarget( pC );

	CDijkstraVertex * pNS = m_pMesh->halfedgeSource( pN );
	CDijkstraVertex * pNT = m_pMesh->halfedgeTarget( pN );

	if( pCT == pNS || pCT == pNT )
	{
		myfile << pCS->ancestor() << " " << pCT->ancestor() << std::endl;
	}
	else if( pCS == pNS || pCS == pNT )
	{
		myfile << pCT->ancestor() << " " << pCS->ancestor() << std::endl;
	}
	else assert(0);

	for( size_t i = 1; i < hes.size(); i ++ )
	{
		CHalfEdge * pC = hes[i-1];
		CHalfEdge * pN = hes[i-0];

		CDijkstraVertex * pCS = m_pMesh->halfedgeSource( pC );
		CDijkstraVertex * pCT = m_pMesh->halfedgeTarget( pC );

		CDijkstraVertex * pNS = m_pMesh->halfedgeSource( pN );
		CDijkstraVertex * pNT = m_pMesh->halfedgeTarget( pN );

		if( pCS == pNS || pCT == pNS )
		{
			myfile << pNS->ancestor() << " " << pNT->ancestor() << std::endl;
			continue;
		}

		if( pCS == pNT || pCT == pNT )
		{
			myfile << pNT->ancestor() << " " << pNS->ancestor() << std::endl;
			continue;
		}
		assert(0);
	}



	myfile.close();


	return +1;
};




}

#endif