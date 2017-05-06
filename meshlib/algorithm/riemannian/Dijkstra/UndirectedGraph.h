/*! \file UndirectedGraph.h
*   \brief Algorithm for handling Undirected Graph
*	\author David Gu
*   \date Documented on 04/01/2011
*/

#ifndef _UNDIRECTED_GRAPH_H_
#define _UNDIRECTED_GRAPH_H_

#include <queue>
#include <set>
//#include "binaryheap.h"

#include "Mesh/boundary.h"
#include "DijkstraMesh.h"

namespace MeshLib
{

template<typename V, typename E> 
class CLink;


template<typename V,typename E>
class CNode
{
public:
	CNode( V& v ) {m_key  = v; m_touched = false; };
	~CNode(){};
	V & key() { return m_key; };
	void insert_link( CLink<V,E> * pL ) { m_links.push_back( pL ); };
	std::vector<CLink<V,E>*> & links()  { return m_links; };
	bool & touched() { return m_touched; };
	int  & label()   { return m_label;   };

protected:
	V m_key;
	int  m_label;
	bool m_touched;
	std::vector<CLink<V,E> *> m_links;
};

template<typename V,typename E>
class CLink
{
public:
	CLink<V,E>( CNode<V,E> * n1, CNode<V,E> * n2, E& e ) {m_node[0]  = n1; m_node[1] = n2; m_key = e;  m_touched = false;};
	~CLink(){};
	
	E & key() { return m_key; };
	CNode<V,E> * node( int k ) { return m_node[k]; };
	
	bool & touched() { return m_touched; };
	
	CNode<V,E>* other( CNode<V,E>* pN ) { return (m_node[0] != pN )? m_node[0]:m_node[1]; };

protected:
	E m_key;
	CNode<V,E> * m_node[2];
	bool m_touched;
};


template<typename V, typename E, typename ECompare, typename VCompare>
class CUndirectedGraph
{
public:
	CUndirectedGraph(){};
	~CUndirectedGraph();
	
	void insert_node( V & v ) { CNode<V,E> * pN = new CNode<V,E>(v); m_nodes.push_back( pN ); };
	void insert_link( CNode<V,E> * pN1, CNode<V,E>* pN2, E & e ) 
	{ 
		CLink<V,E> * pL = new CLink<V,E>( pN1, pN2, e); 
		m_links.push_back( pL ); 

		pL->node(0)->insert_link( pL );
		pL->node(1)->insert_link( pL );
		
	};
	void remove_link( E & e );

	CNode<V,E>* find( V & key )
	{
		for( unsigned int k = 0; k < m_nodes.size() ; k ++ )
		{
			CNode<V,E> * pN = m_nodes[k];
			if( pN->key() == key )
				return pN;
		}
		return NULL;
	};

	CLink<V,E>* find( E & key )
	{
		for( unsigned int k = 0; k < m_links.size() ; k ++ )
		{
			CLink<V,E> * pL = m_links[k];
			if( pL->key() == key )
				return pL;
		}
		return NULL;
	};

	void _prune();

	void _maximal_spanning_tree ( std::vector<E> & seeds );
	void _maximal_spanning_tree2( std::vector<E> & seeds );

	std::vector<CLink<V,E>*> & links() { return m_links; };
	bool _has_cycle();
protected:

	std::vector<CNode<V,E>*> m_nodes;
	std::vector<CLink<V,E>*> m_links;

	//Breadth first search
	void BFS();


};

template<typename V, typename E, typename ECompare, typename VCompare>
CUndirectedGraph<V,E,ECompare,VCompare>::~CUndirectedGraph()
{
	for( unsigned int k = 0; k < m_nodes.size(); k ++ )
	{
		CNode<V,E> * pN = m_nodes[k];
		delete pN;
	}

	for( unsigned int k = 0; k < m_links.size(); k ++ )
	{
		CLink<V,E> * pL = m_links[k];
		delete pL;
	}
};

template<typename V, typename E, typename ECompare, typename VCompare>
void CUndirectedGraph<V,E,ECompare,VCompare>::BFS()
{
	for( unsigned int k = 0; k < m_nodes.size();  k ++ )
	{
		CNode<V,E> * pN = m_nodes[k];
		pN->touched() = false;
	}

	for( unsigned int k = 0; k < m_links.size();  k ++ )
	{
		CLink<V,E> * pL = m_links[k];
		pL->touched() = false;
	}

	while( true )
	{
		CNode<V,E> * root = NULL;

		for( unsigned int k = 0; k < m_nodes.size();  k ++ )
		{
			CNode<V,E> * pN = m_nodes[k];
			if( pN->touched() ) continue;
			root = pN;
			break;
		}

		if( root == NULL ) break;

		root->touched() = true;

		std::queue<CNode<V,E>*> queue;

		queue.push( root );

		while( !queue.empty() )
		{
			CNode<V,E>* pV = queue.front();
			queue.pop();

			std::vector<CLink<V,E>*> & links = pV->links();

			for( unsigned int k = 0; k < links.size(); k ++ )
			{
				CLink<V,E> * pL = links[k];
				CNode<V,E> * pW = pL->other( pV );
				if( pW->touched() ) continue;
				pW->touched() = true;
				pL->touched() = true;
				queue.push( pW );
			}
		}
	}
};


template<typename V, typename E, typename ECompare, typename VCompare>
void CUndirectedGraph<V,E,ECompare,VCompare>::_prune()
{

	while( true )
	{
		
		std::queue<CNode<V,E>*> leaves;

		for( unsigned int k = 0; k < m_nodes.size();  k ++ )
		{
			CNode<V,E> * pN = m_nodes[k];
			if( pN->links().size() == 1 )
			{
				leaves.push( pN );
			}
		}

		if( leaves.empty()) break;

		while( !leaves.empty() )
		{
			CNode<V,E> * pN = leaves.front();
			leaves.pop();

			std::queue<CNode<V,E>*> branch;
			branch.push( pN );
			
			while( !branch.empty() )
			{

				CNode<V,E> * pN = branch.front();
				branch.pop();
				if( pN->links().empty() ) continue;
				CLink<V,E> * pL = pN->links().front();
				CNode<V,E> * pW = pL->other( pN );

				std::vector<CNode<V,E>*>::iterator niter = std::find( m_nodes.begin(), m_nodes.end(), pN );
				if( niter != m_nodes.end() )
				{
					m_nodes.erase( niter );
				}

				std::vector<CLink<V,E>*>::iterator liter = std::find( m_links.begin(), m_links.end(), pL );
				if( liter != m_links.end() )
				{
					m_links.erase( liter );
				}

				liter = std::find( pW->links().begin(), pW->links().end(), pL );
				if( liter != pW->links().end() )
				{
					pW->links().erase( liter );
				}
				
				delete pL;
				delete pN;

				if( pW->links().size() == 1 )
				{
					branch.push( pW );
				}
				
			}
		}	

	}

};



template<typename V, typename E, typename ECompare>
class CompareLink
{
public:
	bool operator()( CLink<V,E>* pL1, CLink<V,E>*pL2 )
	{
		ECompare ec;

		E  pE1 = pL1->key();
		E  pE2 = pL2->key();

		return ec( pE1, pE2 );
	}
};

template<typename V, typename E, typename ECompare, typename VCompare>
void CUndirectedGraph<V,E,ECompare,VCompare>::_maximal_spanning_tree( std::vector<E> & seeds )
{

	std::priority_queue<CLink<V,E>*, std::vector<CLink<V,E>*>, CompareLink<V,E, ECompare> > edges;

	for( unsigned int k = 0; k < m_links.size() ; k ++ )
	{
		CLink<V,E>  * pL = m_links[k];

		edges.push( pL );

	}

	m_links.clear();
	while( !edges.empty() )
	{
		CLink<V,E>  * pL = edges.top();
		edges.pop();
		m_links.push_back( pL );

	}
	std::reverse( m_links.begin(), m_links.end() );

	std::vector<CLink<V,E>*> links = m_links;

	std::set<CLink<V,E>*> active_edges;

	CUndirectedGraph<CNode<V,E>*, CLink<V,E>*,ECompare,VCompare> sub_graph;

	//insert nodes

	for( unsigned int k = 0; k < m_nodes.size(); k ++ )
	{
		CNode<V,E> * pN = m_nodes[k];
		sub_graph.insert_node( pN );
	}


	for( unsigned int k = 0; k < m_links.size(); k ++ )
	{
		CLink<V,E>* pL = m_links[k];

		CNode<V,E>* pN1 = pL->node(0);
		CNode<V,E>* pN2 = pL->node(1);
		
		CNode<CNode<V,E>*, CLink<V,E>*> * pW1 = sub_graph.find( pN1 );
		CNode<CNode<V,E>*, CLink<V,E>*> * pW2 = sub_graph.find( pN2 );

		sub_graph.insert_link( pW1, pW2, pL );

		if( sub_graph._has_cycle() )
		{
			sub_graph.remove_link( pL );
			seeds.push_back( pL->key() );
		}

	}

	std::cout << "Seeds " << seeds.size() << std::endl;
};


template<typename V, typename E, typename ECompare, typename VCompare>
void CUndirectedGraph<V,E,ECompare,VCompare>::remove_link( E & e )
{
	for( std::vector<CLink<V,E>*>::iterator liter = m_links.begin(); liter != m_links.end(); liter ++ )
	{
		CLink<V,E>* pL = *liter;
		if( pL->key() == e )
		{
			CNode<V,E> * pN1 = pL->node(0);
			std::vector<CLink<V,E>*>::iterator titer1 = std::find(pN1->links().begin(), pN1->links().end(), pL);
			pN1->links().erase( titer1 );

			CNode<V,E> * pN2 = pL->node(1);
			std::vector<CLink<V,E>*>::iterator titer2 = std::find(pN2->links().begin(), pN2->links().end(), pL);
			pN2->links().erase( titer2 );
		
			m_links.erase( liter );
			delete pL;
			return;
		}
	}
};

template<typename V, typename E, typename ECompare, typename VCompare>
bool CUndirectedGraph<V,E,ECompare,VCompare>::_has_cycle()
{
	BFS();

	for( unsigned int k = 0; k < m_links.size(); k ++ )
	{
		CLink<V,E> * pL = m_links[k];
		if( !pL->touched() ) return true;
	}
	return false;
}




template<typename V, typename E, typename ECompare, typename VCompare>
void CUndirectedGraph<V,E,ECompare,VCompare>::_maximal_spanning_tree2( std::vector<E> & seeds )
{

	std::priority_queue<CLink<V,E>*, std::vector<CLink<V,E>*>, CompareLink<V,E, ECompare> > edges;

	for( unsigned int k = 0; k < m_links.size() ; k ++ )
	{
		CLink<V,E>  * pL = m_links[k];

		edges.push( pL );

	}

	m_links.clear();
	while( !edges.empty() )
	{
		CLink<V,E>  * pL = edges.top();
		edges.pop();
		m_links.push_back( pL );

	}
	std::reverse( m_links.begin(), m_links.end() );


	for( unsigned int k = 0; k < m_nodes.size() ; k ++ )
	{
		CNode<V,E>  * pN = m_nodes[k];
		pN->label() = k;
	}

	for( unsigned int k = 0; k < m_links.size(); k ++ )
	{
		CLink<V,E>* pL = m_links[k];

		CNode<V,E>* pN1 = pL->node(0);
		CNode<V,E>* pN2 = pL->node(1);
		
		if( pN1->label() != pN2->label() )
		{
			int label_1 = pN1->label();
			int label_2 = pN2->label();

			for( unsigned int j = 0; j < m_nodes.size() ; j ++ )
			{
				CNode<V,E>  * pN = m_nodes[j];
				if( pN->label() == label_1 )
				{
					pN->label() = label_2;
				}
			}
			continue;
		}

		seeds.push_back( pL->key() );
	}

	std::cout << "Seeds " << seeds.size() << std::endl;
};

};
#endif