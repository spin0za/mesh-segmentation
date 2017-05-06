/*! \file GraphEmbedding.h
 *  \brief Euclidean Ricci flow algorithm for planar graph embedding
 *  \author David Gu
 *  \date   documented on 01/10/2014
 *
 *	Algorithm for graph embedding
 */

#ifndef _GRAPH_EMBEDDING_H_
#define _GRAPH_EMBEDDING_H_

#include <map>
#include <vector>
#include <Eigen/Sparse>

#include "Mesh/BaseMesh.h"
#include "Mesh/Vertex.h"
#include "Mesh/HalfEdge.h"
#include "Mesh/Edge.h"
#include "Mesh/Face.h"
#include "mesh/iterators.h"
#include "mesh/boundary.h"
#include "Parser/parser.h"
#include "RicciFlowMesh.h"
#include "BaseRicciFlow.h"


namespace MeshLib
{

namespace RicciFlow
{

/*! \brief Class CGraphEmbedding
*
*	Algorithm for computing graph embedding
*/
template<typename M>
class CGraphEmbedding
  {
  public:
    /*! \brief CGraphEmbedding
	 *  \param pMesh the input mesh
	 */
	  CGraphEmbedding( M * pMesh );
    /*! \brief CTangentialRicciFlow destructor
	 */
	  ~CGraphEmbedding(){};
	/*!	union the graph with the dual graph
	 */
	void graph_union_with_dual_graph( const char * output );
	/*! remove north pole, which is an edge vertex
	 */
	void remove_northpole( const char * output );

  protected:

	  /*!
	   *	pointer to the mesh
	   */
	  M * m_pMesh;
	  /*!
	   *	root vertex
	   */
	  typename M::CVertex * m_root;
	  /*!
	   *	two face-vertex adjacent to the root
	   */
	  std::vector<typename M::CVertex*> m_face_verticies;
	  /*!
	   *	two vertex-vertex adjacent to the root
	   */
	  std::vector<typename M::CVertex*> m_vertex_verticies;
	  /*!
	   *	four edge-vertex adjacent to the root
	   */
	  std::vector<typename M::CVertex*> m_edge_verticies;

	  /*!
	   *	find neighbors to the root
	   */
	  void _find_neighbors( typename M::CVertex * root );

	  /*!
	   *	stereo-graphic projection
	   */
	  void _stereo_graphic_projection( const CPoint & p, CPoint2 & q );

	  /*!
	   *	inverse stereo-graphic projection
	   */
	  void _inverse_stereo_graphic_projection( const CPoint2 & q, CPoint & p );
	
  };

template<typename M>
void CGraphEmbedding<M>::_stereo_graphic_projection( const CPoint & p, CPoint2 & q )
{
	q = CPoint2( p[0]/(1-p[2]); p[1]/(1-p[2]) );
};

template<typename M>
void CGraphEmbedding<M>::_inverse_stereo_graphic_projection( const CPoint2 & q, CPoint & p )
{
	double u = q[0];
	double v = q[1];

	p = CPoint( (2*u)/(1+u*u+v*v), (2*v)/(1+u*u+v*v), (1-u*u-v*v)/(1+u*u+v*v) );
};

template<typename M>
CGraphEmbedding<M>::CGraphEmbedding( M * pMesh )
{
	m_pMesh = pMesh;
};


template<typename M>
void CGraphEmbedding<M>::graph_union_with_dual_graph( const char * output )
{
	int idx = 1;

	for( M::MeshVertexIterator viter( m_pMesh ); !viter.end(); viter ++ )
	{
		M::CVertex * pV = *viter;
		pV->idx() = idx ++;
	}

	for( M::MeshEdgeIterator eiter( m_pMesh ); !eiter.end(); eiter ++ )
	{
		M::CEdge * pE  = *eiter;
		pE->idx() = idx ++;
	}

	for( M::MeshFaceIterator fiter( m_pMesh ); !fiter.end(); fiter ++ )
	{
		M::CFace * pF = *fiter;
		pF->idx() = idx ++;
	}

	std::fstream _os( output, std::fstream::out );
	if( _os.fail() )
	{
		std::cout << "Error is opening file " << output;
		return;
	}

	for( M::MeshVertexIterator viter( m_pMesh ); !viter.end(); viter ++ )
	{
		M::CVertex * v = *viter;
		_os << "Vertex " << v->idx();
		
		for( int i = 0; i < 3; i ++ )
		{
			_os << " " << v->point()[i];
		}
		_os << " {type=(0) K=(0)}" << std::endl;
	}

	for( M::MeshEdgeIterator eiter( m_pMesh ); !eiter.end(); eiter ++ )
	{
		M::CEdge * e = *eiter;
		_os << "Vertex " << e->idx();
		
		M::CVertex * pV1 = m_pMesh->edgeVertex1( e );
		M::CVertex * pV2 = m_pMesh->edgeVertex2( e );

		MeshLib::CPoint p = pV1->point() + pV2->point();

		p /= 2.0;

		for( int i = 0; i < 3; i ++ )
		{
			_os << " " << p[i];
		}	
		_os << " {type=(1) K=(0)}" << std::endl;
	}

	for( M::MeshFaceIterator fiter( m_pMesh ); !fiter.end(); fiter ++ )
	{
		M::CFace * f = *fiter;
		_os << "Vertex " << f->idx();

		MeshLib::CPoint p(0,0,0);
		for( M::FaceVertexIterator fviter( f ); !fviter.end(); fviter ++ )
		{
			M::CVertex * pV = *fviter;
			p += pV->point();
		}
		p /= 3.0;

		for( int i = 0; i < 3; i ++ )
		{
			_os << " " << p[i];
		}	
		_os << " {type=(2) K=(0)}" << std::endl;
	}

	int fid = 1;

	for( M::MeshFaceIterator fiter( m_pMesh ); !fiter.end(); fiter ++ )
	{
		M::CFace * f = *fiter;

		for( M::FaceHalfedgeIterator fhiter( f); !fhiter.end(); fhiter ++ )
		{
			M::CHalfEdge * pH = *fhiter;
			M::CEdge     * pE = m_pMesh->halfedgeEdge( pH );

			M::CVertex   * pS = m_pMesh->halfedgeSource( pH );
			M::CVertex   * pT = m_pMesh->halfedgeTarget( pH );
			
			_os << "Face " << fid++;
			_os << " " << f->idx()<< " " << pS->idx() << " " << pE->idx() << std::endl;

			_os << "Face " << fid++;
			_os << " " << f->idx()<< " " << pE->idx() << " " << pT->idx() << std::endl;
		}
	}

	_os.close();

};

template<typename M>
void CGraphEmbedding<M>::_find_neighbors( typename M::CVertex * root )
{

	for( M::VertexVertexIterator vviter( root ); !vviter.end(); vviter ++ )
	{
		M::CVertex * pW = *vviter;
		
		if( pW->type() == 0 )
		{
			m_vertex_verticies.push_back( pW );
			continue;
		}
		if( pW->type() == 2 )
		{
			m_face_verticies.push_back( pW );
			continue;
		}
	}

	for( M::VertexFaceIterator vfiter( root ); !vfiter.end(); vfiter ++ )
	{
		M::CFace * pF = *vfiter;
		for( M::FaceHalfedgeIterator fhiter( pF ); !fhiter.end(); fhiter ++ )
		{
			M::CHalfEdge * pH = *fhiter;
			M::CVertex * pS = m_pMesh->halfedgeSource( pH );
			M::CVertex * pT = m_pMesh->halfedgeTarget( pH );
			
			if( pS != root && pT != root )
			{
				M::CHalfEdge * pSH = m_pMesh->halfedgeSym( pH );
				M::CHalfEdge * pNH = m_pMesh->halfedgeNext( pSH );
				M::CVertex   * pSV = m_pMesh->halfedgeTarget( pNH );
				std::cout << pSV->id() << std::endl;
				m_edge_verticies.push_back( pSV );
				break;
			}
		}
	}
}

template<typename M>
void CGraphEmbedding<M>::remove_northpole( const char * output )
{

	//find root
	for( M::MeshVertexIterator viter( m_pMesh ); !viter.end(); viter ++ )
	{
		M::CVertex * pV = *viter;
		if( pV->type() == 1 )
		{
			m_root = pV;
			//break;
		}
	}

	assert( m_root != NULL );

	_find_neighbors( m_root );

	//there are two vertex-verticies and face-vertices attaching to the root
	//remove all the faces adjacent to them
	for( M::MeshFaceIterator fiter( m_pMesh ); !fiter.end(); fiter ++ )
	{
		M::CFace * pF = *fiter;
		pF->touched() = false;
	}

	for( size_t i = 0; i < m_face_verticies.size(); i ++ )
	{
		M::CVertex * pW = m_face_verticies[i];
		
		for( M::VertexFaceIterator vfiter( pW ); !vfiter.end(); vfiter ++ )
		{
			M::CFace * pF = *vfiter;
			pF->touched() = true;
		}
	}

	for( size_t i = 0; i < m_vertex_verticies.size(); i ++ )
	{
		M::CVertex * pW = m_vertex_verticies[i];
		
		for( M::VertexFaceIterator vfiter( pW ); !vfiter.end(); vfiter ++ )
		{
			M::CFace * pF = *vfiter;
			pF->touched() = true;
		}
	}

	//find the four edge-vertices, which are against the root
	for( M::MeshVertexIterator viter( m_pMesh ); !viter.end(); viter ++ )
	{
		M::CVertex * pV = *viter;
		pV->target_k() = 0.0;
	}

	for( size_t i= 0; i < m_edge_verticies.size(); i++ )
	{
		m_edge_verticies[i]->target_k()= PI/2.0;
	}

	m_root->type() = -1;

	for( size_t i = 0; i < m_face_verticies.size(); i ++ )
	{
		m_face_verticies[i]->type() = -4;
	}

	for( size_t i = 0; i < m_face_verticies.size(); i ++ )
	{
		m_vertex_verticies[i]->type() = -2;
	}

	std::fstream _os( output, std::fstream::out );
	if( _os.fail() )
	{
		std::cout << "Error is opening file " << output;
		return;
	}

	for( M::MeshVertexIterator viter( m_pMesh ); !viter.end(); viter ++ )
	{
		M::CVertex * v = *viter;
		if( v->type() < 0 )
		{
			_os << "#Vertex " << v->id();
		}
		else
		{
			_os << "Vertex " << v->id();
		}

		for( int i = 0; i < 3; i ++ )
		{
			_os << " " << v->point()[i];
		}
		_os << " {type=(" << v->type() << ") K=(" << v->target_k()<<")}" << std::endl;
	}



	for( M::MeshFaceIterator fiter( m_pMesh ); !fiter.end(); fiter ++ )
	{
		M::CFace * f = *fiter;
		
		if( !f->touched() )
		{
			_os << "Face " << f->id();
		}
		else
		{
			_os << "#Face "<< f->id();
		}
		
		for( M::FaceHalfedgeIterator fhiter( f); !fhiter.end(); fhiter ++ )
		{
			M::CHalfEdge * pH = *fhiter;
			M::CVertex   * pV = m_pMesh->halfedgeTarget( pH );
			_os << " " << pV->id();
		}
					
		_os << std::endl;
	}

	_os.close();

};


} //namespace RicciFlow

} //namespace MeshLib	

#endif  _GRAPH_EMBEDDING_H_