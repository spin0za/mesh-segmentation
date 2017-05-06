/*! \file BaseEmbed.h
 *  \brief Isometrically embed a mesh with canonical metric onto the canonical domain
 *  \author David Gu
 *  \date  documented on 10/17/2010
 *
 *  Isometrically embed a mesh with canonical metric onto the canonical domain
 */
#ifndef  _BASE_EMBED_H_
#define  _BASE_EMBED_H_

#include <vector>
#include <queue>

#include "RicciFlowMesh.h"

namespace MeshLib
{
namespace RicciFlow
{

/*! \brief CBaseEmbed class
 *
 *   embed a mesh with canonical metric onto the canonical domain
 */
template<typename M>
  class CBaseEmbed
  {
  public: 
   /*! \brief CEmbed constructor */
    CBaseEmbed( M * pMesh );
   /*! \brief CEmbed destructor */
    ~CBaseEmbed(){};
	/*! _embed the mesh */
    void _embed();
	/*! embed without normalization */
	void _embed_without_normalization();

  protected:
   
	/*! embed the first face 
	 * \param head the first face
	 */
	  virtual void _embed_first_face( typename M::CFace * head ) = 0;
	/*!
	 *
	 */
	  virtual void _embed_face( typename M::CFace * head ) = 0;

	/*! initialization */
	void _initialize();

  protected:

	/*! the mesh to be embedded */
    M * m_pMesh;

	/*! queue of faces */
	std::queue<typename M::CFace*> m_queue;
  };


//constructor
template<typename M>
CBaseEmbed<M>::CBaseEmbed( M * pMesh )
{
    m_pMesh = pMesh;
};




template<typename M>
void CBaseEmbed<M>::_initialize()
{
	for( M::MeshVertexIterator viter(  m_pMesh ); !viter.end(); viter ++)
	{
	  M::CVertex * v = *viter;
		v->huv() = CPoint2( 0, 0);
		v->touched() = false;
	}

	M::CFace * root_face = NULL;

	for( M::MeshFaceIterator fiter( m_pMesh ); !fiter.end(); fiter ++ )
	{
		M::CFace * f = *fiter;
		f->touched() = false;
		root_face = f;
	}

	_embed_first_face( root_face );
	root_face->touched() = true;

	M::CHalfEdge * he = m_pMesh->faceMostCcwHalfEdge( root_face );
	do{
		M::CHalfEdge  * sh = m_pMesh->halfedgeSym( he );
		if( sh != NULL )
		{
			M::CFace * df = m_pMesh->halfedgeFace( sh );
			m_queue.push( df );
			df->touched() = true;
		}
		he = m_pMesh->faceNextCcwHalfEdge( he );
	}while ( he != m_pMesh->faceMostCcwHalfEdge( root_face ) );

 };


template<typename M>
void CBaseEmbed<M>::_embed()
{

	_initialize();

	int count = 0;

	while( !m_queue.empty() )
	{
		M::CFace *  head = m_queue.front();
		m_queue.pop();
		assert( head->touched() );

		M::CHalfEdge * he = m_pMesh->faceMostCcwHalfEdge( head );
		do{
			M::CHalfEdge * sh = m_pMesh->halfedgeSym( he );
			if( sh != NULL )
			{	  
				M::CFace * df = m_pMesh->halfedgeFace( sh );
				if( !df->touched() )
				{
					m_queue.push( df );
					df->touched() = true;
				}
			}
			he = m_pMesh->halfedgeNext( he );
		}while( he != m_pMesh->faceMostCcwHalfEdge( head ) );

		_embed_face( head );
	}

	//normalize the uv coordinates

	double u_min, u_max, v_min, v_max;
	u_min = v_min = 1e30;
	u_max = v_max = -1e30;


	for( M::MeshVertexIterator viter( m_pMesh ); !viter.end(); viter++ )
	{
		M::CVertex * pV = *viter;
		CPoint2 uv = pV->huv();
		double u_param = uv[0];
		double v_param = uv[1];
		if( u_min>u_param )
			u_min = u_param;
		if( u_max<u_param )
			u_max = u_param;
		if( v_min>v_param )
			v_min = v_param;
		if( v_max<v_param )
			v_max = v_param;
	}

	double range = ( u_max-u_min )>( v_max-v_min ) ? u_max-u_min : v_max-v_min;

	printf( "range = %lf\n", range );

	if( range>1e-6 ) {
		for( M::MeshVertexIterator viter( m_pMesh ); !viter.end(); viter++ )
		{
			M::CVertex * pV = *viter;
			pV->u() -= log( range );
			CPoint2 uv =  pV->huv();
			double u_param = uv[0];
			double v_param = uv[1];
			uv[0] = (u_param - u_min)/range;
			uv[1] = (v_param - v_min)/range;
			pV->huv() = uv;
		}
	}
};

template<typename M>
void CBaseEmbed<M>::_embed_without_normalization()
{

	_initialize();

	int count = 0;

	while( !m_queue.empty() )
	{
		M::CFace *  head = m_queue.front();
		m_queue.pop();
		assert( head->touched() );

		M::CHalfEdge * he = m_pMesh->faceMostCcwHalfEdge( head );
		do{
			M::CHalfEdge * sh = m_pMesh->halfedgeSym( he );
			if( sh != NULL )
			{	  
				M::CFace * df = m_pMesh->halfedgeFace( sh );
				if( !df->touched() )
				{
					m_queue.push( df );
					df->touched() = true;
				}
			}
			he = m_pMesh->halfedgeNext( he );
		}while( he != m_pMesh->faceMostCcwHalfEdge( head ) );

		_embed_face( head );
	}

};

} //namespace RicciFlow
} //namespace MeshLib

#endif 