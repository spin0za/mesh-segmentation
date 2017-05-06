/*!
*      \file Fuchs.h
*      \brief Compute Fuchs transformation
*	   \author David Gu
*      \date Documented 04/11/2011
*
*/

#ifndef  _FUCHS_H_
#define  _FUCHS_H_

#include <map>
#include <vector>
#include <algorithm>

#include "Mobius/Mobius.h"
#include "Riemannian/RicciFlow/RicciFlowMesh.h"
#include "DeckTransformation.h"
#include "TorusDeckTransformation.h"

namespace MeshLib
{

namespace Topology
{

  template<typename M>
  class CSegment
  {
  public:
	  CSegment( M * pMesh ) { m_pMesh = pMesh; m_dual = NULL; };
	  ~CSegment(){};

	  void _insert( typename M::CHalfEdge * pH ) { m_hedges.push_back( pH ); };
	
	  CSegment<M> * & dual() { return m_dual; };
		
	  //verify if pS  is dual to the current segment 
	  bool _dual( CSegment<M> * pS );
		
	  std::string & name() { return m_name; };

	  typename M::CVertex * _source();
	  typename M::CVertex * _target();

	  std::list<typename M::CHalfEdge*> & halfedges() { return m_hedges; };

  protected:

	  std::list<typename M::CHalfEdge*> m_hedges;
	  CSegment<M> * m_dual;
	  M * m_pMesh;
	  std::string m_name;
  };

template<typename M>
inline bool CSegment<M>::_dual( CSegment<M> * pSeg )
{
	M::CHalfEdge * pH = m_hedges.front();

	M::CHalfEdge * pD = pSeg->m_hedges.back();
	

	M::CVertex * pS = m_pMesh->halfedgeSource( pH );
	M::CVertex * pT = m_pMesh->halfedgeTarget( pH );

	int       fS = pS->father();
	int       fT = pT->father();

	M::CVertex * dS = m_pMesh->halfedgeSource( pD );
	M::CVertex * dT = m_pMesh->halfedgeTarget( pD );
	
	int       dfS = dS->father();
	int       dfT = dT->father();

	return( fS == dfT && fT == dfS );
};

template<typename M>
inline typename M::CVertex* CSegment<M>::_source()
{
  if( m_hedges.empty() ) return NULL;
  M::CHalfEdge * pH = m_hedges.front();		
  return m_pMesh->halfedgeSource( pH );
};

template<typename M>
inline typename M::CVertex * CSegment<M>::_target()
{
  if( m_hedges.empty() ) return NULL;
  M::CHalfEdge * pL = m_hedges.back();		
  return m_pMesh->halfedgeTarget( pL );
};

/*!
*	\brief CTorusDeckTransformation
*	
*	Deck Transformation Group for Torus
* 
*/
template<typename M>
class CFuchs :public CDeckTransformation<M>
{

public:
	/*!
	 *	CFuchs constructor
	 *  \param pMesh    closed mesh
	 *  \param pDomain  foundamental domain
	 */
	CFuchs( M * pMesh, M * pDomain ): CDeckTransformation<M>( pMesh, pDomain) {};
	/*!
	 *	CFuchs destructor
	 */
	~CFuchs() 
	{
		for( unsigned int i = 0; i < m_segments.size(); i ++ )
		{
			CSegment<M> * pS = m_segments[i];
			delete pS;
		}	
	};

	/*!
	 *	Compute canoical Fuchs group generators
	 */

	void _calculate_generators();
	/*!
	 *	Output Fuchs group generators
	 */

	void _output_generators( const char * prefix );
	/*!
	 *	Generate fundamental domains transformed by Fuchs generators
	 */

	void HyperbolicTessellation(  const char * prefix  );

protected:
	/*!
	 *	Compute canonical fundamental group generators from the cut graph (boundary of pDomain )
	 */
	void _canonical_basis();
	/*!
	 *	Find the 2g corners on the boundary, whose fathers are the base point on pMesh
	 */
	int _locate_corners( std::vector<typename M::CVertex*> & corners );
	/*!
	 *	Segment the boundary loop to 4g segments, starting and ending at corners	
	 */
	void _segment( int father_of_center );
	/*!
	 *	Compute the Mobius transformation, which maps A to the origin, B to a real number
	 */
	CMobius MT( CPoint A, CPoint B);
	/*!
	 *	Compute the Mobius transformation, which maps (A,B) to (IA,IB)
	 */
	CMobius _transformation( CPoint A, CPoint B, CPoint IA, CPoint IB );
	/*!
	 *	Transform the fundamental domain mesh by the Mobius transformation mob, and color all the vertices by rgb
	 */
	void _transform( CMobius mob, CPoint rgb );

	/*!
	 *	Array of segments
	 */
	std::vector<CSegment<M>*> m_segments;
	/*!
	 *	Array of Fuchs group generators
	 */
	std::vector<CMobius>   m_mobius;
	/*!
	 *	Genus of the surface
	 */
	int m_genus;

};

template<typename M>
int CFuchs<M>::_locate_corners( std::vector<typename M::CVertex*> & corners )
{

	std::priority_queue<M::CVertex*, std::vector<M::CVertex*>, CompareRicciFlowVertex<M::CVertex>> vqueue;

	for( unsigned int i = 0; i < m_boundary.loops().size(); i ++ )
	{
		M::CLoop * pL = m_boundary.loops()[i];

		for( std::list<M::CHalfEdge*>::iterator hiter = pL->halfedges().begin(); hiter != pL->halfedges().end(); hiter ++ )
		{
			M::CHalfEdge * pH = *hiter;
			M::CVertex * pV = m_pDomain->halfedgeTarget( pH );
			vqueue.push( pV );
		}
	}
	
	while( !vqueue.empty() )
	{
		M::CVertex * pV0 = vqueue.top();
		vqueue.pop();
		corners.push_back( pV0 );

		M::CVertex * pV1 = vqueue.top();
		vqueue.pop();
		corners.push_back( pV1 );

		if( pV0->father() != pV1->father() )
		{
			printf("The input mesh doesn't have canonical cut graph. Check the input ! \n" );
			return -1;
		}
		
		int count = 2;
		while( pV0->father() == vqueue.top()->father() )
		{
			corners.push_back( vqueue.top() );
			vqueue.pop();
			count ++;
		}

		if( count > 2 )
		{
			printf("Locate the center, there are %d corners.\n", count );
			return pV0->father();
		}

		corners.clear();

	}

	return -1;
};

//divide the boundary loop to segments
template<typename M>
void CFuchs<M>::_segment( int father_of_center )
{
	M::CLoop * pL = m_boundary.loops().front();

	std::list<M::CHalfEdge*> hedges;
	for( std::list<M::CHalfEdge*>::iterator hiter = pL->halfedges().begin(); hiter != pL->halfedges().end(); hiter ++ )
	{
		M::CHalfEdge * pH = *hiter;
		hedges.push_back( pH );
	}

	//move the head half_edge starts from the center
	while( true )
	{
		std::list<M::CHalfEdge*>::iterator hiter = hedges.begin();
		M::CHalfEdge * pH = *hiter;
		
		if( m_pDomain->halfedgeSource( pH )->father()  == father_of_center ) break;
		
		hedges.erase( hiter );
		hedges.push_back( pH );
	}

	//divide the loop to segments

	while( !hedges.empty() )
	{
		CSegment<M> * pS = new CSegment<M>( m_pMesh );
		assert( pS );

		M::CHalfEdge * pH = NULL;
		do{
			std::list<M::CHalfEdge*>::iterator hiter = hedges.begin();
			pH = *hiter;
			hedges.erase( hiter );
			pS->_insert( pH );
		}while( m_pDomain->halfedgeTarget( pH )->father() != father_of_center && !hedges.empty() );
		m_segments.push_back( pS );
	}

	printf("There are %d segments\n", m_segments.size() );
	
	//match segments

	for( unsigned int i = 0; i < m_segments.size(); i ++ )
	{
		CSegment<M> * pS = m_segments[i];
	
		if( pS->dual() != NULL ) continue;
		
		char line[256];
		sprintf(line, "%d", i );
		pS->name() = line;

		for( unsigned int j = i + 1; j < m_segments.size(); j ++ )
		{
			CSegment<M> * pD = m_segments[j];
			
			if( pS->_dual( pD ) )
			{
				pS->dual() = pD;
				pD->dual() = pS;
				pD->name() = "-" + pS->name();
			}
		}
	}


	unsigned int head = 0;

	for( unsigned int i = 0; i < m_segments.size(); i ++ )
	{
		CSegment<M> * pS = m_segments[(i+0)%m_segments.size()];
		CSegment<M> * pD = m_segments[(i+2)%m_segments.size()];

		CSegment<M> * nS = m_segments[(i+1)%m_segments.size()];
		CSegment<M> * nD = m_segments[(i+3)%m_segments.size()];
		
		if( pS->dual() == pD && nS->dual() == nD )
		{
			head = i;
			break;
		}
	}

	std::vector<CSegment<M>*> tsegments;

	for( unsigned int i = 0; i < m_segments.size(); i ++ )
	{
		tsegments.push_back( m_segments[(i+head)%m_segments.size()] );
	}

	m_segments.clear();

	for( unsigned int i = 0; i < tsegments.size(); i ++ )
	{
		CSegment<M> * pS = tsegments[i];
		m_segments.push_back( pS );
		printf("%s ", pS->name().c_str() );
	}

	m_genus = m_segments.size()/4;
	

}

template<typename M>
inline void CFuchs<M>::_canonical_basis()
{
	std::vector<M::CVertex*> corners;
	int father_of_center = _locate_corners( corners );
	_segment( father_of_center );
};

template<typename M>
inline void CFuchs<M>::_calculate_generators()
{
	_canonical_basis();

	for( unsigned int i = 0; i < m_segments.size(); i ++ )
	{
		CSegment<M> * pS = m_segments[i];
		CSegment<M> * pD = pS->dual();
		
		M::CVertex *   S = pS->_source();
		M::CVertex *   T = pS->_target();

		M::CVertex *  dS = pD->_source();
		M::CVertex *  dT = pD->_target();

		CPoint     vS =  S->point();
		CPoint     vT =  T->point();

		CPoint    vdS = dS->point();
		CPoint    vdT = dT->point();

		CMobius mob = _transformation( vS, vT, vdT, vdS );
		
		m_mobius.push_back( inverse(mob) );
	};
};

template<typename M>
inline CMobius CFuchs<M>::MT( CPoint A, CPoint B)
{
	Complex  za( A[0], A[1] );
	CMobius  m(za,0); 
	Complex  zb( B[0], B[1]);
	zb = m * zb;
	Complex theta= std::conj( zb );

	return CMobius( za, std::arg( theta) );

}

template<typename M>
inline CMobius CFuchs<M>::_transformation( CPoint A, CPoint B, CPoint IA, CPoint IB )
{
	CMobius m0 = MT(A,B);
	CMobius m1 = MT(IA,IB);
	CMobius im1 = inverse( m1 );
	//CMobius m = m0 * im1;
	CMobius m = im1 * m0;

	return m;
}

template<typename M>
void CFuchs<M>::_transform( CMobius mob, CPoint rgb )
{
	for( M::MeshVertexIterator viter( m_pDomain ); !viter.end(); viter ++ )
	{
		M::CVertex * pV = *viter;
		CPoint2		p  = pV->huv();
		Complex   pc( p[0], p[1] );
		pc = mob * pc;
		pV->huv() = CPoint2( pc.real(), pc.imag() );
		pV->rgb() = rgb;
	}
};

template<typename M>
void CFuchs<M>::HyperbolicTessellation( const char * prefix )
{

	float color[][3]={{1,0,0},{0,1,0},{0,0,1},{0.5,0,0},{1,0,1},{0,1,1},{1,0,0.5},{0.5,1,0.5}};
	
	int id = 0;

	int * order = new int[4*m_genus];

	{
		for( int i = m_genus; i > 0; i -- )
		{
			int A  = (i%m_genus) * 4 + 0;
			int B  = (i%m_genus) * 4 + 1;
			int NA = (i%m_genus) * 4 + 2;
			int NB = (i%m_genus) * 4 + 3;

			order[ id++ %( 4 * m_genus ) ] =  NB;
			order[ id++ %( 4 * m_genus ) ] =  A;
			order[ id++ %( 4 * m_genus ) ] =  B;
			order[ id++ %( 4 * m_genus ) ] =  NA;
		}

		int head = order[0];

		for( int i = 0; i < m_genus * 4 - 1; i ++ )
		{
			order[i] = order[i+1];
		}
		order[4*m_genus-1] = head;

		printf("\n");
		for( int i = 0; i < m_genus * 4; i ++ )
		{
			printf("%d ", order[i] );
		}
		printf("\n");
	}

	//for genus 2
	//int order[8] = {0,1,2,7,4,5,6,3};
    //for genus=3
	//int order[12]={ 0,1,2,11,8,9,10,7,4,5,6,3}
	int color_idx = 0;
	for( int k = 0; k < m_genus * 4; k ++ )
	{
		for( int i = 0; i < m_genus * 4 - 2; i ++ )
		{

			for( M::MeshVertexIterator viter( m_pDomain ); !viter.end(); viter ++ )
			{
				M::CVertex * pV = *viter;
				CPoint    p  = pV->point();
				pV->huv() = CPoint2( p[0], p[1] );
			}

			CPoint rgb;
			for( int idx = 0; idx < 3; idx ++ )
			{
				rgb[idx] = color[color_idx%8][idx];
			}
			color_idx ++;

			for( int j = i; j >=0; j-- )
			{
				int idx = (k + j + m_genus * 4 )% ( m_genus * 4);
				idx = order[idx];
				_transform( m_mobius[idx], rgb );
			}
			/*
			for( MeshVertexIterator viter( m_pMesh ); !viter.end(); viter ++ )
			{
				CVertex * pV = *viter;
				CPoint2   p  = v_uv( pV );
				CPoint    n  = v_normal( pV );
				char line[1024];

				sprintf( line, "normal=(%g %g %g) uv=(%g %g) rgb=(%g %g %g)", n[0], n[1], n[2], p[0], p[1], rgb[0], rgb[1], rgb[2]); 
				pV->string() = line;
			}
			*/
			char name[256];
			sprintf(name,"%s_%d_%d.m",prefix, k,i );
			m_pDomain->write_m( name );
		}
	}

	delete []order;
};

template<typename M>
inline void CFuchs<M>::_output_generators( const char * name )
{
	FILE * fp = fopen( name, "w");
	assert( fp != NULL );

	fprintf(fp, "genus=%d\n",m_genus );
	for( unsigned int i = 0; i < m_mobius.size(); i ++ )
	{
		CMobius & mob = m_mobius[i];
		fprintf(fp, "%d=(%.12g %.12g %.12g)\n", i, mob.z().real(), mob.z().imag(), std::arg( mob.theta() ));
	}
	fclose( fp );
}

} //namespace Topology
} //namespace MeshLib
#endif  _FUCHS_H_