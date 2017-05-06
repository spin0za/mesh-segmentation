#ifndef _COVERING_MESH_H_
#define _COVERING_MESH_H_

#include <math.h>
#include <assert.h>
#include <iostream>
#include <fstream>
#include <list>
#include <vector>
#include <map>

#include "Mobius/Mobius.h"
#include "Mesh/BaseMesh.h"
#include "Mesh/Vertex.h"
#include "Mesh/HalfEdge.h"
#include "Mesh/Edge.h"
#include "Mesh/Face.h"
#include "mesh/iterators.h"
#include "mesh/boundary.h"
#include "Parser/parser.h"
#include "Parser/traits_io.h"
#include "Geometry/Line2D.h"

namespace MeshLib{

template<typename T>
void __swap( T & a, T & b ) { T c = a; a = b; b = c; };

namespace Topology
{

template<typename V, typename E, typename F, typename H, typename M>
class CoveringSpace
{
public:

	class CFace
	{
	public:
		CFace(){};
		~CFace(){};
		F * & base() { return m_base; };
		CHalfEdge * & halfedge() { return m_halfedge; };

	protected:
		F * m_base;
		CHalfEdge * m_halfedge;
	};

	class CEdge
	{
	public:
		CEdge()
		{
			for( int i = 0; i < 2; i ++ )
			{
				m_halfedges[i] = NULL;
			}
		};
		~CEdge(){};

		E * & base() { return m_base; };	

		CHalfEdge * & halfedge( int k ) { return m_halfedges[k]; };

		double      & length() { return m_length; };

	protected:
		E * m_base;
		CHalfEdge * m_halfedges[2];
		double m_length;
	};


	class CHalfEdge
	{
	public:
		CHalfEdge(){ m_edge = NULL; m_face = NULL; m_next = NULL; m_prev = NULL; m_vert = NULL; };
		~CHalfEdge(){};

		H *& base(){ return m_base; };

		CHalfEdge *& prev(){ return m_prev; };
		CHalfEdge *& next(){ return m_next; };
		CEdge	  *& edge(){ return m_edge; };
		CFace	  *& face(){ return m_face; };

		CVertex	  *& target(){ return m_vert; };
		CVertex	  *& source(){ return m_prev->target(); };

		CHalfEdge *  dual() 
		{
			if( m_edge == NULL ) return NULL;

			return (m_edge->halfedge(0) != this )? m_edge->halfedge(0): m_edge->halfedge(1);
		};

		double & angle() { return m_angle; };

	protected:
		
		CVertex   * m_vert;
		CHalfEdge * m_prev;
		CHalfEdge * m_next;
		CFace     * m_face;
		CEdge     * m_edge;
		H * m_base;

		double      m_angle;
	};

	class CVertex
	{
	public:
		 CVertex(){};
		 ~CVertex(){};
		 int & id(){ return m_id; };

		 V *& base(){ return m_base; };
		 CHalfEdge * & halfeddge() { return m_halfedge; };
		 std::complex<double> & uv() { return m_uv; };

	protected:
		
		V * m_base;
		int m_id;

		CHalfEdge *			 m_halfedge;
		std::complex<double> m_uv;
		
	};

	//hyperbolic cosine law
	class hyperbolic_cosine_law
	{
	public:
		double operator()( double a, double b, double c )
		{	
			return acos( (cosh(a) * cosh(b)-cosh(c) )/( sinh(a)*sinh(b)) );
		};
	};

	class compareEdgeMeshPair
	{
	public:
		bool operator()( std::pair<E*,M*> e1, std::pair<E*,M*> e2 )
		{
			M * pMesh = e1.second;
			E * pE1   = e1.first;
			E * pE2   = e2.first;

			V * s1 = pMesh->edgeVertex1( pE1 );
			V * t1 = pMesh->edgeVertex2( pE1 );
			
			int sf1 = (s1->father()< t1->father() )? s1->father():t1->father();
			int tf1 = (s1->father()< t1->father() )? t1->father():s1->father();

			V * s2 = pMesh->edgeVertex1( pE2 );
			V * t2 = pMesh->edgeVertex2( pE2 );

			int sf2 = (s2->father()< t2->father() )? s2->father():t2->father();
			int tf2 = (s2->father()< t2->father() )? t2->father():s2->father();


			if( sf1 < sf2 ) return true;
			if( sf1 > sf2 ) return false;

			return tf1 < tf2;
		}
	};

	class compareHalfEdge
	{
	public:
		bool operator()(CHalfEdge * pH1, CHalfEdge * pH2)
		{
			CVertex * s1 = pH1->source();
			CVertex * t1 = pH1->target();
			
			int sid1 = (s1->base()->id() < t1->base()->id())?s1->base()->id():t1->base()->id();
			int tid1 = (s1->base()->id() > t1->base()->id())?s1->base()->id():t1->base()->id();

			CVertex * s2 = pH2->source();
			CVertex * t2 = pH2->target();

			int sid2 = (s2->base()->id() < t2->base()->id())?s2->base()->id():t2->base()->id();
			int tid2 = (s2->base()->id() > t2->base()->id())?s2->base()->id():t2->base()->id();
			
			if( sid1 < sid2 ) return true;
			if( sid1 > sid2 ) return false;

			return tid1 < tid2;
		}
	};

	class CNeighborhood
	{
	public:
		CNeighborhood( V * pV, M * pM );
		~CNeighborhood();
		
		std::vector<CFace*> & faces() { return m_faces; };
		std::vector<CVertex*> & vertices(){ return m_verts; };
		void operator*=( CMobius & mob );
		CVertex *& center() { return m_center; };
		void _embed();

	CFace* createFace( std::vector<CVertex*> & verts )
	{
		CFace * pF = new CFace;
		assert( pF );
		m_faces.push_back(pF);

		std::vector<CHalfEdge*> hes;

		for( size_t i = 0; i < 3; i ++ )
		{
			CHalfEdge * pH = new CHalfEdge;
			assert( pH );
			m_halfedges.push_back( pH );
			pH->target() = verts[i];
			pH->face()   = pF;
			hes.push_back( pH );
		}

		for( size_t i = 0; i < 3; i ++ )
		{
			hes[i]->next() = hes[(i+1)%3];
			hes[i]->prev() = hes[(i+2)%3];
		}

		pF->halfedge() = hes[0];

		return pF;
	};

	protected:

		CVertex * m_center;
		std::vector<CVertex*> m_verts;
		std::vector<CFace*>   m_faces;
		std::list<CHalfEdge*> m_halfedges;
		std::list<CEdge*>     m_edges;

		M*					  m_pMesh;
	};

	class compareFaceInNeighborhood
	{
	public:
		bool operator()(std::pair<CFace*,CNeighborhood*> pF1, std::pair<CFace*,CNeighborhood*> pF2)
		{
			return pF1.first->base()->id() < pF2.first->base()->id();
		}
	};


	class CFundamentalDomain
	{
	public:
		CFundamentalDomain( M * pM, M * pD );
		~CFundamentalDomain(){};
		
		F* dual_face( F * pF ) { return m_face_map[pF]; };
		E* dual_edge( E * pE ) { return m_edge_map[pE]; };
		V* base()   { return m_base;   };
		V* origin() { return m_origin; };
	protected:
		
		M * m_pMesh;
		M * m_pDomain;
		//find covering face for m_pDomain, m_pMesh
		std::map<F*,F*> m_face_map;
		std::map<E*,E*> m_edge_map;
		
		V * m_base;	//base point on the closed mesh
		V * m_origin; //origin point on the domain

	protected:
		//construct the face map
		void _construct_face_map();
		//construct the edge map
		void _construct_edge_map();
		//locate base
		void _locate_base();
	};


	class CLiftedPath
	{
	public:
		CLiftedPath( M* pM, M * pD, const char * input );
		~CLiftedPath();

		void _embed(  CFundamentalDomain * pD );

		//output the path to a mesh file
		void _output( const char * output_file );

		CMobius _fuchsian_transformation();

		void _straighten(  CFundamentalDomain * pD, const char * name, const char * plane_name );
		void _ray(  CFundamentalDomain * pD, CPoint2 S, CPoint2 T, std::vector<std::pair<CPoint,CPoint2> > & pts );


	protected:
		std::vector<CNeighborhood*> m_neighborhoods;
		
		//find the overlapping faces, each pair, first face is in pN1, 2nd face is in pN2
		void _overlap( CNeighborhood * pN1, CNeighborhood * pN2, std::vector<std::pair<CFace*,CFace*>> & pairs );

		CMobius _transform( CFace * pF0, CFace * pF1 );
		CMobius MT( std::complex<double> za, std::complex<double> zb );
		CMobius _transformation( std::complex<double> A, std::complex<double> B, std::complex<double> IA, std::complex<double> IB );

		//if the line intersecting the edge, true - intersect, false- non-intersect
		//q    - intersection point, lambda is the ratio in the line segment from pS to pT

		bool _intersect( CPoint2 pS, CPoint2 pT, E * pE, CPoint2 & q, double & lambda );
		
		//funchsian transform the whole fundamental domain
		void _Fuchsian_Transform(  E* eS, E * eT );

		M * m_pMesh;
		M * m_pDomain;
	};
};


#include "FundamentalDomain.hpp"
#include "Neighborhood.hpp"
#include "LiftedPath.hpp"

} //namespace Topology
} //namespace MeshLib

#endif
