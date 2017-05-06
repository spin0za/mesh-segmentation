/*! \file DelaunayMesh.h
*   \brief Mesh Generator using Delaunay Triangulation
*   \date  Documented on 10/14/2010
*
*   Generate planar triangular mesh with input PLSG (Piecewise Linear Segment Graph )
*/

#ifndef  _DELAUNAY_MESH_H_
#define  _DELAUNAY_MESH_H_

#include <map>
#include <vector>
#include <queue>

#include "DelaunayHalfEdge.h"
#include "Mesh/BaseMesh.h"
#include "Mesh/boundary.h"
#include "Mesh/iterators.h"
#include "Mesh/DynamicMesh.h"
#include "predicates.h"
#include "poly.h"


#define VERYLARGE 1e+10

#ifndef PI
#define PI 3.14159265358979323846
#endif
#ifndef TWOPI
#define TWOPI (2.0*PI)
#endif

namespace MeshLib
{
  
/*! \brief CDeluanayMesh class
*   
*   Delaunay triangulation Algorithm and Mesh structure
*/
template<typename V, typename E, typename F, typename H>
class CDelaunayMesh : public CDynamicMesh<V,E,F,H>
{
public:
	typedef CBoundary<V,E,F,H> CBoundary;
	typedef CLoop<V,E,F,H> CLoop;
	typedef MeshVertexIterator<V,E,F,H>        MeshVertexIterator;
	typedef MeshFaceIterator<V,E,F,H>          MeshFaceIterator;
	typedef MeshEdgeIterator<V,E,F,H>          MeshEdgeIterator;
	typedef VertexVertexIterator<V,E,F,H>      VertexVertexIterator;
	typedef FaceVertexIterator<V,E,F,H>        FaceVertexIterator;
	typedef VertexEdgeIterator<V,E,F,H>        VertexEdgeIterator;
	typedef FaceEdgeIterator<V,E,F,H>          FaceEdgeIterator;
	typedef VertexFaceIterator<V,E,F,H>        VertexFaceIterator;
	typedef FaceHalfedgeIterator<V,E,F,H>      FaceHalfedgeIterator;
	typedef VertexOutHalfedgeIterator<V,E,F,H> VertexOutHalfedgeIterator;

public:
	/*! CDeluanayMesh constructor */
	CDelaunayMesh();
	/*! CDeluanayMesh destructor  */
	~CDelaunayMesh();
	
	/*!
	* process poly file
	* \param input input mesh file name
	* \param output output mesh file name
	*/
	void ProcessPoly( const char * input, const char * output );
	/*!
	* process poly file
	* \param input CPoly
	* \param output mesh
	*/
	void ProcessPoly( CPoly & poly, CDelaunayMesh<V,E,F,H> & omesh );

	/*!
	* Delaunay triangulate the input points
	* \param input input mesh file name
	* \param output output mesh file name
	 */
	void   DelaunayTriangulate( const char * input, const char * output );

	/*!
	* Delaunay triangulate the input points
	* \param input input points
	* \param output output mesh file name
	 */
	void   DelaunayTriangulate( std::vector<CPoint2> & points, CDelaunayMesh<V,E,F,H> & omesh );

protected:
	/*!	Construct initial triangle
	 *  \param grading false means the area is considered
	 */
    void   _initialize( bool grading = false );
	/*!	Insert one point to the current Delaunay triangulation
	 *  update the Delaunay Triangulation
	 *  \param pt the newly insterted point on the plane
	 *  \param pV the newly generated vertex
	 */
	virtual bool   _insert_vertex(CPoint2 & pt, V * & pV );
	/*!	Check if any segment is missing in the current Delaunay triangulation
	*   if so, split the segment to 2 halves, and insert the midpoint to the current DT (Delaunay Triangulation)
	*   update the DT.
	*/
	void   _insert_missing_segments();
	/*! Find the triangle with the worst quality 
     *	compute its Steiner point
	 *  if the Steiner point encroaches any segment, subdivide the segment
	 *  otherwise, insert the Steiner point to the current DT
	 *  \param lower_angle threshould for angle bound
	 */
	void   RemoveBadTriangles( double lowerAngle = 20 );
	/*! Find the triangle with the worst quality 
     *	compute its Steiner point
	 *	if the Steiner point doest encroach any segment
	 *  insert the Steiner point to the current DT.
	 *  if the Steiner point encroaches any segment, 
	 *  if the segment is not frozen then subdivide the segment, and continue.
	 *  if the segment is frozen, insert the Steinter point to the current DT.
	 *  \param lower_angle threshould for angle bound
	 */
	void   RemoveBadTrianglesWithFrozenSegments( double lowerAngle = 20 );
	/*!	Classify all the faces to be inside or outside.
	 *  a. all faces attaching to the initial bounding vertices are outside
	 *  b. the faces adjacent to the center of holes are outside
	 *  c. flood the outside label till reach the segments
	 */
	void   _classify_inside_outside();

public:
	/*!  Save inside faces to another mesh
	 *   \param output_mesh the output mesh
	 */
	void   Convert( CDelaunayMesh<V,E,F,H> & output_mesh );

protected:
	/*!	Insert a new vertex on an edge
	 *  \param e the edge to be split
	 *  \param pt the position of the new vertex
	 */
	virtual V*      _insert_vertex_on_edge( E * e, CPoint2 & pt );

	/*!	Split a segment and the edge attaching to it
	 *	\param pS the segment to be split
	 */
	void	_split_edge_on_segment( CSegment<V> * pS );

	/*!	Test if an edge is legal, if not, swap it
	 *  and test whether the new neighbors are legal
	 *  \param pV the new vertex
	 *  \param pE the edge to be tested
	 *  \param pF the face bound by pV and pE
	 *  \param level recursion level
	 */
	bool	_legalize_edge( const V * pV, E * pE, F * f, int level );

protected:
	/*! Locate a point in DT, find the face containing it
	*	\param pt the point to be located
	*   \param q  
	*/
	F *     _locate_point(  CPoint2 & pt, CPoint & q);
	/*! check if the circle through v0, v1 and centered at their midpoint contains v
	 *  \param v0, v1 the end vertices of an edge
	 *  \param v the point to be tested
	 */
	bool    _encroach( V* v0, V* v1, CPoint2 v);
	/*! Test if the circum circle of e includes the vertices of adjacent faces against e
	 *  \param e edge to be tested
	 */
	bool    _encroach( E *e );

	/*! Find an edge with end vertices v0, v1
	 * \param v0, v1 end vertices
	 */
	E *	          _find_edge( V* v0, V* v1);
	/*! Find a segment with end vertices v0, v1
	 *  \param v0, v1 end vertices
	 */
	CSegment<V> * _find_segment( V * v0, V * v1 );
	/*! create a new segment with end vertices v0, v1
	 *  \param v0, v1 end vertices
	 */
	CSegment<V> *  _create_segment( V* v0, V* v1 );


	/*!	Compute the statistics of triangle of face
	*   \param face input triangle
	*/
	void    _statistics( F * face );
	/*!	Compute the statistics of triangle of [v0,v1,v2]
	*   \param v0,v1,v2 input triangle
	*   \param minAngle, quality output minimal angle and quality of the triangle
	*/
	virtual void    _statistics( V* v0, V* v1, V* v2, double & minAngle, double & quality);
	/*!
	 *	freeze edges on segments, they can not be flipped any more
	 */
	void   _freeze_edges_on_segments();
	
protected:
	/*! oriented area bounded by face [a,b,c]
	* \param a,b,c points of a triangle
	*/
	double  __orient2d(  CPoint2 & a,  CPoint2 &b,  CPoint2 &c );
	/*! circumcenter of face
	 *	\param face input triangle
	 */
	CPoint2 __circumcenter( F* face );
	/*!	angle difference from -PI to PI
	 *	\param a1, a0 rotate a0 to a1
	 */
	double  __angle_difference( double a1, double a0 );
	/*! wheither the triangle [v0,v1,v2] contains bounding vertices
	 *  \param v0,v1,v2 triangle vertices 
	 */
	bool    __contains_bounding_vertices(V* v0, V* v1, V* v2);
	/*!	whether the corner angle at v1 in face [v0,v1,v2] is bounded by segments
	 *  \param v0,v1,v2 three vertices of a triangle
	 */
	bool	__is_segment_bounded_angle(V *v0, V *v1, V *v2 );

protected:
	/*! maximal number of vertices to be allowed in the meshing result */
	size_t m_max_number_vertices;
	/*! input vertices of PLSG */
	int    m_inputVertexNum;
	/*! the angle Lower bound */
	double m_angleLowerbound;
	/*! the quliaty lower bound */
    double m_qualityLowerbound;
	/*! starting face in point location procedure*/
	F *	   m_starting_face;
	/*! grading is false, consider the area; grading is true, consider the angle only*/
	bool   m_grading;
	/*! vector of segments */
	std::vector<CSegment<V>*>  m_segments;
	/*! vector of centers of holes */
	std::vector<CPoint2>       m_hole_centers;

};


//Compute the statistics of the face

template<typename V, typename E, typename F, typename H>
void CDelaunayMesh<V,E,F,H>::_statistics( F * face)
{
	int i = 0;
	V *v[3];
	for( FaceVertexIterator fviter( face ); !fviter.end(); fviter ++ )
	{
		v[i++] = *fviter;
	}
	_statistics( v[0], v[1], v[2], face->minAngle(), face->quality() );
	
	if( face->minAngle() < -1e-4 && face->inside() && ! __contains_bounding_vertices( v[0], v[1], v[2] ) ) 
	{
		//if this happends, the input segments may not be valid
		face->string() = "rgb=(0 0 1)";
		write_m("debug.m" );
	}

}

/* This function returns the minimum angle of a triangle (in degrees),
 * and the quality of the triangle, as required by the refinement algorithm */

template<typename V, typename E, typename F, typename H>
void CDelaunayMesh<V,E,F,H>::_statistics( V* v0, V* v1, V* v2, double & minAngle, double & quality)


{
        
    /* compute edge vectors */
    CPoint2 d[3];

	d[0] = v1->uv()-v0->uv();
	d[1] = v2->uv()-v1->uv();
	d[2] = v0->uv()-v2->uv();
    
    
    /* compute angle of each transformed edge vector */
	double t[3];
	for( int i  = 0; i < 3; i ++ )
	{
		t[i] = atan2(d[i][1], d[i][0]);
	}
    
    /* compute angle differences (angle at each vertex). This way of
     * calculating is not particularly good for small angles */
	double a[3];
	for( int i = 0; i < 3; i ++ )
	{
		a[i] = __angle_difference(t[(i+2)%3] + PI, t[i]);
	}
    
    /* compute minimum angle */
    double amin = (a[0] < a[1]) ? a[0] : a[1];
    if (a[2] < amin) amin = a[2];
    
    minAngle = amin * (360.0/TWOPI);   // return minimum angle

    /* fix for segment-bounded small angles */
    if (__is_segment_bounded_angle(v2, v0, v1))
        a[0] = (PI/3);    // assume a nice angle
    if (__is_segment_bounded_angle(v0, v1, v2))
        a[1] = (PI/3);    // assume a nice angle
    if (__is_segment_bounded_angle(v1, v2, v0))
        a[2] = (PI/3);    // assume a nice angle

    /* recompute minimum angle */
	amin = (a[0] < a[1]) ? a[0] : a[1];
    if (a[2] < amin) amin = a[2];

	if ( m_grading )
        quality = amin * (360.0/TWOPI);    // return triangle quality
    else {
        /* compute perimeter */
        double perimeter = d[0].norm() + d[1].norm() + d[2].norm();        
        quality = amin/(perimeter*perimeter);      // return triangle quality
    }

};


/*! \brief CDelaunayMesh destructor
*/

template<typename V, typename E, typename F, typename H>
CDelaunayMesh<V,E,F,H>::~CDelaunayMesh()
{
	for( size_t i = 0; i <  m_segments.size(); i ++ )
	{
		CSegment<V> * pS = m_segments[i];
		delete pS;
	}

}

//create a new segment

template<typename V, typename E, typename F, typename H>
CSegment<V> *  CDelaunayMesh<V,E,F,H>::_create_segment(V *v0, V *v1)
{
		CSegment<V> * pNS = new CSegment<V>();
		assert( pNS );

		m_segments.push_back( pNS );

		pNS->vertex(0) = v0;
		pNS->vertex(1) = v1;
				
		return pNS;
}

//find a segment

template<typename V, typename E, typename F, typename H>
CSegment<V> *  CDelaunayMesh<V,E,F,H>::_find_segment(V *v0, V *v1)
{    
	for ( size_t i = 0; i < m_segments.size() ; i++) {      // brute force!
        CSegment<V> * pS = m_segments[i];
		V * w0 = pS->vertex(0);
		V * w1 = pS->vertex(1);
        if (v0 == w0 && v1 == w1 || v0 == w1 && v1 == w0)
            return pS;
    }
    
    return NULL;
}

//oriented area of triangle [a,b,c]
template<typename V, typename E, typename F, typename H>
double CDelaunayMesh<V,E,F,H>::__orient2d(  CPoint2 & a,  CPoint2 &b,  CPoint2 &c )
{
	//adaptive exact arithmetic

	double pa[2], pb[2], pc[2];

	pa[0] = a[0]; pa[1] = a[1];
	pb[0] = b[0]; pb[1] = b[1];
	pc[0] = c[0]; pc[1] = c[1];

	return orient2d( pa, pb, pc);

}

/* returns the angle difference a1-a0 in the range -Pi..Pi */
template<typename V, typename E, typename F, typename H>
inline double CDelaunayMesh<V,E,F,H>::__angle_difference( double a1, double a0 )

{
    double d = a1 - a0;
    
    while (d > PI)   d -= TWOPI;
    while (d <= -PI) d += TWOPI;
    
    return d;
}

//locate a point in the current DT
template<typename V, typename E, typename F, typename H>
F* CDelaunayMesh<V,E,F,H>::_locate_point(  CPoint2 & pt, CPoint & q )
{

	assert( m_faces.size() > 0 );

	F * pF = m_starting_face;

	H * he[3];
	
	while( true )
	{
		he[0] = faceHalfedge( pF );
		he[1] = faceNextCcwHalfEdge( he[0] );
		he[2] = faceNextCcwHalfEdge( he[1] );

		double a[3];
		for( int i = 0; i < 3; i ++ )
		{
			a[i] = __orient2d( halfedgeSource( he[i])->uv(), halfedgeTarget( he[i] )->uv(), pt  );
		}

		if( a[0] >= 0 && a[1] >=0 && a[2] >= 0 )
		{
			 V* v[3];
			 CPoint p[3];

			 for( int i = 0; i < 3; i ++ )
			 {
				 v[i] = halfedgeSource( he[i]);
				 p[i] = v[i]->point();
			 }

			double s = __orient2d( v[0]->uv(), v[1]->uv(), v[2]->uv() );
			
			 q = CPoint(0, 0, 0);

			 for( int i = 0; i < 3; i ++ )
			 {
				 q = q + (p[i] * a[(i+1)%3]/s);
			 }


			return pF;
		}

		for( int i = 0; i < 3; i ++ )
		{
			if( a[i] < 0 )
			{
				H * pS = halfedgeSym( he[i] );
				if( pS == NULL ) return NULL;

				pF = halfedgeFace( pS );
				m_starting_face = pF;
				break;
			}
		}
	}
};



/*------------------------------------------------------------------------------------------------*/

template<typename V, typename E, typename F, typename H>
void CDelaunayMesh<V,E,F,H>::_initialize( bool grading )
{
	double M = VERYLARGE;

	V * v[3];

	std::vector<CPoint2> uvs;

	uvs.push_back( CPoint2( 3 * M,0) );
	uvs.push_back( CPoint2( 0, 3 * M) );
	uvs.push_back( CPoint2(-3 * M, -3 * M) );

	 /* read the vertices, starting at index #3 */
    for (int i = 0; i < 3; i++) 
	{
		v[i] = createVertex( ++ m_vertex_id );		
		v[i]->uv() = uvs[i];
		v[i]->point() = CPoint( uvs[i][0], uvs[i][1], 0 );
    }

	V * pV = v[0];
	v[0] = v[1]; v[1] = pV;

	F * pF = createFace( v, ++ m_face_id ); //the front face is the one with infinity point inside
	pF->inside() = false;
	_statistics( pF );

	pV = v[0];
	v[0] = v[1]; v[1] = pV;

	pF = createFace( v, ++ m_face_id );
	pF->inside() = true;
	_statistics( pF );

	for( FaceEdgeIterator feiter( pF ); !feiter.end(); feiter ++ )
	{
		E * pE = *feiter;
		pE->frozen() = true;
	}

	m_grading  = grading;
	m_inputVertexNum = 3;
	m_starting_face = pF;
};


/*---------------------------------------------------------------------------*/

template<typename V, typename E, typename F, typename H>
bool CDelaunayMesh<V,E,F,H>::__is_segment_bounded_angle(V *v0, V *v1, V *v2)

/* is "v1" the apex of a segment bounded angle? */

{
    return (v1->id() < m_inputVertexNum && _find_segment(v0, v1) != NULL && _find_segment(v1, v2) != NULL);
}


template<typename V, typename E, typename F, typename H>
bool CDelaunayMesh<V,E,F,H>::__contains_bounding_vertices(V* v0, V* v1, V* v2)
{
    return (v0->id() < 4) || (v1->id() < 4) || (v2->id() < 4);
}


/*---------------------------------------------------------------------------*/

template<typename V, typename E, typename F, typename H>
bool CDelaunayMesh<V,E,F,H>::_insert_vertex(CPoint2 & pt, V * & pV )
{
    
    /* find triangle "t0" containing "v" */
	CPoint q;

    F * pT = _locate_point( pt, q );

    if (pT == NULL)
	{
        printf("unable to locate a vertex");
		return false;
	}
	

	H * pH[3];

	pH[0] = faceHalfedge( pT );
	pH[1] = faceNextCcwHalfEdge( pH[0] );
	pH[2] = faceNextCcwHalfEdge( pH[1] );

	E * e[3];

	for( int i = 0; i < 3; i ++ )
	{
		e[i] = halfedgeEdge( pH[i] );
	}

	/* check if vertex is on some triangle edge */
    for (int i = 0; i < 3; i++) {

        V * v0 = halfedgeTarget( pH[i] );
        V * v1 = halfedgeSource( pH[i] );
		
        if (fabs( __orient2d(v0->uv(), v1->uv(), pt ) ) < 1e-8 ) 
		{
            pV = _insert_vertex_on_edge( halfedgeEdge( pH[i] ), pt );
            return true;
        }
    }
	bool inside = pT->inside();

	pV = splitFace( pT );
	pV->uv() = pt;
	pV->point() = q;


	for( VertexFaceIterator vfiter( pV ); !vfiter.end(); vfiter ++ )
	{
		F * f = *vfiter;
		_statistics( f );
		f->inside() = inside;
	}




	std::vector<E*> against_edges;
	std::vector<F*> against_faces;

	for( VertexOutHalfedgeIterator vhiter( this, pV ); !vhiter.end(); vhiter ++ )
	{
		H * pH = *vhiter;
		assert( halfedgeSource( pH ) == pV );
		against_edges.push_back( halfedgeEdge( faceNextCcwHalfEdge( pH ) ) );
		against_faces.push_back( halfedgeFace( pH ) );
	}

	for( size_t i = 0; i < against_edges.size(); i ++ )
	{
		_legalize_edge( pV, against_edges[i], against_faces[i], 0 );
	}

	return true;
}

/*---------------------------------------------------------------------------*/

template<typename V, typename E, typename F, typename H>
V* CDelaunayMesh<V,E,F,H>::_insert_vertex_on_edge( E * e, CPoint2 & pt )
{
	H * h0 = edgeHalfedge( e, 0 );
	H * h1 = edgeHalfedge( e, 1 );


	F * t0 = halfedgeFace( h0 );
	F * t1 = halfedgeFace( h1 );

    
    /* preserve inside flags and frozen state */
    bool inside0 = t0->inside();
    bool inside1 = t1->inside();
    bool frozen = e->frozen();

	E * e0 = halfedgeEdge( faceNextCcwHalfEdge( h0 ) );
	E * e1 = halfedgeEdge( faceNextClwHalfEdge( h0 ) );

	E * e2 = halfedgeEdge( faceNextCcwHalfEdge( h1 ) );
	E * e3 = halfedgeEdge( faceNextClwHalfEdge( h1 ) );

	V * v1 = edgeVertex1( e );
	V * v2 = edgeVertex2( e );

	V * pV = splitEdge( e );
	
	pV->uv() = pt;
	double d1 = (pt-v1->uv()).norm();
	double d2 = (pt-v2->uv()).norm();
	pV->point() = ( v1->point() * d2 + v2->point() * d1)/(d1+d2);

	
	F * t2 = halfedgeFace( halfedgeSym( faceNextClwHalfEdge( h0 ) ) );
	F * t3 = halfedgeFace( halfedgeSym( faceNextCcwHalfEdge( h1 ) ) ); 



    /* restore inside flags */
    t0->inside() = inside0;
    t2->inside() = inside0;
    t1->inside() = inside1;
    t3->inside() = inside1;
    
	for( VertexFaceIterator vfiter( pV ); !vfiter.end(); vfiter ++ )
	{
		F * pF = *vfiter;
		_statistics( pF );
	}


    /* set frozen states */
    e->frozen() = frozen;
    
	E * pE = halfedgeEdge( faceNextCcwHalfEdge( halfedgeSym( faceNextCcwHalfEdge( h1 ) ) ) );
	pE->frozen() = frozen;
    

	std::vector<E*> against_edges;
	std::vector<F*> against_faces;

	for( VertexOutHalfedgeIterator vhiter(this, pV ); !vhiter.end(); vhiter ++ )
	{
		H * pH = *vhiter;
		against_faces.push_back( halfedgeFace( pH ) );
		against_edges.push_back( halfedgeEdge( faceNextCcwHalfEdge(pH) ) );
	}

	for( size_t i =  0; i < against_faces.size(); i ++ )
	{
		F * pF = against_faces[i];
		E * pE = against_edges[i];

		if( !pF->inside() ) continue;

		_legalize_edge( pV, pE, pF, 0 );
	}

	return pV;
}

/*---------------------------------------------------------------------------*/

template<typename V, typename E, typename F, typename H>
bool CDelaunayMesh<V,E,F,H>::_encroach(V *v0, V *v1, CPoint2 pt)

/* Does vertex "v" encroach upon segment "v0-v1"? */

{
	CPoint2 p0 = v0->uv();
	CPoint2 p1 = v1->uv();

	double r = (p0-p1).norm2()/4.0;
	double d = (pt - (p0+p1)/2.0).norm2();
    
    return (d < r);
}

/*---------------------------------------------------------------------------*/
template<typename V, typename E, typename F, typename H>
bool CDelaunayMesh<V,E,F,H>::_encroach( E *e )

/* is edge "e" encroached? */
{
	H * h = edgeHalfedge( e, 0 );
	H * s = edgeHalfedge( e, 1 );

    V *v0 = halfedgeTarget( h );
	V *v1 = halfedgeTarget( s );
    
	V * v2 = halfedgeTarget( faceNextCcwHalfEdge( h )  );
	V * v3 = halfedgeTarget( faceNextCcwHalfEdge( s )  );

    return _encroach(v0, v1, v2->uv() )|| _encroach(v0, v1, v3->uv() );
}


/*---------------------------------------------------------------------------*/
template<typename V, typename E, typename F, typename H>
void CDelaunayMesh<V,E,F,H>::_insert_missing_segments()
{
	int id = 0;
    bool change;
    
    do {
        change = false;

		for( size_t i = 0; i < m_segments.size(); i ++ )
		{
			CSegment<V> * pS = m_segments[i];
			
			V * v0 = pS->vertex(0);
			V * v1 = pS->vertex(1);
			E * pE = _find_edge( v0, v1 );


			if( pE == NULL )
			{
				CPoint2 c = (v0->uv() + v1->uv())/2.0;
				V * pV = NULL;
				_insert_vertex( c, pV );
				
				pS->vertex(0)  = v0;
				pS->vertex(1)  = pV;
				CSegment<V> * pNS = _create_segment( pV, v1);

				change = true;
			}
			else if ( _encroach(pE)) 
			{
				_split_edge_on_segment( pS );
                change = true;
            }
		}
	} while (change);

}

/*---------------------------------------------------------------------------*/

template<typename V, typename E, typename F, typename H>
void CDelaunayMesh<V,E,F,H>::_freeze_edges_on_segments()
{
    
    for ( size_t i = 0; i < m_segments.size() ; i++) {
		CSegment<V> * pS = m_segments[i];
		
		V * v0 = pS->vertex(0);
		V * v1 = pS->vertex(1);        
        E * e = _find_edge(v0, v1);
        if (e != NULL) e->frozen() = true;
    }
}

/*---------------------------------------------------------------------------*/

template<typename V, typename E, typename F, typename H>
E * CDelaunayMesh<V,E,F,H>::_find_edge( V * v0, V * v1 )
{
	for( VertexOutHalfedgeIterator vhiter( this, v0 ); !vhiter.end(); vhiter ++ )
	{
		H * h = *vhiter;
		if ( halfedgeTarget( h ) == v1 )
		{
			return halfedgeEdge( h );
		}
	}
	return NULL;
}

/*---------------------------------------------------------------------------*/

template<typename V, typename E, typename F, typename H>
void CDelaunayMesh<V,E,F,H>::_split_edge_on_segment( CSegment<V> *s )
{
    V * v0 = s->vertex(0);
    V * v1 = s->vertex(1);
    
	E * pE = _find_edge( v0, v1 );

	CPoint2 pt = (v0->uv() + v1->uv())/2.0;
    V * v = _insert_vertex_on_edge( pE, pt );
	
	s->vertex(0) = v0;
	s->vertex(1) = v;
	
	CSegment<V> * pS = _create_segment( v, v1 );
}

/*---------------------------------------------------------------------------*/

template<typename V, typename E, typename F, typename H>
void CDelaunayMesh<V,E,F,H>::_classify_inside_outside()
{
	F * head = NULL;

	for( MeshFaceIterator fiter( this ); !fiter.end(); fiter ++ )
	{
		F * pF = *fiter;
		pF->inside() = true;
			
		for( FaceVertexIterator fviter( pF ); !fviter.end(); fviter ++ )
		{
			V * pV = *fviter;
			if( pV->id() <=3 )
			{
				pF->inside() = false;
				head = pF;
				continue;
			}
		}
	}

	for( MeshFaceIterator fiter( this ); !fiter.end(); fiter ++ )
	{
		F * pF = *fiter;
		pF->touched() = false;
	}

	std::queue<F*> queue;
	queue.push( head );
	head->touched() = true;

	while( !queue.empty() )
	{
		F * pF = queue.front();
		queue.pop();

		for( FaceHalfedgeIterator fhiter( pF ); !fhiter.end(); fhiter ++ )
		{
			H * pH = *fhiter;
			E * pE = halfedgeEdge( pH );
			if( pE->frozen() ) continue;

			H * pS = halfedgeSym( pH );
			
			F * pSF = halfedgeFace( pS );
			if( pSF->touched() ) continue;
			pSF->inside() = false;
			pSF->touched() = true;
			queue.push( pSF );
		}
	}


	for( size_t i = 0 ; i <  m_hole_centers.size(); i ++ )
	{
		CPoint2 center = m_hole_centers[i];

		for( MeshFaceIterator fiter( this ); !fiter.end(); fiter ++ )
		{
			F * pF = *fiter;
			pF->touched() = false;
		}
		
		CPoint q;
		head = _locate_point( center, q );
		head->inside() = false;

		std::queue<F*> queue;
		queue.push( head );
		head->touched() = true;

		while( !queue.empty() )
		{
			F * pF = queue.front();
			queue.pop();

			for( FaceHalfedgeIterator fhiter( pF ); !fhiter.end(); fhiter ++ )
			{
				H * pH = *fhiter;
				E * pE = halfedgeEdge( pH );
				if( pE->frozen() ) continue;

				H * pS = halfedgeSym( pH );
				
				F * pSF = halfedgeFace( pS );
				if( pSF->touched() ) continue;
				pSF->inside() = head->inside();
				pSF->touched() = true;
				queue.push( pSF );
			}
		}
	}
}

/*---------------------------------------------------------------------------*/

template<typename V, typename E, typename F, typename H>
void CDelaunayMesh<V,E,F,H>::RemoveBadTriangles( double angle )
{
   double m_angleLowerbound = angle;
   double worstQuality;
    F *worstTriangle;
    bool encroach;
    
	int id = 0;

   
    while ( m_verts.size() < m_max_number_vertices ) {
        /* find the worst triangle */
        worstQuality = 1e30;

        for ( MeshFaceIterator fiter( this ); !fiter.end(); fiter ++ ) {
			F * f = *fiter;
			if( !f->inside() ) continue;

			if( m_grading ) //only consider angles, no area
			{
				if (f->minAngle() < worstQuality) {
					worstQuality = f->minAngle();
					worstTriangle = f;
				}
            }
			else if ( f->quality() < worstQuality) {
                worstQuality = f->quality();
                worstTriangle = f;
            }
        }

		printf("Worst Quality is %f in face %d\n", worstQuality, worstTriangle->id() );
		
		/* check if we're done */
        if ( m_grading &&  worstQuality >= m_angleLowerbound  ) return;
        if ( !m_grading && worstQuality >= m_qualityLowerbound) return;
        
        /* compute worst triangle circumcenter */
        CPoint2 c = __circumcenter(worstTriangle);
        
        /* does it encroach upon any segment? */
		encroach = false;
		
		for( size_t i  = 0; i < m_segments.size() ; i ++ )
		{
			CSegment<V> * pS = m_segments[i];
			V * v0 = pS->vertex(0);
			V * v1 = pS->vertex(1);

			if( _encroach( v0, v1, c ) )
			{
				_split_edge_on_segment( pS );
				encroach = true;
			}
		}
        
        if (!encroach) {        // did not encroach
			V * pV = NULL;
            _insert_vertex( c, pV );
		}
	}
}
/*---------------------------------------------------------------------------*/

template<typename V, typename E, typename F, typename H>
void CDelaunayMesh<V,E,F,H>::RemoveBadTrianglesWithFrozenSegments( double angle )
{
   double m_angleLowerbound = angle;
   double worstQuality;
    F *worstTriangle;
    bool encroach;
    
	int id = 0;

   
    while ( m_verts.size() < m_max_number_vertices ) {
        /* find the worst triangle */
        worstQuality = 1e30;

        for ( MeshFaceIterator fiter( this ); !fiter.end(); fiter ++ ) {
			F * f = *fiter;
			if( !f->inside() ) continue;

			if( m_grading ) //only consider angles, no area
			{
				if (f->minAngle() < worstQuality) {
					worstQuality = f->minAngle();
					worstTriangle = f;
				}
            }
			else if ( f->quality() < worstQuality) {
                worstQuality = f->quality();
                worstTriangle = f;
            }
        }

		printf("Worst Quality is %f in face %d\n", worstQuality, worstTriangle->id() );
		
		/* check if we're done */
        if ( m_grading &&  worstQuality >= m_angleLowerbound  ) return;
        if ( !m_grading && worstQuality >= m_qualityLowerbound) return;
        
        /* compute worst triangle circumcenter */
        CPoint2 c = __circumcenter(worstTriangle);
        
        /* does it encroach upon any segment? */
		encroach = false;
		
		for( size_t i  = 0; i < m_segments.size() ; i ++ )
		{
			CSegment<V> * pS = m_segments[i];
			V * v0 = pS->vertex(0);
			V * v1 = pS->vertex(1);

			/* for frozen edges, do not split, 
			   add Steiner point directly
			 */

			E * e = _find_edge(v0, v1);
			if (e != NULL &&  e->frozen())
			{
				encroach = false;
			}
			else 
			if( _encroach( v0, v1, c ) )
			{
				_split_edge_on_segment( pS );
				encroach = true;
			}
		}
        
        if (!encroach) {        // did not encroach
			V * pV = NULL;
            _insert_vertex( c, pV );
		}
	}
}

/*---------------------------------------------------------------------------*/

template<typename V, typename E, typename F, typename H>
CPoint2 CDelaunayMesh<V,E,F,H>::__circumcenter( F* face )
{

	std::vector<CPoint2> vs;
	for( FaceVertexIterator fviter( face ); !fviter.end(); fviter ++ )
	{
		V* pV = *fviter;
		vs.push_back( pV->uv() );
	}

	CPoint2 a = vs[0];
	CPoint2 b = vs[1];
	CPoint2 c = vs[2];


	double xba, yba, xca, yca;
	double balength, calength;
	double denominator;
	double xcirca, ycirca;

	/* Use coordinates relative to point `a' of the triangle. */
	xba = b[0] - a[0];
	yba = b[1] - a[1];
	xca = c[0] - a[0];
	yca = c[1] - a[1];
	/* Squares of lengths of the edges incident to `a'. */
	balength = xba * xba + yba * yba;
	calength = xca * xca + yca * yca;

	/* Take your chances with floating-point roundoff. */
	denominator = 0.5 / (xba * yca - yba * xca);

	/* Calculate offset (from `a') of circumcenter. */
	xcirca = (yca * balength - yba * calength) * denominator;  
	ycirca = (xba * calength - xca * balength) * denominator;  

	return CPoint2( xcirca + a[0], ycirca + a[1] );
};



/*! \brief Convert Delaunay Mesh to common mesh
 *  \param output_mesh the output mesh
 *
 *  Remove all the faces to the exterior vertices
 */
	
template<typename V, typename E, typename F, typename H>
void CDelaunayMesh<V,E,F,H>::Convert( CDelaunayMesh<V,E,F,H> & output_mesh )
{
	for( MeshVertexIterator viter( this ); !viter.end(); viter ++ )
	{
		V * pV = *viter;
		pV->touched() = false;
		for( VertexFaceIterator vfiter( pV ); !vfiter.end(); vfiter ++ )
		{
			F * pF = *vfiter;
			if( pF->inside() ){
				pV->touched() = true;
				break;
			}
		}
	}
	
	int vid = 1;

	for( MeshVertexIterator viter( this ); !viter.end(); viter ++ )
	{
		V * pV = *viter;
		if( !pV->touched() ) continue;
		V * pW = output_mesh.createVertex( pV->id() );
		pW->point() = pV->point();
		pW->uv() = pV->uv();
		pW->father() = pV->father();
	}
			
	for( MeshFaceIterator fiter( this ); !fiter.end(); fiter ++ )
	{
		F* pF = *fiter;
		if( !pF->inside() ) continue;

		V * v[3];
		int i = 0;
		for( FaceVertexIterator fviter( pF ); !fviter.end(); fviter ++ )
		{
			V * pV = *fviter;
			V * pW = output_mesh.idVertex( pV->id() );
			assert( pW );
			v[i++] = pW;
		}
		F * pWF = output_mesh.createFace( v, pF->id() );
		pWF->inside() = pF->inside();
	}

	for( MeshEdgeIterator eiter( this );  !eiter.end(); eiter ++ )
	{
		E * pE = *eiter;
		if( !pE->frozen() ) continue;

		V * v1 = edgeVertex1( pE );
		V * v2 = edgeVertex2( pE );
		if( !v1->touched() || !v2->touched() ) continue;

		V * w1 = output_mesh.idVertex( v1->id() );
		V * w2 = output_mesh.idVertex( v2->id() );
		E * pWE = output_mesh.vertexEdge( w1, w2 );
		if( pWE )
			pWE->frozen() = true;
	}

	output_mesh.labelBoundary();
};



/*---------------------------------------------------------------------------*/

//input poly file, mesh it

template<typename V, typename E, typename F, typename H>
void CDelaunayMesh<V,E,F,H>::ProcessPoly( const char * input, const char * output )
{

	_initialize( false ); //grading false = consider area; grading true, don't care area
	
	CPoly poly(input);
	std::vector<V*> ws;
	for( size_t i = 0; i < poly.points().size(); i ++ )
	{ 
		V * pV = NULL;
		_insert_vertex( poly.points()[i], pV );
		ws.push_back( pV );
		m_inputVertexNum ++;
	}

	for( size_t i = 0; i < poly.segments().size(); i ++ )
	{
		CPolySegment & s = poly.segments()[i];
		_create_segment( ws[s.vertex(0)-1], ws[s.vertex(1)-1]);
	}

	for( size_t i = 0; i < poly.holes_centers().size(); i ++ )
	{
		m_hole_centers.push_back( poly.holes_centers()[i] );
	}

	//insert points
	_insert_missing_segments();
	_freeze_edges_on_segments();	
	_classify_inside_outside();
	
	//m_qualityLowerbound = 180;
	m_qualityLowerbound = 30;

	RemoveBadTriangles();

	CDelaunayMesh<V,E,F,H> omesh;
	Convert( omesh );
 	omesh.write_m( output );
};	
//input poly file, mesh it

template<typename V, typename E, typename F, typename H>
void CDelaunayMesh<V,E,F,H>::ProcessPoly( CPoly & poly, CDelaunayMesh<V,E,F,H> & omesh )
{

	_initialize( false ); //grading false = consider area; grading true, don't care area
	
	std::vector<V*> ws;
	for( size_t i = 0; i < poly.points().size(); i ++ )
	{ 
		V * pV = NULL;
		_insert_vertex( poly.points()[i], pV );
		ws.push_back( pV );
		m_inputVertexNum ++;
	}

	for( size_t i = 0; i < poly.segments().size(); i ++ )
	{
		CPolySegment & s = poly.segments()[i];
		_create_segment( ws[s.vertex(0)-1], ws[s.vertex(1)-1]);
	}

	for( size_t i = 0; i < poly.holes_centers().size(); i ++ )
	{
		m_hole_centers.push_back( poly.holes_centers()[i] );
	}

	//insert points
	_insert_missing_segments();
	_freeze_edges_on_segments();	
	_classify_inside_outside();
	
	m_qualityLowerbound = 30.0;
	//RemoveBadTrianglesWithFrozenSegments( 20 );
	RemoveBadTrianglesWithFrozenSegments( 30 );

	Convert( omesh );
};	

/*---------------------------------------------------------------------------*/

//pV is the new generated vertex, pH is the halfedge in pF, against pV

template<typename V, typename E, typename F, typename H>
bool CDelaunayMesh<V,E,F,H>::_legalize_edge( const V * pV, E * pE, F * pF, int level )
{

	H * pH = edgeHalfedge( pE, 0 );
	H * pS = edgeHalfedge( pE, 1 );

	if( pV != halfedgeTarget( faceNextCcwHalfEdge( pH ) ) )
	{
		assert( pV == halfedgeTarget( faceNextCcwHalfEdge( pS ) ) );
		H * tmp = pH; pH = pS; pS = tmp;
	}


    V *v0, *v1, *v2, *v3;
    E *e0, *e1, *e2, *e3;
    F *t0, *t1;
    int currentCount;
    double currentMinAngle;
    double a0, a1, q0, q1;
    int proposedCount;
    double proposedMinAngle;
    
    if ( pE->frozen() )              // we don't flip frozen edges
	{
        return false;
	}

    t0 = halfedgeFace( pH );
    t1 = halfedgeFace( pS );


    if (t0->inside() != t1->inside() )
        printf("inconsistency between inside and outside");
    

	v0 = halfedgeSource( pH );
	v2 = halfedgeTarget( pH );

	pH = faceNextCcwHalfEdge( pH );
	e2 = halfedgeEdge( pH );
	v3 = halfedgeTarget( pH );
	assert( v3 == pV );
	pH = faceNextCcwHalfEdge( pH );
	e3 = halfedgeEdge( pH );

	pS = faceNextCcwHalfEdge( pS );
	v1 = halfedgeTarget( pS );
	e0 = halfedgeEdge( pS );
	pS = faceNextCcwHalfEdge( pS );
	e1 = halfedgeEdge( pS );

    
    /* if edge cannot be flipped, do not consider any further */
    if ( __orient2d(v1->uv(), v2->uv(), v3->uv() ) <= 0 || __orient2d(v1->uv() , v3->uv(), v0->uv() ) <= 0)
	{
        return false;
	}

    /* compute the current count and minimum angle */
    currentCount = 0;
    if ( __contains_bounding_vertices(v0, v2, v3)) currentCount++;
    if ( __contains_bounding_vertices(v2, v0, v1)) currentCount++;

	
    a0 = t0->minAngle();
	a1 = t1->minAngle();
    currentMinAngle = (a0 < a1) ? a0 : a1;
    
    /* compute the proposed count and minimum angle */
    proposedCount = 0;
    if (__contains_bounding_vertices(v1, v3, v0)) proposedCount++;
    if (__contains_bounding_vertices(v3, v1, v2)) proposedCount++;
    
	_statistics(v1, v3, v0, a0, q0);
    _statistics(v3, v1, v2, a1, q1);
    
	proposedMinAngle = (a0 < a1) ? a0 : a1;
    
    /* flip if this locally improves the count or minimum angle */
    if (proposedCount < currentCount || ( proposedCount == currentCount && proposedMinAngle > currentMinAngle)) 
	{
		swapEdge( pE );
		
		F * f0 = edgeFace1( pE );
		_statistics( f0 );
		F * f1 = edgeFace2( pE );
		_statistics( f1 );


		pH = edgeHalfedge( pE, 0 );
		pS = edgeHalfedge( pE, 1 );

		if( halfedgeSource( pH ) != pV )
		{
			assert( halfedgeSource( pS ) == pV );
			H * th = pH; pH = pS; pS = th;
		}
		
		E * pHE = halfedgeEdge( faceNextCcwHalfEdge( pH ) );
		E * pSE = halfedgeEdge( faceNextClwHalfEdge( pS ) );

		F * pHF = halfedgeFace( pH );
		F * pSF = halfedgeFace( pS );

		_legalize_edge( pV, pHE, pHF, level + 1 );
		_legalize_edge( pV, pSE, pSF, level + 1 );

		return true;
    }
	return false;

}

/*---------------------------------------------------------------------------*/
/*! \brief CDelaunayMesh constructor
*/

template<typename V, typename E, typename F, typename H>
CDelaunayMesh<V,E,F,H>::CDelaunayMesh()
{

	m_angleLowerbound   = 25; //30 degree
    m_qualityLowerbound = 15;//15.0;
	m_max_number_vertices = 8192;

}

//input poly file, compute the Delaunay triangulation

template<typename V, typename E, typename F, typename H>
void CDelaunayMesh<V,E,F,H>::DelaunayTriangulate( const char * input, const char * output )
{

	_initialize( false ); //grading false = consider area; grading true, don't care area
	
	CPoly poly(input);
	std::vector<V*> ws;
	for( size_t i = 0; i < poly.points().size(); i ++ )
	{ 
		V * pV = NULL;
		_insert_vertex( poly.points()[i], pV );
		ws.push_back( pV );
		m_inputVertexNum ++;
	}

	for( size_t i = 0; i < poly.segments().size(); i ++ )
	{
		CPolySegment & s = poly.segments()[i];
		_create_segment( ws[s.vertex(0)-1], ws[s.vertex(1)-1]);
	}

	for( size_t i = 0; i < poly.holes_centers().size(); i ++ )
	{
		m_hole_centers.push_back( poly.holes_centers()[i] );
	}

	//insert points
	_insert_missing_segments();
	_freeze_edges_on_segments();	
	_classify_inside_outside();

	CDelaunayMesh<V,E,F,H> omesh;
	Convert( omesh );
 	omesh.write_m( output );
};	

//input points

template<typename V, typename E, typename F, typename H>
void CDelaunayMesh<V,E,F,H>::DelaunayTriangulate( std::vector<CPoint2> & points, CDelaunayMesh<V,E,F,H> & omesh  )
{

	_initialize( false ); //grading false = consider area; grading true, don't care area
	
	for( size_t i = 0; i < points.size(); i ++ )
	{ 
		V * pV = NULL;
		_insert_vertex( points[i], pV );
		m_inputVertexNum ++;
	}

	for( MeshFaceIterator fiter( this ); !fiter.end(); fiter ++ )
	{
		F * pF = *fiter;
		pF->inside() = true;
			
		for( FaceVertexIterator fviter( pF ); !fviter.end(); fviter ++ )
		{
			V * pV = *fviter;
			if( pV->id() <=3 )
			{
				pF->inside() = false;
				continue;
			}
		}
	}

	Convert( omesh );
	//omesh.write_m("test.m");
};	


typedef CDelaunayMesh<CDelaunayVertex, CDelaunayEdge, CDelaunayFace, CDelaunayHalfEdge> CDTMesh;


}
#endif  _DELAUNAY_MESH_H_