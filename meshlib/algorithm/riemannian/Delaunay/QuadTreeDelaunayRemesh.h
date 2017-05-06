/*! \file  QuadTreeDelaunayRemesh.h
*   \brief Algorithm: Remesh a given mesh using Ruppert's Delaunay refinement algorithm and Using QuadTree 
*   \date  Documented on 10/15/2010
*
*   Adding Remesh a given mesh using Ruppert's Delaunay refinement algorithm. The input mesh must not have
*   flipped faces. In practice, we prefer to use harmonic mapping result as the input mesh, instead of Riemann mapping result.
*/

#ifndef  _QUAD_TREE_DELAUNAY_REMESH_H_
#define  _QUAD_TREE_DELAUNAY_REMESH_H_

#include <map>
#include <vector>
#include <queue>

#include "DelaunayMesh.h"
#include "Geometry/QuadTree.h"

namespace MeshLib
{

/*!	\brief CQTDelaunayRemeshVertex class
	 *  
	 *	Vertex class for Remeshing using Delaunay triangulation
	 *	Traits, 
	 *  vertex uv coordinates m_uv
	 *  vertex rgb, m_rgb
	 *  vertex conformal factor, m_lambda
	 *  vertex father, m_father
	 *  static with rgb trait m_with_rgb
	 */
  class CQTDelaunayRemeshVertex : public  CVertex
  {
  protected:
    CPoint2  m_uv;
	CPoint   m_rgb;
	int      m_father;
	bool     m_touched;
	double   m_k;
	double   m_lambda;
	static bool m_with_rgb;
  public:
	  
	  int     & father()  { return m_father; };
	  CPoint2 & uv()      { return m_uv; };
	  bool    & touched() { return m_touched; };
	  double  & k()       { return m_k;       };
	  CPoint  & rgb()     { return m_rgb;     };
	  double  & lambda()  { return m_lambda;  };

	  void	  _from_string();
	  void    _to_string();
  };


/*! Read vertex uv, rgb traits from string */

 inline void CQTDelaunayRemeshVertex::_from_string()
 {
		  CParser parser( m_string );
		
		  for( std::list<CToken*>::iterator iter = parser.tokens().begin() ; iter != parser.tokens().end(); ++ iter )
		  {
			  CToken * token = *iter;
			  if( token->m_key == "uv" )
			  {
				  token->m_value >> m_uv;
			  }
			  if( token->m_key == "rgb" )
			  {
				  token->m_value >> m_rgb;
				  m_with_rgb = true;
			  }
		  }
  };
/*! write vertex uv, rgb traits to string */

  inline void CQTDelaunayRemeshVertex::_to_string()
  {
		  CParser parser( m_string );
		  parser._removeToken( "uv" );
		  parser._removeToken( "rgb" );
			
		  parser._toString( m_string );

		  std::stringstream iss;
		  iss << "uv=(" << m_uv[0] << " " << m_uv[1] << ")";

		  if( m_string.length() > 0 )
			  m_string += " ";
		  m_string += iss.str();

		  if( m_with_rgb )
		  {
			  std::stringstream iss;
			  iss << "rgb=(" << m_rgb[0] << " " << m_rgb[1] << " " << m_rgb[2] << ")";

			  if( m_string.length() > 0 )
				  m_string += " ";
			  m_string += iss.str();
		}
  };

  /*!	\brief CQTDelaunayRemeshEdge class
	 *  
	 *	Edge class for Remeshing using Delaunay triangulation
	 *	Traits, 
	 *  frozenn, dihedral angle
	 */

 class CQTDelaunayRemeshEdge : public  CEdge
 {
 public:
	 CQTDelaunayRemeshEdge() { m_frozen = false; };
	 bool & frozen() { return m_frozen; };
	 void _to_string();
	 double & angle() { return m_angle; };

 protected:

	 bool m_frozen;
	 double m_angle;	//dihedral angle
 };

 inline void CQTDelaunayRemeshEdge::_to_string()
 {
	 if( m_frozen ) m_string="sharp";
 }
 /*!	\brief CQTDelaunayRemeshHalfEdge class
	 *  
	 *	HalfEdge class for Remeshing using Delaunay triangulation
	 *	Traits, 
	 *  angle
	 */

 class CQTDelaunayRemeshHalfEdge : public CHalfEdge
 {
 public:
	 double & angle() { return m_angle; };
 public:
	 double m_angle;
 };

  /*!	\brief CQTDelaunayRemeshFace class
	 *  
	 *	Face class for Remeshing using Delaunay triangulation
	 *	Traits, area, uv_area, minimal angle, quality
	 */

class CQTDelaunayRemeshFace : public CFace
 {
 public:
	 
	 CQTDelaunayRemeshFace() { m_inside = true;  }; //by default all faces are inside

	 bool   & inside()   { return m_inside;   };
	 double & minAngle() { return m_minAngle; };
	 double & quality()  { return m_quality;  };
	 double & area()     { return m_area;     };
	 double & uv_area()  { return m_uv_area;  };
	 double & weight()   { return m_weight;   };

	 void   _to_string();
	 bool   & touched() { return m_touched; };
 
	/*! Intersection Testing */

	/* Check if the bounding box of current face p is lower-left corner, q is upper right corner
		intersect the current triangle
	*/
	bool intersect(CPoint2 p, CPoint2 q)
	{
		
		CPoint2 max(-1e+10,-1e+10 );
		CPoint2 min( 1e+10, 1e+10 );

		CHalfEdge * ph = halfedge();
		for( int i = 0; i < 3; i ++, ph = ph->he_next() )
		{
			CQTDelaunayRemeshVertex * pV = (CQTDelaunayRemeshVertex * ) ph->target();
			CPoint2 tq = pV->uv();
			
			for( int k = 0; k < 2; k ++ )
			{
				max[k] = ( max[k] > tq[k])?max[k]:tq[k];
				min[k] = ( min[k] < tq[k])?min[k]:tq[k];
			}
		}

		if( max[0] < p[0] || max[1] < p[1] || min[0] > q[0] || min[1] > q[1] ) return false;
		return true;
	}
	/*!
	 *	Verify if the point pt is inside the triangle
	 */
	 bool  inside( CPoint2 pt )   
	 { 
		CPoint2 q[3];
		CHalfEdge * ph = halfedge();
		for( int i = 0; i < 3; i ++, ph = ph->he_next() )
		{
			CQTDelaunayRemeshVertex * pV = (CQTDelaunayRemeshVertex * ) ph->target();
			q[i] = pV->uv();	
		}

		double s = (q[1]-q[0])^(q[2]-q[0]);
		CPoint bary;

		for( int i = 0; i < 3; i ++ )
		{
			bary[i] = (q[(i+1)%3]-pt)^(q[(i+2)%3]-pt)/s;
		}

		if( bary[0] >= 0 && bary[1] >=0 && bary[2] >= 0 ) return true;
		return false;
	 };

 protected:
	 
	 bool m_touched;
	 bool m_inside;
	 double m_minAngle;
	 double m_quality;
	 double m_area;
	 double m_uv_area;
	 double m_weight;

 };

inline void CQTDelaunayRemeshFace::_to_string()
{
	if( !m_inside )
	{
		m_string = "rgb=(1 0 0)";
	}
};

/*!	\brief CQTSample class
 *
 *	sample of vertices
 */
class CQTSample
{ 
public:
	CQTSample() { m_k = 0; m_lambda = 0; };

	CPoint m_pos;
	CPoint m_rgb;
	double m_lambda;
	double m_k;

	CQTSample operator+( CQTSample & s2 )
	{
		CQTSample r;
		r.m_pos = m_pos + s2.m_pos;
		r.m_rgb = m_rgb + s2.m_rgb;
		r.m_lambda = m_lambda + s2.m_lambda;
		r.m_k = m_k + s2.m_k;
		return r;
	};

	CQTSample operator*( double d )
	{
		CQTSample r;
		r.m_pos = m_pos * d;
		r.m_rgb = m_rgb * d;
		r.m_lambda = m_lambda * d;
		r.m_k = m_k * d;
		return r;
	};
};


/*! \brief CQuadTreeDeluanayRemesh class
*   
*   Remesh a given surface using Ruppert's Delaunay refinement algorithm.
*/
template<typename V, typename E, typename F, typename H>
class CQuadTreeDelaunayRemesh : public CDelaunayMesh<V,E,F,H>
{
public:
	
	/*!
	* Remesh a input mesh with uv coordinates
	* \param input_mesh input mesh file name
	*/
	void Remesh( CQuadTreeDelaunayRemesh<V,E,F,H> & input_mesh );

	/*! \brief Compute the conformal factor
	*/
	void _lambda();

	void Convert( CQuadTreeDelaunayRemesh<V,E,F,H> & output_mesh );
	/*! Locate the point, interpolate all the traits
	* \param p point
	* \param q interpolated position 
	* \param rgb interpolated color
	*/
	F * _locate_point(  CPoint2 & pt, CPoint & q, CPoint & rgb );

	/*! 
	 *	build the quad tree
	 */
	void _construct_tree( int n);

protected:
	/*! Locate the point, interpolate all the traits
	* \param input point
	* \param sample interlated sample value
	*/
	F * _locate_point(  CPoint2 & pt,  CQTSample & sample );
	/*!	background mesh to be sampled.
	 *
	 */
	CQuadTreeDelaunayRemesh<V,E,F,H> * m_background;

	/*! invert a point to the current Delaunay triangulation on an edge
	*  \param e the edge 
	*  \param pt the point
	*/

	V* _insert_vertex_on_edge( E * e, CPoint2 & pt );
	/*! invert a point to the current Delaunay triangulation 
	 * \param pt the input point
	 * \param pV the newly generated vertex
	*/
	bool _insert_vertex(CPoint2 & pt, V * & pV );
	/*! Compute the statistics of a triangle 
	 * \param v0,v1,v2 form the triangle
	 * \param minAngle minimal angle
	 * \param quality the overall quaity
	*/
	void _statistics( V* v0, V* v1, V* v2, double & minAngle, double & quality);

	/*!	Quadtree for the purpose of locating a point
	 *
	 */
	CQuadTree<F> m_qtree;

};

//compute conformal factor
template<typename V, typename E, typename F, typename H>
void CQuadTreeDelaunayRemesh<V,E,F,H>::_lambda()
{
	double total_area = 0;
	double total_uv_area = 0;

	for( MeshFaceIterator fiter( this ); !fiter.end(); fiter ++ )
	{
		F * pF = *fiter;
		std::vector<CPoint> ps;
		std::vector<CPoint2> uvs;
		for( FaceVertexIterator fviter( pF ); !fviter.end(); fviter ++ )
		{
			V * pV = *fviter;
			ps.push_back( pV->point() );
			uvs.push_back( pV->uv() );
		}

		CPoint s = (ps[1]-ps[0])^(ps[2]-ps[0]);
		pF->area() = s.norm()/2.0;
		if( pF->area() < 1e-12 )
		{
			std::cerr << "Warning: Face area is too small" << std::endl;
			pF->area() = 1e-12;
		}
		total_area += pF->area();

		double a= (uvs[1]-uvs[0])^(uvs[2]-uvs[0]);
		pF->uv_area() = a/2.0;
		if( pF->uv_area() < 1e-12 )
		{
			std::cerr << "Warning: UV_area is too small" << std::endl;
			pF->uv_area() = 1e-12;
		}
		total_uv_area += pF->uv_area();
	}
	
	for( MeshFaceIterator fiter( this ); !fiter.end(); fiter ++ )
	{
		F * pF = *fiter;
		pF->area() /= total_area;
		pF->uv_area() /= total_uv_area;
	}
	

	for( MeshVertexIterator viter( this ); !viter.end(); viter ++ )
	{
		V * pV = *viter;
		double area = 0;
		double uv_area = 0;

		for( VertexFaceIterator vfiter( pV ); !vfiter.end(); vfiter ++ )
		{
			 F * pF = *vfiter;
			 area += pF->area();
			 uv_area += pF->uv_area();
		}
		
		pV->lambda() = area/uv_area;
	}
	
};	


//refine the mesh triangulation

template<typename V, typename E, typename F, typename H>
void CQuadTreeDelaunayRemesh<V,E,F,H>::Remesh( CQuadTreeDelaunayRemesh<V,E,F,H> & mesh )
{
	int N = 256;
	double r = 0.99;
	//double r = 0.98;
	//double r = 0.95;
	m_qualityLowerbound = 300; //for K 0.6; //128 for r3 area, 300 for uv_area only
	m_max_number_vertices = 8192*4;
	m_background = &mesh;
	m_starting_face = NULL;
	m_background->m_starting_face = NULL;

	_initialize( false ); //grading false = consider area; grading true, don't care area
	

	std::vector<V*> ws;

	for( int i = 0; i < N; i ++ )
	{
		double alpha = TWOPI * i/(double)N;
		CPoint2 uv( r * cos(alpha), r * sin( alpha) );
		V * pV = NULL;
		_insert_vertex( uv, pV );
		ws.push_back( pV );
	}

	for( int i = 0; i < N; i ++ )
	{
		_create_segment( ws[i], ws[(i+1)%N]);
	}

	//insert points
	_insert_missing_segments();
	_freeze_edges_on_segments();	
	_classify_inside_outside();

	RemoveBadTriangles();
}

/*---------------------------------------------------------------------------*/

template<typename V, typename E, typename F, typename H>
F* CQuadTreeDelaunayRemesh<V,E,F,H>::_locate_point(  CPoint2 & pt, CQTSample & sample )
{

	assert( m_faces.size() > 0 );
	F * pF = m_qtree._locate( pt );

	if( !pF ) return NULL;
	
	H * he[3];

	he[0] = faceHalfedge( pF );
	he[1] = faceNextCcwHalfEdge( he[0] );
	he[2] = faceNextCcwHalfEdge( he[1] );

	double a[3];
	for( int i = 0; i < 3; i ++ )
	{
		a[i] = __orient2d( halfedgeSource( he[i])->uv(), halfedgeTarget( he[i] )->uv(), pt  );
	}
	double s = __orient2d( halfedgeSource( he[0] )->uv(), halfedgeSource( he[1] )->uv(), halfedgeSource( he[2] )->uv());

	 V* v[3];
	 CQTSample samples[3];

	 for( int i = 0; i < 3; i ++ )
	 {
		 v[i] = halfedgeSource( he[i]);
		 samples[i].m_rgb = v[i]->rgb();
		 samples[i].m_pos = v[i]->point();
		 samples[i].m_k	  = v[i]->k();
		 samples[i].m_lambda = v[i]->lambda();
	 }

	s = __orient2d( v[0]->uv(), v[1]->uv(), v[2]->uv() );
	if( s < 1e-8 )
	{
		//q=CPoint(1,0,0);
		sample = samples[0];
	}
	else
	{
		for( int i = 0; i < 3; i ++ )
		{
			sample = sample + samples[i] * (a[(i+1)%3]/s);
		}
	}

	return pF;
};


/*---------------------------------------------------------------------------*/

template<typename V, typename E, typename F, typename H>
F* CQuadTreeDelaunayRemesh<V,E,F,H>::_locate_point(  CPoint2 & pt, CPoint & q, CPoint & rgb )
{

	assert( m_faces.size() > 0 );

	if( m_starting_face == NULL )
	{
		double max_area = -1;
		for( MeshFaceIterator ffiter( this ); !ffiter.end(); ffiter ++ )
		{
			F * pF = *ffiter;
			if( pF->uv_area() > max_area )
			{
				m_starting_face = pF;
				max_area = pF->uv_area();
			}
		}
	}

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
			double s = __orient2d( halfedgeSource( he[0] )->uv(), 
			halfedgeSource( he[1] )->uv(), 
			halfedgeSource( he[2] )->uv());

		if( a[0]/s >= 0 && a[1]/s >=0 && a[2]/s >= 0 )
		{
			 V* v[3];
			 CPoint p[3];
			 CPoint c[3];

			 for( int i = 0; i < 3; i ++ )
			 {
				 v[i] = halfedgeSource( he[i]);
				 c[i] = v[i]->rgb();
				 p[i] = v[i]->point();
			 }

			double s = __orient2d( v[0]->uv(), v[1]->uv(), v[2]->uv() );
			if( s < 1e-8 )
			{
				//q=CPoint(1,0,0);
				q = v[0]->point();
				rgb=v[0]->rgb();
			}
			else
			{
				//q = CPoint(a[0]/s, a[1]/s, a[2]/s);
			
				q = CPoint(0,0,0);
				rgb=CPoint(0,0,0);
				for( int i = 0; i < 3; i ++ )
				{
					q = q + (p[i] * a[(i+1)%3]/s);
					rgb = rgb + (c[i] * a[(i+1)%3]/s);
				}
			}

			return pF;
		}

		for( int i = 0; i < 3; i ++ )
		{
			if( a[i]/s < 0 )
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


template<typename V, typename E, typename F, typename H>
bool CQuadTreeDelaunayRemesh<V,E,F,H>::_insert_vertex(CPoint2 & pt, V * & pV )
{
    
    /* find triangle "t0" containing "v" */
	CPoint pos,rgb;

    F * pT = _locate_point( pt, pos, rgb );
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
	V * v[3];

	for( int i = 0; i < 3; i ++ )
	{
		e[i] = halfedgeEdge( pH[i] );
		v[i] = halfedgeSource( pH[i] );
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
	
	CQTSample sample;
	F * bF = m_background->_locate_point( pt, sample );

	pV->point() = sample.m_pos;
	pV->rgb() = sample.m_rgb;
	pV->lambda() = sample.m_lambda;
	pV->k() = sample.m_k;

	for( VertexFaceIterator vfiter( pV ); !vfiter.end(); vfiter ++ )
	{
		F * f = *vfiter;
		CDelaunayMesh::_statistics( f );
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
V* CQuadTreeDelaunayRemesh<V,E,F,H>::_insert_vertex_on_edge( E * e, CPoint2 & pt )
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
	
	CPoint pos, rgb;
	pV->uv() = pt;
	
	CQTSample sample;

	m_background->_locate_point( pt, sample );
	pV->point() = sample.m_pos;
	pV->rgb() = sample.m_rgb;
	pV->lambda() = sample.m_lambda;
	pV->k() = sample.m_k;
	
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
		CDelaunayMesh::_statistics( pF );
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
void CQuadTreeDelaunayRemesh<V,E,F,H>::_statistics( V* v0, V* v1, V* v2, double & minAngle, double & quality)

/* This function returns the minimum angle of a triangle (in degrees),
 * and the quality of the triangle, as required by the refinement algorithm */

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
		double lambda    = (v0->lambda() + v1->lambda() + v2->lambda())/3.0;
        quality = 1.0/lambda * amin/(perimeter*perimeter);      // return triangle quality
    
		//double k = ( fabs(v0->k() ) + fabs( v1->k()) + fabs( v2->k()))/3.0;

		//quality = quality/( 1.0 +  100*k );
	}
};

/*---------------------------------------------------------------------------*/

template<typename V, typename E, typename F, typename H>
void CQuadTreeDelaunayRemesh<V,E,F,H>::Convert( CQuadTreeDelaunayRemesh<V,E,F,H> & output_mesh )
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
	    pW->rgb()  = pV->rgb();
		pW->lambda() = pV->lambda();
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
};

/*---------------------------------------------------------------------------*/

template<typename V, typename E, typename F, typename H>
void CQuadTreeDelaunayRemesh<V,E,F,H>::_construct_tree( int n)
{

	std::vector<F*> fs;
	for( MeshFaceIterator fiter( this ); !fiter.end(); fiter ++ )
	{
		F* pF = *fiter;
		fs.push_back( pF );
	}
	m_qtree._construct( fs, n );
};


typedef CQuadTreeDelaunayRemesh<CQTDelaunayRemeshVertex, CQTDelaunayRemeshEdge, CQTDelaunayRemeshFace, CQTDelaunayRemeshHalfEdge> CQTDRMesh;


}
#endif  _DELAUNAY_REMESH_H_