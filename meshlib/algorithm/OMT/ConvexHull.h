#ifndef _POWER_DIAGRAM_H_
#define _POWER_DIAGRAM_H_

#include "Geometry/point.h"
#include <vector>
#include <iostream>
#include <queue>
#include <Eigen/Sparse>


namespace MeshLib
{
namespace OMT
{
#define SWAP(t,x,y)	{ t = x; x = y; y = t; }

#define EPSILON     1e-10

#define ONHULL   	true
#define REMOVED  	true
#define VISIBLE  	true
#define PROCESSED	true

#define X 0
#define Y 1
#define Z 2


/* Define structures for vertices, edges and faces */
typedef struct tVertexStructure tsVertex;
typedef tsVertex *tVertex;

typedef struct tEdgeStructure tsEdge;
typedef tsEdge *tEdge;

typedef struct tFaceStructure tsFace;
typedef tsFace *tFace;

typedef struct tHalfEdgeStructure tsHalfEdge;
typedef tsHalfEdge *tHalfEdge;

struct tVertexStructure {
   double   v[3];
   int	    vnum;
   tEdge    duplicate;	/* pointer to incident cone edge (or NULL) */
   bool     onhull;		/* T iff point on hull. */
   bool	    mark;		/* T iff point already processed. */
   double   target_dual_area;/* target area */
   double   dual_area;	/* current area*/
   double   weight;
   double   old_height;
   double   delta_h;
   int      m_valence;

   bool     fixed;	
   tHalfEdge halfedge;  
   tVertex  next, prev;
};

struct tEdgeStructure {
   tFace     adjface[2];
   tVertex   endpts[2];
   tHalfEdge hes[2];
   bool      next_swap;
   double    event_time;
	 
   double    a;		//the volume of the tetrahedron is 
   double    b;     //at + b

   double    weight;
   tFace    newface;            /* pointer to incident cone face. */
   bool     remove;		/* T iff edge should be delete. */
   tEdge    next, prev;
   bool     convex;
};

struct tFaceStructure {
   double    v[3];		/* dual point */
   double    n[3];      /* normal     */
   int       fnum;      /* face id    */
   tEdge     edge[3];
   tVertex   vertex[3];
   tHalfEdge halfedge;
   bool	     visible;	        /* T iff face visible from new point. */
   tFace     next, prev;
};

struct tHalfEdgeStructure{
	tFace   face;
	tVertex vert;
	tEdge   edge;
	tHalfEdge he_next, he_prev;
	tHalfEdge next, prev;
};

class CEvent
{
public:
	double m_t;
	tEdge  m_e;
	CEvent( double t, tEdge e )
	{
		m_t = t;
		m_e = e;
	};

	bool operator<( const CEvent & even ) const
	{
		return m_t > even.m_t;
	};
};



/*
 *! add p to head
 */
template<typename T>
inline void ADD(T* & head, T* p )
{
	if ( head )  { 
		p->next = head; 
		p->prev = head->prev; 
		head->prev = p; 
		p->prev->next = p; 
	} 
	else { 
		head = p; 
		head->next = head->prev = p; 
	}
};

/*
 * remove p from head to buffer
 */
template<typename T>
void DELETE(T* & head, T* & buffer, T *p )
{
	if ( head )  { 
		if ( head == head->next ) 
					head = NULL;  
		else if ( p == head ) 
			head = head->next; 
		p->next->prev = p->prev;  
		p->prev->next = p->next;  
		ADD<T>( buffer, p); 
	} 
};

/*!
 *	allocate a new object, if the buffer is not empty, get it from the buffer
 */
template<typename T>
inline T* NEW(T* & head )
{
	T * p;

	if ( head )  { 
		p = head;
		if( head->next == head )
		{
			head = NULL;
			return p;
		}
		head->next->prev = head->prev;
		head->prev->next = head->next;
		head = head->next;
		return p;
	} 
	
	p = new T;
	return p;
};


class CConvexHull
{
public:

	CConvexHull()
	{
		vertices = NULL;
		edges    = NULL;
		faces    = NULL;

		vertex_buffer = NULL;
		edge_buffer   = NULL;
		face_buffer   = NULL;

		//check = true;
		check = false;
	};
	bool    ConstructHull( void ); //if all vertices are on the hull, return false, otherwise stops and return true
	void    AddOneVertex( double x, double y, double z, double weight, bool fixed, int vid);
	void    Reset();
	void	GradientDescent();
	void    draw_mesh();
	void    draw_edges();
	void    draw_dual_mesh( bool show_next_swap_edge );
	void	draw_dual_edges();
	void    set_target_dual_area();
	void    one_step();
	bool    one_step_Newton();
	void    kinetic_one_step();

	void    hull(); //reset and construuct a convex hull of current vertices

	void	print( const char * name );
	void    print_dual( const char * name );

	void	VertexGradient();
	void	VertexNewton();
	double  VertexDualAreaError();
	void    EdgeEventTime();

public:
	/* Function declarations */
	tVertex MakeNullVertex();
	tEdge   MakeNullEdge( void );
	tFace   MakeNullFace( void );
	tFace   MakeFace( tVertex v0, tVertex v1, tVertex v2, tFace f );
	tFace	MakeConeFace( tEdge e, tVertex p );
	void    MakeCcw( tFace f, tEdge e, tVertex p );

	void    ReadVertices( void );
	void    SubVec( int a[3], int b[3], int c[3]);
	void    DoubleTriangle( void );
	bool	AddOne( tVertex p );
	int     VolumeSign(tFace f, tVertex p);
	
	bool    CleanUp( tVertex *pvnext );
	void    CleanEdges( void );
	void    CleanFaces( void );
	bool    CleanVertices( tVertex *pvnext );
	bool	Collinear( tVertex a, tVertex b, tVertex c );
	
	void	EdgeOrderOnFaces ( void );
	void    AttachHalfEdges( void );
	void    FaceDualPoint( void );
	void    VertexDualArea( void );
	void	EdgeWeight( void ); 
	void    VertexWeight( void );
	void    VertexValence( void );
	
	double  EventTime( tEdge e );
	bool	EdgeConvexity();
	void    FixEdge( double time );

private:
	void   _face_dual_point(  tFace  f );
	void   _vertex_dual_area( tVertex v );
protected:

	tVertex vertices;
	tEdge   edges;
	tFace   faces;
	tHalfEdge halfedges;

	tVertex vertex_buffer;
	tEdge   edge_buffer;
	tFace   face_buffer;
	tHalfEdge halfedge_buffer;

protected:

	void    CheckEuler(int V, int E, int F );
	void    Checks( void );
	void	Consistency( void );
	void	CheckEndpts( void );
	bool    check;

public:

	void ___Compute_Event_Time();
	void ___Update( double delta_t );
	tEdge ___next_swap_edge( double current_time );
	//havn't tested
	void    EdgeSwap( tEdge e, double current_time );
	bool	Convexity( void );

	tEdge   swap_edge;
	double  current_time;
};


/*---------------------------------------------------------------------
MakeNullVertex: Makes a vertex, nulls out fields.
---------------------------------------------------------------------*/
tVertex	CConvexHull::MakeNullVertex()
{
   tVertex  v;

   v = NEW<tsVertex>( vertex_buffer );

   assert( v != NULL );
   v->duplicate = NULL;
   v->onhull = !ONHULL;
   v->mark = !PROCESSED;
   ADD<tsVertex>( vertices, v );
	
   return v;
};

/*---------------------------------------------------------------------
MakeNullEdge creates a new cell and initializes all pointers to NULL
and sets all flags to off.  It returns a pointer to the empty cell.
---------------------------------------------------------------------*/
tEdge CConvexHull::MakeNullEdge()
{
   tEdge  e;

   e = NEW<tsEdge>( edge_buffer );

   e->adjface[0] = e->adjface[1] = e->newface = NULL;
   e->endpts[0] = e->endpts[1] = NULL;
   e->remove = !REMOVED;
   ADD<tsEdge>( edges, e );
   return e;
};

/*--------------------------------------------------------------------
MakeNullFace creates a new face structure and initializes all of its
flags to NULL and sets all the flags to off.  It returns a pointer
to the empty cell.
---------------------------------------------------------------------*/
tFace 	CConvexHull::MakeNullFace( void )
{
   tFace  f;

   f = NEW<tsFace>( face_buffer );

   for (int i=0; i < 3; ++i ) {
      f->edge[i] = NULL;
      f->vertex[i] = NULL;
   }
   f->visible = !VISIBLE;
   ADD<tsFace>( faces, f );
   return f;
};

/*---------------------------------------------------------------------
MakeFace creates a new face structure from three vertices (in ccw
order).  It returns a pointer to the face.
---------------------------------------------------------------------*/
tFace CConvexHull::MakeFace( tVertex v0, tVertex v1, tVertex v2, tFace fold )
{
   tFace  f;
   tEdge  e0, e1, e2;

   /* Create edges of the initial triangle. */
   if( !fold ) {
     e0 = MakeNullEdge();
     e1 = MakeNullEdge();
     e2 = MakeNullEdge();
   }
   else { /* Copy from fold, in reverse order. */
     e0 = fold->edge[2];
     e1 = fold->edge[1];
     e2 = fold->edge[0];
   }
   e0->endpts[0] = v0;              e0->endpts[1] = v1;
   e1->endpts[0] = v1;              e1->endpts[1] = v2;
   e2->endpts[0] = v2;              e2->endpts[1] = v0;
	
   /* Create face for triangle. */
   f = MakeNullFace();
   f->edge[0]   = e0;  f->edge[1]   = e1; f->edge[2]   = e2;
   f->vertex[0] = v0;  f->vertex[1] = v1; f->vertex[2] = v2;
	
   /* Link edges to face. */
   e0->adjface[0] = e1->adjface[0] = e2->adjface[0] = f;
	
   return f;
};

/*---------------------------------------------------------------------
CleanUp goes through each data structure list and clears all
flags and NULLs out some pointers.  The order of processing
(edges, faces, vertices) is important.
---------------------------------------------------------------------*/
bool	CConvexHull::CleanUp( tVertex *pvnext )
{
   CleanEdges();
   CleanFaces();
   return CleanVertices( pvnext );
}

/*---------------------------------------------------------------------
CleanEdges runs through the edge list and cleans up the structure.
If there is a newface then it will put that face in place of the 
visible face and NULL out newface. It also deletes so marked edges.
---------------------------------------------------------------------*/
void	CConvexHull::CleanEdges( void )
{
   tEdge  e;	/* Primary index into edge list. */
   tEdge  t;	/* Temporary edge pointer. */
		
   /* Integrate the newface's into the data structure. */
   /* Check every edge. */
   e = edges;
   do {
      if ( e->newface ) { 
	 if ( e->adjface[0]->visible )
	    e->adjface[0] = e->newface; 
	 else	e->adjface[1] = e->newface;
	 e->newface = NULL;
      }
      e = e->next;
   } while ( e != edges );

   /* Delete any edges marked for deletion. */
   while ( edges && edges->remove ) { 
      e = edges;
      DELETE<tsEdge>( edges, edge_buffer, e );
   }
   e = edges->next;
   do {
      if ( e->remove ) {
	 t = e;
	 e = e->next;
	 DELETE<tsEdge>( edges, edge_buffer, t );
      }
      else e = e->next;
   } while ( e != edges );
}

/*---------------------------------------------------------------------
CleanFaces runs through the face list and deletes any face marked visible.
---------------------------------------------------------------------*/
void	CConvexHull::CleanFaces()
{
   tFace  f;	/* Primary pointer into face list. */
   tFace  t;	/* Temporary pointer, for deleting. */
	

   while ( faces && faces->visible ) { 
      f = faces;
      DELETE<tsFace>( faces, face_buffer, f );
   }
   f = faces->next;
   do {
      if ( f->visible ) {
	 t = f;
	 f = f->next;
	 DELETE<tsFace>( faces, face_buffer, t );
      }
      else f = f->next;
   } while ( f != faces );
};

/*---------------------------------------------------------------------
CleanVertices runs through the vertex list and deletes the 
vertices that are marked as processed but are not incident to any 
undeleted edges. 
The pointer to vnext, pvnext, is used to alter vnext in
ConstructHull() if we are about to delete vnext.
---------------------------------------------------------------------*/
bool CConvexHull::CleanVertices( tVertex *pvnext )
{
   bool hidden_vertex = false;
   tEdge    e;
   tVertex  v, t;
	
   /* Mark all vertices incident to some undeleted edge as on the hull. */
   e = edges;
   do {
      e->endpts[0]->onhull = e->endpts[1]->onhull = ONHULL;
      e = e->next;
   } while (e != edges);
	
   /* Delete all vertices that have been processed but
      are not on the hull. */
   while ( vertices && vertices->mark && !vertices->onhull ) { 
      /* If about to delete vnext, advance it first. */
      v = vertices;
      if( v == *pvnext )
         *pvnext = v->next;
      DELETE<tsVertex>( vertices, vertex_buffer, v );
	  hidden_vertex = true;
   }
   v = vertices->next;
   do {
      if ( v->mark && !v->onhull ) {    
	 t = v; 
	 v = v->next;
	 if( t == *pvnext )
         *pvnext = t->next;
	 DELETE<tsVertex>( vertices, vertex_buffer, t );
	 hidden_vertex = true;
      }
      else v = v->next;
   } while ( v != vertices );
	
   /* Reset flags. */
   v = vertices;
   do {
      v->duplicate = NULL; 
      v->onhull = !ONHULL; 
      v = v->next;
   } while ( v != vertices );

   return hidden_vertex;
};

/*---------------------------------------------------------------------
Collinear checks to see if the three points given are collinear,
by checking to see if each element of the cross product is zero.
---------------------------------------------------------------------*/
bool CConvexHull::Collinear( tVertex a, tVertex b, tVertex c )
{
   return 
         ( c->v[Z] - a->v[Z] ) * ( b->v[Y] - a->v[Y] ) -
         ( b->v[Z] - a->v[Z] ) * ( c->v[Y] - a->v[Y] ) == 0
      && ( b->v[Z] - a->v[Z] ) * ( c->v[X] - a->v[X] ) -
         ( b->v[X] - a->v[X] ) * ( c->v[Z] - a->v[Z] ) == 0
      && ( b->v[X] - a->v[X] ) * ( c->v[Y] - a->v[Y] ) -
         ( b->v[Y] - a->v[Y] ) * ( c->v[X] - a->v[X] ) == 0  ;
};

/*---------------------------------------------------------------------
AddOne is passed a vertex.  It first determines all faces visible from 
that point.  If none are visible then the point is marked as not 
onhull.  Next is a loop over edges.  If both faces adjacent to an edge
are visible, then the edge is marked for deletion.  If just one of the
adjacent faces is visible then a new face is constructed.
---------------------------------------------------------------------*/
bool CConvexHull::AddOne( tVertex p )
{
  
   tFace  f; 
   tEdge  e, temp;
   int 	  vol;
   bool	  vis = false;


   /* Mark faces visible from p. */
   f = faces;
   do {
      vol = VolumeSign( f, p );
      if ( vol < 0 ) {
		f->visible = VISIBLE;  
		vis = true;                      
      }
      f = f->next;
   } while ( f != faces );

   /* If no faces are visible from p, then p is inside the hull. */
   if ( !vis ) {
      p->onhull = !ONHULL;  
      return false; 
   }

   /* Mark edges in interior of visible region for deletion.
      Erect a newface based on each border edge. */
   e = edges;
   do {
      temp = e->next;
      if ( e->adjface[0]->visible && e->adjface[1]->visible )
	 /* e interior: mark for deletion. */
	 e->remove = REMOVED;
      else if ( e->adjface[0]->visible || e->adjface[1]->visible ) 
	 /* e border: make a new face. */
	 e->newface = MakeConeFace( e, p );
      e = temp;
   } while ( e != edges );
   return true;
};

/*---------------------------------------------------------------------
VolumeSign returns the sign of the volume of the tetrahedron determined by f
and p.  VolumeSign is +1 iff p is on the negative side of f,
where the positive side is determined by the rh-rule.  So the volume 
is positive if the ccw normal to f points outside the tetrahedron.
The final fewer-multiplications form is due to Bob Williamson.
---------------------------------------------------------------------*/
int  CConvexHull::VolumeSign( tFace f, tVertex p )
{
   double  vol;
   double  ax, ay, az, bx, by, bz, cx, cy, cz;

   ax = f->vertex[0]->v[X] - p->v[X];
   ay = f->vertex[0]->v[Y] - p->v[Y];
   az = f->vertex[0]->v[Z] - p->v[Z];
   bx = f->vertex[1]->v[X] - p->v[X];
   by = f->vertex[1]->v[Y] - p->v[Y];
   bz = f->vertex[1]->v[Z] - p->v[Z];
   cx = f->vertex[2]->v[X] - p->v[X];
   cy = f->vertex[2]->v[Y] - p->v[Y];
   cz = f->vertex[2]->v[Z] - p->v[Z];

   vol =   ax * (by*cz - bz*cy)
         + ay * (bz*cx - bx*cz)
         + az * (bx*cy - by*cx);


   /* The volume should be an integer. */
/*
   if      ( vol >  0.0 )  return  1;
   else if ( vol <  0.0 )  return -1;
   else                     return  0;
*/
   /* The volume should be an integer. */
   if      ( vol >  EPSILON )  return  1;
   else if ( vol < -EPSILON )  return -1;
   else                     return  0;
};

/*---------------------------------------------------------------------
MakeConeFace makes a new face and two new edges between the 
edge and the point that are passed to it. It returns a pointer to
the new face.
---------------------------------------------------------------------*/
tFace	CConvexHull::MakeConeFace( tEdge e, tVertex p )
{
   tEdge  new_edge[2];
   tFace  new_face;
   int 	  i, j;

   /* Make two new edges (if don't already exist). */
   for ( i=0; i < 2; ++i ) 
      /* If the edge exists, copy it into new_edge. */
      if ( !( new_edge[i] = e->endpts[i]->duplicate) ) {
	 /* Otherwise (duplicate is NULL), MakeNullEdge. */
	 new_edge[i] = MakeNullEdge();
	 new_edge[i]->endpts[0] = e->endpts[i];
	 new_edge[i]->endpts[1] = p;
	 e->endpts[i]->duplicate = new_edge[i];
      }

   /* Make the new face. */
   new_face = MakeNullFace();   
   new_face->edge[0] = e;
   new_face->edge[1] = new_edge[0];
   new_face->edge[2] = new_edge[1];
   MakeCcw( new_face, e, p ); 
        
   /* Set the adjacent face pointers. */
   for ( i=0; i < 2; ++i )
      for ( j=0; j < 2; ++j )  
	 /* Only one NULL link should be set to new_face. */
	 if ( !new_edge[i]->adjface[j] ) {
	    new_edge[i]->adjface[j] = new_face;
	    break;
	 }
        
   return new_face;
};

/*---------------------------------------------------------------------
MakeCcw puts the vertices in the face structure in counterclock wise 
order.  We want to store the vertices in the same 
order as in the visible face.  The third vertex is always p.

Although no specific ordering of the edges of a face are used
by the code, the following condition is maintained for each face f:
one of the two endpoints of f->edge[i] matches f->vertex[i]. 
But note that this does not imply that f->edge[i] is between
f->vertex[i] and f->vertex[(i+1)%3].  (Thanks to Bob Williamson.)
---------------------------------------------------------------------*/
void CConvexHull::MakeCcw( tFace f, tEdge e, tVertex p )
{
   tFace  fv;   /* The visible face adjacent to e */
   int    i;    /* Index of e->endpoint[0] in fv. */
   tEdge  s;	/* Temporary, for swapping */
      
   if  ( e->adjface[0]->visible )      
        fv = e->adjface[0];
   else fv = e->adjface[1];
       
   /* Set vertex[0] & [1] of f to have the same orientation
      as do the corresponding vertices of fv. */ 
   for ( i=0; fv->vertex[i] != e->endpts[0]; ++i )
      ;
   /* Orient f the same as fv. */
   if ( fv->vertex[ (i+1) % 3 ] != e->endpts[1] ) {
      f->vertex[0] = e->endpts[1];  
      f->vertex[1] = e->endpts[0];    
   }
   else {                               
      f->vertex[0] = e->endpts[0];   
      f->vertex[1] = e->endpts[1];      
      SWAP( s, f->edge[1], f->edge[2] );
   }
   /* This swap is tricky. e is edge[0]. edge[1] is based on endpt[0],
      edge[2] on endpt[1].  So if e is oriented "forwards," we
      need to move edge[1] to follow [0], because it precedes. */
   
   f->vertex[2] = p;
};

/*---------------------------------------------------------------------
ConstructHull adds the vertices to the hull one at a time.  The hull
vertices are those in the list marked as onhull.
---------------------------------------------------------------------*/
bool CConvexHull::ConstructHull()
{
   DoubleTriangle();

   tVertex  v, vnext;
   bool	    changed;	/* T if addition changes hull; not used. */

   v = vertices;
   do {
      vnext = v->next;
      if ( !v->mark ) {
        v->mark = PROCESSED;
		changed = AddOne( v );
		if( CleanUp( &vnext ) )
		{
			return true;
		}
		/* Pass down vnext in case it gets deleted. */
		//fprintf(stderr,"%d\n", v->vnum );
	  }
/*
	  if( check )
	  {
		  fprintf( stderr, "\n\nVertex %d\n", v->vnum );
		  Checks();
		  
	  }
	if( v->next == vertices )
	  {
		  char name[256];
		  sprintf(name,"debug_%d.m", v->vnum );
		  print(name);
	  }
*/
		v = vnext;
   } while ( v != vertices );

   EdgeOrderOnFaces();
   AttachHalfEdges();
   FaceDualPoint();
   VertexDualArea();

   return false;
};


/*---------------------------------------------------------------------
 DoubleTriangle builds the initial double triangle.  It first finds 3 
 noncollinear points and makes two faces out of them, in opposite order.
 It then finds a fourth point that is not coplanar with that face.  The  
 vertices are stored in the face structure in counterclockwise order so 
 that the volume between the face and the point is negative. Lastly, the
 3 newfaces to the fourth point are constructed and the data structures
 are cleaned up. 
---------------------------------------------------------------------*/
void  CConvexHull::DoubleTriangle()
{
   tVertex  v0, v1, v2, v3;
   tFace    f0, f1 = NULL;
   int      vol;
	
   /* Find 3 noncollinear points. */
   v0 = vertices;
   while ( Collinear( v0, v0->next, v0->next->next ) )
      if ( ( v0 = v0->next ) == vertices )
         printf("DoubleTriangle:  All points are Collinear!\n"), exit(0);
   v1 = v0->next;
   v2 = v1->next;
	
   /* Mark the vertices as processed. */
   v0->mark = PROCESSED;
   v1->mark = PROCESSED;
   v2->mark = PROCESSED;
   
   /* Create the two "twin" faces. */
   f0 = MakeFace( v0, v1, v2, f1 );
   f1 = MakeFace( v2, v1, v0, f0 );

   /* Link adjacent face fields. */
   f0->edge[0]->adjface[1] = f1;
   f0->edge[1]->adjface[1] = f1;
   f0->edge[2]->adjface[1] = f1;
   f1->edge[0]->adjface[1] = f0;
   f1->edge[1]->adjface[1] = f0;
   f1->edge[2]->adjface[1] = f0;
	
   /* Find a fourth, noncoplanar point to form tetrahedron. */
   v3 = v2->next;
   vol = VolumeSign( f0, v3 );
   while ( !vol )   {
      if ( ( v3 = v3->next ) == v0 ) 
         printf("DoubleTriangle:  All points are coplanar!\n"), exit(0);
      vol = VolumeSign( f0, v3 );
   }
	
   /* Insure that v3 will be the first added. */
   vertices = v3;
};

/*------------------------------------------------------------------
  EdgeOrderOnFaces: puts e0 between v0 and v1, e1 between v1 and v2,
  e2 between v2 and v0 on each face.  This should be unnecessary, alas.
  Not used in code, but useful for other purposes.
------------------------------------------------------------------*/
void CConvexHull::EdgeOrderOnFaces () 
{
  tFace f = faces;
  tEdge new_edge;
  int i,j;

  do {
    for (i = 0; i < 3; i++) {
      if (!(((f->edge[i]->endpts[0] == f->vertex[i]) &&
             (f->edge[i]->endpts[1] == f->vertex[(i+1)%3])) ||
            ((f->edge[i]->endpts[1] == f->vertex[i]) &&
             (f->edge[i]->endpts[0] == f->vertex[(i+1)%3])))) {
        /* Change the order of the edges on the face: */
        for (j = 0; j < 3; j ++) {
          /* find the edge that should be there */
          if (((f->edge[j]->endpts[0] == f->vertex[i]) &&
               (f->edge[j]->endpts[1] == f->vertex[(i+1)%3])) ||
              ((f->edge[j]->endpts[1] == f->vertex[i]) &&
               (f->edge[j]->endpts[0] == f->vertex[(i+1)%3]))) {
            /* Swap it with the one erroneously put into its place: */
            new_edge = f->edge[i];
            f->edge[i] = f->edge[j];
            f->edge[j] = new_edge;
          }
        }
      }
    }
    f = f->next;
  } while (f != faces);

};

/*---------------------------------------------------------------------
ReadVertices: Reads in the vertices, and links them into a circular
list with MakeNullVertex.  There is no need for the # of vertices to be
the first line: the function looks for EOF instead.  Sets the global
variable vertices via the ADD macro.
---------------------------------------------------------------------*/
void CConvexHull::AddOneVertex( double x, double y, double z, double weight, bool fixed, int vid )
{
   tVertex  v;

	v = MakeNullVertex();
	v->v[X] = x;
	v->v[Y] = y;
	v->v[Z] = z;
	v->fixed = fixed;
	v->target_dual_area = weight;
	v->vnum = vid;
};




/*---------------------------------------------------------------------
Consistency runs through the edge list and checks that all
adjacent faces have their endpoints in opposite order.  This verifies
that the vertices are in counterclockwise order.
---------------------------------------------------------------------*/
void	CConvexHull::Consistency( void )
{
   register tEdge  e;
   register int    i, j;

   e = edges;

   do {
      /* find index of endpoint[0] in adjacent face[0] */
      for ( i = 0; e->adjface[0]->vertex[i] != e->endpts[0]; ++i )
	 ;
   
      /* find index of endpoint[0] in adjacent face[1] */
      for ( j = 0; e->adjface[1]->vertex[j] != e->endpts[0]; ++j )
	 ;

      /* check if the endpoints occur in opposite order */
      if ( !( e->adjface[0]->vertex[ (i+1) % 3 ] ==
	      e->adjface[1]->vertex[ (j+2) % 3 ] ||
	      e->adjface[0]->vertex[ (i+2) % 3 ] ==
	      e->adjface[1]->vertex[ (j+1) % 3 ] )  )
	 break;
      e = e->next;

   } while ( e != edges );

   if ( e != edges )
      fprintf( stderr, "Checks: edges are NOT consistent.\n");
   else
      fprintf( stderr, "Checks: edges consistent.\n");

}

/*---------------------------------------------------------------------
Convexity checks that the volume between every face and every
point is negative.  This shows that each point is inside every face
and therefore the hull is convex.
---------------------------------------------------------------------*/
bool	CConvexHull::Convexity( void )
{
   register tFace    f;
   register tVertex  v;
   int               vol;

   f = faces;
   
   do {
      v = vertices;
      do {
	 if ( v->mark ) 
	 {
	    vol = VolumeSign( f, v );
	    if ( vol < 0 )
		{
		   printf("Degenerated v %d\n", v->vnum );
		   tHalfEdge h = f->halfedge;
		   for(int i = 0; i < 3; i ++ )
		   {
			   printf("Degenerated v %d\n", h->vert->vnum );
			   h = h->he_next;
		   }
	       return false;
		}
	 }
	 v = v->next;
      } while ( v != vertices );

      f = f->next;

   } while ( f != faces );

   if ( f != faces )
      fprintf( stderr, "Checks: NOT convex.\n");
   else if ( check ) 
      fprintf( stderr, "Checks: convex.\n");

   return true;
}

/*---------------------------------------------------------------------
CheckEuler checks Euler's relation, as well as its implications when
all faces are known to be triangles.  Only prints positive information
when debug is true, but always prints negative information.
---------------------------------------------------------------------*/
void	CConvexHull::CheckEuler( int V, int E, int F )
{
   if ( check )
      fprintf( stderr, "Checks: V, E, F = %d %d %d:\t", V, E, F);

   if ( (V - E + F) != 2 )
      fprintf( stderr, "Checks: V-E+F != 2\n");
   else if ( check )
      fprintf( stderr, "V-E+F = 2\t");


   if ( F != (2 * V - 4) )
      fprintf( stderr, "Checks: F=%d != 2V-4=%d; V=%d\n",
	      F, 2*V-4, V);
   else if ( check ) 
      fprintf( stderr, "F = 2V-4\t");
   
   if ( (2 * E) != (3 * F) )
      fprintf( stderr, "Checks: 2E=%d != 3F=%d; E=%d, F=%d\n",
	      2*E, 3*F, E, F );
   else if ( check ) 
      fprintf( stderr, "2E = 3F\n");
}

/*-------------------------------------------------------------------*/
void	CConvexHull::Checks( void )
{
   tVertex  v;
   tEdge    e;
   tFace    f;
   int 	   V = 0, E = 0 , F = 0;

   Consistency();
   Convexity();
   if ( v = vertices )
      do {
         if (v->mark) V++;
	 v = v->next;
      } while ( v != vertices );
   if ( e = edges )
      do {
         E++;
	 e = e->next;
      } while ( e != edges );
   if ( f = faces )
      do {
         F++;
	 f  = f ->next;
      } while ( f  != faces );
   CheckEuler( V, E, F );
   CheckEndpts();
}

/*-------------------------------------------------------------------
Checks that, for each face, for each i={0,1,2}, the [i]th vertex of
that face is either the [0]th or [1]st endpoint of the [ith] edge of
the face.
-------------------------------------------------------------------*/
void CConvexHull::CheckEndpts ( void )
{
   int 	   i;
   tFace   fstart;
   tEdge   e;
   tVertex v;
   bool error = false;

   fstart = faces;
   if (faces) do {
      for( i=0; i<3; ++i ) {
         v = faces->vertex[i];
         e = faces->edge[i];
         if ( v != e->endpts[0] && v != e->endpts[1] ) {
            error = true;
            fprintf(stderr,"CheckEndpts: Error!\n");
            fprintf(stderr,"  addr: %8x;", faces );
            fprintf(stderr,"  edges:");
            fprintf(stderr,"(%3d,%3d)", 
               e->endpts[0]->vnum,
               e->endpts[1]->vnum);
            fprintf(stderr,"\n");
         }
      }
      faces= faces->next;
   } while ( faces != fstart );

   if ( error )
     fprintf(stderr,"Checks: ERROR found and reported above.\n");
   else
     fprintf(stderr,"Checks: All endpts of all edges of all faces check.\n");

}


void CConvexHull::print( const char * name )
{
	FILE * fp = fopen( name, "w" );

	tVertex v = vertices;
	do{
		if( v->mark != PROCESSED )
		{
			v = v->next;
			continue;
		}
		fprintf(fp, "Vertex %d %f %f %f\n", 
			v->vnum,
		v->v[0], v->v[1],v->v[2]);
		v = v->next;
	}while( v!= vertices );

	tFace  f = faces;
	int fid = 1;
	do{
		fprintf(fp, "Face %d %d %d %d\n",
			fid++,
			f->vertex[0]->vnum,
			f->vertex[1]->vnum,
			f->vertex[2]->vnum);
		f = f->next;
	}while( f!= faces );

	fclose(fp);
};

void CConvexHull::AttachHalfEdges()
{
	tFace f = faces;
	do{
		tHalfEdge he[3];
		for( int i = 0; i < 3; i ++ )
		{
			he[i] = new tsHalfEdge;
			ADD<tsHalfEdge>( halfedges, he[i] );
		}
		for( int i = 0; i < 3; i ++ )
		{
			he[i]->he_next = he[(i+1)%3];
			he[i]->he_prev = he[(i+2)%3];
		}
		for( int i = 0; i < 3; i ++ )
		{
			he[i]->face = f;
		}
		f->halfedge = he[0];
		for( int i = 0; i < 3; i ++ )
		{
			he[i]->edge = f->edge[i];
		}
		for( int i = 0; i < 3; i ++ )
		{
			he[i]->vert = f->vertex[(i+1)%3];
			he[i]->vert->halfedge = he[i];
		}
		
		f = f->next;
	}while( f!= faces );

	tEdge e = edges;
	do{
		tFace f0 = e->adjface[0];
		tHalfEdge h = f0->halfedge;
		for( int i = 0; i < 3; i ++ )
		{
			if( h->edge == e )
			{
				e->hes[0] = h;
				break;
			}
			h = h->next;
		}
		tFace f1 = e->adjface[1];
		h = f1->halfedge;
		for( int i = 0; i < 3; i ++ )
		{
			if( h->edge == e )
			{
				e->hes[1] = h;
				break;
			}
			h = h->next;
		}
		e = e->next;
	}while( e!= edges );
};

void CConvexHull::_face_dual_point( tFace f )
{
		CPoint p[3];			
		for( int i = 0; i < 3; i ++ )
		{
			tVertex v = f->vertex[i];
			p[i] = CPoint( v->v[0], v->v[1], v->v[2] );
		}

		CPoint n = (p[1]-p[0])^(p[2]-p[0]);

		double A =  -n[X]/n[Z];
		double B =  -n[Y]/n[Z];
		double C =  p[0]*n/n[Z];

		f->v[0] = A;
		f->v[1] = B;
		f->v[2] = -C;

		n /= n.norm();

		f->n[X] = n[0];
		f->n[Y] = n[1];
		f->n[Z] = n[2];

};

void CConvexHull::FaceDualPoint()
{
	int fnum = 1;
	tFace f = faces;
	do{
		f->fnum = fnum ++;
		_face_dual_point( f );
		f = f->next;
	}while( f!= faces );

};


void CConvexHull::print_dual( const char * name )
{
	FILE * fp = fopen( name, "w" );

	tFace   f = faces;
	do{
		fprintf(fp, "Vertex %d %f %f %f\n", 
			f->fnum,
		f->v[0], f->v[1],f->v[2]);
		f = f->next;
	}while( f!= faces );

	tVertex v = vertices;
	do{

		tHalfEdge h = v->halfedge;
		std::vector<int> ids;
		do{
			tFace f = h->face;
			ids.push_back( f->fnum );
	
			tEdge e = h->edge;
			h = (e->hes[0]!= h)? e->hes[0]:e->hes[1];
			h = h->he_prev;

		}while( h!= v->halfedge );

		


		fprintf(fp, "Face %d", v->vnum);
		for( size_t i = 0; i < ids.size(); i ++)
		{
			fprintf(fp, " %d", ids[i]);
		}
		fprintf(fp,"\n");
		
		v = v->next;
	}while( v!= vertices );

	fclose(fp);
};

void CConvexHull::Reset()
{
	while( faces )
	{
		tFace f = faces;
		DELETE( faces, face_buffer, f );
	}

	while( edges )
	{
		tEdge e = edges;
		DELETE( edges, edge_buffer, e );
	}

	while( halfedges )
	{
		tHalfEdge h = halfedges;
		DELETE( halfedges, halfedge_buffer, h );
	}

	tVertex v = vertices;
	do{
		v->duplicate = NULL;
		v->onhull = !ONHULL;
		v->mark = !PROCESSED;
		v = v->next;
	}while( v!= vertices );
	if( vertex_buffer )
	{
		while(  vertex_buffer )
		{
			v = vertex_buffer;
			DELETE( vertex_buffer, vertices, v);
			v->duplicate = NULL;
			v->onhull = !ONHULL;
			v->mark = !PROCESSED;
		}
	}
};

void CConvexHull::_vertex_dual_area( tVertex v )
{
		tHalfEdge h = v->halfedge;
		std::vector<CPoint2> pts;
		do{
			tFace f = h->face;
			CPoint2 p(f->v[0], f->v[1]);
			pts.push_back( p );
	
			tEdge e = h->edge;
			h = (e->hes[0]!= h)? e->hes[0]:e->hes[1];
			h = h->he_prev;

		}while( h!= v->halfedge );

		double s = 0;
		for( size_t i = 1; i < pts.size() -1 ; i ++ )
		{
			s += (pts[i]-pts[0])^(pts[i+1]-pts[0]);
		}	
		v->dual_area = s/2.0;
};

void CConvexHull::VertexDualArea()
{
	tVertex v = vertices;
	do{
		_vertex_dual_area( v );
		v = v->next;
	}while( v!= vertices );
};

void CConvexHull::GradientDescent()
{
	int i = 0;

	tVertex v = vertices;
	do{
		v->target_dual_area = 0.1;
		v = v->next;
	}while( v!= vertices );

	int k = 0;
	do{

	ConstructHull();
	char name[256];
	sprintf(name, "debug_%d.m", k ++ );
	print_dual( name );

	v = vertices;
	double error = 0;
	do{
		double d = fabs(v->dual_area - v->target_dual_area);
		error  = (error > d )?error:d;

		v = v->next;
	}while( v!= vertices );
	fprintf(stderr, "Error is %f\n", error );
	if( error < 1e-3 ) break;

	v = vertices;
	do{
		if( v->fixed ) 
		{
			v = v->next;
			continue;
		}
		double d = v->target_dual_area - v->dual_area;
		v->v[2] += d * 0.2;
		v = v->next;
	}while( v!= vertices );
	Reset();



	}while( true );
};


/*! draw mesh */
void CConvexHull::draw_mesh()
{
	glColor3d(  229.0/255.0, 162.0/255.0, 141.0/255.0 );
	 glEnable( GL_POLYGON_OFFSET_FILL );      
     glPolygonOffset( 1.0, 1.0 );

	  glBegin(GL_TRIANGLES);

	  tFace f = faces;

	  do{
		  glNormal3d( f->n[X], f->n[Y], f->n[Z] );
		  for( int i = 0; i < 3; i ++ )
		  {
			  tVertex v = f->vertex[i];
			  glVertex3f( v->v[X], v->v[Y], v->v[Z] );
		  }
		  for( int i = 0; i < 3; i ++ )
		  {
			  tVertex v = f->vertex[i];
			  glVertex3f( v->v[X], v->v[Y], -2 );
		  }

		  f = f->next;
	  }while( f != faces );

	glEnd();
	glDisable( GL_POLYGON_OFFSET_FILL );      

}

void CConvexHull::draw_edges()
{
	glDisable(GL_LIGHTING );
	glColor3f(0,0,0);
	glBegin(GL_LINES);

	tEdge e = edges;
	do{
		if( e->next_swap )
		{
			glColor3f(1,0,0);
		}
		else
		{
			glColor3f(0,0,0);
		}
		if( !e->convex )
		{
			glColor3f(0,0,1);
		}

		tVertex v0 = e->endpts[0];
        glVertex3d ( v0->v[X], v0->v[Y], v0->v[Z]);
		tVertex v1 = e->endpts[1];
        glVertex3d ( v1->v[X], v1->v[Y], v1->v[Z]);
		e = e->next;
	}while( e != edges );

	 glEnd();
	 glEnable(GL_LIGHTING );
	 
};


/*! draw mesh */
void CConvexHull::draw_dual_mesh( bool show_next_swap_edge )
{
	glColor3d(  229.0/255.0, 162.0/255.0, 141.0/255.0 );
	 glEnable( GL_POLYGON_OFFSET_FILL );      
     glPolygonOffset( 1.0, 1.0 );

	  glBegin(GL_TRIANGLES);


	 tVertex v = vertices;
	 do{
		 if( v->fixed ) 
		 {
			 v = v->next;
			 continue;
		 }
		tHalfEdge h = v->halfedge;
		std::vector<CPoint> pts;
		do{
			tFace f = h->face;
			CPoint q(f->v[X],f->v[Y],f->v[Z] );
			pts.push_back( q );
	
			tEdge e = h->edge;
			h = (e->hes[0]!= h)? e->hes[0]:e->hes[1];
			h = h->he_prev;

		}while( h!= v->halfedge );

		CPoint n = (pts[1]-pts[0])^(pts[2]-pts[0]);

		for( size_t i = 2; i < pts.size(); i ++ )
		{
			CPoint d = (pts[1]-pts[0])^(pts[i]-pts[0]);
			if( d.norm() > n.norm() ) n = d;
		}

		n /= n.norm();
		glNormal3d( n[X], n[Y], n[Z] );

		  for( size_t i = 1; i < pts.size()-1; i ++ )
		  {
			  glVertex3f( pts[0][X], pts[0][Y], pts[0][Z] );
			  glVertex3f( pts[i][X], pts[i][Y], pts[i][Z] );
			  glVertex3f( pts[i+1][X], pts[i+1][Y], pts[i+1][Z] );
		  }
		  v = v->next;
	  }while( v != vertices );

	glEnd();

	if( show_next_swap_edge && swap_edge != NULL )
	{
		tFace f0 = swap_edge->hes[0]->face;

		glPushMatrix();
		glTranslatef( f0->v[0], f0->v[1], f0->v[2] );
		glColor3f(0,1,0);
		glutSolidSphere( 0.1, 8, 8 );
		glPopMatrix();
	}	

	//draw back faces
	{
	glFrontFace(GL_CW); 
		//glCullFace(GL_FRONT);
	glColor3f(0,1,0);
	 glEnable( GL_POLYGON_OFFSET_FILL );      
     glPolygonOffset( 1.0, 1.0 );

	  glBegin(GL_TRIANGLES);


	 tVertex v = vertices;
	 do{
		 if( v->fixed ) 
		 {
			 v = v->next;
			 continue;
		 }
		tHalfEdge h = v->halfedge;
		std::vector<CPoint> pts;
		do{
			tFace f = h->face;
			CPoint q(f->v[X],f->v[Y],f->v[Z] );
			pts.push_back( q );
	
			tEdge e = h->edge;
			h = (e->hes[0]!= h)? e->hes[0]:e->hes[1];
			h = h->he_prev;

		}while( h!= v->halfedge );

		CPoint n = (pts[1]-pts[0])^(pts[2]-pts[0]);

		for( size_t i = 2; i < pts.size(); i ++ )
		{
			CPoint d = (pts[1]-pts[0])^(pts[i]-pts[0]);
			if( d.norm() > n.norm() ) n = d;
		}

		n /= n.norm();
		glNormal3d( -n[X], -n[Y], -n[Z] );

		  for( size_t i = 1; i < pts.size()-1; i ++ )
		  {
			  glVertex3f( pts[0][X], pts[0][Y], pts[0][Z] );
			  glVertex3f( pts[i][X], pts[i][Y], pts[i][Z] );
			  glVertex3f( pts[i+1][X], pts[i+1][Y], pts[i+1][Z] );
		  }
		  v = v->next;
	  }while( v != vertices );

	glEnd();
	glFrontFace(GL_CCW); 
	}	

	//draw planar projection
	{
	glFrontFace(GL_CW); 
		//glCullFace(GL_FRONT);
	glColor3d(  229.0/255.0, 162.0/255.0, 141.0/255.0 );
	 glEnable( GL_POLYGON_OFFSET_FILL );      
     glPolygonOffset( 1.0, 1.0 );

	  glBegin(GL_TRIANGLES);


	 tVertex v = vertices;
	 do{
		 if( v->fixed ) 
		 {
			 v = v->next;
			 continue;
		 }
		tHalfEdge h = v->halfedge;
		std::vector<CPoint> pts;
		do{
			tFace f = h->face;
			CPoint q(f->v[X],f->v[Y],f->v[Z] );
			pts.push_back( q );
	
			tEdge e = h->edge;
			h = (e->hes[0]!= h)? e->hes[0]:e->hes[1];
			h = h->he_prev;

		}while( h!= v->halfedge );

		CPoint n = (pts[1]-pts[0])^(pts[2]-pts[0]);

		for( size_t i = 2; i < pts.size(); i ++ )
		{
			CPoint d = (pts[1]-pts[0])^(pts[i]-pts[0]);
			if( d.norm() > n.norm() ) n = d;
		}

		n /= n.norm();
		glNormal3d( n[X], n[Y], n[Z] );

		  for( size_t i = 1; i < pts.size()-1; i ++ )
		  {
			  glVertex3f( pts[0][X], pts[0][Y], -2 );
			  glVertex3f( pts[i][X], pts[i][Y], -2 );
			  glVertex3f( pts[i+1][X], pts[i+1][Y], -2 );
		  }
		  v = v->next;
	  }while( v != vertices );

	glEnd();
	glFrontFace(GL_CCW); 
	}	

}

void CConvexHull::draw_dual_edges()
{
	glDisable(GL_LIGHTING );
	glColor3f(0,0,0);
	glBegin(GL_LINES);

	tEdge e = edges;
	
	do{
		if( e->next_swap )
		{
			glLineWidth( 2.0 );
			glColor3f(1,0,0);
		}
		else
		{
			glLineWidth( 1.0 );
			glColor3f(0,0,0);
		}

		if( !e->convex )
		{
			//printf("Bad Edge (a,b) (%f,%f) event time %f\n", e->a, e->b, e->event_time );
			glColor3f(0,0,1);
		}
		tFace f0 = e->adjface[0];
		tFace f1 = e->adjface[1];

		if( f0->halfedge->vert->fixed && f0->halfedge->he_next->vert->fixed && f0->halfedge->he_prev->vert->fixed )
		{
			e = e->next;
			continue;
		}

		if( f1->halfedge->vert->fixed && f1->halfedge->he_next->vert->fixed && f1->halfedge->he_prev->vert->fixed )
		{
			e = e->next;
			continue;
		}

        glVertex3d ( f0->v[X], f0->v[Y], f0->v[Z]);
        glVertex3d ( f1->v[X], f1->v[Y], f1->v[Z]);

		glVertex3d ( f0->v[X], f0->v[Y], -2);
		glVertex3d ( f1->v[X], f1->v[Y], -2);
		e = e->next;
	}while( e != edges );

	 glEnd();
	 glEnable(GL_LIGHTING );
	 
};

void CConvexHull::set_target_dual_area()
{
	tVertex v = vertices;
	double s = 0;
	double target_s = 0;

	int vnum = 0;
	do{
		if( !v->fixed )
		{
			s += v->dual_area;
			target_s += v->target_dual_area;
			vnum ++;
		}
		v = v->next;
	}while( v != vertices );

	v = vertices;
	do{
		if( !v->fixed )
		{
			v->target_dual_area *= (s/target_s);
		}
		v = v->next;
	}while( v != vertices );

};


void CConvexHull::one_step()
{
	tVertex v = vertices;
	double error = 0;
	do{
		if( v->fixed  )
		{
			v = v->next;
			continue;
		}
		double d = fabs(v->dual_area - v->target_dual_area);
		error  = (error > d )?error:d;
		v = v->next;
	}while( v!= vertices );
	fprintf(stderr, "Error is %f\n", error );


	v = vertices;
	do{
		if( v->fixed ) 
		{
			v = v->next;
			continue;
		}
		double d = v->target_dual_area - v->dual_area;
		v->v[2] += d * 0.02;
		v = v->next;
	}while( v!= vertices );
	Reset();
	ConstructHull();

};

void CConvexHull::EdgeWeight()
{
	tEdge e = edges;
	do{
		CPoint2 p0( e->adjface[0]->v[X], e->adjface[0]->v[Y] );
		CPoint2 p1( e->adjface[1]->v[X], e->adjface[1]->v[Y] );
		
		CPoint2 q0( e->endpts[0]->v[X], e->endpts[0]->v[Y] );
		CPoint2 q1( e->endpts[1]->v[X], e->endpts[1]->v[Y] );

		e->weight = (p0-p1).norm()/(q0-q1).norm();

		e = e->next;
	}while( e != edges );
};

void CConvexHull::VertexWeight()
{
	tVertex v = vertices;
	do{
		v->weight = 0;
		v = v->next;
	}while( v !=  vertices );

	tEdge e = edges;
	do{
		
		e->endpts[0]->weight += e->weight;
		e->endpts[1]->weight += e->weight;

		e = e->next;
	}while( e != edges );
};


bool CConvexHull::one_step_Newton()
{
	EdgeWeight();
	VertexWeight();

	tVertex v = vertices;
	double error = 0;
	int vnum = 0;
	int vid = 0;
	do{
		if( v->fixed  )
		{
			v = v->next;
			continue;
		}
		double d = fabs(v->dual_area - v->target_dual_area);
		if ( d > error  )
		{
			error = d;
			vid = v->vnum;
		}
		vnum ++;
		v = v->next;
	}while( v!= vertices );
	fprintf(stderr, "Face %d Error is %f\n", vid, error );

	Eigen::SparseMatrix<double> M(vnum,vnum);
	M.setZero();
	std::vector<Eigen::Triplet<double> > M_coefficients;

	v = vertices;
	do{
		if( v->fixed ) 
		{
			v = v->next;
			continue;
		}
		M_coefficients.push_back( Eigen::Triplet<double>( v->vnum-1,v->vnum-1, v->weight ));
		v = v->next;
	}while( v!= vertices );

	tEdge e = edges;
	do{
		tVertex v0 = e->endpts[0];
		tVertex v1 = e->endpts[1];
		if( v0->fixed || v1->fixed )
		{
			e = e->next;
			continue;
		}
		
		M_coefficients.push_back( Eigen::Triplet<double>( v0->vnum-1,v1->vnum-1, -e->weight ));
		e = e->next;
	}while( e!= edges );

	M.setFromTriplets(M_coefficients.begin(), M_coefficients.end());

	Eigen::VectorXd b(vnum);

	v = vertices;
	do{
		if( v->fixed )
		{
			v = v->next;
			continue;
		}
		b(v->vnum-1) = v->target_dual_area - v->dual_area;
		v = v->next;
	}while( v != vertices );




		Eigen::SimplicialLDLT<Eigen::SparseMatrix<double>> solver;
		std::cerr << "Eigen Decomposition" << std::endl;
		solver.compute(M);
		std::cerr << "Eigen Decomposition Finished" << std::endl;
	
		if( solver.info() != Eigen::Success )
		{
			std::cerr << "Waring: Eigen decomposition failed" << std::endl;
			return false;
		}

		Eigen::VectorXd x = solver.solve(b);
		if( solver.info() != Eigen::Success )
		{
			std::cerr << "Waring: Eigen decomposition failed" << std::endl;
			return false;
		}

		v = vertices;
		do{
			if( v->fixed )
			{
				v = v->next;
				continue;
			}
			v->old_height = v->v[2];
			v = v->next;
		}while( v != vertices );

		double step = 0.2;
		while( true )
		{
			v = vertices;
			do{
				if( v->fixed )
				{
					v = v->next;
					continue;
				}
				v->v[2] = v->old_height + x(v->vnum-1)*step;
				v = v->next;
			}while( v != vertices );

			if( vertex_buffer )
			{
				v = vertex_buffer;
				do{
					if( v->fixed )
					{
						v = v->next;
						continue;
					}
					v->v[2] = v->old_height + x(v->vnum-1)*step;
					v = v->next;
				}while( v != vertex_buffer );
			}
			Reset();
			if( !ConstructHull() ) break;
			fprintf(stderr, "Hidden vertex\n");
			step /= 2.0;
		}
		return true;
};

double det( double A[3][3] )
{
	return A[0][0]*(A[1][1]*A[2][2]-A[2][1]*A[1][2])-A[0][1]*(A[1][0]*A[2][2]-A[1][2]*A[2][0])+A[0][2]*(A[1][0]*A[2][1]-A[1][1]*A[2][0]);
}

void vertex_neighbor( tVertex v )
{
	printf("Neighborhood of %d\n", v->vnum );

	tHalfEdge h = v->halfedge;

	do{
		printf("Face %d\n", h->face->fnum );
		tEdge e = h->edge;
		tHalfEdge sh = (e->hes[0] != h )?e->hes[0]:e->hes[1];
		h = sh->he_prev;
	}while( h!= v->halfedge );


};

void vertex_clw_neighbor( tVertex v )
{
	printf("Neighborhood of %d\n", v->vnum );

	tHalfEdge h = v->halfedge;

	do{
		printf("Face %d\n", h->face->fnum );
		h = h->he_next;
		tEdge e = h->edge;
		h = (e->hes[0] != h )?e->hes[0]:e->hes[1];
		
	}while( h!= v->halfedge );


};

void CConvexHull::EdgeSwap( tEdge e, double current_time )
{
	tHalfEdge h0 = e->hes[0];
	tHalfEdge h1 = h0->he_next;
	tHalfEdge h2 = h1->he_next;

	tHalfEdge h3 = e->hes[1];
	tHalfEdge h4 = h3->he_next;
	tHalfEdge h5 = h4->he_next;

	tVertex   v0 = h2->vert;
	tVertex   v1 = h0->vert;
	tVertex   v2 = h1->vert;
	tVertex   v3 = h4->vert;

	tEdge     e1 = h1->edge;
	tEdge     e2 = h2->edge;
	tEdge     e4 = h4->edge;
	tEdge     e5 = h5->edge;

	tFace     f0 = h0->face;
	tFace     f1 = h3->face;

	tHalfEdge sh1 = (e1->hes[0] != h1 )? e1->hes[0]: e1->hes[1];
	tHalfEdge sh2 = (e2->hes[0] != h2 )? e2->hes[0]: e2->hes[1];
	tHalfEdge sh4 = (e4->hes[0] != h4 )? e4->hes[0]: e4->hes[1];
	tHalfEdge sh5 = (e5->hes[0] != h5 )? e5->hes[0]: e5->hes[1];


	//vertex_neighbor( v3 );


	//change pointers

	v1->halfedge = sh1;
	v2->halfedge = sh2;
	v3->halfedge = sh5;
	v0->halfedge = sh4;


	



	tHalfEdge th;
	if( e1->hes[1]  != h1 ) SWAP( th, e1->hes[0], e1->hes[1] ); 
	if( e2->hes[1]  != h2 ) SWAP( th, e2->hes[0], e2->hes[1] ); 
	if( e4->hes[1]  != h4 ) SWAP( th, e4->hes[0], e4->hes[1] ); 
	if( e5->hes[1]  != h5 ) SWAP( th, e5->hes[0], e5->hes[1] ); 

	e1->hes[1] = h2;
	e2->hes[1] = h4;
	e4->hes[1] = h5;
	e5->hes[1] = h1;
	
	h2->edge = e1;
	h4->edge = e2;
	h5->edge = e4;
	h1->edge = e5;

	h0->vert = v3;
	h1->vert = v1;
	h2->vert = v2;
	h3->vert = v2;
	h4->vert = v0;
	h5->vert = v3;

	f0->vertex[0] = h2->vert;
	f0->vertex[1] = h0->vert;
	f0->vertex[2] = h1->vert;

	f0->edge[0] = h0->edge;
	f0->edge[1] = h1->edge;
	f0->edge[2] = h2->edge;

	f1->vertex[0] = h5->vert;
	f1->vertex[1] = h3->vert;
	f1->vertex[2] = h4->vert;

	f1->edge[0] = h3->edge;
	f1->edge[1] = h4->edge;
	f1->edge[2] = h5->edge;

	e->endpts[0] = v3;
	e->endpts[1] = v2;

	for( int i = 0; i < 2; i ++ )
	{
		e1->adjface[i] = e1->hes[i]->face;
		e2->adjface[i] = e2->hes[i]->face;
		e4->adjface[i] = e4->hes[i]->face;
		e5->adjface[i] = e5->hes[i]->face;
	}

	//vertex_clw_neighbor( v3 );


	_face_dual_point( f0 );
	_face_dual_point( f1 );

	_vertex_dual_area( v0 );
	_vertex_dual_area( v1 );
	_vertex_dual_area( v2 );
	_vertex_dual_area( v3 );

/*
	v0->delta_h = v0->target_dual_area - v0->dual_area;
	v1->delta_h = v1->target_dual_area - v1->dual_area;
	v2->delta_h = v2->target_dual_area - v2->dual_area;
	v3->delta_h = v3->target_dual_area - v3->dual_area;
*/
	double t;
	/*
	t = EventTime( e );
	t += current_time;
	e->event_time = t;
	*/

	e->event_time = current_time - EPSILON;

	t = EventTime( e1 );
	t += current_time;
	e1->event_time = t;
	
	t = EventTime( e2 );
	t += current_time;
	e2->event_time = t;
	
	t = EventTime( e4 );
	t += current_time;
	e4->event_time = t;
	
	t = EventTime( e5 );
	t += current_time;
	e5->event_time = t;
	
	v0->m_valence --;
	v1->m_valence --;
	v2->m_valence ++;
	v3->m_valence ++;
};


double CConvexHull::EventTime( tEdge e )
{
	tHalfEdge h0 = e->hes[0];
	tHalfEdge h1 = h0->he_next;
	tHalfEdge h2 = h1->he_next;

	tHalfEdge h3 = e->hes[1];
	tHalfEdge h4 = h3->he_next;
	tHalfEdge h5 = h4->he_next;

	tVertex   v0 = h2->vert;
	tVertex   v1 = h0->vert;
	tVertex   v2 = h1->vert;
	tVertex   v3 = h4->vert;

	if( v0->fixed && v1->fixed && v2->fixed && v3->fixed )
	{
		e->a = 0;
		e->b = 1;
		return 1e+20;
	}

	CPoint    p0( v0->v[0], v0->v[1], v0->v[2] );
	CPoint    p1( v1->v[0], v1->v[1], v1->v[2] );
	CPoint    p2( v2->v[0], v2->v[1], v2->v[2] );
	CPoint    p3( v3->v[0], v3->v[1], v3->v[2] );
/*
	double    H0 = (v0->fixed)?0:v0->target_dual_area - v0->dual_area;
	double    H1 = (v1->fixed)?0:v1->target_dual_area - v1->dual_area;
	double    H2 = (v2->fixed)?0:v2->target_dual_area - v2->dual_area;
	double    H3 = (v3->fixed)?0:v3->target_dual_area - v3->dual_area;
*/
	double    H0 = (v0->fixed)?0:v0->delta_h;
	double    H1 = (v1->fixed)?0:v1->delta_h;
	double    H2 = (v2->fixed)?0:v2->delta_h;
	double    H3 = (v3->fixed)?0:v3->delta_h;

	CPoint    q0 = p0 - p3;
	CPoint    q1 = p1 - p3;
	CPoint    q2 = p2 - p3;

	double    d0 = H0 - H3;
	double    d1 = H1 - H3;
	double    d2 = H2 - H3;


	double A[3][3];

	A[0][0] = q0[0]; A[0][1] = q0[1]; A[0][2] = q0[2];
	A[1][0] = q1[0]; A[1][1] = q1[1]; A[1][2] = q1[2];
	A[2][0] = q2[0]; A[2][1] = q2[1]; A[2][2] = q2[2];

	double S = det( A );
	A[0][2] = d0; A[1][2] = d1; A[2][2] = d2;

	double  D = det( A );

	if( fabs( D) < 1e-8 ) return 1e+20;

	double t = -S/D;

	e->a = D;
	e->b = S;

/* debug */
/*
	A[0][2] = q0[2] + t * d0;
	A[1][2] = q1[2] + t * d1;
	A[2][2] = q2[2] + t * d2;

	D = det(A);

	v0->v[2] += t * H0;
	v1->v[2] += t * H1;
	v2->v[2] += t * H2;
	v3->v[2] += t * H3;

	this->VolumeSign( h0->face, v3 );
*/
	return t;
};

void CConvexHull::VertexGradient()
{
	tVertex v = vertices;
	do{
		if( v->fixed  )
		{
			v->delta_h = 0;
			v = v->next;
			continue;
		}
		v->delta_h = v->target_dual_area - v->dual_area;
		v = v->next;
	}while( v!= vertices );

}

double CConvexHull::VertexDualAreaError()
{
	tVertex v = vertices;
	double error = 0;
	do{
		if( v->fixed  )
		{
			v->delta_h = 0;
			v = v->next;
			continue;
		}
		double d = fabs( v->delta_h );
		error  = (error > d )?error:d;
		v = v->next;
	}while( v!= vertices );
	fprintf(stderr, "Error is %f\n", error );
	return error;
}

void CConvexHull::___Compute_Event_Time()
{
	FaceDualPoint();
	VertexDualArea();

	double current_time = 0;

	tVertex v = vertices;
	double error = 0;
	do{
		if( v->fixed  )
		{
			v->delta_h = 0;
			v = v->next;
			continue;
		}
		v->delta_h = v->target_dual_area - v->dual_area;
		double d = fabs( v->delta_h );
		error  = (error > d )?error:d;
		v = v->next;
	}while( v!= vertices );
	fprintf(stderr, "Error is %f\n", error );


	tEdge e = edges;
	do{
		e->next_swap = false;
		e->event_time = EventTime( e );
		//printf("Event time %f\n", e->event_time );
		e = e->next;
	}while( e != edges );

}

tEdge CConvexHull::___next_swap_edge( double current_time )
{
	tEdge swap_edge = NULL;
	double min_t = 1e+20;

	tEdge e = edges;
	do{
		if( e == this->swap_edge )
		{
			e->next_swap = false;
			e = e->next;
			continue;
		}
		e->next_swap = false;
		
		if( e->event_time > current_time )
		{
			if( e->event_time < min_t ) 
			{
				min_t = e->event_time;
				swap_edge = e;
			}
		}
		e = e->next;
	}while( e != edges );

	swap_edge->next_swap = true;
	return swap_edge;
}


void CConvexHull::___Update( double delta_t )
{
	tVertex v = vertices;
	do{
		if( !v->fixed ) 
		{
			v->v[2] += v->delta_h * delta_t;
		}
		v = v->next;
	}while( v != vertices );
}


void CConvexHull::EdgeEventTime()
{
	tEdge e = edges;
	do{
		e->next_swap = false;
		e->event_time = EventTime( e );
		//printf("Event time %f\n", e->event_time );
		e = e->next;
	}while( e != edges );

}

void CConvexHull::kinetic_one_step()
{
/* step 0 */
	FaceDualPoint();
	VertexDualArea();
	VertexGradient();
	EdgeEventTime();
	
	VertexDualAreaError();

	current_time  = 0;

	do{
		/*! step 1 */
		swap_edge = ___next_swap_edge( current_time );

		tHalfEdge h0 = swap_edge->hes[0];
		tHalfEdge h1 = swap_edge->hes[1];

		int i = VolumeSign(h0->face, h1->he_next->vert);

		if( !Convexity() )
		{
			printf("Before edge swap, it is not convex!\n" );
		}

		___Update( swap_edge->event_time - current_time );


		i = VolumeSign(h0->face, h1->he_next->vert);


		tEdge pe = edges;
		do{
			tHalfEdge h0 = pe->hes[0];
			tHalfEdge h1 = pe->hes[1];
			
			int v = VolumeSign( h0->face, h1->he_next->vert );
			if( v < 0 )
			{
				printf("After Update, it is not convex!\n");
			}
			pe = pe->next;
		}while( pe != edges );

		current_time = swap_edge->event_time;
/*! step 3 */
		EdgeSwap( swap_edge, current_time );
		
		if( !Convexity() )
		{
			printf("After Edge Swap, it is not convex!\n");
		}


	}while( current_time < 1.0 );

	FaceDualPoint();


};


void CConvexHull::hull()
{
	Reset();
	ConstructHull();
};


void CConvexHull::VertexValence()
{

	 tVertex v = vertices;
	 do{
		 v->m_valence = 0;
		tHalfEdge h = v->halfedge;
		do{
			v->m_valence ++;	
			tEdge e = h->edge;
			h = (e->hes[0]!= h)? e->hes[0]:e->hes[1];
			h = h->he_prev;

		}while( h!= v->halfedge );

		v = v->next;
	  }while( v != vertices );


};





void CConvexHull::VertexNewton()
{
	EdgeWeight();
	VertexWeight();

	tVertex v = vertices;
	double error = 0;
	int vnum = 0;
	int vid = 0;
	do{
		if( v->fixed  )
		{
			v = v->next;
			continue;
		}
		double d = fabs(v->dual_area - v->target_dual_area);
		if ( d > error  )
		{
			error = d;
			vid = v->vnum;
		}
		vnum ++;
		v = v->next;
	}while( v!= vertices );
	fprintf(stderr, "Face %d Error is %f\n", vid, error );

	Eigen::SparseMatrix<double> M(vnum,vnum);
	M.setZero();
	std::vector<Eigen::Triplet<double> > M_coefficients;

	v = vertices;
	do{
		if( v->fixed ) 
		{
			v = v->next;
			continue;
		}
		M_coefficients.push_back( Eigen::Triplet<double>( v->vnum-1,v->vnum-1, v->weight ));
		v = v->next;
	}while( v!= vertices );

	tEdge e = edges;
	do{
		tVertex v0 = e->endpts[0];
		tVertex v1 = e->endpts[1];
		if( v0->fixed || v1->fixed )
		{
			e = e->next;
			continue;
		}
		
		M_coefficients.push_back( Eigen::Triplet<double>( v0->vnum-1,v1->vnum-1, -e->weight ));
		e = e->next;
	}while( e!= edges );

	M.setFromTriplets(M_coefficients.begin(), M_coefficients.end());

	Eigen::VectorXd b(vnum);

	v = vertices;
	do{
		if( v->fixed )
		{
			v = v->next;
			continue;
		}
		b(v->vnum-1) = v->target_dual_area - v->dual_area;
		v = v->next;
	}while( v != vertices );




		Eigen::SimplicialLDLT<Eigen::SparseMatrix<double>> solver;
		std::cerr << "Eigen Decomposition" << std::endl;
		solver.compute(M);
		std::cerr << "Eigen Decomposition Finished" << std::endl;
	
		if( solver.info() != Eigen::Success )
		{
			std::cerr << "Waring: Eigen decomposition failed" << std::endl;
			return;
		}

		Eigen::VectorXd x = solver.solve(b);
		if( solver.info() != Eigen::Success )
		{
			std::cerr << "Waring: Eigen decomposition failed" << std::endl;
			return;
		}

		double sum = 0;
		for( size_t i = 0; i < vnum; i ++ )
		{
			sum += x(i);
		}
		sum /= vnum;
		for( size_t i = 0; i < vnum; i ++ )
		{
			x(i) -= sum;
		}
		

		v = vertices;
		do{
			if( v->fixed )
			{
				v = v->next;
				continue;
			}
			v->delta_h = x(v->vnum-1);
			v = v->next;
		}while( v != vertices );

};


bool CConvexHull::EdgeConvexity()
{
	bool convex = true;
	tEdge e = edges;
	do{
		tHalfEdge h0 = e->hes[0];
		tHalfEdge h1 = e->hes[1];
		tFace   f = h0->face;
		tVertex v = h1->he_next->vert;
	    int vol = VolumeSign( f, v );

		e->convex = ( vol >= 0 );
		if( vol < 0 ) convex = false;
		e = e->next;
	}while( e!= edges );

	return convex;
};

void CConvexHull::FixEdge( double time )
{
	tEdge e = edges;
	do{
		if( !e->convex )
			EdgeSwap( e, time );
		e = e->next;
	}while( e != edges );
}


};
};

#endif