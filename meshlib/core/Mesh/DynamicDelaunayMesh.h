/*!
*      \file DelaunayMesh.h
*      \brief Delaunay Mesh
*
*      Dynamic Mesh Class for Delaunay mesh in Generalized Ricci Flow
*	   \author David Gu
*      \date 09/26/2013
*
*/

#ifndef  _DYNAMIC_DELAUNAY_MESH_H_ 
#define  _DYNAMIC_DELAUNAY_MESH_H_

#include <map>
#include <vector>
#include <queue>

#include "Mesh/BaseMesh.h"
#include "Mesh/boundary.h"
#include "Mesh/iterators.h"
#include "Mesh/DynamicMesh.h"


namespace MeshLib
{
#define EPSILON 1e-8

/*-------------------------------------------------------------------------------------------------------------------------------------

	Delaunay Mesh Class

--------------------------------------------------------------------------------------------------------------------------------------*/
/*! \brief CDynamicDelaunayMesh class : Delaunay mesh
* 
*  Delaunay Mesh for Generalized Ricci Flow
*/
template<typename V, typename E, typename F, typename H>
class CDynamicDelaunayMesh : public CDynamicMesh<V,E,F,H>
{
public:
	/*! CDynamicDelaunayMesh constructor */
	CDynamicDelaunayMesh(){};
	/*! CDynamicDelaunayMesh destructor */
	~CDynamicDelaunayMesh(){};
	/*! preserve Delaunay */
	void _preserve_Delaunay();
	/*! verify if a mesh is Delaunay */
	bool _is_Delaunay();
	/*! verify if the edge lengths form a metric, satisfying triangle inequality */
	bool _is_metric();

protected:

	/*! cosine law
	* \param edge lengths a,b,c
	* \param return angles
	*/
	double _cosine_law( double a, double b, double c );
	/*! 
	* from edge length to corner angles for a face f
	*/
	void   _metric_2_angle( F * f );
	/*! 
	* from edge length to corner angles for the whole mesh
	*/
	void   _metric_2_angle();
	/*!
	*	swap edge
	*/
	bool _swap_edge( E * e );
	
	/*! verify if an edge is Delaunay */
	bool _is_Delaunay( E * e, double tolerance = 1e-8);

	/*! verify if the edge lengths form a metric for the current face */
	bool _is_metric( F * f );

};

/*---------------------------------------------------------------------------*/

//Euclidean cosine law
template<typename V, typename E, typename F, typename H>
double CDynamicDelaunayMesh<V,E,F,H>::_cosine_law( double a, double b, double c )
{
          double cs =  ( a*a + b * b  - c * c )/( 2.0 * a * b );
          assert( cs <= 1.0 && cs >= -1.0 );
          return acos( cs );    
};

//convert metric to angle
template<typename V, typename E, typename F, typename H>
void CDynamicDelaunayMesh<V,E,F,H>::_metric_2_angle( F * f )
{

		H *    he[3];
		double lg[3];

		H * h = faceHalfedge( f );
		for( int i = 0; i < 3; i ++ )
		{
			he[i] = h;
			h = faceNextCcwHalfEdge( h );
		}

		for( int i = 0; i < 3; i ++ )
		{
			E * e = halfedgeEdge( he[i] );
			lg[i] = e->length();
		}

		for(int i = 0; i < 3; i ++ )
		{
          he[(i+1)%3]->angle() = _cosine_law( lg[(i+1)%3] ,lg[(i+2)%3] ,lg[i] );
		}

};

//convert metric to angle
template<typename V, typename E, typename F, typename H>
void CDynamicDelaunayMesh<V,E,F,H>::_metric_2_angle()
{

	for( std::list<F*>::iterator fiter = faces().begin(); fiter != faces().end(); fiter ++ )
	{
		F * f = *fiter;
		_metric_2_angle( f );		
	}
};

//convert metric to angle
template<typename V, typename E, typename F, typename H>
bool CDynamicDelaunayMesh<V,E,F,H>::_swap_edge( E * e )
{
		H *    lh[3];
		double ll[3];
		double la[3];

		lh[0] = edgeHalfedge( e, 0 );

		for( int i = 0; i < 2; i ++ )
		{
			lh[i+1] = faceNextCcwHalfEdge( lh[i] );
		}

		for( int i = 0; i < 3; i ++ )
		{
			E * e = halfedgeEdge( lh[i] );
			ll[i] = e->length();
			la[i] = lh[i]->angle();
		}


		H *    rh[3];
		double rl[3];
		double ra[3];

		rh[0] = edgeHalfedge( e, 1 );

		for( int i = 0; i < 2; i ++ )
		{
			rh[i+1] = faceNextCcwHalfEdge( rh[i] );
		}

		for( int i = 0; i < 3; i ++ )
		{
			E * e = halfedgeEdge( rh[i] );
			rl[i] = e->length();
			ra[i] = rh[i]->angle();
		}

		double diagonal = sqrt( ll[1]*ll[1] + rl[2] * rl[2] - 2 * ll[1] * rl[2] * cos( la[0] + ra[2] ) );

		//if( ll[0] > diagonal )
		if( la[1]+ ra[1] > PI + EPSILON)
		{
			std::cerr << "Edge Swap" << std::endl;
			swapEdge( e );
			e->length() = diagonal;
			
			F * f = halfedgeFace( lh[0] );
			_metric_2_angle( f );
			
			f = halfedgeFace( rh[0] );
			_metric_2_angle( f );

			return true;
		}

		return false;
};

//convert metric to angle
template<typename V, typename E, typename F, typename H>
void CDynamicDelaunayMesh<V,E,F,H>::_preserve_Delaunay()
{
	_metric_2_angle();

	while( true )
	{
		bool swapped = false;
		for( std::list<E*>::iterator eiter = edges().begin(); eiter != edges().end(); eiter ++ )
		{
			E * e = *eiter;
			if( _swap_edge(e) ) swapped = true;
		}
		if( !swapped ) break;
	}

	if( !_is_Delaunay() ) 
	{
		std::cerr << "Mesh is not Delaunay " << std::endl;
	}
};


//verify if an edge is Delaunay with tolerance

template<typename V, typename E, typename F, typename H>
bool CDynamicDelaunayMesh<V,E,F,H>::_is_Delaunay( E * e, double tolerance )
{
		H *    lh;
		lh = edgeHalfedge( e, 0 );
		lh = faceNextCcwHalfEdge( lh );

		H *    rh;
		rh = edgeHalfedge( e, 1 );
		rh = faceNextCcwHalfEdge( rh );
		
		return (lh->angle() + rh->angle() < (PI + tolerance));
};

//verify if the whole mesh is Delaunay

template<typename V, typename E, typename F, typename H>
bool CDynamicDelaunayMesh<V,E,F,H>::_is_Delaunay()
{
	_metric_2_angle();

	for( std::list<E*>::iterator eiter = edges().begin(); eiter != edges().end(); eiter ++ )
	{
		E * e = *eiter;
		if( !_is_Delaunay(e  ) ) return false;
	}
	return true;
};

//verify if the triangle inequality holds for the whole mesh

template<typename V, typename E, typename F, typename H>
bool CDynamicDelaunayMesh<V,E,F,H>::_is_metric( F * f)
{
	H * h[3];
	
	h[0] = faceHalfedge( f );
	h[1] = faceNextCcwHalfEdge( h[0] );
	h[2] = faceNextCcwHalfEdge( h[1] );

	double l[3];

	for( int i = 0; i < 3; i ++ )
	{
		E * e = halfedgeEdge( h[i] );
		l[i] = e->length();
	};

	for( int i = 0; i < 3; i ++ )
	{
		if( l[(i+0)%3] + l[(i+1)%3] <= l[(i+2)%3] ) return false;
	}

	return true;
};

//verify if the triangle inequality holds for the whole mesh

template<typename V, typename E, typename F, typename H>
bool CDynamicDelaunayMesh<V,E,F,H>::_is_metric()
{
	for( std::list<F*>::iterator fiter = faces().begin(); fiter != faces().end(); fiter ++ )
	{
		F * f = *fiter;
		if( !_is_metric( f  ) ) return false;
	}
	return true;
};

}
#endif  _DYNAMIC_DELAUNAY_MESH_H_