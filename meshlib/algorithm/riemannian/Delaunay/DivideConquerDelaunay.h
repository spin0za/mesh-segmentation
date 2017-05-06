/*! \file  DelaunayMesh.h
*   \brief Mesh Generator using Delaunay Triangulation
*   \date  Documented on 10/14/2010
*
*   Generate planar triangular mesh with input PLSG (Piecewise Linear Segment Graph )
*/

#ifndef  _DIVIDE_CONQUEER_CONVEX_HULL_H_
#define  _DIVIDE_CONQUEER_CONVEX_HULL_H_

#include <map>
#include <vector>
#include <queue>
#include <iostream>
#include "Geometry/Point.h"

#define INF 1e99

namespace MeshLib
{

struct Point {
  double x, y, z;
  Point *prev, *next;
  void act() {  
    if (prev->next != this) prev->next = next->prev = this;  // insert
    else { prev->next = next; next->prev = prev; }  // delete
  }
};

/*!
 *	Divide Conquer Algorithm to compute Delaunay Triangulation
 *  based on "A Minimalist's Implementation of the 3-d Divide-and-Conquer Convex Hull Algorithm"
 *  by	Timothy M. Chan
 */

class CHull
{
public:
	CHull(){ nil.x = INF; nil.y = INF; nil.z = INF; nil.prev = 0; nil.next = 0; NIL = &nil;};
	~CHull(){};
	
	/*!
	 *	Compute the Delaunay triangulation
	 *  \param pts input points {(x_i,y_i,x_i^2+y_i^2)}
	 *  \param tris, output triangles {(index1, index2, index3)}
	 *  note that the orientation of the output triangles may not be consistently positive
	 */
	void _delaunay( std::vector<CPoint> & pts, std::vector<int> & triangles );

protected:
	Point *sort(Point P[], int n);
	double turn(Point *p, Point *q, Point *r);
	double time(Point *p, Point *q, Point *r);
	void hull(Point *list, int n, Point **A, Point **B);
protected:
	Point * NIL;
	Point   nil;
};






}
#endif  _DIVIDE_CONQUEER_CONVEX_HULL_H_