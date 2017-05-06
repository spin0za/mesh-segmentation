/*!
*      \file HyperbolicCircle.h
*      \brief Planar circle
*	   \author David Gu
*      \date 10/19/2010
*
*/

#ifndef _MESHLIB_HYPERBOLIC_CIRCLE_H_
#define _MESHLIB_HYPERBOLIC_CIRCLE_H_

#include "circle.h"

namespace MeshLib{

/*!
*	\brief CCircle class, circle on the plane
*
*	Circle on the two dimensional plane
*/
class CHyperbolicCircle
{
	public:
	/*!
	*	CHyperbolicCircle default constructor, it is initialized to be ((0,0),1)
	*/
    CHyperbolicCircle(){ m_c = CPoint2(0,0); m_r = 1.0; };
	/*!
	*	CHyperbolicCircle class copy operator
	*/
    CHyperbolicCircle( const CPoint2 & center, const double r ) 
	{
		double norm = center[0]*center[0]+center[1]*center[1];
		double c = (exp(r)-1.0)/(exp(r)+1.0);
		c = c * c;
		double a = 1 - c * norm;

		double B = 2.0 * ( c - 1.0) * center[0]/a;
		double C = 2.0 * ( c - 1.0) * center[1]/a;
		double D = (norm - c)/a;

		m_c[0] = -B/2.0;
		m_c[1] = -C/2.0;
		m_r    = sqrt( B*B/4.0+C*C/4.0-D);
	};	
	
	/*!
	 *	CHyperbolicCircle class destructor
	 */
    ~CHyperbolicCircle(){};
	

	/*! The center of the circle
	 *
	 */
	CPoint2 & c() { return m_c; };

	/*! The radius of the circle
	 *
	 */
	double & r() { return  m_r; };
	/*!
	 *	Intersection between two hyperbolic circles
	 */
	int intersect( CHyperbolicCircle & circle, CPoint2 & p0,  CPoint2 & p1 );

	private:
		/*!
		* Center
		*/
		CPoint2 m_c;
		/*!
		 * radius
		 */
		double m_r;
		
};



inline int CHyperbolicCircle::intersect(CHyperbolicCircle & circle , CPoint2 & p0, CPoint2 & p1 )
{
	CCircle c0( m_c, m_r );
	CCircle c1( circle.c(), circle.r() );
	return _circle_circle_intersection( c0, c1, p0, p1);

};


}; //namespace

#endif
