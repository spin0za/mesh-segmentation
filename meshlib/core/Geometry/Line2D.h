/*!
*      \file Line2D.h
*      \brief Planar Line
*	   \author David Gu
*      \date 01/21/2011
*
*/

#ifndef _MESHLIB_LINE_2D_H_
#define _MESHLIB_LINE_2D_H_


namespace MeshLib{

/*!
*	\brief CLine2D class, line on the plane
*
*	Line on the two dimensional plane <p,n> = d
*/
class CLine2D
{
	public:

	/*!
	*	CLine2D Constructor, 
	*   \param n normal
	*   \param d distance
	*/

    CLine2D( CPoint2 n,  double d){ assert( mag2(n) > 0); m_n = n/mag(n); m_d = d; };

	/*!
	 *	CLine2D class destructor
	 *  \param  p1
	 *  \param  p2
	 */

	CLine2D( const CPoint2 p1, const CPoint2 p2 ) 
	{	
		CPoint2 n = p2-p1; 
		assert( mag2(n) > 0 );

		m_n[0] = -n[1]; 
		m_n[1] =  n[0]; 
		
		m_n = m_n/mag( m_n );
		m_d = p1 * m_n;
	};

	double _evaluate( CPoint2 q )
	{
		return q * m_n - m_d;
	};

	/*!
	*	CLine2D Destructor
	*/
    ~CLine2D(){};
	

	/*! The normal of the line
	 *
	 */
	CPoint2 & n() { return m_n; };

	/*! The distance
	 *
	 */
	double & d() { return  m_d; };

	private:
		/*!
		* normal
		*/
		CPoint2 m_n;
		/*!
		 * distance
		 */
		double m_d;

};

/*!
*	Intersection between 2 CLine2D, 
*   \param L1 first line
*   \param L2 second line
*   \param p  intersection point
*
*   return 0, no intersection
*   return 1, with intersection
*/

inline int _intersect(CLine2D & L1, CLine2D & L2, CPoint2 & p )
{

	double m[2][2];

	m[0][0] = L1.n()[0];
	m[0][1] = L1.n()[1];

	m[1][0] = L2.n()[0];
	m[1][1] = L2.n()[1];
	
	double h[2];
	h[0] = L1.d();
	h[1] = L2.d();

	double det = m[0][0] * m[1][1] - m[1][0] * m[0][1];

	if( fabs( det ) < 1e-18 )
	{
		fprintf(stderr, "parallel lines\n");
		return 0;
	}

	double inv[2][2];

	inv[0][0] =  m[1][1]/det;
	inv[1][1] =  m[0][0]/det;
	inv[0][1] = -m[0][1]/det;
	inv[1][0] = -m[1][0]/det;

	double x[2];
	x[0] = inv[0][0] * h[0] + inv[0][1] * h[1];
	x[1] = inv[1][0] * h[0] + inv[1][1] * h[1];

	p[0] = x[0];
	p[1] = x[1];

	return 1;
};


}; //namespace

#endif
