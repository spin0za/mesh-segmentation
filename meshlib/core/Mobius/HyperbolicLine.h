/*!
*      \file HyperbolicLine.h
*      \brief Hyperbolic Line
*	   \author David Gu
*      \date 04/10/2011
*
*/

#ifndef _MESHLIB_HYPERBOLIC_LINE_H_
#define _MESHLIB_HYPERBOLIC_LINE_H_

#include <assert.h>
#include <math.h>
#include <complex>

namespace MeshLib{

	enum HyperbolicLineType {CIRCLE, LINE};

/*!
*	\brief CHyperbolicLine  class, Hyperbolic Lines
*/
class CHyperbolicLine
{
public:
	/*!
	 *	CHyperbolic Line constructor
	 */
	CHyperbolicLine( std::complex<double> z1, std::complex<double> z2 )
	{
		//rotate m_p[0] ccwly to reach m_p[1]

		if( z1.real() * z2.imag() - z1.imag() * z2.real() > 0 )
		{
			m_p[0] = z1;
			m_p[1] = z2;
		}
		else
		{
			m_p[0] = z2;
			m_p[1] = z1;
		}

		std::complex<double> mid = (m_p[0]+m_p[1])/2.0;

		if( std::abs( mid ) < 1e-10 ) 
		{
			m_type = LINE;
			return;
		}
		else
		{
			m_type = CIRCLE;
		}

		std::complex<double> D = m_p[1] - m_p[0];
		std::complex<double> d( -D.imag(), D.real() );


		std::complex<double> t = ( 1.0 + m_p[0] * std::conj( m_p[0] ) - mid * std::conj( m_p[0] ) - std::conj( mid ) * m_p[0] )/( d * std::conj( m_p[0] ) + std::conj(d) * m_p[0] );

		m_c = mid + d * t.real();
		
		m_r = std::abs( m_p[0] - m_c );

		printf("%f %f\n", m_r - std::abs( m_p[1] - m_c ), std::norm( m_c ) -1 - m_r * m_r);
	}
	/*!
	 *	CHyperbolicLine destructor
	 */
	~CHyperbolicLine(){};

	/*!
	 *	Center of the circle
	 */
	std::complex<double> center() { return m_c; };
	/*!
	 *	Radius of the circle
	 */
	double radius() { return m_r; };
	/*!
	 *	return starting, end point of the hyperbolic line
	 */

	std::complex<double> operator[](int k ) 
	{
		return m_p[k];
	}
protected:
	/*!
	 *	Two points on the hyperbolic line
	 */
	std::complex<double>  m_p[2];
	/*!
	 *	center of the circle
	 */
	std::complex<double>  m_c;
	/*!
	 *	radius of the circle
	 */
	double m_r;
	/*!
	 *	type of the line, which is either a circle or a line
	 */
	HyperbolicLineType m_type;
};



}//name space MeshLib

#endif //_MESHLIB_HYPERBOLIC_LINE_H_ defined