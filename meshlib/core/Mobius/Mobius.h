/*!
*      \file Mobius.h
*      \brief Mobius transformation
*	   \author David Gu
*      \date 10/07/2010
*
*/

#ifndef _MESHLIB_MOBIUS_H_
#define _MESHLIB_MOBIUS_H_

#include <assert.h>
#include <math.h>
#include <complex>

namespace MeshLib{

/*!
 *	use complex number in STL
 */
typedef std::complex<double> Complex;

/*!
*	\brief CMobius  class, Mobius transformation
*/
class CMobius
{
public:
	/*!
	 *	CMobius default constructor, z0 is the origin, rotation angle is 0
	 */
	CMobius(){ m_z0=Complex(0,0); m_theta = std::polar(1.0,0.0);};
	/*!
	 *	CMobius constructor, copy operator
	 */
	CMobius( const CMobius & mob ){ m_z0= mob.m_z0; m_theta = mob.m_theta;};
	/*!
	 *	CMobius constructor, with prescribed z0 and rotation angle
	 *  \param z0 center of the Mobius transformation
	 *  \param theta rotation angle
	 */
	CMobius( Complex z0, double theta ) { m_z0 = z0; m_theta = std::polar(1.0,theta);};
	/*!
	 *	Mobius destructor
	 */
	~CMobius(){};
	/*!
	 *	Transform a complex number by the current Mobius transformation
	 *  \param z the input complex number
	 */
	Complex operator*( Complex z ) { return m_theta * (z - m_z0)/(Complex(1.0,0.0)-std::conj(m_z0)*z); };
	/*!
	 *	The z0 of the Mobius transformation
	 */
	Complex & z() { return m_z0; };
	/*!
	 *	The rotation part
	 */
	Complex & theta() { return m_theta; };

protected:
	/*!
	 *	Rotation angle of the Mobius transformation, \f$e^{i\theta}\f$
	 */
	Complex  m_theta;
	/*!
	 *	center of the Mobius transformation
	 */
	Complex  m_z0;
	
	friend CMobius inverse( CMobius & mob );
	friend CMobius operator*( CMobius & mob1, CMobius & mob2 );

};

/*! 
 *	Compute the inverse mobius transformation
 */
CMobius inverse( CMobius & mob )
{
	CMobius imb;
	imb.m_theta = std::conj( mob.m_theta );
	imb.m_z0    = - mob.m_theta * mob.m_z0;
	return imb;
};

/*!
 *	Compose Mobius transformation
 */
CMobius operator*( CMobius & mob1, CMobius & mob2 )
{
	Complex m1[2][2];
	Complex m2[2][2];

	m1[0][0] = mob1.m_theta;
	m1[0][1] =-mob1.m_z0 * mob1.m_theta;
	m1[1][0] =-std::conj( mob1.m_z0 );
	m1[1][1] = 1.0;

	m2[0][0] = mob2.m_theta;
	m2[0][1] =-mob2.m_z0 * mob2.m_theta;
	m2[1][0] =-std::conj( mob2.m_z0 );
	m2[1][1] = 1.0;

	Complex m[2][2];
	for( int i = 0; i < 2; i ++ )
	for( int j = 0; j < 2; j ++ )
	{
		m[i][j] = 0;
		for( int k = 0; k < 2; k ++ )
		{
			m[i][j] += m1[i][k] * m2[k][j];
		}
	};
	
	Complex d = m[1][1];

	for( int i = 0; i < 2; i ++ )
	for( int j = 0; j < 2; j ++ )
	{
		m[i][j] /= d;
	}

	CMobius mob;

	mob.m_theta = m[0][0];
	mob.m_z0    = std::conj( -m[1][0] );

	return mob;
};

}//name space MeshLib

#endif //_MESHLIB_POINT_H_ defiined