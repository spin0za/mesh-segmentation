/*!
*      \file Hyperbolic.h
*      \brief Hyperbolic Geometric Functions
*	   \author David Gu
*      \date 11/14/2010
*
*/

#ifndef _MESHLIB_HYPERBOLIC_H_
#define _MESHLIB_HYPEROBLIC_H_


namespace MeshLib{

inline double _sinh( double x )
{
	double e = exp( x );
	return (e - 1.0/e)/2.0;
};

inline double _cosh( double x )
{
	double e = exp( x );
	return (e + 1.0/e)/2.0;
};

inline double _tanh( double x )
{
	double e = exp( x );
	return ( e-1.0/e )/( e + 1.0/e );
};

inline double _asinh( double x )
{
	return log( x + sqrt( 1 + x*x));
};

inline double _acosh( double x )
{
	return 2 * log (sqrt((x+1)/2) + sqrt((x-1)/2));
};

inline double _atanh( double x )
{
	return (log  (1+x) - log (1-x))/2;
};

//Calculate corner angle
inline double _hyperbolic_cosine_law( double a, double b, double c )
{
	double C;
	C = acos( (cosh(a) * cosh(b)-cosh(c))/(sinh(a)*sinh(b)) );
	return C;
};


inline double _hyperbolic_inverse_cosine_law( double a, double b, double C )
{
  double c;
  c = cosh(a) * cosh(b) - sinh(a)*sinh(b) * cos(C) ;
  c = log( c+sqrt(c*c-1) );
  return c;
}


}; //namespace

#endif
