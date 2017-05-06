#include "functions.h"

CPoint rgb2yuv(CPoint rgb)
{
	double y = 0.299 * rgb[0] + 0.587 * rgb[1] + 0.114 * rgb[2];
	double u = -0.169 * rgb[0] - 0.331 * rgb[1] + 0.5 * rgb[2] + 0.5;
	double v = 0.5 * rgb[0] - 0.419 * rgb[1] - 0.081 * rgb[2] + 0.5;

	return CPoint(y, u, v);
}

CPoint yuv2rgb(CPoint yuv)
{
	double r = yuv[0] + 1.13983 * (yuv[2] - 0.5);
	double g = yuv[0] - 0.39465 * (yuv[1] - 0.5) - 0.58060 * (yuv[2] - 0.5);
	double b = yuv[0] + 2.03221 * (yuv[1] - 0.5);

	return CPoint(r, g, b);
}

int power(int x, int y)
{
	int z = 1;
	for (int i = 0; i < y; i++)
	{
		z *= x;
	}
	return z;
}

double distance(QPoint p, QPoint q)
{
	return sqrt((p.rx() - q.rx()) * (p.rx() - q.rx()) + (p.ry() - q.ry()) * (p.ry() - q.ry()));
}

double distance(CPoint p, CPoint q)
{
	return (p - q).norm();
}