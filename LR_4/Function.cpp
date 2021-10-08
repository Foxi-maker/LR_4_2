#include "Function.h"

double ExFirstRight(double, double)
{
	return 0.0;
}

double ExFirstGamma_1(double)
{
	return 1.;
}

double ExFirstGamma_2(double)
{
	return 1.;
}

double ExFirstGamma_3(double)
{
	return 1.;
}

double ExFirstGamma_4(double)
{
	return 1.;
}

double ExSecondRight(double x1, double x2)
{
	return 0.0;
}

double ExSecondGamma_1(double)
{
	return -1.;
}

double ExSecondGamma_2(double x2)
{
	return 1. + x2;
}

double ExSecondGamma_3(double x2)
{
	return 1.+x2;
}

double ExSecondGamma_4(double)
{
	return 1.;
}

double ExThirdRight(double, double)
{
	return 4.;
}

double ExThirdGamma_1(double x1)
{
	return x1*x1;
}

double ExThirdGamma_2(double)
{
	return 2.;
}

double ExThirdGamma_3(double)
{
	return 0.0;
}

double ExThirdGamma_4(double x1)
{
	return 1+ x1*x1;
}

double ExNineRight(double x1, double x2)
{
	return (x1*x1+x2*x2)*sin(x1*x2);
}

double ExNineGamma_1(double x1)
{
	return x1;
}

double ExNineGamma_2(double x2)
{
	return -sin(x2);
}

double ExNineGamma_3(double)
{
	return 0.0;
}

double ExNineGamma_4(double x1)
{
	return -x1 * cos(x1);
}
