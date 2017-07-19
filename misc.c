/*!
* misc.h
* date: 2017/03/21 15:02
*
* author: dav1sNJU
*
* brief:  random stuff
*
*
*/
#include "misc.h"

//inline float fabs(const float x) { return (x > 0 ? x : -x); };
double dabs(const double x) { return (x > 0 ? x : -x); };

int sign(const int x) { return (x < 0 ? -1 : 1); };
int fsign(const float x) { return (x < 0 ? -1 : 1); };
int dsign(const double x) { return (x < 0 ? -1 : 1); };

int square(const int x) { return x*x; };
float fsquare(const float x) { return x*x; };
double dsquare(const double x) { return x*x; };

int bound(const int x, const int min, const int max)
{
    return (x < min ? min : (x > max ? max : x));
};
float fbound(const float x, const float min, const float max)
{
    return (x < min ? min : (x > max ? max : x));
};
double dbound(const double x, const double min, const double max)
{
    return (x < min ? min : (x > max ? max : x));
};

bool check_bound(const int x, const int min, const int max)
{
    return ((x < min) || (x > max));
};
bool fcheck_bound(const float x, const float min, const float max)
{
    return ((x < min) || (x > max));
};
bool dcheck_bound(const double x, const double min, const double max)
{
    return ((x < min) || (x > max));
};
//
//int fvlib_round(float x) { return (int)(x + 0.5F); };
//
//int dvlib_round(double x) { return (int)(x + 0.5); };
//
//double gaussian(double val, double sigma)
//{
//    return exp(-dsquare(val / sigma) / 2) / (sqrt(2 * M_PI)*sigma);
//};
