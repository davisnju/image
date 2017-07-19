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
#pragma once
#ifndef MISC_H
#define MISC_H

#include <stdint.h>
#include <math.h>
#include <stdbool.h>

#ifndef M_PI
#define M_PI 3.141592653589793
#endif

typedef struct rgb { uint8_t r, g, b; } rgb;
typedef struct yuv { uint8_t  y, u, v; } yuv;

//inline float fabs(const float x) { return (x > 0 ? x : -x); };
double dabs(const double x);

int sign(const int x);
int fsign(const float x) ;
int dsign(const double x);

int square(const int x);
float fsquare(const float x);
double dsquare(const double x);

int bound(const int x, const int min, const int max);
float fbound(const float x, const float min, const float max);
double dbound(const double x, const double min, const double max);

/* return ((x < min) || (x > max)); */
bool check_bound(const int x, const int min, const int max);
bool fcheck_bound(const float x, const float min, const float max);
bool dcheck_bound(const double x, const double min, const double max);
//
//int fvlib_round(float x);
//
//int dvlib_round(double x);
//
//double gaussian(double val, double sigma);

#endif /* !MISC_H */
