/*!
 * convolve.h
 * date: 2017/03/21 14:45
 *
 * author: dav1sNJU
 *
 * brief: convolution
 *
*/
#pragma once
#ifndef CONVOLVE_H
#define CONVOLVE_H
#include "misc.h"
#include "image.h"

static void convolve_even_s(Image_f *dst, Image_f *src,
                            double *mask, int len);
static void convolve_odd_s(Image_f *dst, Image_f *src,
                           double *mask, int len);

void convolve_even(Image_f *dst, Image_f *src,
                   double *mask, int len);
void convolve_odd(Image_f *dst, Image_f *src,
                  double *mask, int len);

#endif /* !CONVOLVE_H */
