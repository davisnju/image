/*!
 * filter.h
 * date: 2017/03/21 14:48
 *
 * author: dav1sNJU
 *
 * brief:  simple filters 
 *
 *
*/
#pragma once

#ifndef FILTER_H
#define FILTER_H

#include "image.h"
#include <math.h>

// 滤波器带宽
#define WIDTH 4.0

/* normalize mask so it integrates to one */
static void normalize(double *mask, int len);
void normalizeWrapper(double *mask, int len);

/* make filters */
//#define MAKE_FILTER(name, fun)                                \
//static std::vector<float> make_ ## name (float sigma) {       \
//  sigma = max(sigma, 0.01F);                                \
//  int len = (int)ceil(sigma * WIDTH) + 1;                     \
//  std::vector<float> mask(len);                               \
//  for (int i = 0; i < len; i++) {                             \
//    mask[i] = fun;                                            \
//  }                                                           \
//  return mask;                                                \
//}
//
//MAKE_FILTER(fgauss, exp(-0.5*square(i / sigma)));

/* convolve image with gaussian filter */
// 默认FLOATtoFLOAT, Image_f ---> Image_f
static void smooth(Image_f *dst, Image_f *src, float sigma);
void smoothWrapper(Image_f *dst, Image_f *src, float sigma);

/* convolve image with gaussian filter */
static void smoothUINT8_TtoFLOAT(Image *dst, Image_f *src, float sigma);
void smoothUINT8_TtoFLOATWrapper(Image *dst, Image_f *src, float sigma);

/* compute laplacian */
static void laplacian(Image_f *dst, Image_f *src);
void laplacianWrapper(Image_f *dst, Image_f *src);

#endif /* !FILTER_H */