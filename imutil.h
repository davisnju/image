/*!
 * imutil.h
 * date: 2017/03/21 15:02
 *
 * author: dav1sNJU
 *
 * brief: some image utilities
 *
 *
*/
#pragma once
#ifndef IMUTIL_H
#define IMUTIL_H

#include "image.h"
#include "misc.h"

/* compute minimum and maximum value in an image */
void min_max(Image *im, uint8_t *ret_min, uint8_t *ret_max);

void float_min_max(Image_f *im, float *ret_min, float *ret_max);

void long_min_max(Image_long *im, long *ret_min, long *ret_max);

void short_min_max(Image_short *im, short *ret_min, short *ret_max);

/* threshold image */
Image *thresholdImage(Image *dst, Image *src, int t);

#endif /* !IMUTIL_H */
