/*!
 * imconv.h
 * date: 2017/03/21 15:02
 *
 * author: dav1sNJU
 *
 * brief: image conversion
 *
 *
*/
#pragma once
#ifndef IMCONV_H
#define IMCONV_H

#include "image.h"
#include "imutil.h"
#include "misc.h"

#include <limits.h>

// RGBtoGRAY权重
#define RED_WEIGHT  0.299
#define GREEN_WEIGHT    0.587
#define BLUE_WEIGHT 0.114

// 接口
void imageRGBtoYUV(Image_rgb *input, Image_yuv *output);

void imageYUVtoRGB(Image_yuv *input, Image_rgb *output);

void imageRGBtoGRAY(Image_rgb *input, Image *output);

void imageGRAYtoRGB(Image *input, Image_rgb *output);

void imageUINT8_TtoFLOAT(Image *input, Image_f *output);

void imageINTtoFLOAT(Image_int *input, Image_f *output);

void imageFLOATtoUINT8_T(Image_f *input, Image *output
                         , float min, float max);

void imageFLOATtoUINT8_T2(Image_f *input, Image *output);

void imageUINT8_TtoLONG(Image *input, Image_long *output);

void imageLONGtoUINT8_T(Image_long *input, Image *output, long min, long max);

void imageLONGtoUINT8_T2(Image_long *input, Image *output);

void imageSHORTtoUINT8_T(Image_short *input, Image *output,
                         short min, short max);

void imageSHORTtoUINT8_T2(Image_short *input, Image *output);

// 静态函数只能在该头文件对应的c中调用
static void s_imageRGBtoYUV(Image_rgb *input, Image_yuv *output);

static void s_imageYUVtoRGB(Image_yuv *input, Image_rgb *output);

static void s_imageRGBtoGRAY(Image_rgb *input, Image *output);

static void s_imageGRAYtoRGB(Image *input, Image_rgb *output);

static void s_imageUINT8_TtoFLOAT(Image *input, Image_f *output);

static void s_imageINTtoFLOAT(Image_int *input, Image_f *output);

static void s_imageFLOATtoUINT8_T(Image_f *input, Image *output
                                  , float min, float max);

static void s_imageFLOATtoUINT8_T2(Image_f *input, Image *output);

static void s_imageUINT8_TtoLONG(Image *input, Image_long *output);

static void s_imageLONGtoUINT8_T(Image_long *input, Image *output, long min, long max);

static void s_imageLONGtoUINT8_T2(Image_long *input, Image *output);

static void s_imageSHORTtoUINT8_T(Image_short *input, Image *output,
                                  short min, short max);

static void s_imageSHORTtoUINT8_T2(Image_short *input, Image *output);

#endif /* !IMCONV_H */
