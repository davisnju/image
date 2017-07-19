/*!
 * image.h
 * date: 2017/03/21 15:01
 *
 * author: dav1sNJU
 *
 * brief: a simple image class
 *
 *
*/
#pragma once
#ifndef IMAGE_H
#define IMAGE_H

#include "misc.h"

//#include "common.h"


typedef int point_val_type; //二维坐标点值类型
typedef point_val_type pv_t;     //二维坐标点值类型

/*  二维坐标点Point p, p.x表示列号, p.y表示行号, 值类型为pv_t;              */
typedef struct point {
    pv_t x; // col
    pv_t y; // row
}Point;

/*
*	Point 比较仿函数
*  比较顺序：先比较行号，再比较列号
*/
int point_cmpfunc(const void *a, const void *b);


/* uint8_t 图像 */
typedef struct image{

    /* image data. */
    uint8_t *data;

    /* row pointers. */
    uint8_t **access;

    /* width, height. */
    int w, h;
}Image;

/* use imRef to access image data. 访问第y(<height)行第x(<width)列元素.  */
#define imRef(im, x, y) (im->access[y][x])

/* use imPtr to get pointer to image data. 取第y(<height)行第x(<width)列元素地址 */
#define imPtr(im, x, y) &(im->access[y][x])

/* create a new uint8_t image */
Image *create_new_image( const int width, const int height, const int init);

/* create an uint8_t image */
void create_image(Image *dst, const int width, const int height, const int init);

/* delete an image */
void destroy_image(Image *im);

/* init an image */
void init_image(Image *im, const uint8_t val);

/* copy an image */
void copy_image(Image *dst, Image *src);

/* get the width of an image. */
int image_width(Image *im);

/* get the height of an image. */
int image_height(Image *im);

/* int 图像 */
typedef struct image_int{

    /* image data. */
    int *data;

    /* row pointers. */
    int **access;

    /* width, height. */
    int w, h;
}Image_int;

/* create a new int image */
Image_int *create_new_image_int(const int width, const int height, const int init);

/* create an image */
void create_image_int(Image_int *dst, const int width, const int height, const int init);

/* delete an image */
void destroy_image_int(Image_int *im);


/* float 图像 */
typedef struct image_f{

    /* image data. */
    float *data;

    /* row pointers. */
    float **access;

    /* width, height. */
    int w, h;
}Image_f;

/* create a new image and return ptr */
Image_f * create_new_imagef(const int width, const int height, const int init);

/* create an image */
void create_imagef(Image_f *dst, const int width, const int height, const int init);

/* delete an image */
void destroy_imagef(Image_f *im);

/* init an image */
void init_imagef(Image_f *im, const float val);

/* copy an image */
void copy_imagef(Image_f *dst, Image_f *src);

/* get the width of an image. */
int imagef_width(Image_f *im);

/* get the height of an image. */
int imagef_height(Image_f *im);

/* rgb 图像 */
typedef struct image_rgb{

    /* image data. */
    rgb *data;

    /* row pointers. */
    rgb **access;

    /* width, height. */
    int w, h;
}Image_rgb;

/* create an image */
void create_image_rgb(Image_rgb *dst, const int width, const int height, const int init);

/* delete an image */
void destroy_image_rgb(Image_rgb *im);

/* init an image */
void init_image_rgb(Image_rgb *im, const rgb val);

/* copy an image */
void copy_image_rgb(Image_rgb *dst, Image_rgb *src);

/* get the width of an image. */
int image_rgb_width(Image_rgb *im);

/* get the height of an image. */
int image_rgb_height(Image_rgb *im);

/* yuv 图像 */
typedef struct image_yuv{

    /* image data. */
    yuv *data;

    /* row pointers. */
    yuv **access;

    /* width, height. */
    int w, h;
}Image_yuv;

/* create a yuv image */
void create_image_yuv(Image_yuv *dst, const int width, const int height, const int init);

/* delete a yuv image */
void destroy_image_yuv(Image_yuv *im);

/* init an image */
void init_image_yuv(Image_yuv *im, const yuv val);

/* copy an image */
void copy_image_yuv(Image_yuv *dst, Image_yuv *src);

// image<long>
typedef struct image_long{

    /* image data. */
    long *data;

    /* row pointers. */
    long **access;

    /* width, height. */
    int w, h;
}Image_long;

/* create a new long image */
Image_long *create_new_image_long(const int width, const int height, const int init);

/* create a long image */
void create_image_long(Image_long *dst, const int width, const int height, const int init);

/* delete a long image */
void destroy_image_long(Image_long *im);

// image<short>
typedef struct image_short{

    /* image data. */
    short *data;

    /* row pointers. */
    short **access;

    /* width, height. */
    int w, h;
}Image_short;

/* create a new short image */
Image_short *create_new_image_short(const int width, const int height, const int init);

/* create a short image */
void create_image_short(Image_short *dst, const int width, const int height, const int init);

/* delete a short image */
void destroy_image_short(Image_short *im);

typedef Image Mat8U;
typedef Image_f MatF;
typedef Image_int Mat;
typedef Image_long MatL;
typedef Image_short MatS;

#endif /* !IMAGE_H */

