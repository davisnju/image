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
#include "image.h"
#include <stdlib.h>
#include <string.h>
#include <assert.h>

int point_cmpfunc(const void *a, const void *b)
{
    const Point *pa = a, *pb = b;

    if (pa->y < pb->y || (pa->y == pb->y && pa->x < pb->x)) return -1;
    else if (pa->y > pb->y || (pa->y == pb->y && pa->x > pb->x)) return 1;
    else return 0;
};

/* create a new uint8_t image */
Image *create_new_image(const int width, const int height, const int init)
{
    Image *ret = (Image *)malloc(sizeof(Image));

    create_image(ret, width, height, init);

    return ret;
};

/* create an uint8_t image */
void create_image(Image * dst, const int width, const int height, const int init)
{
    assert(dst);
    dst->w = width;
    dst->h = height;
    dst->data = (uint8_t *)malloc(sizeof(uint8_t) * width * height);  // allocate space for image data
    dst->access = (uint8_t **)malloc(sizeof(uint8_t *) * height);   // allocate space for row pointers

    // initialize row pointers
    for (int i = 0; i < height; i++)
    {
        dst->access[i] = dst->data + (i * width /** sizeof(uint8_t)*/);
    }

    if (init && dst->data)
    {
        memset(dst->data, 0, width * height * sizeof(uint8_t));
    }
};

/* delete an uint8_t image */
void destroy_image(Image *im)
{
    if (!im)
    {
        return;
    }
    free(im->data);
    im->data = NULL;
    free(im->access);
    im->access = NULL;
};

/* init an uint8_t image */
void init_image(Image *im, const uint8_t val)
{
    uint8_t *ptr = &(im->access[0][0])/*imPtr(this, 0, 0)*/;
    uint8_t *end = &(im->access[im->h - 1][im->w - 1])/*imPtr(this, im->w - 1, im->h - 1)*/;
    while (ptr <= end)
    {
        *ptr++ = val;
    }
};

/* copy an uint8_t image */
void copy_image(Image *dst, Image *src)
{
    destroy_image(dst);
    create_image(dst, src->w, src->h, 0);
    memcpy(dst->data, src->data, dst->w * dst->h * sizeof(uint8_t));
};

/* get the width of an image. */
int image_width(Image *im) { return im->w; }

/* get the height of an image. */
int image_height(Image *im) { return im->h; }

/* create a new int image */
Image_int *create_new_image_int(const int width, const int height, const int init)
{
    Image_int *ret = (Image_int *)malloc(sizeof(Image_int));

    create_image_int(ret, width, height, init);

    return ret;
};
/* create an int image */
void create_image_int(Image_int * dst, const int width, const int height, const int init)
{
    assert(dst);
    dst->w = width;
    dst->h = height;
    dst->data = (int *)malloc(sizeof(int) * width * height);  // allocate space for image data
    dst->access = (int **)malloc(sizeof(int *) * height);   // allocate space for row pointers

    // initialize row pointers
    for (int i = 0; i < height; i++)
    {
        dst->access[i] = dst->data + (i * width/* * sizeof(int)*/);
    }

    if (init && dst->data)
    {
        memset(dst->data, 0, width * height * sizeof(int));
    }
};

/* delete an int image */
void destroy_image_int(Image_int *im)
{
    if (!im)
    {
        return;
    }
    free(im->data);
    im->data = NULL;
    free(im->access);
    im->access = NULL;
};

// image<float>
/* create an new image and return ptr */
Image_f * create_new_imagef(const int width, const int height, const int init)
{
    Image_f *res = (Image_f *)malloc(sizeof(Image_f));
    create_imagef(res, width, height, init);
    return res;
};

/* create an float image */
void create_imagef(Image_f * dst, const int width, const int height, const int init)
{
    assert(dst);
    dst->w = width;
    dst->h = height;
    dst->data = (float *)malloc(sizeof(float) * width * height);  // allocate space for image data
    dst->access = (float **)malloc(sizeof(float *) * height);   // allocate space for row pointers

    // initialize row pointers
    for (int i = 0; i < height; i++)
    {
        dst->access[i] = dst->data + i * width;
    }

    if (init && dst->data)
    {
        memset(dst->data, 0, width * height * sizeof(float));
    }
};

/* delete an float image */
void destroy_imagef(Image_f *im)
{
    if (!im)
    {
        return;
    }
    free(im->access);
    free(im->data);
    im->data = NULL;
    im->access = NULL;
};

/* init an float image */
void init_imagef(Image_f *im, const float val)
{
    float *ptr = &(im->access[0][0])/*imPtr(this, 0, 0)*/;
    float *end = &(im->access[im->h - 1][im->w - 1])/*imPtr(this, im->w - 1, im->h - 1)*/;
    while (ptr <= end)
    {
        *ptr++ = val;
    }
};

/* copy an float image */
void copy_imagef(Image_f *dst, Image_f *src)
{
    destroy_imagef(dst);
    create_imagef(dst, src->w, src->h, 0);
    memcpy(dst->data, src->data, dst->w * dst->h * sizeof(float));
};

/* get the width of an image. */
int imagef_width(Image_f *im) { return im->w; }

/* get the height of an image. */
int imagef_height(Image_f *im) { return im->h; }

// image<rgb>
/* create an rgb image */
void create_image_rgb(Image_rgb * dst, const int width, const int height, const int init)
{
    assert(dst);
    dst->w = width;
    dst->h = height;
    dst->data = (rgb *)malloc(sizeof(rgb) * width * height);  // allocate space for image data
    dst->access = (rgb **)malloc(sizeof(rgb *) * height);   // allocate space for row pointers

    // initialize row pointers
    for (int i = 0; i < height; i++)
    {
        dst->access[i] = dst->data + i * width;
    }

    if (init && dst->data)
    {
        memset(dst->data, 0, width * height * sizeof(rgb));
    }
};

/* delete an rgb image */
void destroy_image_rgb(Image_rgb *im)
{
    free(im->access);
    free(im->data);
    im->data = NULL;
    im->access = NULL;
};

/* init an rgb image */
void init_image_rgb(Image_rgb *im, const rgb val)
{
    rgb *ptr = &(im->access[0][0])/*imPtr(this, 0, 0)*/;
    rgb *end = &(im->access[im->h - 1][im->w - 1])/*imPtr(this, im->w - 1, im->h - 1)*/;
    while (ptr <= end)
    {
        *ptr++ = val;
    }
};

/* copy an rgb image */
void copy_image_rgb(Image_rgb *dst, Image_rgb *src)
{
    destroy_image_rgb(dst);
    create_image_rgb(dst, src->w, src->h, 0);
    memcpy(dst->data, src->data, dst->w * dst->h * sizeof(rgb));
};

/* get the width of an image. */
int image_rgb_width(Image_rgb *im) { return im->w; }

/* get the height of an image. */
int image_rgb_height(Image_rgb *im) { return im->h; }

// image<yuv>
/* create an yuv image */
void create_image_yuv(Image_yuv * dst, const int width, const int height, const int init)
{
    assert(dst);
    dst->w = width;
    dst->h = height;
    dst->data = (yuv *)malloc(sizeof(yuv) * width * height);  // allocate space for image data
    dst->access = (yuv **)malloc(sizeof(yuv *) * height);   // allocate space for row pointers

    // initialize row pointers
    for (int i = 0; i < height; i++)
    {
        dst->access[i] = dst->data + i * width;
    }

    if (init && dst->data)
    {
        memset(dst->data, 0, width * height * sizeof(yuv));
    }
};

/* delete an yuv image */
void destroy_image_yuv(Image_yuv *im)
{
    if (!im)
    {
        return;
    }
    free(im->access);
    free(im->data);
    im->data = NULL;
    im->access = NULL;
};

/* init an image */
void init_image_yuv(Image_yuv *im, const yuv val)
{
    yuv *ptr = &(im->access[0][0])/*imPtr(this, 0, 0)*/;
    yuv *end = &(im->access[im->h - 1][im->w - 1])/*imPtr(this, im->w - 1, im->h - 1)*/;
    while (ptr <= end)
    {
        *ptr++ = val;
    }
};

/* copy an yuv image */
void copy_image_yuv(Image_yuv *dst, Image_yuv *src)
{
    destroy_image_yuv(dst);
    create_image_yuv(dst, src->w, src->h, 0);
    memcpy(dst->data, src->data, dst->w * dst->h * sizeof(yuv));
};

/* create a new long image */
Image_long *create_new_image_long(const int width, const int height, const int init)
{
    Image_long *imL = (Image_long *)malloc(sizeof(Image_long));
    create_image_long(imL, width, height, init);

    return imL;
};
/* create an long image */
void create_image_long(Image_long * dst, const int width, const int height, const int init)
{
    assert(dst);
    dst->w = width;
    dst->h = height;
    dst->data = (long *)malloc(sizeof(long) * width * height);  // allocate space for image data
    dst->access = (long **)malloc(sizeof(long *) * height);   // allocate space for row pointers

    // initialize row pointers
    for (int i = 0; i < height; i++)
    {
        dst->access[i] = dst->data + (i * width/* * sizeof(long)*/);
    }

    if (init && dst->data)
    {
        memset(dst->data, 0, width * height * sizeof(long));
    }
};

/* delete an long image */
void destroy_image_long(Image_long *im)
{
    if (!im)
    {
        return;
    }
    free(im->data);
    im->data = NULL;
    free(im->access);
    im->access = NULL;
};

/* create a new short image */
Image_short *create_new_image_short(const int width, const int height, const int init)
{
    Image_short *imL = (Image_short *)malloc(sizeof(Image_short));
    create_image_short(imL, width, height, init);

    return imL;
};

/* create a short image */
void create_image_short(Image_short *dst, const int width, const int height, const int init)
{
    assert(dst);
    dst->w = width;
    dst->h = height;
    dst->data = (short *)malloc(sizeof(short) * width * height);  // allocate space for image data
    dst->access = (short **)malloc(sizeof(short *) * height);   // allocate space for row pointers

                                                              // initialize row pointers
    for (int i = 0; i < height; i++)
    {
        dst->access[i] = dst->data + (i * width);
    }

    if (init && dst->data)
    {
        memset(dst->data, 0, width * height * sizeof(short));
    }
}

/* delete a short image */
void destroy_image_short(Image_short *im)
{
    if (!im) 
    {
        return; 
    }
    free(im->data);
    im->data = NULL;
    free(im->access);
    im->access = NULL;
}