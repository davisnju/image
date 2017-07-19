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
#include "imutil.h"


/* compute minimum and maximum value in an uint8_t image */
void min_max(Image *im, uint8_t *ret_min, uint8_t *ret_max)
{
    int width = im->w;
    int height = im->h;

    uint8_t min = imRef(im, 0, 0);
    uint8_t max = imRef(im, 0, 0);
    for (int y = 0; y < height; y++)
    {
        for (int x = 0; x < width; x++)
        {
            uint8_t val = imRef(im, x, y);
            if (min > val)
                min = val;
            if (max < val)
                max = val;
        }
    }

    *ret_min = min;
    *ret_max = max;
};

void float_min_max(Image_f *im, float *ret_min, float *ret_max)
{
    int width = im->w;
    int height = im->h;

    float min = imRef(im, 0, 0);
    float max = imRef(im, 0, 0);
    for (int y = 0; y < height; y++)
    {
        for (int x = 0; x < width; x++)
        {
            float val = imRef(im, x, y);
            if (min > val)
                min = val;
            if (max < val)
                max = val;
        }
    }

    *ret_min = min;
    *ret_max = max;
};

void long_min_max(Image_long *im, long *ret_min, long *ret_max)
{
    int width = im->w;
    int height = im->h;

    long min = imRef(im, 0, 0);
    long max = imRef(im, 0, 0);
    for (int y = 0; y < height; y++)
    {
        for (int x = 0; x < width; x++)
        {
            long val = imRef(im, x, y);
            if (min > val)
                min = val;
            if (max < val)
                max = val;
        }
    }

    *ret_min = min;
    *ret_max = max;
};

void short_min_max(Image_short *im, short *ret_min, short *ret_max)
{
    int width = im->w;
    int height = im->h;

    short min = imRef(im, 0, 0);
    short max = imRef(im, 0, 0);
    for (int y = 0; y < height; y++)
    {
        for (int x = 0; x < width; x++)
        {
            short val = imRef(im, x, y);
            if (min > val)
                min = val;
            if (max < val)
                max = val;
        }
    }

    *ret_min = min;
    *ret_max = max;
};

/* threshold image */
Image *thresholdImage(Image *dst, Image *src, int t)
{
    int width = src->w;
    int height = src->h;

    for (int y = 0; y < height; y++)
    {
        for (int x = 0; x < width; x++)
        {
            imRef(dst, x, y) = (imRef(src, x, y) >= t);
        }
    }

    return dst;
};