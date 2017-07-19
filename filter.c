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
#include "filter.h"
#include "convolve.h" //
#include "imconv.h"   //
#include <stdlib.h>

/* normalize mask so it integrates to one */
static void normalize(double *mask, int len)
{
    //int len = mask.size();
    double sum = 0;
    for (int i = 1; i < len; i++)
    {
        sum += fabs(mask[i]);
    }
    sum = 2 * sum + fabs(mask[0]);
    for (int i = 0; i < len; i++)
    {
        mask[i] /= sum;
    }
};
void normalizeWrapper(double *mask, int len)
{
    normalize(mask, len);
};

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
static void smooth(Image_f *dst, Image_f *src, float sigma)
{
    // make_fgauss
    sigma = (sigma > 0.01F) ? sigma : 0.01F;
    int len = (int)ceil(sigma * WIDTH) + 1;
    double *mask = (double *)malloc(len * sizeof(double));
    for (int i = 0; i < len; i++)
    {
        mask[i] = exp(-0.5*dsquare((double)i / sigma));
    }
    normalize(mask, len);

    Image_f *tmp = (Image_f *)malloc(sizeof(Image_f));
    create_imagef(tmp, src->h, src->w, 1);

    convolve_even(tmp, src, mask, len);
    convolve_even(dst, tmp, mask, len);

    free(mask);
    destroy_imagef(tmp);
    free(tmp);
};

void smoothWrapper(Image_f *dst, Image_f *src, float sigma)
{
    smooth(dst, src, sigma);
};


/* compute laplacian */
static void laplacian(Image_f *dst, Image_f *src)
{
    int width = src->w;
    int height = src->h;

    for (int y = 1; y < height - 1; y++)
    {
        for (int x = 1; x < width - 1; x++)
        {
            float d2x = *(src->data + y * width + x - 1) + *(src->data + y * width + x + 1) -
                2 * *(src->data + y * width + x);
            float d2y = *(src->data + (y - 1) * width + x) + *(src->data + (y + 1) * width + x) -
                2 * *(src->data + y * width + x);
            *(dst->data + y * width + x) = d2x + d2y;
        }
    }
};
void laplacianWrapper(Image_f *dst, Image_f *src)
{
    laplacian(dst, src);
};