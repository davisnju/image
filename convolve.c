/*!
* convolve.h
* date: 2017/03/21 14:45
*
* author: dav1sNJU
*
* brief: convolution
*
*/
#include "convolve.h"
#include "common.h"
#include <stdlib.h>

/* convolve src with mask.  dst is flipped! 
 *
 * @param dst 输出图像
 * @param src 原图像
 * @param mask 数组 长度为len
 */
static void convolve_even_s(Image_f *dst, Image_f *src,
                            double *mask, int len)
{
    int width = src->w;
    int height = src->h;

    for (int y = 0; y < height; y++)
    {
        for (int x = 0; x < width; x++)
        {
            double sum = mask[0] * imRef(src, x, y);
            for (int i = 1; i < len; i++)
            {
                int xli = max(x - i, 0),
                    xgi = min(x + i, width - 1);
                sum += mask[i] *
                    (imRef(src, xli, y)
                    + imRef(src, xgi, y));
            }
            imRef(dst, y, x) = (float)sum;
        }
    }
};

/* convolve src with mask.  dst is flipped! */
static void convolve_odd_s(Image_f *dst, Image_f *src,
                           double *mask, int len)
{
    int width = src->w;
    int height = src->h;

    for (int y = 0; y < height; y++)
    {
        for (int x = 0; x < width; x++)
        {
            double sum = mask[0] * imRef(src, x, y);
            for (int i = 1; i < len; i++)
            {
                int xli = max(x - i, 0),
                    xgi = min(x + i, width - 1);
                sum += mask[i] *
                    (imRef(src, xli, y)
                    - imRef(src, xgi, y));
            }
            imRef(dst, y, x) = (float)sum;
        }
    }
};

void convolve_even(Image_f *dst, Image_f *src,
                   double *mask, int len)
{
    convolve_even_s(dst, src, mask, len);
};
void convolve_odd(Image_f *dst, Image_f *src,
                  double *mask, int len)
{
    convolve_odd_s(dst, src, mask, len);
};
