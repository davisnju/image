/*!
 * segment-image.h
 * date: 2017/03/21 15:03
 *
 * author: dav1sNJU
 *
 * brief: 
 *
 *
 */
#pragma once
#ifndef SEGMENT_IMAGE
#define SEGMENT_IMAGE

#include <stdlib.h>
#include "image.h"
#include "misc.h"
#include "filter.h"
#include "segment_graph.h"
#include "edge.h"

#include "element_lib.h"

typedef struct All_element_context Context;

// random color
rgb random_rgb();
yuv random_yuv();

// dissimilarity measure between pixels
static float diff(Image_f *r, Image_f *g, Image_f *b,
                         int x1, int y1, int x2, int y2);

/*
* Segment an image
*
* Input:
* im: image to segment.
* sigma: to smooth the image.
* c: constant for treshold function.
* min_size: minimum component size (enforced by post-processing stage).
*
* Returns:
* dst: a color image representing the segmentation.
* labelImg: a float image representing the segmentation label.
* num_ccs: number of connected components in the segmentation.
*/
Universe * segment_image_cntxt(Context *cntxt, float sigma, float c, int min_size,
                             Image_yuv *dst, Image_f *labelImg, int *num_ccs);

void gradient(Image_f * imLabel, Image *segOut, float thd);

// ============================ 暂时未使用 =================================
/*
* Segment an image
*
* Input:
* im: image to segment.
* sigma: to smooth the image.
* c: constant for treshold function.
* min_size: minimum component size (enforced by post-processing stage).
*
* Returns:
* dst: a color image representing the segmentation.
* labelImg: a float image representing the segmentation label.
* num_ccs: number of connected components in the segmentation.
*/
Universe * segment_image_rgb(Image_rgb *im, float sigma, float c, int min_size,
                             Image_rgb *dst, Image_f *labelImg, int *num_ccs);

void gradientYUV(Image_yuv * im, Image_f * ed1, Image_f * mag);


#endif /* !SEGMENT_IMAGE */
