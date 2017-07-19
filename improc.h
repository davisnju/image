/*!
 * improc.h
 * date: 2017/03/21 15:02
 *
 * author: dav1sNJU
 *
 * brief: 图像处理
 *
 *
*/
#pragma once
#ifndef IMPROC_H
#define IMPROC_H

#include "image.h"

#include "element_lib.h"

typedef struct All_element_context Context;
typedef struct Pixel Pixel;

/*
 *  计算单通道梯度
 *  设置offset来确定梯度计算像素间隔
 *  offset >= 0
 */
void gradientOneChannelOffset(Image_f * imLabel, Image *segOut, float thd, int offset);

/*
 *	分别输入三个通道的图像，执行梯度边缘检测
 *  设置offset来确定梯度计算像素间隔
 *  offset >= 0
 */
void gradientThreeChannelOffset(Image_f * Y, Image_f * Cr, Image_f * Cb, 
                                Image *gradLine, float gradThred, int offset);

/*
 *	根据context计算每个像素点的对比度值，存入mag(Image_f类型)
 */
void gradientCntxt(Context * cntxt, Image_f * mag);


// dissimilarity measure between pixels
float yuvdiff(Image_yuv *ima, int x1, int y1, int x2, int y2);

/* 图像二值化处理 */
void binarizeImage(Image * im, const int width, const int height);

/* 清除孤立像素点 */
void cleanImage(Image * im, const int width, const int height);

/* 腐蚀 */
void thinImage(Image * im, const int width, const int height);

/* 图像膨胀算法 */
void imdilate(Image *binIm, const int width, const int height);

/* 腐蚀函数的三条规则1 */
int G1(int *p);
/* 腐蚀函数的三条规则2 */
int G2(int *p);
/* 腐蚀函数的三条规则31 */
int G3(int *p);
/* 腐蚀函数的三条规则32 */
int G32(int *p); 

/* 打印二值图像 用来测试 */
//void printImage(Image *im, char * name);

#endif /* !IMPROC_H */
