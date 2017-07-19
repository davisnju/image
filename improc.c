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
#include "improc.h"
#include <stdio.h>
#include <string.h>
#include <math.h>
#include <assert.h>

/*
*	计算单通道梯度
*  设置offset来确定梯度计算像素间隔
*  offset >= 0
*/
void gradientOneChannelOffset(Image_f *imLabel, Image *segOut, float thd, int offset)
{
    assert(imLabel);
    assert(segOut);
    assert(offset >= 0);
    int w = imLabel->w, h = imLabel->h;
    float fx = 0, fy = 0;
    for (int y = 0; y < h; ++y)
    {
        for (int x = 0; x < w; ++x)
        {
            fx = 0.0f;
            fy = 0.0f;
            // dx
            if (x <= offset)
                fx = (float)(imRef(imLabel, 1 + offset, y) - imRef(imLabel, 0, y));
            else if (x >= w - 1 - offset)
                fx = (float)(imRef(imLabel, x, y) - imRef(imLabel, x - 1 - offset, y));
            else
                fx = 0.5f * (imRef(imLabel, x + 1 + offset, y) - imRef(imLabel, x - 1 - offset, y));

            // dy
            if (y <= offset)
                fy = (float)(imRef(imLabel, x, 1 + offset) - imRef(imLabel, x, 0));
            else if (y >= h - 1 - offset)
                fy = (float)(imRef(imLabel, x, y) - imRef(imLabel, x, y - 1 - offset));
            else
                fy = 0.5f * (imRef(imLabel, x, y + 1 + offset) - imRef(imLabel, x, y - 1 - offset));

            /* 提取边界线 */
            if ((float)sqrt(fsquare(fx) + fsquare(fy)) > thd)// 计算2-范数
                imRef(segOut, x, y) = 1u;
        }
    }

};


static void getThreeChannelOffsetFxFy(Image_f * fxPy, Image_f * fxPcr, Image_f * fxPcb,
    Image_f * fyPy, Image_f * fyPcr, Image_f * fyPcb,
    Image_f * Y, Image_f * Cr, Image_f * Cb, int offset)
{
    int w = Y->w, h = Y->h;
    for (int y = 0; y < h; y++)
    {
        int roffset = y * w;
        for (int x = 0; x < w; x++)
        {
            //fx
            if (x <= offset)
            {
                imRef(fxPy, x, y) = imRef(Y, 1 + offset, y) - imRef(Y, 0, y);
                imRef(fxPcr, x, y) = imRef(Cr, 1 + offset, y) - imRef(Cr, 0, y);
                imRef(fxPcb, x, y) = imRef(Cb, 1 + offset, y) - imRef(Cb, 0, y);
            }
            else if (x >= w - 1 - offset)
            {
                imRef(fxPy, x, y) = imRef(Y, x, y) - imRef(Y, x - 1 - offset, y);
                imRef(fxPcr, x, y) = imRef(Cr, x, y) - imRef(Cr, x - 1 - offset, y);
                imRef(fxPcb, x, y) = imRef(Cb, x, y) - imRef(Cb, x - 1 - offset, y);
            }
            else
            {
                imRef(fxPy, x, y) = 0.5f * (imRef(Y, x + 1 + offset, y) - imRef(Y, x - 1 - offset, y));
                imRef(fxPcr, x, y) = 0.5f * (imRef(Cr, x + 1 + offset, y) - imRef(Cr, x - 1 - offset, y));
                imRef(fxPcb, x, y) = 0.5f * (imRef(Cb, x + 1 + offset, y) - imRef(Cb, x - 1 - offset, y));
            }
            //fy
            if (y <= offset)
            {
                imRef(fyPy, x, y) = imRef(Y, x, 1 + offset) - imRef(Y, x, 0);
                imRef(fyPcr, x, y) = imRef(Cr, x, 1 + offset) - imRef(Cr, x, 0);
                imRef(fyPcb, x, y) = imRef(Cb, x, 1 + offset) - imRef(Cb, x, 0);
            }
            else if (y >= h - 1 - offset)
            {
                imRef(fyPy, x, y) = imRef(Y, x, y) - imRef(Y, x, y - 1 - offset);
                imRef(fyPcr, x, y) = imRef(Cr, x, y) - imRef(Cr, x, y - 1 - offset);
                imRef(fyPcb, x, y) = imRef(Cb, x, y) - imRef(Cb, x, y - 1 - offset);
            }
            else
            {
                imRef(fyPy, x, y) = 0.5f * (imRef(Y, x, y + 1 + offset) - imRef(Y, x, y - 1 - offset));
                imRef(fyPcr, x, y) = 0.5f * (imRef(Cr, x, y + 1 + offset) - imRef(Cr, x, y - 1 - offset));
                imRef(fyPcb, x, y) = 0.5f * (imRef(Cb, x, y + 1 + offset) - imRef(Cb, x, y - 1 - offset));
            }
        }
    }
};
/*
*	分别输入三个通道的图像，执行梯度边缘检测
*  设置offset来确定梯度计算像素间隔
*  offset >= 0
*/
void gradientThreeChannelOffset(Image_f * Y, Image_f * Cr, Image_f * Cb,
    Image *gradLine, float gradThred, int offset)
{
    assert(Y);
    assert(Cr);
    assert(Cb);
    assert(gradLine);
    assert(offset >= 0);
    int w = Y->w, h = Y->h;
    Image_f *fxPy = create_new_imagef(w, h, 0);
    Image_f *fxPcr = create_new_imagef(w, h, 0);
    Image_f *fxPcb = create_new_imagef(w, h, 0);
    Image_f *fyPy = create_new_imagef(w, h, 0);
    Image_f *fyPcr = create_new_imagef(w, h, 0);
    Image_f *fyPcb = create_new_imagef(w, h, 0);
    Image_f * ed1 = create_new_imagef(w, h, 0);

    getThreeChannelOffsetFxFy(fxPy, fxPcr, fxPcb, fyPy, fyPcr, fyPcb, Y, Cr, Cb, offset);
    float maxed1 = -1;
    for (int y = 0; y < h; y++)
        for (int x = 0; x < w; x++)
        {
            float fx = (float)sqrt(fsquare(imRef(fxPy, x, y)) + fsquare(imRef(fxPcr, x, y)) + fsquare(imRef(fxPcb, x, y)));
            float fy = (float)sqrt(fsquare(imRef(fyPy, x, y)) + fsquare(imRef(fyPcr, x, y)) + fsquare(imRef(fyPcb, x, y)));
            float ed1xy = (float)sqrt(fsquare(fx) + fsquare(fy)); // 计算2-范数
            if (maxed1 < ed1xy) maxed1 = ed1xy;
            imRef(ed1, x, y) = ed1xy;
        }
    /* 提取边界线 */
    for (int y = 0; y < h; y++)
        for (int x = 0; x < w; x++)
            if (imRef(ed1, x, y) / maxed1 > gradThred) 
                imRef(gradLine, x, y) = 1u;

    destroy_imagef(fxPy);
    destroy_imagef(fxPcr);
    destroy_imagef(fxPcb);
    free(fxPy);
    free(fxPcr);
    free(fxPcb);
    destroy_imagef(fyPy);
    destroy_imagef(fyPcr);
    destroy_imagef(fyPcb);
    free(fyPy);
    free(fyPcr);
    free(fyPcb);
    destroy_imagef(ed1);
    free(ed1);
};

void getCntxtFxFy(Context * cntxt, Image_f * fxPy, Image_f * fxPcr, Image_f * fxPcb,
    Image_f * fyPy, Image_f * fyPcr, Image_f * fyPcb)
{
    int w = cntxt->pic_col, h = cntxt->pic_row;
    Pixel *pp1 = NULL, *pp2 = NULL;
    for (int y = 0; y < h; y++)
    {
        int roffset = y * w;
        for (int x = 0; x < w; x++)
        {
            if (x == 0)
                pp1 = &(cntxt->pic[1 + roffset]), pp2 = &(cntxt->pic[0 + roffset]);
            else if (x == w - 1)
                pp1 = &(cntxt->pic[x + roffset]), pp2 = &(cntxt->pic[x - 1 + roffset]);
            else
                pp1 = &(cntxt->pic[x + 1 + roffset]), pp2 = &(cntxt->pic[x - 1 + roffset]);
            if (x == 0 || x == w - 1)
            {
                imRef(fxPy, x, y) = (float)(pp1->Y - pp2->Y);
                imRef(fxPcr, x, y) = (float)(pp1->Cr - pp2->Cr);
                imRef(fxPcb, x, y) = (float)(pp1->Cb - pp2->Cb);
            }
            else 
            {
                imRef(fxPy, x, y) = 0.5f *((float)(pp1->Y - pp2->Y));
                imRef(fxPcr, x, y) = 0.5f *((float)(pp1->Cr - pp2->Cr));
                imRef(fxPcb, x, y) = 0.5f *((float)(pp1->Cb - pp2->Cb));
            }
            if (y == 0)
                pp1 = &(cntxt->pic[x + w]), pp2 = &(cntxt->pic[x + roffset]);
            else if (y == h - 1)
                pp1 = &(cntxt->pic[x + roffset]), pp2 = &(cntxt->pic[x + roffset - w]);
            else
                pp1 = &(cntxt->pic[x + w + roffset]), pp2 = &(cntxt->pic[x - w + roffset]);
            if (y == 0 || y == h - 1)
            {
                imRef(fyPy, x, y) = (float)(pp1->Y - pp2->Y);
                imRef(fyPcr, x, y) = (float)(pp1->Cr - pp2->Cr);
                imRef(fyPcb, x, y) = (float)(pp1->Cb - pp2->Cb);
            }
            else
            {
                imRef(fyPy, x, y) = 0.5f *((float)(pp1->Y - pp2->Y));
                imRef(fyPcr, x, y) = 0.5f *((float)(pp1->Cr - pp2->Cr));
                imRef(fyPcb, x, y) = 0.5f *((float)(pp1->Cb - pp2->Cb));
            }
        }
    }
};
/*
 *	根据context计算每个像素点的对比度值
 */
void gradientCntxt(Context * cntxt, Image_f * mag)
{
    assert(cntxt);
    assert(mag);
    int w = cntxt->pic_col, h = cntxt->pic_row;
    Image_f *fxPy = create_new_imagef(w, h, 1);
    Image_f *fxPcr = create_new_imagef(w, h, 1);
    Image_f *fxPcb = create_new_imagef(w, h, 1);
    Image_f *fyPy = create_new_imagef(w, h, 1);
    Image_f *fyPcr = create_new_imagef(w, h, 1);
    Image_f *fyPcb = create_new_imagef(w, h, 1);
    Image_f * ed1 = create_new_imagef(w, h, 1);

    getCntxtFxFy(cntxt, fxPy, fxPcr, fxPcb, fyPy, fyPcr, fyPcb);
    //printf("\n");
    //for (int y = 0; y < h; y++)
    //{
    //    for (int x = 0; x < w; x++)
    //    {
    //        printf("%4.1f ", imRef(fxPy, x, y));
    //    }
    //    printf("\n");
    //}

    float maxed1 = -1.0f;
    for (int x = 0; x < w; x++)
    {
        for (int y = 0; y < h; y++)
        {
            float ed1fx = (float)sqrt(fsquare(imRef(fxPy, x, y)) + fsquare(imRef(fxPcr, x, y)) + fsquare(imRef(fxPcb, x, y)));
            float ed1fy = (float)sqrt(fsquare(imRef(fyPy, x, y)) + fsquare(imRef(fyPcr, x, y)) + fsquare(imRef(fyPcb, x, y)));
            float ed1xy = (float)sqrt(fsquare(ed1fx) + fsquare(ed1fy));
            imRef(ed1, x, y) = ed1xy;
            if (ed1xy > maxed1)  maxed1 = ed1xy;
        }
    }
    if (fabs(maxed1) > 0.0f)
    {
        for (int x = 0; x < w; x++)
            for (int y = 0; y < h; y++)
                imRef(mag, x, y) = imRef(ed1, x, y) / maxed1;
    }
    else
    {
        for (int x = 0; x < w; x++)
            for (int y = 0; y < h; y++)
                imRef(mag, x, y) = 0;
    }

    destroy_imagef(fxPy);
    destroy_imagef(fxPcr);
    destroy_imagef(fxPcb);
    free(fxPy);
    free(fxPcr);
    free(fxPcb);
    destroy_imagef(fyPy);
    destroy_imagef(fyPcr);
    destroy_imagef(fyPcb);
    free(fyPy);
    free(fyPcr);
    free(fyPcb);
    destroy_imagef(ed1);
    free(ed1);
};

/* 图像二值化处理 */
void binarizeImage(Image * im, const int w, const int h)
{
    for (int y = 0; y < h; y++)
    {
        for (int x = 0; x < w; x++)
        {
            if (imRef(im, x, y) > 0u)
                imRef(im, x, y) = 1u;
        }
    }
};


void cleanImage(Image * im, const int w, const int h)
{
    int roff[] = { -1, 0, 1, -1, 1, -1, 0, 1 },
        coff[] = { -1, -1, -1, 0, 0, 1, 1, 1 };
    for (int y = 0; y < h; y++)
    {
        for (int x = 0; x < w; x++)
        {
            if (imRef(im, x, y) == 0u)
                continue;

            int isalated = 1;
            for (int i = 0; i < 8; i++)
            {
                int c = x + coff[i], r = y + roff[i];
                if (r >= 0 && r < h && c >= 0 && c < w && imRef(im, c, r))
                {
                    isalated = 0;
                    break;
                }
            }
            if (isalated == 1)
                imRef(im, x, y) = 0u;
        }
    }
};


/* 
 * 图像腐蚀算法
 * 像素邻居编号如下:
 *     6  7  8
 *     5  x  1
 *     4  3  2
 */
void thinImage(Image *im, const int w, const int h)
{
    int pixels[8] = { 0 };
    int roff[] = { 0, 1, 1, 1, 0, -1, -1, -1 },
        coff[] = { 1, 1, 0, -1, -1, -1, 0, 1 };
    // 迭代 直到图像不再变化
    int oddchanged = 1, evenchanged = 1;
    while (oddchanged || evenchanged)
    {
        oddchanged = 0;        
        for (int y = 0; y < h; y++)// the first subiteration
            for (int x = 0; x < w; x++)
            {
                if (imRef(im, x, y) == 0u) continue;
                for (int i = 0; i < 8; i++)
                {
                    int c = x + coff[i], r = y + roff[i];
                    if (r >= 0 && r < h && c >= 0 && c < w)
                        pixels[i] = imRef(im, c, r);
                    else
                        pixels[i] = 0;
                }
                if (G1(pixels) && G2(pixels) && G3(pixels))
                {
                    imRef(im, x, y) = 0u;
                    oddchanged = 1;
                }
            }// end of the first subiteration
        evenchanged = 0;
        for (int y = 0; y < h; y++)// the second subiteration
            for (int x = 0; x < w; x++)
            {
                if (imRef(im, x, y) == 0u) continue;
                for (int i = 0; i < 8; i++)
                {
                    int c = x + coff[i], r = y + roff[i];
                    if (r >= 0 && r < h && c >= 0 && c < w)
                        pixels[i] = imRef(im, c, r);
                    else
                        pixels[i] = 0;
                }
                if (G1(pixels) && G2(pixels) && G32(pixels))
                {
                    imRef(im, x, y) = 0u;
                    evenchanged = 1;
                }
            }// end of the second subiteration
    }// end of while()
};

int G1(int *p)
{
    int b0, b1, b2, b3;
    b0 = p[0] == 0 && (p[1] == 1 || p[2] == 1);
    b1 = p[2] == 0 && (p[3] == 1 || p[4] == 1);
    b2 = p[4] == 0 && (p[5] == 1 || p[6] == 1);
    b3 = p[6] == 0 && (p[7] == 1 || p[0] == 1);
    int ret = b0 + b1 + b2 + b3;
    return ret == 1;
};
int G2(int *p)
{
    int minnp = 0;
    int np1 = (p[0] || p[1]) + (p[2] || p[3]) + (p[4] || p[5]) + (p[6] || p[7]),
        np2 = (p[1] || p[2]) + (p[3] || p[4]) + (p[5] || p[6]) + (p[0] || p[7]);
    minnp = min(np1, np2);
    return minnp >= 2 && minnp <= 3;
};
int G3(int *p)
{
    int ret = ((p[6] || p[7] || (!p[1])) && p[0] == 0) || ((p[1] || p[2] || (!p[7])) && p[0] == 0);
    return ret;
};
int G32(int *p)
{
    int ret = ((p[2] || p[3] || (!p[5])) && p[4] == 0) || ((p[5] || p[6] || (!p[3])) && p[4] == 0);
    return ret;
};


// dissimilarity measure between pixels
float yuvdiff(Image_yuv *ima,
              int x1, int y1, int x2, int y2)
{
    return (float)sqrt(square((int)imRef(ima, x1, y1).y - (int)imRef(ima, x2, y2).y) +
                       square((int)imRef(ima, x1, y1).u - (int)imRef(ima, x2, y2).u) +
                       square((int)imRef(ima, x1, y1).v - (int)imRef(ima, x2, y2).v));

};

/* 膨胀 sd_disk(2) */
void imdilate(Image *binIm, const int width, const int height)
{
    int se_disk_r[] = { -2, -1, -1, -1, 0, 0, 0, 0, 1, 1, 1, 2 },
        se_disk_c[] = { 0, -1, 0, 1, -2, -1, 1, 2, -1, 0, 1, 0 };
    Image *copy = (Image *)malloc(sizeof(Image));
    create_image(copy, width, height, 1);
    copy_image(copy, binIm);

    int ri = 0, ci = 0;
    for (int y = 0; y < height; y++)
    {
        for (int x = 0; x < width; x++)
        {
            if (imRef(copy, x, y) == 1u)
            {
                for (int i = 0; i < 12; i++)
                {
                    ri = y + se_disk_r[i];
                    ci = x + se_disk_c[i];
                    if (ri >= 0 && ri < height && ci >= 0 && ci < width)
                        imRef(binIm, ci, ri) = 1u;
                }
            }
        }
    }

    destroy_image(copy);
    free(copy);
};

//
//void printImage(Image *im, char * name)
//{
//    int w = im->w, h = im->h;
//    printf("%s\n", name);
//    for (int y = 0; y < h; y++)
//    {
//        for (int x = 0; x < w; x++)
//        {
//            printf("%d", imRef(im, x, y));
//            if (x < w - 1)putchar(' ');
//        }
//        putchar('\n');
//    }
//
//};
