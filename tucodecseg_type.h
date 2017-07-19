/*!
 * tucodecseg_type.h
 * date: 2017/03/28 10:47
 *
 * author: dav1sNJU
 *
 * brief: 
 *
 * TODO: long description
 *
*/
#pragma once

#ifndef TUCODECSEG_TYPE_H
#define TUCODECSEG_TYPE_H

#include "common.h"
#include "edge.h"
#include "border_band.h"

typedef int LabelVal_t;
typedef int16_t DiffVal_t;

struct Border_seg      //边界信息
{
    int height;    //高度
    int width;     //宽度
    uint8_t direction;            // 1:横向扩展 非1:纵向扩展
    int* r;
    int* c;
    int8_t * delta;   //1 3 4type时，为下一行边界最左像素y坐标减去上一行像素 2则为下一列
    uint8_t * predict_info;    // 参考边带线的像素值
};
typedef struct Border_seg Border_seg;
// 预测信息
typedef struct PredictInfo
{
    int width;
    int col;
    uint8_t * YData;
    uint8_t * CbData;
    uint8_t * CrData;
}PredictInfo;
// 块的边带信息
typedef struct BorderContext
{
    uint8_t block_size;                       //block_size = 64 >> lvl_i; lvl_i=0,1,2,3,4

    LabelVal_t *borderLabelMap;  //块边界标签  // elts_count 4096
    LabelVal_t *flatLabelMap;    //块平坦标签  // elts_count 4096
    
    DiffVal_t *pic_border_diffY;            // elts_count 4096
    DiffVal_t *pic_border_diffCb;           // elts_count 4096
    DiffVal_t *pic_border_diffCr;           // elts_count 4096

    Border_seg **border;                    // elts_count = max_border_num
    LabelVal_t borderLabel0;            // 边带起始标号
    int border_num;
    int max_border_num;          //max_border_num = 1 << (lvl_i << 1);
}BorderContext;

typedef struct SegContext
{
    uint8_t border_processed;    //5张边带图是否已存在
    int block_r;
    int block_c;
    BorderContext border_context[5];

    EdgeLine *edgeList;
    EdgeLine *longEdgeList;
    Image_f *edge_prob_map;
    Image_int *edge_label_map;
    float *edgesContrast;
    int edgesNum;
    int longEdgesNum;
}SegContext;

enum depthMap
{
    BLOCK_64,
    BLOCK_32,
    BLOCK_16,
    BLOCK_8,
    BLOCK_4
};
//初始化分割数据内存
void initSegContext(SegContext *seg_context);
void destroySegContext(SegContext *seg_context);

#endif /* !TUCODECSEG_TYPE_H */
