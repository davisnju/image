/*!
 * tucodecseg_type.c
 * date: 2017/03/28 10:56
 *
 * author: dav1sNJU
 *
 * brief: 
 *
 * TODO: long description
 *
*/

#include "tucodecseg_type.h"

void initSegContext(SegContext *seg_context)
{
    seg_context->border_processed = 0;
    for (uint8_t lvl_i = 0; lvl_i < 5; ++lvl_i)
    {
        int block_size = 64 >> lvl_i;
        seg_context->border_context[lvl_i].block_size = block_size;
        seg_context->border_context[lvl_i].border_num = 0;
        seg_context->border_context[lvl_i].borderLabelMap = (LabelVal_t *)calloc(4096, sizeof(LabelVal_t));
        seg_context->border_context[lvl_i].flatLabelMap = (LabelVal_t *)calloc(4096, sizeof(LabelVal_t));

        seg_context->border_context[lvl_i].pic_border_diffY = (DiffVal_t *)calloc(4096, sizeof(DiffVal_t));
        seg_context->border_context[lvl_i].pic_border_diffCr = (DiffVal_t *)calloc(4096, sizeof(DiffVal_t));
        seg_context->border_context[lvl_i].pic_border_diffCb = (DiffVal_t *)calloc(4096, sizeof(DiffVal_t));

        int max_border_num = (1 << lvl_i) * (1 << lvl_i);
        seg_context->border_context[lvl_i].border = (Border_seg **)malloc(sizeof(Border_seg *) * max_border_num);
        seg_context->border_context[lvl_i].max_border_num = max_border_num;
        seg_context->border_context[lvl_i].border_num=0;
        Border_seg **pBorder = seg_context->border_context[lvl_i].border;
        for (int i = 0; i < max_border_num; i++)
        {
            pBorder[i] = (Border_seg *)malloc(sizeof(Border_seg));
            pBorder[i]->direction = 0;
            pBorder[i]->height = 0;
            pBorder[i]->width = 0;
            pBorder[i]->delta = (int8_t *)calloc(block_size, sizeof(int8_t));
            pBorder[i]->predict_info = (uint8_t *)calloc(7, sizeof(uint8_t));
            pBorder[i]->r = (int *)calloc(block_size, sizeof(int));
            pBorder[i]->c = (int *)calloc(block_size, sizeof(int));
        }
    }

    seg_context->edge_prob_map = create_new_imagef(64, 64, 1); // 图像对比度信息
    seg_context->edge_label_map = create_new_image_int(64, 64, 1);// 64x64块边界标签图
    seg_context->edgesNum = 0;
    seg_context->longEdgesNum = 0;
    seg_context->edgeList = (EdgeLine *)malloc(sizeof(EdgeLine) * 4096); //边界线列表
    seg_context->longEdgeList = (EdgeLine *)malloc(sizeof(EdgeLine) * 4096); //边界线列表
    seg_context->edgesContrast = (float *)calloc(4096, sizeof(float));

}

void destroySegContext(SegContext *seg_context)
{
    seg_context->border_processed = 0;
    for (uint8_t lvl_i = 0; lvl_i < 5; ++lvl_i)
    {

        int max_border_num = 4096;// 1 << (lvl_i << 1);

        Border_seg **pBorder = seg_context->border_context[lvl_i].border;
        for (int i = 0; i < max_border_num; i++)
        {
            free(pBorder[i]->c);
            free(pBorder[i]->r);
            free(pBorder[i]->predict_info);
            free(pBorder[i]->delta);
            free(pBorder[i]);
        }

        free(seg_context->border_context[lvl_i].border);
        free(seg_context->border_context[lvl_i].pic_border_diffY);
        free(seg_context->border_context[lvl_i].pic_border_diffCr);
        free(seg_context->border_context[lvl_i].pic_border_diffCb);

        free(seg_context->border_context[lvl_i].borderLabelMap);
        free(seg_context->border_context[lvl_i].flatLabelMap);
    }
    free(seg_context->edgesContrast);
    free(seg_context->longEdgeList);
    free(seg_context->edgeList);
    destroy_imagef(seg_context->edge_prob_map);
    free(seg_context->edge_prob_map);
    
    free(seg_context);
}