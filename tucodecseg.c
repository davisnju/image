/************************************************************************/
/*                     tucodecseg                                       */
/************************************************************************/
/*!
* tucodecseg.h
* date: 2017/03/21 14:50
*
* author: dav1sNJU
*
*
*/
#include "tucodecseg.h"
#include "common.h"
#include "segment_image.h"
#include "improc.h"
#include "edgeproc.h"
#include "edge_band_cluster_v4_rewrite.h"
#include "all_element_init.h"
#include <assert.h>


#include "advance_tool_for_edge_band_cluster.h"
#include <string.h>

// 预先在contexxt打上边带线标签
//#define EDGELINE_LABEL

// 使用长边界生成边带
//#define LONG_EDGE_LINE

//显示边界线信息
//#define SHOW_LONGEDGE_INFO

// 时间复杂度测试
//#define TIMECOST_TEST

#ifdef TIMECOST_TEST
#include <time.h>
#endif /* TIMECOST_TEST */

//显示分割线segline
//#define SHOW_SED_SEG_LINE
// 显示长边界 
//#define SHOW_EDGE_LINE
//显示边界线
//#define SHOW_EDGE_LINE
//显示边带
//#define SHOW_BAND

/* typedef */
typedef struct edge_point Edge_point;
typedef struct Border_line Border_line;

float pixelDiff(Pixel *p1, Pixel *p2)
{
    return (float)square(p1->Y - p2->Y) + square(p1->Cb - p2->Cb) + square(p1->Cr - p2->Cr);
}

float predictUseRefAndDeltaAtLevel(Context * context, int block_r, int block_c, size_t lvl)
{
    BorderContext *pBC = &context->seg_context->border_context[lvl];
    SegContext *pSC = context->seg_context;
    Pixel *pPic = context->pic;
    Image_rgb *predictIm = (Image_rgb *)malloc(sizeof(Image_rgb));
    create_image_rgb(predictIm, 64, 64, 1);
    Image_rgb *Im = (Image_rgb *)malloc(sizeof(Image_rgb));
    create_image_rgb(Im, 64, 64, 1);
    Image_rgb *predictImShow = (Image_rgb *)malloc(sizeof(Image_rgb));
    create_image_rgb(predictImShow, 64, 64, 1);

    for (int r = 0; r < 64; ++r)
    {
        for (int c = 0; c < 64; ++c)
        {
            int abs_pos = (c + block_c) + (r + block_r) * context->pic_col;
            imRef(predictIm, c, r).r = pPic[abs_pos].Y;
            imRef(predictIm, c, r).g = pPic[abs_pos].Cb;
            imRef(predictIm, c, r).b = pPic[abs_pos].Cr;
            //imRef(predictImShow, c, r).r = pPic[abs_pos].Y;
            imRef(Im, c, r).r = pPic[abs_pos].Y;
            imRef(Im, c, r).g = pPic[abs_pos].Cb;
            imRef(Im, c, r).b = pPic[abs_pos].Cr;
        }
    }
    int search_width = _WIDTH_FOR_SEARCH_WIDTH_;  //3
    int ref_len = 2 * search_width + 1,           //7
        ref_extend_len = 4 * search_width + 1,    //13
        ref_extend_inter_len = _MULTIPLE_ * ref_extend_len + _MULTIPLE_ - 1,  //55
        ref_inter_len = _MULTIPLE_ * ref_len + _MULTIPLE_ - 1;
    Pixel ref[_ARRAY_SIZE_BOUNDARY_SEARCH_ + 5] = { 0 },
        ref_inter[4 * _WIDTH_FOR_SEARCH_WIDTH_ * _MULTIPLE_ + 2 * _MULTIPLE_ - 1 + 5] = { 0 };
    
    // 预测
    float mse = 0.0f;
    int bandPixelCnt = 0;
    if (lvl >= 0 && pSC->border_processed) // 32x32 16x16 8x8 4x4
    {
        int subBlockNum = 1 << (lvl << 1);
        int subBlockSize = 64 >> lvl;
        int subblock_c = 0, subblock_r = 0;
        printf("level %d, ", lvl);
        for (size_t bi = 0; bi < subBlockNum; bi++)
        {
            subblock_c = subBlockSize * (bi % (1 << lvl));   // 子块起始点相对坐标
            subblock_r = subBlockSize * (bi / (1 << lvl));  // 子块起始点相对坐标
            Border_seg *pBand = pBC->border[bi];
            int borderLineHeight = pBand->height;
            // 边带存在
            if (borderLineHeight > 0)
            {
                printf("block %d, band height %d:\n", bi, borderLineHeight);
                // 从预测信息里读取参考线
                for (int i = 0; i < ref_len; ++i)
                {
                    ref[i].Y = pBand->predict_info[i];
                }
                // 参考线插值
                interpolation(ref_inter, ref, ref_len, _MULTIPLE_ - 1);
                // 遍历一条边带的所有边带线
                int r0 = pBand->r[0], c0 = pBand->c[0];
                for (int borderlineIdx = 0; borderlineIdx < borderLineHeight; borderlineIdx++)
                {
                    int r = r0, c = c0;  // 边带起始点绝对坐标
                    int ref_delta = pBand->delta[borderlineIdx];
                    int8_t shape_delta = get_delta_in_encode(pBand->r[borderlineIdx], pBand->c[borderlineIdx], r0, c0, pBand->direction, ref_delta);
                    // 遍历一条边带线的像素
                    int j0 = -1;
                    for (int j = 0; j < 2 * pBand->width + 1; j++)
                    {
                        if (pBand->direction == 1) // 边带线水平扩展
                        {
                            c = c0 + shape_delta / 4 - pBand->width + j;
                            r = r0 + borderlineIdx;
                        }
                        else
                        {
                            c = c0 + borderlineIdx;
                            r = r0 + shape_delta / 4 - pBand->width + j;
                        }
                        int ref_intern_idx =(j + 1) * _MULTIPLE_ - 1;
                        ref_intern_idx += ref_delta;
                        ref_intern_idx = bound(ref_intern_idx, 0, ref_inter_len - 1);
                        if (c >= block_c + subblock_c && r >= block_r + subblock_r
                            && c < block_c + subblock_c + pBC->block_size&&r < block_r + subblock_r + pBC->block_size)
                        {
                            // shape ok, Y have some problem
                            imRef(predictIm, c, r).r = ref_inter[ref_intern_idx].Y;
                            imRef(predictIm, c, r).g = ref_inter[ref_intern_idx].Y;
                            imRef(predictIm, c, r).b = ref_inter[ref_intern_idx].Y;
                            imRef(predictImShow, c, r).r = 255;
                            mse += square(imRef(predictIm, c, r).r - imRef(Im, c, r).r);
                            if (j0<0)
                            {
                                j0 = j;
                                printf("bli:%2d, shape delta:%3d, ref delta:%3d |", borderlineIdx, shape_delta, ref_delta);
                            }
                            ++bandPixelCnt;
                            printf("%3d(%3d)", ref_inter[ref_intern_idx].Y, 
                                                imRef(Im,c,r).r - ref_inter[ref_intern_idx].Y);
                            if (abs(imRef(Im, c, r).r - ref_inter[ref_intern_idx].Y)>25)
                            {
                                printf("*");
                            }
                            else
                            {
                                printf(" ");
                            }
                        }
                        if (j0 >= 0 && j == 2 * pBand->width)
                        {
                            printf("\n");
                        }
                    }// for (int j = 0; j < 2 * pBand->width + 1; j++) 遍历一条边带线的像素
                } // for (int borderlineIdx = 0; borderlineIdx < borderLineHeight; borderlineIdx++) 遍历一条边带的所有边带线
            } // if (borderLineHeight > 0)
        } // for (size_t bi = 0; bi < subBlockNum; bi++)
    } // if (lvl >= 0 && pSC->border_processed)
    if (bandPixelCnt)
    {
        mse /= bandPixelCnt;
    }
    printf("mse at level %d = %f, band pixel count = %d\n", lvl, mse, bandPixelCnt);

    for (int r = 0; r < 64; ++r)
    {
        for (int c = 0; c < 64; ++c)if(!imRef(predictImShow, c, r).r)
        {
            int abs_pos = (c + block_c) + (r + block_r) * context->pic_col;
            imRef(predictIm, c, r).r = pPic[abs_pos].Y;
            imRef(predictIm, c, r).g = pPic[abs_pos].Cb;
            imRef(predictIm, c, r).b = pPic[abs_pos].Cr;
        }
    }

#if defined _MSC_VER && defined USE_OPENCV
    //显示
    char windowName[256];
    memset(windowName, 0, sizeof(windowName));
    sprintf(windowName, "image", lvl);
    if (!lvl)
    {
        showRGBImage(Im, windowName, 10 + 64 * 5, 10 + 400);
    }

    memset(windowName, 0, sizeof(windowName));
    sprintf(windowName, "predict at level %d", lvl);
    showRGBImage(predictIm, windowName, 10 + 64 * lvl, 10 + 400);
    memset(windowName, 0, sizeof(windowName));
    sprintf(windowName, "band at level %d", lvl);
    showRGBImage(predictImShow, windowName, 10 + 64 * lvl, 10 + 300);
    //cvWaitKey(0);
#endif

    return mse;
}

void predictUseRefAndDelta(Context *context, int block_r, int block_c)
{
    float mse_lvl_avg = 0.0f, mse_by_level[5] = { 0 };
    for (size_t lvl = 0; lvl < 5; lvl++)
    {
        float mse_at_lvl = predictUseRefAndDeltaAtLevel(context, block_r, block_c, lvl);
        mse_lvl_avg += mse_at_lvl;
        mse_by_level[lvl] = mse_at_lvl;
    }
    mse_lvl_avg /= 5;
    printf("avg mse = %f\n", mse_lvl_avg); 
    for (size_t lvl = 0; lvl < 5; lvl++)
    {
        printf("mse at lvl %d:%f\n", lvl, mse_by_level[lvl]);
    }

#if defined _MSC_VER && defined USE_OPENCV
    cvWaitKey(0);
#endif
}

//接口：生成64x64块边带信息, 返回边带数目0-1
int makeBandInfo(Context *context, int block_r, int block_c)
{
    int ret = 0;
    //int w = context->pic_col, h = context->pic_row;

    // 从文件读取中间结果
    FILE *fp;
    if (fp = fopen("edgeQuadtree.bin", "rb"))
    {
        uint8_t dataBuf[5][64][64] = { 0 };            // [lvli(0-4)][c][r]
        for (int i = 0; i < 5; i++)
        {
            uint8_t *pLine = (uint8_t*)(dataBuf[i]);
            fread(pLine, sizeof(uint8_t), 64 * 64, fp);
        }
        fclose(fp);

        for (int lvl = 0; lvl < 5; ++lvl)
        {
            Image *imb = create_new_image(64, 64, 1);
            for (int i = 0; i < 64; i++)
            {
                for (int j = 0; j < 64; j++)
                {
                    imRef(imb, j, i) = dataBuf[lvl][j][i] > 0;  //注意matlab中列优先的性质
                }
            }

            makeBandInfoAtLevel(context, block_r, block_c, imb, lvl);
        }



        predictUseRefAndDelta(context, block_r, block_c);

        ret = 1;
    }

#if defined _MSC_VER && defined USE_OPENCV
    cvWaitKey(0);
#endif
    return ret;
};

Border_line * getSubblockEdgeline(Context * context, int block_r, int block_c, int subblock_c, int subblock_r, int lvl, SegContext * pSC)
{
    Border_line * pBL = (Border_line *)malloc(sizeof(Border_line));

    int blockSize = 64 >> lvl;
    vector_int *edge_in_block = NULL;
    create_vector_int(&edge_in_block, 16);
    //vector_int *edge_id_ref = NULL;
    //create_vector_int(&edge_id_ref, pSC->longEdgesNum + 1);

    //for (int i = 0; i < pSC->longEdgesNum; ++i)
    //{
    //    push_back_vector_int(edge_id_ref, i + 1);
    //}

    for (int y = 0; y < blockSize; ++y)
    {
        for (int x = 0; x < blockSize; ++x)
        {
            if (imRef(pSC->edge_label_map, x + subblock_c, y + subblock_r) > 0)
            {
                push_back_vector_int(edge_in_block, imRef(pSC->edge_label_map, x + subblock_c, y + subblock_r));
            }
        }
    }
    if (edge_in_block->ele_num == 0)
    {
        pBL->pixel_count = -1;
        return pBL;
    }
    vector_int *edge_in_block_set = unique_vector_int(edge_in_block);
    if (edge_in_block_set->ele_num == 0)
    {
        pBL->pixel_count = -1;
        return pBL;
    }
    int maxProbEdgeId = -1;
    float maxEdgeProb = -0.1f;
    for (size_t i = 0; i < edge_in_block_set->ele_num; i++)
    {
        int edgeId = edge_in_block_set->head[i] - 1;    // edgeId = edgeLabel-1
        float edgeProb = pSC->edgesContrast[edgeId];
        if (edgeProb > maxEdgeProb)
        {
            maxEdgeProb = edgeProb;
            maxProbEdgeId = edgeId;
        }
    }
    EdgeLine *piEdgeLine = &pSC->longEdgeList[maxProbEdgeId];
    int in_block_edge_length = 0;
    for (size_t i = 0; i < edge_in_block->ele_num; i++)
    {
        if (edge_in_block->head[i] == (maxProbEdgeId + 1))
        {
            ++in_block_edge_length;
        }
    }

    pBL->edge_point = (Edge_point *)malloc(sizeof(Edge_point) * in_block_edge_length);
    pBL->pixel_count = 0;// in_block_edge_length;
    pBL->label = -1;
    Edgeline_elt *pElt = piEdgeLine->head;
    Edgeline_elt *pLastElt = piEdgeLine->end;
    int el_dir = getEdgeLineDirection(pElt, pLastElt);
    unsigned int edgePointIdx = 0;
    while (pElt)
    {
        Point *pep = &pElt->pos;
        pElt = pElt->next;

        if (edgePointIdx 
            && ((el_dir == 1 && pBL->edge_point[edgePointIdx - 1].x == pep->y) 
                || (el_dir != 1 && pBL->edge_point[edgePointIdx - 1].y == pep->x)))
        {
            break;
        }
        // not in this block
        if (((pep->x < subblock_c) || (pep->x >= subblock_c + blockSize))
            || ((pep->y < subblock_r) || (pep->y >= subblock_r + blockSize)))
        {
            continue;
        }
        // in this block
        Edge_point *pjEdgaePoint = &(pBL->edge_point[edgePointIdx]);
        pjEdgaePoint->x = pep->y;      // r a y(dw) h  使用相对64x64块的坐标
        pjEdgaePoint->y = pep->x;      // c b x(dw) w  使用相对64x64块的坐标
        pjEdgaePoint->label = -1;
        ++edgePointIdx;
    }
    pBL->pixel_count = edgePointIdx;

    destroy_vector_int(edge_in_block);
    //destroy_vector_int(edge_id_ref);
    destroy_vector_int(edge_in_block_set);
    return pBL;
}

// 采用绝对坐标比较边带线是否落在块内
int check_bandline_in_block(Border_seg * pPBand, int idx, int subblock_abs_r, int subblock_abs_c, int subBlockSize)
{
    // 边带线方向水平
    if (pPBand->direction == 1)
    {
        int bandline_col_scope[2] = { pPBand->c[idx] - pPBand->width, pPBand->c[idx] + pPBand->width };
        if ((
            !check_bound(bandline_col_scope[0], subblock_abs_c, subblock_abs_c + subBlockSize - 1)
            || !check_bound(bandline_col_scope[1], subblock_abs_c, subblock_abs_c + subBlockSize - 1)
             || !check_bound(subblock_abs_c,bandline_col_scope[0],bandline_col_scope[1])
             || !check_bound(subblock_abs_c + subBlockSize - 1, bandline_col_scope[0], bandline_col_scope[1])
            )
            && !check_bound(pPBand->r[idx], subblock_abs_r, subblock_abs_r + subBlockSize - 1))
        {
            return 1;
        }
    }
    else
    {
        int bandline_row_scope[2] = { pPBand->r[idx] - pPBand->width, pPBand->r[idx] + pPBand->width };
        if ((
            !check_bound(bandline_row_scope[0], subblock_abs_r, subblock_abs_r + subBlockSize - 1)
            || !check_bound(bandline_row_scope[1], subblock_abs_r, subblock_abs_r + subBlockSize - 1)
            || !check_bound(subblock_abs_r, bandline_row_scope[0], bandline_row_scope[1])
            || !check_bound(subblock_abs_r + subBlockSize - 1, bandline_row_scope[0], bandline_row_scope[1])
            )
            && !check_bound(pPBand->c[idx], subblock_abs_c, subblock_abs_c + subBlockSize - 1))
        {
            return 1;
        }
    }
    return 0;
}

// 重新计算delta、diff
void recompute_delta_diff(Context * context, BorderContext * pBC, Border_seg * pBand, Border_seg * pPBand, int start_borderline_idx
    , int block_abs_r, int block_abs_c, int subblock_abs_r, int subblock_abs_c, int subblock_size)
{
    int8_t delta = 0, prevDelta = 0;
    DiffVal_t line_diff_y[_ARRAY_SIZE_BOUNDARY_SEARCH_ + 5];
    DiffVal_t line_diff_u[_ARRAY_SIZE_BOUNDARY_SEARCH_ + 5];
    DiffVal_t line_diff_v[_ARRAY_SIZE_BOUNDARY_SEARCH_ + 5];
    int search_width = _WIDTH_FOR_SEARCH_WIDTH_;  //3
    int ref_len = 2 * search_width + 1,           //7
        ref_extend_len = 4 * search_width + 1,    //13
        ref_extend_inter_len = _MULTIPLE_ * ref_extend_len + _MULTIPLE_ - 1,  //55
        tar_inter_len = 0;
    Pixel ref[_ARRAY_SIZE_BOUNDARY_SEARCH_ + 5],
        tar[_ARRAY_SIZE_BOUNDARY_SEARCH_ * 2 + 5],
        tar_extend_inter[4 * _WIDTH_FOR_SEARCH_WIDTH_ * _MULTIPLE_ + 2 * _MULTIPLE_ - 1 + 5] = { 0 },
        ref_extend_inter[4 * _WIDTH_FOR_SEARCH_WIDTH_ * _MULTIPLE_ + 2 * _MULTIPLE_ - 1 + 5] = { 0 };
    // 参考线
    for (int i = 0; i < ref_len; ++i)
    {
        //int pos = pBand->direction == 1 ? (pBand->c[0] - search_width + i) + pBand->r[0] * context->pic_col
        //                           : (pBand->c[0]) + (pBand->r[0] - search_width + i) * context->pic_col;
        //ref[i] = context->pic[pos];
        ref[i].Y = pBand->predict_info[i];
    }
    // 计算参考线扩展-插值结果
    extend_and_interpolation(ref, ref_len, ref_extend_inter, ref_extend_len, _MULTIPLE_ - 1);
    for (int subbli = 0; subbli < pBand->height; ++subbli)
    {
        //pBand->delta[subbli] = pPBand->delta[start_borderline_idx + subbli] - pPBand->delta[start_borderline_idx];
        // 目标边带线
        for (int i = 0; i < ref_extend_len; ++i)
        {
            int pos = pBand->direction == 1 ? (pBand->c[subbli] - search_width * 2 + i) + pBand->r[subbli] * context->pic_col
                                       : (pBand->c[subbli]) + (pBand->r[subbli] - search_width * 2 + i) * context->pic_col;
            tar[i] = context->pic[pos];
        }
        // 计算目标线扩展-插值结果
        extend_and_interpolation(tar, ref_extend_len, tar_extend_inter, ref_extend_len, _MULTIPLE_ - 1);

        // 重新计算delta和diff
        compare_line_for_block(ref_extend_inter, ref_extend_inter_len, tar_extend_inter, ref_extend_inter_len, &delta,
            line_diff_y, line_diff_u, line_diff_v, _MULTIPLE_, prevDelta);

        prevDelta = delta;
        pBand->delta[subbli] = delta;
        int r, c;
        for (int i = 0; i < ref_len; ++i)
        {
            if (pBand->direction == 1)
            {
                c = pBand->c[subbli] - search_width + i;
                r = pBand->r[subbli];
            }
            else
            {
                c = pBand->c[subbli];
                r = pBand->r[subbli] - search_width + i;
            }
            if (!check_bound(c, subblock_abs_c, subblock_abs_c + subblock_size - 1)
                && !check_bound(r, subblock_abs_r, subblock_abs_r + subblock_size - 1))
            {
                int pos = (c - block_abs_c) + (r - block_abs_r) * 64;    // 绝对坐标(pic)转相对坐标(64x64)
                pBC->pic_border_diffY[pos] = line_diff_y[i];
                pBC->pic_border_diffCb[pos] = line_diff_u[i];
                pBC->pic_border_diffCr[pos] = line_diff_v[i];
            }
        }
    }
}

void copy_parent_band(Context *context, SegContext *pSC, int block_r, int block_c, int lvl
    , int subblock_rel_r, int subblock_rel_c, int subBlockSize)
{
    BorderContext *pBC = &pSC->border_context[lvl], *pPBC = &pSC->border_context[lvl - 1];
    int subBlockIdx = subblock_rel_c / subBlockSize + subblock_rel_r / subBlockSize * (1 << lvl);
    // 边带标签
    LabelVal_t inBlockLabel = -1;
    for (int r = 0; r < subBlockSize; r++)
    {
        for (int c = 0; c < subBlockSize; c++)
        {
            if (pPBC->borderLabelMap[(r + subblock_rel_r) * 64 + (c + subblock_rel_c)] > 0)
            {
                pBC->borderLabelMap[(r + subblock_rel_r) * 64 + (c + subblock_rel_c)] =
                    pBC->borderLabel0 + subBlockIdx;
                pBC->flatLabelMap[(r + subblock_rel_r) * 64 + (c + subblock_rel_c)] = -1;
                if (inBlockLabel < 0)
                {
                    inBlockLabel = pPBC->borderLabelMap[(r + subblock_rel_r) * 64 + (c + subblock_rel_c)];
                }
            }
            else
            {
                pBC->borderLabelMap[(r + subblock_rel_r) * 64 + (c + subblock_rel_c)] = 0;
            }
        }
    }
    int parentBorderIdx = inBlockLabel - pPBC->borderLabel0;
    assert(parentBorderIdx < pPBC->max_border_num);
    Border_seg *pPBand = pPBC->border[parentBorderIdx];     // 父节点边带
    assert(subBlockIdx < pBC->max_border_num);
    Border_seg *pBand = pBC->border[subBlockIdx];           // 待求边带

    // 其他边带信息:高度、宽度、边带线方向、中心点坐标（r,c）、delta、参考线像素值
    pBand->direction = pPBand->direction;                   // 边带线方向
    pBand->width = pPBand->width;                           // 扩展宽度，边带线范围[w-0+w]

    Point start_point_in_sub_block = { -1,-1 };
    int start_borderline_idx = -1;
    //printf("\nabs_pos(%d,%d),size=%d:\n", subblock_rel_r + block_r, subblock_rel_c + block_c, subBlockSize);
    for (int idx = 0; idx < pPBand->height; ++idx)
    {
        //printf("(%2d, %2d)", pPBand->r[idx], pPBand->c[idx]);  //绝对坐标 in Context Pic
        int isInBlock = check_bandline_in_block(pPBand, idx, subblock_rel_r + block_r, 
            subblock_rel_c + block_c,
            subBlockSize);// 采用绝对坐标比较边带线是否落在块内
        if (isInBlock)
        {
            if (start_point_in_sub_block.x < 0)
            {
                start_point_in_sub_block.x = pPBand->c[idx];
                start_point_in_sub_block.y = pPBand->r[idx];
                start_borderline_idx = idx;                  // 落在子块里的边带线在父块边带的index
                pBand->height = 0;
            }
            pBand->r[pBand->height] = pPBand->r[idx];     // 中心点坐标（r,c）
            pBand->c[pBand->height] = pPBand->c[idx];
            ++(pBand->height);                            // 高度
            //printf("*");
        }
        //if (idx < pPBand->height - 1)
        //{
        //    printf("->");
        //}
        //else
        //{
        //    printf("\n");
        //}

    }
    //printf("height = %d\n", pBand->height);
    // 更新ref
    //for (int  d = 0; d < pPBand->width * 2 + 1; d++)
    //{
    //    pBand->predict_info[d] = pPBand->predict_info[d];
    //}
    int pic_pos;
    for (int d = -pPBand->width; d < pPBand->width + 1; d++)
    {
        // 边带线方向水平
        if (pBand->direction == 1)
        {
            pic_pos = (start_point_in_sub_block.x + d) + (start_point_in_sub_block.y) * context->pic_col;
        }
        else // 边带线方向垂直
        {
            pic_pos = (start_point_in_sub_block.x) + (start_point_in_sub_block.y + d)* context->pic_col;
        }
        pBand->predict_info[pPBand->width + d] = context->pic[pic_pos].Y;
    }

    // 重新计算delta、diff
    recompute_delta_diff(context, pBC, pBand, pPBand, start_borderline_idx, block_r, block_c, subblock_rel_r + block_r, subblock_rel_c + block_c, subBlockSize);

    pBC->border_num += 1;
};

void copy_border_line(Border_line * dst, Border_line * src)
{
    if (src && dst)
    {
        memcpy(dst, src, sizeof(Border_line));
        memcpy(dst->edge_point, src->edge_point, sizeof(edge_point) * src->pixel_count);
    }
};

void makeBandInfoAtLevel(Context *context, int block_r, int block_c, Image *imb, int lvl)
{
    SegContext *pSC = context->seg_context;
    if (lvl == 0 && pSC->border_processed)
    {
        // 释放上一个64块的链表
        for (int i = 0; i < pSC->longEdgesNum; i++)
        {
            clearEdgeLine(&pSC->longEdgeList[i]);
        }
        for (int i = 0; i < pSC->edgesNum; i++)
        {
            destroyEdgeLine(&pSC->edgeList[i]);
        }

        pSC->edgesNum = 0;
        pSC->longEdgesNum = 0;
        pSC->border_processed = 0;
    }
    else if (lvl == 0 && !pSC->border_processed)
    {
        // initialize 64x64 block SegContext

        //从文件读取边界概率图
        FILE *fp;
        assert(fp = fopen("edge_prob_map.bin", "rb"));
        fread(pSC->edge_prob_map->data, sizeof(float), 64 * 64, fp);
        fclose(fp);

        Image *im = create_new_image(64, 64, 1);
        for (int y = 0; y < 64; y++)
        {
            for (int x = 0; x < 64; x++)
            {
                imRef(im, x, y) = imRef(pSC->edge_prob_map, x, y) > 0.3f;
            }
        }
        int minlength = 5;
        int linkStatus = edgelink(im, pSC->edgeList, &pSC->edgesNum, NULL/*edgeIm*/, NULL/*etype*/,
            pSC->longEdgeList, &pSC->longEdgesNum, minlength);

        // 计算每条边界的概率并保存
        int maxBorderSize = 0;
        for (int i = 0; i < pSC->longEdgesNum; i++)
        {
            EdgeLine *piEdgeLine = &pSC->longEdgeList[i];
            Edgeline_elt *pElt = piEdgeLine->head;
            if (piEdgeLine->ele_num >=2 && point_cmpfunc(&pElt->pos, &pElt->next->pos) > 0)
            {
                reverseEdgeLine(piEdgeLine);
            }
            unsigned int edgePointIdx = 0;
            pSC->edgesContrast[i] = 0;
            while (pElt)
            {
                Point *pep = &pElt->pos;
                // 记录到SegContext中
                imRef(pSC->edge_label_map, pep->x, pep->y) = i + 1;
                pSC->edgesContrast[i] += imRef(pSC->edge_prob_map, pep->x, pep->y);
                pElt = pElt->next;
                ++edgePointIdx;
                ++maxBorderSize;
            }
            pSC->edgesContrast[i] /= piEdgeLine->ele_num;
        }

        // 选概率最高的边界生成边带
        Border_line *borderLineList = getSubblockEdgeline(context, block_r, block_c, 0, 0, 0, pSC);
        pSC->border_context[0].borderLabel0 = 1;
        edge_band_cluster_for_block(context, block_r, block_c, lvl, 0, 0, pSC, borderLineList, 1);
        pSC->border_processed = 1;
    }

    rgb color[4096];
    for (int i = 0; i < 4096; ++i)
    {
        color[i] = random_rgb();
    }
    color[0].r = 0;
    color[0].g = 0;
    color[0].b = 0;
    Image_rgb *bandStartPointIm = (Image_rgb *)malloc(sizeof(Image_rgb));
    create_image_rgb(bandStartPointIm, 64, 64, 1);

    if (lvl >= 1 && pSC->border_processed) // 32x32 16x16 8x8 4x4
    {
        int subBlockNum = 1 << (lvl << 1);
        int subBlockSize = 64 >> lvl;
        int borderInSubBlocks = 0, subblock_c = 0, subblock_r = 0;
        Border_line * borderLineList = (Border_line *)malloc(subBlockNum * sizeof(Border_line));
        for (int i = 0; i < subBlockNum; ++i)
        {
            subblock_c = subBlockSize * (i % (1 << lvl));   // 子块起始点相对坐标
            subblock_r = subBlockSize * (i / (1 << lvl));  // 子块起始点相对坐标
            pSC->border_context[lvl].borderLabel0 = pSC->border_context[lvl - 1].borderLabel0 + square(1 << (lvl - 1));  //子块边带的起始标号

            // 检查是否存在父节点边带
            int existParentBand = 0;
            for (int r = 0; !existParentBand && r < subBlockSize; r++)
            {
                for (int c = 0; !existParentBand && c < subBlockSize; c++)
                {
                    existParentBand = pSC->border_context[lvl - 1].borderLabelMap[(r + subblock_r) * 64 + (c + subblock_c)] > 0;
                }
            }
            // 存在父节点边带, 则保持父节点信息
            if (existParentBand)
            {
                copy_parent_band(context, pSC, block_r, block_c, lvl, subblock_r, subblock_c, subBlockSize);
                continue;
            }
            // 不存在父节点边带
            // 获取子块中最高概率的边界线
            Border_line *pBL = getSubblockEdgeline(context, block_r, block_c, subblock_c, subblock_r, lvl, pSC);
            // 若上述边界线存在
            if (pBL->pixel_count > 0)
            {
                copy_border_line(&borderLineList[borderInSubBlocks], pBL);
                ++borderInSubBlocks;
                edge_band_cluster_for_block(context, block_r, block_c, lvl, subblock_r, subblock_c, pSC, pBL/*borderLineList*/, 1/*borderInSubBlocks*/);
            }

        }

        //
        for (int bli = 0; bli < borderInSubBlocks; ++bli)
        {
            for (int blpi = 0; blpi < borderLineList[bli].pixel_count; blpi++)
            {
                int r = borderLineList[bli].edge_point[blpi].x,
                    c = borderLineList[bli].edge_point[blpi].y;
                imRef(bandStartPointIm, c, r) = color[bli + 1];
            }
        }


    }
    //printf("\n%d:\n", lvl);
    //for (int y = 0; y < 64; y++)
    //{
        //printf("%2d: ", y);
        //for (int x = 0; x < 64; x++)
        //{
            //printf("%2d ", (pSC->border_context[lvl].borderLabelMap[x + y * 64]));
            //printf("%2d", imRef(pSC->edge_label_map, x, y));
        //}
        //printf("\n");
    //}

    //TODO:利用并查集对平坦区域进行标号
    // 

    //Image_rgb *colormap = (Image_rgb *)malloc(sizeof(Image_rgb));
    //create_image_rgb(colormap, 64*5*10, 10, 1);
    //for (int c = 0; c < 64 * 5; ++c)
    //{
    //    for (int i = 0; i < 10; ++i)
    //    {
    //        imRef(colormap, c, 0) = color[c];
    //        imRef(colormap, c, 0 + i) = color[c];
    //        imRef(colormap, c + i, 0) = color[c];
    //        imRef(colormap, c + i, 0 + i) = color[c];
    //    }
    //}

    Image_rgb *bandIm = (Image_rgb *)malloc(sizeof(Image_rgb));
    create_image_rgb(bandIm, 64, 64, 1);
    Image_rgb *edgeLabelmap = (Image_rgb *)malloc(sizeof(Image_rgb));
    create_image_rgb(edgeLabelmap, 64, 64, 1);
    //printf("\n%d:\n", lvl);
    for (int y = 0; y < 64; y++)
    {
        //printf("%2d: ", y);
        for (int x = 0; x < 64; x++)
        {
            //printf("%2d ", );
            int borderlabel = bound(pSC->border_context[lvl].borderLabelMap[x + y * 64], 0, 4095);
            imRef(bandIm, x, y) = color[borderlabel];
            int edgelabel = imRef(pSC->edge_label_map, x, y);
            imRef(edgeLabelmap, x, y) = color[edgelabel];
        }
        //printf("\n");
    }

#if defined _MSC_VER && defined USE_OPENCV

    //showRGBImage(edgeLabelmap, "elm", 300, 10);
    //cvWaitKey(0);
    char windowName[256];
    memset(windowName, 0, sizeof(windowName));
    sprintf(windowName, "%dBand", lvl);
    showRGBImage(bandIm, windowName, 10 + lvl * 64, 10 + 84 * 2);

    memset(windowName, 0, sizeof(windowName));
    sprintf(windowName, "%dBorderLine", lvl);
    showRGBImage(bandStartPointIm, windowName, 10 + lvl * 64, 10 + 84);

    memset(windowName, 0, sizeof(windowName));
    sprintf(windowName, "%dline at", lvl);
    showImage(imb, windowName, 10 + lvl * 64, 10);
    //showRGBImage(colormap, "colormap", 10, 30+64);
    if(lvl == 4)cvWaitKey(0);
#endif
}


//int getBandInfo(Context *context, int block_r, int block_c, int block_depth)
//{
//    BorderContext *pBC = &(context->seg_context->border_context[block_depth - 1]);
//    int ret = 0;
//    int w = context->pic_col, h = context->pic_row;
//
//
//
//
//    return ret;
//}

//////////////////////////////////////////////////////////////////////////

void gradientEdgeDetectionInterface(Context *context, float sigma, float gradThred, int offset)
{
    int w = context->pic_col, h = context->pic_row;
    Image *gradLine = (Image *)malloc(sizeof(Image));//边界线  二值图//malloc free ok
    create_image(gradLine, w, h, 1);
    gradientEdgeDetection(context, gradLine, sigma, gradThred, offset);
};

void gradientEdgeDetection(Context * cntxt, Image *gradLine, float sigma, float gradThred, int offset)
{
    int width = cntxt->pic_col, height = cntxt->pic_row;
    Image_f Y, Cr, Cb, *Py = &Y, *Pcr = &Cr, *Pcb = &Cb;
    create_imagef(Py, width, height, 1);
    create_imagef(Pcr, width, height, 1);
    create_imagef(Pcb, width, height, 1);

    // 预处理 smooth each color channel  
    int roffset = 0;
    for (int y = 0; y < height; y++)
    {
        for (int x = 0; x < width; x++)
        {
            int pos = x + roffset;
            imRef(Py, x, y) = cntxt->pic[pos].Y;
            imRef(Pcr, x, y) = cntxt->pic[pos].Cr;
            imRef(Pcb, x, y) = cntxt->pic[pos].Cb;
        }
        roffset += width;
    }

    smoothWrapper(Py, Py, sigma);
    smoothWrapper(Pcr, Pcr, sigma);
    smoothWrapper(Pcb, Pcb, sigma);

    gradientThreeChannelOffset(Py, Pcr, Pcb, gradLine, gradThred, offset);

    int edgesNum = 0;  //边界线数目
    uint32_t maxEdgeNo = width * height;
    EdgeLine *edgeList = (EdgeLine *)malloc(sizeof(EdgeLine) * maxEdgeNo); //边界线列表
    int longEdgeNo = 0; //长边界线数目
    EdgeLine *longEdgeList = (EdgeLine *)malloc(sizeof(EdgeLine) * maxEdgeNo); //长边界线列表
    int linkStatus = edgelink(gradLine, edgeList, &edgesNum, NULL/*edgeIm*/, NULL/*etype*/,
        longEdgeList, &longEdgeNo, 0);
    Image_f *mag = create_new_imagef(width, height, 0); // 图像对比度信息//malloc free ok

    gradientCntxt(cntxt, mag); //图像求梯度fx fy 归一化对比度

#if defined _MSC_VER  && defined SHOW_EDGE_LINE
    showImage(gradLine, "grad line", width, 10);//梯度边界检测
    Image_rgb *edgelIm = (Image_rgb *)malloc(sizeof(Image_rgb));
    create_image_rgb(edgelIm, width, height, 1);
#endif

#ifdef SHOW_LONGEDGE_INFO
    unsigned int minEdgeSize = 100000;
#endif // SHOW_LONGEDGE_INFO    

    int maxBorderSize = 0;
    Border_line *borderLineList = (Border_line *)malloc(sizeof(Border_line) * longEdgeNo); //边界线列表  //malloc free ok
    for (int i = 0; i < longEdgeNo; i++)
    {
#if defined _MSC_VER  && defined SHOW_EDGE_LINE
        rgb color = random_rgb();
#endif
        EdgeLine *piEdgeLine = &longEdgeList[i];
        Border_line *piBorderLine = &borderLineList[i];
        piBorderLine->edge_point = (Edge_point *)malloc(sizeof(Edge_point) * piEdgeLine->ele_num);//malloc free ok
        piBorderLine->pixel_count = piEdgeLine->ele_num;
        piBorderLine->label = -1;
        Edgeline_elt *pElt = piEdgeLine->head;
        unsigned int edgePointIdx = 0;
        while (edgePointIdx < piEdgeLine->ele_num)//for (unsigned int j = 0; j < longEdgeiSize; j++)
        {
            Point *pep = &pElt->pos;/* atEdgeVec(piEdgeVec, j);*/
            Edge_point *pjEdgaePoint = &(piBorderLine->edge_point[edgePointIdx]);
            pjEdgaePoint->x = pep->y;      // r a y(dw) h
            pjEdgaePoint->y = pep->x;      // c b x(dw) w
            pjEdgaePoint->label = -1;
            pElt = pElt->next;
            ++edgePointIdx;
            ++maxBorderSize;

#if defined _MSC_VER  && defined SHOW_EDGE_LINE
            imRef(edgelIm, pep->x, pep->y).r = color.r;
            imRef(edgelIm, pep->x, pep->y).g = color.g;
            imRef(edgelIm, pep->x, pep->y).b = color.b;
#endif /* _MSC_VER */
        }

#ifdef SHOW_LONGEDGE_INFO
        if (minEdgeSize > longEdgeList[i].ele_num) minEdgeSize = longEdgeList[i].ele_num;
#endif // SHOW_LONGEDGE_INFO   
    }

#ifdef SHOW_LONGEDGE_INFO
    printf("sigma %.2f, grad thred %.2f, offset %d\n", sigma, gradThred, offset);
    printf("edge %d, long edge %d, min length %d\n", edgesNum, longEdgeNo, minEdgeSize);
#endif // SHOW_LONGEDGE_INFO    

#if defined _MSC_VER  && defined SHOW_EDGE_LINE && defined USE_OPENCV
    showRGBImage(edgelIm, "Long Edge", 0, 10 + height);//显示长边界rgb图像
    //cvWaitKey(0);
#endif /* _MSC_VER */
    cntxt->border_num = 0;
    init_border(cntxt, maxBorderSize);
    for (int k = 0; k < maxBorderSize; k++)
        cntxt->border[k]->height = 0;

    edge_band_cluster(cntxt, borderLineList, longEdgeNo, mag/* 对比度 */);

#if defined _MSC_VER  && defined SHOW_BAND && defined USE_OPENCV
    // 显示边带
    Image *borderBinMap = (Image *)malloc(sizeof(Image));
    create_image(borderBinMap, cntxt->pic_col, cntxt->pic_row, 1);
    for (int y = 0; y < cntxt->pic_row; y++)
    {
        int roffset = y * cntxt->pic_col;
        for (int x = 0; x < cntxt->pic_col; x++)
        {
            if (cntxt->pic[roffset + x].border_label > 0)
                imRef(borderBinMap, x, y) = 1u;
}
    }
    showImage(borderBinMap, "Band BinIm", width, 10 + height);//显示长边界rgb图像
    //cvWaitKey(0);
#endif /* _MSC_VER */

    // 释放内存
    int borderLineNo = edgesNum;
    for (int i = 0; i < borderLineNo; i++)
        free(borderLineList[i].edge_point);
    free(borderLineList);

    for (int i = 0; i < longEdgeNo; i++)
        clearEdgeLine(&longEdgeList[i]);
    free(longEdgeList);

    for (int i = 0; i < edgesNum; i++)
        destroyEdgeLine(&edgeList[i]);
    free(edgeList);
    destroy_imagef(mag);
    free(mag);
};


/*
*
* minlength for edgelink
* lenthThred for edgeband
*/
void tucodecSeg(Context *cntxt, float sigma, float k, int min_size,
    const int minlength, const float lenthThred)
{
#ifdef TIMECOST_TEST
    time_t start = clock();
#endif // TIMECOST_TEST

    int w = cntxt->pic_col, h = cntxt->pic_row;
    // ============================= 图像分割 =========================================
    int num_ccs = 0;
    Image_f *blockLabelIm = create_new_imagef(w, h, 0);
    Universe *pu = segment_image_cntxt(cntxt, sigma, k, min_size, NULL, blockLabelIm, &num_ccs); // 分区并查集 //malloc free ok
    if (pu == NULL) assert(pu);
    destroy_universe(pu);
    //// =========================== 向context添加平坦区域标识信息 =========================== 
    for (int y = 0; y < h; y++)
    {
        int roffset = y * w;
        for (int x = 0; x < w; x++)
        {
            int pos = x + roffset;
            cntxt->pic[pos].border_label = -1;
            cntxt->pic[pos].block_label = (int)(imRef(blockLabelIm, x, y) + 0.5);
        }
        }

#ifdef TIMECOST_TEST
    time_t end = clock();
    printf("got %d componets, segment_image cost:%dms\n", num_ccs, (long)end - start);
#endif // TIMECOST_TEST
    // ============================= 边界提取 =========================================
    Image *segOutline = (Image *)malloc(sizeof(Image));//边界线  二值图//malloc free ok
    create_image(segOutline, w, h, 1);
    gradient(blockLabelIm, segOutline, 0.0f);//根据不同平坦区域标签的梯度变化确定边界
    destroy_imagef(blockLabelIm);
    free(blockLabelIm);
    int edgesNum = 0;  //边界线数目
    uint32_t maxEdgeNo = w * h;
    EdgeLine *edgeList = (EdgeLine *)malloc(sizeof(EdgeLine) * maxEdgeNo); //边界线列表//malloc free ok
    int longEdgeNo = 0; //长边界线数目
    EdgeLine *longEdgeList = (EdgeLine *)malloc(sizeof(EdgeLine) * maxEdgeNo); //长边界线列表//malloc free ok

#ifdef TIMECOST_TEST
    printf("edgelink\n");
    start = clock();
#endif // TIMECOST_TEST

    int linkStatus = edgelink(segOutline, edgeList, &edgesNum, NULL/*edgeIm*/, NULL/*etype*/,
        longEdgeList, &longEdgeNo, minlength);    // 根据边界图来对边界进行解析

#ifdef TIMECOST_TEST
    end = clock();
    printf("edge_link cost:%dms\n", end - start);
#endif // TIMECOST_TEST

    //assert(edgesNum == 1);
    //assert(longEdgeNo == 1);

    destroy_image(segOutline);
    free(segOutline);

    Image_f *mag = (Image_f *)malloc(sizeof(Image_f)); // 图像对比度信息//malloc free ok
    create_imagef(mag, w, h, 1);
    gradientCntxt(cntxt, mag); //图像求梯度fx fy 归一化对比度
#if defined _MSC_VER  && defined SHOW_EDGE_LINE
    Image_rgb *edgelIm = (Image_rgb *)malloc(sizeof(Image_rgb));
    create_image_rgb(edgelIm, w, h, 1);
#endif

#ifdef SHOW_LONGEDGE_INFO
    unsigned int minEdgeSize = 100000;
#endif // SHOW_LONGEDGE_INFO    

    // 长边界
    int maxBorderSize = 0;
    Border_line *borderLineList = (Border_line *)malloc(sizeof(Border_line) * longEdgeNo); //边界线列表  //malloc free ok
    for (int i = 0; i < longEdgeNo; i++)
    {
#if defined _MSC_VER  && defined SHOW_EDGE_LINE
        rgb color = random_rgb();
#endif
        EdgeLine *piEdgeLine = &longEdgeList[i];
        Border_line *piBorderLine = &borderLineList[i];
        piBorderLine->edge_point = (Edge_point *)malloc(sizeof(Edge_point) * piEdgeLine->ele_num);//malloc free ok
        piBorderLine->pixel_count = piEdgeLine->ele_num;
        piBorderLine->label = -1;
        Edgeline_elt *pElt = piEdgeLine->head;
        unsigned int edgePointIdx = 0;
        while (edgePointIdx < piEdgeLine->ele_num)//for (unsigned int j = 0; j < longEdgeiSize; j++)
        {
            Point *pep = &pElt->pos;/* atEdgeVec(piEdgeVec, j);*/
            Edge_point *pjEdgaePoint = &(piBorderLine->edge_point[edgePointIdx]);
            pjEdgaePoint->x = pep->y;      // r a y(dw) h
            pjEdgaePoint->y = pep->x;      // c b x(dw) w
            pjEdgaePoint->label = -1;
            pElt = pElt->next;
            ++edgePointIdx;
            ++maxBorderSize;

#if defined _MSC_VER  && defined SHOW_EDGE_LINE
            imRef(edgelIm, pep->x, pep->y).r = color.r;
            imRef(edgelIm, pep->x, pep->y).g = color.g;
            imRef(edgelIm, pep->x, pep->y).b = color.b;
#endif /* _MSC_VER */
        }

#ifdef SHOW_LONGEDGE_INFO
        if (minEdgeSize > longEdgeList[i].ele_num)minEdgeSize = longEdgeList[i].ele_num;
#endif // SHOW_LONGEDGE_INFO    

    }

#ifdef SHOW_LONGEDGE_INFO
    printf("edge %d, long edge %d, min length %d\n", edgesNum, longEdgeNo, minEdgeSize);
#endif // SHOW_LONGEDGE_INFO    

#if defined _MSC_VER  && defined SHOW_EDGE_LINE && defined USE_OPENCV
    showRGBImage(edgelIm, "Long Edge", 10, 100);//显示长边界rgb图像
    //cvWaitKey(0);
#endif /* _MSC_VER */

#ifdef TIMECOST_TEST
    printf("edge band clustering...\n");
    start = clock();
#endif // TIMECOST_TEST

    // ============================== 生成边带 ==============================
    cntxt->border_num = 0;
    init_border(cntxt, maxBorderSize);
    for (int k = 0; k < maxBorderSize; k++)
        cntxt->border[k]->height = 0;

    edge_band_cluster(cntxt, borderLineList, longEdgeNo, mag);

#ifdef TIMECOST_TEST
    end = clock();
    printf("\nedge band cluster cost:%dms\n", end - start);
#endif // TIMECOST_TEST

#if defined _MSC_VER  && defined SHOW_BAND && defined USE_OPENCV
    // 显示边带
    Image *borderBinMap = (Image *)malloc(sizeof(Image));
    create_image(borderBinMap, cntxt->pic_col, cntxt->pic_row, 1);
    for (int y = 0; y < cntxt->pic_row; y++)
    {
        int roffset = y * cntxt->pic_col;
        for (int x = 0; x < cntxt->pic_col; x++)
        {
            if (cntxt->pic[roffset + x].border_label > 0)
                imRef(borderBinMap, x, y) = 1u;
    }
}
    showImage(borderBinMap, "Band BinIm", 10 + w, 100);//显示长边界rgb图像
    //cvWaitKey(0);
#endif /* _MSC_VER */

    // ============================== 释放内存 ==============================
    for (int i = 0; i < longEdgeNo; i++)
        free(borderLineList[i].edge_point);
    free(borderLineList);
    for (int i = 0; i < longEdgeNo; i++)
        clearEdgeLine(&longEdgeList[i]);
    free(longEdgeList);
    for (int i = 0; i < edgesNum; i++)
        destroyEdgeLine(&edgeList[i]);
    free(edgeList);
    destroy_imagef(mag);
    free(mag);
};


void tucodecSegInterface(Context *cntxt)
{

    // 提取边界线， 并生成边带，添加到cntxt里
    tucodecSeg(cntxt, SIGMA, SEG_PARAM_C, MINSIZE, MINLENGTH, LENGTHTHRED);

};


// ============================ 暂时未使用 =================================

void tucodecSegInterfaceUseSED(Context *context, char *SED_Segline)
{
    if (!SED_Segline)
    {
        assert("can't find segline file!");
        return;
    }
    tucodecSegUseSED(context, SED_Segline);
};


void getSeglineFromFile(char * SED_Segline_Dir, Image * segOutline)
{
    int pixelNum = segOutline->w * segOutline->h;
    float *fbuf = (float *)malloc(sizeof(float) * pixelNum);
    FILE * inputfp = fopen(SED_Segline_Dir, "rb");
    fread((char *)fbuf, pixelNum, sizeof(float), inputfp);
    fclose(inputfp);

    for (int i = 0; i < pixelNum; i++)
    {
        segOutline->data[i] = fbuf[i] > 0.1f ? 1u : 0u;
    }

};

void tucodecSegUseSED(Context *cntxt, char *SED_Segline)
{

    int w = cntxt->pic_col, h = cntxt->pic_row;
    int num_ccs = 0;
    Image_f *blockLabelIm = (Image_f *)malloc(sizeof(Image_f));//平坦区域标签图//malloc free ok
    create_imagef(blockLabelIm, w, h, 1);

    // 图像分割
    Universe *pu = segment_image_cntxt(cntxt, SIGMA, SEG_PARAM_C, MINSIZE, NULL, blockLabelIm, &num_ccs); // 分区并查集 //malloc free ok
    if (pu == NULL)
    {
        //printf("segment_image 划分失败！");
        return;
    }
    destroy_universe(pu);

    //////////////////////////////////////////////////////////////////////////
    //读取SED边界
    Image *segOutlineSED = (Image *)malloc(sizeof(Image));//边界线  二值图//malloc free ok
    create_image(segOutlineSED, cntxt->real_pic_col, cntxt->real_pic_row, 1);

    getSeglineFromFile(SED_Segline, segOutlineSED);
#if defined SHOW_SED_SEG_LINE && defined USE_OPENCV
    showSEDSegline(segOutlineSED, "SED segline", 10, 100);
    //cvWaitKey(0);
#endif
    //////////////////////////////////////////////////////////////////////////
    //根据不同平坦区域标签的梯度变化确定边界
    Image *segOutline = (Image *)malloc(sizeof(Image));//边界线  二值图//malloc free ok
    create_image(segOutline, w, h, 1);
    gradient(blockLabelIm, segOutline, 0.0f);
    //////////////////////////////////////////////////////////////////////////

    int edgesNum = 0;  //边界线数目
    uint32_t maxEdgeNo = w * h;
    EdgeLine *edgeList = (EdgeLine *)malloc(sizeof(EdgeLine) * maxEdgeNo); //边界线列表//malloc free ok
    int longEdgeNo = 0; //长边界线数目
    EdgeLine *longEdgeList = (EdgeLine *)malloc(sizeof(EdgeLine) * maxEdgeNo); //长边界线列表//malloc free ok
                                                                               //Image *edgeIm = (Image *)malloc(sizeof(Image));  // 边界二值图
                                                                               //create_image(edgeIm, w, h, 1);
        //
    int linkStatus = edgelink(segOutlineSED/*segOutline*/, edgeList, &edgesNum, NULL/*edgeIm*/, NULL/*etype*/,
        longEdgeList, &longEdgeNo,
#ifndef LONG_EDGE_LINE
        0
#else //#ifdef LONG_EDGE_LINE
        MINLENGTH
#endif /* LONG_EDGE_LINE */                                  
        );    // 根据边界图来对边界进行解析


    destroy_image(segOutline);
    free(segOutline);
    destroy_image(segOutlineSED);
    free(segOutlineSED);


    // ============================== 生成边带 ==============================
    Image_f *mag = (Image_f *)malloc(sizeof(Image_f)); // 图像对比度信息//malloc free ok
    create_imagef(mag, w, h, 1);

    gradientCntxt(cntxt, mag); //图像求梯度fx fy 归一化对比度
    // 长边界 time-consuming
#if defined _MSC_VER  && defined SHOW_EDGE_LINE
    Image_rgb *edgelIm = (Image_rgb *)malloc(sizeof(Image_rgb));
    create_image_rgb(edgelIm, w, h, 1);
#endif
    int maxBorderSize = 0;
    Border_line *borderLineList = (Border_line *)malloc(sizeof(Border_line) * longEdgeNo); //边界线列表  //malloc free ok
    for (int i = 0; i < longEdgeNo; i++)
    {
#if defined _MSC_VER  && defined SHOW_EDGE_LINE
        rgb color = random_rgb();
#endif /* _MSC_VER */
        EdgeLine *piEdgeLine = &longEdgeList[i];
        Border_line *piBorderLine = &borderLineList[i];
        piBorderLine->edge_point = (Edge_point *)malloc(sizeof(Edge_point) * piEdgeLine->ele_num);//malloc free ok
        piBorderLine->pixel_count = piEdgeLine->ele_num;
        piBorderLine->label = -1;
        Edgeline_elt *pElt = piEdgeLine->head;
        unsigned int edgePointIdx = 0;
        while (edgePointIdx < piEdgeLine->ele_num)//for (unsigned int j = 0; j < longEdgeiSize; j++)
        {
            Point *pep = &pElt->pos;/* atEdgeVec(piEdgeVec, j);*/
            Edge_point *pjEdgaePoint = &(piBorderLine->edge_point[edgePointIdx]);
            pjEdgaePoint->x = pep->y;      // r a y(dw) h
            pjEdgaePoint->y = pep->x;      // c b x(dw) w
            pjEdgaePoint->label = -1;
            pElt = pElt->next;
            ++edgePointIdx;
            ++maxBorderSize;

#if defined _MSC_VER  && defined SHOW_EDGE_LINE
            imRef(edgelIm, pep->x, pep->y).r = color.r;
            imRef(edgelIm, pep->x, pep->y).g = color.g;
            imRef(edgelIm, pep->x, pep->y).b = color.b;
#endif /* _MSC_VER */
        }
        }

#if defined _MSC_VER  && defined SHOW_EDGE_LINE && defined USE_OPENCV
    showRGBImage(edgelIm, "Long Edge", 100, 100);//显示长边界rgb图像
                                                //cvWaitKey(0);
#endif /* _MSC_VER */
    //// =========================== 向context添加边界/平坦区域标识信息 =========================== 
    for (int y = 0; y < h; y++)
    {
        int roffset = y * w;
        for (int x = 0; x < w; x++)
        {
            int pos = x + roffset;
#ifdef EDGELINE_LABEL
            //imRef(edgeIm,x,y)>0

            int edgeLabel = imRef(edgeIm, x, y); //像素的长边带标签
            if (!edgeLabel)    // 像素在平坦区域
            {
                cntxt->pic[pos].border_label = -1;
                cntxt->pic[pos].block_label = (int)round(imRef(blockLabelIm, x, y) + 0.5);
                continue;
            }
            else
            {
                //imRef(edgeb_map_bd, x, y) = 1u;
                cntxt->pic[pos].border_label = edgeLabel;
                cntxt->pic[pos].block_label = -1;
            }
#else
            cntxt->pic[pos].border_label = -1;
            cntxt->pic[pos].block_label = (int)(imRef(blockLabelIm, x, y) + 0.5);
#endif /* EDGELINE_LABEL */
    }
}
    //printf("%d\n", maxLabel);
    destroy_imagef(blockLabelIm);
    free(blockLabelIm);

    // 边带生成
    cntxt->border_num = 0;
    init_border(cntxt, maxBorderSize);
    for (int k = 0; k < maxBorderSize; k++)
    {
        cntxt->border[k]->height = 0;
    }

    edge_band_cluster(cntxt,
        borderLineList,
#ifndef LONG_EDGE_LINE
        edgesNum
#else //#ifdef LONG_EDGE_LINE
        longEdgeNo
#endif /* LONG_EDGE_LINE */       
        , mag // 对比度
        );

    int borderLineNo = 0;
#ifndef LONG_EDGE_LINE
    borderLineNo = edgesNum;
#else //#ifdef LONG_EDGE_LINE
    borderLineNo = longEdgeNo;
#endif /* !LONG_EDGE_LINE */

    for (int i = 0; i < borderLineNo; i++)
    {
        free(borderLineList[i].edge_point);
    }
    free(borderLineList);

    for (int i = 0; i < longEdgeNo; i++)
    {
        clearEdgeLine(&longEdgeList[i]);
    }
    free(longEdgeList);

    for (int i = 0; i < edgesNum; i++)
    {
        destroyEdgeLine(&edgeList[i]);
    }
    free(edgeList);
    destroy_imagef(mag);
    free(mag);
};

void ContextToYUV(Context *context, Image_yuv *im)
{
    //assert(!"The method or operation is not implemented.");
    int W = context->pic_col, H = context->pic_row;
    for (int y = 0; y < H; y++)
    {
        int roffset = y * W;
        for (int x = 0; x < W; x++)
        {
            imRef(im, x, y).y = context->pic[x + roffset].Y;
            imRef(im, x, y).u = context->pic[x + roffset].Cr;
            imRef(im, x, y).v = context->pic[x + roffset].Cb;
        }
    }
};


#if defined _MSC_VER && defined USE_OPENCV
static void showRGBImage(Image_rgb *im, char *name, int x, int y)
{
    int w = im->w, h = im->h;
    // 显示和保存图像
    int channels = 3;
    IplImage* ipl = cvCreateImage(cvSize(w, h), IPL_DEPTH_8U, channels);//cvCloneImage(pImg);//
    uchar * iplData = (uchar *)ipl->imageData;
    int step = ipl->widthStep;
    // 转换边带图像
    for (int x = 0; x < w; ++x)
    {
        for (int y = 0; y < h; ++y)
        {
            iplData[y * step + x * channels + 0] = imRef(im, x, y).b;     //B
            iplData[y * step + x * channels + 1] = imRef(im, x, y).g;     //G
            iplData[y * step + x * channels + 2] = imRef(im, x, y).r;     //R
        }
    }
    /*if (!cvSaveImage(name, ipl, 0))
    printf("Could not save: %s\n", name);*/
    cvNamedWindow(name, 1);//创建窗口
    cvMoveWindow(name, x, y);// 移动窗口
    cvShowImage(name, ipl);//显示图像
};

static void showImage(Image *im, char *name, int x, int y)
{
    int w = im->w, h = im->h;
    // 显示和保存图像
    int channels = 3;
    IplImage* ipl = cvCreateImage(cvSize(w, h), IPL_DEPTH_8U, channels);//cvCloneImage(pImg);//
    uchar * iplData = (uchar *)ipl->imageData;
    int step = ipl->widthStep;
    // 转换边带图像
    for (int x = 0; x < w; ++x)
    {
        for (int y = 0; y < h; ++y)
        {
            iplData[y * step + x * channels + 0] = imRef(im, x, y) > 0 ? 255 : 0;     //B
            iplData[y * step + x * channels + 1] = imRef(im, x, y) > 0 ? 255 : 0;     //G
            iplData[y * step + x * channels + 2] = imRef(im, x, y) > 0 ? 255 : 0;     //R
        }
    }
    /*if (!cvSaveImage(name, ipl, 0))
    printf("Could not save: %s\n", name);*/
    cvNamedWindow(name, 1);//创建窗口
    cvMoveWindow(name, x, y);// 移动窗口
    cvShowImage(name, ipl);//显示图像
};
static void showSEDSegline(Image *im, char *name, int x, int y)
{
    int w = im->w, h = im->h;
    // 显示和保存图像
    int channels = 3;
    IplImage* ipl = cvCreateImage(cvSize(w, h), IPL_DEPTH_8U, channels);//cvCloneImage(pImg);//
    uchar * iplData = (uchar *)ipl->imageData;
    int step = ipl->widthStep;
    // 转换边带图像
    for (int x = 0; x < w; ++x)
    {
        for (int y = 0; y < h; ++y)
        {
            iplData[y * step + x * channels + 0] = imRef(im, x, y)>0 ? 255 : 0 /*(int)(imRef(im, x, y) * 255 + 0.5f)*/;     //B
            iplData[y * step + x * channels + 1] = imRef(im, x, y)>0 ? 255 : 0 /*(int)(imRef(im, x, y) * 255 + 0.5f)*/;     //G
            iplData[y * step + x * channels + 2] = imRef(im, x, y) > 0 ? 255 : 0 /*(int)(imRef(im, x, y) * 255 + 0.5f)*/;     //R
        }
    }
    /*if (!cvSaveImage(name, ipl, 0))
    printf("Could not save: %s\n", name);*/
    cvNamedWindow(name, 1);//创建窗口
    cvMoveWindow(name, x, y);// 移动窗口
    cvShowImage(name, ipl);//显示图像
};
#endif // _MSC_VER  &&  USE_OPENCV
