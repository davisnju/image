/*!
* edge.h
* date: 2017/03/21 14:46
*
* author: dav1sNJU
*
* brief: 边界元、边界线、边界线链表元素的定义及操作
*
*
*/
#include "edge.h"
#include <assert.h>

int edge_compare(const void *a, const void *b)
{
    const Edge *ea = a, *eb = b;
    if (ea->w < eb->w) return -1;
    else if (ea->w > eb->w) return 1;
    else return 0;
};

void initEdgeLine(EdgeLine * line)
{
    line->head = NULL;
    line->end = NULL;
    line->ele_num = 0;
};

EdgeLine *creatEdgeLine() 
{
    EdgeLine *ret = (EdgeLine *)malloc(sizeof(EdgeLine));
    ret->head = NULL;
    ret->end = NULL;
    ret->ele_num = 0;
    return ret;
};

int destroyEdgeLine(EdgeLine * line)
{
    if (!line)
        return 0;
    Edgeline_elt *pElt = line->head;
    Edgeline_elt *pNextElt = NULL;
    while (pElt != line->end)
    {
        pNextElt = pElt->next;
        free(pElt);
        pElt = pNextElt;
    }
    free(line->end);//释放end
    line->head = NULL;
    line->end = NULL;
    line->ele_num = 0;
    return 0;
};

int clearEdgeLine(EdgeLine * line)
{
    line->head = NULL;
    line->end = NULL;
    line->ele_num = 0;
    return 0;
};

int appendPointToEdgeLine(EdgeLine * line, Point *point)
{
    Edgeline_elt *pElt = (Edgeline_elt *)malloc(sizeof(Edgeline_elt));
    if (!pElt)
        return -1;
    pElt->pos = *point;
    pElt->prev = NULL;
    pElt->next = NULL;
    return appendEltToEdgeLine(line, pElt);
};

int appendEltToEdgeLine(EdgeLine * line, Edgeline_elt * pElt)
{
    if (!line || !pElt)
        return -1;
    //empty
    if (line->ele_num == 0)
    {
        line->head = pElt;
        line->end = pElt;
        line->ele_num = 1;
    }
    else
    {
        pElt->prev = line->end;
        line->end->next = pElt;
        line->end = pElt;
        line->ele_num += 1;
    }
    return 0;
};

/* 将链表元素添加到边界线指向的链表尾 */
void appendEdgeLineToEdgeLine(EdgeLine * lineDst, EdgeLine * lineSrc) 
{
    assert(lineDst);
    assert(lineSrc); 
    if (lineSrc->ele_num == 0)
    {
        return;
    }
    else if (lineDst->ele_num == 0)
    {
        lineDst->head = lineSrc->head;
        lineDst->end = lineSrc->end;
        lineDst->ele_num = lineSrc->ele_num;
    }
    else
    {
        lineDst->end->next = lineSrc->head;
        lineSrc->head->prev = lineDst->end;
        lineDst->end = lineSrc->end;
        lineDst->ele_num += lineSrc->ele_num;
    }
    clearEdgeLine(lineSrc);
};

int reverseEdgeLine(EdgeLine * line)
{
    if (!line)
        return 0;

    Edgeline_elt *pElt = line->head;
    if (!pElt)
        return 0;

    Edgeline_elt *pTmp = NULL;

    //交换每个元素的prev-next指针
    while (pElt)
    {
        pTmp = pElt->prev;
        pElt->prev = pElt->next;
        pElt->next = pTmp;
        pElt = pElt->prev;    // 反转之后的下一个
    }
    //交换边界线head-end指针
    Edgeline_elt *pEnd = line->end;
    line->end = line->head;
    line->head = pEnd;

    return 0;
};

int getEdgeLineDirection(Edgeline_elt * pElt, Edgeline_elt * pLastElt)
{
    if (pElt->pos.x == pLastElt->pos.x)
    {
        return 1;
    }
    double tan = abs((pElt->pos.y - pLastElt->pos.y) / (pElt->pos.x - pLastElt->pos.x));
    if (tan > 1)
    {
        return 1;
    }
    else
    {
        return 2;
    }
    return 0;
}
// 以下代码暂时没用到
/*
int createEdgeVec(EdgeVec **vec, const unsigned int len)
{
    int vec_creat_ret = -1;
    if (len >= 0)
    {
        *vec = (EdgeVec *)malloc(sizeof(EdgeVec));
        (*vec)->type_len = sizeof(Point);
        (*vec)->ele_num = 0;
        (*vec)->user_len = 0;
        (*vec)->mem_len = (len * sizeof(Point));
        Point *vec_creat_ptr1 = (Point *)malloc((*vec)->mem_len);
        if (vec_creat_ptr1 == NULL)
            return -1;
        (*vec)->head = vec_creat_ptr1;
        (*vec)->end = vec_creat_ptr1;
        vec_creat_ret = 0;
    }
    return vec_creat_ret;
};

int isFullEdgeVec(EdgeVec *vec)
{
    int ret = -1;
    if ((vec->mem_len - vec->user_len) <= (2 * vec->type_len))
        ret = 0;
    return ret;
};
int resizeEdgeVec(EdgeVec *vec)
{
    int ret = -1;
    void* ptr = (void*)malloc(2 * vec->mem_len);
    if (ptr != NULL)
    {
        memcpy(ptr, (void*)(vec->head), vec->user_len);
        free(vec->head);
        vec->head = ptr;
        vec->end = (Point *)(vec->head) + vec->user_len / vec->type_len;
        vec->mem_len = 2 * vec->mem_len;
        ret = 0;
    }
    return ret;
};
int push_back_EdgeVec(EdgeVec *vec, Point *value)
{
    int vec_psbk_ret = -1;
    if (vec != NULL)
    {
        if (isFullEdgeVec(vec) == 0)
        {
            resizeEdgeVec(vec);
        }
        //memcpy(vec->end, (void *)value, vec->type_len);
        vec->end[0] = *value;
        vec->end++;
        vec->user_len += vec->type_len;
        vec->ele_num++;
        vec_psbk_ret = 0;
    }
    return vec_psbk_ret;
};
int pop_back_EdgeVec(EdgeVec *vec)
{
    int vec_ppbk_ret = -1;
    if (vec->ele_num != 0)
    {
        vec->end -= 1;
        vec->user_len -= vec->type_len;
        vec->ele_num--;
        vec_ppbk_ret = 0;
    }
    return  vec_ppbk_ret;
};
Point atEdgeVec(EdgeVec *vec, unsigned int i)
{
    Point vec_at_ret;
    vec_at_ret.x = -1;
    vec_at_ret.y = -1;
    if (!vec || !vec->head)
    {
        return vec_at_ret;
    }
    if (i < vec->ele_num && i >= 0)
    {
        Point *vec_at_ptr = vec->head + i;
        vec_at_ret = *vec_at_ptr;
    }
    if (vec_at_ret.x < 0)
    {
        return vec_at_ret;
    }
    return vec_at_ret;
};
int insertToEdgeVec(EdgeVec *vec, unsigned int i, Point *value)
{
    int vec_inst_ret = -1;
    if (isFullEdgeVec(vec) == 0)
    {
        resizeEdgeVec(vec);
    }
    if (i <= vec->ele_num && i >= 0)
    {
        Point* vec_inst_ptr1 = vec->end + 1;
        Point* vec_inst_ptr2 = vec->end;
        int vec_inst_num = vec->ele_num - i;  //
        Point* vec_inst_ptr3 = vec->head + i;
        while (vec_inst_num--)
        {
            memcpy(vec_inst_ptr1, vec_inst_ptr2, vec->type_len);
            vec_inst_ptr1 -= 1;
            vec_inst_ptr2 -= 1;
        }
        memcpy(vec_inst_ptr3, (void*)value, vec->type_len);
        vec->user_len += vec->type_len;
        vec->end += 1;
        vec->ele_num++;
        vec_inst_ret = 0;
    }
    return vec_inst_ret;
};
int eraseInEdgeVec(EdgeVec *vec, unsigned int i)
{
    int vec_erase_ret = -1;
    if (vec->ele_num != 0)
    {
        if (i >= 0 && i < vec->ele_num)
        {
            Point* vec_erase_ptr1 = vec->head + i;
            Point* vec_erase_ptr2 = vec_erase_ptr1 + 1;
            int vec_erase_num = vec->ele_num - i - 1;  //
            while (vec_erase_num--)
            {
                memcpy(vec_erase_ptr1, vec_erase_ptr2, vec->type_len);
                vec_erase_ptr1 += 1;
                vec_erase_ptr2 += 1;
            }
            vec->ele_num--;
            vec->user_len -= vec->type_len;
            vec->end -= 1;
            vec_erase_ret = 0;
        }
    }
    return vec_erase_ret;
};
Point frontInEdgeVec(EdgeVec *vec)
{
    Point vec_front_ret;
    if (vec->ele_num != 0)
    {
        vec_front_ret = atEdgeVec(vec, 0);
    }
    return vec_front_ret;
};
Point backInEdgeVec(EdgeVec*vec)
{
    Point vec_back_ret;
    if (vec->ele_num != 0)
    {
        vec_back_ret = atEdgeVec(vec, (vec->ele_num - 1));
    }
    return vec_back_ret;
};
unsigned int getEdgeVecCapacity(EdgeVec *vec)
{
    int ret = 0;
    if (vec != NULL)
        ret = vec->mem_len;
    return ret;
};
unsigned int getEdgeVecSize(EdgeVec* vec)
{
    int ret = 0;
    if (vec != NULL)
    {
        if (vec->type_len == 0)  // error
            return ret;
        ret = vec->user_len / vec->type_len;
    }
    return ret;
};
int clearEdgeVec(EdgeVec* vec)
{
    int ret = -1;
    if (vec != NULL)
    {
        free(vec->head);
        vec->head = NULL;//(Edge *)malloc(vec->mem_len);   //
        vec->end = vec->head;
        vec->ele_num = 0;
        vec->mem_len = 0;
        vec->user_len = 0;
    }
    return ret;
};   
unsigned int beginOfEdgeVec(EdgeVec* vec)
{
    return 0;
};
unsigned int endOfEdgeVec(EdgeVec* vec)
{
    return (vec->ele_num - 1);
};
*/