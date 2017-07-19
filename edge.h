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
#pragma once
#ifndef EDGELIST_H
#define EDGELIST_H

#include "vector_c.h"
#include <stdlib.h>
#include <string.h>

/*
 * 边界元
 * segment-graph 用到的结构体
 */
typedef struct edge {
    float w;
    int a, b;
} Edge;
/* 边界元比较函数，比较w值 */
int edge_compare(const void *a, const void *b);

//边界线链表元素
typedef struct edgeline_elt
{
    Point pos;
    struct edgeline_elt * prev;
    struct edgeline_elt * next;
}Edgeline_elt;

//边界线结构体
typedef struct edgeline
{
    unsigned int ele_num;/*当前元素个数*/
    Edgeline_elt *head;  /*数据头指针*/
    Edgeline_elt *end;   /*数据尾指针*/
}EdgeLine;

/* 新建边界线，并初始化成员 */
EdgeLine *creatEdgeLine();

/* 初始化边界线 初始化其成员 */
void initEdgeLine(EdgeLine * line);

/* 清除边界线所指向的链表中所有元素内存空间 */
int destroyEdgeLine(EdgeLine * line);

/* 清除边界线 不清除元素
 * 因为可能这些元素可能在其它边界线
 * 程序确保最后能正确释放这些元素
 */
int clearEdgeLine(EdgeLine * line);

/* 将Point添加到边界线指向的链表尾 */
int appendPointToEdgeLine(EdgeLine * line, Point *point);

/* 将链表元素添加到边界线指向的链表尾 */
int appendEltToEdgeLine(EdgeLine * line, Edgeline_elt * elelt);

/* 将链表元素添加到边界线指向的链表尾 */
void appendEdgeLineToEdgeLine(EdgeLine * lineDst, EdgeLine * lineSrc);

/* 反转双向链表 */
int reverseEdgeLine(EdgeLine * line);

int getEdgeLineDirection(Edgeline_elt * pElt, Edgeline_elt * pLastElt);

// ============================ 以下代码暂时没用到 ======================================
/*
typedef vector_point EdgeVec;

int createEdgeVec(EdgeVec **vec, const unsigned int len);

int isFullEdgeVec(EdgeVec *vec);
int resizeEdgeVec(EdgeVec *vec);
int push_back_EdgeVec(EdgeVec *vec, Point *value);
int pop_back_EdgeVec(EdgeVec *vec);
Point atEdgeVec(EdgeVec *vec, unsigned int i);
int insertToEdgeVec(EdgeVec *vec, unsigned int i, Point *value);
int eraseInEdgeVec(EdgeVec *vec, unsigned int i);
Point frontInEdgeVec(EdgeVec *vec);
Point backInEdgeVec(EdgeVec *vec);
unsigned int getEdgeVecCapacity(EdgeVec* vec);
unsigned int getEdgeVecSize(EdgeVec* vec);
int clearEdgeVec(EdgeVec* vec);
unsigned int beginOfEdgeVec(EdgeVec* vec);
unsigned int endOfEdgeVec(EdgeVec* vec);
*/


#endif /* !EDGELIST_H */