/*!
 * vector_c.h
 * date: 2017/03/21 14:50
 *
 * author: dav1sNJU
 *
 * brief: vector ctype
 *
 *
*/
#pragma once
#ifndef VECTOR_C_H
#define VECTOR_C_H

#include "image.h"


typedef struct vector_int
{
    unsigned int mem_len;/*申请的内存长度，单位字节*/
    unsigned int user_len;/*已经使用的内存长度，单位字节*/
    unsigned int ele_num;/*当前元素个数*/
    unsigned int type_len;/*存放数据类型的长度，单位字节*/
    int *head;/*数据头指针*/
    int *end;/*数据尾指针*/
}vector_int;


/*
*函数功能:创建vector
*函数参数:
*       vec: vec结构指针
*       len: 初次开辟的元素个数
*/
int create_vector_int(vector_int **vec, const unsigned int len);

/*
*函数功能:测试容器是不是已经占满内存
*函数参数:
*       vec: vec结构指针
*/
int isFull_vector_int(vector_int *vec);
/*
*函数功能:重新开辟更大的空间来存放元素
*函数参数:
*       vec: vec结构指针
*/
int resize_vector_int(vector_int *vec);
/*
*函数功能:在容器尾部添加元素
*函数参数:
*       vec: vec结构指针
*       value:所要添加元素，务必要保证元素类型一致
*/
int push_back_vector_int(vector_int *vec, int value);
/*
*函数功能:删除容器尾部的元素
*函数参数:
*       vec: vec结构指针
*/
int pop_back_vector_int(vector_int *vec);
/*
*函数功能:定位位置为i的元素值
*函数参数:
*       vec: vec结构指针
*       i: 在容器中的位置，注意位置是从0开始的
*/
int at_vector_int(vector_int *vec, unsigned int i);
/*
*函数功能:在位置i处插入元素
*函数参数:
*       vec: vec结构指针
*       i:在容器中的位置，注意位置是从0开始的
*       value: 插入的元素
*/
//int insert_vector_int(vector_int *vec, unsigned int i, int value);
/*
*函数功能:在位置i处删除元素
*函数参数:
*       vec: vec结构指针
*       i:删除元素在容器中的位置
*/
//int erase_vector_int(vector_int *vec, unsigned int i);
/*
*函数功能:在位置i处修改元素
*函数参数:
*       vec: vec结构指针
*       i:修改元素在容器中的位置
*/
//int set_vector_int(vector_int *vec, unsigned int i, int value);
/*
*函数功能:返回容器当前第一个元素的值
*函数参数:
*       vec: vec结构指针
* 返回值:该类型的元素
*/
int front_vector_int(vector_int *vec);
/*
*函数功能:返回容器尾部的元素
*函数参数:
*       vec: vec结构指针
* 返回值:该类型的元素
*/
int back_vector_int(vector_int *vec);
/*
*函数功能:返回容器内申请的最大元素个数
*函数参数:
*       vec: vec结构指针
*/
unsigned int getCapacity_vector_int(vector_int* vec);
/*
*函数功能:返回当前容器内元素个数
*函数参数:
*       vec: vec结构指针
*/
unsigned int getSize_vector_int(vector_int* vec);
/*
*函数功能:清除vector中所有元素
*函数参数:
*       vec: vec结构指针
*/
int clear_vector_int(vector_int* vec);
/*
*函数功能:返回迭代容器首位置
*函数参数:
*       vec: vec结构指针
*/
int begin_vector_int(vector_int* vec);
/*
*函数功能:返回迭代器容器中尾位置
*函数参数:
*       vec: vec结构指针
*/
int end_vector_int(vector_int* vec);
/*
*函数功能:返回迭代器容器是否为空
*函数参数:
*       vec: vec结构指针
*/
int empty_vector_int(vector_int* vec);

/* 销毁向量 */
int destroy_vector_int(vector_int* vec);

/* int 比较函数 */
int int_cmpfunc(const void * a, const void * b);

/* 在vector_int中查找x 
* 不存在，则返回 -1
* 否则返回位置  0-elenum-1
*/
int find_in_vector_int(vector_int* vec, const int x);

//类似于unique(v1.begin,v1.end)
vector_int * unique_vector_int(vector_int * v1);

//////////////////////////////////////////////////////////////////////////
// vector<Point>
typedef struct vector_point
{
    unsigned int mem_len;/*申请的内存长度，单位字节*/
    unsigned int user_len;/*已经使用的内存长度，单位字节*/
    unsigned int ele_num;/*当前元素个数*/
    unsigned int type_len;/*存放数据类型的长度，单位字节*/
    Point *head;/*数据头指针*/
    Point *end;/*数据尾指针*/
}vector_point;
int create_vector_point(vector_point **vec, const unsigned int len);
int isFull_vector_point(vector_point *vec);
int resize_vector_point(vector_point *vec);
int push_back_vector_point(vector_point *vec, Point *value);
int pop_back_vector_point(vector_point *vec);
Point at_vector_point(vector_point *vec, unsigned int i);
Point front_vector_point(vector_point *vec);
Point back_vector_point(vector_point *vec);
unsigned int getCapacity_vector_point(vector_point* vec);
unsigned int getSize_vector_point(vector_point* vec);
int clear_vector_point(vector_point* vec);
int begin_vector_point(vector_point* vec);
int end_vector_point(vector_point* vec);
int empty_vector_point(vector_point* vec);
int destroy_vector_point(vector_point* vec);

/* 在vector_point中查找point x
* 不存在，则返回 -1
* 否则返回位置  0-elenum-1
*/
int find_in_vector_point(vector_point* vec, const Point *x);

#endif /* !VECTOR_C_H */