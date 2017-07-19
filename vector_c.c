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
#include "vector_c.h"
#include <stdlib.h>
#include <string.h>
#include <stdint.h>
#include <assert.h>

int create_vector_int(vector_int **vec, const unsigned int len)
{
    int vec_creat_ret = -1;

    *vec = (vector_int *)malloc(sizeof(vector_int));
    (*vec)->type_len = sizeof(int);
    (*vec)->ele_num = 0;
    (*vec)->user_len = 0;
    (*vec)->mem_len = (len * sizeof(int));
    void *vec_creat_ptr1 = (void *)malloc((*vec)->mem_len);
    if (vec_creat_ptr1 == NULL)
        return -1;
    (*vec)->head = vec_creat_ptr1;
    (*vec)->end = vec_creat_ptr1;
    vec_creat_ret = 0;

    return vec_creat_ret;
}

/*
*函数功能:测试容器是不是已经占满内存
*函数参数:
*       vec: vec结构指针
*/
int isFull_vector_int(vector_int *vec)
{
    int ret = 0;
    if (vec && (vec->mem_len - vec->user_len) <= (vec->type_len))
    {
        ret = 1;
    }
    return ret;
}
/*
*函数功能:重新开辟更大的空间来存放元素
*函数参数:
*       vec: vec结构指针
*/
int resize_vector_int(vector_int *vec)
{
    int ret = -1;
    if (!vec)
    {
        return ret;
    }
    void* ptr = (void*)malloc(2 * vec->mem_len);
    if (ptr != NULL)
    {
        memcpy(ptr, vec->head, vec->user_len);
        free(vec->head);
        vec->head = ptr;
        vec->end = (int *)(vec->head) + vec->user_len / vec->type_len;
        vec->mem_len = 2 * vec->mem_len;
        ret = 0;
    }
    return ret;
}
/*
*函数功能:在容器尾部添加元素
*函数参数:
*       vec: vec结构指针
*       value:所要添加元素，务必要保证元素类型一致
*/
int push_back_vector_int(vector_int *vec, int value)
{
    int vec_psbk_ret = -1;
    if (vec != NULL)
    {
        if (isFull_vector_int(vec))
        {
            resize_vector_int(vec);
        }
        memcpy(vec->end, (void *)&value, vec->type_len);
        vec->end += 1;
        vec->user_len += vec->type_len;
        vec->ele_num++;
        vec_psbk_ret = 0;
    }
    return vec_psbk_ret;
}
/*
*函数功能:删除容器尾部的元素
*函数参数:
*       vec: vec结构指针
*/
int pop_back_vector_int(vector_int *vec)
{
    int vec_ppbk_ret = -1;
    if (vec && vec->ele_num != 0)
    {
        vec->end -= 1;
        vec->user_len -= vec->type_len;
        vec->ele_num--;
        vec_ppbk_ret = 0;
    }
    return  vec_ppbk_ret;
}
/*
*函数功能:定位位置为i的元素值
*函数参数:
*       vec: vec结构指针
*       i: 在容器中的位置，注意位置是从0开始的
*/
int at_vector_int(vector_int *vec, unsigned int i)
{
    int vec_at_ret = -1;
    if (i < vec->ele_num)
    {
        int *vec_at_ptr = vec->head + i;
        vec_at_ret = *vec_at_ptr;
    }
    if (vec_at_ret < 0)
    {
        return -1;
    }
    return vec_at_ret;
}
/*
*函数功能:在位置i处插入元素
*函数参数:
*       vec: vec结构指针
*       i:在容器中的位置，注意位置是从0开始的
*       value: 插入的元素
*/
//int insert_vector_int(vector_int *vec, unsigned int i, int value)
//{
//    int vec_inst_ret = -1;
//    if (isFull_vector_int(vec) == 0)
//    {
//        resize_vector_int(vec);
//    }
//    if (i <= vec->ele_num)
//    {
//        int* vec_inst_ptr1 = vec->end + 1;
//        int* vec_inst_ptr2 = vec->end;
//        int vec_inst_num = vec->ele_num - i;  //
//        int* vec_inst_ptr3 = vec->head + i;
//        while (vec_inst_num--)
//        {
//            memcpy(vec_inst_ptr1, vec_inst_ptr2, vec->type_len);
//            vec_inst_ptr1 -= 1;
//            vec_inst_ptr2 -= 1;
//        }
//        memcpy(vec_inst_ptr3, (void*)&value, vec->type_len);
//        vec->user_len += vec->type_len;
//        vec->end += 1;
//        vec->ele_num++;
//        vec_inst_ret = 0;
//    }
//    return vec_inst_ret;
//}
/*
*函数功能:在位置i处删除元素
*函数参数:
*       vec: vec结构指针
*       i:删除元素在容器中的位置
*/
//int erase_vector_int(vector_int *vec, unsigned int i)
//{
//    int vec_erase_ret = -1;
//    if (vec->ele_num != 0)
//    {
//        if (i < vec->ele_num)
//        {
//            int* vec_erase_ptr1 = vec->head + i;
//            int* vec_erase_ptr2 = vec_erase_ptr1 + 1;
//            int vec_erase_num = vec->ele_num - i - 1;  //
//            while (vec_erase_num--)
//            {
//                memcpy(vec_erase_ptr1, vec_erase_ptr2, vec->type_len);
//                vec_erase_ptr1 += 1;
//                vec_erase_ptr2 += 1;
//            }
//            vec->ele_num--;
//            vec->user_len -= vec->type_len;
//            vec->end -= 1;
//            vec_erase_ret = 0;
//        }
//    }
//    return vec_erase_ret;
//}
/*
*函数功能:在位置i处修改元素
*函数参数:
*       vec: vec结构指针
*       i:修改元素在容器中的位置
*/
//int set_vector_int(vector_int *vec, unsigned int i, int value)
//{
//    int vec_at_ret = -1;
//    if (i < vec->ele_num)
//    {
//        int *vec_at_ptr = vec->head + i;
//        *vec_at_ptr = value;
//        vec_at_ret = 0;
//    }
//    if (i >= vec->ele_num)
//    {
//        insert_vector_int(vec, i, value);
//    }
//    return vec_at_ret;
//}
/*
*函数功能:返回容器当前第一个元素的值
*函数参数:
*       vec: vec结构指针
* 返回值:该类型的元素
*/
int front_vector_int(vector_int *vec)
{
    int vec_front_ret = -1;
    if (vec && vec->ele_num != 0)
    {
        vec_front_ret = at_vector_int(vec, 0);
    }
    return vec_front_ret;
}
/*
*函数功能:返回容器尾部的元素
*函数参数:
*       vec: vec结构指针
* 返回值:该类型的元素
*/
int back_vector_int(vector_int *vec)
{
    int vec_back_ret = -1;
    if (vec && vec->ele_num != 0)
    {
        vec_back_ret = at_vector_int(vec, (vec->ele_num - 1));
    }
    return vec_back_ret;
}
/*
*函数功能:返回容器内申请的最大元素个数
*函数参数:
*       vec: vec结构指针
*/
unsigned int getCapacity_vector_int(vector_int* vec)
{
    int ret = 0;
    if (vec != NULL)
        ret = vec->mem_len / vec->type_len;
    return ret;
}
/*
*函数功能:返回当前容器内元素个数
*函数参数:
*       vec: vec结构指针
*/
unsigned int getSize_vector_int(vector_int* vec)
{
    unsigned int ret = 0u;
    if (vec != NULL)
        ret = vec->user_len / vec->type_len;
    return ret;
}
/*
*函数功能:清除vector中所有元素
*函数参数:
*       vec: vec结构指针
*/
int clear_vector_int(vector_int* vec)
{
    int ret = -1;
    if (vec != NULL)
    {
        vec->ele_num = 0;
        //free(vec->head);
        //vec->head = (int *)malloc(vec->mem_len);   //
        memset((void *)vec->head, 0, vec->mem_len);
        vec->end = vec->head;
        vec->user_len = 0;
        ret = 0;
    }
    return ret;
}
/*
*函数功能:返回迭代容器首位置
*函数参数:
*       vec: vec结构指针
*/
int begin_vector_int(vector_int* vec)
{
    if (!vec)
    {
        return -1;
    }
    return 0;
}
/*
*函数功能:返回迭代器容器中尾位置
*函数参数:
*       vec: vec结构指针
*/
int end_vector_int(vector_int* vec)
{
    if (!vec)
    {
        return -1;
    }
    return (vec->ele_num - 1);
}
/*
*函数功能:返回迭代器容器是否为空
*函数参数:
*       vec: vec结构指针
*/
int empty_vector_int(vector_int* vec)
{
    if (!vec)
    {
        return -1;
    }
    return vec->ele_num == 0;
}

/* 销毁向量 free ok */
int destroy_vector_int(vector_int* vec)
{
    int ret = -1;
    if (vec != NULL)
    {
        vec->end = NULL;
        free(vec->head);
        //vec->ele_num = 0;
        //vec->head = NULL;   //
        //vec->user_len = 0;
        free(vec);
        ret = 0;
    }
    return ret;
}

int int_cmpfunc(const void * a, const void * b)
{
    return (*(int*)a - *(int*)b);
}

/* 二分查找 */
int bsearch_t(int array[], int low, int high, int target)
{
    if (low > high) return -1;

    int mid = (high - low) / 2 + low;
    if (array[mid] > target)
        return    bsearch_t(array, low, mid - 1, target);
    if (array[mid] < target)
        return    bsearch_t(array, mid + 1, high, target);

    //if (midValue == target)
    return mid;
}

/* 在vector_int中查找x
* 不存在，则返回 -1
* 否则返回位置  0-elenum-1
*/
int find_in_vector_int(vector_int* vec, const int x)
{
    int ret = -1;

    ret = bsearch_t(vec->head, 0, vec->ele_num - 1, x);

    return ret;
};

//类似于unique(v1.begin,v1.end)
vector_int * unique_vector_int(vector_int * v1)
{
    unsigned int v1Len = v1->ele_num;
    vector_int *ret = NULL;
    create_vector_int(&ret, v1Len + 5);

    qsort(v1->head, v1Len, sizeof(int), int_cmpfunc);

    for (unsigned int i = 0; i < v1Len; i++)
    {
        int v1e = v1->head[i];
        //if (find_in_vector_int(v2, v1e) > -1 && find_in_vector_int(ret, v1e) == -1)
        if (!ret->ele_num || v1e != back_vector_int(ret)/*ret->head[ret->ele_num - 1]*/)
        {
            push_back_vector_int(ret, v1e);
        }
    }
    return ret;
}

//////////////////////////////////////////////////////////////////////////
int create_vector_point(vector_point **vec, const unsigned int len)
{
    int vec_creat_ret = -1;
    *vec = (vector_point *)malloc(sizeof(vector_point));
    (*vec)->type_len = sizeof(Point);
    (*vec)->ele_num = 0;
    (*vec)->user_len = 0;
    (*vec)->mem_len = (len * sizeof(Point));
    void *vec_creat_ptr1 = (void *)malloc((*vec)->mem_len);
    if (vec_creat_ptr1 == NULL)
        return -1;
    (*vec)->head = vec_creat_ptr1;
    (*vec)->end = vec_creat_ptr1;
    vec_creat_ret = 0;

    return vec_creat_ret;
};
int isFull_vector_point(vector_point *vec)
{
    int ret = 0;
    if (vec && (vec->mem_len - vec->user_len) <= (vec->type_len))
    {
        ret = 1;
    }
    return ret;
};
int resize_vector_point(vector_point *vec)
{
    int ret = -1;
    if (!vec)
    {
        return ret;
    }
    void *ptr = (void *)malloc(2 * vec->mem_len);
    if (ptr)
    {
        memcpy(ptr, vec->head, vec->user_len);
        free(vec->head);
        vec->head = ptr;
        vec->end = (Point *)(vec->head) + vec->user_len / vec->type_len;
        vec->mem_len = 2 * vec->mem_len;
        ret = 0;
    }
    return ret;
};
int push_back_vector_point(vector_point *vec, Point *value)
{
    int vec_psbk_ret = -1;
    if (vec != NULL)
    {
        if (isFull_vector_point(vec))
        {
            resize_vector_point(vec);
        }
        memcpy(vec->end, (void *)value, vec->type_len);
        vec->end += 1;
        vec->user_len += vec->type_len;
        vec->ele_num++;
        vec_psbk_ret = 0;
    }
    return vec_psbk_ret;
};
int pop_back_vector_point(vector_point *vec)
{
    int vec_ppbk_ret = -1;
    if (vec && vec->ele_num != 0)
    {
        vec->end -= 1;
        vec->user_len -= vec->type_len;
        vec->ele_num--;
        vec_ppbk_ret = 0;
    }
    return  vec_ppbk_ret;
};

Point at_vector_point(vector_point *vec, unsigned int i)
{
    Point vec_at_ret = { -1, -1 };
    if (vec && i < vec->ele_num)
    {
        Point *vec_at_ptr = vec->head + i;
        vec_at_ret = *vec_at_ptr;
    }
    return vec_at_ret;
};
Point front_vector_point(vector_point *vec)
{
    Point vec_front_ret;
    vec_front_ret.x = -1;
    vec_front_ret.y = -1;
    if (vec && vec->ele_num > 0)
    {
        vec_front_ret = at_vector_point(vec, 0);
    }
    return vec_front_ret;
};
Point back_vector_point(vector_point *vec)
{
    Point vec_back_ret = { -1, -1 };
    if (vec && vec->ele_num > 0)
    {
        vec_back_ret = at_vector_point(vec, (vec->ele_num - 1));
    }
    return vec_back_ret;
};
unsigned int getCapacity_vector_point(vector_point* vec)
{
    unsigned int ret = 0;
    if (vec != NULL)
    {
        ret = vec->mem_len / vec->type_len;
    }
    return ret;
};
unsigned int getSize_vector_point(vector_point* vec)
{
    unsigned int ret = 0;
    if (vec != NULL)
    {
        ret = vec->user_len / vec->type_len;
    }
    return ret;
};
int clear_vector_point(vector_point* vec)
{
    int ret = -1;
    if (vec != NULL)
    {
        free(vec->head);
        vec->ele_num = 0;
        vec->head = (Point *)malloc(vec->mem_len);   //
        vec->end = vec->head;
        vec->user_len = 0;
        ret = 0;
    }
    return ret;
};
int begin_vector_point(vector_point* vec)
{
    if (!vec)
    {
        return -1;
    }
    return 0;
};
int end_vector_point(vector_point* vec)
{
    if (!vec)
    {
        return -1;
    }
    return (vec->ele_num - 1);
};
int empty_vector_point(vector_point* vec)
{
    if (!vec)
    {
        return -1;
    }
    return vec->ele_num == 0;
};
int destroy_vector_point(vector_point* vec)
{
    int ret = -1;
    if (vec != NULL)
    {
        free(vec->head);
        vec->ele_num = 0;
        vec->head = NULL;   //
        vec->end = vec->head;
        vec->user_len = 0;
        free(vec);
        ret = 0;
    }
    return ret;
};

/* 二分查找 */
static int bsearch_point_vec(Point array[], int low, int high, const Point *target)
{
    if (low > high) return -1;

    int mid = (high - low) / 2 + low;
    if (point_cmpfunc(&array[mid],target) > 0)
        return    bsearch_point_vec(array, low, mid - 1, target);
    if (point_cmpfunc(&array[mid], target) < 0)
        return    bsearch_point_vec(array, mid + 1, high, target);

    //if (midValue == target)
    return mid;
};

/* 在vector_point中查找point x
* 不存在，则返回 -1
* 否则返回位置  0-elenum-1
*/
int find_in_vector_point(vector_point* vec, const Point *x)
{
    int ret = -1;

    ret = bsearch_point_vec(vec->head, 0, vec->ele_num - 1, x);

    return ret;
};
