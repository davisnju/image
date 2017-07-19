/*!
* disjoint-set.h
* date: 2017/03/21 14:45
*
* author: dav1sNJU
*
* brief: 并查集算法，用于segment-graph
*
*/
#include "disjoint_set.h"
#include <stdlib.h>

int get_universe_size(Universe *u, int x)
{
    return u->elts[x].size;
};

int get_universe_num_sets(Universe *u)
{
    return u->num;
};

/*
* 初始化结构体 Universe
* 令Universe->num=元素数目
* 令Universe->elts中的所有元素
* rank=0, size=1, p = i;
* @param pu 指向目标结构体
* @param elements 表示有元素数目
*/
Universe * init_universe(const int elements)
{
    Universe *pu = (Universe *)malloc(sizeof(Universe));
    Uni_elt *elts = (Uni_elt *)malloc(sizeof(Uni_elt) * elements);
    pu->num = elements;
    pu->size = elements;
    pu->elts = elts;
    for (int i = 0; i < elements; i++)
    {
        elts[i].rank = 0;
        elts[i].size = 1;
        elts[i].p = i;
    }
    return pu;
};

/*
* 销毁Universe结构体，释放内存
*/
int destroy_universe(Universe *pu)
{
    if (!pu) return 1;

    free(pu->elts);
    pu->elts = NULL;
    free(pu);            //
    return 0;
};

/*
* 查找节点x的根节点id
*/
int find_in_universe(Universe *pu, int x)
{
    int y = x;
    Uni_elt *pEltX = &pu->elts[x];
    while (y != pu->elts[y].p)  //如果和 parent 的 node id 不一样，即节点y不是根节点
        y = pu->elts[y].p;      //赋值为父亲的ID号
    pEltX->p = y;          //更新根节点，加快下次查询速度
    return y;
};

/*
* 合并节点x, y, 并查集连通元素数目-1
*/
void join_universe(Universe *pu, int x, int y)
{
    Uni_elt *pEltX = &pu->elts[x], *pEltY = &pu->elts[y];
    if (pEltX->rank > pEltY->rank)
    {
        pEltY->p = x;
        pEltX->size += pEltY->size;
    }
    else
    {
        pEltX->p = y;
        pEltY->size += pEltX->size;
        if (pEltX->rank == pEltY->rank)
            pEltY->rank++;
    }
    pu->num--;
};
