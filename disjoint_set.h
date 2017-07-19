/*!
 * disjoint-set.h
 * date: 2017/03/21 14:45
 *
 * author: dav1sNJU
 *
 * brief: 并查集算法，用于segment-graph
 *
*/
#pragma once
#ifndef DISJOINT_SET
#define DISJOINT_SET

// disjoint-set forests using union-by-rank and path compression (sort of).
typedef struct {
    int rank;
    int p;        // 该节点对应的根节点--并查集
    int size;
} Uni_elt;

typedef struct universe {
    Uni_elt *elts;
    int num;
    int size;     // Universe 所含元素数目
} Universe;

int get_universe_size(Universe *u, int x);

int get_universe_num_sets(Universe *u);

/*
* 初始化结构体 Universe
* 令Universe->num=元素数目
* 令Universe->elts中的所有元素
* rank=0, size=1, p = i;
* @param elements 表示有元素数目
*/
Universe * init_universe(const int elements);

/*
* 销毁Universe结构体，释放内存
*/
int destroy_universe(Universe *pu);

/*
* 查找节点x的根节点id
*/
int find_in_universe(Universe *pu, int x);

void join_universe(Universe *pu, int x, int y);

#endif /* !DISJOINT_SET */
