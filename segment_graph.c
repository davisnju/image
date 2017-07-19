/*!
* segment-graph.h
* date: 2017/03/21 15:02
*
* author: dav1sNJU
*
* brief:  Segment a graph
*
*
*/
#include "image.h"
#include "imconv.h"
#include "segment_graph.h"

/*
* Segment a graph
*
* Returns a disjoint-set forest representing the segmentation.
*
* num_vertices: number of vertices in graph.
* num_Edges: number of Edges in graph
* Edges: array of Edges.
* c: constant for treshold function.
*/
void segment_graph(Universe *u, int num_vertices, int num_Edges
                   , Edge *edges, float c)
{
    // sort Edges by weight
    qsort((void *)edges, (unsigned int)num_Edges, sizeof(Edge), edge_compare); //对边界的权重进行排序
       
    // init thresholds
    float *threshold = (float *)malloc(sizeof(float) * num_vertices);//在每个像素点设计一个threshold
    for (int i = 0; i < num_vertices; i++)
        threshold[i] = THRESHOLD(1, c);  // to set the threshold as the default value.

    // for each Edge, in non-decreasing weight order...
    //int ra = 0;// debugging 
    for (int i = 0; i < num_Edges; i++)
    {
        Edge *pEdge = &edges[i];

        // components conected by this Edge, 在这里可以调整REGION的大小。
        int roota = find_in_universe(u, pEdge->a);   //查找节点a的根节点
        int rootb = find_in_universe(u, pEdge->b);   //查找节点a的根节点
        if (roota != rootb)
        {
            if ((pEdge->w <= threshold[roota]) &&
                (pEdge->w <= threshold[rootb]))
            { // &&u->size(a)+u->size(b)<5000
                join_universe(u, roota, rootb); //合并
                roota = find_in_universe(u, roota); 
                threshold[roota] = pEdge->w + THRESHOLD(get_universe_size(u, roota), c); //更新threshold，如何做的。
            }
        }
    }

    // free up
    free(threshold);

    //return u;
};

