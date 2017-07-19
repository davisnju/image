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
#pragma once
#ifndef SEGMENT_GRAPH
#define SEGMENT_GRAPH

#include <stdlib.h>
#include <math.h>
#include "disjoint_set.h"
#include "edge.h"
#include "vector_c.h"

// threshold function
#define THRESHOLD(size, c) (c/size)

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
                   , Edge *Edges, float c);

#endif /* !SEGMENT_GRAPH */
