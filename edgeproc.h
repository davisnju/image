/*!
 * edgeproc.h
 * date: 2017/03/21 14:47
 *
 * author: dav1sNJU
 *
 * brief: 边界检测及相关处理
 *
 *
*/
#pragma once
#ifndef EDGEPROC_H
#define EDGEPROC_H

#include "element_lib.h"

#include "edge.h"
#include "border_band.h"

#define EPS 0.000002f
#define DEFAULT_EDGELIST_SIZE 10

/* 
 * 边界的连接信息
 * spurdegree 分支数
 * spurnode 分支节点（首或尾）
 * startnode 首节点号
 * sconns 首节点连接数
 * endnode 尾节点号
 * econns 尾节点连接数
 */
typedef struct connectionInfo
{
    int spurdegree, spurnode, startnode, sconns, endnode, econns;
}ConnectionInfo;

/*
% EDGELINK - Link edge points in an image into lists
% Usage: [edgelist edgeim, etypr] = edgelink(im, minlength, location)
% Arguments:  im         - Binary edge image, it is assumed that edges have been thinned (or are nearly thin).
%             minlength  - Optional minimum edge length of interest, defaults to 1 if omitted or specified as []. Ignored at the moment.
% Returns:    edgelist - a cell array of edge lists in row, column coords in the form
%                     { [r1 c1   [r1 c1   etc }
%                        r2 c2    ...
%                        ...
%                        rN cN]   ....]
%             edgeim   - Image with pixels labeled with edge number.
%             etype   - Array of values, one for each edge segment indicating its type
%                      0  - Start free, end free
%                      1  - Start free, end junction
%                      2  - Start junction, end free (should not happen)
%                      3  - Start junction, end junction
%                      4  - Loop
%
% This function links edge points together into lists of coordinate pairs.
% Where an edge junction is encountered the list is terminated and a separate
% list is generated for each of the branches.
*/
int edgelink(Image *im, EdgeLine *edgeList, int *edgesNum,
             Image_int *edgeIm, vector_int *etype, EdgeLine *longEdgeList, int *longEdgeNo,
             const int minlength);

int cleanedgelist(EdgeLine *edgeList, EdgeLine *longEdgeList, int edgeNo, const int minlength, int *longEdgeNo);


void findEndsJunctions(Image * im, vector_int *Rj, vector_int *Cj, vector_int *Re, vector_int *Ce);


void mergenodes(Mat8U *A, Mat8U *B, Image_int *Econnections, Image_int *Nconnections,
                EdgeLine *edgeList, int n1, int n2);

void removeedge(Mat8U *A, Mat8U *B, Image_int *Econnections, Image_int *Nconnections,
                EdgeLine *edgeList,
                int n);

void initNode(Image_long *node, EdgeLine *edgeList, int Nedges);
void initAdjacencyMatrix(Mat8U *A, Mat8U *B, Image_int *Econnections, Image_int *Nconnections,
    int Nnodes, int Nedges, Image_long *node);
//% TRACKEDGE
//%
//% Function to track all the edge points starting from an end point or junction.
//% As it tracks it stores the coords of the edge points in an array and labels the
//% pixels in the edge image with the - ve of their edge number.This continues
//% until no more connected points are found, or a junction point is encountered.
//%
//% Usage:   endType = trackedge(Image_f *imf, Image *junction, rstart, cstart, edgeNo
//                              ,EdgeLine *edgeList)
//%
//% Arguments : rstart, cstart - Row and column No of starting point.
//%              edgeNo - The current edge number.
//%                       should not be immediately connected to a junction(if possible).
//%
//% Returns:    edgePoints  - Nx2 array of row and col values for each edge point. (edgeList[edgeNo])
//%               endType - 0 for a free end
//%                         1 for a junction
//%                         5 for a loop
int trackedge(Image_f *imf, Image *junction, int rstart, int cstart, int edgeNo,
              int r2, int c2, int avoidJunction,
              EdgeLine *edgeList);


//% AVAILABLEPIXELS
//%
//% Find all the pixels that could be linked to point r, c
//%
//% Arguments:  rp, cp - Row, col coordinates of pixel of interest.
//%             edgeNo - The edge number of the edge we are seeking to
//%                      track.If not supplied its value defaults to 0
//% resulting in all adjacent junctions being returned,
//% (see note below)
//%
//% Returns:    ra, ca - Row and column coordinates of available non - junction
//%                      pixels.
//%             rj, cj - Row and column coordinates of available junction
//%                      pixels.
//%
//% A pixel is avalable for linking if it is :
//% 1) Adjacent, that is it is 8 - connected.
//% 2) Its value is 1 indicating it has not already been assigned to an edge
//% 3) or it is a junction that has not been labeled - edgeNo indicating we have
//%    not already assigned it to the current edge being tracked.If edgeNo is
//% 0 all adjacent junctions will be returned
//If edgeNo not supplied set to 0 to allow all adjacent junctions to be returned
void availablepixels(Image_f *imf, Image *junction, int rp, int cp, int edgeNo,
                     vector_int *ra, vector_int *ca, vector_int *rj, vector_int *cj);

/* find(B(n,:)); */
vector_int *findLink(Mat8U * B, int n);

/* find(A(n,:)); */
vector_int *findNodes(Mat8U * A, int node);
int isJunction(int *x);
int isEnding(int *x);
ConnectionInfo *getConnectionInfo(Image_int *Nconnections, EdgeLine *edgeList, int n);

// Function to compute the path length of an edgelist
float getEdgeLineLength(EdgeLine * edgelist);

vector_int * intersect_row(vector_int * ra, vector_int * ca, vector_int * rak, vector_int * cak);

float norm_t(int r, int c);

void unitvector(float r, float c, float *or, float *oc);



// ============================ 暂时未使用 =================================
/*
Image *edge_band(Image_yuv * rim, Image_f *lab_im, Image_f *mag, EdgeVec **edgeList, unsigned int edgesNum,
    float *edgeDiff, int *edgeDiffSize,
    const float lenthThred);
*/
#endif /* !EDGEPROC_H */
