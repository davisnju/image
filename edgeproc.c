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
#include "edgeproc.h"
#include "imconv.h"
#include "improc.h"
#include <stdlib.h>
#include <stdint.h>
#include <math.h>
#include <float.h>
#include <assert.h>

//#define TIMECOST_TEST 

#ifdef TIMECOST_TEST
#include <time.h>
#include <stdio.h>
#endif /* TIMECOST_TEST */

int edgelink(Image *im, EdgeLine *edgeList, int *edgesNum,
             Image_int *edgeIm, vector_int *etype, EdgeLine *longEdgeList, int *longEdgeNo,
             const int minlength)
{
    int w = im->w, h = im->h;
    binarizeImage(im, w, h); // 确保输入im为二值图像
    cleanImage(im, w, h);    // 去除孤立像素
    thinImage(im, w, h);    // 腐蚀，thins objects to lines.

    // Find endings and junctions in edge data 
    vector_int *Rj, *Cj, *Re, *Ce;  
    int vecsize = min(1000, max(20, max(w, h)));   //20 <= vecsize <= 1000
    create_vector_int(&Rj, vecsize);
    create_vector_int(&Cj, vecsize);
    create_vector_int(&Re, vecsize);
    create_vector_int(&Ce, vecsize);

    findEndsJunctions(im, Rj, Cj, Re, Ce);

    int Njunct = getSize_vector_int(Rj), Nends = getSize_vector_int(Re);
    Image *junction = create_new_image(w, h, 1);  
    int ri = 0, ci = 0, Rjj, Cjj, edgeNo = 0;;
    for (int i = 0; i < Njunct; ++i)
    {
        ri = Rj->head[i]/*at_vector_int(Rj, i)*/;
        ci = Cj->head[i]/*at_vector_int(Cj, i)*/;
        imRef(junction, ci, ri) = 1u;
    }

    Image_f * imf = create_new_imagef(w, h, 0);
    imageUINT8_TtoFLOAT(im, imf);

    for (int i = 0; i < Nends; ++i)
    {
        ri = Re->head[i]/*at_vector_int(Re, i)*/;
        ci = Ce->head[i]/*at_vector_int(Ce, i)*/;
        if (fabs(imRef(imf, ci, ri) - 1) < EPS)
        {
            ++edgeNo;
            int endType = trackedge(imf, junction, ri, ci, edgeNo, -1, -1, 0, edgeList);  // 返回[edgelist{edgeNo} endType]
        }
    }
    for (int j = 0; j < Njunct; ++j)
    {
        Rjj = at_vector_int(Rj, j);
        Cjj = at_vector_int(Cj, j);
        if (imRef(junction, Cjj, Rjj) != 2u) // We have not visited this junction
        {
            imRef(junction, Cjj, Rjj) = 2u;
            vector_int *rj, *cj, *ra, *ca;   // free ok
            create_vector_int(&rj, 10);
            create_vector_int(&cj, 10);
            create_vector_int(&ra, 10);
            create_vector_int(&ca, 10);
            availablepixels(imf, junction, Rjj, Cjj, 0, ra, ca, rj, cj);

            int rjLen = getSize_vector_int(rj);
            for (int k = 0; k < rjLen; ++k)// For all adjacent junctions...
            {
                int rjk = at_vector_int(rj, k), cjk = at_vector_int(cj, k);
                ++edgeNo;
                Point ejj = { Rjj, Cjj };
                Point ejk = { rjk, cjk };
                initEdgeLine(&edgeList[edgeNo - 1]);
                appendPointToEdgeLine(&edgeList[edgeNo - 1], &ejj);
                appendPointToEdgeLine(&edgeList[edgeNo - 1], &ejk);

                imRef(imf, Cjj, Rjj) = (float)-edgeNo;
                imRef(imf, cjk, rjk) = (float)-edgeNo;
                vector_int *rak = NULL, *cak = NULL;   
                create_vector_int(&rak, 10);
                create_vector_int(&cak, 10);
                availablepixels(imf, junction, rjk, cjk, 0, rak, cak, NULL, NULL);
                if (getSize_vector_int(ra) > 0 && getSize_vector_int(rak) > 0)
                {
                    // 相同点[r0,c0,r1,c1,...]
                    vector_int *commonrc = intersect_row(ra, ca, rak, cak);  
                    int commonrcLen = getSize_vector_int(commonrc);
                    for (int n = 0; n < commonrcLen; n += 2)
                    {
                        edgeNo = edgeNo + 1;
                        int crcn_r = at_vector_int(commonrc, n),
                            crcn_c = at_vector_int(commonrc, n + 1);
                        float distj = norm_t(crcn_r - Rjj, crcn_c - Cjj);
                        float distk = norm_t(crcn_r - rjk, crcn_c - cjk);
                        if (distj < distk)
                            trackedge(imf, junction, Rjj, Cjj, edgeNo, crcn_r, crcn_c, 1, edgeList);
                        else
                            trackedge(imf, junction, rjk, cjk, edgeNo, crcn_r, crcn_c, 1, edgeList);
                    }
                    destroy_vector_int(commonrc);
                }

                int rakLen = getSize_vector_int(rak);
                for (int m = 0; m < rakLen; ++m) //m = 1:length(rak)
                {
                    int rakm = at_vector_int(rak, m), cakm = at_vector_int(cak, m);
                    if (fabs(imRef(imf, cakm, rakm) - 1) < EPS)//EDGEIM(rak(m), cak(m)) == 1
                    {
                        ++edgeNo;
                        trackedge(imf, junction, Rjj, Cjj, edgeNo, rakm, cakm, 0, edgeList);
                    }
                }
                // Mark that we have visited junction(rj(k) cj(k))
                imRef(junction, cjk, rjk) = 2u;

                destroy_vector_int(rak);
                destroy_vector_int(cak);
            }

            // Finally track any remaining unlabeled pixels adjacent to original junction j
            int raLen = getSize_vector_int(ra);
            for (int m = 0; m < raLen; ++m) //m = 1:length(ra)
            {
                int ram = at_vector_int(ra, m), cam = at_vector_int(ca, m);
                if (fabs(imRef(imf, cam, ram) - 1) < EPS)//EDGEIM(ra(m), ca(m)) == 1
                {
                    ++edgeNo;
                    trackedge(imf, junction, Rjj, Cjj, edgeNo, ram, cam, 0, edgeList);
                }
            }
            destroy_vector_int(rj);
            destroy_vector_int(cj);
            destroy_vector_int(ra);
            destroy_vector_int(ca);
        }// If we have not visited this junction
    }//  For each junction

    for (int x = 0; x < w; x++)
    {
        for (int y = 0; y < h; y++)
        {
            if (fabs(imRef(imf, x, y) - 1) < EPS) // We have an unlabeled edge
            {
                ++edgeNo;
                int endType = trackedge(imf, junction, y, x, edgeNo, -1, -1, 0, edgeList);  // 返回[edgelist{edgeNo} endType]
            }
        }
    }
    if (edgeIm)
    {
        for (int x = 0; x < w; x++)
            for (int y = 0; y < h; y++)
                imRef(edgeIm, x, y) = (int)round(-imRef(imf, x, y));
    }

    *edgesNum = edgeNo;

    if (minlength >= 0 && edgeNo > 0)
        cleanedgelist(edgeList/*tmpEdgeLine*/, longEdgeList, edgeNo, minlength, longEdgeNo);

    destroy_vector_int(Rj);
    destroy_vector_int(Cj);
    destroy_vector_int(Re);
    destroy_vector_int(Ce);
    destroy_image(junction);
    free(junction);
    destroy_imagef(imf);
    free(imf);

    if (edgeNo < 0)  return -1;
    return 0;
};

//% TRACKEDGE
//%
//% Function to track all the edge points starting from an end point or junction.
//% As it tracks it stores the coords of the edge points in an array and labels the
//% pixels in the edge image with the - ve of their edge number.This continues
//% until no more connected points are found, or a junction point is encountered.
//%
//% Usage:   endType = trackedge(Image_f *imf, Image *junction, rstart, cstart, edgeNo
//                              ,EdgeVec *edgeList)
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
              EdgeLine *edgeList)
{
    int endType = 0;
    if (avoidJunction < 0) avoidJunction = 0;

    Point edgep = { 0,0 };
    Point *edgePoint = &edgep;

    int vecsize = max(imf->w + imf->h, 20);

    initEdgeLine(&edgeList[edgeNo - 1]);
    EdgeLine *edgePoints = &edgeList[edgeNo - 1]; // = (vector_int *)malloc(sizeof(vector_int));

    edgePoint->x = cstart;
    edgePoint->y = rstart;
    appendPointToEdgeLine(edgePoints, edgePoint);

    imRef(imf, cstart, rstart) = (float)-edgeNo;

    int preferredDirection = 0;

    float dirn[] = { 0, 0 }, dirna[] = { 0, 0 }, dirnbest[] = { 0, 0 };
    int rbest, cbest;
    int r = rstart;
    int c = cstart;

    if (r2 >= 0 && c2 >= 0)
    {
        edgePoint->x = c2;
        edgePoint->y = r2;
        appendPointToEdgeLine(edgePoints, edgePoint);
        imRef(imf, c2, r2) = (float)-edgeNo;
        //Initialise direction vector of path and set the current point on the path
        unitvector((float)(r2 - rstart), (float)(c2 - cstart), &dirn[0], &dirn[1]);
        r = r2;
        c = c2;
        preferredDirection = 1;
    }

    vector_int *rj = NULL, *cj = NULL, *ra = NULL, *ca = NULL;
    create_vector_int(&rj, 10);
    create_vector_int(&cj, 10);
    create_vector_int(&ra, 10);
    create_vector_int(&ca, 10);
    availablepixels(imf, junction, r, c, edgeNo, ra, ca, rj, cj);

    float dotp = 0.0f;
    while ((!empty_vector_int(ra)) || (!empty_vector_int(rj)))
    {
        if ((!empty_vector_int(rj) && !avoidJunction)
            || (!empty_vector_int(rj) && empty_vector_int(ra)))
        {
            if (preferredDirection)
            {
                dotp = -FLT_MAX;
                int rjLen = rj->ele_num;/* getSize_vector_int(rj);*/
                for (int n = 0; n < rjLen; ++n) // n = 1:length(rj)
                {
                    int rjn = rj->head[n]/*at_vector_int(rj, n)*/, cjn = cj->head[n];/* at_vector_int(cj, n);*/
                    unitvector((float)(rjn - r), (float)(cjn - c), &dirna[0], &dirna[1]);
                    float dp = dirn[0] * dirna[0] + dirn[1] * dirna[1];//dirn*dirna';
                    if (dp > dotp)
                    {
                        dotp = dp;
                        rbest = rjn;
                        cbest = cjn;
                        dirnbest[0] = dirna[0];
                        dirnbest[1] = dirna[1];
                    }
                }
            }
            else
            {
                float dirnbest2 = FLT_MAX;
                int rjLen = rj->ele_num;/* getSize_vector_int(rj);*/
                for (int n = 0; n < rjLen; ++n) // n = 1:length(rj)
                {
                    int rjn = rj->head[n]/*at_vector_int(rj, n)*/, cjn = cj->head[n];/* at_vector_int(cj, n);*/
                    float dist = (float)rjn - r + cjn - c;//sum([rj(n)-r;  cj(n)-c]); 
                    if (dist < dirnbest2)
                    {
                        rbest = rjn;
                        cbest = cjn;
                        dirnbest2 = dist;
                        unitvector((float)rjn - r, (float)cjn - c, &dirnbest[0], &dirnbest[1]);
                    }
                }
                preferredDirection = 1;
            }
        }
        else
        {
            dotp = -FLT_MAX;
            int raLen = ra->ele_num;
            for (int m = 0; m < raLen; ++m) //m = 1:length(ra)
            {
                int ram = ra->head[m], cam = ca->head[m];
                unitvector((float)ram - r, (float)cam - c, &dirna[0], &dirna[1]);
                float dp = dirn[0] * dirna[0] + dirn[1] * dirna[1];//dirn*dirna';
                if (dp > dotp)
                {
                    dotp = dp;
                    rbest = ram;
                    cbest = cam;
                    dirnbest[0] = dirna[0];
                    dirnbest[1] = dirna[1];
                }
            }
            avoidJunction = 0;  //Clear the avoidJunction flag if it had been set
        }
        r = rbest;
        c = cbest;
        edgePoint->x = c;
        edgePoint->y = r;
        appendPointToEdgeLine(edgePoints, edgePoint);
        dirn[0] = dirnbest[0];
        dirn[1] = dirnbest[1];
        imRef(imf, c, r) = (float)-edgeNo;

        // If this point is a junction exit here
        if (imRef(junction, c, r) > 0u)
        {
            endType = 1;  // Mark end as being a junction
            break;
        }
        else
        {    // Get the next set of available pixels to link.
            availablepixels(imf, junction, r, c, edgeNo, ra, ca, rj, cj);
        }
    } // while()
    
    int edgesPointsSize = edgePoints->ele_num;
    if (edgesPointsSize / 2 >= 4)
    {
        int rs = edgePoints->head->pos.y,
            cs = edgePoints->head->pos.x,
            re = edgePoints->end->pos.y, 
            ce = edgePoints->end->pos.x;
        if (abs(rs - re) <= 1 && abs(cs - ce) <= 1)
        {
            edgePoint->x = c;
            edgePoint->y = r;
            appendPointToEdgeLine(edgePoints, edgePoint);
            endType = 5; // Mark end as being a loop
        }
    }

    destroy_vector_int(rj);
    destroy_vector_int(cj);
    destroy_vector_int(ra);
    destroy_vector_int(ca);
    return endType;
};

void initNode(Image_long *node, EdgeLine *edgeList, int Nedges)
{
    Edgeline_elt *pes = NULL, *pee = NULL;
    int pos = 0;
    for (int i = 0; i < Nedges; i++, pos += 2)
    {
        int eilen = edgeList[i].ele_num;/* getEdgeVecSize(edgeList[i]);*/
        if (eilen == 0)
            continue;
        pes = edgeList[i].head;/* atEdgeVec(edgeList[i], 0);*/
        imRef(node, 0, pos) = pes->pos.y;
        imRef(node, 1, pos) = pes->pos.x;
        pee = edgeList[i].end; /*backInEdgelist(edgeList[i]);*/
        imRef(node, 0, pos + 1) = pee->pos.y;
        imRef(node, 1, pos + 1) = pee->pos.x;
    }
};

void initAdjacencyMatrix(Mat8U *A, Mat8U *B, Image_int *Econnections, Image_int *Nconnections,
    int Nnodes, int Nedges, Image_long *node)
{
    for (int n = 0; n < Nnodes - 1; n++)
    {
        int edgen = n >> 1;
        for (int m = n + 1; m < Nnodes; m++)
        {
            // If nodes m & n are connected
            imRef(A, m, n) = imRef(node, 0, n) == imRef(node, 0, m) && imRef(node, 1, n) == imRef(node, 1, m);
            imRef(A, n, m) = imRef(A, m, n);
            if (imRef(A, m, n) != 0)
            {
                int edgem = m >> 1;
                imRef(B, edgem, edgen) = 1;
                imRef(B, edgen, edgem) = 1;
            }
        }
    }

    int sum = 0;
    for (int n = 0; n < Nnodes; n++)
    {
        sum = 0;
        for (int m = 0; m < Nnodes; m++)
            sum += imRef(A, m, n);
        imRef(Nconnections, 0, n) = sum;
    }
    for (int edgen = 0; edgen < Nedges; edgen++)
    {
        sum = 0;
        for (int m = 0; m < Nedges; m++)
            sum += imRef(B, m, edgen);
        imRef(Econnections, 0, edgen) = sum;
    }
};

static void getLongEdgeLine(EdgeLine *edgeList, EdgeLine *longEdgeList, int edgeNo, const int minlength, int *longEdgeNo)
{
    int m = 0;
    for (int n = 0; n < edgeNo; n++)
    {
        int l = edgeList[n].ele_num;
        if (l > 0)
        {
            float llength = (float)getEdgeLineLength(&edgeList[n]);
            if (llength >= minlength)
            {
                initEdgeLine(&longEdgeList[m]);
                longEdgeList[m].head = edgeList[n].head;
                longEdgeList[m].end = edgeList[n].end;
                longEdgeList[m].ele_num = edgeList[n].ele_num;
                ++m;
            }
        }
    }    
    *longEdgeNo = m;
};

static void handleOneConns(Mat8U *A, Mat8U *B, Image_int *Econnections, Image_int *Nconnections,
    int Nnodes, int Nedges, EdgeLine *edgeList)
{
    for (int n = 0; n < Nedges; n++)
    {
        int edgenSize = edgeList[n].ele_num;
        if (!(imRef(B, n, n)) && edgenSize > 0)
        {
            ConnectionInfo *ci = getConnectionInfo(Nconnections, edgeList, n);
            if (ci->sconns == 1)
            {
                vector_int *node2merge = findNodes(A, ci->startnode);//find(A(startnode,:));
                int node2Num = node2merge->ele_num;

                for (int i = 0; i < node2Num; i++)
                    mergenodes(A, B, Econnections, Nconnections, edgeList, node2merge->head[i], ci->startnode);
                destroy_vector_int(node2merge);
            }
            edgenSize = edgeList[n].ele_num;
            if (edgenSize > 0 && ci->econns == 1)
            {
                vector_int *node2merge = findNodes(A, ci->endnode);// find(A(endnode, :));
                int node2Num = node2merge->ele_num;
                for (int i = 0; i < node2Num; i++)
                    mergenodes(A, B, Econnections, Nconnections, edgeList, node2merge->head[i], ci->endnode);
                destroy_vector_int(node2merge);
            }
            free(ci);
        }
    }
};
static void handleUnvisitedEdges(Mat8U *A, Mat8U *B, Image_int *Econnections, Image_int *Nconnections,
    int Nnodes, int Nedges, EdgeLine *edgeList, const int minlength)
{
    for (int n = 0; n < Nedges; n++)
    {
        int l = edgeList[n].ele_num;
        if (l > 0 && getEdgeLineLength(&edgeList[n]) < minlength
            && (!imRef(Econnections, 0, n)
                || (imRef(Econnections, 0, n) == 1 && imRef(B, n, n) == 1)))
        {
            removeedge(A, B, Econnections, Nconnections, edgeList, n);
        }
    }
};

static void mergeTwoSpurDegreeNodes(Mat8U *A, Mat8U *B, Image_int *Econnections, Image_int *Nconnections,
     EdgeLine *edgeList,int n, vector_int *linkingedges, const unsigned int ll)
{
    vector_int *spurs = NULL;
    create_vector_int(&spurs, 10);
    push_back_vector_int(spurs, n);

    float edgeLength = getEdgeLineLength(&edgeList[n]);//
    int edgeSize = edgeList[n].ele_num;
    vector_int *len = NULL;
    create_vector_int(&len, 10);
    push_back_vector_int(len, (int)(edgeLength + 0.5));
    for (unsigned int i = 0; i < ll; i++)
    {
        int linkingedgei = at_vector_int(linkingedges, i);//linkingedges(i)
        ConnectionInfo *ci2 = getConnectionInfo(Nconnections, edgeList, linkingedgei);
        if (ci2->spurdegree != 0)
        {
            push_back_vector_int(spurs, linkingedgei);
            push_back_vector_int(len, (int)(getEdgeLineLength(&edgeList[n]) + 0.5));
        }
        free(ci2);
    }
    push_back_vector_int(linkingedges, n);
    int idx = -1, minlen = -1;
    int lenLen = getSize_vector_int(len);
    for (int i = 0; i < lenLen; i++)
    {
        if (minlen < 0 || at_vector_int(len, i) < minlen)
        {
            idx = i;
            minlen = at_vector_int(len, i);
        }
    }
    int edge2delete = at_vector_int(spurs, idx);
    ConnectionInfo *ci3 = getConnectionInfo(Nconnections, edgeList, edge2delete);
    vector_int *nodes2merge = findNodes(A, ci3->spurnode);
    
    assert((getSize_vector_int(nodes2merge) == 2u));// error('attempt to merge other than 2 nodes');

    removeedge(A, B, Econnections, Nconnections, edgeList, edge2delete);
    mergenodes(A, B, Econnections, Nconnections, edgeList, nodes2merge->head[0], nodes2merge->head[1]);

    destroy_vector_int(spurs);
    destroy_vector_int(len);
    destroy_vector_int(nodes2merge);
};

int cleanedgelist(EdgeLine *edgeList, EdgeLine *longEdgeList, int edgeNo, const int minlength, int *longEdgeNo)
{    
    assert(edgeNo > 0);
    assert(edgeList);
    assert(longEdgeList);
    int Nnodes = edgeNo << 1, Nedges = edgeNo;
    Image_long *node = create_new_image_long(2, Nnodes, 0);
    initNode(node, edgeList, Nedges);

    Mat8U *A = create_new_image(Nnodes, Nnodes, 1),  // Adjacency matrix for nodes
        *B = create_new_image(Nedges, Nedges, 1);     // Adjacency matrix for edges

    Image_int *Nconnections = create_new_image_int(1, Nnodes, 1),  // Connection count array for nodes
        *Econnections = create_new_image_int(1, Nedges, 1);     // Connection count array for edges

    initAdjacencyMatrix(A, B, Econnections, Nconnections, Nnodes, Nedges, node);
    destroy_image_long(node);
    free(node);


    handleOneConns(A, B, Econnections, Nconnections, Nnodes, Nedges, edgeList);

    if (minlength >= 0)
    {
        for (int n = 0; n < Nedges; n++)
        {
            ConnectionInfo *ci = getConnectionInfo(Nconnections, edgeList, n);
            int l = edgeList[n].ele_num;
            if (l > 0 && getEdgeLineLength(&edgeList[n]) < minlength)
            {
                if (!imRef(Econnections, 0, n) ||
                    (imRef(Econnections, 0, n) == 1 && imRef(B, n, n) == 1))
                    removeedge(A, B, Econnections, Nconnections, edgeList, n);
                else if (ci->spurdegree == 2)
                {
                    vector_int *linkingedges = findLink(B, n);
                    unsigned int ll = getSize_vector_int(linkingedges);
                    if (ll == 1)
                        removeedge(A, B, Econnections, Nconnections, edgeList, n);
                    else
                        mergeTwoSpurDegreeNodes(A, B, Econnections, Nconnections, edgeList, n, linkingedges, ll);
                    destroy_vector_int(linkingedges);
                }
                else if (ci->spurdegree == 3)
                    removeedge(A, B, Econnections, Nconnections, edgeList, n);

            }
            free(ci);
        } //for (int n = 0; n < Nedges; n++)

        handleUnvisitedEdges(A, B, Econnections, Nconnections, Nnodes, Nedges, edgeList, minlength);
    } // if (minlength >= 0)

    getLongEdgeLine(edgeList, longEdgeList, edgeNo, minlength, longEdgeNo);

    destroy_image(A);
    destroy_image(B);
    destroy_image_int(Nconnections);
    destroy_image_int(Econnections);
    free(A);
    free(B);
    free(Nconnections);
    free(Econnections);

    return 0;
};

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
                     vector_int *ra, vector_int *ca, vector_int *rj, vector_int *cj)
{
    // row and column offsets for the eight neighbours of a point
    clear_vector_int(ra);
    clear_vector_int(ca);
    clear_vector_int(rj);
    clear_vector_int(cj);

    if (edgeNo < 0) edgeNo = 0;

    int roff[] = { -1, 0, 1, 1, 1, 0, -1, -1 },
        coff[] = { -1, -1, -1, 0, 1, 1, 1, 0 },
        rarr[8] = { 0 }, carr[8] = { 0 },
        w = imf->w, h = imf->h;
    vector_int *ind = NULL;
    create_vector_int(&ind, 10);
    for (int i = 0; i < 8; i++)
    {
        int r = rp + roff[i], c = cp + coff[i];
        rarr[i] = r;
        carr[i] = c;
        if (r >= 0 && r < h && c >= 0 && c < w)
            push_back_vector_int(ind, i);
    }

    int indLen = getSize_vector_int(ind);
    for (int i = 0; i < indLen; ++i)
    {
        int pos = at_vector_int(ind, i);
        if (pos < 0 || pos > 7) //something is wrong
            continue;

        int ci = carr[pos], ri = rarr[pos];
        if ((fabs(imRef(imf, ci, ri) - 1) < EPS) && !imRef(junction, ci, ri))
        {
            push_back_vector_int(ra, ri);//ra = [ra; r(i)];
            push_back_vector_int(ca, ci);//ca = [ca; c(i)];
        }
        else if (rj && cj
                 && (fabs(imRef(imf, ci, ri) + edgeNo) > EPS)
                 && imRef(junction, ci, ri))
        {
            push_back_vector_int(rj, ri);//rj = [rj; r(i)];
            push_back_vector_int(cj, ci);//cj = [cj; c(i)];            
        }
    }
    destroy_vector_int(ind);
};

//Recompute connection counts
static void recompute_connection_counts(Mat8U *A, Mat8U *B, Image_int *Econnections, Image_int *Nconnections,
                                int Nnodes, int Nedges, int edge1, int edge2)
{
    for (int i = 0; i < Nedges; i++)
    {
        imRef(B, edge1, i) = imRef(B, edge1, i) | imRef(B, edge2, i);
        imRef(B, i, edge1) = imRef(B, i, edge1) | imRef(B, i, edge2);
    }
    imRef(B, edge1, edge1) = 0;

    int sum = 0;
    for (int n = 0; n < Nnodes; n++)
    {
        sum = 0;
        for (int m = 0; m < Nnodes; m++)
            sum += imRef(A, m, n);
        imRef(Nconnections, 0, n) = sum;
    }
    for (int edgen = 0; edgen < Nedges; edgen++)
    {
        sum = 0;
        for (int m = 0; m < Nedges; m++)
            sum += imRef(B, m, edgen);
        imRef(Econnections, 0, edgen) = sum;
    }
};

static void setA(Mat8U *A, int Nnodes, int e1, int e2)
{
    for (int i = 0; i < Nnodes; i++)
    {
        imRef(A, e1, i) = imRef(A, e2, i);
        imRef(A, i, e1) = imRef(A, i, e2);
    }
};

//% Internal function to merge 2 edgelists together at the specified nodes and
//% perform the necessary updates to the edge adjacency and node adjacency
//% matrices and the connection count arrays
void mergenodes(Mat8U *A, Mat8U *B, Image_int *Econnections, Image_int *Nconnections,
                EdgeLine *edgeList, int n1, int n2)
{
    int edge1 = n1 >> 1, edge2 = n2 >> 1;// Indices of the edges associated with the nodes
    if (edge1 == edge2) return;
    int Nnodes = A->w, Nedges = B->w;
    // Get indices of nodes at each end of the two edges
    int s1 = 2 * edge1, e1 = 2 * edge1 + 1,
        s2 = 2 * edge2, e2 = 2 * edge2 + 1;

    if (!imRef(A, n1, n2))  return;// error('Attempt to merge nodes that are not connected');
    int flipedge1 = 0, flipedge2 = 1;
    if (n1 % 2 == 0)    flipedge1 = 1; // node n1 is the start of edge1 // edge1 will need to be reversed in order to join edge2
    if (n2 % 2 == 0)  flipedge2 = 0;  //node n2 is the start of edge2
    
    if (!flipedge1 && !flipedge2)
    {
        appendEdgeLineToEdgeLine(&edgeList[edge1], &edgeList[edge2]);
        setA(A, Nnodes, e1, e2);
        imRef(Nconnections, 0, e1) = imRef(Nconnections, 0, e2);
    }
    else if (!flipedge1 && flipedge2)
    {
        reverseEdgeLine(&edgeList[edge2]);
        appendEdgeLineToEdgeLine(&edgeList[edge1], &edgeList[edge2]);
        setA(A, Nnodes, e1, s2);
        imRef(Nconnections, 0, e1) = imRef(Nconnections, 0, s2);
    }
    else if (flipedge1 && !flipedge2)
    {
        reverseEdgeLine(&edgeList[edge1]);
        appendEdgeLineToEdgeLine(&edgeList[edge1], &edgeList[edge2]);
        setA(A, Nnodes, s1, e1);
        setA(A, Nnodes, e1, e2);
        imRef(Nconnections, 0, s1) = imRef(Nconnections, 0, e1);
        imRef(Nconnections, 0, e1) = imRef(Nconnections, 0, e2);
    }
    else if (flipedge1 && flipedge2)
    {
        reverseEdgeLine(&edgeList[edge1]);
        reverseEdgeLine(&edgeList[edge2]);
        appendEdgeLineToEdgeLine(&edgeList[edge1], &edgeList[edge2]);
        setA(A, Nnodes, s1, e1);
        setA(A, Nnodes, e1, s2);
        imRef(Nconnections, 0, s1) = imRef(Nconnections, 0, e1);
        imRef(Nconnections, 0, e1) = imRef(Nconnections, 0, s2);
    }
    else return;  //error('We should not have got here - edgelists cannot be merged');  
    recompute_connection_counts(A, B, Econnections, Nconnections, Nnodes, Nedges, edge1, edge2);
    removeedge(A, B, Econnections, Nconnections, edgeList, edge2);  // Finally discard edge2
};

void removeedge(Mat8U *A, Mat8U *B, Image_int *Econnections, Image_int *Nconnections,
                EdgeLine *edgeList,
                int n)
{
    //clearEdgeLine(&edgeList[n]);
    int Nnodes = A->w, Nedges = B->w;
    for (int i = 0; i < Nedges; i++)
    {
        imRef(Econnections, 0, i) -= imRef(B, n, i);
        imRef(B, n, i) = 0;
        imRef(B, i, n) = 0;
    }
    imRef(Econnections, 0, n) = 0;

    int nodes2delete[] = { 2 * n, 2 * n + 1 };

    for (int i = 0; i < Nnodes; i++)
    {
        imRef(Nconnections, 0, i) -= imRef(A, nodes2delete[0], i);
        imRef(Nconnections, 0, i) -= imRef(A, nodes2delete[1], i);
        imRef(A, i, nodes2delete[0]) = 0;
        imRef(A, nodes2delete[0], i) = 0;
        imRef(A, i, nodes2delete[1]) = 0;
        imRef(A, nodes2delete[1], i) = 0;
    }
}

/*
% FINDENDSJUNCTIONS - find junctions and endings in a line/edge image
%
% Usage: [rj, cj, re, ce] = findendsjunctions(edgeim, disp)
%
% Arguments:  edgeim - A binary image marking lines/edges in an image.  It is
%                      assumed that this is a thinned or skeleton image
% Returns:    rj, cj - Row and column coordinates of junction points in the image.
%             re, ce - Row and column coordinates of end points in the image.
*/
void findEndsJunctions(Image * im, vector_int *Rj, vector_int *Cj, vector_int *Re, vector_int *Ce)
{
    int w = im->w, h = im->h;
    int pixels[9] = { 0 };
    int roff[] = { -1, 0, 1, -1, 0, 1, -1, 0, 1 },
        coff[] = { -1, -1, -1, 0, 0, 0, 1, 1, 1 };
    for (int y = 0; y < h; y++)
    {
        for (int x = 0; x < w; x++)
        {

            if (imRef(im, x, y) == 0u)
                continue;

            for (int i = 0; i < 9; i++)
            {
                int c = x + coff[i], r = y + roff[i];
                if (r >= 0 && r < h && c >= 0 && c < w)
                    pixels[i] = imRef(im, c, r);
                else
                    pixels[i] = 0;
            }
            if (isJunction(pixels)) // find Junctions
            {
                push_back_vector_int(Rj, y);
                push_back_vector_int(Cj, x);
            }
            else if (isEnding(pixels)) // find Ends
            {
                push_back_vector_int(Re, y);
                push_back_vector_int(Ce, x);
            }

        }
    }
};

/* find(B(n,:)); */
vector_int * findLink(Mat8U * B, int n)
{

    int Nedges = B->w;
    vector_int *link = NULL;
    create_vector_int(&link, Nedges / 2 + 10);
    for (int i = 0; i < Nedges; i++)
        if (imRef(B, n, i) > 0)
            push_back_vector_int(link, i);

    return link;
};

/* find(A(startnode,:)); */
vector_int * findNodes(Mat8U * A, int node)
{
    int Nnodes = A->w;
    vector_int *nodes = NULL;
    create_vector_int(&nodes, Nnodes / 2 + 10);
    for (int i = 0; i < Nnodes; i++)
        if (imRef(A, node, i) > 0)
            push_back_vector_int(nodes, i);

    return nodes;
};

// Function to compute the path length of an edgelist
float getEdgeLineLength(EdgeLine * edgelist)
{
    unsigned int size = edgelist->ele_num;
    float len = 0;
    Edgeline_elt *pElt = edgelist->head;
    while (pElt != edgelist->end)
    {
        len += (float)sqrt(square(pElt->pos.x - pElt->next->pos.x) + square(pElt->pos.y - pElt->next->pos.y));
        pElt = pElt->next;
    }


    return len;
};

ConnectionInfo *getConnectionInfo(Image_int *Nconnections, EdgeLine *edgeList, int n)
{
    ConnectionInfo *ci = (ConnectionInfo *)malloc(1 * sizeof(ConnectionInfo));
    ci->spurdegree = 0;
    ci->spurnode = 0;
    ci->startnode = 0;
    ci->sconns = 0;
    ci->endnode = 0;
    ci->econns = 0;

    if (edgeList[n].ele_num == 0u)
        return ci;
    ci->startnode = 2 * n;
    ci->endnode = 2 * n + 1;
    ci->sconns = imRef(Nconnections, 0, ci->startnode);  // No of connections to start node
    ci->econns = imRef(Nconnections, 0, ci->endnode);    //% No of connections to end node

    if (ci->sconns == 0 && ci->econns >= 1)
    {
        ci->spurdegree = ci->econns;
        ci->spurnode = ci->endnode;
    }
    else if (ci->sconns >= 1 && ci->econns == 0)
    {
        ci->spurdegree = ci->sconns;
        ci->spurnode = ci->startnode;
    }
    else // 首尾都存在分支
    {
        ci->spurdegree = 0;
        ci->spurnode = -1;
    }

    return ci;
}

int isEnding(int *x)
{
    int a[] = { 0, 1, 2, 5, 8, 7, 6, 3, 0 };
    int crossings = 0;
    for (int i = 0; i < 8; i++)
        crossings += abs(x[a[i]] - x[a[i + 1]]);
    return x[4] && crossings == 2;
};
/*
% Pixels in the 3x3 region are numbered as follows
%
%       1 4 7
%       2 5 8
%       3 6 9
*/
int isJunction(int *x)
{
    int a[] = { 0, 1, 2, 5, 8, 7, 6, 3, 0 };
    int crossings = 0;
    for (int i = 0; i < 8; i++)
        crossings += abs(x[a[i]] - x[a[i + 1]]);
    return x[4] && crossings >= 6;
};

/* 对排好序的向量求相同子序列 */
vector_int *intersect_row(vector_int *ra, vector_int *ca, vector_int *rak, vector_int *cak)
{
    vector_int *commonrc = NULL;
    create_vector_int(&commonrc, 10);

    unsigned int raLen = getSize_vector_int(ra), caLen = getSize_vector_int(ca),
        rakLen = getSize_vector_int(rak), cakLen = getSize_vector_int(cak);
    raLen = min(raLen, caLen);
    rakLen = min(rakLen, cakLen);
    vector_point *pointsVecA = NULL, *pointsVecB = NULL;   //free ok
    create_vector_point(&pointsVecA, raLen + 5);
    create_vector_point(&pointsVecB, rakLen + 5);

    for (unsigned int i = 0; i < raLen; i++)
    {
        int rai = ra->head[i], cai = ca->head[i];
        Point pti = { cai, rai };
        push_back_vector_point(pointsVecA, &pti);
    }
    for (unsigned int i = 0; i < rakLen; i++)
    {
        int raki = rak->head[i], caki = cak->head[i];
        Point pti = { caki, raki };
        push_back_vector_point(pointsVecB, &pti);
    }
    qsort(pointsVecA->head, raLen, sizeof(Point), point_cmpfunc);
    qsort(pointsVecB->head, rakLen, sizeof(Point), point_cmpfunc);

    for (unsigned int i = 0; i < raLen; i++)
    {
        Point *ptAi = &pointsVecA->head[i];
        if (find_in_vector_point(pointsVecB, ptAi) > -1)
        {
            push_back_vector_int(commonrc, ptAi->y); // r
            push_back_vector_int(commonrc, ptAi->x); // c
        }
    }
    destroy_vector_point(pointsVecA);
    destroy_vector_point(pointsVecB);
    return commonrc;
};

float norm_t(int r, int c)
{
    return (float)sqrt(square(r) + square(c));
};

void unitvector(float r, float c, float *or, float *oc)
{
    double sqrtrc = sqrt(fsquare(r) + fsquare(c));
    *or = (float)(r / sqrtrc);
    *oc = (float)(c / sqrtrc);
};



// ============================ 暂时未使用 =================================
//
//Image *edge_band(Image_yuv * rim, Image_f *lab_im, Image_f *mag, EdgeVec **edgeList, unsigned int edgesNum,
//                 float *edgeDiff, int *edgeDiffSize,
//                 const float lenthThred)
//{
//    return NULL;
//};