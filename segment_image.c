/*!
* segment-image.h
* date: 2017/03/21 15:03
*
* author: dav1sNJU
*
* brief:
*
*
*/
#include <math.h>
#include "image.h"
#include "segment_image.h"
#include <assert.h>

// dissimilarity measure between pixels
static float diff(Image_f *r, Image_f *g, Image_f *b,
                  int x1, int y1, int x2, int y2)
{
    return (float)sqrt(fsquare(imRef(r, x1, y1) - imRef(r, x2, y2)) +
                       fsquare(imRef(g, x1, y1) - imRef(g, x2, y2)) +
                       fsquare(imRef(b, x1, y1) - imRef(b, x2, y2)));
};

static void getContextYUV(Context *cntxt, Image_f *Py, Image_f *Pcr, Image_f *Pcb,
                            int width, int height)
{
    int roffset = 0;
    for (int y = 0; y < height; y++)
    {
        for (int x = 0; x < width; x++)
        {
            int pos = x + roffset;
            imRef(Py, x, y) = cntxt->pic[pos].Y;
            imRef(Pcr, x, y) = cntxt->pic[pos].Cr;
            imRef(Pcb, x, y) = cntxt->pic[pos].Cb;
        }
        roffset += width;
    }
};

static void buildGraph(Image_f *Py, Image_f *Pcr, Image_f *Pcb, 
                    Edge *edges, int width, int height, int *n)
{
    int num = 0;
    for (int y = 0; y < height; y++)
        for (int x = 0; x < width; x++)
        {
            if (x < width - 1)
            {
                edges[num].a = y * width + x;
                edges[num].b = y * width + (x + 1);
                edges[num].w = diff(Py, Pcr, Pcb, x, y, x + 1, y);// 计算权重
                num++;
            }

            if (y < height - 1)
            {
                edges[num].a = y * width + x;
                edges[num].b = (y + 1) * width + x;
                edges[num].w = diff(Py, Pcr, Pcb, x, y, x, y + 1);
                num++;
            }

            if ((x < width - 1) && (y < height - 1))
            {
                edges[num].a = y * width + x;
                edges[num].b = (y + 1) * width + (x + 1);
                edges[num].w = diff(Py, Pcr, Pcb, x, y, x + 1, y + 1);
                num++;
            }

            if ((x < width - 1) && (y > 0))
            {
                edges[num].a = y * width + x;
                edges[num].b = (y - 1) * width + (x + 1);
                edges[num].w = diff(Py, Pcr, Pcb, x, y, x + 1, y - 1);
                num++;
            }
        }
    *n = num;
};

Universe * segment_image_cntxt(Context *cntxt, float sigma, float c, int min_size,
                             Image_yuv *dst, Image_f *labelImg, int *num_ccs)
{
    int width = cntxt->pic_col, height = cntxt->pic_row;
    Image_f *Py = create_new_imagef(width, height, 0);
    Image_f *Pcr = create_new_imagef(width, height, 0);
    Image_f *Pcb = create_new_imagef(width, height, 0);    
    getContextYUV(cntxt, Py, Pcr, Pcb, width, height);

    smoothWrapper(Py, Py, sigma);
    smoothWrapper(Pcr, Pcr, sigma);
    smoothWrapper(Pcb, Pcb, sigma);

    // build graph
    Edge *edges = (Edge *)malloc(sizeof(Edge) * width * height * 4);
    if (!edges) return NULL;
    int num = 0;//edges num    
    buildGraph(Py, Pcr, Pcb, edges, width, height, &num);

    destroy_imagef(Py);
    destroy_imagef(Pcr);
    destroy_imagef(Pcb);

    // make a disjoint-set forest   
    int max_vertices_num = width * height;
    Universe *u = init_universe(max_vertices_num);

    segment_graph(u, max_vertices_num, num, edges, c);

    // post process small components
    for (int i = 0; i < num; i++)
    {
        Edge *piEdge = &edges[i];
        int a = find_in_universe(u, piEdge->a); //找到发现的标号类别。
        int b = find_in_universe(u, piEdge->b);
        if ((a != b) && ((get_universe_size(u, a) < min_size) || (get_universe_size(u, b) < min_size)))
            join_universe(u, a, b);
    }
    free(edges);
    *num_ccs = u->num;
    for (int y = 0; y < height; y++)
    {
        int roffset = y * width;
        for (int x = 0; x < width; x++)
        {
            int comp = find_in_universe(u, x + roffset);  //分析属于哪个类别;
            imRef(labelImg, x, y) = (float)comp;
        }
    }
    return u;
};


void gradient(Image_f * imLabel, Image *segOut, float thd) 
{
    //G(:, 1) = A(:, 2) - A(:, 1);
    //G(:, j) = 0.5*(A(:, j + 1) - A(:, j - 1));
    //G(:, N) = A(:, N) - A(:, N - 1);
    int w = imLabel->w, h = imLabel->h;
    float fx = 0, fy = 0;
    for (int y = 0; y < h; ++y)
    {
        for (int x = 0; x < w; ++x)
        {
            fx = 0.0f;
            fy = 0.0f;
            // dx
            if (x == 0)
                fx = (float)(imRef(imLabel, 1, y) - imRef(imLabel, 0, y));
            else if (x == w - 1)
                fx = (float)(imRef(imLabel, x, y) - imRef(imLabel, x - 1, y));
            else
                fx = 0.5f * (imRef(imLabel, x + 1, y) - imRef(imLabel, x - 1, y));

            // dy
            if (y == 0)
                fy = (float)(imRef(imLabel, x, 1) - imRef(imLabel, x, 0));
            else if (y == h - 1)
                fy = (float)(imRef(imLabel, x, y) - imRef(imLabel, x, y - 1));
            else
                fy = 0.5f * (imRef(imLabel, x, y + 1) - imRef(imLabel, x, y - 1));

            /* 提取边界线 */
            if ((float)sqrt(fsquare(fx) + fsquare(fy)) > thd)
            {
                imRef(segOut, x, y) = 1u;
            }
        }
    }
};

// random color
rgb random_rgb()
{
    rgb c;
    //double r;
    c.r = (uint8_t)rand();
    c.g = (uint8_t)rand();
    c.b = (uint8_t)rand();
    return c;
};
yuv random_yuv()
{
    yuv c;
    //double r;
    c.y = (uint8_t)rand();
    c.u = (uint8_t)rand();
    c.v = (uint8_t)rand();
    return c;
};

// ============================ 暂时未使用 =================================
/*
static void getYUVFxFy(Image_yuv * im, Image_f * fxPy, Image_f * fxPcr, Image_f * fxPcb,
    Image_f * fyPy, Image_f * fyPcr, Image_f * fyPcb)
{
    int w = im->w, h = im->h;
    for (int x = 0; x < w; x++)
    {
        for (int y = 0; y < h; y++)
        {
            //fx
            if (x == 0)
            {
                imRef(fxPy, x, y) = (float)(imRef(im, 1, y).y - imRef(im, 0, y).y);
                imRef(fxPcr, x, y) = (float)(imRef(im, 1, y).u - imRef(im, 0, y).u);
                imRef(fxPcb, x, y) = (float)(imRef(im, 1, y).v - imRef(im, 0, y).v);
            }
            else if (x == w - 1)
            {
                imRef(fxPy, x, y) = (float)(imRef(im, x, y).y - imRef(im, x - 1, y).y);
                imRef(fxPcr, x, y) = (float)(imRef(im, x, y).u - imRef(im, x - 1, y).u);
                imRef(fxPcb, x, y) = (float)(imRef(im, x, y).v - imRef(im, x - 1, y).v);
            }
            else
            {
                imRef(fxPy, x, y) = 0.5f *((float)(imRef(im, x + 1, y).y - imRef(im, x - 1, y).y));
                imRef(fxPcr, x, y) = 0.5f *((float)(imRef(im, x + 1, y).u - imRef(im, x - 1, y).u));
                imRef(fxPcb, x, y) = 0.5f *((float)(imRef(im, x + 1, y).v - imRef(im, x - 1, y).v));
            }
            //fy
            if (y == 0)
            {
                imRef(fyPy, x, y) = (float)(imRef(im, x, 1).y - imRef(im, x, y).y);
                imRef(fyPcr, x, y) = (float)(imRef(im, x, 1).u - imRef(im, x, y).u);
                imRef(fyPcb, x, y) = (float)(imRef(im, x, 1).v - imRef(im, x, y).v);
            }
            else if (y == h - 1)
            {
                imRef(fyPy, x, y) = (float)(imRef(im, x, y).y - imRef(im, x, y - 1).y);
                imRef(fyPcr, x, y) = (float)(imRef(im, x, y).u - imRef(im, x, y - 1).u);
                imRef(fyPcb, x, y) = (float)(imRef(im, x, y).v - imRef(im, x, y - 1).v);
            }
            else
            {
                imRef(fyPy, x, y) = 0.5f *((float)(imRef(im, x, y + 1).y - imRef(im, x, y - 1).y));
                imRef(fyPcr, x, y) = 0.5f *((float)(imRef(im, x, y + 1).u - imRef(im, x, y - 1).u));
                imRef(fyPcb, x, y) = 0.5f *((float)(imRef(im, x, y + 1).v - imRef(im, x, y - 1).v));
            }
        }
    }
};

void gradientYUV(Image_yuv * im, Image_f * ed1, Image_f * mag)
{
    assert(im);
    assert(ed1);
    assert(mag);
    int w = im->w, h = im->h;
    Image_f *fxPy = create_new_imagef(w, h, 0);
    Image_f *fxPcr = create_new_imagef(w, h, 0);
    Image_f *fxPcb = create_new_imagef(w, h, 0);
    Image_f *fyPy = create_new_imagef(w, h, 0);
    Image_f *fyPcr = create_new_imagef(w, h, 0);
    Image_f *fyPcb = create_new_imagef(w, h, 0);

    getYUVFxFy(im, fxPy, fxPcr, fxPcb, fyPy, fyPcr, fyPcb);

    float ed1err = 0.0f;
    float maxed1 = 1 - .0f;
    for (int x = 0; x < w; x++)
        for (int y = 0; y < h; y++)
        {
            float ed1fx = (float)sqrt(fsquare(imRef(fxPy, x, y)) + fsquare(imRef(fxPcr, x, y)) + fsquare(imRef(fxPcb, x, y)));
            float ed1fy = (float)sqrt(fsquare(imRef(fyPy, x, y)) + fsquare(imRef(fyPcr, x, y)) + fsquare(imRef(fyPcb, x, y)));
            float ed1xy = (float)sqrt(fsquare(ed1fx) + fsquare(ed1fy));
            imRef(ed1, x, y) = ed1xy;
            if (ed1xy > maxed1)maxed1 = ed1xy;
        }

    for (int x = 0; x < w; x++)
        for (int y = 0; y < h; y++)
            imRef(mag, x, y) = imRef(ed1, x, y) / maxed1;

    destroy_imagef(fxPy);
    destroy_imagef(fxPcr);
    destroy_imagef(fxPcb);
    free(fxPy);
    free(fxPcr);
    free(fxPcb);
    destroy_imagef(fyPy);
    destroy_imagef(fyPcr);
    destroy_imagef(fyPcb);
    free(fyPy);
    free(fyPcr);
    free(fyPcb);
};

Universe * segment_image_rgb(Image_rgb *im, float sigma, float c, int min_size,
                             Image_rgb *dst, Image_f *labelImg, int *num_ccs)
{
    int width = im->w;
    int height = im->h;
    Image_f *Pr = create_new_imagef(width, height, 0);
    Image_f *Pg = create_new_imagef(width, height, 0);
    Image_f *Pb = create_new_imagef(width, height, 0);
    
    // smooth each color channel  
    for (int y = 0; y < height; y++)
        for (int x = 0; x < width; x++)
        {
            imRef(Pr, x, y) = imRef(im, x, y).r;
            imRef(Pg, x, y) = imRef(im, x, y).g;
            imRef(Pb, x, y) = imRef(im, x, y).b;
        }

    smoothWrapper(Pr, Pr, sigma);
    smoothWrapper(Pg, Pg, sigma);
    smoothWrapper(Pb, Pb, sigma);

    // build graph
    Edge *edges = (Edge *)malloc(sizeof(Edge) * width * height * 4);
    if (!edges) return NULL;
    int num = 0;
    // 在这里定义几个neighbor 之间的关系，共4个邻域。
    for (int y = 0; y < height; y++)
        for (int x = 0; x < width; x++)
        {
            if (x < width - 1)
            {
                edges[num].a = y * width + x;
                edges[num].b = y * width + (x + 1);
                edges[num].w = diff(Pr, Pg, Pb, x, y, x + 1, y);// 计算权重
                num++;
            }

            if (y < height - 1)
            {
                edges[num].a = y * width + x;
                edges[num].b = (y + 1) * width + x;
                edges[num].w = diff(Pr, Pg, Pb, x, y, x, y + 1);
                num++;
            }

            if ((x < width - 1) && (y < height - 1))
            {
                edges[num].a = y * width + x;
                edges[num].b = (y + 1) * width + (x + 1);
                edges[num].w = diff(Pr, Pg, Pb, x, y, x + 1, y + 1);
                num++;
            }

            if ((x < width - 1) && (y > 0))
            {
                edges[num].a = y * width + x;
                edges[num].b = (y - 1) * width + (x + 1);
                edges[num].w = diff(Pr, Pg, Pb, x, y, x + 1, y - 1);
                num++;
            }
        }

    destroy_imagef(Pr); free(Pr);
    destroy_imagef(Pg); free(Pg);
    destroy_imagef(Pb); free(Pb);

    // segment
    // make a disjoint-set forest    ====================
    int max_vertices_num = width * height;
    Universe *u = init_universe(max_vertices_num);

    segment_graph(u, max_vertices_num, num, edges, c);

    // post process small components
    for (int i = 0; i < num; i++)
    {
        int a = find_in_universe(u, edges[i].a); //找到发现的标号类别。
        int b = find_in_universe(u, edges[i].b);

        if ((a != b)
            && ((get_universe_size(u, a) < min_size) || (get_universe_size(u, b) < min_size)))
            join_universe(u, a, b);
    }

    free(edges);
    *num_ccs = u->num;

    //根据得到的轮廓线
    // pick random colors for each component  
    rgb *colors = (rgb *)malloc(sizeof(rgb) * width * height);
    for (int i = 0; i < width * height; i++)
        colors[i] = random_rgb();

    for (int y = 0; y < height; y++)
        for (int x = 0; x < width; x++)
        {
            int comp = find_in_universe(u, y * width + x);  //分析属于哪个类别;
            imRef(labelImg, x, y) = (float)comp;
            //((rgb *)(data + j*w + i))->b
            imRef(dst, x, y).r = colors[comp].r; //对每个类别分配颜色。
            imRef(dst, x, y).g = colors[comp].g; //对每个类别分配颜色。
            imRef(dst, x, y).b = colors[comp].b; //对每个类别分配颜色。
        }

    free(colors);
    return u;
};
*/