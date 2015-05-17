//
//  rendering.h - the magic happens here!
//  Renderer
//
//  Created by Raphael Kargon on 6/8/14.
//  Copyright (c) 2014 Raphael Kargon. All rights reserved.
//

#ifndef __Renderer__rendering__
#define __Renderer__rendering__

#include "camera.h"
#include "kdtree.h"
#include "mesh.h"
#include "raster.h"
#include "scene.h"

#define AMB_OCC_SAMPLES 100
#define RAY_DEPTH 5

extern int num_rays_traced; //keeps track of how many rays have been traced, to measure performance

color calcLighting(const vertex& v, const vertex& n, const material& mat, scene* sc);
color traceRay(const ray& viewray, int depth, scene* sc);
color tracePath(const ray& viewray, int depth, scene* sc);
double ambientOcclusion(const ray& viewray, kdtree* kdt, int samples = AMB_OCC_SAMPLES);


/* Rasterization */
void generate_maps_vector(int mapflags, raster *imgrasters, scene *sc);
void zBufferDraw_vector(raster *imgrasters, scene *sc);
void paintNormalMap(raster *imgrasters, scene *sc);
void SSAO(raster *imgrasters, scene *sc);
void rayTraceUnthreaded(raster *imgrasters, scene *sc, int tilesize);
void pathTraceUnthreaded(raster *imgrasters, scene *sc, int tilesize, int pathTracingSamples);

//stores increment values of barycentric coordinates for 4 pixel values at once
typedef struct EdgeVect{
    static const int stepX = 4;
    static const int stepY = 1;
    
    __v4si oneStepX;
    __v4si oneStepY;
    
    __v4si init(const point2D<int>& v0, const point2D<int>& v1, const point2D<int>& origin);
} EdgeVect;

#endif /* defined(__Renderer__rendering__) */
