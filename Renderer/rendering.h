//
//  rendering.h - the magic happens here!
//  Renderer
//
//  Functions for rendering 3D scenes.
//
//  Rasterization:  Projection of triangles onto the image plane, does not calculate shadows
//
//  Raytracing:     Sends ray from the camera, and checks for shadows / reflection by
//                  bouncing rays through the 3D scene. Does not do indirect lighting.
//
//  Path Tracing:   Sends rays off in random direction, and averages the results to sample all light
//                  falling onto an object. More photorealistic, much slower.
//
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
#define PATH_TRACE_SAMPLES 100
#define RAY_DEPTH 5

extern int num_rays_traced; //keeps track of how many rays have been traced, to measure performance

color calc_lighting(const vertex& v, const vertex& n, const material& mat, const scene* sc);
color trace_ray(const ray& viewray, scene* sc, int depth=1);
color trace_path(const ray& viewray, scene* sc, int depth=1);
double ambient_occlusion(const ray& viewray, kdtree* kdt);
color rayTraceDistanceField(const ray& viewray, scene *sc, int num_iterations, int depth = 1);

/* Rasterization */
void generate_maps_vector(int mapflags, raster *imgrasters, scene *sc);
void zBufferDraw(raster *imgrasters, scene *sc);
void paintNormalMap(raster *imgrasters, scene *sc);
void SSAO(raster *imgrasters, scene *sc);

/* Used in multithreaded rendering */
color ray_trace_pixel(double x, double y, int w, int h, scene *sc);
color path_trace_pixel(double x, double y, int w, int h, scene *sc);
color amb_occ_pixel(double x, double y, int w, int h, scene *sc);
color ray_march_pixel(double x, double y, int w, int h, scene *sc);

//stores increment values of barycentric coordinates for 4 pixel values at once
typedef struct EdgeVect{
    static const int stepX = 4;
    static const int stepY = 1;
    
    __v4si oneStepX;
    __v4si oneStepY;
    
    __v4si init(const point_2d<int>& v0, const point_2d<int>& v1, const point_2d<int>& origin);
} EdgeVect;

#endif /* defined(__Renderer__rendering__) */
