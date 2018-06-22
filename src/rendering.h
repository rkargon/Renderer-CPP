//
//  rendering.h - the magic happens here!
//  Renderer
//
//  Functions for rendering 3D scenes.
//
//  Rasterization:  Projection of triangles onto the image plane, does not
//  calculate shadows
//
//  Raytracing:     Sends ray from the camera, and checks for shadows /
//  reflection by
//                  bouncing rays through the 3D scene. Does not do indirect
//                  lighting.
//
//  Path Tracing:   Sends rays off in random direction, and averages the results
//  to sample all light
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

#include <string>

extern int num_rays_traced; // keeps track of how many rays have been traced, to
                            // measure performance

struct render_options {
  unsigned int threads;
  unsigned int samples;
  unsigned int ray_depth;
  unsigned int tilesize;
  unsigned int antialiasing_per_side;
};

color calc_lighting(const vertex &v, const vertex &n, const material &mat,
                    const scene &sc);
color trace_ray(const ray &viewray, const scene &sc, unsigned int depth = 1,
                unsigned int max_depth = 10);
color trace_path(const ray &viewray, const scene &sc, unsigned int depth = 1,
                 unsigned int max_depth = 10);
double ambient_occlusion(const ray &viewray, kdtree *kdt, int samples);
color raytrace_distance_field(const ray &viewray, const scene &sc,
                              int num_iterations, int depth = 1);

/* Rasterization */
void generate_maps_vector(int mapflags, raster &imgrasters, const scene &sc);
void zbuffer_draw(raster &imgrasters, const scene &sc);
void paint_normal_map(raster &imgrasters, const scene &sc);
void SSAO(raster &imgrasters, const scene &sc);

/* Used in multithreaded rendering */
color ray_trace_pixel(double x, double y, int w, int h, const scene &sc,
                      const render_options &opts);
color path_trace_pixel(double x, double y, int w, int h, const scene &sc,
                       const render_options &opts);
color amb_occ_pixel(double x, double y, int w, int h, const scene &sc,
                    const render_options &opts);
color ray_march_pixel(double x, double y, int w, int h, const scene &sc,
                      const render_options &opts);

// stores increment values of barycentric coordinates for 4 pixel values at once
typedef struct edge_vect {
  static const int step_x = 4;
  static const int step_y = 1;

  __v4si one_step_x;
  __v4si one_step_y;

  __v4si init(const point_2d<int> &v0, const point_2d<int> &v1,
              const point_2d<int> &origin);
} edge_vect;

#endif /* defined(__Renderer__rendering__) */
