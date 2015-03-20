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
#include "mesh.h"
#include "kdtree.h"
#include "scene.h"

#define AMB_OCC_SAMPLES 100
#define RAY_DEPTH 5

extern int num_rays_traced; //keeps track of how many rays have been traced, to measure performance

color calcLighting(const vertex& v, const vertex& n, const material& mat, scene* sc);
color traceRay(const ray& viewray, int depth, scene* sc);
color tracePath(const ray& viewray, int depth, scene* sc);
real ambientOcclusion(const ray& viewray, kdtree* kdt, int samples = AMB_OCC_SAMPLES);


#endif /* defined(__Renderer__rendering__) */
