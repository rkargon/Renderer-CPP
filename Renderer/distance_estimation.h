//
//  fractal.hpp
//  Renderer
//
//  Created by Raphael Kargon on 1/5/17.
//  Copyright © 2017 Raphael Kargon. All rights reserved.
//

#ifndef distance_estimation_hpp
#define distance_estimation_hpp

#include <complex>
#include <functional>
#include <stdio.h>

#include "geom.h"
#include "raster.h"


// These functions return the min or max of two numbers, but with an adjustable amount of polynomial smoothing to make the transition smooth.
double smooth_min(double a, double b, double k);
double smooth_max(double a, double b, double k);

// Represents a 3D scene as a distance function - a function that, given a vertex, returns the distance from the vertex to the closest point in the scene.
using distance_estimator=std::function<double (const vertex&)>;


/* Distance Estimator Primitives */
distance_estimator de_plane(vertex normal);
distance_estimator de_sphere(vertex center, double radius);
//distance_estimator de_box(const bounds& bounding_box);
distance_estimator de_torus(double r1, double r2);
distance_estimator de_sierpinski_tetrahedron(int num_iterations);
distance_estimator de_menger(double scale, vertex center, int num_iterations);
distance_estimator de_mandelbulb(double power, int num_iterations, double bailout = 2);
// A fast implementation of the 8th-power mandelbulb
// NOTE: bailout here is for r^2
distance_estimator de_mandelbulb_8_fast(int num_iterations, double bailout = 4);
distance_estimator de_mandelbox(double scale, double folding_limit, double min_radius_2, double fixed_radius_2, int num_iterations);

/* Distance Estimator Operations */
distance_estimator de_union(distance_estimator de1, distance_estimator de2);
distance_estimator operator||(distance_estimator de1, distance_estimator de2);
// k determines how much blending occurs. 
distance_estimator de_blend(distance_estimator de1, distance_estimator de2, double k);
distance_estimator de_intersect(distance_estimator de1, distance_estimator de2);
distance_estimator operator&&(distance_estimator de1, distance_estimator de2);
distance_estimator de_intersect_blend(distance_estimator de1, distance_estimator de2, double k);
// Subtracts the second object from the first
distance_estimator de_subtract(distance_estimator de1, distance_estimator de2);
distance_estimator de_twist(distance_estimator de);


/**
 *  Estimates the normal of a distance field at a given point using a numerical gradient.
 *
 *  @param v        The vertex at which to find the normal
 *  @param obj      The given object as a distance field
 *  @param epsilon  Distance between individual points for sampling the gradient
 *
 *  @return A normalized vector representing the normal at this point.
 */
vertex estimate_normal(const vertex& v, const distance_estimator& obj, double epsilon = 0.00001);

/**
 *  Intersects a ray with a 3D scene given as a distance function using ray marching.
 *
 *  @param r             A ray to intersect with the scene
 *  @param obj           The 3D scene to render, as a distance function
 *  @param t             If the ray intersects, this stored the distance of the intersection along the ray from the ray's origin.
 *  @param steps         Stores number of steps ray was marched before intersecting.
 *  @param max_ray_steps Maximum number of steps to march the ray before assuming it does not intersect.
 *  @param .001          Tolerance for a ray to intersect with the scene.
 *
 *  @return Whether or not the ray intersects with the scene.
 */
bool ray_march(const ray& r, const distance_estimator& obj, double *t, int *steps, int max_ray_steps, double epsilon = 0.001);

    
#endif /* distance_estimation_hpp */
