//
//  geom.h
//  Renderer
//
// Functions and data structures relating to geometry. 
//
//  Created by Raphael Kargon on 6/5/14.
//  Copyright (c) 2014 Raphael Kargon. All rights reserved.
//

#ifndef __Renderer__geom__
#define __Renderer__geom__

#include <iostream>
#include <vector>
#include "common.h"

class vertex;
class meshvertex;
class face;
class edge;
class mesh;
class bounds;
typedef vertex color; //In the 'color' type, values range from 0 to 1 (not 0 to 255)
typedef bounds ray;

template <typename T> struct point_2d{
    T x,y;
    point_2d(){}
    point_2d(T a, T b) :x(a),y(b){}
};
//cross product of AB, and AC
template<typename T> T orient_2d(point_2d<T> a, point_2d<T> b, point_2d<T> c){
    return ((b.x - a.x) * (c.y - a.y) - (b.y - a.y) * (c.x - a.x));
}

//Represents a point in three-dimensional space in cartesian coordinates. Also used to represent color values.
class vertex
{
public:
    vertex();
    vertex(double a, double b, double c);
    
    inline double len() const {return sqrt(x*x+y*y+z*z);}
    inline double lensqr() const {return (x*x+y*y+z*z);}
    void normalize();
    void set(const double a, const double b, const double c);
    vertex reflection(const vertex& n) const;
    vertex unitvect() const; //return a normalized copy of this vector
    void to_polar(double &radius, double& theta, double& phi) const;
    // Convenience function for getting polar coordinates, given a pre-calculated radius. This saves some computation time.
    void to_polar_angles(const double radius_squared, double &theta, double &phi) const;
    vertex box_fold(const double l) const;
    // Returns the scaling factor for a sphere fold 
    double sphere_fold_ratio(const double min_radius_2, const double fixed_radius_2) const;
    
    void operator+= (const vertex& v2);
    void operator+= (const vertex *v2);
    void operator-= (const vertex& v2);
    void operator-= (const vertex *v2);
    void operator*= (const double r);
    void operator*= (const vertex& v2);
    void operator*= (const vertex *v2);
    vertex operator-();
    friend vertex operator+ (const vertex& v1, const vertex& v2);
    friend vertex operator+ (const vertex& v1, const vertex *v2);
    friend vertex operator- (const vertex& v1, const vertex& v2);
    friend vertex operator- (const vertex& v1, const vertex *v2);
    friend vertex operator* (const vertex& v1, double r);
    friend vertex operator* (const double r, const vertex& v1);
    friend vertex operator* (const vertex& v1, const vertex& v2);
    friend vertex operator/ (const double r, const vertex& v1);
    friend vertex operator/ (const vertex& v1, const double r);
    friend bool operator== (const vertex& v1, const vertex& v2);
    friend bool operator!= (const vertex& v1, const vertex& v2);
    friend std::ostream& operator<<(std::ostream& os, const vertex& v);
    
    //data
    union{
        struct{double r, g, b;};
        struct{double t, u, v;};
        struct{double x, y, z;};
        struct{double vec[3];};
        //ok so maybe I'm lazy
    };
};
struct vertex_ptr_hasher{
    std::size_t operator()(const vertex* v) const;
};
struct vertex_hasher{
    std::size_t operator()(const vertex v) const;
};

//stores a vertex, as well as texture coordinates and adjacent faces.
class meshvertex : public vertex
{
public:
    std::vector<face*> faces;
    double tex_u,tex_v;
    
    meshvertex();
    meshvertex(double x, double y, double z);
    vertex vertex_normal();
    void operator=(const meshvertex& v);
    friend bool operator==(const meshvertex& v1, const meshvertex& v2);
};

//Represents an axis-aligned bounding box. Also represents rays, with an origin and direction
class bounds{
public:
    union{
        struct{vertex min, max;};
        struct{vertex org, dir;};
        struct{vertex varr[2];};
    };
    
    bounds();
    bounds(const vertex& a, const vertex& b);
    
    double area() const;
    inline double volume() const{return (max.x-min.x)*(max.y-min.y)*(max.z-min.z);}
    inline double d(const int k) const{return max.vec[k]-min.vec[k];} //doesn't check axis indices for validity. Be careful.
};

//Represents a triangular face on a 3D mesh
class face
{
public:
    vertex normal;
    meshvertex *vertices[3];
    mesh *obj;
    
    face();
    face(const vertex& norm, meshvertex * const v1, meshvertex * const v2, meshvertex * const v3, mesh * const object);
    bounds bounding_box() const;
    vertex center() const;
    vertex generate_normal();
    bool intersect_ray_triangle(const ray& r, vertex *tuv);
    bool is_perpendicular(const vertex& normal);
    inline double min_coord(int axis) const{
        return fmin(fmin(vertices[0]->vec[axis],vertices[1]->vec[axis]),vertices[2]->vec[axis]);
    }
    inline double max_coord(int axis) const{
        return fmax(fmax(vertices[0]->vec[axis],vertices[1]->vec[axis]),vertices[2]->vec[axis]);
    }
};

class edge
{
public:
    vertex *v1, *v2;
    edge();
    edge(vertex *vert1, vertex *vert2);
    vertex get_vector();
    double length();
    
    void operator=(const edge& e);
    friend bool operator==(const edge& e1, const edge& e2);
};
struct edge_hasher{
    std::size_t operator()(const edge& e) const;
};
struct edge_ptr_hasher{
    std::size_t operator()(edge* const& e) const;
};
struct edge_ptr_equality{
    template<typename X, typename Y>
    bool operator()(X const &lhs, Y const &rhs) const;
};

/* vertex functions */
inline vertex cross(const vertex& A, const vertex& B){ return vertex(A.y*B.z - A.z*B.y, A.z*B.x - A.x*B.z, A.x*B.y - A.y*B.x); }
inline double dot(const vertex& A, const vertex& B){ return (A.x*B.x + A.y*B.y + A.z*B.z); }
vertex rotate(const vertex& v, const vertex& axis, const double dtheta);
vertex abs(const vertex& v);
vertex lerp(const vertex& v1, const vertex& v2, const double r);
//linearly interpolate between 3 vertices using barycentric coordinates
vertex lerp(const vertex& v1, const vertex& v2, const vertex& v3, double w1, double w2, double w3);
vertex lerp(const vertex& v1, const vertex& v2, const vertex& v3, const vertex& w);
//return a vertex containig the min/max of each coordinate from the given two inputs.
vertex min(const vertex& v1, const vertex& v2);
vertex max(const vertex& v1, const vertex& v2);
vertex clamp(const vertex& v, const vertex& min, const vertex& max);
vertex random_direction(); //picks a random point on the unit sphere (uniformly)
vertex from_polar(const double radius, const double theta, const double phi);

/* color functions */
uint color_to_rgb(const color& c);
uint normal_to_rgb(const vertex& n);
color rgb_to_color(const uint rgb);
vertex rgb_to_normal(const uint n);
color hsv_to_rgb(const int hue, const double saturation, const double value);

/* bounding box related functions */
bounds calc_bounding_box(const std::vector<face*>& faces);
void intersect_bounding_boxes(const bounds& b1, const bounds& b2, bounds& newbounds);
void list_bounds(bounds& newbounds, int nverts, const vertex* vertices ...);

/* geometry intersections */
bool ray_AABB_intersect(const bounds& AABB, const ray& r);
face *ray_faces_intersect(const std::vector<face*>& faces, const ray& r, bool lazy, vertex *tuv);
bool ray_sphere_intersect(const ray& r, const double rad, double& t); //intersects a ray with a sphere centered on the origin. Assumes ray direction is normalized

/* misc geometry functions */
point_2d<double> random_point_in_unit_circle();

#endif /* defined(__Renderer__geom__) */