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

template <typename T> struct point2D{
    T x,y;
    point2D(){}
    point2D(T a, T b) :x(a),y(b){}
};
//cross product of AB, and AC
template<typename T> T orient2D(point2D<T> a, point2D<T> b, point2D<T> c){
    return ((b.x - a.x) * (c.y - a.y) - (b.y - a.y) * (c.x - a.x));
}

//Represents a point in three-dimensional space. Also used to represent color values.
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
struct vertexPtrHasher{
    std::size_t operator()(const vertex* v) const;
};
struct vertexHasher{
    std::size_t operator()(const vertex v) const;
};
struct ptrEquality{
    template< typename X, typename Y >
    bool operator() ( X const &lhs, Y const &rhs ) const
    { return * lhs == * rhs; }
};

//stores a vertex, as well as texture coordinates and adjacent faces.
class meshvertex : public vertex
{
public:
    std::vector<face*> faces;
    double tex_u,tex_v;
    
    meshvertex();
    meshvertex(double x, double y, double z);
    vertex vertexNormal();
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
    bounds boundingBox() const;
    vertex center() const;
    vertex generateNormal();
    bool intersectRayTriangle(const ray& r, vertex *tuv);
    bool isPerpendicular(const vertex& normal);
    inline double minCoord(int axis) const{
        return fmin(fmin(vertices[0]->vec[axis],vertices[1]->vec[axis]),vertices[2]->vec[axis]);
    }
    inline double maxCoord(int axis) const{
        return fmax(fmax(vertices[0]->vec[axis],vertices[1]->vec[axis]),vertices[2]->vec[axis]);
    }
};

class edge
{
public:
    vertex *v1, *v2;
    edge();
    edge(vertex *vert1, vertex *vert2);
    vertex getVector();
    double length();
    
    void operator=(const edge& e);
    friend bool operator==(const edge& e1, const edge& e2);
};
struct edgeHasher{
    std::size_t operator()(const edge& e) const;
};
struct edgePtrHasher{
    std::size_t operator()(edge* const& e) const;
};
struct edgePtrEquality{
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
vertex randomDirection(); //picks a random point on the unit sphere (uniformly)

/* color functions */
uint colorToRGB(const color& c);
uint normalToRGB(const vertex& n);
color RGBToColor(const uint rgb);
vertex RGBToNormal(const uint n);

/* bounding box related functions */
bounds calcBoundingBox(const std::vector<face*>& faces);
void intersectBoundingBoxes(const bounds& b1, const bounds& b2, bounds& newbounds);
void listBounds(bounds& newbounds, int nverts, const vertex* vertices ...);

/* geometry intersections */
bool rayAABBIntersect(const bounds& AABB, const ray& r);
face *rayFacesIntersect(const std::vector<face*>& faces, const ray& r, bool lazy, vertex *tuv);
bool raySphereIntersect(const ray& r, const double rad, double& t); //intersects a ray with a sphere centered on the origin. Assumes ray direction is normalized

#endif /* defined(__Renderer__geom__) */
