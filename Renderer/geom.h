//
//  geom.h
//  Renderer
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
typedef vertex color;
typedef struct bounds bounds;
typedef bounds ray;

template <typename T> struct point2D{
    T x,y;
    point2D(){}
    point2D(T a, T b) :x(a),y(b){}
};
template<typename T> T orient2D(point2D<T> a, point2D<T> b, point2D<T> c){
    return ((b.x - a.x) * (c.y - a.y) - (b.y - a.y) * (c.x - a.x));
}

class vertex
{
public:
    vertex();
    vertex(real a, real b, real c);
    
    inline real len() const {return _sqrt(x*x+y*y+z*z);}
    inline real lensqr() const {return (x*x+y*y+z*z);}
    void normalize();
    void set(const real a, const real b, const real c);
    vertex reflection(const vertex& n) const;
    vertex unitvect() const;
    
    void operator+= (const vertex& v2);
    void operator+= (const vertex *v2);
    void operator-= (const vertex& v2);
    void operator-= (const vertex *v2);
    void operator*= (const real r);
    void operator*= (const vertex& v2);
    void operator*= (const vertex *v2);
    vertex operator-();
    friend vertex operator+ (const vertex& v1, const vertex& v2);
    friend vertex operator+ (const vertex& v1, const vertex *v2);
    friend vertex operator- (const vertex& v1, const vertex& v2);
    friend vertex operator- (const vertex& v1, const vertex *v2);
    friend vertex operator* (const vertex& v1, real r);
    friend vertex operator* (const real r, vertex& v1);
    friend vertex operator* (const vertex& v1, const vertex& v2);
    friend bool operator== (const vertex& v1, const vertex& v2);
    friend bool operator!= (const vertex& v1, const vertex& v2);
    friend std::ostream& operator<<(std::ostream& os, const vertex& v);
    
    //data
    union{
        struct{real r, g, b;};
        struct{real t, u, v;};
        struct{real x, y, z;};
        struct{real vec[3];};
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
    real tex_u,tex_v;
    
    meshvertex();
    meshvertex(real x, real y, real z);
    vertex vertexNormal();
    void operator=(const meshvertex& v);
    friend bool operator==(const meshvertex& v1, const meshvertex& v2);
};

struct bounds{
    union{
        struct{vertex min, max;};
        struct{vertex org, dir;};
        struct{vertex varr[2];};
    };
    
    bounds();
    bounds(const vertex& a, const vertex& b);
    
    real area() const;
    inline real volume() const{return (max.x-min.x)*(max.y-min.y)*(max.z-min.z);}
    inline real d(const int k) const{return max.vec[k]-min.vec[k];} //doesn't check axis indices for validity. Be careful.
};

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
    inline real minCoord(int axis) const{
        return _min(_min(vertices[0]->vec[axis],vertices[1]->vec[axis]),vertices[2]->vec[axis]);
    }
    inline real maxCoord(int axis) const{
        return _max(_max(vertices[0]->vec[axis],vertices[1]->vec[axis]),vertices[2]->vec[axis]);
    }
};

class edge
{
public:
    vertex *v1, *v2;
    edge();
    edge(vertex *vert1, vertex *vert2);
    vertex getVector();
    real length();
    
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
//NOTE: be careful with arguments, as you should know, cross(A,B) != cross(B,A)
inline vertex cross(const vertex& A, const vertex& B){ return vertex(A.y*B.z - A.z*B.y, A.z*B.x - A.x*B.z, A.x*B.y - A.y*B.x); }
inline real dot(const vertex& A, const vertex& B){ return (A.x*B.x + A.y*B.y + A.z*B.z); }
vertex lerp(const vertex& v1, const vertex& v2, const real r);
vertex lerp(const vertex& v1, const vertex& v2, const vertex& v3, real w1, real w2, real w3);
vertex min3(const vertex& v1, const vertex& v2);
vertex max3(const vertex& v1, const vertex& v2);
vertex randomDirection(); //picks a random point on the unit sphere (uniformly)

//color functions
uint colorToRGB(const color& c);
uint normalToRGB(const vertex& n);
color RGBToColor(const uint rgb);
vertex RGBToNormal(const uint n);

//bounding box related functions
bounds calcBoundingBox(const std::vector<face*>& faces);
void intersectBoundingBoxes(const bounds& b1, const bounds& b2, bounds& newbounds);
void listBounds(bounds& newbounds, int nverts, const vertex* vertices ...);

//geometry intersections
bool rayAABBIntersect(const bounds& AABB, const ray& r);
face *rayFacesIntersect(const std::vector<face*>& faces, const ray& r, bool lazy, vertex *tuv);
bool raySphereIntersect(const ray& r, const real rad, real& t); //intersects a ray with a sphere centered on the origin. Assumes ray direction is normalized

#endif /* defined(__Renderer__geom__) */
