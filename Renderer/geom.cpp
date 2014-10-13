//
//  geom.cpp
//  Renderer
//
//  Created by Raphael Kargon on 6/5/14.
//  Copyright (c) 2014 Raphael Kargon. All rights reserved.
//

#include "geom.h"

/* VERTEX */
vertex::vertex() :x(0),y(0),z(0){};
vertex::vertex(real a, real b, real c) :x(a),y(b),z(c) {};

//Scales the vector such that its magnitude is one, while keeping its orientation the same.
//If the vector is 0, it is unchanged.
void vertex::normalize(){
    real len = x * x + y * y + z * z;
    if (len != 0 && len != 1) {
        len = 1.0/_sqrt(len);
        x *= len;
        y *= len;
        z *= len;
    }
}
void vertex::set(const real a, const real b, const real c){x=a;y=b;z=c;}

//Calculates the refelction of this vertex across a normal vector n.
//REFL = V - N*(2(V.N))
vertex vertex::reflection(const vertex& n) const {return *this - n*(2*dot(*this,n));}

//Returns the corresponding unit vector of this vertex, but this vertex is not modified.
vertex vertex::unitvect() const{
    real len = x*x+y*y+z*z;
    if(len!=0 && len!=1){
        len = 1.0/_sqrt(len);
        return vertex(x*len, y*len, z*len);
    }
    return vertex(x,y,z);
}

void vertex::operator+= (const vertex& v2){x+=v2.x; y+=v2.y; z+=v2.z;}
void vertex::operator+= (const vertex *v2){x+=v2->x; y+=v2->y; z+=v2->z;}
void vertex::operator-= (const vertex& v2){x-=v2.x; y-=v2.y; z-=v2.z;}
void vertex::operator-= (const vertex *v2){x-=v2->x; y-=v2->y; z-=v2->z;}
void vertex::operator*= (const real r){x*=r; y*=r; z*=r;} //scalar multiplication
void vertex::operator*= (const vertex& v2){x*=v2.x; y*=v2.y; z*=v2.z;} //coordinate-wise multiplication
void vertex::operator*= (const vertex *v2){x*=v2->x; y*=v2->y; z*=v2->z;}
vertex vertex::operator-() {return vertex(-x, -y, -z);}

vertex operator+ (const vertex& v1, const vertex& v2) {return vertex(v1.x+v2.x, v1.y+v2.y, v1.z+v2.z);}
vertex operator+ (const vertex& v1, const vertex *v2) {return vertex(v1.x+v2->x, v1.y+v2->y, v1.z+v2->z);}
vertex operator- (const vertex& v1, const vertex& v2) {return vertex(v1.x-v2.x, v1.y-v2.y, v1.z-v2.z);}
vertex operator- (const vertex& v1, const vertex *v2) {return vertex(v1.x-v2->x, v1.y-v2->y, v1.z-v2->z);}
vertex operator* (const vertex& v1, real r){return vertex(v1.x*r, v1.y*r, v1.z*r);}
vertex operator* (const real r, vertex& v1){return vertex(v1.x*r, v1.y*r, v1.z*r);}
vertex operator* (const vertex& v1, const vertex& v2) {return vertex(v1.x*v2.x, v1.y*v2.y, v1.z*v2.z);} //coordinate-wise multiplication
bool operator== (const vertex& v1, const vertex& v2){return v1.x==v2.x && v1.y==v2.y && v1.z==v2.z;}
bool operator!= (const vertex& v1, const vertex& v2){return !(v1==v2);}
std::ostream& operator<< (std::ostream &os, const vertex &v){return os << "(" << v.x << ", " << v.y << ", " << v.z << ")";}
std::ostream& operator<< (std::ostream &os, const vertex *&v){return os << "(" << v->x << ", " << v->y << ", " << v->z << ")";}

//two vertex*'s or two vertex's are equal in these hash functions if they represent the same point in space.
std::size_t vertexPtrHasher::operator()(const vertex* v) const{
    size_t result = 1;
    unsigned realsize bits = *(realsize *)&v->x;
    result = 31 * result + (size_t) (bits ^ (bits >> realbytes*4));
    bits  = *(realsize *)&v->y;
    result = 31 * result + (size_t) (bits ^ (bits >> realbytes*4));
    bits = *(realsize *)&v->z;
    result = 31 * result + (size_t) (bits ^ (bits >> realbytes*4));
    return result;
}

std::size_t vertexHasher::operator()(const vertex v) const{
    size_t result = 1;
    unsigned realsize bits = *(realsize *)&v.x;
    result = 31 * result + (size_t) (bits ^ (bits >> realbytes*4));
    bits  = *(realsize *)&v.y;
    result = 31 * result + (size_t) (bits ^ (bits >> realbytes*4));
    bits = *(realsize *)&v.z;
    result = 31 * result + (size_t) (bits ^ (bits >> realbytes*4));
    return result;
}

/* FACE */
face::face() = default;
face::face(const vertex& norm, meshvertex * const v1, meshvertex * const v2, meshvertex * const v3, mesh * const object){
    this->obj = object;
    this->normal = norm;
    vertices[0] = v1;
    vertices[1] = v2;
    vertices[2] = v3;
    
    if(normal.lensqr() && isPerpendicular(normal)){
        normal.normalize();
    }
    else{
        normal = generateNormal(); //if normal is invalid, generate it
    }
}
bool face::isPerpendicular(const vertex &normal){
    vertex v12 = *vertices[1] - *vertices[0];
    vertex v23 = *vertices[2] - *vertices[1];
    //normal must be perpendicular to two of the edges, and must have nonzero length
    return (normal.lensqr() && !dot(v12, normal) && !dot(v23, normal));
}
//use cross product of edges to create normal. Orientation depends on order of vertices
vertex face::generateNormal(){
    vertex v12 = *vertices[1] - *vertices[0];
    vertex v23 = *vertices[2] - *vertices[1];
    vertex n = cross(v12, v23);
    n.normalize();
    return n;
}
//the center of the face, the arithmetic mean of vertices.
vertex face::center() const{
    return (*vertices[0]+*vertices[1]+*vertices[2])*(1/3.0);
}

//tuv stores:
//t: distance from origin to intersect point
//u, v: barycentric coordinates of intersect based on edge1, edge2
//Muller-Trumbore algorithm
bool face::intersectRayTriangle(const ray& r, vertex *tuv){
    vertex edge1, edge2, tvec, pvec, qvec;
    real det, inv_det;
    real t, u, v;
    
    edge1 = (*vertices[1]-*vertices[0]);
    edge2 = (*vertices[2]-*vertices[0]);
    
    pvec = cross(r.dir, edge2);
    det = dot(edge1, pvec);
    
    tvec = r.org - *vertices[0];
    inv_det = 1.0/det;
    qvec = cross(tvec, edge1);
    
    if(det > EPSILON){
        u = dot(tvec, pvec);
        if(u<0.0 || u>det) return false;
        v = dot(r.dir, qvec);
        if(v<0.0 || u + v > det) return false;
    }
    else if (det<-EPSILON){
        u = dot(tvec, pvec);
        if(u>0.0 || u<det) return false;
        v = dot(r.dir, qvec);
        if(v>0.0 || u + v < det) return false;
        
    }
    else return false;
    
    u*=inv_det;
    v*=inv_det;
    t = dot(edge2, qvec) * inv_det;
    if(tuv!=nullptr) tuv->set(t, u, v);
    return (t>EPSILON);
}

bounds face::boundingBox() const{
    bounds b;
    b.min = vertex(minCoord(0), minCoord(1), minCoord(2));
    b.max = vertex(maxCoord(0), maxCoord(1), maxCoord(2));
    return b;
}


bounds::bounds() :min(0,0,0),max(0,0,0){}
bounds::bounds(const vertex& a, const vertex& b) :min(a),max(b){}

//surface area
real bounds::area() const{
    real dx = _abs(max.x - min.x);
    real dy = _abs(max.y - min.y);
    real dz = _abs(max.z - min.z);
    return 2*(dx*dy + dy*dz + dz*dx);
}

/* MESHVERTEX */
meshvertex::meshvertex(){
    x=0;
    y=0;
    z=0;
}
meshvertex::meshvertex(real a, real b, real c){
    x=a; y=b; z=c;
}
//sum the normal of each adjacent face, and normalize.
vertex meshvertex::vertexNormal(){
    vertex v=vertex();
    for(face *f : faces) v+=f->normal;
    v.normalize();
    return v;
};

void meshvertex::operator=(const meshvertex& v){
    x=v.x;
    y=v.y;
    z=v.z;
    faces = v.faces;
    tex_u = v.tex_u;
    tex_v = v.tex_v;
}

bool operator==(const meshvertex& v1, const meshvertex& v2){return v1.x==v2.x && v1.y==v2.y && v1.z==v2.z;}

/* EDGE */
edge::edge(){}
edge::edge(vertex *vert1, vertex *vert2)  :v1(vert1), v2(vert2) {}
vertex edge::getVector() {return (*v2) - (*v1);}
real edge::length() {return (*v2-*v1).len();}
void edge::operator=(const edge& e){
    v1 = e.v1;
    v2 = e.v2;
}
bool operator==(const edge& e1, const edge& e2){
    if(*e1.v1 == *e2.v1) return (*e1.v2 == *e2.v2);
    else if(*e1.v1 == *e2.v2) return (*e1.v2 == *e2.v1);
    else return false;
}

//Two edges are equal in this hash operation if their vertices represent the same points in space
// (but in any order, edge(v1,v2) and edge(v2, v1) have same hash
std::size_t edgeHasher::operator()(const edge &e) const{
    size_t result1 = 1;
    unsigned realsize bits = *(realsize *)&e.v1->x;
    result1 = 31 * result1 + (size_t) (bits ^ (bits >> realbytes*4));
    bits = *(realsize *)&e.v1->y;
    result1 = 31 * result1 + (size_t) (bits ^ (bits >> realbytes*4));
    bits = *(realsize *)&e.v1->z;
    result1 = 31 * result1 + (size_t) (bits ^ (bits >> realbytes*4));
    size_t result2 = 1;
    bits = *(realsize *)&e.v2->x;
    result2 = 31 * result2 + (size_t) (bits ^ (bits >> realbytes*4));
    bits = *(realsize *)&e.v2->y;
    result2 = 31 * result2 + (size_t) (bits ^ (bits >> realbytes*4));
    bits = *(realsize *)&e.v2->z;
    result2 = 31 * result2 + (size_t) (bits ^ (bits >> realbytes*4));
    size_t result;
    if(result1>result2) result = 31 * (31+result1) +result2;
    else result = 31 * (31+result2) + result1;
    return result;
}

std::size_t edgePtrHasher::operator()(edge* const& e) const{
    size_t result1 = 1;
    unsigned realsize bits = *(realsize *)&e->v1->x;
    result1 = 31 * result1 + (size_t) (bits ^ (bits >> realbytes*4));
    bits = *(realsize *)&e->v1->y;
    result1 = 31 * result1 + (size_t) (bits ^ (bits >> realbytes*4));
    bits = *(realsize *)&e->v1->z;
    result1 = 31 * result1 + (size_t) (bits ^ (bits >> realbytes*4));
    size_t result2 = 1;
    bits = *(realsize *)&e->v2->x;
    result2 = 31 * result2 + (size_t) (bits ^ (bits >> realbytes*4));
    bits = *(realsize *)&e->v2->y;
    result2 = 31 * result2 + (size_t) (bits ^ (bits >> realbytes*4));
    bits = *(realsize *)&e->v2->z;
    result2 = 31 * result2 + (size_t) (bits ^ (bits >> realbytes*4));
    size_t result;
    if(result1>result2) result = 31 * (31+result1) +result2;
    else result = 31 * (31+result2) + result1;
    return result;
}
template< typename X, typename Y >
bool edgePtrEquality::operator()(X const &lhs, Y const &rhs) const{
    return *lhs==*rhs;
}

//linearly interpolate two 3d vertices
vertex lerp(const vertex& v1, const vertex& v2, const real r){return vertex(v1.x+r*(v2.x-v1.x), v1.y+r*(v2.y-v1.y), v1.z+r*(v2.z-v1.z));}
//linearly interpolate 3 vertices based on barycentric coordinates.
vertex lerp(const vertex& v1, const vertex& v2, const vertex& v3, real w1, real w2, real w3){return v1*w1 + v2*w2 + v3*w3;}
//The maximum of each coordainte of the two vertices
vertex max3(const vertex& v1, const vertex& v2){return vertex(_max(v1.x,v2.x), _max(v1.y,v2.y), _max(v1.z,v2.z));}
vertex min3(const vertex& v1, const vertex& v2){return vertex(_min(v1.x,v2.x), _min(v1.y,v2.y), _min(v1.z,v2.z));}
//A random unit vector.
vertex randomDirection(){
    real theta = (real)rand()/RAND_MAX;
    theta *= M_PI*2;
    real z = (real)rand()/RAND_MAX;
    z = z*2 - 1;
    real w = _sqrt(1-z*z);
    return vertex(w*_cos(theta), w*_sin(theta), z);
}

//color functions
//Take a vertex of three values in range [0,1] and convert to 32-bit RGB
uint colorToRGB(const color& c){
    return (unsigned char)clamp(c.r*255, 0, 255)<<16 | (unsigned char)clamp(c.g*255, 0, 255)<<8 | (unsigned char)clamp(c.b*255, 0, 255);
}
//Convert a normal vector to 32-bit RGB. FOr each coordinate, the range [-1, 1] is mapped to [0, 255]
uint normalToRGB(const vertex& n){
    return (unsigned char)clamp(n.x*128+128, 0, 255)<<16 | (unsigned char)clamp(n.y*128+128, 0, 255)<<8 | (unsigned char)clamp(n.z*128+128, 0, 255);
}
//Convert 32-bit RGB to vertex
color RGBToColor(const uint rgb){
    return color(((rgb>>16)&0xff)/255.0,
                 ((rgb>>8)&0xff)/255.0,
                 (rgb&0xff)/255.0);
}
//Convert 32-bit RGB to normal vector
color RGBToNormal(const uint n){
    return vertex(((n>>16)&0xff)/255.0 - 0.5,
                  ((n>>8)&0xff)/255.0 - 0.5,
                  (n&0xff)/255.0 - 0.5);
}

//bounding box related functions

//THe bounding box of a set of faces.
bounds calcBoundingBox(const std::vector<face *>& faces){
    real minx, miny, minz, maxx, maxy, maxz, tmp;
    minx = miny = minz = maxx = maxy = maxz = _nan("");
    for(face *f: faces){
        tmp = f->minCoord(0);
        if(isnan(minx) || minx>tmp) minx = tmp;
        tmp = f->minCoord(1);
        if(isnan(miny) || miny>tmp) miny = tmp;
        tmp = f->minCoord(2);
        if(isnan(minz) || minz>tmp) minz = tmp;
        
        tmp = f->maxCoord(0);
        if(isnan(maxx) || maxx<tmp) maxx = tmp;
        tmp = f->maxCoord(1);
        if(isnan(maxy) || maxy<tmp) maxy = tmp;
        tmp = f->maxCoord(2);
        if(isnan(maxz) || maxz<tmp) maxz = tmp;
    }
    bounds b;
    b.min = vertex(minx, miny, minz);
    b.max = vertex(maxx, maxy, maxz);
    return b;
}
//Stores intersection of b1 and b2 in newbounds
void intersectBoundingBoxes(const bounds& b1, const bounds& b2, bounds& newbounds){
    newbounds.min = max3(b1.min, b2.min);
    newbounds.max = min3(b1.max, b2.max);
}
void listBounds(bounds& newbounds, int nverts, const vertex *v1 ...){
    newbounds.min = vertex(NAN,NAN,NAN);
    newbounds.max = vertex(NAN,NAN,NAN);
    const vertex *vtmp;
    
    va_list vertices;
    va_start(vertices, v1);
    vtmp = v1;
    for(int i=0; i<nverts; i++){
        if(newbounds.min!=newbounds.min || newbounds.min.x>vtmp->x) newbounds.min.x = vtmp->x;
        if(newbounds.min!=newbounds.min || newbounds.min.y>vtmp->y) newbounds.min.y = vtmp->y;
        if(newbounds.min!=newbounds.min || newbounds.min.z>vtmp->z) newbounds.min.z = vtmp->z;
        if(newbounds.max!=newbounds.max || newbounds.max.x<vtmp->x) newbounds.max.x = vtmp->x;
        if(newbounds.max!=newbounds.max || newbounds.max.y<vtmp->y) newbounds.max.y = vtmp->y;
        if(newbounds.max!=newbounds.max || newbounds.max.z<vtmp->z) newbounds.max.z = vtmp->z;
        vtmp = va_arg(vertices, vertex*);
    }
}

//geometry intersection

//Intersect a ray with a bounding box (AABB = Axis Aligned Bounding Box)
bool rayAABBIntersect(const bounds& AABB, const ray& r){
    real tmp;
    real tmin = (AABB.min.x-r.org.x) / r.dir.x;
    real tmax = (AABB.max.x-r.org.x) / r.dir.x;
    if(tmin>tmax){
        tmp = tmin;
        tmin = tmax;
        tmax = tmp;
    }
    
    real tymin = (AABB.min.y-r.org.y) / r.dir.y;
    real tymax = (AABB.max.y-r.org.y) / r.dir.y;
    if (tymin > tymax) {
        tmp = tymin;
        tymin = tymax;
        tymax = tmp;
    }
    if(tmin > tymax || tmax < tymin) return false;
    if (tymin > tmin) tmin = tymin;
    if (tymax < tmax) tmax = tymax;
    
    real tzmin = (AABB.min.z-r.org.z) / r.dir.z;
    real tzmax = (AABB.max.z-r.org.z) / r.dir.z;
    if (tzmin > tzmax) {
        tmp = tzmin;
        tzmin = tzmax;
        tzmax = tmp;
    }
    if (tmin > tzmax || tmax < tzmin) return false;
    if (tzmin > tmin) tmin = tzmin;
    if (tzmax < tmax) tmax = tzmax;
    
    if(tmin<=0 && tmax<=0) return false; //ray should not intersect bounding box if box is behind the origin of the ray
    return true;
}
//Intersect a ray with a setof faces. If lazy==false, return nearest face. Otherwise, return first intersection.
//(t,u,v) stores (intersection distance, barycentric coordinate from v12, barycentric coordiante from v13)
face *rayFacesIntersect(const std::vector<face*>& faces, const ray& r, bool lazy, vertex *tuv){
    face *f = nullptr;
    vertex tuvtmp;
    real zmin = _nan("");
    for(face *ftmp : faces){
        if(ftmp->intersectRayTriangle(r, &tuvtmp) && (tuvtmp.t < zmin || isnan(zmin))){
            zmin = tuvtmp.t;
            f = ftmp;
            if(tuv!=nullptr) *tuv = tuvtmp;
            if(lazy) return f;
        }
    }
    return f;
}
//Intersect a ray with a sphere, t stores distance from ray origin to intersect
bool raySphereIntersect(const ray& r, const real rad, real& t){
    real DdotO = dot(r.dir, r.org);
    real lensq = r.org.lensqr();
    real disc = DdotO*DdotO - (lensq-rad*rad);
    if(disc<0) return false; //no intersect
    else if (disc==0){ //one intersect
        t = -DdotO;
        return (t>=0);
    }
    else{
        //return nearest intersection, that is not behind ray origin
        disc = _sqrt(disc);
        //two intersection points.
        real t1=disc-DdotO, t2=(-disc)-DdotO;
        if(t1>=0 && t2>=0){
            t = _min(t1,t2);
            return true;
        }
        else if (t1>=0 || t2>=0){
            t=_max(t1,t2);
            return true;
        }
        else return false;
    }
}
