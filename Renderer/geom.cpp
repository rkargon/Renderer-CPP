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
vertex::vertex(double a, double b, double c) :x(a),y(b),z(c) {};

//Scales the vector such that its magnitude is one, while keeping its orientation the same.
//If the vector is 0, it is unchanged.
void vertex::normalize(){
    double len = x * x + y * y + z * z;
    if (len != 0 && len != 1) {
        len = 1.0/sqrt(len);
        x *= len;
        y *= len;
        z *= len;

    }
}
void vertex::set(const double a, const double b, const double c){x=a;y=b;z=c;}

//Calculates the refelction of this vertex across a normal vector n.
//REFL = V - N*(2(V.N))
vertex vertex::reflection(const vertex& n) const {return *this - n*(2*dot(*this,n));}

//Returns the corresponding unit vector of this vertex, but this vertex is not modified.
vertex vertex::unitvect() const{
    double len = x*x+y*y+z*z;
    if(len!=0 && len!=1){
        len = 1.0/sqrt(len);
        return vertex(x*len, y*len, z*len);
    }
    return vertex(x,y,z);
}

void vertex::to_polar(double &radius, double& theta, double& phi) const {
    radius = this->len();
    theta = acos(z/radius);
    phi = atan2(y, x);
}

void vertex::to_polar_angles(const double radius, double &theta, double &phi) const {
    theta = acos(z/radius);
    phi = atan2(y, x);
}

vertex vertex::box_fold(const double l) const{
    return 2 * clamp(*this, {-l, -l, -l}, {l,l,l}) - *this;
}

// TODO test if optimization makes a difference
double vertex::sphere_fold_ratio(const double min_radius_2, const double fixed_radius_2) const{
    double r2 = this->lensqr();
    if (r2 < min_radius_2){
        return fixed_radius_2 / min_radius_2;
    } else if (r2 < fixed_radius_2){
        return fixed_radius_2 / r2;
    } else {
        return 1;
    }
    

//    double rad = this->len();
//     branchless optimization
//    double f = (rad < r) * ((r * r) / (rad * rad)) + (rad >= r) * ((rad < 1) / rad + (rad >= 1));
//    return *this * f;
}

void vertex::operator+= (const vertex& v2){x+=v2.x; y+=v2.y; z+=v2.z;}
void vertex::operator+= (const vertex *v2){x+=v2->x; y+=v2->y; z+=v2->z;}
void vertex::operator-= (const vertex& v2){x-=v2.x; y-=v2.y; z-=v2.z;}
void vertex::operator-= (const vertex *v2){x-=v2->x; y-=v2->y; z-=v2->z;}
void vertex::operator*= (const double r){x*=r; y*=r; z*=r;} //scalar multiplication
void vertex::operator*= (const vertex& v2){x*=v2.x; y*=v2.y; z*=v2.z;} //coordinate-wise multiplication
void vertex::operator*= (const vertex *v2){x*=v2->x; y*=v2->y; z*=v2->z;}
vertex vertex::operator-() {return vertex(-x, -y, -z);}

vertex operator+ (const vertex& v1, const vertex& v2) {return vertex(v1.x+v2.x, v1.y+v2.y, v1.z+v2.z);}
vertex operator- (const vertex& v1, const vertex& v2) {return vertex(v1.x-v2.x, v1.y-v2.y, v1.z-v2.z);}
vertex operator* (const vertex& v1, const double r){return vertex(v1.x*r, v1.y*r, v1.z*r);}
vertex operator* (const double r, const vertex& v1){return vertex(v1.x*r, v1.y*r, v1.z*r);}
vertex operator* (const vertex& v1, const vertex& v2) {return vertex(v1.x*v2.x, v1.y*v2.y, v1.z*v2.z);} //coordinate-wise multiplication
vertex operator/ (const double r, const vertex& v1){return vertex(v1.x/r, v1.y/r, v1.z/r);}
vertex operator/ (const vertex& v1, const double r){return vertex(v1.x/r, v1.y/r, v1.z/r);}
bool operator== (const vertex& v1, const vertex& v2){return v1.x==v2.x && v1.y==v2.y && v1.z==v2.z;}
bool operator!= (const vertex& v1, const vertex& v2){return !(v1==v2);}
std::ostream& operator<< (std::ostream &os, const vertex &v){return os << "(" << v.x << ", " << v.y << ", " << v.z << ")";}
std::ostream& operator<< (std::ostream &os, const vertex *&v){return os << "(" << v->x << ", " << v->y << ", " << v->z << ")";}

//two vertex*'s or two vertex's are equal in these hash functions if they represent the same point in space.
std::size_t vertex_ptr_hasher::operator()(const vertex* v) const{
    size_t result = 1;
    unsigned long bits = *(long *)&v->x;
    result = 31 * result + (size_t) (bits ^ (bits >> 8*4));
    bits  = *(long *)&v->y;
    result = 31 * result + (size_t) (bits ^ (bits >> 8*4));
    bits = *(long *)&v->z;
    result = 31 * result + (size_t) (bits ^ (bits >> 8*4));
    return result;
}

std::size_t vertex_hasher::operator()(const vertex v) const{
    size_t result = 1;
    unsigned long bits = *(long *)&v.x;
    result = 31 * result + (size_t) (bits ^ (bits >> 8*4));
    bits  = *(long *)&v.y;
    result = 31 * result + (size_t) (bits ^ (bits >> 8*4));
    bits = *(long *)&v.z;
    result = 31 * result + (size_t) (bits ^ (bits >> 8*4));
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

    if(normal.lensqr() && is_perpendicular(normal)){
        normal.normalize();
    }
    else{
        normal = generate_normal(); //if normal is invalid, generate it
    }
}
bool face::is_perpendicular(const vertex &normal){
    vertex v12 = *vertices[1] - *vertices[0];
    vertex v23 = *vertices[2] - *vertices[1];
    //normal must be perpendicular to two of the edges, and must have nonzero length
    return (normal.lensqr() && !dot(v12, normal) && !dot(v23, normal));
}
//use cross product of edges to create normal. Orientation depends on order of vertices
vertex face::generate_normal(){
    vertex v12 = *vertices[1] - *vertices[0];
    vertex v23 = *vertices[2] - *vertices[1];
    vertex n = cross(v12, v23);
    return n.unitvect();
}
//the center of the face, the arithmetic mean of vertices.
vertex face::center() const{
    return (*vertices[0]+*vertices[1]+*vertices[2])*(1/3.0);
}

//tuv stores:
//t: distance from origin to intersect point
//u, v: barycentric coordinates of intersect based on edge1, edge2
//Muller-Trumbore algorithm
bool face::intersect_ray_triangle(const ray& r, vertex *tuv){
    vertex edge1, edge2, tvec, pvec, qvec;
    double det, inv_det;
    double t, u, v;

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

bounds face::bounding_box() const{
    bounds b;
    b.min = vertex(min_coord(0), min_coord(1), min_coord(2));
    b.max = vertex(max_coord(0), max_coord(1), max_coord(2));
    return b;
}


bounds::bounds() :min(0,0,0),max(0,0,0){}
bounds::bounds(const vertex& a, const vertex& b) :min(a),max(b){}

//surface area
double bounds::area() const{
    double dx = fabs(max.x - min.x);
    double dy = fabs(max.y - min.y);
    double dz = fabs(max.z - min.z);
    return 2*(dx*dy + dy*dz + dz*dx);
}

/* MESHVERTEX */
meshvertex::meshvertex(){
    x=0;
    y=0;
    z=0;
}
meshvertex::meshvertex(double a, double b, double c){
    x=a; y=b; z=c;
}
//sum the normal of each adjacent face, and normalize.
vertex meshvertex::vertex_normal(){
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
vertex edge::get_vector() {return (*v2) - (*v1);}
double edge::length() {return (*v2-*v1).len();}
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
std::size_t edge_hasher::operator()(const edge &e) const{
    size_t result1 = 1;
    unsigned long bits = *(long *)&e.v1->x;
    result1 = 31 * result1 + (size_t) (bits ^ (bits >> 8*4));
    bits = *(long *)&e.v1->y;
    result1 = 31 * result1 + (size_t) (bits ^ (bits >> 8*4));
    bits = *(long *)&e.v1->z;
    result1 = 31 * result1 + (size_t) (bits ^ (bits >> 8*4));
    size_t result2 = 1;
    bits = *(long *)&e.v2->x;
    result2 = 31 * result2 + (size_t) (bits ^ (bits >> 8*4));
    bits = *(long *)&e.v2->y;
    result2 = 31 * result2 + (size_t) (bits ^ (bits >> 8*4));
    bits = *(long *)&e.v2->z;
    result2 = 31 * result2 + (size_t) (bits ^ (bits >> 8*4));
    size_t result;
    if(result1>result2) result = 31 * (31+result1) +result2;
    else result = 31 * (31+result2) + result1;
    return result;
}

std::size_t edge_ptr_hasher::operator()(edge* const& e) const{
    size_t result1 = 1;
    unsigned long bits = *(long *)&e->v1->x;
    result1 = 31 * result1 + (size_t) (bits ^ (bits >> 8*4));
    bits = *(long *)&e->v1->y;
    result1 = 31 * result1 + (size_t) (bits ^ (bits >> 8*4));
    bits = *(long *)&e->v1->z;
    result1 = 31 * result1 + (size_t) (bits ^ (bits >> 8*4));
    size_t result2 = 1;
    bits = *(long *)&e->v2->x;
    result2 = 31 * result2 + (size_t) (bits ^ (bits >> 8*4));
    bits = *(long *)&e->v2->y;
    result2 = 31 * result2 + (size_t) (bits ^ (bits >> 8*4));
    bits = *(long *)&e->v2->z;
    result2 = 31 * result2 + (size_t) (bits ^ (bits >> 8*4));
    size_t result;
    if(result1>result2) result = 31 * (31+result1) +result2;
    else result = 31 * (31+result2) + result1;
    return result;
}
template< typename X, typename Y >
bool edge_ptr_equality::operator()(X const &lhs, Y const &rhs) const{
    return *lhs==*rhs;
}

// Rotates the given vertex about an axis, the given number of degrees
vertex rotate(const vertex& v, const vertex& axis, const double dtheta){
    vertex a = axis.unitvect();
    double l = a.x, m = a.y, n= a.z;
    double s=sin(dtheta), c=cos(dtheta);
    
    //columns of rotation matrix
    vertex col1(l*l*(1-c) + c, l*m*(1-c) + n*s, l*n*(1-c) - m*s);
    vertex col2(m*l*(1-c) - n*s, m*m*(1-c) + c, m*n*(1-c) + l*s);
    vertex col3(n*l*(1-c) + m*s, n*m*(1-c) - l*s, n*n*(1-c) + c);
    
    return {dot(v, col1), dot(v, col2), dot(v, col3)};
}

vertex abs(const vertex& v){
    return {fabs(v.x), fabs(v.y), fabs(v.z)};
}

//linearly interpolate two 3d vertices
vertex lerp(const vertex& v1, const vertex& v2, const double r){return vertex(v1.x+r*(v2.x-v1.x), v1.y+r*(v2.y-v1.y), v1.z+r*(v2.z-v1.z));}
//linearly interpolate 3 vertices based on barycentric coordinates.
vertex lerp(const vertex& v1, const vertex& v2, const vertex& v3, double w1, double w2, double w3){return v1*w1 + v2*w2 + v3*w3;}
vertex lerp(const vertex& v1, const vertex& v2, const vertex& v3, const vertex& w){return v1*w.x + v2*w.y + v3*w.z;}
vertex min(const vertex& v1, const vertex& v2){return vertex(fmin(v1.x,v2.x), fmin(v1.y,v2.y), fmin(v1.z,v2.z));}
vertex max(const vertex& v1, const vertex& v2){return vertex(fmax(v1.x,v2.x), fmax(v1.y,v2.y), fmax(v1.z,v2.z));}
vertex clamp(const vertex& v, const vertex& min_v, const vertex& max_v){return max(min(v, max_v), min_v);}

//A random unit vector.
vertex random_direction(){
    double theta = (double)rand()/RAND_MAX;
    theta *= M_PI*2;
    double z = (double)rand()/RAND_MAX;
    z = z*2 - 1;
    double w = sqrt(1-z*z);
    return vertex(w*cos(theta), w*sin(theta), z);
}

vertex from_polar(const double radius, const double theta, const double phi){
    return radius * vertex{cos(phi) * sin(theta), sin(phi)*sin(theta), cos(theta)};
}

//color functions
//Take a vertex of three values in range [0,1] and convert to 32-bit RGB
uint color_to_rgb(const color& c){
    return (unsigned char)clamp(c.r*255, 0, 255)<<16 | (unsigned char)clamp(c.g*255, 0, 255)<<8 | (unsigned char)clamp(c.b*255, 0, 255);
}
//Convert a normal vector to 32-bit RGB. FOr each coordinate, the range [-1, 1] is mapped to [0, 255]
uint normal_to_rgb(const vertex& n){
    return (unsigned char)clamp(n.x*128+128, 0, 255)<<16 | (unsigned char)clamp(n.y*128+128, 0, 255)<<8 | (unsigned char)clamp(n.z*128+128, 0, 255);
}
//Convert 32-bit RGB to vertex
color rgb_to_color(const uint rgb){
    return color(((rgb>>16)&0xff)/255.0,
                 ((rgb>>8)&0xff)/255.0,
                 (rgb&0xff)/255.0);
}
//Convert 32-bit RGB to normal vector
color rgb_to_normal(const uint n){
    return vertex(((n>>16)&0xff)/255.0 - 0.5,
                  ((n>>8)&0xff)/255.0 - 0.5,
                  (n&0xff)/255.0 - 0.5);
}

color hsv_to_rgb(const int hue, const double saturation, const double value){
    double chroma = value * saturation;
    double hue_mod = hue / 60.0;
    double x = chroma * (1 - fabs(fmod(hue_mod, 2) - 1));
    switch (int(hue_mod)) {
        case 0:
            return {chroma, x, 0};
        case 1:
            return {x, chroma, 0};
        case 2:
            return {0, chroma, x};
        case 3:
            return {0, x, chroma};
        case 4:
            return {x, 0, chroma};
        case 5:
            return {chroma, 0, x};
        default:
            return {0,0,0};
    }
}

/* bounding box related functions */

//THe bounding box of a set of faces.
bounds calc_bounding_box(const std::vector<face *>& faces){
    double minx, miny, minz, maxx, maxy, maxz, tmp;
    minx = miny = minz = maxx = maxy = maxz = nan("");
    for(face *f: faces){
        tmp = f->min_coord(0);
        if(isnan(minx) || minx>tmp) minx = tmp;
        tmp = f->min_coord(1);
        if(isnan(miny) || miny>tmp) miny = tmp;
        tmp = f->min_coord(2);
        if(isnan(minz) || minz>tmp) minz = tmp;

        tmp = f->max_coord(0);
        if(isnan(maxx) || maxx<tmp) maxx = tmp;
        tmp = f->max_coord(1);
        if(isnan(maxy) || maxy<tmp) maxy = tmp;
        tmp = f->max_coord(2);
        if(isnan(maxz) || maxz<tmp) maxz = tmp;
    }
    bounds b;
    b.min = vertex(minx, miny, minz);
    b.max = vertex(maxx, maxy, maxz);
    return b;
}
//Stores intersection of b1 and b2 in newbounds
void intersect_bounding_boxes(const bounds& b1, const bounds& b2, bounds& newbounds){
    newbounds.min = max(b1.min, b2.min);
    newbounds.max = min(b1.max, b2.max);
}
void list_bounds(bounds& newbounds, int nverts, const vertex *v1 ...){
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
bool ray_AABB_intersect(const bounds& AABB, const ray& r){
    double tmp;
    double tmin = (AABB.min.x-r.org.x) / r.dir.x;
    double tmax = (AABB.max.x-r.org.x) / r.dir.x;
    if(tmin>tmax){
        tmp = tmin;
        tmin = tmax;
        tmax = tmp;
    }

    double tymin = (AABB.min.y-r.org.y) / r.dir.y;
    double tymax = (AABB.max.y-r.org.y) / r.dir.y;
    if (tymin > tymax) {
        tmp = tymin;
        tymin = tymax;
        tymax = tmp;
    }
    if(tmin > tymax || tmax < tymin) return false;
    if (tymin > tmin) tmin = tymin;
    if (tymax < tmax) tmax = tymax;

    double tzmin = (AABB.min.z-r.org.z) / r.dir.z;
    double tzmax = (AABB.max.z-r.org.z) / r.dir.z;
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
//Intersect a ray with a set of faces. If lazy==false, return nearest face. Otherwise, return first intersection.
//(t,u,v) stores (intersection distance, barycentric coordinate from v12, barycentric coordiante from v13)
face *ray_faces_intersect(const std::vector<face*>& faces, const ray& r, bool lazy, vertex *tuv){
    face *f = nullptr;
    vertex tuvtmp;
    double zmin = nan("");
    for(face *ftmp : faces){
        if(ftmp->intersect_ray_triangle(r, &tuvtmp) && (tuvtmp.t < zmin || isnan(zmin))){
            zmin = tuvtmp.t;
            f = ftmp;
            if(tuv!=nullptr) *tuv = tuvtmp;
            if(lazy) return f;
        }
    }
    return f;
}
//Intersect a ray with a sphere, t stores distance from ray origin to intersect
bool ray_sphere_intersect(const ray& r, const double rad, double& t){
    double d_dot_o = dot(r.dir, r.org);
    double lensq = r.org.lensqr();
    double disc = d_dot_o*d_dot_o - (lensq-rad*rad);
    if(disc<0) return false; //no intersect
    else if (disc==0){ //one intersect
        t = -d_dot_o;
        return (t>=0);
    }
    else{
        //return nearest intersection, that is not behind ray origin
        disc = sqrt(disc);
        //two intersection points.
        double t1=disc-d_dot_o, t2=(-disc)-d_dot_o;
        if(t1>=0 && t2>=0){
            t = fmin(t1,t2);
            return true;
        }
        else if (t1>=0 || t2>=0){
            t=fmax(t1,t2);
            return true;
        }
        else return false;
    }
}

point_2d<double> random_point_in_unit_circle(){
    double theta = double(rand())/RAND_MAX * M_PI * 2;
    double r = double(rand())/RAND_MAX;
    return {r * cos(theta), r * sin(theta)};
}
