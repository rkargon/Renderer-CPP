//
//  camera.cpp
//  Renderer
//
//  Created by Raphael Kargon on 6/8/14.
//  Copyright (c) 2014 Raphael Kargon. All rights reserved.
//

#include "camera.h"

camera::camera() :center(vertex(-5,0,0)), focus(vertex()), normal(vertex(1,0,0)), vert(vertex(0,0,1)), fov(0.75), mindist(0.01),  maxdist(100), ortho(false) {
    calcImageVectors();
}
camera::camera(const vertex& center, const vertex& focus, const vertex& normal, const vertex& vert, const real fov, const real mindist, const real maxdist, const bool ortho){
    this->center = center;
    this->focus = focus;
    this->normal = normal;
    this->vert = vert;
    this->fov = fov;
    this->mindist = mindist;
    this->maxdist = maxdist > 0 ? maxdist : _INFINITY;
    this->ortho = ortho;
    calcImageVectors();
}

//projections

/*  
 ray camera::castRay(const real px, const real py, const int w, const int h) const
 
 Casts a ray from a given point on a screen with given dimensions.
 
 Parameters:
 px, py - real - The coordinates of the point the ray should be cast from
 w, h - real - the dimensions of the screen from which the ray is being cast.
 
 Returns:
 A normalized ray pointing into the 3D scene from the camera.
 */
ray camera::castRay(const real px, const real py, const real w, const real h) const{
    real img_w = 2*tan(fov/2);
    real img_h = (img_w * h / w);
    real x= px/w - 0.5;
    real y = py/h - 0.5;
    ray castray;
    
    if(ortho){
        castray.dir = normal;
        castray.org = center + cx*(x*img_w) + cy*(-y*img_h); //-y because of screen coordinates
        return castray;
    }
    else{
        castray.org = center;
        castray.dir = normal + cx*(x*img_w) + cy*(-y*img_h);
        castray.dir.normalize();
        return castray;
    }
}

/*
 point2D<real> camera::projectVertex(const vertex& v, const int w, const int h) const
 
 Projects a 3D vertex onto a 2D screen.
 
 Parameters:
 v - vertex - the point being projected.
 w,h - real - the dimensions of the screen being projected onto.
 
 Returns:
 point2D<real> - The x and y screen coordiantes of the projected point. The coordinates are not integers to allow for sub-pixel precision.
 */
point2D<real> camera::projectVertex(const vertex& v, const real w, const real h) const{
    vertex dv = v-center;
    real x,y;
    if(!ortho) dv.normalize();
    real dvdotnorm = dot(dv, normal);
    if(dvdotnorm<=mindist || dvdotnorm > maxdist) return point2D<real>(_nan(""), _nan(""));
    x = dot(dv, cx);
    y = dot(dv, cy);
    real px = (0.5+x/fov)*w;
    real py = (0.5-y*w/(h*fov))*h;
    return point2D<real>(px,py);
}

/*
 real camera::vertexDepth(const vertex& v)const
 
 Finds the distance from the camera to the vertex. When camera is orthographic, the distance is from the vertex to the *image plane*, not the camera's center.
 */
real camera::vertexDepth(const vertex& v)const{
    vertex dv = v - center;
    if(ortho) return dot(dv, normal);
    else return dv.len() * (dot(dv, normal)<0 ? -1 : 1);
}

/*
 vertex camera::viewVector(const vertex& v) const
 
 Returns the vector fromm the camera to a vertex. If the projection is orthographic, then the vector is the camera's normal vector, scaled appropriately.
 */
vertex camera::viewVector(const vertex& v) const{
    return ortho ? normal * dot(v, normal): v - center;
}


/*
 real camera::faceDepth(const face& f) const
 
 An approximation of the distance of from a face to the camera, using the vertexDepth of the center of the face.
 */
real camera::faceDepth(const face& f) const{
    return vertexDepth(f.center());
}

//manipulation of camera

/*
 void camera::centerFocus()
 
 Moves the camera such that it's current orientation points directly to the focus point, while not changing its distance from the focus.
 */
void camera::centerFocus(){
    vertex camdist = focus - center;
    center = focus - normal*camdist.len();
}

/*
 void camera::shiftFocus(const real dx, const real dy)
 
 Moves the focus of the camera the given dx and dy, relative to the image plane. 
 This is used to pan the camera with the mouse. 
 
 Parameters:
 dx, dy - real - the x and y distances, relative to the image plane (ie perpendicular to this->normal) ti move the focus of the camera.
 */
void camera::shiftFocus(const real dx, const real dy){
    focus += cx*dx + cy*dy;
}

/*
 void camera::zoom(const real zoomfactor)
 
 Moves the camera away or towards the focus while preserving it's orientation.
 Basically zooms in on the focus.
 
 Parameters:
 zoomfactor - real - the amount to scale to distance form the focus of the camera. i.e. how much to zoom in.
 */
void camera::zoom(const real zoomfactor){
    center = focus + (center - focus)*zoomfactor;
}

/*
 vertex camera::getImagePlaneVector(const real dx, const real dy)
 
 Returns the 3D vector (relative to the camera's center) that corresponds to the given 2D coordaintes relative to the image plane.
 
 Parameters:
 dx, dy - real - The 2D coordinates relative to the iamge plane
 
 Returns:
 vertex - The 3D vector corresponding to the point (dx, dy) relative to the image plane.
 */
vertex camera::getImagePlaneVector(const real dx, const real dy){
    return cx*dx + cy*dy;
}

/* void camera::rotateAxis(const vertex& axis, const real dtheta)
 
 Rotates the camera by a given amount (in radians) along a given axis.
 
 Parameters:
 axis - vertex - The axis around which the camera should be rotated
 dtheta - real - the amount, in radians, to be rotated
 */
void camera::rotateAxis(const vertex& axis, const real dtheta){
    vertex a = axis.unitvect();
    real l = a.x, m = a.y, n= a.z;
    real s=_sin(dtheta), c=_cos(dtheta);
    
    //columns of rotation matrix
    vertex col1(l*l*(1-c) + c, l*m*(1-c) + n*s, l*n*(1-c) - m*s);
    vertex col2(m*l*(1-c) - n*s, m*m*(1-c) + c, m*n*(1-c) + l*s);
    vertex col3(n*l*(1-c) + m*s, n*m*(1-c) - l*s, n*n*(1-c) + c);
    
    //matrix multiplication
    vert = vertex(dot(vert,col1), dot(vert,col2), dot(vert,col3));
    normal = vertex(dot(normal,col1), dot(normal,col2), dot(normal,col3));
    calcImageVectors();
}

// These functions call rotateAxis with the given dtheta and the corresponding local axes.
void camera::rotateLocalX(const real dtheta){
    rotateAxis(cross(normal, vert), dtheta);
}
void camera::rotateLocalY(const real dtheta){
    rotateAxis(vert, dtheta);
}
void camera::rotateLocalZ(const real dtheta){
    rotateAxis(normal, dtheta);
}

/* void camera::setGlobalRotation(const real theta, const real rho, const real psi)

 Uses a rotation matrix with given global angles to set up axis vectors

 Parameters:
 theta, rho, psi - real - the global angles, in radians, of the local x,y,z vectors (cx, vert, normal)
*/
void camera::setGlobalRotation(const real theta, const real rho, const real psi){
    real st = _sin(theta), ct = _cos(theta), sr = _sin(rho), cr=_cos(rho), sp=_sin(psi), cp=_cos(psi);
    vert = vertex(-sr, cr*sp, cr*cp);
    normal = vertex(ct*cr, ct*sr*sp-cp*st, st*sp+ct*cp*sr);
    calcImageVectors();
}

/* zFaceComparator::zFaceComparator(camera *cam, bool nearestFirst)
 Constructs a 'zFaceComparator' object that uses a camera and an ordering type (nearest or farthest first)
 */
zFaceComparator::zFaceComparator(camera *cam, bool nearestFirst)
{
    this->cam = cam;
    this->nearestFirst = nearestFirst;
}

/* bool zFaceComparator::operator()(face *f1, face *f2)
 
 Given two faces, determines which one is first in an ordering (specified by this->nearestFirst)
 by using their approximate distances from the camera.
 The ordering is strict ,such that if the distance of two faces is the same, the function returns false
 
 Parameters:
 f1, f2 - face* - two faces to be compared by their distance.
 */
bool zFaceComparator::operator()(face *f1, face *f2)
{
    if(f1==nullptr) return false;
    if(f2==nullptr) return true;
    if(nearestFirst) return (cam->faceDepth(*f1) < cam->faceDepth(*f2));
    else return (cam->faceDepth(*f1) > cam->faceDepth(*f2));
}