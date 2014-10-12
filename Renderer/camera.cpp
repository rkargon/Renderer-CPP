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
ray camera::castRay(const real px, const real py, const int w, const int h) const{
    real img_w = 2*tan(fov/2);
    real img_h = (img_w * h / w);
    real x= px/w - 0.5;
    real y = py/h - 0.5;
    ray castray;
    
    if(ortho){
        castray.dir = normal;
        castray.org = center + cx*(x*img_w) + cy*(-y*img_h); //-y because of java panel image coordinates
        return castray;
    }
    else{
        castray.org = center;
        castray.dir = normal + cx*(x*img_w) + cy*(-y*img_h);
        castray.dir.normalize();
        return castray;
    }
}

point2D<real> camera::projectVertex(const vertex& v, const int w, const int h) const{
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

real camera::vertexDepth(const vertex& v)const{
    vertex dv = v - center;
    if(ortho) return dot(dv, normal);
    else return dv.len() * (dot(dv, normal)<0 ? -1 : 1);
}

vertex camera::viewVector(const vertex& v) const{
    return ortho ? normal * dot(v, normal): v - center;
}

real camera::faceDepth(const face& f) const{
    return vertexDepth(f.center());
}

//manipulation of camera
void camera::centerFocus(){
    vertex camdist = focus - center;
    center = focus - normal*camdist.len();
}

void camera::shiftFocus(const real dx, const real dy){
    focus += cx*dx + cy*dy;
}

void camera::zoom(const real zoomfactor){
    center = focus + (center - focus)*zoomfactor;
}

vertex camera::getImagePlaneVector(const real dx, const real dy){
    return cx*dx + cy*dy;
}

void camera::rotateAxis(const vertex& axis, const real dtheta){
    vertex a = axis.unitvect();
    real l = a.x, m = a.y, n= a.z;
    real s=_sin(dtheta), c=_cos(dtheta);
    //axes={vert,
    //      norm}
    //columns of rotation matrix
    vertex col1(l*l*(1-c) + c, l*m*(1-c) + n*s, l*n*(1-c) - m*s);
    vertex col2(m*l*(1-c) - n*s, m*m*(1-c) + c, m*n*(1-c) + l*s);
    vertex col3(n*l*(1-c) + m*s, n*m*(1-c) - l*s, n*n*(1-c) + c);
    
    //matrix multiplication
    vert = vertex(dot(vert,col1), dot(vert,col2), dot(vert,col3));
    normal = vertex(dot(normal,col1), dot(normal,col2), dot(normal,col3));
    calcImageVectors();
}
void camera::rotateLocalX(const real dtheta){
    rotateAxis(cross(normal, vert), dtheta);
}
void camera::rotateLocalY(const real dtheta){
    rotateAxis(vert, dtheta);
}
void camera::rotateLocalZ(const real dtheta){
    rotateAxis(normal, dtheta);
}

// uses a rotation matrix with given global angles to set up axis vectors
void camera::setGlobalRotation(const real theta, const real rho, const real psi){
    real st = _sin(theta), ct = _cos(theta), sr = _sin(rho), cr=_cos(rho), sp=_sin(psi), cp=_cos(psi);
    vert = vertex(-sr, cr*sp, cr*cp);
    normal = vertex(ct*cr, ct*sr*sp-cp*st, st*sp+ct*cp*sr);
    calcImageVectors();
}

zFaceComparator::zFaceComparator(camera *cam, bool nearestFirst)
{
    this->cam = cam;
    this->nearestFirst = nearestFirst;
}
bool zFaceComparator::operator()(face *f1, face *f2)
{
    //return nearestFirst != (cam->faceDepth(*f1) > cam->faceDepth(*f2));
    if(f1==nullptr) return false;
    if(f2==nullptr) return true;
    if(nearestFirst) return (cam->faceDepth(*f1) < cam->faceDepth(*f2));
    else return (cam->faceDepth(*f1) > cam->faceDepth(*f2));
}