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
camera::camera(const vertex& center, const vertex& focus, const vertex& normal, const vertex& vert, const double fov, const double mindist, const double maxdist, const bool ortho){
    this->center = center;
    this->focus = focus;
    this->normal = normal;
    this->vert = vert;
    this->fov = fov;
    this->mindist = mindist;
    this->maxdist = maxdist > 0 ? maxdist : HUGE_VAL;
    this->ortho = ortho;
    calcImageVectors();
}

/* Projections */

// Cast a ray from the screen into 3D space
ray camera::castRay(const double px, const double py, const double w, const double h) const{
    double img_w = 2*tan(fov/2);
    double img_h = (img_w * h / w);
    double x= px/w - 0.5;
    double y = py/h - 0.5;
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


// Project a vertex into 3D space onto a screen
point2D<double> camera::projectVertex(const vertex& v, const double w, const double h) const{
    vertex dv = v-center;
    double x,y;
    
    //in perspective mode, only the relative angle (ie cosine, ie dot product)
    //  of the vector matters
    if(!ortho) dv.normalize();
    double dvdotnorm = dot(dv, normal);
    if(dvdotnorm<=mindist || dvdotnorm > maxdist){
        return point2D<double>(nan(""), nan(""));
    }
    x = dot(dv, cx);
    y = dot(dv, cy);
    
    double px = (0.5+x/fov)*w;
    double py = (0.5-y*w/(h*fov))*h;
    return point2D<double>(px,py);
}

// The distance of a vertex from the camera.
double camera::vertexDepth(const vertex& v)const{
    vertex dv = v - center;
    if(ortho) return dot(dv, normal);
    else return dv.len() * (dot(dv, normal)<0 ? -1 : 1);
}

// Get vector from camera to given vertex
vertex camera::viewVector(const vertex& v) const{
    return ortho ? normal * dot(v, normal): v - center;
}

double camera::faceDepth(const face& f) const{
    return vertexDepth(f.center());
}

/* manipulation of camera */

void camera::centerFocus(){
    vertex camdist = focus - center;
    center = focus - normal*camdist.len();
}

void camera::shiftFocus(const double dx, const double dy){
    focus += cx*dx + cy*dy;
}

void camera::zoom(const double zoomfactor){
    center = focus + (center - focus)*zoomfactor;
}

vertex camera::getImagePlaneVector(const double dx, const double dy){
    return cx*dx + cy*dy;
}

/* Rotation */

void camera::rotateAxis(const vertex& axis, const double dtheta){
    vertex a = axis.unitvect();
    double l = a.x, m = a.y, n= a.z;
    double s=sin(dtheta), c=cos(dtheta);
    
    //columns of rotation matrix
    vertex col1(l*l*(1-c) + c, l*m*(1-c) + n*s, l*n*(1-c) - m*s);
    vertex col2(m*l*(1-c) - n*s, m*m*(1-c) + c, m*n*(1-c) + l*s);
    vertex col3(n*l*(1-c) + m*s, n*m*(1-c) - l*s, n*n*(1-c) + c);
    
    //matrix multiplication
    vert = vertex(dot(vert,col1), dot(vert,col2), dot(vert,col3));
    normal = vertex(dot(normal,col1), dot(normal,col2), dot(normal,col3));
    calcImageVectors();
}

void camera::rotateLocalX(const double dtheta){
    rotateAxis(cross(normal, vert), dtheta);
}
void camera::rotateLocalY(const double dtheta){
    rotateAxis(vert, dtheta);
}
void camera::rotateLocalZ(const double dtheta){
    rotateAxis(normal, dtheta);
}

void camera::setGlobalRotation(const double theta, const double rho, const double psi){
    double st = sin(theta), ct = cos(theta), sr = sin(rho), cr=cos(rho), sp=sin(psi), cp=cos(psi);
    vert = vertex(-sr, cr*sp, cr*cp);
    normal = vertex(ct*cr, ct*sr*sp-cp*st, st*sp+ct*cp*sr);
    calcImageVectors();
}