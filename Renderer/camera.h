//
//  camera.h
//  Renderer
//
//  Created by Raphael Kargon on 6/8/14.
//  Copyright (c) 2014 Raphael Kargon. All rights reserved.
//

#ifndef __Renderer__camera__
#define __Renderer__camera__

#include "geom.h"

class camera{
public:
    vertex center;
    vertex focus;
    vertex normal;
    vertex vert;
    real fov;
    real mindist, maxdist;
    bool ortho;
    
    camera();
    camera(const vertex& center, const vertex& focus, const vertex& normal, const vertex& vert, const real fov, const real mindist, const real maxdist, const bool ortho);
    
    //projections
    ray castRay(const real px, const real py, const int w, const int h) const;
    point2D<real> projectVertex(const vertex& v, int w, int h) const;
    real vertexDepth(const vertex& v) const;
    vertex viewVector(const vertex& v) const;
    real faceDepth(const face& f) const;
    
    //manipulation of camera
    void centerFocus();
    void shiftFocus(const real dx, const real dy);
    void zoom(const real zoomfactor);
    vertex getImagePlaneVector(const real dx, const real dy);
    
    //rotation
    void rotateAxis(const vertex& axis, const real dtheta);
    void rotateLocalX(const real dtheta);
    void rotateLocalY(const real dtheta);
    void rotateLocalZ(const real dtheta);
    void setGlobalRotation(const real theta, const real rho, const real psi);
    
private:
    vertex cx, cy;
    void calcImageVectors(){
		cx = cross(normal, vert);
		cy = cross(cx, normal);
		cx.normalize();
		cy.normalize();
    }
    
};

//sorts faces with nearest one first.
struct zFaceComparator{
    camera *cam;
    bool nearestFirst;
    bool operator()(face* f1, face* f2);
    
    zFaceComparator(camera *cam, bool nearestFirst);
};

#endif /* defined(__Renderer__camera__) */
