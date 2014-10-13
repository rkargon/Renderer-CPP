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

// Represents the 'camera' from which a 3D scene is projected and rendered.
// Contains a position vector, and orientation vectors, as well as a focus point around which the camera may be rotated.
// The minimum and maximum clipping distance, the field of view, and whether the projection is orthographic or perspective may be changed.
class camera{
public:
    vertex center; //location of the camera
    vertex focus; //What the camera is focused on, can be used to center the view on an object
    vertex normal; //The normal vector pointing straight out of the camera
    vertex vert; //The vector pointing 'up' from the camera, used along with 'normal' to determine orientation
    real fov; //The field of view of the camera
    real mindist, maxdist; //Clipping distances
    bool ortho; //Whether the projection is orthographic or perspective
    
    camera();
    camera(const vertex& center, const vertex& focus, const vertex& normal, const vertex& vert, const real fov, const real mindist, const real maxdist, const bool ortho);
    
    /* projections */
    ray castRay(const real px, const real py, const real w, const real h) const;
    point2D<real> projectVertex(const vertex& v, real w, real h) const;
    real vertexDepth(const vertex& v) const;
    vertex viewVector(const vertex& v) const;
    real faceDepth(const face& f) const;
    
    /* manipulation of camera */
    void centerFocus();
    void shiftFocus(const real dx, const real dy);
    void zoom(const real zoomfactor);
    vertex getImagePlaneVector(const real dx, const real dy);
    
    /* rotation */
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

//sorts faces with increasing distance from camera
struct zFaceComparator{
    camera *cam;
    bool nearestFirst;
    bool operator()(face* f1, face* f2);
    
    zFaceComparator(camera *cam, bool nearestFirst);
};

#endif /* defined(__Renderer__camera__) */
