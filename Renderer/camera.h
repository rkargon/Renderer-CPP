//
//  camera.h
//  Renderer
//
//  Represents a camera in 3D space, through which one views and renders a 3D scene.
//
//  Created by Raphael Kargon on 6/8/14.
//  Copyright (c) 2014 Raphael Kargon. All rights reserved.
//

#ifndef __Renderer__camera__
#define __Renderer__camera__

#include "geom.h"

class camera{
public:
    vertex center; //location of the camera
    vertex focus; //What the camera is focused on, can be used to center the view on an object
    vertex normal; //The normal vector pointing straight out of the camera
    vertex vert; //The vector pointing 'up' from the camera, used along with 'normal' to determine orientation
    double fov; //The field of view of the camera
    double mindist, maxdist; //Clipping distances
    bool ortho; //Whether the projection is orthographic or perspective
    
    camera();
    camera(const vertex& center, const vertex& focus, const vertex& normal, const vertex& vert, const double fov, const double mindist, const double maxdist, const bool ortho);
    
    /* projections */
    
    ray castRay(const double px, const double py, const double w, const double h) const;
    point2D<double> projectVertex(const vertex& v, double w, double h) const;
    double vertexDepth(const vertex& v) const;
    vertex viewVector(const vertex& v) const; //get vector from camera to given vertex
    double faceDepth(const face& f) const;
    
    /* manipulation of camera */
    
    //Moves the camera so that is is facing the focus point
    void centerFocus();
    void shiftFocus(const double dx, const double dy);
    void zoom(const double zoomfactor);
    vertex getImagePlaneVector(const double dx, const double dy);
    
    /* rotation */
    void rotateAxis(const vertex& axis, const double dtheta); //rotate along arbitrary axis
    void rotateLocalX(const double dtheta); //rotate along local X axis
    void rotateLocalY(const double dtheta); //rotate along local y axis
    void rotateLocalZ(const double dtheta); //rotate along local z axis
    void setGlobalRotation(const double theta, const double rho, const double psi); //set global rotation angles
    
private:
    //keeps track of image plane vectors
    vertex cx, cy;
    //recalculates image plane vectors
    void calcImageVectors(){
		cx = cross(normal, vert);
		cy = cross(cx, normal);
		cx.normalize();
		cy.normalize();
    }
    
};

#endif /* defined(__Renderer__camera__) */
