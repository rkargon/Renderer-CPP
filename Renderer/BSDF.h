//
//  BRDF.h
//  Renderer
//
//  Created by Raphael Kargon on 11/9/14.
//  Copyright (c) 2014 Raphael Kargon. All rights reserved.
//

#ifndef __Renderer__BSDF__
#define __Renderer__BSDF__

#include "geom.h"

class BSDF
{
public:
    virtual color getLight(color incidentColor, vertex incidentDirection, vertex normal, vertex returningDirection) = 0;
    virtual vertex getIncidentDirection(vertex normal, vertex viewDirection) = 0;
};

class DiffuseBSDF : public BSDF
{
public:
    color col;
    
    DiffuseBSDF(color c = color(1,1,1));
    virtual color getLight(color incidentColor, vertex incidentDirection, vertex normal, vertex returningDirection);
    virtual vertex getIncidentDirection(vertex normal, vertex viewDirection);
};

class EmissionBSDF : public BSDF
{
public:
    color col;
    double intensity;
    
    EmissionBSDF(color c=color(1,1,1), double i=1);
    virtual color getLight(color incidentColor, vertex incidentDirection, vertex normal, vertex returningDirection);
    virtual vertex getIncidentDirection(vertex normal, vertex viewDirection);
};

class GlossyBSDF : public BSDF
{
public:
    color col;
    double roughness;
    
    GlossyBSDF(color c=color(1,1,1), double r=0.8);
    virtual color getLight(color incidentColor, vertex incidentDirection, vertex normal, vertex returningDirection);
    virtual vertex getIncidentDirection(vertex normal, vertex viewDirection);
};

class MixBSDF : public BSDF
{
public:
    double fac;
    BSDF *mat1;
    BSDF *mat2;
    
    //whether current sample should come from material 1
    bool materialCounter;
    
    MixBSDF(double f = 0.5, BSDF *m1=nullptr, BSDF *m2=nullptr);
    virtual color getLight(color incidentColor, vertex incidentDirection, vertex normal, vertex returningDirection);
    virtual vertex getIncidentDirection(vertex normal, vertex viewDirection);
};

//A test material.
class TestBSDF : public BSDF
{
    virtual color getLight(color incidentColor, vertex incidentDirection, vertex normal, vertex returningDirection);
    virtual vertex getIncidentDirection(vertex normal, vertex viewDirection);
};

#endif /* defined(__Renderer__BSDF__) */
