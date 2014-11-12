//
//  BRDF.h
//  Renderer
//
//  Created by Raphael Kargon on 11/9/14.
//  Copyright (c) 2014 Raphael Kargon. All rights reserved.
//

#ifndef __Renderer__BRDF__
#define __Renderer__BRDF__

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
    real intensity;
    
    EmissionBSDF(color c=color(1,1,1), real i=1);
    virtual color getLight(color incidentColor, vertex incidentDirection, vertex normal, vertex returningDirection);
    virtual vertex getIncidentDirection(vertex normal, vertex viewDirection);
};

class GlossyBSDF : public BSDF
{
public:
    color col;
    real roughness;
    
    GlossyBSDF(color c=color(1,1,1), real r=0.8);
    virtual color getLight(color incidentColor, vertex incidentDirection, vertex normal, vertex returningDirection);
    virtual vertex getIncidentDirection(vertex normal, vertex viewDirection);
};

class MixBSDF : public BSDF
{
public:
    real fac;
    //uh oh. You know some very nasty people will end up causing infinite loops :-(
    BSDF *mat1;
    BSDF *mat2;
    
    //whether current sample should come from material 1
    bool materialCounter;
    
    MixBSDF(real f = 0.5, BSDF *m1=nullptr, BSDF *m2=nullptr);
    virtual color getLight(color incidentColor, vertex incidentDirection, vertex normal, vertex returningDirection);
    virtual vertex getIncidentDirection(vertex normal, vertex viewDirection);
};

//A test material.
class TestBSDF : public BSDF
{
    virtual color getLight(color incidentColor, vertex incidentDirection, vertex normal, vertex returningDirection);
    virtual vertex getIncidentDirection(vertex normal, vertex viewDirection);
};

#endif /* defined(__Renderer__BRDF__) */
