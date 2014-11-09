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
    virtual color getLight(color incidentColor, vertex incidentDirection, vertex normal, vertex returningDirection) =0;
};

class EmissionBSDF : BSDF
{
public:
    color col;
    real intensity;
    
    EmissionBSDF(color col=color(1,1,1), real intensity=1);
    virtual color getLight(color incidentColor, vertex incidentDirection, vertex normal, vertex returningDirection);
};

#endif /* defined(__Renderer__BRDF__) */
