//
//  BRDF.cpp
//  Renderer
//
//  Created by Raphael Kargon on 11/9/14.
//  Copyright (c) 2014 Raphael Kargon. All rights reserved.
//

#include "BSDF.h"

EmissionBSDF::EmissionBSDF(color c, real i)
:col(c),
intensity(i){};

color EmissionBSDF::getLight(color incidentColor, vertex incidentDirection, vertex normal, vertex returningDirection){
    return col * intensity;
}