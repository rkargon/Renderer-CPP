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

vertex EmissionBSDF::getIncidentDirection(vertex normal, vertex viewDirection){
    //does not actually matter, material does not reflect any incident light
    return viewDirection.reflection(normal);
}

GlossyBSDF::GlossyBSDF(color c, real r)
:col(c),
roughness(clamp(r, 0, 1)){};

color GlossyBSDF::getLight(color incidentColor, vertex incidentDirection, vertex normal, vertex returningDirection){
    return col*incidentColor;
}

vertex GlossyBSDF::getIncidentDirection(vertex normal, vertex viewDirection){
    //returns a random ray in a cone centered around the reflection of viewDirection
    //the cone's width depends on the roughness
    return viewDirection.reflection(normal);
}

color TestBSDF::getLight(color incidentColor, vertex incidentDirection, vertex normal, vertex returningDirection){
    return incidentColor;
}

vertex TestBSDF::getIncidentDirection(vertex normal, vertex viewDirection){
    vertex d = randomDirection();
    d += viewDirection.reflection(normal)*5;
    d.normalize();
    //this keeps it in the proper hemisphere
    if (dot(d, normal) * dot(viewDirection, normal) > 0) d *= -1;
    return d;
}

