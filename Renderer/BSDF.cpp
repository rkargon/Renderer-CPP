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

DiffuseBSDF::DiffuseBSDF(color c)
:col(c){};

color DiffuseBSDF::getLight(color incidentColor, vertex incidentDirection, vertex normal, vertex returningDirection){
    return col*incidentColor*dot(normal, incidentDirection);
}

vertex DiffuseBSDF::getIncidentDirection(vertex normal, vertex viewDirection){
    vertex d;
    real cos=1, r=0;
    int i=0;
    while(cos > r){
        i++;
        d = randomDirection();
        cos = dot(normal, d);
        r = real(rand())/real(RAND_MAX);
    }
    if (dot(d, normal) * dot(viewDirection, normal) > 0) d *= -1;
    return d;
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
    return (viewDirection.reflection(normal) + randomDirection()*roughness).unitvect();
}

MixBSDF::MixBSDF(real f, BSDF *m1, BSDF *m2)
:fac(f),
mat1(m1),
mat2(m2){
    materialCounter = true;
};

color MixBSDF::getLight(color incidentColor, vertex incidentDirection, vertex normal, vertex returningDirection){
    materialCounter = !materialCounter;
    return (materialCounter ? mat1 : mat2)->getLight(incidentColor, incidentDirection, normal, returningDirection);
}

vertex MixBSDF::getIncidentDirection(vertex normal, vertex viewDirection){
    return (materialCounter ? mat1 : mat2)->getIncidentDirection(normal, viewDirection);
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

