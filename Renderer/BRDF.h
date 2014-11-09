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

class BRDF
{
    virtual color getLight(color incidentColor, vertex incidentDirection, vertex normal, vertex returningDirection) =0;
};

#endif /* defined(__Renderer__BRDF__) */
