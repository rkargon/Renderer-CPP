//
//  raster.h
//  Renderer
//
//  Created by Raphael Kargon on 4/25/15.
//  Copyright (c) 2015 Raphael Kargon. All rights reserved.
//

#ifndef __Renderer__raster__
#define __Renderer__raster__

#include <stdio.h>

class raster {
public:
    int *colbuffer;
    double *zbuffer;
    int *normbuffer;
    
    raster(int w, int h);
    ~raster();
    
    void resize(const int w, const int h);
private:
    int w, h;
};

#endif /* defined(__Renderer__raster__) */
