//
//  raster.cpp
//  Renderer
//
//  Created by Raphael Kargon on 4/25/15.
//  Copyright (c) 2015 Raphael Kargon. All rights reserved.
//

#include "raster.h"

raster::raster(const int w, const int h){
    this->w = w;
    this->h = h;
    colbuffer = new int[w*h];
    zbuffer = new double[w*h];
    normbuffer = new int[w*h];
}

raster::~raster(){
    delete[] colbuffer;
    delete[] zbuffer;
    delete[] normbuffer;
}

void raster::resize(const int w, const int h){
    this->w = w;
    this->h = h;
    
    delete[] colbuffer;
    delete[] zbuffer;
    delete[] normbuffer;
    
    colbuffer = new int[w*h];
    zbuffer = new double[w*h];
    normbuffer = new int[w*h];
}