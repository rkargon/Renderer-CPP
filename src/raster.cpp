//
//  raster.cpp
//  Renderer
//
//  Created by Raphael Kargon on 4/25/15.
//  Copyright (c) 2015 Raphael Kargon. All rights reserved.
//

#include "raster.h"

raster::raster(const int w, const int h)
    : _w(w), _h(h), colbuffer(std::make_unique<int[]>(_w * _h)),
      zbuffer(std::make_unique<double[]>(_w * _h)),
      normbuffer(std::make_unique<int[]>(_w * _h)) {}

void raster::resize(const int w, const int h) {
  _w = w;
  _h = h;
  colbuffer = std::make_unique<int[]>(_w * _h);
  zbuffer = std::make_unique<double[]>(_w * _h);
  normbuffer = std::make_unique<int[]>(_w * _h);
}
