//
//  raster.h
//  Renderer
//
//  Created by Raphael Kargon on 4/25/15.
//  Copyright (c) 2015 Raphael Kargon. All rights reserved.
//

#ifndef __Renderer__raster__
#define __Renderer__raster__

#include <memory>

class raster {
public:
  std::unique_ptr<int[]> colbuffer;
  std::unique_ptr<double[]> zbuffer;
  std::unique_ptr<int[]> normbuffer;

  raster(int w, int h);

  void resize(const int w, const int h);
  int width() const { return _w; }
  int height() const { return _h; }
  std::size_t size() const { return _w * _h; }
  std::size_t datasize() const {
    return size() * (sizeof(int) + sizeof(double) + sizeof(int));
  }

private:
  int _w, _h;
};

#endif /* defined(__Renderer__raster__) */
