//
//  BSDF.h
//  Renderer
//
//  Created by Raphael Kargon on 11/9/14.
//  Copyright (c) 2014 Raphael Kargon. All rights reserved.
//

#ifndef __Renderer__BSDF__
#define __Renderer__BSDF__

#include "geom.h"

#include "tiny_obj_loader/tiny_obj_loader.h"

#include <memory>

class BSDF {
public:
  virtual color getLight(color incidentColor, vertex incidentDirection,
                         vertex normal, vertex returningDirection) const = 0;
  virtual vertex getIncidentDirection(vertex normal,
                                      vertex viewDirection) const = 0;

  static std::unique_ptr<BSDF> from_obj_mat(const tinyobj::material_t &mat);
};

class DiffuseBSDF : public BSDF {
public:
  color col;

  DiffuseBSDF(color c = color(1, 1, 1));
  virtual color getLight(color incidentColor, vertex incidentDirection,
                         vertex normal, vertex returningDirection) const;
  virtual vertex getIncidentDirection(vertex normal,
                                      vertex viewDirection) const;
};

class EmissionBSDF : public BSDF {
public:
  color col;

  EmissionBSDF(color c = color(10, 10, 10));
  virtual color getLight(color incidentColor, vertex incidentDirection,
                         vertex normal, vertex returningDirection) const;
  virtual vertex getIncidentDirection(vertex normal,
                                      vertex viewDirection) const;
};

class GlossyBSDF : public BSDF {
public:
  color col;
  double roughness;

  GlossyBSDF(color c = color(1, 1, 1), double r = 0.2);
  virtual color getLight(color incidentColor, vertex incidentDirection,
                         vertex normal, vertex returningDirection) const;
  virtual vertex getIncidentDirection(vertex normal,
                                      vertex viewDirection) const;
};

class MixBSDF : public BSDF {
public:
  double fac;
  BSDF *mat1;
  BSDF *mat2;

  // whether current sample should come from material 1
  mutable bool materialCounter;

  MixBSDF(double f = 0.5, BSDF *m1 = nullptr, BSDF *m2 = nullptr);
  virtual color getLight(color incidentColor, vertex incidentDirection,
                         vertex normal, vertex returningDirection) const;
  virtual vertex getIncidentDirection(vertex normal,
                                      vertex viewDirection) const;
};
#endif /* defined(__Renderer__BSDF__) */
