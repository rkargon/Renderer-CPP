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

#include "sampling.h"

#include "tiny_obj_loader/tiny_obj_loader.h"

#include <memory>
#include <utility>

class BSDF {
public:
  virtual direction_sample
  sample_direction(const vertex &normal,
                   const vertex &outgoing_direction) const = 0;

  virtual color bsdf(const vertex &incoming_dir,
                     const vertex &outgoing_dir) const = 0;

  static std::unique_ptr<BSDF> from_obj_mat(const tinyobj::material_t &mat);
};

class DiffuseBSDF : public BSDF {
public:
  color col;

  DiffuseBSDF(color c = color(1, 1, 1));

  virtual direction_sample
  sample_direction(const vertex &normal,
                   const vertex &outgoing_direction) const override;

  virtual color bsdf(const vertex &incoming_dir,
                     const vertex &outgoing_dir) const override;
};

class EmissionBSDF : public BSDF {
public:
  color col;

  EmissionBSDF(color c = color(10, 10, 10));

  virtual direction_sample
  sample_direction(const vertex &normal,
                   const vertex &outgoing_direction) const override;

  virtual color bsdf(const vertex &incoming_dir,
                     const vertex &outgoing_dir) const override;
};

class GlossyBSDF : public BSDF {
public:
  color col;
  double roughness;

  GlossyBSDF(color c = color(1, 1, 1), double r = 0.2);

  virtual direction_sample
  sample_direction(const vertex &normal,
                   const vertex &outgoing_direction) const override;

  virtual color bsdf(const vertex &incoming_dir,
                     const vertex &outgoing_dir) const override;
};

class MixBSDF : public BSDF {
public:
  double fac;
  BSDF *mat1;
  BSDF *mat2;

  // whether current sample should come from material 1
  mutable bool materialCounter;

  MixBSDF(double f = 0.5, BSDF *m1 = nullptr, BSDF *m2 = nullptr);

  virtual direction_sample
  sample_direction(const vertex &normal,
                   const vertex &outgoing_direction) const override;

  virtual color bsdf(const vertex &incoming_dir,
                     const vertex &outgoing_dir) const override;
};
#endif /* defined(__Renderer__BSDF__) */
