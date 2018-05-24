//
//  mesh.h
//  Renderer
//
//  Created by Raphael Kargon on 6/5/14.
//  Copyright (c) 2014 Raphael Kargon. All rights reserved.
//

#ifndef __Renderer__mesh__
#define __Renderer__mesh__

#include <fstream>
#include <iostream>
#include <set>
#include <string>

#include "BSDF.h"
#include "geom.h"
#include "scene.h"

/* Texture projections:
 Determine how a texture is mapped onto an object's faces.

TEX_PROJ_SPHERICAL
    Projects vertices onto a sphere and uses spherical coordinates as texture
coordinates.
    Each vertex is converted from cartesian to spherical coordinates which are
normalized to [0, 1]
TEX_PROJ_CUBIC
TEXT_PROJ_CYLINDRICAL
 */
enum tex_projection_t {
  TEX_PROJ_SPHERICAL,
  TEX_PROJ_CUBIC,
  TEXT_PROJ_CYLINDRICAL
};

class mesh {
public:
  std::string name;
  vertex origin;
  // TODOP multiple materials
  material *mat;
  BSDF *bsdf;
  bool smooth;
  // TODO change these to unordered_maps if we want to modify
  std::vector<vertex> vertices; // List of *unique* vertices in mesh
  std::vector<edge> edges;      // List of *unique* edges in mesh
  std::vector<face> faces;      // List of faces in mesh
  std::vector<std::set<face_id>> face_adjacencies;
  std::vector<vertex> vertex_normals;

  // loads object from STL file
  // smooth off by default, materials uninitialized.
  mesh(std::ifstream &infile, std::string objname);

  void project_texture(tex_projection_t proj);

  /* Mesh manipulation */
  // NOTE: These invalidate KD Trees :-(
  void move(const vertex &dv);
  void scale(const vertex &ds, const vertex &center);
  void scale_centered(const vertex &ds); // scale around a given origin

  // calculates the centroid, or arithmetic mean, of all the vertices in an
  // object.
  vertex centroid() const;
};

#endif
