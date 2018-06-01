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
#include <type_traits>

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
  material_id mat_id;
  const BSDF *bsdf;
  bool smooth;
  // TODO change these to unordered_maps if we want to modify
  std::vector<vertex> vertices; // List of *unique* vertices in mesh
  std::vector<edge> edges;      // List of *unique* edges in mesh
  std::vector<face> faces;      // List of faces in mesh
  std::vector<std::set<face_id>> face_adjacencies;
  std::vector<vertex> vertex_normals;

  mesh();
  // loads object from STL file
  // smooth off by default, materials uninitialized.
  // TODO exceptions
  mesh(std::ifstream &infile, std::string objname);

  // TODO perfect forwarding?
  // mesh(mesh &&m)
  //     : name(std::move(m.name)), origin(m.origin), mat(m.mat), bsdf(m.bsdf),
  //       smooth(m.smooth), vertices(std::move(m.vertices)),
  //       edges(std::move(m.edges)), faces(std::move(m.faces)),
  //       face_adjacencies(std::move(m.faces)),
  //       vertex_normals(std::move(m.vertex_normals)) {
  //   reset_face_obj_ptrs();
  // }

  // TODO consolidate forwarding fuckshit
  // template <typename M,
  //           typename = std::enable_if_t<
  //               std::is_same<mesh, std::remove_reference<M>>::value>>
  template <typename M>
  mesh(M &&m)
      : name(std::forward<M>(m).name), origin(m.origin), mat_id(m.mat_id),
        bsdf(m.bsdf), smooth(m.smooth), vertices(std::forward<M>(m).vertices),
        edges(std::forward<M>(m).edges), faces(std::forward<M>(m).faces),
        face_adjacencies(std::forward<M>(m).face_adjacencies),
        vertex_normals(std::forward<M>(m).vertex_normals) {
    reset_face_obj_ptrs();
  }

  mesh &operator=(mesh &m) {
    std::cout << "COPY ASSIGN" << std::endl;
    name = m.name;
    origin = m.origin;
    mat_id = m.mat_id;
    bsdf = m.bsdf;
    smooth = m.smooth;

    vertices = m.vertices;
    edges = m.edges;
    faces = m.faces;
    face_adjacencies = m.face_adjacencies;
    vertex_normals = m.vertex_normals;

    reset_face_obj_ptrs();
    return *this;
  }

  mesh &operator=(mesh &&m) {
    std::cout << "MOVE ASSIGN" << std::endl;
    name = std::move(m).name;
    origin = m.origin;
    mat_id = m.mat_id;
    bsdf = m.bsdf;
    smooth = m.smooth;

    vertices = std::move(m).vertices;
    edges = std::move(m).edges;
    faces = std::move(m).faces;
    face_adjacencies = std::move(m).face_adjacencies;
    vertex_normals = std::move(m).vertex_normals;

    reset_face_obj_ptrs();
    return *this;
  }

  // template <typename M> mesh &operator=(M &&m) {
  //   std::cout << "ASSIGN" << std::endl;
  //   name = std::forward<M>(m).name;
  //   origin = m.origin;
  //   bsdf = m.bsdf;
  //   smooth = m.smooth;

  //   vertices = std::forward<M>(m).vertices;
  //   edges = std::forward<M>(m).edges;
  //   faces = std::forward<M>(m).faces;
  //   face_adjacencies = std::forward<M>(m).face_adjacencies;
  //   vertex_normals = std::forward<M>(m).vertex_normals;

  //   reset_face_obj_ptrs();
  //   return *this;
  // }

  void project_texture(tex_projection_t proj);

  /* Mesh manipulation */
  // NOTE: These invalidate KD Trees :-(
  void move(const vertex &dv);
  void scale(const vertex &ds, const vertex &center);
  void scale_centered(const vertex &ds); // scale around a given origin

  // calculates the centroid, or arithmetic mean, of all the vertices in an
  // object.
  vertex centroid() const;

private:
  void reset_face_obj_ptrs();
};

std::ostream &operator<<(std::ostream &os, const mesh &m);

#endif
