//
//  mesh.cpp
//  Renderer
//
//  Created by Raphael Kargon on 6/5/14.
//  Copyright (c) 2014 Raphael Kargon. All rights reserved.
//

#include "mesh.h"

#include <glm/gtc/type_ptr.hpp>
#include <glm/gtx/io.hpp>

#include <set>
#include <unordered_map>
#include <unordered_set>

struct stl_face {
  float normal[3], vs[3][3];
  std::uint16_t attr;
};

vertex vertex_from_floats(const float *fs) {
  return vertex(fs[0], fs[1], fs[2]);
}

mesh::mesh() : origin(0), mat_id(0), bsdf(nullptr), smooth(false) {}

mesh::mesh(std::ifstream &infile, std::string objname)
    : name(objname), smooth(false) {
  std::cout << "Loading STL object \"" << name << "\"" << std::endl;
  // read header
  std::unique_ptr<char[]> header = std::make_unique<char[]>(80);
  infile.read(header.get(), 80);
  std::cout << header.get() << std::endl;

  stl_face facebuffer;
  char *fb_ptr = reinterpret_cast<char *>(&facebuffer);
  infile.read(fb_ptr, 4);

  // read face data
  std::unordered_map<vertex, vertex_id> vertices_hash;
  std::unordered_set<edge> edges_hash;
  while (!infile.read(fb_ptr, 50).eof()) {
    // load normal and vertices
    auto norm = vertex_from_floats(facebuffer.normal);
    vertex_id vids[3];
    for (int vi = 0; vi < 3; ++vi) {
      auto vtmp = vertex_from_floats(facebuffer.vs[vi]);
      auto res = vertices_hash.emplace(vtmp, vertices_hash.size());
      vids[vi] = res.first->second;
      if (vids[vi] == this->vertices.size()) {
        this->vertices.push_back(vtmp);
        this->face_adjacencies.emplace_back();
      }
    }

    // Store edges of triangle
    edges_hash.emplace(vids[0], vids[1]);
    edges_hash.emplace(vids[1], vids[2]);
    edges_hash.emplace(vids[2], vids[0]);

    // Load face with normal, vertices, and link to this object
    face f(norm, vids[0], vids[1], vids[2], this);
    this->faces.push_back(f);
    for (int vi = 0; vi < 3; ++vi) {
      this->face_adjacencies[vids[vi]].insert(this->faces.size() - 1);
    }
  }

  // TODO unnecesary?
  // Load unique vertices into this->vertices
  // vertices.reserve(vertices_hash.size());
  // for (const auto &kv_iter : vertices_hash) {
  //   vertices[kv_iter.second] = kv_iter.first;
  // }

  // load unique edges into this->edges
  edges.insert(edges.begin(), edges_hash.begin(), edges_hash.end());

  // calc vertex normals
  vertex_normals.reserve(vertices.size());
  for (vertex_id i = 0; i < vertices.size(); ++i) {
    vertex n(0, 0, 0);
    for (face_id f_id : face_adjacencies[i]) {
      n += faces[f_id].normal;
    }
    vertex_normals.push_back(glm::normalize(n));
  }

  origin = centroid();

  std::cout << vertices.size() << " vertices." << std::endl;
  std::cout << edges.size() << " edges." << std::endl;
  std::cout << faces.size() << " faces." << std::endl;
  if (infile.bad()) {
    std::cout << "Reading of " << name << " failed, creating an empty object."
              << std::endl;
    return;
  }
}

// void mesh::project_texture(tex_projection_t proj) {
//   vertex vn;
//   for (meshvertex *v : vertices) {
//     vn = v->unitvect();
//     switch (proj) {
//     case TEX_PROJ_SPHERICAL:
//       v->tex_u = std::atan2(vn.x, vn.y);
//       if (v->tex_u < 0)
//         v->tex_u += M_PI * 2;
//       v->tex_u /= (2 * M_PI);
//       v->tex_v = std::acos(vn.z) / M_PI;
//       break;

//     default:
//       break;
//     }
//   }
// }

void mesh::move(const vertex &dv) {
  for (auto &v : vertices) {
    v += dv;
  }
  origin += dv;
}

void mesh::scale(const vertex &ds, const vertex &scale_center) {
  vertex dv;
  for (auto &v : vertices) {
    dv = v - scale_center; // get vertex relative to center
    dv *= ds;              // scale vertex relative to center
    v = scale_center + dv; // return new vertex
  }
  dv = origin - scale_center;
  dv *= ds;
  origin = scale_center + dv;
}
void mesh::scale_centered(const vertex &ds) { scale(ds, origin); }

vertex mesh::centroid() const {
  vertex mean;
  for (const auto &v : vertices) {
    mean += v;
  }
  return mean * (1.0 / vertices.size());
}

void mesh::reset_face_obj_ptrs() {
  std::cout << "UPDATING FACE OBJ PTRS! (" << name << ", " << faces.size()
            << " faces)" << std::endl;
  for (auto &f : faces) {
    f.obj = this;
  }
}

std::ostream &operator<<(std::ostream &os, const mesh &m) {
  return os << "mesh: (" << &m << ") " << m.name << " [" << m.vertices.size()
            << " vertices, " << m.edges.size() << " edges, " << m.faces.size()
            << " faces]." << std::endl;
}
