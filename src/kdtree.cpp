//
//  kdtree.cpp
//  Renderer
//
//  Created by Raphael Kargon on 6/13/14.
//  Copyright (c) 2014 Raphael Kargon. All rights reserved.
//

#include "kdtree.h"

#include "mesh.h"

#include <iostream>

double kdtree::completion = 0;

kdtree kdtree::build_tree(const std::vector<const face *> &faces) {
  std::cout << "build kdtree faces: " << faces.size() << std::endl;
  if (faces.empty()) {
    return {};
  }
  // std::cout << faces.front() << std::endl;
  // std::cout << *(faces.front()->obj) << std::endl;
  // std::cout << faces.back() << std::endl;
  // std::cout << *(faces.back()->obj) << std::endl;

  bounds facebounds = calc_bounding_box(faces);
  std::vector<face_wrapper *> faceswrapper =
      face_wrapper::to_wrapper_list(faces);
  std::vector<planar_event *> events = build_event_list(faceswrapper);
  return kdtree(events, facebounds, 0);
}

kdtree::kdtree(std::vector<planar_event *> &events, const bounds &facebounds,
               int depth) {
  if (depth == 0) {
    completion = 0;
  }
  this->boundingbox = facebounds;
  int nfaces = 0;

  // face classifications are initially not 0 (see
  // face_wrapper::to_wrapper_list), and are not 0 when passed from parent
  // nodes. So an initial pass to reset classification values is not necessary
  for (planar_event *e : events) {
    if (e->fw->classification) {
      nfaces++;
      e->fw->classification = 0;
    }
  }

#define QUICKBUILD
#ifdef QUICKBUILD

  // exit if max recursion depth is reached or too few faces remaining
  if (depth > KDTfMAX_DEPTH || nfaces <= KDTfMIN_FACES ||
      facebounds.volume() < EPSILON) {
    // generate face list
    for (planar_event *e : events) {
      if (e->fw->classification == 0) {
        this->faces.push_back(e->fw->f);
        e->fw->classification = 1;
      }
    }
    completion += 1.0 / (1 << depth);
    return;
  }
#endif

  double cost = KDT_CTRAVERSAL + KDT_CINTERSECT * nfaces;
  double area = boundingbox.area();
  plane_data optimalplane{};
  optimalplane = find_optimal_plane(nfaces, events, boundingbox, area);

  // only split if SAH finds it would reduce cost
  if (cost <= optimalplane.split_cost) {
    // generate face list
    for (planar_event *e : events) {
      if (e->fw->classification == 0) {
        this->faces.push_back(e->fw->f);
        e->fw->classification = 1;
      }
    }
    completion += 1.0 / (1 << depth);
    return;
  }

  this->pos = optimalplane.plane_pos;
  this->axis = optimalplane.plane_axis;
  this->planarside = optimalplane.planar_side;
  bounds lowerbounds = boundingbox;
  bounds upperbounds = boundingbox;
  lowerbounds.max[axis] = pos;
  upperbounds.min[axis] = pos;

  std::vector<planar_event *> totalevents_left, totalevents_right,
      oldevents_left, oldevents_right, newevents_left, newevents_right;
  std::vector<face_wrapper *> faces_bothsides;
  for (planar_event *e : events)
    e->fw->classification = 3;
  for (planar_event *e : events) {
    if (e->axis == axis) {
      if (e->type == 0 && e->pos <= pos) {
        e->fw->classification = 2;
      } // left only
      else if (e->type == 2 && e->pos >= pos) {
        e->fw->classification = 1;
      } // right only
      else if (e->type == 1) {
        if (e->pos == pos) {
          e->fw->classification = planarside ? 1 : 2;
        } else if (e->pos > pos) {
          e->fw->classification = 1;
        } else {
          e->fw->classification = 2;
        }
      }
    }
  }
  for (planar_event *e : events) {
    if (e->fw->classification == 1) {
      oldevents_right.push_back(e);
    } else if (e->fw->classification == 2) {
      oldevents_left.push_back(e);
    } else if (e->fw->classification == 3) {
      faces_bothsides.push_back(e->fw);
      e->fw->classification = 4;
    }
  }

  for (face_wrapper *fw : faces_bothsides) {
    generate_clipped_events(fw, pos, axis, boundingbox, newevents_left,
                            newevents_right);
  }
  std::sort(newevents_left.begin(), newevents_left.end(),
            kdtree::planar_event::planar_event_comparator());
  std::sort(newevents_right.begin(), newevents_right.end(),
            kdtree::planar_event::planar_event_comparator());
  std::merge(oldevents_left.begin(), oldevents_left.end(),
             newevents_left.begin(), newevents_left.end(),
             std::back_inserter(totalevents_left),
             kdtree::planar_event::planar_event_comparator());
  std::merge(oldevents_right.begin(), oldevents_right.end(),
             newevents_right.begin(), newevents_right.end(),
             std::back_inserter(totalevents_right),
             kdtree::planar_event::planar_event_comparator());

  //    double lowerboundsarea = lowerbounds.area();
  //    double upperboundsarea = upperbounds.area();

  if (!totalevents_left.empty()) {
    this->lower = std::unique_ptr<kdtree>(
        new kdtree(totalevents_left, lowerbounds, depth + 1));
  }
  if (!totalevents_right.empty()) {
    this->upper = std::unique_ptr<kdtree>(
        new kdtree(totalevents_right, upperbounds, depth + 1));
  }

  // clean up events and face wrappers
  if (depth == 0) {
    std::vector<face_wrapper *> allfaces;
    for (planar_event *e : events) {
      if (e->fw->classification != 5) {
        e->fw->classification = 5;
        allfaces.push_back(e->fw);
      }
      delete e;
    }
    for (face_wrapper *fw : allfaces)
      delete fw;
  }
  for (planar_event *e : newevents_left)
    delete e;
  for (planar_event *e : newevents_right)
    delete e;
}

/**
 Generates splitting plane events at the 6 planes of the bounding box of each
 face. These represent candidates for a splitting plane of the kdtere. The SAH
 algorithm finds the optimal plane out of these to split the tree
 */
std::vector<kdtree::planar_event *>
kdtree::build_event_list(const std::vector<face_wrapper *> &faces) {
  std::vector<planar_event *> events;
  double minx, maxx, miny, maxy, minz, maxz;
  for (face_wrapper *fw : faces) {
    minx = fw->f->min_coord(0);
    maxx = fw->f->max_coord(0);
    miny = fw->f->min_coord(1);
    maxy = fw->f->max_coord(1);
    minz = fw->f->min_coord(2);
    maxz = fw->f->max_coord(2);

    if (minx == maxx)
      events.push_back(new planar_event(fw, minx, 0, 1));
    else {
      events.push_back(new planar_event(fw, minx, 0, 2));
      events.push_back(new planar_event(fw, maxx, 0, 0));
    }
    if (miny == maxy)
      events.push_back(new planar_event(fw, miny, 1, 1));
    else {
      events.push_back(new planar_event(fw, miny, 1, 2));
      events.push_back(new planar_event(fw, maxy, 1, 0));
    }
    if (minz == maxz)
      events.push_back(new planar_event(fw, minz, 2, 1));
    else {
      events.push_back(new planar_event(fw, minz, 2, 2));
      events.push_back(new planar_event(fw, maxz, 2, 0));
    }
  }
  std::sort(events.begin(), events.end(),
            planar_event::planar_event_comparator());
  return events;
}
/**
 * Finds optimal splitting plane in O(N) time.
 * From: http://www.eng.utah.edu/~cs6965/papers/kdtree.pdf
 *
 * @param nfaces
 *            The number of faces in the current bounding box
 * @param events
 *            A list of splitting plane candidates
 * @param boundingbox
 *            The current bounding box
 * @return A plane_data struct containing data for the optimal splitting plane
 */

kdtree::plane_data
kdtree::find_optimal_plane(int nfaces,
                           const std::vector<planar_event *> &events,
                           bounds &boundingbox, double area) {

  int Nleft[3] = {0}, Nplanar[3] = {0}, Nright[3] = {nfaces, nfaces, nfaces};
  int axis_best = 0;
  double pos_best = 0, cost_best = HUGE_VAL;
  bool p_side = 0;

  // temp variables
  int plane_axis, start_count, planar_count, end_count;
  double plane_pos;
  planar_event *e;
  std::pair<double, bool> planecostdata;

  for (int i = 0; i < events.size();) {
    e = events[i];
    plane_pos = e->pos;
    plane_axis = e->axis;
    start_count = planar_count = end_count = 0;
    while (i < events.size() && e->axis == plane_axis && e->pos == plane_pos &&
           e->type == 0) {
      end_count++;
      i++;
      e = events[i];
    }
    while (i < events.size() && e->axis == plane_axis && e->pos == plane_pos &&
           e->type == 1) {
      planar_count++;
      i++;
      e = events[i];
    }
    while (i < events.size() && e->axis == plane_axis && e->pos == plane_pos &&
           e->type == 2) {
      start_count++;
      i++;
      e = events[i];
    }

    Nplanar[plane_axis] = planar_count;
    Nright[plane_axis] -= planar_count;
    Nright[plane_axis] -= end_count;

    planecostdata =
        SAH(boundingbox, plane_pos, plane_axis, area, Nleft[plane_axis],
            Nplanar[plane_axis], Nright[plane_axis]);
    if (planecostdata.first < cost_best) {
      cost_best = planecostdata.first;
      pos_best = plane_pos;
      axis_best = plane_axis;
      p_side = planecostdata.second;
    }
    Nleft[plane_axis] += start_count;
    Nleft[plane_axis] += planar_count;
    Nplanar[plane_axis] = 0;
  }
  return plane_data{pos_best, axis_best, cost_best, p_side};
}

/**
 * Computes approximate cost of traversing a ray through a voxel with a
 * given split plane, using the Surface Area Heuristic
 *
 * @param bounds
 *            The bounds of the whole voxel
 * @param position
 *            The position of the splitting plane
 * @param axis
 *            The axis of the splitting plane
 * @param area
 *            The area of the voxel
 * @param Nl
 *            The number of faces to the left (below) the splitting plane
 * @param Np
 *            The number of faces exactly (flat) on the splitting plane
 * @param Nr
 *            The number of faces to the right (above) the splitting plane
 * @return a double array {cost, pside}. pside = 0 or 1, depending if planar
 *         faces should be grouped in the left or right voxel, respectively
 */

std::pair<double, bool> kdtree::SAH(const bounds &boundingbox, double position,
                                    int axis, double area, int Nl, int Np,
                                    int Nr) {
  bounds upperbounds = boundingbox;
  bounds lowerbounds = boundingbox;
  lowerbounds.max[axis] = position;
  upperbounds.min[axis] = position;

  // probability of a ray hitting each subvoxel, given that boundingbox is
  // already intersected.
  double lower_area = lowerbounds.area() / area;
  double upper_area = upperbounds.area() / area;
  // if bounding box has no area, or split cuts off just a plane,
  if (area == 0 || boundingbox.d(axis) == 0 || lower_area == 0 ||
      upper_area == 0)
    return std::pair<double, bool>(HUGE_VAL, false);

  double costleftbias = (lower_area * (Nl + Np) + upper_area * Nr) *
                        lambda(lower_area, upper_area, Nl + Np, Nr);
  double costrightbias =
      lower_area * Nl +
      upper_area * (Np + Nr) * lambda(lower_area, upper_area, Nl, Np + Nr);
  // planar faces go into left or right depending on which results in lowest
  // cost
  return std::pair<double, bool>(
      fmin(costleftbias, costrightbias) * KDT_CINTERSECT + KDT_CTRAVERSAL,
      costrightbias < costleftbias);
}

#define ROUGH_FACE_SPLIT
#ifdef ROUGH_FACE_SPLIT
void kdtree::generate_clipped_events(
    kdtree::face_wrapper *fw, double pos, int axis, const bounds &boundingbox,
    std::vector<planar_event *> &newevents_left,
    std::vector<planar_event *> &newevents_right) {
  bounds upperbounds = boundingbox;
  bounds lowerbounds = boundingbox;
  bounds facebound = fw->f->bounding_box();
  lowerbounds.max[axis] = pos;
  upperbounds.min[axis] = pos;
  intersect_bounding_boxes(lowerbounds, facebound, lowerbounds);
  intersect_bounding_boxes(upperbounds, facebound, upperbounds);

  for (int k = 0; k < 3; k++) {
    if (lowerbounds.min[k] == lowerbounds.max[k]) {
      newevents_left.push_back(new planar_event(fw, lowerbounds.min[k], k, 1));
    } else {
      newevents_left.push_back(new planar_event(fw, lowerbounds.min[k], k, 2));
      newevents_left.push_back(new planar_event(fw, lowerbounds.max[k], k, 0));
    }

    if (upperbounds.min[k] == upperbounds.max[k]) {
      newevents_right.push_back(new planar_event(fw, upperbounds.min[k], k, 1));
    } else {
      newevents_right.push_back(new planar_event(fw, upperbounds.min[k], k, 2));
      newevents_right.push_back(new planar_event(fw, upperbounds.max[k], k, 0));
    }
  }
}
#else
void kdtree::generate_clipped_events(
    kdtree::face_wrapper *fw, double pos, int axis, const bounds &boundingbox,
    std::vector<planar_event *> &newevents_left,
    std::vector<planar_event *> &newevents_right) {
  int sideflags = (fw->f->vertices[0]->vec[axis] < pos ? 0 : 1) |
                  (fw->f->vertices[1]->vec[axis] < pos ? 0 : 2) |
                  (fw->f->vertices[2]->vec[axis] < pos ? 0 : 4);

  // checks if sideflags is a power of two, ie there is only one vertex on the
  // right
  bool oneUpperSide = (sideflags & (sideflags - 1)) == 0;

  // whichever side has only one vertex, assign 'lonely' vertex to v1
  vertex *v1 = nullptr, *v2 = nullptr, *v3 = nullptr;
  switch (oneUpperSide ? sideflags : (~sideflags) & 7) {
  case 1:
    v1 = fw->f->vertices[0];
    v2 = fw->f->vertices[1];
    v3 = fw->f->vertices[2];
    break;
  case 2:
    v1 = fw->f->vertices[1];
    v2 = fw->f->vertices[0];
    v3 = fw->f->vertices[2];
    break;
  case 4:
    v1 = fw->f->vertices[2];
    v2 = fw->f->vertices[0];
    v3 = fw->f->vertices[1];
    break;
  }
  double r12 = (pos - v1->vec[axis]) / (v2->vec[axis] - v1->vec[axis]);
  double r13 = (pos - v1->vec[axis]) / (v3->vec[axis] - v1->vec[axis]);
  // points on split plane along triangle edges
  vertex v12cut = vertex::lerp(*v1, *v2, r12);
  vertex v13cut = vertex::lerp(*v1, *v3, r13);

  int ax1 = (axis + 2) % 3, ax2 = (axis + 1) % 3; // other axes
  // clamp values of intersecting points to the bounding box
  clamp(&v12cut[ax1], boundingbox.min[ax1], boundingbox.max[ax1]);
  clamp(&v12cut[ax2], boundingbox.min[ax2], boundingbox.max[ax2]);
  clamp(&v13cut[ax1], boundingbox.min[ax1], boundingbox.max[ax1]);
  clamp(&v13cut[ax2], boundingbox.min[ax2], boundingbox.max[ax2]);

  bounds lowerbounds, upperbounds;
  if (oneUpperSide) {
    vertex::list_bounds(upperbounds, 3, v1, &v12cut, &v13cut);
    vertex::list_bounds(lowerbounds, 4, v2, v3, &v12cut, &v13cut);
  } else {
    vertex::list_bounds(lowerbounds, 3, v1, &v12cut, &v13cut);
    vertex::list_bounds(upperbounds, 4, v2, v3, &v12cut, &v13cut);
  }
  // intersecting triangle bounding boxes with bounding box of volume
  vertex::intersect_bounding_boxes(lowerbounds, boundingbox, lowerbounds);
  vertex::intersect_bounding_boxes(upperbounds, boundingbox, upperbounds);

  for (int k = 0; k < 3; k++) {
    if (lowerbounds.min[k] == lowerbounds.max[k]) {
      newevents_left.push_back(new planar_event(fw, lowerbounds.min[k], k, 1));
    } else {
      newevents_left.push_back(new planar_event(fw, lowerbounds.min[k], k, 2));
      newevents_left.push_back(new planar_event(fw, lowerbounds.max[k], k, 0));
    }
    if (upperbounds.min[k] == upperbounds.max[k]) {
      newevents_right.push_back(new planar_event(fw, upperbounds.min[k], k, 1));
    } else {
      newevents_right.push_back(new planar_event(fw, upperbounds.min[k], k, 2));
      newevents_right.push_back(new planar_event(fw, upperbounds.max[k], k, 0));
    }
  }
}
#endif

int kdtree::face_wrapper::i = 0;
kdtree::face_wrapper::face_wrapper(const face *f, int classification)
    : f(f), classification(classification) {
  i++;
}
kdtree::face_wrapper::face_wrapper(const face *f) : f(f) { i++; }
kdtree::face_wrapper::~face_wrapper() { i--; }

std::vector<const face *> kdtree::face_wrapper::to_face_list(
    const std::vector<face_wrapper *> &wrapped_faces) {
  std::vector<const face *> faces;
  for (face_wrapper *fw : wrapped_faces)
    faces.push_back(fw->f);
  return faces;
}
std::vector<kdtree::face_wrapper *>
kdtree::face_wrapper::to_wrapper_list(const std::vector<const face *> &faces) {
  std::vector<face_wrapper *> wrapped_faces;
  for (const face *f : faces) {
    wrapped_faces.push_back(new face_wrapper(f, 1));
  }
  return wrapped_faces;
}

int kdtree::planar_event::i = 0;
kdtree::planar_event::planar_event(face_wrapper *fw, double p, int k, int type)
    : fw(fw), pos(p), axis(k), type(type) {
  i++;
}
kdtree::planar_event::~planar_event() { i--; }

bool kdtree::planar_event::planar_event_comparator::
operator()(const kdtree::planar_event *pe1, const kdtree::planar_event *pe2) {
  if (pe1->pos == pe2->pos) {
    if (pe1->axis == pe2->axis)
      return pe1->type < pe2->type;
    return pe1->axis < pe2->axis;
  }
  return pe1->pos < pe2->pos;
}

const face *kdtree::ray_tree_intersect(const kdtree *kdt, const ray &r,
                                       bool lazy, vertex *tuv, bool isRoot) {
  if (kdt == nullptr) {
    return nullptr;
  }
  vertex tuv1, tuv2;
  bool tuvvalid = (tuv != nullptr); // in case no variable to store intersection
                                    // coordinates is given
  const face *f1, *f2;

  // does ray intersect this bounding box?
  if (ray_AABB_intersect(kdt->boundingbox, r)) {
    // leaf node or nah?
    if (kdt->lower != nullptr || kdt->upper != nullptr) {
      f1 = ray_tree_intersect(kdt->lower.get(), r, lazy, &tuv1, true);
      if (lazy && f1 != nullptr) {
        if (tuvvalid) {
          *tuv = tuv1;
        }
        return f1;
      }
      f2 = ray_tree_intersect(kdt->upper.get(), r, lazy, &tuv2, true);
      if (lazy && f2 != nullptr) {
        if (tuvvalid) {
          *tuv = tuv2;
        }
        return f2;
      }

      // non-lazy compare return face from each child
      if (f1 == nullptr) {
        if (tuvvalid) {
          *tuv = tuv2;
        }
        return f2;
      } else if (f2 == nullptr) {
        if (tuvvalid) {
          *tuv = tuv1;
        }
        return f1;
      } else {
        if (tuv1[0] < tuv2[0]) {
          if (tuvvalid) {
            *tuv = tuv1;
          }
          return f1;
        } else {
          if (tuvvalid) {
            *tuv = tuv2;
          }
          return f2;
        }
      }
    } else {
      // std::cout << "intersect " << kdt->faces.size() << " faces" <<
      // std::endl;
      f1 = ray_faces_intersect(kdt->faces, r, lazy, &tuv1);
      if (tuvvalid) {
        *tuv = tuv1;
      }
      return f1;
    }
  }
  return nullptr;
}

std::vector<std::pair<vertex, vertex>> kdtree::wireframe() {
  return wireframe(true, 0);
}

std::vector<std::pair<vertex, vertex>> kdtree::wireframe(bool isRoot,
                                                         int depth) {
  std::vector<std::pair<vertex, vertex>> edges;
  if (isRoot) {
    // create corner vertices (capital letters indicate greater magnitude along
    // given axis)
    vertex xyz(boundingbox.min);
    vertex XYZ(boundingbox.max);
    vertex xyZ(xyz.x, xyz.y, XYZ.z);
    vertex xYz(xyz.x, XYZ.y, xyz.z);
    vertex xYZ(xyz.x, XYZ.y, XYZ.z);
    vertex Xyz(XYZ.x, xyz.y, xyz.z);
    vertex XyZ(XYZ.x, xyz.y, XYZ.z);
    vertex XYz(XYZ.x, XYZ.y, xyz.z);

    // from bottom corner
    edges.emplace_back(xyz, xyZ);
    edges.emplace_back(xyz, xYz);
    edges.emplace_back(xyz, Xyz);

    // from top corner
    edges.emplace_back(XYZ, XYz);
    edges.emplace_back(XYZ, XyZ);
    edges.emplace_back(XYZ, xYZ);

    // connect remaining edges
    edges.emplace_back(xyZ, xYZ);
    edges.emplace_back(xyZ, XyZ);
    edges.emplace_back(xYz, XYz);
    edges.emplace_back(xYz, xYZ);
    edges.emplace_back(Xyz, XYz);
    edges.emplace_back(Xyz, XyZ);
  }

  if (faces.empty()) {
    switch (axis) {
    case 0:
      edges.emplace_back(vertex(pos, boundingbox.min.y, boundingbox.min.z),
                         vertex(pos, boundingbox.min.y, boundingbox.max.z));
      edges.emplace_back(vertex(pos, boundingbox.min.y, boundingbox.min.z),
                         vertex(pos, boundingbox.max.y, boundingbox.min.z));
      edges.emplace_back(vertex(pos, boundingbox.min.y, boundingbox.max.z),
                         vertex(pos, boundingbox.max.y, boundingbox.max.z));
      edges.emplace_back(vertex(pos, boundingbox.max.y, boundingbox.min.z),
                         vertex(pos, boundingbox.max.y, boundingbox.max.z));
      break;
    case 1:
      edges.emplace_back(vertex(boundingbox.min.x, pos, boundingbox.min.z),
                         vertex(boundingbox.min.x, pos, boundingbox.max.z));
      edges.emplace_back(vertex(boundingbox.min.x, pos, boundingbox.min.z),
                         vertex(boundingbox.max.x, pos, boundingbox.min.z));
      edges.emplace_back(vertex(boundingbox.min.x, pos, boundingbox.max.z),
                         vertex(boundingbox.max.x, pos, boundingbox.max.z));
      edges.emplace_back(vertex(boundingbox.max.x, pos, boundingbox.min.z),
                         vertex(boundingbox.max.x, pos, boundingbox.max.z));
      break;
    case 2:
      edges.emplace_back(vertex(boundingbox.min.x, boundingbox.min.y, pos),
                         vertex(boundingbox.min.x, boundingbox.max.y, pos));
      edges.emplace_back(vertex(boundingbox.min.x, boundingbox.min.y, pos),
                         vertex(boundingbox.max.x, boundingbox.min.y, pos));
      edges.emplace_back(vertex(boundingbox.max.x, boundingbox.min.y, pos),
                         vertex(boundingbox.max.x, boundingbox.max.y, pos));
      edges.emplace_back(vertex(boundingbox.min.x, boundingbox.max.y, pos),
                         vertex(boundingbox.max.x, boundingbox.max.y, pos));
      break;
    }
    if (lower) {
      auto loweredges = lower->wireframe(false, depth + 1);
      edges.insert(edges.end(), loweredges.begin(), loweredges.end());
    }
    if (upper) {
      auto upperedges = upper->wireframe(false, depth + 1);
      edges.insert(edges.end(), upperedges.begin(), upperedges.end());
    }
  }
  return edges;
}

void kdtree::calc_stats(int &nodes, int &leaf_nodes, int &non_empty_leaf_nodes,
                        int &triangles_in_leaves, double &est_traversals,
                        double &est_leaves_visited,
                        double &est_tris_intersected, double &estcost,
                        double root_area) {
  if (root_area < 0)
    root_area = this->boundingbox.area();
  double area_ratio = this->boundingbox.area() / root_area;

  nodes++;
  est_traversals += area_ratio;
  estcost += area_ratio * 15;

  if (!this->faces.empty()) {
    // every leaf node is non empty;
    leaf_nodes++;
    non_empty_leaf_nodes++;
    triangles_in_leaves += this->faces.size();
    est_leaves_visited += area_ratio;
    est_tris_intersected += area_ratio * this->faces.size();
    estcost += area_ratio * 20;
    // constants used for estcost (15, 20) are used to match up to Ingo Wald
    // 2006 paper, to compare for debugging purposes
  }

  if (this->lower) {
    this->lower->calc_stats(nodes, leaf_nodes, non_empty_leaf_nodes,
                            triangles_in_leaves, est_traversals,
                            est_leaves_visited, est_tris_intersected, estcost,
                            root_area);
  }
  if (this->upper) {
    this->upper->calc_stats(nodes, leaf_nodes, non_empty_leaf_nodes,
                            triangles_in_leaves, est_traversals,
                            est_leaves_visited, est_tris_intersected, estcost,
                            root_area);
  }
}

void kdtree::print_stats() {
  clock_t begin = clock();
  for (int i = 1; i <= 100000; i++) {
    ray_tree_intersect(this, ray(vertex(), vertex((double)rand() / RAND_MAX,
                                                  (double)rand() / RAND_MAX,
                                                  (double)rand() / RAND_MAX)),
                       false, nullptr);
  }
  clock_t end = clock();
  std::cout << "100,000 kdt intersects: "
            << double(end - begin) / CLOCKS_PER_SEC << std::endl;

  int n = 0, ln = 0, ln_ne = 0, t_tot = 0;
  double est_t = 0, est_l = 0, est_tr = 0, est_c = 0;
  calc_stats(n, ln, ln_ne, t_tot, est_t, est_l, est_tr, est_c);
  std::cout << std::endl << "total nodes: " << n << std::endl;
  std::cout << "leaf nodes : " << ln << std::endl;
  std::cout << "non-empty leaf nodes: " << ln_ne << std::endl;
  std::cout << "avg. faces per leaf node: " << double(t_tot) / ln_ne
            << std::endl;
  std::cout << std::endl << "Average, for each ray:" << std::endl;
  std::cout << "est. num of traversals: " << est_t << std::endl;
  std::cout << "est. num of leaf nodes visited: " << est_l << std::endl;
  std::cout << "est. num of triangles intersected: " << est_tr << std::endl;
  std::cout << "est. cost (using SAH): " << est_c << std::endl;
}
